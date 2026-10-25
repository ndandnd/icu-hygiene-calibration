#!/usr/bin/env python3
"""Build a bounded new study, preserving all earlier selections and failures."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]

def write(path, value):
    data = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n'
    if path.exists() and path.read_text() != data:
        raise ValueError('Refusing design drift: '+str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(data)

def collect(value, used):
    if isinstance(value, dict):
        for key, child in value.items():
            if key in ('seed', 'bootstrap_seed', 'base_seed') and type(child) is int:
                used.add(child)
            elif key == 'seeds' and isinstance(child, list):
                used.update(x for x in child if type(x) is int)
            collect(child, used)
    elif isinstance(value, list):
        for child in value: collect(child, used)

def main():
    old = ABM/'research/full_icu_confirmation_2026-10-08'
    config = json.loads((old/'configs/holdout_all.json').read_text())
    blocks = json.loads((old/'seedbook.json').read_text())['seed_blocks']['holdout']
    assert config['seeds'] == blocks['A']+blocks['B']
    upper_seeds = blocks['A'][:128]+blocks['B'][:128]
    used = set(config['seeds']); scanned = 0
    for path in (ABM/'research').rglob('*.json'):
        if HERE in path.parents or path.stat().st_size > 8_000_000: continue
        try: value=json.loads(path.read_text())
        except (ValueError, UnicodeError): continue
        collect(value, used); scanned += 1
    old_count=len(used)
    def seeds(module,n):
        values=[]; counter=0
        while len(values)<n:
            s=int.from_bytes(hashlib.sha256(f'{HERE.name}/{module}/{counter}'.encode()).digest()[:4],'big');counter+=1
            if s not in used: values.append(s); used.add(s)
        return values
    grid_seeds=seeds('clock_grid',128)
    policy_seeds=seeds('clock_policy',512)
    boundary_seeds=seeds('clock_boundary',512)
    def scenario(name,p,rho,e,clock,h1=.9,env=1.):
        return dict(id=name,p=p,rho=rho,e=e,h0=.9,h1=h1,opportunities=clock,environment_multiplier=env)
    grid=[]
    for i in range(9):
        p=.15+.15*i/8
        grid.append(scenario(f'clustered_grid_{i:02d}',p,.5,1.,'separate'))
    for i,p in enumerate((.30,.3125,.325,.3375,.35,.375,.40,.45,.50,.60)):
        grid.append(scenario(f'imperfect_grid_{i:02d}',p,0.,.8,'separate'))
    surface=[]
    for original in config['scenarios']:
        s=dict(original);s['id']=s['id'].replace('_target_candidate','_surfaces_off');s['environment_multiplier']=0.;surface.append(s)
    upper=dict(config['scenarios'][0]);upper['p']*=1.25;upper['id']='independent_upper'
    common={k:config[k] for k in ('warmup_days','observation_days','block_days','model_overrides')}
    modules={
        'clock_grid':(grid,grid_seeds),
        'surfaces_off':(surface,config['seeds']),
        'independent_upper':([upper],upper_seeds),
        'clock_boundary':([scenario('independent_boundary_p1',1.,0.,1.,'separate',.95)],boundary_seeds)}
    for name,(ss,values) in modules.items():
        write(HERE/'configs'/f'{name}.json',dict(common,scenarios=ss,seeds=values,description='Separate bounded follow-up after known earlier results; never changes old selections.'))
    write(HERE/'seedbook.json',dict(schema_version=1,
        clock_grid=dict(A=grid_seeds[:64],B=grid_seeds[64:]),
        clock_policy=dict(A=policy_seeds[:256],B=policy_seeds[256:]),
        clock_boundary=dict(A=boundary_seeds[:256],B=boundary_seeds[256:]),
        surfaces_off=blocks,
        independent_upper=dict(A=blocks['A'][:128],B=blocks['B'][:128]),
        audit=dict(accessible_public_json_files=scanned,unique_old_explicit_seeds=old_count,
                   excluded_json_above_bytes=8_000_000,absent_external_archives_not_certified=True,
                   reuse_is_deliberate_for_matched_ablation_and_upper_comparisons=True)))
    write(HERE/'output_scopes.json',dict(banks={name:f'results/{HERE.name}/{name}' for name in (*modules,'clock_policy')}))
    write(HERE/'design.json',dict(schema_version=1,study=HERE.name,
        planned_jobs=dict(clock_grid=2432,surfaces_off=1536,independent_upper=256,clock_boundary=512,clock_policy_maximum=1024),
        bootstrap_draws=50000,bootstrap_seed=20261011,
        interval_levels=dict(descriptive=.95,primary_six_contrast_family=1-.05/6),
        target=3.,target_band=[2.7,3.3],
        baseline_selection=dict(exact_hit_priority=True,lowest_increasing_adjacent_bracket=True,
            interpolate_p_against_pooled_rate=True,no_fallback=True,no_grid_extension=True,
            policy_outcomes_never_select_p=True),
        primary_contrasts=['clock_minus_zone_clustered_F','clock_minus_zone_imperfect_F',
            'surfaces_off_minus_full_independent_F','surfaces_off_minus_full_persistent_F',
            'surfaces_off_minus_full_imperfect_F','upper_minus_center_independent_F'],
        missing_primary_contrast_preserves_family_size=True,
        old_holdout_config_sha256=hashlib.sha256((old/'configs/holdout_all.json').read_bytes()).hexdigest(),
        old_holdout_seedbook_sha256=hashlib.sha256((old/'seedbook.json').read_bytes()).hexdigest(),
        old_matching_failures_are_unchanged=True,
        old_resampling_strata=['A_core128','A_remainder128','B_core128','B_remainder128']))
    print('Declared maximum5760 records/pairs; selection uses only2432 new baseline records.')

if __name__=='__main__': main()
