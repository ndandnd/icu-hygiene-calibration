#!/usr/bin/env python3
"""Baseline-only selection on the complete, frozen nineteen-condition grid."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
import analyze

HERE=Path(__file__).resolve().parent
ABM=HERE.parents[1]


def choose(points,target=3.):
    """Exact hit first; otherwise lowest increasing adjacent bracket. No fallback."""
    analyze.require(points and all(type(p) in (int,float) and type(r) in (int,float) and
        math.isfinite(p) and math.isfinite(r) and 0<=p<=1 and r>=0 for p,r in points),'Finite baseline points required')
    analyze.require(all(a[0]<b[0] for a,b in zip(points,points[1:])),'Strictly ordered unique probabilities required')
    hits=[dict(p=float(p),rate=float(r)) for p,r in points if r==target]
    brackets=[];nonincreasing=[]
    for a,b in zip(points,points[1:]):
        if b[1]<=a[1]:nonincreasing.append(dict(p_low=a[0],p_high=b[0],rate_low=a[1],rate_high=b[1]))
        if a[1]<target<b[1]:
            p=a[0]+(target-a[1])*(b[0]-a[0])/(b[1]-a[1])
            brackets.append(dict(p_low=a[0],p_high=b[0],rate_low=a[1],rate_high=b[1],selected_p=float(p)))
    common=dict(exact_hits=hits,increasing_brackets=brackets,nonincreasing_adjacent_segments=nonincreasing)
    if hits:return dict(common,status='selected',method='lowest_exact_hit',selected_p=hits[0]['p'])
    if brackets:return dict(common,status='selected',method='lowest_increasing_adjacent_bracket',selected_p=brackets[0]['selected_p'])
    return dict(common,status='failed',method=None,selected_p=None,reason='No exact hit or increasing adjacent target bracket; no fallback/extension.')


def policy_config(cfg,book,profiles):
    scenarios=[]
    for profile,rho,e in (('clustered',.5,1.),('imperfect',0.,.8)):
        p=profiles[profile]['selected_p']
        if p is not None:scenarios.append(dict(id=profile+'_selected_clock',p=p,rho=rho,e=e,
            h0=.9,h1=.95,opportunities='separate',environment_multiplier=1.))
    policy=dict({k:cfg[k] for k in ('warmup_days','observation_days','block_days','model_overrides')},
        scenarios=scenarios,seeds=book['clock_policy']['A']+book['clock_policy']['B'],
        description='Baseline-only nineteen-point selection fixed before selected-clock policy; no old fitting decision changed.')
    analyze.require(len(policy['seeds'])==512 and len(set(policy['seeds']))==512 and
        len(book['clock_policy']['A'])==len(book['clock_policy']['B'])==256 and
        all(type(s) is int and s>=0 for s in policy['seeds']),'Wrong frozen policy banks/seeds')
    return policy


def selection_draws(values,scenarios,groups,draws=50000,seed=20261011):
    rng=np.random.default_rng(seed);result={p:np.full(draws,np.nan) for p in ('clustered','imperfect')}
    failures={p:0 for p in result};nonincreasing={p:0 for p in result};weight_hash=analyze.hashlib.sha256()
    indices={p:[i for i,s in enumerate(scenarios) if s['id'].startswith(p+'_')] for p in result}
    for first in range(0,draws,500):
        last=min(first+500,draws);w=np.zeros((last-first,len(values)),dtype=np.int64)
        for group in groups:w[:,group]=rng.multinomial(len(group),np.full(len(group),1/len(group)),size=last-first)
        weight_hash.update(w.astype('<i8',copy=False).tobytes(order='C'))
        totals=np.einsum('bn,nck->bck',w,values)
        rates=1000*totals[:,:,0]/totals[:,:,1]
        for profile,ii in indices.items():
            for local,row in enumerate(rates):
                selected=choose([(scenarios[i]['p'],float(row[i])) for i in ii])
                if selected['nonincreasing_adjacent_segments']:nonincreasing[profile]+=1
                if selected['selected_p'] is None:failures[profile]+=1
                else:result[profile][first+local]=selected['selected_p']
    intervals={}
    for profile,values in result.items():
        valid=values[np.isfinite(values)];intervals[profile]=dict(draws=draws,failed_draws=failures[profile],
            valid_draws=len(valid),failure_fraction=failures[profile]/draws,draws_with_nonincreasing_adjacent_segments=nonincreasing[profile],
            unconditional95=analyze.interval(values,.95),conditional_on_valid_selection95=None if not len(valid) else np.quantile(valid,[.025,.975]).tolist(),
            no_fallback=True)
    return intervals,weight_hash.hexdigest()


def write_new(path,value):
    text=json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if path.exists():
        analyze.require(path.read_text()==text,'Existing selection output differs; no reselection')
        return
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.tmp')
    with temp.open('x') as stream:stream.write(text)
    temp.replace(path)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid',type=Path,required=True)
    parser.add_argument('--output-study',type=Path,default=HERE)
    args=parser.parse_args(argv)
    analyze.require(args.output_study.resolve()==HERE,'Output study must be this frozen study directory')
    from freeze import verify
    frozen=verify();cfg=analyze.read(HERE/'configs/clock_grid.json');book=analyze.read(HERE/'seedbook.json')
    analyze.require(len(cfg['scenarios'])==19 and len(cfg['seeds'])==128 and
        book['clock_grid']['A']+book['clock_grid']['B']==cfg['seeds'] and
        len(book['clock_grid']['A'])==len(book['clock_grid']['B'])==64,'Wrong fixed grid/banks')
    bank=analyze.load_bank(args.grid,cfg,frozen['identity_sha256'],'clock_grid')
    total=bank['values'].sum(axis=0);rates=1000*total[:,0]/total[:,1]
    distributions,weights=selection_draws(bank['values'],cfg['scenarios'],[list(range(64)),list(range(64,128))])
    profiles={};selected=[]
    for profile,rho,e in (('clustered',.5,1.),('imperfect',0.,.8)):
        ii=[i for i,s in enumerate(cfg['scenarios']) if s['id'].startswith(profile+'_')]
        profiles[profile]=choose([(cfg['scenarios'][i]['p'],float(rates[i])) for i in ii])
        profiles[profile]['grid_points']=[dict(scenario_id=cfg['scenarios'][i]['id'],p=cfg['scenarios'][i]['p'],
            acquired_total=int(total[i,0]),patient_days_total=float(total[i,1]),rate=float(rates[i])) for i in ii]
        profiles[profile]['selection_uncertainty']=distributions[profile]
        profiles[profile]['banks']={}
        for label,rows in (('A',range(64)),('B',range(64,128))):
            summed=bank['values'][list(rows)].sum(axis=0);rr=1000*summed[:,0]/summed[:,1]
            profiles[profile]['banks'][label]=choose([(cfg['scenarios'][i]['p'],float(rr[i])) for i in ii])
        if profiles[profile]['selected_p'] is not None:
            selected.append(profile)
    policy=policy_config(cfg,book,profiles)
    receipt=dict(schema_version=1,study=HERE.name,status='complete_baseline_selection',parent_freeze_identity=frozen['identity_sha256'],
        grid_manifest_sha256=bank['manifest_sha256'],grid_binding_sha256=bank['binding_sha256'],grid_input_hashes=bank['input_sha256'],
        selected_profiles=selected,profiles=profiles,bootstrap_draws=50000,bootstrap_seed=20261011,
        weight_stream_sha256_int64LE=weights,source_sha256=analyze.digest(Path(__file__)),
        seed_banks=dict(grid=book['clock_grid'],policy=book['clock_policy']),
        policy_outcomes_used=False,selection_scope='Only complete baseline outcomes. Conditional intervals retain valid selections and explicitly count failures; no fallback, extension, or holdout reselection.')
    verify();analyze.recheck_bank(bank)
    write_new(HERE/'selection.json',receipt)
    files={f'research/{HERE.name}/selection.json':analyze.digest(HERE/'selection.json')}
    if selected:
        write_new(HERE/'configs/clock_policy.json',policy)
        files[f'research/{HERE.name}/configs/clock_policy.json']=analyze.digest(HERE/'configs/clock_policy.json')
    else:analyze.require(not (HERE/'configs/clock_policy.json').exists(),'Policy config exists after failed selection')
    secondary=dict(schema_version=1,study=HERE.name,parent_freeze_identity=frozen['identity_sha256'],files=files,
        grid_manifest_sha256=bank['manifest_sha256'],grid_input_hashes=bank['input_sha256'])
    secondary['identity_sha256']=analyze.identity(secondary)
    write_new(HERE/'clock_policy_freeze.json',secondary)
    print(json.dumps(dict(status=receipt['status'],selected_profiles=selected,secondary_freeze_identity=secondary['identity_sha256'])))


if __name__=='__main__':main()
