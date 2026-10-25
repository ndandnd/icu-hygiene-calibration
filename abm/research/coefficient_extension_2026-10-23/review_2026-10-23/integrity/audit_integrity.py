#!/usr/bin/env python3
"""Read-only raw-batch audit independent of study runner/gate imports. No simulations."""
import sys
sys.dont_write_bytecode = True
import json, hashlib, re, platform, importlib.metadata, math, statistics, fcntl
from pathlib import Path
from decimal import Decimal, localcontext
from datetime import datetime, timezone
from collections import Counter
HERE=Path(__file__).resolve().parent
STUDY=HERE.parents[1]
ABM=STUDY.parents[1]
ARCHIVE=ABM.parent/'research/author_notes_2026-10-10/collections/bridge_20261011_061144_utc/results/bridge'

def strict_object(items):
    d={}
    for k,v in items:
        if k in d: raise ValueError('Duplicate JSON key '+k)
        d[k]=v
    return d

def read(p): return json.loads(p.read_text(), object_pairs_hook=strict_object, parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
def canonical(x): return json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def digest(x): return hashlib.sha256(canonical(x).encode()).hexdigest()
def eq(a,b): return canonical(a)==canonical(b)
checks=Counter()
def check(c, label):
    checks[label]+=1
    if not c: raise AssertionError(label)

review=(STUDY/'CODE_REVIEW.md').read_text()
pinned=dict(re.findall(r'\| `([^`]+)` \| `([0-9a-f]{64})` \|',review))
for rel,h in pinned.items(): check(sha(STUDY/rel)==h,'review_pinned_file_hash')
manifest=read(STUDY/'records/manifest.json'); science=manifest['science']
check(digest({k:v for k,v in manifest.items() if k!='identity_sha256'})==manifest['identity_sha256'],'manifest_identity')
check(digest(science)==manifest['science_sha256'],'science_hash')
for rel,h in manifest['source_sha256'].items():
    p=ARCHIVE/'part_00/manifest.json' if rel.startswith('extended_manifest') else ABM/rel
    check(sha(p)==h,'frozen_source_hash')
archmanifest=read(ARCHIVE/'part_00/manifest.json')
archived=archmanifest['resolved_config']['science']['pairs']
seedbook=read(ABM/'research/full_icu_confirmation_2026-10-08/bridge/seedbook.json')
banks={seed:bank for bank,ss in seedbook['seed_blocks'].items() for seed in ss}
order=[s for ab in zip(seedbook['seed_blocks']['A'],seedbook['seed_blocks']['B']) for s in ab]
check(len(order)==len(set(order))==256 and set(order)==set(seedbook['all_seeds']),'seedbook_coverage')
check(science['resolved_config']['seeds']==seedbook['all_seeds'],'config_seed_order')
expected=[{k:p[k] for k in ('id','h0','h1','rho','e')} for p in archived]

def pair(prefix,a,b,r,e): return {'id':f'{prefix}h{a:g}_to_{b:g}_rho{r:g}_e{e:g}','h0':float(a),'h1':float(b),'rho':float(r),'e':float(e)}
for r in (0.,.5):
    for e in (.91,.92,.93,.94,.95,.96,.98): expected.append(pair('x_cross_',.9,.95,r,e))
for a,b in ((.6,.7),(.7,.8)):
    for r,e in ((0.,1.),(.5,1.),(0.,.8),(.5,.8)): expected.append(pair('x_range_',a,b,r,e))
for a,b in ((.8,.9),(.95,.99)): expected.append(pair('x_range_',a,b,.5,.8))
for r in (0.,.5):
    with localcontext() as ctx:
        ctx.prec=100
        rd=Decimal(str(r))
        roots=[]
        for target in ('0.90','0.95'):
            # Positive solution x=1-h to (1-rho)x^2+rho*x=1-target.
            x=(-rd+(rd*rd+4*(1-rd)*(1-Decimal(target))).sqrt())/(2*(1-rd))
            roots.append(float(1-x))
    for e in (1.,.8): expected.append(dict(id=f'x_trans_c0.9_to_0.95_rho{r:g}_e{e:g}',h0=roots[0],h1=roots[1],rho=r,e=e))
for r in (0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95):
    for e in (.5,.6,.7,.75,.78,.8,.82,.85,.9,.95,.98,1.): expected.append(pair('x_fine_',.9,.95,r,e))
check(len(expected)==182,'expected_arm_count')
check(eq([{k:p[k] for k in ('id','h0','h1','rho','e')} for p in science['pairs']],expected),'ordered_arm_parameters')
ids=[p['id'] for p in expected]
check(science['kinds']['all']==ids,'ordered_arm_ids')
check(len(set(ids))==182,'unique_arm_ids')
byid={p['id']:p for p in expected}
primaries={p['legacy_id']:p['id'] for p in archived if p.get('legacy_id')}
key=lambda p: tuple(p[k] for k in ('h0','h1','rho','e'))
arch_bykey={key(p):p['id'] for p in archived}
dups={p['id']:arch_bykey[key(p)] for p in expected[22:] if key(p) in arch_bykey}
internal={p['id']:q['id'] for p in expected if p['id'].startswith('x_cross_') for q in expected if q['id'].startswith('x_fine_') and key(p)==key(q)}
check(len(dups)==8 and len(internal)==4,'duplicate_counts')
benchmark=read(STUDY/'BENCHMARK.json')
check(benchmark['decision']['surface']=='full','benchmark_full_plan')
projection=max(x['process_cpu_seconds_total'] for x in benchmark['benchmark']['seeds'])*256/3600
check(projection==benchmark['projection']['projected_cpu_hours_used_for_decision'] and projection<=60,'prospective_budget_rule')
codehashes=benchmark['sources']['code']
for rel,h in codehashes.items(): check(sha(ABM/rel)==h,'runtime_code_hash')
actual_runtime=dict(python=sys.version,implementation=platform.python_implementation(),system=platform.system(),machine=platform.machine(),**{n:importlib.metadata.version(n) for n in ('Mesa','numpy','pandas','networkx')})
check(eq(actual_runtime,manifest['runtime']),'review_matches_production_runtime')
records=STUDY/'records'
paths=list(records.glob('[0-9]*.json'))
check(len(paths)==256 and {int(p.stem) for p in paths}==set(order),'exact_record_seed_set')
check(not list(records.glob('failed_*.json')),'no_failed_seed_files')
check(not list(records.glob('.*.tmp')),'no_partial_temporary_files')
with (records/'.run.lock').open('r') as handle:
    fcntl.flock(handle.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    fcntl.flock(handle.fileno(),fcntl.LOCK_UN)
checks['runner_lock_free']+=1
exposure=('event_counts','executed_ticks','model_kwargs','observation_patient_days','observation_susceptible_patient_days','observation_admissions','observation_imported_colonizations','total_admissions','total_imported_colonizations')
histkeys=('K_by_opportunities','raw_pairs_without_reset_by_opportunities')
non_science={'cpu_seconds','process_cpu_seconds_total','wall_seconds','elapsed_seconds','load_average_start','load_average_end','completed_at_utc','pid','payload_sha256','runtime','code_sha256','started_at_utc'}
receipts=[]; cpu=[]; processcpu=[]; begins=[]; ends=[]
for seed in order:
    p=records/f'{seed}.json'; row=read(p)
    check(row['payload_sha256']==digest({k:v for k,v in row.items() if k!='payload_sha256'}),'record_payload')
    check(row['seed']==seed and row['bank']==banks[seed] and row['status']=='complete' and row['study']==STUDY.name,'record_identity')
    check(row['surface']=='all' and set(row['pairs'])==set(ids),'record_full_arm_set')
    for fld,desired in (('runtime',manifest['runtime']),('code_sha256',codehashes),('source_sha256',manifest['source_sha256']),('science_sha256',manifest['science_sha256']),('identity_sha256',manifest['identity_sha256'])): check(eq(row[fld],desired),'record_'+fld)
    check(row['executed_ticks']==378432 and not row['workflow']['alternative'] and row['workflow']['name']=='original_swap_isolation','record_workflow_horizon')
    check(set(row['profiles'])==set(primaries),'primary_profile_identity')
    oldpaths=list(ARCHIVE.glob(f'part_*/completed/{seed}.json'))
    check(len(oldpaths)==1,'unique_archived_record')
    oldpath=oldpaths[0]; old=read(oldpath)
    check(old['payload_sha256']==digest({k:v for k,v in old.items() if k!='payload_sha256'}),'archived_payload')
    check(row['ordered_bridge_events_sha256']==old['ordered_bridge_events_sha256'],'ordered_event_hash')
    check(row['gate']['archived_extended_record']['file_sha256']==sha(oldpath),'stored_archived_hash')
    for field in exposure: check(eq(row[field],old[field]),'exposure_field')
    for pdef in expected:
        pr=row['pairs'][pdef['id']]
        desired_params={k:pdef[k] for k in ('h0','h1','rho','e')}
        for legacy,primary_id in primaries.items():
            if pdef['id']==primary_id: desired_params['legacy_id']=legacy
        check(eq(pr['parameters'],desired_params),'per_arm_parameters')
        check(set(pr)=={'baseline','policy','parameters'},'per_pair_schema')
        for arm in ('baseline','policy'):
            a=pr[arm]
            check(math.isfinite(a['K']) and a['K']>=0,'finite_nonnegative_coefficient')
            check(math.isclose(a['K'],math.fsum(a['K_by_worker'].values()),rel_tol=2e-12,abs_tol=1e-10),'worker_coefficient_sum')
            ev=a['event_counts']
            check(ev['opportunities_observed']+ev['opportunities_warmup']==row['event_counts']['opportunity'],'opportunity_counter_reconciliation')
            check(ev['imported_carrier_contacts']+ev['susceptible_contacts']==row['event_counts']['patient'],'patient_counter_reconciliation')
            check(ev['implicit_resets']==row['event_counts']['reset'],'reset_counter_reconciliation')
    for pid in [p['id'] for p in archived]:
        for arm in ('baseline','policy'): check(eq(row['pairs'][pid][arm],old['pairs'][pid][arm]),'archived_arm_bitwise')
    for legacy,pid in primaries.items():
        for arm in ('baseline','policy'):
            check(eq(row['profiles'][legacy][arm],old['pairs'][pid][arm]),'primary_bitwise')
            check(eq(row['profiles'][legacy][arm],row['pairs'][pid][arm]),'primary_internal')
    for newid,oldid in dups.items():
        for arm in ('baseline','policy'):
            a,b=row['pairs'][newid][arm],old['pairs'][oldid][arm]
            if oldid in primaries.values():
                check(all(a[k] is None for k in histkeys),'duplicate_histogram_exclusion')
                a={k:v for k,v in a.items() if k not in histkeys}; b={k:v for k,v in b.items() if k not in histkeys}
            check(eq(a,b),'duplicate_archived_arm_bitwise')
    for cross,fine in internal.items():
        for arm in ('baseline','policy'): check(eq(row['pairs'][cross][arm],row['pairs'][fine][arm]),'duplicate_internal_arm_bitwise')
    g=row['gate']
    check(g['passed'] and g['counts']==dict(primary_profiles_equal=3,archived_pairs_equal=22,archived_duplicates_checked=8,archived_duplicates_equal=8,internal_duplicates_checked=4,internal_duplicates_equal=4),'stored_gate_counts')
    check(all(g['informational_exposure_fields_equal'].values()),'stored_exposure_pass')
    for group in ('primary_profiles','archived_pairs','archived_duplicates','internal_duplicates'):
        check(all(v['equal'] and all(t['equal'] for t in v['arms'].values()) for v in g[group].values()),'stored_gate_details')
    for fixture in ('benchmark','determinism'):
        fp=STUDY/'gate_checks'/fixture/'records'/f'{seed}.json'
        if fp.exists():
            fr=read(fp)
            check(eq({k:v for k,v in row.items() if k not in non_science},{k:v for k,v in fr.items() if k not in non_science}),'fixture_'+fixture+'_science_bitwise')
    cpu.append(row['cpu_seconds']);processcpu.append(row['process_cpu_seconds_total']);begins.append(row['started_at_utc']);ends.append(row['completed_at_utc'])
    receipts.append(dict(seed=seed,bank=banks[seed],record_sha256=sha(p),record_payload_sha256=row['payload_sha256'],archived_file_sha256=sha(oldpath),event_sha256=row['ordered_bridge_events_sha256']))
check('stop done=256 already_complete=0 failed=0 gate_failed=0 interrupted_running=0 not_started=0 signal=None' in (records/'batch.log').read_text(),'batch_terminal_log')
start,end=min(begins),max(ends)
result=dict(verdict='PASS',audit_time_utc=datetime.now(timezone.utc).isoformat(),checks=dict(checks),record_count=256,pair_count=182,coefficient_arms=256*182*2,banks=dict(Counter(banks.values())),runtime=actual_runtime,manifest_identity=manifest['identity_sha256'],science_sha256=manifest['science_sha256'],manifest_file_sha256=sha(records/'manifest.json'),protocol_sha256=sha(STUDY/'PROTOCOL.md'),launch_notes_sha256=sha(STUDY/'LAUNCH_NOTES.md'),code_review_sha256=sha(STUDY/'CODE_REVIEW.md'),benchmark_sha256=sha(STUDY/'BENCHMARK.json'),benchmark_projected_cpu_hours=projection,worker_measured_cpu_hours=sum(cpu)/3600,whole_process_cpu_hours=sum(processcpu)/3600,first_started_at_utc=start,last_completed_at_utc=end,elapsed_hours=(datetime.fromisoformat(end)-datetime.fromisoformat(start)).total_seconds()/3600,failed_records=0,gate_failures=0,missing_records=0,audit_script_sha256=sha(Path(__file__)),records=receipts)
(HERE/'receipt.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))
