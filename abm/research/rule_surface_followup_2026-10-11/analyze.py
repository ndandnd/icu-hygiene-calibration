#!/usr/bin/env python3
"""Strict complete-bank analysis of the separately frozen rule/surface study.

Period and cross-counter validators are retained from the independently reviewed
complete sensitivity audit v2 (source SHA69381043dfcd6fda4aca70b3f8e5cdf4cb31a24209c5242fd2aff0598defb000).
No statistical routine is imported from an earlier producer or auditor.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ABM=HERE.parents[1]
TOLERANCE=dict(rel_tol=2e-10,abs_tol=2e-10)
COUNTERS=('acquired_colonizations','imported_colonizations','cumul_patients','discharged_patients',
 'bed_visits','patient_contacts','environment_contacts','patient_days','susceptible_patient_days',
 'hygiene_opportunity_count','hygiene_completion_count','hygiene_effective_clear_count','hygiene_decontamination_count')
EXPOSURES={'patient_days','susceptible_patient_days','environment_contaminated_site_days'}
HYGIENE=COUNTERS[9:]
RESPONSE_KINDS=('pickups','pickups_by_completion_state','acquisitions_by_patient_source_status',
 'acquisitions_by_route','acquisitions_by_opportunities','acquisitions_by_route_and_opportunities',
 'hand_at_risk_opportunities','hand_survived_opportunities','hand_cleared_opportunities',
 'hand_implicit_reset_opportunities','environment_events')
MODEL=dict(num_HCWs=10,num_Patients=30,num_Goo=30,prob_patient_sick=.03,prob_new_patient=.003,
 isolation_factor=.6,cleaningDay=1.,isolate_sick=True,height=11,width=32,stay_extension_days=0,
 los_mode='synthetic',susceptible_los_mean_days=5.,colonized_los_mean_days=14.,
 initial_environment_contaminated=False,data_collection_interval=2592)
SOURCE_PATHS=('research/full_icu_response_2026-10-08/run.py',
 'research/full_icu_response_2026-10-08/instrumented_model.py','src/CPEmodel.py','src/agents.py','src/routes.py',
 'research/full_icu_calibration_2026-10-07/run.py','research/full_icu_response_2026-10-08/PROTOCOL.md')
PACKAGES=dict(Mesa='1.2.1',numpy='1.26.4',pandas='3.0.6',networkx='3.7')

def require(value,message):
    if not value:raise ValueError(message)

def digest(path):
    return hashlib.sha256(safe(path).read_bytes()).hexdigest()

def identity(value,key='identity_sha256'):
    return hashlib.sha256(json.dumps({k:v for k,v in value.items() if k!=key},
        sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def parse(raw):
    def pairs(items):
        result={}
        for key,value in items:
            require(key not in result,'Duplicate JSON key')
            result[key]=value
        return result
    def invalid(value):raise ValueError('Nonfinite JSON: '+value)
    return json.loads(raw,object_pairs_hook=pairs,parse_constant=invalid)

def read(path):return parse(safe(path).read_bytes())

def bound_read(path):
    raw=safe(path).read_bytes()
    return parse(raw),hashlib.sha256(raw).hexdigest()

def safe(path):
    path=Path(os.path.abspath(path))
    require(not any(p.is_symlink() for p in (path,*path.parents)),'Symlink input/output path refused')
    return path

def compare(actual,expected,where=''):
    if isinstance(expected,dict):
        require(isinstance(actual,dict) and set(actual)==set(expected),'Object fields differ: '+where)
        for key,value in expected.items():compare(actual[key],value,where+'/'+key)
    elif isinstance(expected,list):
        require(isinstance(actual,list) and len(actual)==len(expected),'Array length differs: '+where)
        for i,value in enumerate(expected):compare(actual[i],value,where+'/'+str(i))
    elif type(expected) is float:
        if where.endswith('/level'):
            require(type(actual) is float and actual==expected,'Declared interval level differs: '+where)
            return
        require(type(actual) in (int,float) and math.isfinite(actual) and math.isclose(actual,expected,**TOLERANCE),
                'Numeric mismatch: '+where)
    else:
        require(type(actual) is type(expected) and actual==expected,'Exact value/type mismatch: '+where)

def exact(actual,expected,where=''):
    """Declarations and identities have no numerical-comparison tolerance."""
    require(type(actual) is type(expected),'Exact declaration type differs: '+where)
    if isinstance(expected,dict):
        require(set(actual)==set(expected),'Exact declaration fields differ: '+where)
        for key,value in expected.items():exact(actual[key],value,where+'/'+key)
    elif isinstance(expected,list):
        require(len(actual)==len(expected),'Exact declaration length differs: '+where)
        for index,value in enumerate(expected):exact(actual[index],value,where+'/'+str(index))
    else:require(actual==expected,'Exact declaration value differs: '+where)

def flat(value,prefix=''):
    require(isinstance(value,dict),'Nested counter dictionary required')
    result={}
    for key,item in value.items():
        name=prefix+'/'+key
        if isinstance(item,dict):result.update(flat(item,name))
        else:
            require(type(item) is int and item>=0,'Nonnegative integer role/response counter required')
            result[name]=item
    return result

def reconcile(part,snapshot=False):
    """Independent cross-counter invariants, in addition to temporal sums."""
    roles=part['hygiene_counts_by_role']
    require(isinstance(roles,dict) and roles,'Hygiene role counters required')
    require(set(roles)=={'Nurse','Dr','XrayDr'},'Exact fixed worker roles required')
    for counts in roles.values():
        require(isinstance(counts,dict) and set(counts)==set(HYGIENE) and
                all(type(counts[k]) is int and counts[k]>=0 for k in HYGIENE),
                'Exact nonnegative hygiene role counter schema required')
        require(all(counts[lo]<=counts[hi] for hi,lo in zip(HYGIENE,HYGIENE[1:])),
                'Hygiene role counter ordering failed')
    for key in HYGIENE:
        require(type(part[key]) is int and part[key]>=0 and sum(c[key] for c in roles.values())==part[key],
                'Hygiene role/aggregate reconciliation failed: '+key)
    require(all(part[lo]<=part[hi] for hi,lo in zip(HYGIENE,HYGIENE[1:])),
            'Aggregate hygiene counter ordering failed')
    response=part['response_counts'];flat(response)
    require(set(response)==set(RESPONSE_KINDS),'Exact response histogram categories required')
    for kind in ('acquisitions_by_route','acquisitions_by_opportunities','acquisitions_by_route_and_opportunities'):
        require(sum(response[kind].values())==part['acquired_colonizations'],
                'Acquisition histogram/aggregate reconciliation failed')
    routes=response['acquisitions_by_route']
    require(routes.get('untracked_hand',0)==0,'Untracked acquired route violates frozen instrumentation')
    require(sum(response['acquisitions_by_patient_source_status'].values())==routes.get('patient_pickup',0),
            'Patient-source acquisition histogram reconciliation failed')
    require(sum(response['pickups'].values())==sum(response['pickups_by_completion_state'].values()),
            'Pickup histogram reconciliation failed')
    for gap in set().union(*(set(response[k]) for k in
            ('hand_at_risk_opportunities','hand_survived_opportunities','hand_cleared_opportunities'))):
        require(response['hand_at_risk_opportunities'].get(gap,0)==
                response['hand_survived_opportunities'].get(gap,0)+response['hand_cleared_opportunities'].get(gap,0),
                'Hand-at-risk histogram reconciliation failed')
    environment=response['environment_events']
    require(environment.get('dirty_scheduled_clear',0)<=environment.get('scheduled_clear',0),
            'Environmental clear histogram ordering failed')
    if snapshot:
        require(part['cumul_patients']-part['discharged_patients']==part['census'],
                'Snapshot patient admission/discharge reconciliation failed')
    else:
        require(part['cumul_patients']-part['discharged_patients']==part['census_end']-part['census_start'],
                'Period patient admission/discharge reconciliation failed')

def period(total,blocks,start,end,n):
    reconcile(start,True);reconcile(end,True);reconcile(total)
    for block in blocks:reconcile(block)
    require(len(blocks)==n and type(total['ticks']) is int and total['ticks']==2592*n and
            type(total['days']) in (int,float) and total['days']==n,'Wrong complete period')
    require(type(start['total_ticks']) is type(end['total_ticks']) is int,'Integer cumulative ticks required')
    require(end['total_ticks']-start['total_ticks']==2592*n,'Wrong cumulative horizon')
    for key in COUNTERS+('environment_contaminated_site_days',):
        before,after=start[key],end[key];target=after-before
        require(type(before) in (int,float) and type(after) in (int,float) and
                math.isfinite(before) and math.isfinite(after) and before>=0 and after>=0 and target>=-1e-9,
                'Invalid/decreasing cumulative counter')
        values=[b[key] for b in blocks]
        if key in EXPOSURES:
            require(all(type(v) in (int,float) and math.isfinite(v) and v>=-1e-9 for v in values),
                    'Invalid daily exposure')
            require(type(total[key]) in (int,float) and math.isfinite(total[key]) and
                    math.isclose(total[key],target,**TOLERANCE) and
                    math.isclose(math.fsum(values),target,**TOLERANCE),'Exposure conservation failed: '+key)
        else:
            require(type(before) is type(after) is int and all(type(v) is int and v>=0 for v in values) and
                    type(total[key]) is int and total[key]==target==sum(values),
                    'Integer counter conservation failed: '+key)
    for key in ('hygiene_counts_by_role','response_counts'):
        before,after=flat(start[key]),flat(end[key]);summed=flat(total[key]);daily=[flat(b[key]) for b in blocks]
        keys=set(before)|set(after)|set(summed)|set().union(*(set(b) for b in daily))
        for name in keys:
            target=after.get(name,0)-before.get(name,0)
            require(target>=0 and summed.get(name,0)==target==sum(b.get(name,0) for b in daily),
                    'Nested counter conservation failed: '+key+name)
    for day,b in enumerate(blocks):
        require(type(b['ticks']) is int and b['ticks']==2592 and type(b['days']) in (int,float) and b['days']==1 and
                type(b['day_start']) in (int,float) and type(b['day_end']) in (int,float) and
                b['day_start']==day and b['day_end']==day+1,
                'Daily time/grid identity differs')
        require(-1e-9<=b['susceptible_patient_days']<=b['patient_days']+1e-9 and b['patient_days']<=30+1e-9,
                'Daily census exposure invalid')
        require(b['environment_contaminated_site_days']<=30+1e-9,'Daily environmental exposure invalid')
    endpoint_keys={'census':'census','colonized_census':'colonized_census',
                   'environment_contaminated_sites_endpoint':'environment_contaminated_sites_endpoint'}
    for base,field in endpoint_keys.items():
        require(type(start[base]) is type(end[base]) is int and 0<=start[base]<=30 and 0<=end[base]<=30,
                'Invalid census/environment snapshot endpoint')
        chain=[start[base]]
        for b in blocks:
            require(type(b[field+'_start']) is type(b[field+'_end']) is int and
                    b[field+'_start']==chain[-1] and 0<=b[field+'_end']<=30,
                    'Daily endpoint continuity differs')
            chain.append(b[field+'_end'])
        require(type(total[field+'_start']) is type(total[field+'_end']) is int and
                chain[-1]==end[base] and total[field+'_start']==start[base] and total[field+'_end']==end[base],
                'Period snapshot endpoint differs')
    for part in [total,*blocks]:
        require(part['colonized_census_start']<=part['census_start'] and
                part['colonized_census_end']<=part['census_end'],'Colonized census exceeds total census')
        days=part['days'];patient=part['patient_days'];susceptible=part['susceptible_patient_days']
        derived=dict(mean_census=patient/days,mean_susceptible_census=susceptible/days,
                     colonized_patient_days=patient-susceptible,mean_colonized_census=(patient-susceptible)/days,
                     mean_environment_contaminated_fraction=part['environment_contaminated_site_days']/(30*days))
        for denominator in ('patient_days','susceptible_patient_days'):
            derived['acquired_per_1000_'+denominator]=1000*part['acquired_colonizations']/part[denominator] if part[denominator] else None
        compare({k:part[k] for k in derived},derived,'derived_daily_or_period')
    require(total['patient_days']>0 and
            math.isclose(total['acquired_per_1000_patient_days'],
                         1000*total['acquired_colonizations']/total['patient_days'],**TOLERANCE),
            'Period exposure/rate differs')


def relative(root,name):
    p=Path(name)
    require(isinstance(name,str) and name and not p.is_absolute() and '..' not in p.parts and
            p.as_posix()==name,'Canonical relative path required')
    q=safe(root/p)
    require(q.is_relative_to(safe(root)),'Escaped root')
    return q


def validate_record(row,job,software):
    exact({k:row[k] for k in ('job_id','scenario','seed')},job,'job')
    sc=job['scenario'];shared=sc['h0']==sc['h1']
    require(row['status']=='complete' and set(row['arms'])=={'baseline','policy'} and
            type(row['distinct_observation_trajectories']) is int and
            row['distinct_observation_trajectories']==(1 if shared else 2) and
            type(row['warmup_executions']) is int and row['warmup_executions']==1 and
            type(row['executed_ticks_per_trajectory']) is int and row['executed_ticks_per_trajectory']==378432,
            'Execution/status differs')
    require(set(row['execution'])==set(software)|{'platform'} and
            {k:v for k,v in row['execution'].items() if k!='platform'}==software and
            isinstance(row['execution']['platform'],str) and row['execution']['platform'],
            'Record runtime differs from own manifest')
    exact(row['model_kwargs'],dict(MODEL,seed=job['seed'],prob_transmission=sc['p'],
        icu_hcw_wash_rate=sc['h0'],outside_hcw_wash_rate=sc['h0'],hygiene_opportunities=sc['opportunities'],
        hygiene_autocorrelation=sc['rho'],hygiene_efficacy=sc['e'],environment_multiplier=sc['environment_multiplier']),
        'model_kwargs')
    start,warm=row['initial_snapshot'],row['warmup_snapshot']
    require(type(start['total_ticks']) is type(warm['total_ticks']) is int and
            start['total_ticks']==0 and warm['total_ticks']==233280,'Wrong warmup ticks')
    period(row['warmup'],row['warmup_blocks'],start,warm,90)
    values=[];windows=[]
    if shared:exact(row['arms']['policy'],{'shared_with':'baseline'},'baseline-only policy')
    for name,h in (('baseline',sc['h0']),('policy',sc['h1'])):
        arm=row['arms']['baseline' if shared else name]
        require(isinstance(arm,dict) and set(arm)=={'h','observation','observation_blocks','final_snapshot'},
                'Complete branch required')
        require(arm['h']==h and type(arm['final_snapshot']['total_ticks']) is int and
                arm['final_snapshot']['total_ticks']==378432,'Branch completion/horizon differs')
        period(arm['observation'],arm['observation_blocks'],warm,arm['final_snapshot'],56)
        values.extend((arm['observation']['acquired_colonizations'],arm['observation']['patient_days']))
        blocks=arm['observation_blocks']
        windows.extend((sum(b['acquired_colonizations'] for b in blocks[:28]),math.fsum(b['patient_days'] for b in blocks[:28]),
                        sum(b['acquired_colonizations'] for b in blocks[28:]),math.fsum(b['patient_days'] for b in blocks[28:])))
    require(windows[1]>0 and windows[3]>0 and windows[5]>0 and windows[7]>0,'Positive observation-half exposure required')
    return values,windows


def inventory(folder,count):
    folder=safe(folder);directory=safe(folder/'completed')
    require(directory.is_dir(),'Complete bank required; outcomes were not read')
    expected={f'job_{i:08d}.json' for i in range(count)}
    paths=list(directory.iterdir())
    require(len(paths)==count and {p.name for p in paths}==expected and
            all(not p.is_symlink() and p.is_file() for p in paths),'Exact complete record inventory required; outcomes were not read')
    return sorted(paths)


def load_bank(folder,cfg,freeze_identity,module,repo=ABM):
    seeds,scenarios=cfg['seeds'],cfg['scenarios'];n=len(seeds)
    require(n==len(set(seeds)) and all(type(s) is int and s>=0 for s in seeds),'Unique integer seeds required')
    paths=inventory(folder,n*len(scenarios))
    manifest,msha=bound_read(folder/'manifest.json');binding,bsha=bound_read(folder/'study_binding.json')
    require(type(manifest['schema_version']) is int and manifest['schema_version']==1 and
            manifest['identity']==identity(manifest,'identity'),'Wrong production manifest identity')
    exact(manifest['config'],cfg,'config')
    exact(manifest['resolved_config'],dict(cfg,warmup_ticks=233280,observation_ticks=145152,block_ticks=2592),'resolved config')
    jobs=[dict(job_id=i,scenario=sc,seed=s) for i,(sc,s) in enumerate((sc,s) for sc in scenarios for s in seeds)]
    exact(manifest['jobs'],jobs,'jobs')
    sources={p:digest(relative(repo,p)) for p in SOURCE_PATHS}
    sources['array_worker.py']=digest(relative(repo,'research/full_icu_response_2026-10-08/array_worker.py'))
    exact(manifest['sources'],sources,'producer sources')
    require(type(manifest['task_count']) is int and manifest['task_count']>0,'Partition integer required')
    software=manifest['software']
    require(set(software)==set(PACKAGES)|{'python'} and all(software[k]==v for k,v in PACKAGES.items()) and
            isinstance(software['python'],str) and software['python'].split()[0]=='3.12.3','Pinned production releases required')
    exact(binding,dict(study=HERE.name,freeze_identity=freeze_identity,
        config_path=f'research/{HERE.name}/configs/{module}.json',config_sha256=digest(HERE/'configs'/f'{module}.json'),
        output_scope=f'results/{HERE.name}/{module}',task_count=manifest['task_count']),'study binding')
    values=np.zeros((n,len(scenarios),4));windows=np.zeros((n,len(scenarios),8));hashes={};seen=set()
    for path in paths:
        row,h=bound_read(path);j=row['job_id']
        require(type(j) is int and j not in seen and 0<=j<len(jobs) and path.name==f'job_{j:08d}.json' and
                row['identity']==manifest['identity'],'Wrong/duplicate record identity')
        ci,si=divmod(j,n);values[si,ci],windows[si,ci]=validate_record(row,jobs[j],software)
        hashes['completed/'+path.name]=h;seen.add(j)
    require(seen==set(range(len(jobs))) and np.isfinite(values).all() and (values[:,:,[1,3]]>0).all(),'Complete finite vectors required')
    return dict(config=cfg,manifest=manifest,manifest_sha256=msha,binding_sha256=bsha,
                input_sha256=hashes,values=values,windows=windows,root=folder)


def recheck_bank(bank):
    require(digest(bank['root']/'manifest.json')==bank['manifest_sha256'] and
            digest(bank['root']/'study_binding.json')==bank['binding_sha256'],'Bank metadata changed during analysis')
    inventory(bank['root'],len(bank['input_sha256']))
    for name,h in bank['input_sha256'].items():require(digest(relative(bank['root'],name))==h,'Record changed during analysis')


def load_old(repo,index_path):
    index,ih=bound_read(index_path)
    require(type(index['schema_version']) is int and index['schema_version']==1 and
            index['kind']=='verified_original_holdout_reuse_index' and index['status']=='pass' and
            index['identity_sha256']==identity(index),'Invalid frozen old-input index')
    require(index['n_records']==1536 and index['n_seeds_per_profile']==512,'Wrong old bank counts')
    ledger={};metadata={}
    for key in ('config','seedbook','report','freeze'):
        descriptor=index[key];p=relative(repo,descriptor['path']);require(digest(p)==descriptor['sha256'],'Old descriptor changed')
        metadata[key]=read(p);ledger[descriptor['path']]=descriptor['sha256']
    cfg=metadata['config'];seeds=cfg['seeds'];scenarios=cfg['scenarios'];book=index['ordered_statistical_banks']
    require(book['A']+book['B']==seeds and len(book['A'])==len(book['B'])==256 and len(set(seeds))==512,'Wrong old ordered banks')
    require(index['independent_upper_comparison_seeds']==book['A'][:128]+book['B'][:128],'Upper subset changed')
    require(metadata['seedbook']['seed_blocks']['holdout']==book,'Old bank identity changed')
    freeze=metadata['freeze'];require(freeze['identity_sha256']==index['freeze']['identity_sha256']==identity(freeze),'Old freeze identity changed')
    for name,h in freeze['files'].items():require(digest(relative(repo,name))==h,'Old frozen source/config changed')
    for name,h in index['sources'].items():require(digest(relative(repo,name))==h,'Old scientific source changed');ledger[name]=h
    require(set(index['record_metadata'])==set(index['records']) and len(index['records'])==1536,'Incomplete indexed records')
    manifests={}
    for name,h in index['manifests'].items():
        p=relative(repo,name);require(digest(p)==h,'Old manifest changed');m=read(p)
        require(m['identity']==identity(m,'identity'),'Old manifest identity changed')
        expected_sources={p:index['sources'][p] for p in SOURCE_PATHS}
        expected_sources['array_worker.py']=index['sources']['research/full_icu_response_2026-10-08/array_worker.py']
        exact(m['sources'],expected_sources,'old sources');manifests[name]=m;ledger[name]=h
        inventory(p.parent,len(m['jobs']))
    values=np.zeros((512,3,4));windows=np.zeros((512,3,8));seen=set()
    for name,h in index['records'].items():
        p=relative(repo,name);require(digest(p)==h,'Old record changed');row=read(p);meta=index['record_metadata'][name]
        mname=str(Path(name).parent.parent/'manifest.json');require(mname in manifests,'Missing own old manifest')
        m=manifests[mname];j=row['job_id'];require(type(j) is int and 0<=j<len(m['jobs']) and p.name==f'job_{j:08d}.json','Old job mismatch')
        require(row['identity']==m['identity']==meta['manifest_identity'],'Wrong own original manifest')
        ci=next((i for i,s in enumerate(scenarios) if s==row['scenario']),None)
        require(ci is not None and row['seed'] in seeds and row['scenario']['id']==meta['scenario'] and row['seed']==meta['seed'],'Wrong old profile/seed')
        si=seeds.index(row['seed']);require((si,ci) not in seen,'Repeated old paired seed/profile')
        require(meta['statistical_bank']==('A' if si<256 else 'B') and
                meta['execution_site'] in ('local','cluster') and ('holdout_'+meta['execution_site'])==p.parent.parent.name,'Old bank/site association differs')
        values[si,ci],windows[si,ci]=validate_record(row,m['jobs'][j],m['software']);seen.add((si,ci));ledger[name]=h
    require(len(seen)==1536 and np.isfinite(values).all() and (values[:,:,[1,3]]>0).all(),'Incomplete old vectors')
    return dict(config=cfg,values=values,windows=windows,book=book,index=index,index_sha256=ih,ledger=ledger,repo=repo,index_path=index_path)


def interval(values,level):
    missing=int(np.count_nonzero(~np.isfinite(values)))
    return dict(level=level,undefined_draws=missing,interval=None if missing else
                np.quantile(values,[(1-level)/2,1-(1-level)/2]).tolist())


def difference_or_none(a,b):
    return None if a is None or b is None else float(a-b)


def verify_selection(selection,secondary,grid,book,parent_identity,policy_cfg=None):
    """Recompute the fixed point choice; never use policy/validation outcomes."""
    from select_clock import choose,policy_config
    require(type(selection['schema_version']) is int and selection['schema_version']==1 and
            selection['study']==HERE.name and selection['status']=='complete_baseline_selection' and
            selection['parent_freeze_identity']==parent_identity and
            selection['policy_outcomes_used'] is False and
            selection['bootstrap_draws']==50000 and selection['bootstrap_seed']==20261011 and
            selection['source_sha256']==digest(HERE/'select_clock.py'),'Wrong fixed baseline-selection receipt')
    require(type(secondary['schema_version']) is int and secondary['schema_version']==1 and
            secondary['study']==HERE.name and secondary['parent_freeze_identity']==parent_identity and
            secondary['identity_sha256']==identity(secondary),'Wrong secondary freeze')
    exact(selection['seed_banks'],dict(grid=book['clock_grid'],policy=book['clock_policy']),'selection seed banks')
    require(len(book['clock_grid']['A'])==len(book['clock_grid']['B'])==64,'Fixed grid A/B64 required')
    exact(grid['config']['seeds'],book['clock_grid']['A']+book['clock_grid']['B'],'grid ordered seeds')
    for obj in (selection,secondary):
        exact(obj['grid_manifest_sha256'],grid['manifest_sha256'],'selection grid manifest')
        exact(obj['grid_input_hashes'],grid['input_sha256'],'selection complete grid records')
    exact(selection['grid_binding_sha256'],grid['binding_sha256'],'selection grid binding')
    require(set(selection['profiles'])=={'clustered','imperfect'},'Exactly two planned selection profiles required')
    selected=[];cfg=grid['config'];totals=grid['values'].sum(axis=0)
    rates=1000*totals[:,0]/totals[:,1]
    for profile in ('clustered','imperfect'):
        ii=[i for i,s in enumerate(cfg['scenarios']) if s['id'].startswith(profile+'_')]
        chosen=choose([(cfg['scenarios'][i]['p'],float(rates[i])) for i in ii])
        actual=selection['profiles'][profile]
        exact({k:actual[k] for k in chosen},chosen,'point selection '+profile)
        points=[dict(scenario_id=cfg['scenarios'][i]['id'],p=cfg['scenarios'][i]['p'],
            acquired_total=int(totals[i,0]),patient_days_total=float(totals[i,1]),rate=float(rates[i])) for i in ii]
        exact(actual['grid_points'],points,'selection pooled grid '+profile)
        for label,rows in (('A',range(64)),('B',range(64,128))):
            total=grid['values'][list(rows)].sum(axis=0);rr=1000*total[:,0]/total[:,1]
            exact(actual['banks'][label],choose([(cfg['scenarios'][i]['p'],float(rr[i])) for i in ii]),
                'selection bank '+profile+'/'+label)
        if chosen['selected_p'] is not None:selected.append(profile)
    exact(selection['selected_profiles'],selected,'selected profile order')
    files={f'research/{HERE.name}/selection.json':digest(HERE/'selection.json')}
    if selected:
        require(policy_cfg is not None,'Selected policy config required')
        exact(policy_cfg,policy_config(cfg,book,selection['profiles']),'derived selected policy config')
        files[f'research/{HERE.name}/configs/clock_policy.json']=digest(HERE/'configs/clock_policy.json')
    else:require(policy_cfg is None,'No policy config allowed after both selections fail')
    exact(secondary['files'],files,'secondary file inventory')


def estimates(sums,n):
    r0=np.divide(1000*sums[...,0],sums[...,1],out=np.full(sums.shape[:-1],np.nan),where=sums[...,1]>0)
    r1=np.divide(1000*sums[...,2],sums[...,3],out=np.full(sums.shape[:-1],np.nan),where=sums[...,3]>0)
    f=1-np.divide(sums[...,2]*sums[...,1],sums[...,3]*sums[...,0],out=np.full(sums.shape[:-1],np.nan),where=(sums[...,0]>0)&(sums[...,3]>0))
    return np.stack((r0,r1,f,r0-r1,(sums[...,0]-sums[...,2])/n),axis=-1)


def boot(values,windows,strata,n,rng,draws=50000):
    require(sorted(i for s in strata for i in s)==list(range(len(values))),'Bootstrap strata must partition all rows')
    samples=np.empty((draws,values.shape[1],5));win_samples=np.empty((draws,values.shape[1],2));weight_hash=hashlib.sha256()
    for first in range(0,draws,500):
        last=min(first+500,draws);weights=np.zeros((last-first,len(values)),dtype=np.int64)
        for bank in strata:weights[:,bank]=rng.multinomial(len(bank),np.full(len(bank),1/len(bank)),size=last-first)
        weight_hash.update(weights.astype('<i8',copy=False).tobytes(order='C'))
        samples[first:last]=estimates(np.einsum('bn,nck->bck',weights,values),n)
        w=np.einsum('bn,nck->bck',weights,windows)
        win_samples[first:last,:,0]=1000*(w[:,:,2]/w[:,:,3]-w[:,:,0]/w[:,:,1])
        win_samples[first:last,:,1]=1000*(w[:,:,6]/w[:,:,7]-w[:,:,4]/w[:,:,5])
    return samples,win_samples,weight_hash.hexdigest()


def summarize(name,scenario,values,windows,samples,win_samples,n,bank_rows,matching=False):
    total=values.sum(axis=0);point=estimates(total,n);keys=('control_rate','policy_rate','F','absolute_rate_reduction','mean_carriers_prevented')
    result=dict(id=name,scenario=scenario,n=n,control_count=int(total[0]),control_patient_days=float(total[1]),
                policy_count=int(total[2]),policy_patient_days=float(total[3]))
    for i,k in enumerate(keys):result[k]=None if not np.isfinite(point[i]) else float(point[i]);result[k+'95']=interval(samples[:,i],.95)
    result['banks']={}
    for label,rows in bank_rows.items():
        p=estimates(values[rows].sum(axis=0),len(rows));result['banks'][label]=dict(n=len(rows),**{k:None if not np.isfinite(p[i]) else float(p[i]) for i,k in enumerate(keys)})
    w=windows.sum(axis=0);result['observation_halves']={}
    for arm,j in (('baseline',0),('policy',1)):
        off=4*j;result['observation_halves'][arm]=dict(first28_rate=float(1000*w[off]/w[off+1]),last28_rate=float(1000*w[off+2]/w[off+3]),last_minus_first95=interval(win_samples[:,j],.95))
    if matching:
        ci=interval(samples[:,0],1-.05/6);bounds=ci['interval'];result['baseline_target_interval']=ci
        result['baseline_target_band']=[2.7,3.3];result['sampled_baseline_target_gate']=bool(bounds and bounds[0]>=2.7 and bounds[1]<=3.3)
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--new-root',type=Path,required=True);parser.add_argument('--old-root',type=Path,required=True)
    parser.add_argument('--old-index',type=Path,default=HERE/'old_input_index.json');parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    from freeze import verify,verify_config
    frozen=verify();design=read(HERE/'design.json');book=read(HERE/'seedbook.json')
    require(type(design['schema_version']) is int and design['schema_version']==1 and
            design['bootstrap_draws']==50000 and design['bootstrap_seed']==20261011 and len(design['primary_contrasts'])==6,'Wrong fixed analysis design')
    require(args.old_index.resolve()==(HERE/'old_input_index.json').resolve() and
            frozen['files'].get(f'research/{HERE.name}/old_input_index.json')==digest(args.old_index),'Old reuse index must be master-frozen')
    selection,selection_sha=bound_read(HERE/'selection.json');secondary,secondary_sha=bound_read(HERE/'clock_policy_freeze.json')
    require(secondary['identity_sha256']==identity(secondary) and secondary['parent_freeze_identity']==frozen['identity_sha256'],'Wrong secondary selection freeze')
    for name,h in secondary['files'].items():require(digest(relative(ABM,name))==h,'Secondary selected source changed')
    require(selection['parent_freeze_identity']==frozen['identity_sha256'],'Wrong selection parent')
    old=load_old(args.old_root,args.old_index);banks={}
    modules=['clock_grid','surfaces_off','independent_upper','clock_boundary']
    if selection['selected_profiles']:modules.append('clock_policy')
    else:require(not (args.new_root/'clock_policy'/'completed').exists(),'Unexpected policy output after failed selection')
    # Require every complete inventory before parsing any new outcomes.
    for module in modules:
        cfg=read(HERE/'configs'/f'{module}.json');inventory(args.new_root/module,len(cfg['scenarios'])*len(cfg['seeds']))
    for module in modules:
        cfg=read(HERE/'configs'/f'{module}.json');f=verify_config(HERE/'configs'/f'{module}.json')
        banks[module]=load_bank(args.new_root/module,cfg,f['identity_sha256'],module)
    grid=banks['clock_grid']
    verify_selection(selection,secondary,grid,book,frozen['identity_sha256'],
        None if 'clock_policy' not in banks else banks['clock_policy']['config'])
    oldn=512;data=np.zeros((512,8,4));windows=np.zeros((512,8,8));data[:,:3]=old['values'];windows[:,:3]=old['windows']
    data[:,3:6]=banks['surfaces_off']['values'];windows[:,3:6]=banks['surfaces_off']['windows']
    exact(banks['surfaces_off']['config']['seeds'],old['config']['seeds'],'Reused surface seed order')
    core=list(range(128))+list(range(256,384));expected=[old['config']['seeds'][i] for i in core]
    exact(banks['independent_upper']['config']['seeds'],expected,'Upper seed subset')
    data[core,6]=banks['independent_upper']['values'][:,0];windows[core,6]=banks['independent_upper']['windows'][:,0]
    data[core,7]=old['values'][core,0];windows[core,7]=old['windows'][core,0]
    ns=np.array([512]*6+[256]*2)
    rng=np.random.default_rng(20261011);samples,win_samples,wh=boot(data,windows,[list(range(a,a+128)) for a in (0,128,256,384)],ns,rng)
    entries=[];full_scenarios=old['config']['scenarios'];off_scenarios=banks['surfaces_off']['config']['scenarios']
    names=['old_full_'+s['id'] for s in full_scenarios]+[s['id'] for s in off_scenarios]+['independent_upper','independent_center_matched256']
    scenarios=full_scenarios+off_scenarios+[banks['independent_upper']['config']['scenarios'][0],full_scenarios[0]]
    for i,name in enumerate(names):
        rows={'A':list(range(256)),'B':list(range(256,512))} if i<6 else {'A':list(range(128)),'B':list(range(256,384))}
        entries.append(summarize(name,scenarios[i],data[:,i],windows[:,i],samples[:,i],win_samples[:,i],int(ns[i]),rows))
    contrasts=[];family=1-.05/6
    for i,profile in enumerate(('independent','persistent','imperfect')):
        contrasts.append(dict(name=f'surfaces_off_minus_full_{profile}_F',estimate=difference_or_none(entries[i+3]['F'],entries[i]['F']),uncertainty=interval(samples[:,i+3,2]-samples[:,i,2],family)))
    contrasts.append(dict(name='upper_minus_center_independent_F',estimate=difference_or_none(entries[6]['F'],entries[7]['F']),uncertainty=interval(samples[:,6,2]-samples[:,7,2],family)))
    weights=dict(old_four_strata=wh)
    for module in ('clock_policy','clock_boundary'):
        if module not in banks:continue
        bank=banks[module];values=bank['values'];win=bank['windows'];n=len(values)
        b=book[module];require(b['A']+b['B']==bank['config']['seeds'] and len(b['A'])==len(b['B'])==256,'Wrong fresh clock banks')
        samples_c,windows_c,weights[module]=boot(values,win,[list(range(256)),list(range(256,512))],n,rng)
        for i,sc in enumerate(bank['config']['scenarios']):
            entry=summarize(sc['id'],sc,values[:,i],win[:,i],samples_c[:,i],windows_c[:,i],n,{'A':list(range(256)),'B':list(range(256,512))},matching=module=='clock_policy')
            if module=='clock_boundary':entry['interpretation']='Unfitted p=1 boundary comparison; not a global feasibility proof.'
            entries.append(entry)
            if module=='clock_policy':
                profile='clustered' if sc['id'].startswith('clustered_') else 'imperfect';old_i=1 if profile=='clustered' else 2
                contrasts.append(dict(name=f'clock_minus_zone_{profile}_F',estimate=difference_or_none(entry['F'],entries[old_i]['F']),uncertainty=interval(samples_c[:,i,2]-samples[:,old_i,2],family)))
    for profile in ('clustered','imperfect'):
        name=f'clock_minus_zone_{profile}_F'
        if not any(c['name']==name for c in contrasts):contrasts.append(dict(name=name,estimate=None,uncertainty=dict(level=family,undefined_draws=None,interval=None),unavailable_reason='Prespecified baseline-only selection failed; family size remains six.'))
    require({c['name'] for c in contrasts}==set(design['primary_contrasts']),'Wrong primary family')
    contrasts.sort(key=lambda c:design['primary_contrasts'].index(c['name']))
    report=dict(schema_version=1,study=HERE.name,status='complete',bootstrap_draws=50000,bootstrap_seed=20261011,
        primary_family_size=6,primary_interval_level=family,descriptive_interval_level=.95,
        bootstrap_method='Synchronous whole-seed vectors: old A-core/A-remainder/B-core/B-remainder128 each; fresh clock A/B256 independently. Fixed chunk500; groups old,selected-clock,boundary.',
        weight_stream_sha256_int64LE=weights,scenarios=entries,primary_contrasts=contrasts,
        selected_profiles=selection['selected_profiles'],selection_sha256=selection_sha,
        secondary_freeze_sha256=secondary_sha,
        parent_freeze_identity=frozen['identity_sha256'],secondary_freeze_identity=secondary['identity_sha256'],
        old_input_index_sha256=old['index_sha256'],source_sha256=digest(Path(__file__)),
        new_bank_bindings={m:dict(manifest_sha256=b['manifest_sha256'],binding_sha256=b['binding_sha256'],input_sha256=b['input_sha256']) for m,b in banks.items()},
        interpretation='Separate fixed follow-up. Policy intervals condition on each once-selected p; baseline selection uncertainty is recorded separately and is not used to reselect or propagated as hypothetical policy outcomes. Old point forecasts are unchanged; their newly recomputed four-stratum intervals do not replace the archived two-bank original intervals. Old matching failures unchanged; clock-versus-zone fixed configurations; surface ablation includes feedback; null intervals do not establish equivalence; boundary is not a fit or global maximum.')
    verify()
    for b in banks.values():recheck_bank(b)
    require(digest(args.old_index)==old['index_sha256'],'Old index changed')
    for name,h in old['ledger'].items():require(digest(relative(args.old_root,name))==h,'Old input changed during analysis')
    require(digest(HERE/'selection.json')==selection_sha and digest(HERE/'clock_policy_freeze.json')==secondary_sha,
            'Selection bytes changed during analysis')
    safe(args.output);require(not args.output.exists(),'Fresh analysis output required');args.output.mkdir(parents=True)
    text=json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n'
    (args.output/'report.json').write_text(text)
    print('Complete frozen rule/surface study analyzed')


if __name__=='__main__':main()

