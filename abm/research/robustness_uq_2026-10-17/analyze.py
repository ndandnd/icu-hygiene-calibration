#!/usr/bin/env python3
"""Prospective bounded workflow-study selection and fixed-p policy inference.

CLI: --stage select|analyze --bank DIR --config FILE --output FILE --probe FILE
Selection additionally requires --policy-seeds seedbook.json and may write
--policy-config FILE. The seedbook contains a policy seed list. Analysis requires
--selection FILE and --derived-freeze FILE. Outputs are immutable; CLI is pinned.
Public interfaces: select(records,cfg,policy_seeds), analyze_bank(records,cfg).
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import sys
import uuid
import numpy as np

HERE=Path(__file__).resolve().parent
ABM=HERE.parents[1]
PROFILE_ORDER=('independent','clustered','imperfect')
PROFILE_PARAMETERS={'independent':(0.,1.),'clustered':(.5,1.),'imperfect':(0.,.8)}
TARGET=3.
BAND=(2.7,3.3)
DRAWS=20000
BOOTSTRAP_SEED=20261017
FAMILY_LEVEL=1-.05/3
OLD_ANALYZER=ABM/'research/rule_surface_followup_2026-10-11/analyze.py'
spec=importlib.util.spec_from_file_location('_workflow_accounting_validator',OLD_ANALYZER)
accounting=importlib.util.module_from_spec(spec)
spec.loader.exec_module(accounting)
accounting.RESPONSE_KINDS=accounting.RESPONSE_KINDS+('workflow_events',)
require=accounting.require
exact=accounting.exact
read=accounting.read
digest=accounting.digest
safe=accounting.safe
canonical=lambda x:json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)
identity=lambda x,key='identity':hashlib.sha256(canonical({k:v for k,v in x.items() if k!=key}).encode()).hexdigest()
_original_reconcile=accounting.reconcile


def reconcile(part,snapshot=False):
    _original_reconcile(part,snapshot)
    response=part['response_counts'];env=response['environment_events'];flow=response['workflow_events']
    require(set(env)<={'scheduled_clear','dirty_scheduled_clear','deposition','discharge_clear',
        'dirty_discharge_clear','relocation_clear','dirty_relocation_clear'},'Unknown environment event')
    require(set(flow)<= {'isolation_relocation','isolation_deferred_no_empty_bed','occupied_destination_rejected'},
            'Unknown workflow counter')
    require(flow.get('occupied_destination_rejected',0)==0,'Occupied destination attempted in production')
    for kind in ('discharge','relocation'):
        require(env.get('dirty_'+kind+'_clear',0)<=env.get(kind+'_clear',0),'Dirty terminal clear exceeds clears')
    require(env.get('discharge_clear',0)==part['discharged_patients'],'Discharge terminal-clear reconciliation')
    require(env.get('relocation_clear',0)==flow.get('isolation_relocation',0),'Relocation terminal-clear reconciliation')
    net=env.get('deposition',0)-sum(env.get(k,0) for k in
        ('dirty_scheduled_clear','dirty_discharge_clear','dirty_relocation_clear'))
    endpoint=part['environment_contaminated_sites_endpoint'] if snapshot else (
        part['environment_contaminated_sites_endpoint_end']-part['environment_contaminated_sites_endpoint_start'])
    require(net==endpoint,'Clean-initial environmental mass balance failed')


accounting.reconcile=reconcile


def profile(sc):
    name=sc['id'].split('_',1)[0]
    require(name in PROFILE_ORDER,'Unknown new-study profile')
    return name


def jobs_for(cfg):
    return [dict(job_id=i,scenario=s,seed=seed) for i,(s,seed) in
            enumerate((s,seed) for s in cfg['scenarios'] for seed in cfg['seeds'])]


def validate_config(cfg):
    require(isinstance(cfg,dict) and set(cfg)<= {'scenarios','seeds','warmup_days','observation_days',
        'block_days','description','model_overrides'},'Unknown study config field')
    require(cfg['warmup_days']==90 and cfg['observation_days']==56 and cfg['block_days']==1,'Wrong fixed horizons')
    exact(cfg['model_overrides'],dict(prob_patient_sick=.03,initial_environment_contaminated=False),'overrides')
    seeds=cfg['seeds'];ss=cfg['scenarios']
    require(isinstance(seeds,list) and seeds and all(type(x) is int and 0<=x<2**32 for x in seeds)
        and len(set(seeds))==len(seeds),'Unique 32-bit seeds required')
    require(isinstance(ss,list) and ss and len({s['id'] for s in ss})==len(ss),'Unique scenarios required')
    for sc in ss:
        require(set(sc)=={'id','p','h0','h1','rho','e','opportunities','environment_multiplier'},'Wrong scenario schema')
        require(isinstance(sc['id'],str) and sc['id'],'Scenario id required')
        require(type(sc['p']) in (int,float) and math.isfinite(sc['p']) and 0<sc['p']<=1,'Positive finite p required')
        rho,e=PROFILE_PARAMETERS[profile(sc)]
        require(sc['rho']==rho and sc['e']==e and sc['h0']==.9 and sc['h1'] in (.9,.95)
            and sc['opportunities']=='zone_separate' and sc['environment_multiplier']==1.,'Profile/settings drift')
    require(len({s['h1'] for s in ss})==1,'Do not mix grid and policy records')


def validate_records(records,cfg,software=None,manifest_identity=None):
    """Validate complete Cartesian bank and independently checked daily accounting."""
    validate_config(cfg);jobs=jobs_for(cfg);n=len(cfg['seeds'])
    require(isinstance(records,list) and len(records)==len(jobs),'Complete config-by-seed bank required')
    if software is None:
        software={k:v for k,v in records[0]['execution'].items() if k!='platform'}
    validate_software(software)
    mid=manifest_identity if manifest_identity is not None else records[0]['identity']
    require(isinstance(mid,str) and mid,'Record identity required')
    values=np.empty((n,len(cfg['scenarios']),4));seen=set()
    for row in records:
        j=row['job_id']
        require(type(j) is int and 0<=j<len(jobs) and j not in seen and row['identity']==mid,'Duplicate/wrong record identity')
        require(row['execution']['platform'].startswith('Linux'),'Production bank must use Linux')
        point,_=accounting.validate_record(row,jobs[j],software)
        ci,si=divmod(j,n);values[si,ci]=point;seen.add(j)
    require(seen==set(range(len(jobs))) and np.isfinite(values).all(),'Complete finite bank required')
    return values


def validate_software(software):
    require(set(software)==set(accounting.PACKAGES)|{'python'} and
        all(software[k]==v for k,v in accounting.PACKAGES.items()) and
        isinstance(software['python'],str) and software['python'].split()[0]=='3.12.3','Pinned production runtime required')

def production_anchor(path,parent_identity):
    probe,sha=accounting.bound_read(path)
    require(probe.get('status')=='passed' and probe.get('production') is True and
        probe.get('reference_verified') is True and probe.get('parent_freeze_identity')==parent_identity,
        'Successful parent-bound production probe required')
    versions=probe.get('versions',{})
    require(versions.get('platform','').startswith('Linux'),'Production probe must use Linux')
    software={k:v for k,v in versions.items() if k!='platform'}
    validate_software(software)
    return software,sha


def estimates(total,n):
    r0=np.divide(1000*total[...,0],total[...,1],out=np.full(total.shape[:-1],np.nan),where=total[...,1]>0)
    r1=np.divide(1000*total[...,2],total[...,3],out=np.full(total.shape[:-1],np.nan),where=total[...,3]>0)
    f=1-np.divide(r1,r0,out=np.full_like(r0,np.nan),where=r0>0)
    return np.stack((r0,r1,f,r0-r1,(total[...,0]-total[...,2])/n),axis=-1)


def interval(x,level=.95):
    x=np.asarray(x);valid=x[np.isfinite(x)];bad=int(x.size-valid.size)
    quantiles=None if not len(valid) else np.quantile(valid,[(1-level)/2,1-(1-level)/2]).tolist()
    return dict(level=level,undefined_draws=bad,defined_draws=len(valid),interval=None if bad else quantiles,
                conditional_on_defined_interval=quantiles if bad else None)


def finite(x):
    return float(x) if np.isfinite(x) else None


def bootstrap_totals(values,stream,draws=DRAWS):
    """Resample complete seed vectors synchronously across all conditions."""
    n=len(values);require(n>0,'Nonempty seed bank required')
    rng=np.random.default_rng(np.random.SeedSequence([BOOTSTRAP_SEED,stream]))
    sha=hashlib.sha256();samples=np.empty((draws,values.shape[1],4))
    for first in range(0,draws,500):
        size=min(500,draws-first)
        weights=rng.multinomial(n,np.full(n,1/n),size=size)
        sha.update(weights.astype('<i8',copy=False).tobytes(order='C'))
        samples[first:first+size]=np.einsum('bn,nck->bck',weights,values)
    return samples,sha.hexdigest()


def choose(points,target=TARGET):
    require(points and all(type(p) in (int,float) and type(r) in (int,float) and
        math.isfinite(p) and math.isfinite(r) and 0<p<=1 and r>=0 for p,r in points),'Finite positive-p nonnegative-rate grid required')
    require(all(a[0]<b[0] for a,b in zip(points,points[1:])),'Strictly increasing p grid required')
    hits=[dict(p=float(p),rate=float(r)) for p,r in points if r==target]
    brackets=[];nonincreasing=[];zero_brackets=[]
    for low,high in zip(points,points[1:]):
        if high[1]<=low[1]:nonincreasing.append([list(low),list(high)])
        if low[1]<target<high[1]:
            if low[1]==0:
                zero_brackets.append([list(low),list(high)]);continue
            w=(math.log(target)-math.log(low[1]))/(math.log(high[1])-math.log(low[1]))
            p=math.exp((1-w)*math.log(low[0])+w*math.log(high[0]))
            brackets.append(dict(p_low=low[0],p_high=high[0],rate_low=low[1],rate_high=high[1],
                                 log_weight=w,selected_p=p))
    details=dict(exact_hits=hits,increasing_positive_adjacent_brackets=brackets,
        nonincreasing_adjacent_segments=nonincreasing,zero_lower_rate_brackets=zero_brackets)
    if hits:return dict(details,status='selected',method='lowest_exact_target_hit',selected_p=hits[0]['p'])
    if brackets:return dict(details,status='selected',method='lowest_positive_increasing_adjacent_log_log_bracket',selected_p=brackets[0]['selected_p'])
    return dict(details,status='failed',method=None,selected_p=None,reason='No exact hit or positive increasing adjacent bracket; no fallback.')


def select_values(values,cfg,policy_seeds):
    require(len(cfg['seeds'])==128 and len(cfg['scenarios'])==12,'Fixed 128-seed twelve-condition grid required')
    require(all(s['h1']==.9 for s in cfg['scenarios']),'Only baseline outcomes may select')
    require(isinstance(policy_seeds,list) and len(policy_seeds)==512 and len(set(policy_seeds))==512 and
        all(type(s) is int and 0<=s<2**32 for s in policy_seeds) and not set(policy_seeds)&set(cfg['seeds']),
        '512 distinct fresh policy seeds required')
    book=read(HERE/'seedbook.json')
    exact(cfg['seeds'],book['grid'],'frozen grid seeds')
    exact(policy_seeds,book['policy'],'frozen policy seeds')
    totals=values.sum(axis=0);rates=1000*totals[:,0]/totals[:,1]
    boot,weight_hash=bootstrap_totals(values,0);draw_rates=1000*boot[:,:,0]/boot[:,:,1]
    profiles={};policy=[]
    old=read(ABM/'research/full_icu_confirmation_2026-10-08/configs/holdout_all.json')
    multipliers=(.5,1.,2.,4.);suffixes=('m0p5','m1','m2','m4')
    for index,name in enumerate(PROFILE_ORDER):
        ii=[i for i,s in enumerate(cfg['scenarios']) if profile(s)==name]
        require(len(ii)==4,'Four fixed points per profile required')
        ss=[cfg['scenarios'][i] for i in ii]
        require([s['id'] for s in ss]==[name+'_'+suffix for suffix in suffixes] and
            [s['p'] for s in ss]==[old['scenarios'][index]['p']*x for x in multipliers],'Grid identity drift')
        selected=choose([(s['p'],float(rates[i])) for s,i in zip(ss,ii)])
        selected['grid_points']=[dict(scenario_id=s['id'],p=s['p'],control_count=int(totals[i,0]),
            control_patient_days=float(totals[i,1]),rate=float(rates[i])) for s,i in zip(ss,ii)]
        draws=[];brackets={};failed=0;nonincreasing=0;zero=0
        for row in draw_rates:
            result=choose([(s['p'],float(row[i])) for s,i in zip(ss,ii)])
            nonincreasing+=bool(result['nonincreasing_adjacent_segments']);zero+=bool(result['zero_lower_rate_brackets'])
            if result['selected_p'] is None:failed+=1;continue
            draws.append(result['selected_p'])
            key='exact:'+str(result['selected_p']) if result['method']=='lowest_exact_target_hit' else str(
                result['increasing_positive_adjacent_brackets'][0]['p_low'])+':'+str(result['increasing_positive_adjacent_brackets'][0]['p_high'])
            brackets[key]=brackets.get(key,0)+1
        selected['selection_uncertainty']=dict(draws=DRAWS,valid_draws=len(draws),failed_draws=failed,
            failure_fraction=failed/DRAWS,conditional_on_success95=None if not draws else np.quantile(draws,[.025,.975]).tolist(),
            conditional_interval=True,selected_bracket_frequencies=brackets,
            draws_with_nonincreasing_segments=nonincreasing,draws_with_zero_lower_brackets=zero)
        profiles[name]=selected
        if selected['selected_p'] is not None:
            sc=dict(ss[0],id=name,p=selected['selected_p'],h1=.95);policy.append(sc)
    policy_cfg=None if not policy else dict({k:cfg[k] for k in ('warmup_days','observation_days','block_days','model_overrides')},
        scenarios=policy,seeds=policy_seeds,description='Once-selected alternative workflow; fixed p and fresh paired policy holdout.')
    return dict(schema_version=1,study=HERE.name,status='complete_baseline_selection',profiles=profiles,
        policy_config=policy_cfg,policy_outcomes_used=False,bootstrap_draws=DRAWS,bootstrap_seed=BOOTSTRAP_SEED,
        bootstrap_stream=0,weight_stream_sha256_int64LE=weight_hash,
        interpretation='p uncertainty conditions on positive valid selections; failures retained. No policy response is inferred at bootstrap-selected p.')


def select(records,cfg,policy_seeds):
    return select_values(validate_records(records,cfg),cfg,policy_seeds)


def summarize_values(values,cfg):
    n=len(values);require(n==512 and len(cfg['scenarios'])<=3,'Fixed 512-seed selected policy bank required')
    names=[s['id'] for s in cfg['scenarios']]
    require(names==[p for p in PROFILE_ORDER if p in names] and all(s['h1']==.95 for s in cfg['scenarios']),
        'Selected profile order/settings drift')
    point=estimates(values.sum(axis=0),n);totals,weight_hash=bootstrap_totals(values,1)
    samples=estimates(totals,n);keys=('control_rate','policy_rate','F','absolute_rate_reduction','mean_carriers_prevented')
    entries=[]
    for i,s in enumerate(cfg['scenarios']):
        summary=dict(scenario=s,n=n)
        for j,key in enumerate(keys):summary[key]=finite(point[i,j]);summary[key+'95']=interval(samples[:,i,j])
        ci=summary['control_rate95']['interval']
        summary['descriptive_baseline_target_band']=list(BAND)
        summary['descriptive_baseline95_inside_target_band']=bool(ci is not None and ci[0]>=BAND[0] and ci[1]<=BAND[1])
        entries.append(summary)
    contrasts=[]
    for left,right in (('independent','clustered'),('independent','imperfect'),('clustered','imperfect')):
        item=dict(name=left+'_minus_'+right+'_F')
        if left in names and right in names:
            i,j=names.index(left),names.index(right)
            item.update(estimate=finite(point[i,2]-point[j,2]),uncertainty=interval(samples[:,i,2]-samples[:,j,2],FAMILY_LEVEL))
        else:item.update(estimate=None,uncertainty=dict(level=FAMILY_LEVEL,undefined_draws=None,interval=None),
                         unavailable_reason='Prespecified baseline selection failed; family remains three.')
        contrasts.append(item)
    return dict(schema_version=1,study=HERE.name,status='complete_fixed_p_policy_analysis',scenarios=entries,
        primary_contrasts=contrasts,bootstrap_draws=DRAWS,bootstrap_seed=BOOTSTRAP_SEED,bootstrap_stream=1,
        weight_stream_sha256_int64LE=weight_hash,bootstrap_method='Synchronous whole-seed vectors over all512 seeds; policy stream independent of grid stream.',
        interpretation='Policy intervals condition on once-selected fixed p; selection uncertainty is not propagated to unobserved policy outcomes. Target-band checks are descriptive, not profile baseline equivalence. Original failures unchanged.')


def analyze_bank(records,cfg):
    exact(cfg['seeds'],read(HERE/'seedbook.json')['policy'],'frozen policy seeds')
    return summarize_values(validate_records(records,cfg),cfg)


def verify_freeze(path):
    frozen=read(path)
    require(frozen['study']==HERE.name and frozen['identity_sha256']==identity(frozen,'identity_sha256'),'Invalid study freeze')
    for name,sha in frozen['files'].items():
        p=accounting.relative(ABM,name);require(digest(p)==sha,'Frozen file changed: '+name)
    return frozen


def load_bank(bank,cfg,frozen_sources=None,freeze_identity=None,parent_freeze_identity=None,
              probe_software=None,probe_sha=None):
    """Stream full records into a compact array; never retain all daily histories."""
    validate_config(cfg)
    bank=safe(bank);manifest,manifest_sha=accounting.bound_read(bank/'manifest.json');mid=manifest['identity']
    require(type(manifest['schema_version']) is int and manifest['schema_version']==1 and mid==identity(manifest),'Manifest identity mismatch')
    require(manifest['study']==HERE.name and manifest['model_factory']=='WorkflowRobustnessModel' and
        type(manifest['task_count']) is int and manifest['task_count']>0,'Wrong study/model/partition identity')
    if frozen_sources is not None:exact(manifest['sources'],frozen_sources,'frozen source closure')
    if freeze_identity is not None:require(manifest['freeze_identity']==freeze_identity,'Wrong execution freeze')
    if parent_freeze_identity is not None:require(manifest['parent_freeze_identity']==parent_freeze_identity,'Wrong execution parent freeze')
    if probe_software is not None:
        require(probe_sha is not None and manifest.get('production_probe_sha256')==probe_sha,
                'Bank is not bound to production probe')
        exact(manifest['software'],probe_software,'production probe runtime including full Python build')
    exact(manifest['config'],cfg,'manifest config')
    exact(manifest['resolved_config'],dict(cfg,warmup_ticks=233280,observation_ticks=145152,block_ticks=2592),'resolved config')
    jobs=jobs_for(cfg);exact(manifest['jobs'],jobs,'manifest jobs')
    validate_software(manifest['software'])
    required=set(accounting.SOURCE_PATHS)|{f'research/{HERE.name}/workflow_variant.py',f'research/{HERE.name}/worker.py'}
    require(required<=set(manifest['sources']),'Scientific manifest omits inherited/variant source')
    for name,sha in manifest['sources'].items():
        p=accounting.relative(ABM,name)
        require(digest(p)==sha,'Manifest source changed: '+name)
    paths=accounting.inventory(bank,len(jobs));hashes={};seen=set();n=len(cfg['seeds'])
    values=np.empty((n,len(cfg['scenarios']),4))
    for path in paths:
        row,sha=accounting.bound_read(path)
        j=row['job_id']
        require(type(j) is int and 0<=j<len(jobs) and j not in seen and path.name==f"job_{j:08d}.json" and
            row['identity']==mid,'Filename/job/identity mismatch')
        require(row['execution']['platform'].startswith('Linux'),'Production bank must use Linux')
        point,_=accounting.validate_record(row,jobs[j],manifest['software'])
        ci,si=divmod(j,n);values[si,ci]=point;seen.add(j);hashes['completed/'+path.name]=sha
    require(seen==set(range(len(jobs))) and np.isfinite(values).all(),'Complete finite bank required')
    receipt=dict(manifest_identity=mid,manifest_sha256=manifest_sha,complete_input_hashes=hashes)
    return values,receipt


def write_new(path,value):
    path=safe(path);payload=json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if path.exists():require(path.read_text()==payload,'Existing output differs; refusing overwrite');return
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
    try:
        with temp.open('x') as stream:stream.write(payload);stream.flush();os.fsync(stream.fileno())
        try:os.link(temp,path)
        except FileExistsError:require(path.read_text()==payload,'Concurrent output differs; refusing overwrite')
    finally:
        if temp.exists():temp.unlink()


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('select','analyze'),required=True)
    parser.add_argument('--bank',type=Path,required=True);parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--policy-seeds',type=Path)
    parser.add_argument('--policy-config',type=Path);parser.add_argument('--freeze',type=Path,default=HERE/'freeze.json')
    parser.add_argument('--derived-freeze',type=Path)
    parser.add_argument('--selection',type=Path)
    parser.add_argument('--probe',type=Path,required=True)
    args=parser.parse_args(argv)
    require(sys.version.split()[0]=='3.12.3' and np.__version__=='1.26.4','Pinned analysis Python3.12.3/NumPy1.26.4 required')
    frozen=verify_freeze(args.freeze);execution_freeze=frozen;cfg=read(args.config);name=str(safe(args.config).relative_to(ABM))
    if frozen['files'].get(name)!=digest(args.config):
        require(args.derived_freeze is not None,'Policy config needs derived freeze')
        derived=verify_freeze(args.derived_freeze)
        require(derived['parent_freeze_identity']==frozen['identity_sha256'] and derived['files'].get(name)==digest(args.config),'Unbound derived config')
        execution_freeze=derived
    if args.stage=='analyze':
        require(args.selection is not None and execution_freeze is not frozen,'Selected-policy analysis requires derived selection binding')
        selection=read(args.selection);selection_name=str(safe(args.selection).relative_to(ABM))
        require(execution_freeze['files'].get(selection_name)==digest(args.selection) and
            selection['parent_freeze_identity']==frozen['identity_sha256'],'Selection receipt is not bound to derived freeze')
        exact(cfg,selection['policy_config'],'selected policy config')
    sources=dict(frozen['files']);sources.update(execution_freeze['files'])
    probe_software,probe_sha=production_anchor(args.probe,frozen['identity_sha256'])
    require(sys.version==probe_software['python'],'Analysis Python build differs from production probe')
    values,receipt=load_bank(args.bank,cfg,sources,execution_freeze['identity_sha256'],frozen['identity_sha256'],
                             probe_software,probe_sha)
    if args.stage=='select':
        require(args.policy_seeds is not None,'Policy seedbook required')
        book=read(args.policy_seeds);seeds=book['policy']
        require(safe(args.policy_seeds)==safe(HERE/'seedbook.json'),'Use this study frozen seedbook')
        require(frozen['files'].get(str(safe(args.policy_seeds).relative_to(ABM)))==digest(args.policy_seeds),'Policy seedbook not frozen')
        report=select_values(values,cfg,seeds)
    else:
        exact(cfg['seeds'],read(HERE/'seedbook.json')['policy'],'frozen policy seeds')
        report=summarize_values(values,cfg)
    report.update(input_receipt=receipt,config_sha256=digest(args.config),parent_freeze_identity=frozen['identity_sha256'],
        production_probe_sha256=probe_sha,production_probe_software=probe_software,
        analysis_runtime=dict(python=sys.version,numpy=np.__version__,platform=platform.platform()),source_sha256=digest(Path(__file__)))
    verify_freeze(args.freeze)
    if execution_freeze is not frozen:verify_freeze(args.derived_freeze)
    require(digest(args.bank/'manifest.json')==receipt['manifest_sha256'],'Manifest changed during analysis')
    require(digest(args.probe)==probe_sha,'Production probe changed during analysis')
    accounting.inventory(args.bank,len(receipt['complete_input_hashes']))
    for name,sha in receipt['complete_input_hashes'].items():require(digest(args.bank/name)==sha,'Input changed during analysis')
    write_new(args.output,report)
    if args.stage=='select' and args.policy_config is not None and report['policy_config'] is not None:write_new(args.policy_config,report['policy_config'])
    print(canonical(dict(status=report['status'],output=str(args.output))))


if __name__=='__main__':main()
