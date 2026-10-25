#!/usr/bin/env python3
"""Declared pooled-rate/paired-contrast analysis of complete new banks only."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from freeze import verify

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--module', choices=['clock','sensitivity'], required=True)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    frozen=verify()
    binding=json.loads((args.results/'study_binding.json').read_text())
    if binding['freeze_identity']!=frozen['identity_sha256'] or binding['study']!=HERE.name or binding['config_sha256']!=digest(HERE/'configs'/(args.module+'.json')):
        raise ValueError('Wrong frozen study binding')
    cfg = json.loads((HERE / 'configs' / (args.module+'.json')).read_text())
    seedbook = json.loads((HERE / 'seedbook.json').read_text())[args.module]
    manifest = json.loads((args.results / 'manifest.json').read_text())
    scenarios, seeds = cfg['scenarios'], cfg['seeds']
    if manifest['config'] != cfg:
        raise ValueError('Wrong configuration')
    files = sorted((args.results / 'completed').glob('job_*.json'))
    if len(files) != len(scenarios)*len(seeds):
        raise ValueError('Incomplete bank; analysis refused')
    data = np.zeros((len(seeds),len(scenarios),4))
    windows = np.zeros((len(seeds),len(scenarios),4))
    seen = set()
    file_hashes = {}
    for path in files:
        row = json.loads(path.read_text())
        job = row['job_id']; ci, si = divmod(job,len(seeds))
        if job in seen or not 0 <= ci < len(scenarios) or row['scenario'] != scenarios[ci] or row['seed'] != seeds[si]:
            raise ValueError('Wrong or duplicate record')
        if row['status'] != 'complete' or row['identity'] != manifest['identity']:
            raise ValueError('Wrong record identity')
        seen.add(job);file_hashes[path.name]=digest(path)
        b = row['arms']['baseline']['observation']
        q = b if args.module == 'clock' else row['arms']['policy']['observation']
        data[si,ci] = [b['acquired_colonizations'], b['patient_days'], q['acquired_colonizations'], q['patient_days']]
        blocks=row['arms']['baseline']['observation_blocks']
        if len(blocks)!=56:raise ValueError('Unexpected daily observation blocks')
        windows[si,ci]=[sum(x['acquired_colonizations'] for x in blocks[:28]),sum(x['patient_days'] for x in blocks[:28]),
                       sum(x['acquired_colonizations'] for x in blocks[28:]),sum(x['patient_days'] for x in blocks[28:])]
    if len(seen) != len(files) or (data[:,:,[1,3]]<=0).any():
        raise ValueError('Incomplete/invalid exposure')
    total=data.sum(axis=0);r0=1000*total[:,0]/total[:,1];r1=1000*total[:,2]/total[:,3]
    fraction=np.divide(total[:,2]*total[:,1],total[:,3]*total[:,0],out=np.full(len(scenarios),np.nan),where=total[:,0]>0)
    fraction=1-fraction
    rng=np.random.default_rng(20261010);draws=50000
    samples=np.empty((draws,len(scenarios),3))
    window_samples=np.empty((draws,len(scenarios)))
    banks=[[seeds.index(s) for s in seedbook[name]] for name in ('A','B')]
    for start in range(0,draws,500):
        stop=min(start+500,draws);weights=np.zeros((stop-start,len(seeds)))
        for bank in banks:
            weights[:,bank]=rng.multinomial(len(bank),np.full(len(bank),1/len(bank)),size=stop-start)
        sums=np.einsum('bn,nck->bck',weights,data)
        samples[start:stop,:,0]=1000*sums[:,:,0]/sums[:,:,1]
        samples[start:stop,:,1]=1000*sums[:,:,2]/sums[:,:,3]
        ratio=np.divide(sums[:,:,2]*sums[:,:,1],sums[:,:,3]*sums[:,:,0],
                        out=np.full((stop-start,len(scenarios)),np.nan),where=sums[:,:,0]>0)
        samples[start:stop,:,2]=1-ratio
        win=np.einsum('bn,nck->bck',weights,windows)
        window_samples[start:stop]=1000*(win[:,:,2]/win[:,:,3]-win[:,:,0]/win[:,:,1])
    def interval(values,level):
        count=int((~np.isfinite(values)).sum())
        return dict(level=level,undefined_draws=count,interval=None if count else
                    np.quantile(values,[(1-level)/2,1-(1-level)/2]).tolist())
    family=1-.05/12 if args.module=='clock' else 1-.05/4
    entries=[]
    for i,sc in enumerate(scenarios):
        entry=dict(scenario=sc,control_rate=float(r0[i]),policy_rate=float(r1[i]),
                   F=None if not np.isfinite(fraction[i]) else float(fraction[i]),
                   mean_carriers_prevented=float((total[i,0]-total[i,2])/len(seeds)),
                   control_rate95=interval(samples[:,i,0],.95))
        sums=windows[:,i,:].sum(axis=0)
        entry['control_windows']=dict(first28_rate=float(1000*sums[0]/sums[1]),last28_rate=float(1000*sums[2]/sums[3]),
                                       last_minus_first95=interval(window_samples[:,i],.95))
        entry['banks']={}
        for name,bank in zip(('A','B'),banks):
            arm=data[bank,i,:].sum(axis=0)
            entry['banks'][name]=dict(n=len(bank),control_rate=float(1000*arm[0]/arm[1]),
                policy_rate=float(1000*arm[2]/arm[3]),F=None if arm[0]==0 else float(1-arm[2]*arm[1]/(arm[3]*arm[0])))
        if args.module=='clock':
            ci=interval(samples[:,i,0],family);bounds=ci['interval']
            entry.update(primary_rate_interval=ci,
                         sampled_point_approximately_matches=bool(bounds and bounds[0]>=2.7 and bounds[1]<=3.3))
            entry['target_classification']='unresolved' if bounds is None else ('below_target' if bounds[1]<3 else
                                         'above_target' if bounds[0]>3 else 'contains_target')
        else:
            entry.update(policy_rate95=interval(samples[:,i,1],.95),F95=interval(samples[:,i,2],.95))
        entries.append(entry)
    contrasts=[]
    brackets=[]
    if args.module=='clock':
        for profile in ('independent','clustered','imperfect'):
            selected=[e for e in entries if e['scenario']['id'].startswith(profile+'_')]
            for low,high in zip(selected,selected[1:]):
                if (low['control_rate']-3)*(high['control_rate']-3)<=0:
                    brackets.append(dict(profile=profile,p_low=low['scenario']['p'],p_high=high['scenario']['p'],
                                         interpretation='Adjacent sampled point estimates bracket target; no probability selected.'))
    if args.module=='sensitivity':
        for profile in ('clustered','imperfect'):
            center=next(i for i,s in enumerate(scenarios) if s['id']==profile+'_center')
            for label in ('low','high'):
                idx=next(i for i,s in enumerate(scenarios) if s['id']==profile+'_'+label)
                value=fraction[idx]-fraction[center]
                contrasts.append(dict(name=profile+'_'+label+'_minus_center_F',
                    estimate=None if not np.isfinite(value) else float(value),
                    uncertainty=interval(samples[:,idx,2]-samples[:,center,2],family)))
    report=dict(module=args.module,seeds=len(seeds),records=len(files),bootstrap_draws=draws,
                bootstrap_seed=20261010,method='Whole-seed vector resampling independently within A/B, synchronized across conditions.',
                scenarios=entries,primary_contrasts=contrasts,descriptive_adjacent_brackets=brackets,
                interpretation='Fixed synthetic follow-up; no refit, equilibrium, or global probability-infeasibility claim.',
                source_sha256=digest(Path(__file__)),config_sha256=digest(HERE/'configs'/(args.module+'.json')),
                manifest_sha256=digest(args.results/'manifest.json'),input_sha256=file_hashes)
    verify()
    args.output.mkdir(parents=True,exist_ok=True);target=args.output/'report.json'
    if target.exists():raise ValueError('Analysis output already exists; do not overwrite')
    temp=target.with_suffix('.json.tmp');temp.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n');temp.replace(target)
    print(args.module+' complete bank analyzed')


if __name__=='__main__':main()
