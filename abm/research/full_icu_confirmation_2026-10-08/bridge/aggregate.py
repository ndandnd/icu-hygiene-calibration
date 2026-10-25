#!/usr/bin/env python3
"""Validate all frozen bridge shards, then aggregate individual background rows.

No simulations. Refuse incomplete banks before opening any outcome record.
Analysis is separately hashed and does not modify frozen simulation sources.
"""
from copy import deepcopy
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parent
NEW = HERE.parent
ABM = HERE.parents[2]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def backend():
    sys.path.insert(0, str(HERE))
    return load_module('_bridge_reviewed_aggregate_backend', HERE/'run.py')


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def regular(path):
    if not path.is_file() or path.is_symlink(): raise ValueError('Missing regular artifact: '+str(path))


def nonnegative(value, label, integer=False):
    if integer:
        if type(value) is not int or value < 0: raise ValueError('Invalid integer '+label)
    elif isinstance(value, bool) or not isinstance(value, (int,float)) or not math.isfinite(value) or value < 0:
        raise ValueError('Invalid nonnegative '+label)
    return value


def close(a,b):
    return math.isclose(a,b,rel_tol=2e-10,abs_tol=1e-9)


def equal_payload(actual, expected):
    """Allow insignificant float rounding between recorded runtime builds."""
    if isinstance(expected,dict):
        return isinstance(actual,dict) and set(actual)==set(expected) and all(equal_payload(actual[k],v) for k,v in expected.items())
    if isinstance(expected,list):
        return isinstance(actual,list) and len(actual)==len(expected) and all(equal_payload(a,b) for a,b in zip(actual,expected))
    if isinstance(expected,float):
        return type(actual) in (int,float) and math.isfinite(actual) and close(actual,expected)
    return type(actual) is type(expected) and actual==expected


def validate_row(row, seed, config, identity, fixed_model):
    if row.get('seed') != seed or row.get('status')!='complete' or row.get('identity_sha256')!=identity:
        raise ValueError('Invalid completed row identity/seed/status')
    first=config['profiles'][0]
    kwargs=dict(fixed_model,**config['model_overrides'],seed=seed,
        hygiene_opportunities=config['opportunities'],hygiene_autocorrelation=first['rho'],hygiene_efficacy=first['e'])
    if row['model_kwargs']!=kwargs:raise ValueError('Row model kwargs mismatch')
    if row['executed_ticks']!=config['warmup_ticks']+config['observation_ticks']:
        raise ValueError('Row horizon mismatch')
    nonnegative(row['elapsed_seconds'],'elapsed seconds')
    pd=nonnegative(row['observation_patient_days'],'patient-days')
    sd=nonnegative(row['observation_susceptible_patient_days'],'susceptible-days')
    maximum=30*config['observation_ticks']/2592
    if pd>maximum+1e-8 or sd>pd+1e-8:raise ValueError('Exposure bounds failed')
    for value in (pd,sd):
        if abs(value*2592-round(value*2592))>1e-6:raise ValueError('Exposure is not a tick count')
    for key in ('observation_admissions','observation_imported_colonizations','total_admissions','total_imported_colonizations'):
        nonnegative(row[key],key,integer=True)
    if row['total_admissions'] < 30: raise ValueError('Initial admissions missing')
    if row['observation_imported_colonizations']>row['observation_admissions'] or row['total_imported_colonizations']>row['total_admissions'] or row['observation_admissions']>row['total_admissions'] or row['observation_imported_colonizations']>row['total_imported_colonizations']:
        raise ValueError('Admission/import conservation failed')
    event_counts=row['event_counts']
    if set(event_counts)-{'patient','environment','opportunity','reset'}:raise ValueError('Unknown event counter')
    for value in event_counts.values():nonnegative(value,'event counter',True)
    ordered=row['ordered_bridge_events_sha256']
    if not isinstance(ordered,str) or len(ordered)!=64 or any(c not in '0123456789abcdef' for c in ordered):
        raise ValueError('Invalid ordered-event digest')
    if set(row['profiles'])!={p['id'] for p in config['profiles']}:raise ValueError('Profile membership mismatch')
    labels={'0','1','2','>=3'}
    previous_counts=None
    for profile in config['profiles']:
        arms=row['profiles'][profile['id']]
        if set(arms)!={'baseline','policy'}:raise ValueError('Arm membership mismatch')
        for name in ('baseline','policy'):
            arm=arms[name];k=nonnegative(arm['K'],'K');warm=nonnegative(arm['K_from_warmup_pickups'],'warmup K')
            raw=nonnegative(arm['raw_pairs_without_reset'],'raw pairs',True)
            if warm>k+1e-8 or k>raw+1e-8:raise ValueError('Coefficient/warmup/raw-pair bounds failed')
            hist=arm['K_by_opportunities'];rh=arm['raw_pairs_without_reset_by_opportunities']
            if not isinstance(hist,dict) or not isinstance(rh,dict) or set(hist)!=labels or set(rh)!=labels:
                raise ValueError('Expected cap2 histograms')
            for key in labels:
                nonnegative(hist[key],'K bin');nonnegative(rh[key],'raw bin',True)
                if hist[key]>rh[key]+1e-8:raise ValueError('Coefficient bin exceeds raw pair count')
            if not close(sum(hist.values()),k) or sum(rh.values())!=raw:
                raise ValueError('Histogram conservation failed')
            worker=arm['K_by_worker'];laws=arm['completion_law_endpoint']
            valid_workers={str(i) for i in range(-19,0)}
            if not isinstance(worker,dict) or not isinstance(laws,dict) or set(worker)-valid_workers or set(laws)-valid_workers:
                raise ValueError('Worker IDs invalid')
            for value in worker.values():nonnegative(value,'worker K')
            if not close(sum(worker.values()),k):raise ValueError('Worker coefficient conservation failed')
            for law in laws.values():
                if not isinstance(law,list) or len(law)!=3:raise ValueError('Invalid completion law')
                for mass in law:nonnegative(mass,'completion-law mass')
                if not close(sum(law),1.):raise ValueError('Completion-law mass is not one')
            counts=arm['event_counts']
            if set(counts)-{'opportunities_observed','opportunities_warmup','implicit_resets','imported_carrier_contacts','susceptible_contacts'}:
                raise ValueError('Unknown coefficient counter')
            for value in counts.values():nonnegative(value,'coefficient counter',True)
            if counts.get('imported_carrier_contacts',0)+counts.get('susceptible_contacts',0)!=event_counts.get('patient',0) or counts.get('opportunities_observed',0)+counts.get('opportunities_warmup',0)!=event_counts.get('opportunity',0) or counts.get('implicit_resets',0)!=event_counts.get('reset',0):
                raise ValueError('Event/coefficient counter conservation failed')
            if previous_counts is not None and counts!=previous_counts:raise ValueError('Background differs across profiles/arms')
            previous_counts=counts
        if arms['policy']['K']>arms['baseline']['K']+1e-8:
            raise ValueError('Increasing completion produced larger K')


def preflight(root, book):
    """Check complete file membership for all38 partitions before any row reads."""
    for site,parts in book['partitions'].items():
        base=root/f'results/{NEW.name}/bridge_{site}'
        if not base.is_dir():raise ValueError('Incomplete bridge site: '+site)
        actual={p.name for p in base.glob('part_*') if p.is_dir()}
        expected={f"part_{p['index']:02d}" for p in parts}
        if actual!=expected:raise ValueError('Incomplete/extra partition scopes: '+site)
        for part in parts:
            directory=root/part['output_scope']
            if directory.resolve()!=directory:raise ValueError('Symlinked partition scope')
            for name in ('bridge_binding.json','manifest.json','summary.json'):
                regular(directory/name)
            complete=directory/'completed'
            if not complete.is_dir() or complete.is_symlink():raise ValueError('Missing completed directory')
            filenames={p.name for p in complete.glob('*.json')}
            if filenames!={f'{s}.json' for s in part['seeds']}:
                raise ValueError('Incomplete/extra completed seeds')
            for seed in part['seeds']:regular(complete/f'{seed}.json')
    # All immutable completion summaries must confirm the fixed sample before
    # any individual outcome record is opened. Partial progress files do not count.
    for parts in book['partitions'].values():
        for part in parts:
            marker=json.loads((root/part['output_scope']/'summary.json').read_text())
            if marker.get('status')!='complete' or marker.get('seeds_completed')!=len(part['seeds']) or marker.get('seeds_planned')!=len(part['seeds']):
                raise ValueError('Shard completion gate failed')


def validate_outer(outer,part,site,book,freeze,freeze_sha,root):
    prefix=f'research/{NEW.name}';config_path=root/part['config_path']
    expected={'study':NEW.name,'module':'bridge','freeze_identity':freeze['identity_sha256'],
        'freeze_sha256':freeze_sha,'site':site,'partition_index':part['index'],
        'partition_count':len(book['partitions'][site]),'config_path':part['config_path'],
        'config_sha256':digest(config_path),'protocol_path':f'{prefix}/bridge/PROTOCOL.md',
        'protocol_sha256':digest(root/f'{prefix}/bridge/PROTOCOL.md'),
        'seed_audit_path':f'{prefix}/seed_audit.json','seed_audit_sha256':digest(root/f'{prefix}/seed_audit.json'),
        'bridge_seedbook_path':f'{prefix}/bridge/seedbook.json',
        'bridge_seedbook_sha256':digest(root/f'{prefix}/bridge/seedbook.json'),
        'seeds':part['seeds'],'output_scope':part['output_scope'],
        'backend':f'{prefix}/bridge/run.py',
        'dispatch_config_path':f'{prefix}/bridge/configs/bridge_{site}.json',
        'dispatch_config_sha256':digest(root/f'{prefix}/bridge/configs/bridge_{site}.json'),
        'dispatch_output_scope':f'results/{NEW.name}/bridge_{site}'}
    if any(outer.get(k)!=v for k,v in expected.items()):raise ValueError('Outer binding mismatch')
    science={f'{prefix}/bridge/{n}' for n in ('run.py','coefficient.py','background.py')}
    science|={f'src/{n}' for n in ('CPEmodel.py','agents.py','routes.py')}
    if outer.get('scientific_sources')!={n:freeze['files'][n] for n in science}:
        raise ValueError('Outer scientific source chain mismatch')


def collect(root,freeze_path):
    """Complete-design gate, then validation and loading. Do not call on running banks."""
    root=root.resolve()
    verifier=load_module('_bridge_aggregate_freeze',root/f'research/{NEW.name}/array_worker.py')
    freeze=verifier.verify_freeze(freeze_path,root)
    if freeze.get('phase')!='bridge':raise ValueError('Bridge phase freeze required')
    prefix=f'research/{NEW.name}';book_path=root/f'{prefix}/bridge/seedbook.json'
    book=json.loads(book_path.read_text());audit=json.loads((root/f'{prefix}/seed_audit.json').read_text())
    partition_worker=load_module('_bridge_aggregate_partition',root/f'{prefix}/bridge/array_worker.py')
    partition_worker.validate_partition(book,audit,'local',0)
    for name in (f'{prefix}/bridge/seedbook.json',f'{prefix}/bridge/PROTOCOL.md',f'{prefix}/seed_audit.json',f'{prefix}/seedbook.json',f'{prefix}/runtime_validation.json'):
        if freeze['files'].get(name)!=digest(root/name):raise ValueError('Required aggregate input not frozen')
    primary=json.loads((root/f'{prefix}/seedbook.json').read_text())
    prior_declared={s for banks in primary['seed_blocks'].values() for values in banks.values() for s in values}
    if prior_declared.intersection(book['all_seeds']):raise ValueError('Bridge/primary seed overlap')
    all_config_path=root/f'{prefix}/bridge/configs/bridge_all.json'
    if freeze['files'].get(str(all_config_path.relative_to(root)))!=digest(all_config_path):
        raise ValueError('All-background config not frozen')
    all_config=json.loads(all_config_path.read_text())
    if all_config!=dict(book['scientific_config'],seeds=book['all_seeds']):raise ValueError('All-config membership/science mismatch')
    needed=[f'{prefix}/array_worker.py',f'{prefix}/bridge/array_worker.py',
        *[f'{prefix}/bridge/{n}' for n in ('run.py','coefficient.py','background.py','test_bridge.py')],
        *[f'src/{n}' for n in ('CPEmodel.py','agents.py','routes.py')],
        *[p['config_path'] for parts in book['partitions'].values() for p in parts],
        *[f'{prefix}/bridge/configs/bridge_{site}.json' for site in ('local','cluster')]]
    if any(name not in freeze['files'] for name in needed):raise ValueError('Required scientific/binding/config file omitted from freeze')
    preflight(root,book)
    reviewed=backend();config=reviewed.normalize(all_config)
    expected_sources={n:freeze['files'][f'{prefix}/bridge/{n}'] for n in ('run.py','coefficient.py','background.py','test_bridge.py')}
    expected_sources.update({f'src/{n}':freeze['files'][f'src/{n}'] for n in ('CPEmodel.py','agents.py','routes.py')})
    for name in ('run.py','coefficient.py','background.py','test_bridge.py'):
        if digest(HERE/name)!=expected_sources[name]:raise ValueError('Loaded reviewed analysis backend differs from frozen science')
    runtime_reference=json.loads((root/f'{prefix}/runtime_validation.json').read_text())
    rows_by_seed={};receipts=[];freeze_sha=digest(freeze_path)
    for site,parts in book['partitions'].items():
        for part in parts:
            directory=root/part['output_scope']
            outer=json.loads((directory/'bridge_binding.json').read_text())
            validate_outer(outer,part,site,book,freeze,freeze_sha,root)
            manifest=json.loads((directory/'manifest.json').read_text())
            identity=manifest['identity_sha256'];payload={k:v for k,v in manifest.items() if k!='identity_sha256'}
            if hashlib.sha256(canonical(payload).encode()).hexdigest()!=identity:raise ValueError('Inner manifest identity mismatch')
            raw_config=json.loads((root/part['config_path']).read_text());resolved=reviewed.normalize(raw_config)
            if raw_config!=dict(book['scientific_config'],seeds=part['seeds']):
                raise ValueError('Partition executable config differs from declared membership/science')
            if manifest.get('study')!=NEW.name+'_bridge' or manifest.get('source_sha256')!=expected_sources or manifest.get('fixed_model')!=reviewed.FIXED_MODEL or manifest.get('input_config')!=raw_config or manifest.get('resolved_config')!=resolved:
                raise ValueError('Inner science/config identity mismatch')
            if manifest.get('versions')!=runtime_reference[f'{site}_runtime']:
                raise ValueError('Runtime differs from declared site validation')
            shard=[];hashes={}
            for seed in part['seeds']:
                path=directory/'completed'/f'{seed}.json';row=json.loads(path.read_text())
                validate_row(row,seed,resolved,identity,reviewed.FIXED_MODEL)
                if seed in rows_by_seed:raise ValueError('Duplicate background across shards')
                rows_by_seed[seed]=row;shard.append(row);hashes[path.name]=digest(path)
            stored_summary=json.loads((directory/'summary.json').read_text())
            expected_summary=reviewed.summarize(shard,resolved)
            expected_summary.update(identity_sha256=identity,status='complete')
            if not equal_payload(stored_summary,expected_summary):raise ValueError('Shard summary fails individual-row reconciliation')
            receipts.append({'site':site,'partition_index':part['index'],'output_scope':part['output_scope'],
                'identity':identity,'manifest_sha256':digest(directory/'manifest.json'),'jobs':len(shard),
                'completed_records_sha256':hashes,'bridge_binding_sha256':digest(directory/'bridge_binding.json'),
                'summary_sha256':digest(directory/'summary.json'),'runtime':manifest['versions'],
                'config_path':part['config_path'],'config_sha256':digest(root/part['config_path'])})
    if set(rows_by_seed)!=set(book['all_seeds']) or len(rows_by_seed)!=256:raise ValueError('Final exact256 membership failed')
    sources={f'{prefix}/bridge/{n}':expected_sources[n] for n in ('run.py','coefficient.py','background.py')}
    sources.update({f'src/{n}':expected_sources[f'src/{n}'] for n in ('CPEmodel.py','agents.py','routes.py')})
    provenance={'freeze_sha256':freeze_sha,'freeze_identity':freeze['identity_sha256'],
        'config_path':str(all_config_path.relative_to(root)),'config_sha256':digest(all_config_path),
        'scientific_sources':sources,'source_sha256':sources,'banks':receipts,
        'seedbook_sha256':digest(book_path),'runtime_validation_sha256':digest(root/f'{prefix}/runtime_validation.json')}
    return [rows_by_seed[s] for s in book['all_seeds']],config,book,provenance


# Shared complete-design loader for release tooling; no analysis is published here.
load_complete = collect

def difference(mean_a,se_a,mean_b,se_b):
    estimate=mean_b-mean_a
    se=None if se_a is None or se_b is None else math.hypot(se_a,se_b)
    return {'estimate':estimate,'monte_carlo_se':se,
        'normal_approximation_95_interval':None if se is None else [estimate-1.96*se,estimate+1.96*se]}


def summarize_complete(rows,config,book,provenance,analysis_hash):
    reviewed=backend()
    if len(rows)!=256 or [r['seed'] for r in rows]!=book['all_seeds']:
        raise ValueError('Exact ordered256 backgrounds required')
    result=reviewed.summarize(rows,config)
    result.update(provenance,schema_version=1,study=NEW.name,module='bridge',status='complete',n_backgrounds=256,seeds=book['all_seeds'],
                  analysis_code_sha256=analysis_hash,statistical_banks={},bank_B_minus_A={})
    by_seed={r['seed']:r for r in rows}
    for name,seeds in book['seed_blocks'].items():
        result['statistical_banks'][name]=reviewed.summarize([by_seed[s] for s in seeds],dict(config,seeds=seeds))
        result['statistical_banks'][name].update(status='complete',seeds=seeds)
    for profile in config['profiles']:
        key=profile['id'];a=result['statistical_banks']['A']['profiles'][key];b=result['statistical_banks']['B']['profiles'][key]
        contrasts={name:difference(a[name]['mean'],a[name]['monte_carlo_se'],b[name]['mean'],b[name]['monte_carlo_se'])
            for name in ('K0','K1','paired_K1_minus_K0')}
        if a['limiting_response'] is None or b['limiting_response'] is None:
            contrasts['limiting_response']={'estimate':None,'monte_carlo_se':None,'normal_approximation_95_interval':None}
        else:
            contrasts['limiting_response']=difference(a['limiting_response'],a['limiting_response_delta_monte_carlo_se'],b['limiting_response'],b['limiting_response_delta_monte_carlo_se'])
        result['bank_B_minus_A'][key]=contrasts
    result['bank_comparison_note']='Independent statistical banks; B−A SE combines bank variances. Delta response SE preserves within-background arm covariance. Descriptive diagnostics do not identify a compute-site effect or establish equivalence.'
    return result


def publish(summary,output):
    raw=json.dumps(summary,indent=2,sort_keys=True,allow_nan=False)+'\n'
    binding={'schema_version':1,'module':'bridge','status':'complete','n_backgrounds':256,
        'summary_sha256':hashlib.sha256(raw.encode()).hexdigest(),
        **{key:summary[key] for key in ('config_sha256','seeds','freeze_sha256','freeze_identity','scientific_sources','source_sha256','analysis_code_sha256')},
        'summary_filename':'summary.json','config_path':summary['config_path']}
    if (output/'summary.json').exists() or (output/'binding.json').exists():raise FileExistsError('Aggregate outputs are immutable')
    output.mkdir(parents=True,exist_ok=True)
    with (output/'summary.json').open('x') as stream:stream.write(raw)
    with (output/'binding.json').open('x') as stream:stream.write(json.dumps(binding,indent=2,sort_keys=True,allow_nan=False)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--abm-root',type=Path,default=ABM)
    parser.add_argument('--freeze',type=Path,default=NEW/'freeze_bridge.json')
    parser.add_argument('--output',type=Path,default=HERE/'aggregate_report')
    args=parser.parse_args()
    rows,config,book,provenance=collect(args.abm_root,args.freeze)
    summary=summarize_complete(rows,config,book,provenance,digest(Path(__file__)))
    publish(summary,args.output)
    print(canonical({'status':'complete','n_backgrounds':256,'partitions':len(provenance['banks']),
                     'output':str(args.output),'summary_sha256':digest(args.output/'summary.json')}))


if __name__=='__main__':main()
