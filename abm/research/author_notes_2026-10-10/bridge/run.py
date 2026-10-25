#!/usr/bin/env python3
"""Deterministic original p=0 replay; 22 independent coefficient comparisons."""
import argparse
import fcntl
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import signal
import statistics
import sys
import time

try:
    from .background import ABM, CORE_SHA256, FIXED_MODEL, make_background_model
    from .coefficient import stationary_q2
except ImportError:
    from background import ABM, CORE_SHA256, FIXED_MODEL, make_background_model
    from coefficient import stationary_q2

HERE = Path(__file__).resolve().parent
LEGACY_PROFILES = [('iid_perfect', 0., 1.), ('persistent_perfect', .5, 1.),
                   ('iid_imperfect', 0., .8)]


def coefficient_pairs():
    pairs = []
    for h0, h1, profiles in ((.9, .95, [(None, r, e) for r in (0., .25, .5, .75) for e in (1., .9, .8, .6)]),
                             (.8, .9, LEGACY_PROFILES), (.95, .99, LEGACY_PROFILES)):
        for _, rho, e in profiles:
            legacy = next((name for name, r, eff in LEGACY_PROFILES if (r, eff) == (rho, e)), None)
            pair = dict(id=f'h{h0:g}_to_{h1:g}_rho{rho:g}_e{e:g}', h0=h0, h1=h1, rho=rho, e=e)
            if (h0, h1) == (.9, .95) and legacy:
                pair['legacy_id'] = legacy
            pairs.append(pair)
    return pairs


PAIRS = coefficient_pairs()
SCIENCE = dict(warmup_ticks=90*2592, observation_ticks=56*2592,
               background_h=.90, background_rho=0., background_e=1.,
               opportunities='zone_separate', prob_patient_sick=.03,
               raw_histogram_cap=12, pairs=PAIRS)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_references():
    return json.loads((HERE/'reference_index.json').read_text())


def normalize(raw):
    if not isinstance(raw, dict) or set(raw)-{'seeds', 'description'}:
        raise ValueError('Only seeds and description are configurable; scientific parameters are fixed')
    seeds = raw.get('seeds')
    bank = [row['seed'] for row in read_references()['seeds']]
    if not isinstance(seeds, list) or not seeds or any(type(s) is not int or s not in bank for s in seeds) or len(set(seeds)) != len(seeds):
        raise ValueError('Unique original-bank seed subset required')
    return dict(seeds=seeds, science=SCIENCE)


def source_hashes(abm_root=ABM):
    abm_root = Path(abm_root)
    values = {str((HERE/name).relative_to(ABM)): sha(HERE/name)
              for name in ('run.py','background.py','coefficient.py','test_bridge.py','benchmark.py','array_worker.py','PROTOCOL.md','reference_index.json')}
    refs = read_references()
    for path, expected in refs['old_source_sha256'].items():
        if sha(abm_root/path) != expected:
            raise RuntimeError('Original source drift: '+path)
        values[path] = expected
    for name, expected in CORE_SHA256.items():
        path = 'src/'+name
        if sha(abm_root/path) != expected:
            raise RuntimeError('Frozen core changed: '+path)
        values[path] = expected
    return values


def verify_reference_files(config, abm_root=ABM):
    refs = {row['seed']: row for row in read_references()['seeds']}
    for seed in config['seeds']:
        row = refs[seed]
        if sha(Path(abm_root)/row['path']) != row['sha256']:
            raise RuntimeError('Original completed record changed: '+row['path'])
    return {seed: refs[seed] for seed in config['seeds']}


def build_manifest(raw, config, abm_root=ABM):
    sources = source_hashes(abm_root)
    refs = verify_reference_files(config, abm_root)
    science_identity = hashlib.sha256(canonical(dict(science=SCIENCE, sources=sources)).encode()).hexdigest()
    versions = dict(python=sys.version, platform=platform.platform(),
                    **{name: importlib.metadata.version(name) for name in ('Mesa','numpy','pandas','networkx','matplotlib')})
    manifest = dict(schema_version=1, study='author_notes_2026-10-10_bridge', input_config=raw,
                    resolved_config=config, fixed_model=dict(FIXED_MODEL, prob_patient_sick=.03),
                    source_sha256=sources, reference_records={str(s): dict(path=r['path'],sha256=r['sha256']) for s,r in refs.items()},
                    scientific_identity_sha256=science_identity, versions=versions,
                    background_note='Original literal-seed p=0 law h=.90/rho=0/e=1 for all 22 comparisons; no new epidemiological sampling.',
                    coefficient_note='Each pair uses its own h0 warmup/control and h1 after warmup; completion memory persists across resets.',
                    checkpoint_note='Complete seed files atomic and immutable. Incomplete seeds restart from initialization; no event traces stored.')
    manifest['identity_sha256'] = hashlib.sha256(canonical(manifest).encode()).hexdigest()
    return manifest


def atomic_json(path, value, immutable=False):
    path = Path(path)
    if immutable and path.exists():
        if json.loads(path.read_text()) != value:
            raise ValueError('Immutable output differs: '+str(path))
        return
    temporary = path.with_suffix(path.suffix+'.tmp')
    with temporary.open('w') as handle:
        handle.write(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)
    descriptor = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def seal_seed(row):
    row['payload_sha256'] = hashlib.sha256(canonical(row).encode()).hexdigest()
    return row


def verify_seed_payload(row):
    expected = row.get('payload_sha256')
    payload = {key: value for key, value in row.items() if key != 'payload_sha256'}
    if expected != hashlib.sha256(canonical(payload).encode()).hexdigest():
        raise ValueError('Completed seed payload digest mismatch')


def simulate_seed(seed, factory, deadline=None, *, warmup_ticks=SCIENCE['warmup_ticks'],
                  observation_ticks=SCIENCE['observation_ticks']):
    kwargs = dict(FIXED_MODEL, prob_patient_sick=.03, seed=seed,
                  hygiene_opportunities='zone_separate', hygiene_autocorrelation=0., hygiene_efficacy=1.)
    model = factory(warmup_ticks=warmup_ticks, pairs=PAIRS, **kwargs)
    boundary = None
    started = time.monotonic()
    for tick in range(warmup_ticks+observation_ticks+1):
        if tick == warmup_ticks:
            boundary = (model.patient_days, model.susceptible_patient_days, model.cumul_patients,
                        model.imported_colonizations)
        if tick == warmup_ticks+observation_ticks:
            break
        if tick % 128 == 0 and deadline is not None and time.monotonic() >= deadline:
            raise TimeoutError('Incomplete seed: restart from initialization')
        model.step()
    result = model.bridge_result()
    result.update(seed=seed, status='complete', model_kwargs=kwargs,
                  elapsed_seconds=time.monotonic()-started, executed_ticks=model.total_ticks,
                  observation_patient_days=model.patient_days-boundary[0],
                  observation_susceptible_patient_days=model.susceptible_patient_days-boundary[1],
                  observation_admissions=model.cumul_patients-boundary[2],
                  observation_imported_colonizations=model.imported_colonizations-boundary[3],
                  total_admissions=model.cumul_patients,total_imported_colonizations=model.imported_colonizations)
    return result


def validate_replay(row, reference, abm_root=ABM):
    path = Path(abm_root)/reference['path']
    if sha(path) != reference['sha256']:
        raise RuntimeError('Original reference changed during replay')
    old = json.loads(path.read_text())
    for key in ('seed','status','model_kwargs','event_counts','executed_ticks',
                'ordered_bridge_events_sha256','observation_patient_days',
                'observation_susceptible_patient_days','observation_admissions',
                'observation_imported_colonizations','total_admissions','total_imported_colonizations'):
        if row[key] != old[key]:
            raise RuntimeError(f'Original replay mismatch: {key}, seed {row["seed"]}')
    for pair in PAIRS:
        legacy = pair.get('legacy_id')
        if legacy:
            for arm in ('baseline','policy'):
                # Identical recurrence/order should reproduce K bit for bit.
                if row['pairs'][pair['id']][arm] != old['profiles'][legacy][arm]:
                    raise RuntimeError(f'Original coefficient mismatch: {legacy}/{arm}, seed {row["seed"]}')
    row['original_replay_validation'] = dict(status='pass', reference_sha256=reference['sha256'],
                                              matched_legacy_profiles=3, ordered_event_hash_equal=True,
                                              coefficients_bitwise_equal=True)


def mean_uncertainty(values):
    mean = statistics.fmean(values)
    se = statistics.stdev(values)/math.sqrt(len(values)) if len(values)>1 else None
    return dict(mean=mean, monte_carlo_se=se,
                normal_approximation_95_interval=None if se is None else [mean-1.96*se,mean+1.96*se])


def summarize(rows, config):
    result = dict(seeds_completed=len(rows), seeds_planned=len(config['seeds']),pairs={},
                  uncertainty_note='Descriptive seed-level normal/delta Monte Carlo intervals; conditional survival integrated exactly.')
    if not rows:
        return result
    for pair in PAIRS:
        key=pair['id']
        x=[row['pairs'][key]['baseline']['K'] for row in rows]
        y=[row['pairs'][key]['policy']['K'] for row in rows]
        xm,ym=statistics.fmean(x),statistics.fmean(y)
        ratio=None if xm==0 else ym/xm
        response=None if ratio is None else 1-ratio
        residual=None if ratio is None else [(b-ratio*a)/xm for a,b in zip(x,y)]
        se=None if residual is None or len(rows)<2 else statistics.stdev(residual)/math.sqrt(len(rows))
        out=dict(parameters=pair,K0=mean_uncertainty(x),K1=mean_uncertainty(y),
                 paired_K1_minus_K0=mean_uncertainty([b-a for a,b in zip(x,y)]),
                 K1_over_K0=ratio,limiting_response=response,
                 limiting_response_delta_monte_carlo_se=se,
                 limiting_response_normal_95_interval=None if se is None else [response-1.96*se,response+1.96*se])
        q0,q1=(stationary_q2(pair[h],pair['rho'],pair['e']) for h in ('h0','h1'))
        out['stationary_q2_response_for_comparison_only']=None if q0==0 else 1-q1/q0
        for arm in ('baseline','policy'):
            out[arm+'_K_from_warmup_pickups']=mean_uncertainty([r['pairs'][key][arm]['K_from_warmup_pickups'] for r in rows])
            hist=rows[0]['pairs'][key][arm]['K_by_opportunities']
            if hist is not None:
                out[arm+'_K_by_opportunities_mean']={label:statistics.fmean(r['pairs'][key][arm]['K_by_opportunities'][label] for r in rows) for label in hist}
        result['pairs'][key]=out
    result['raw_pairs_without_reset_by_opportunities_total']={label:sum(r['raw_pairs_without_reset_by_opportunities'][label] for r in rows)
                                                             for label in rows[0]['raw_pairs_without_reset_by_opportunities']}
    result['observation_event_counts_by_role_total']={role:{kind:sum(r['observation_event_counts_by_role'][role][kind] for r in rows)
                                                           for kind in ('opportunity','patient','surface','reset')}
                                                      for role in ('Nurse','Dr','XrayDr')}
    result['observation_patient_days']=mean_uncertainty([r['observation_patient_days'] for r in rows])
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--abm-root',type=Path,default=ABM)
    parser.add_argument('--wall-seconds',type=float,default=3300.)
    parser.add_argument('--freeze-only',action='store_true')
    args=parser.parse_args(argv)
    if not math.isfinite(args.wall_seconds) or not 0 < args.wall_seconds <= 3300:
        parser.error('wall-seconds must be finite in (0,3300]')
    raw=json.loads(args.config.read_text()); config=normalize(raw)
    manifest=build_manifest(raw,config,args.abm_root)
    refs=verify_reference_files(config,args.abm_root)
    args.output.mkdir(parents=True,exist_ok=True)
    completed=args.output/'completed'; completed.mkdir(exist_ok=True)
    def interrupted(signum,frame):
        raise InterruptedError('Incomplete seed interrupted; restart from initialization')
    with (args.output/'.run.lock').open('a+') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        atomic_json(args.output/'manifest.json',manifest,immutable=True)
        if args.freeze_only:
            print('Frozen '+manifest['identity_sha256']);return 0
        previous={sig:signal.signal(sig,interrupted) for sig in (signal.SIGTERM,signal.SIGUSR1,signal.SIGINT)}
        rows=[]
        try:
            factory=make_background_model(args.abm_root)
            deadline=time.monotonic()+args.wall_seconds
            for seed in config['seeds']:
                path=completed/f'{seed}.json'
                if path.exists():
                    row=json.loads(path.read_text())
                    verify_seed_payload(row)
                    if row.get('identity_sha256')!=manifest['identity_sha256'] or row.get('seed')!=seed or row.get('status')!='complete':
                        raise ValueError('Invalid completed seed: '+str(path))
                    validate_replay(row,refs[seed],args.abm_root)
                else:
                    if time.monotonic()>=deadline:break
                    row=simulate_seed(seed,factory,deadline)
                    validate_replay(row,refs[seed],args.abm_root)
                    if source_hashes(args.abm_root)!=manifest['source_sha256']:
                        raise RuntimeError('Source changed during seed')
                    row['identity_sha256']=manifest['identity_sha256']
                    seal_seed(row)
                    atomic_json(path,row,immutable=True)
                rows.append(row)
        except (TimeoutError,InterruptedError):
            pass
        finally:
            for sig,handler in previous.items():signal.signal(sig,handler)
        summary=summarize(rows,config); complete=len(rows)==len(config['seeds'])
        summary.update(identity_sha256=manifest['identity_sha256'],status='complete' if complete else 'partial')
        atomic_json(args.output/('summary.json' if complete else 'progress.json'),summary,immutable=complete)
        print(canonical(summary))
        return 0 if complete else 75


if __name__=='__main__':sys.exit(main())
