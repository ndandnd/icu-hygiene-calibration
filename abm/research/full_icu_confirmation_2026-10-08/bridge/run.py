#!/usr/bin/env python3
"""Portable p=0 coefficient runner with immutable per-seed files and resume.

No positive-p simulations. Different completion profiles are Rao-Blackwellized
on one p=0 background per seed. --freeze-only prepares a reviewable manifest.
"""
import argparse
import fcntl
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import statistics
import sys
import time

from background import ABM, CORE_SHA256, FIXED_MODEL, make_background_model
from coefficient import probability, stationary_q2

HERE = Path(__file__).resolve().parent
DEFAULT_PROFILES = [dict(id='iid_perfect', rho=0., e=1.),
                    dict(id='persistent_perfect', rho=.5, e=1.),
                    dict(id='iid_imperfect', rho=0., e=.8)]


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def normalize(raw):
    allowed = {'seeds', 'warmup_ticks', 'observation_ticks', 'h0', 'h1', 'profiles',
               'opportunities', 'histogram_cap', 'description', 'model_overrides'}
    if not isinstance(raw, dict) or set(raw) - allowed:
        raise ValueError('Unsupported bridge config keys')
    config = dict(raw)
    seeds = config.get('seeds')
    if not isinstance(seeds, list) or not seeds or any(type(s) is not int or s < 0 for s in seeds) or len(seeds) != len(set(seeds)):
        raise ValueError('Unique nonnegative integer seeds required')
    for key, default in (('warmup_ticks', 90*2592), ('observation_ticks', 28*2592)):
        config.setdefault(key, default)
        if type(config[key]) is not int or config[key] < (0 if key == 'warmup_ticks' else 1):
            raise ValueError(key + ' must be a valid integer tick count')
    for key, default in (('h0', .9), ('h1', .95)):
        config[key] = probability(config.get(key, default), key)
    config.setdefault('opportunities', 'zone_separate')
    if config['opportunities'] not in ('separate', 'shared', 'zone_separate'):
        raise ValueError('Unknown opportunity convention')
    config.setdefault('histogram_cap', None)
    cap = config['histogram_cap']
    if cap is not None and (type(cap) is not int or cap < 0):
        raise ValueError('histogram_cap must be nonnegative integer or null')
    profiles = config.setdefault('profiles', DEFAULT_PROFILES)
    if not isinstance(profiles, list) or not profiles:
        raise ValueError('Nonempty profiles required')
    identifiers = set()
    for profile in profiles:
        if not isinstance(profile, dict) or set(profile) != {'id', 'rho', 'e'}:
            raise ValueError('Each profile requires exactly id, rho, e')
        if not isinstance(profile['id'], str) or not profile['id'] or profile['id'] in identifiers:
            raise ValueError('Unique nonempty profile IDs required')
        identifiers.add(profile['id'])
        probability(profile['e'], 'e')
        if probability(profile['rho'], 'rho') == 1:
            raise ValueError('rho must be below 1')
    overrides = config.setdefault('model_overrides', {})
    if not isinstance(overrides, dict) or set(overrides) - {'prob_patient_sick'}:
        raise ValueError('Only prob_patient_sick may override fixed background model')
    if 'prob_patient_sick' in overrides:
        probability(overrides['prob_patient_sick'], 'prob_patient_sick')
    return config


def source_hashes(abm_root):
    paths = [HERE / name for name in ('run.py', 'coefficient.py', 'background.py', 'test_bridge.py')]
    values = {str(path.relative_to(HERE)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    for name, expected in CORE_SHA256.items():
        actual = hashlib.sha256((abm_root / 'src' / name).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError('Frozen core hash changed: ' + name)
        values['src/' + name] = actual
    return values


def build_manifest(raw, config, abm_root):
    runtime = dict(python=sys.version, platform=platform.platform())
    for package in ('Mesa', 'numpy', 'pandas', 'networkx'):
        runtime[package] = importlib.metadata.version(package)
    result = dict(schema_version=1, study='full_icu_confirmation_2026-10-08_bridge',
                  input_config=raw, resolved_config=config, fixed_model=FIXED_MODEL,
                  source_sha256=source_hashes(abm_root), versions=runtime,
                  estimand='E_B sum imported-source same-worker bridges c_i d_j survival_h(i,j)',
                  background_law_note='One clean p=0 background under h0 and first profile. All profiles share its ideal background law; this is not alternative-arm literal-seed matching.',
                  reservoir_note='Clean hands/sites at tick 0; imported patients retained. Environment contacts contribute intervening opportunities only.',
                  checkpoint_note='Completed seed files immutable; incomplete seeds restart from initialization; no biological checkpoint.')
    result['identity_sha256'] = hashlib.sha256(canonical(result).encode()).hexdigest()
    return result


def immutable_json(path, value):
    """Exclusive creation under run lock; reject any conflicting existing file."""
    if path.exists():
        if json.loads(path.read_text()) != value:
            raise ValueError('Immutable output differs: ' + str(path))
        return
    temporary = path.with_suffix(path.suffix + '.tmp')
    with temporary.open('w') as handle:
        handle.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


def simulate_seed(seed, config, factory, deadline=None):
    kwargs = dict(FIXED_MODEL, **config['model_overrides'])
    first = config['profiles'][0]
    kwargs.update(seed=seed, hygiene_opportunities=config['opportunities'],
                  hygiene_autocorrelation=first['rho'], hygiene_efficacy=first['e'])
    model = factory(warmup_ticks=config['warmup_ticks'], h0=config['h0'], h1=config['h1'],
                    profiles=config['profiles'], histogram_cap=config['histogram_cap'], **kwargs)
    boundary = None
    started = time.monotonic()
    for tick in range(config['warmup_ticks'] + config['observation_ticks'] + 1):
        if tick == config['warmup_ticks']:
            boundary = (model.patient_days, model.susceptible_patient_days, model.cumul_patients,
                        model.imported_colonizations)
        if tick == config['warmup_ticks'] + config['observation_ticks']:
            break
        if tick % 128 == 0 and deadline is not None and time.monotonic() >= deadline:
            raise TimeoutError('Incomplete seed: restart from initialization')
        model.step()
    result = model.bridge_result()
    result.update(seed=seed, status='complete', model_kwargs=kwargs,
                  elapsed_seconds=time.monotonic()-started,
                  executed_ticks=model.total_ticks,
                  observation_patient_days=model.patient_days-boundary[0],
                  observation_susceptible_patient_days=model.susceptible_patient_days-boundary[1],
                  observation_admissions=model.cumul_patients-boundary[2],
                  observation_imported_colonizations=model.imported_colonizations-boundary[3],
                  total_admissions=model.cumul_patients, total_imported_colonizations=model.imported_colonizations)
    return result


def mean_uncertainty(values):
    mean = statistics.fmean(values)
    se = statistics.stdev(values)/math.sqrt(len(values)) if len(values) > 1 else None
    return dict(mean=mean, monte_carlo_se=se,
                normal_approximation_95_interval=None if se is None else [mean-1.96*se, mean+1.96*se])


def summarize(rows, config):
    result = dict(seeds_completed=len(rows), seeds_planned=len(config['seeds']),
                  uncertainty_note='Seed-level Monte Carlo uncertainty, conditional survival integrated exactly. Normal/delta intervals are asymptotic diagnostics; very small banks do not provide reliable coverage.',
                  profiles={})
    if not rows:
        return result
    for profile in config['profiles']:
        key = profile['id']
        x = [row['profiles'][key]['baseline']['K'] for row in rows]
        y = [row['profiles'][key]['policy']['K'] for row in rows]
        xmean, ymean = statistics.fmean(x), statistics.fmean(y)
        ratio = None if xmean == 0 else ymean/xmean
        r = None if ratio is None else 1-ratio
        residual = None if ratio is None else [(b-ratio*a)/xmean for a, b in zip(x, y)]
        se = None if residual is None or len(rows) < 2 else statistics.stdev(residual)/math.sqrt(len(rows))
        arm = dict(K0=mean_uncertainty(x), K1=mean_uncertainty(y),
                   paired_K1_minus_K0=mean_uncertainty([b-a for a,b in zip(x,y)]),
                   K1_over_K0=ratio, limiting_response=r,
                   limiting_response_delta_monte_carlo_se=se,
                   limiting_response_normal_95_interval=None if se is None else [r-1.96*se, r+1.96*se],
                   stationary_q2_baseline=stationary_q2(config['h0'], profile['rho'], profile['e']),
                   stationary_q2_policy=stationary_q2(config['h1'], profile['rho'], profile['e']))
        q0, q1 = arm['stationary_q2_baseline'], arm['stationary_q2_policy']
        arm['stationary_q2_response_for_comparison_only'] = None if q0 == 0 else 1-q1/q0
        for name in ('baseline', 'policy'):
            arm[name + '_K_from_warmup_pickups'] = mean_uncertainty([row['profiles'][key][name]['K_from_warmup_pickups'] for row in rows])
            histograms = [row['profiles'][key][name]['K_by_opportunities'] for row in rows]
            if histograms[0] is not None:
                arm[name + '_K_by_opportunities_mean'] = {label: statistics.fmean(hist[label] for hist in histograms) for label in histograms[0]}
        result['profiles'][key] = arm
    result['observation_patient_days'] = mean_uncertainty([row['observation_patient_days'] for row in rows])
    result['total_event_counts'] = {key: sum(row['event_counts'].get(key, 0) for row in rows) for key in ('patient', 'environment', 'opportunity', 'reset')}
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--abm-root', default=ABM, type=Path)
    parser.add_argument('--wall-seconds', default=3300., type=float)
    parser.add_argument('--freeze-only', action='store_true')
    args = parser.parse_args(argv)
    if not math.isfinite(args.wall_seconds) or args.wall_seconds <= 0:
        parser.error('wall-seconds must be finite and positive')
    raw = json.loads(args.config.read_text())
    config = normalize(raw)
    manifest = build_manifest(raw, config, args.abm_root)
    args.output.mkdir(parents=True, exist_ok=True)
    completed = args.output / 'completed'
    completed.mkdir(exist_ok=True)
    with (args.output / '.run.lock').open('a+') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        immutable_json(args.output / 'manifest.json', manifest)
        if args.freeze_only:
            print('Frozen ' + manifest['identity_sha256'])
            return 0
        factory = make_background_model(args.abm_root)
        deadline = time.monotonic() + args.wall_seconds
        rows = []
        for seed in config['seeds']:
            path = completed / f'{seed}.json'
            if path.exists():
                row = json.loads(path.read_text())
                if row.get('identity_sha256') != manifest['identity_sha256'] or row.get('seed') != seed or row.get('status') != 'complete':
                    raise ValueError('Invalid completed seed: ' + str(path))
            else:
                if time.monotonic() >= deadline:
                    break
                try:
                    row = simulate_seed(seed, config, factory, deadline)
                except TimeoutError:
                    break
                if source_hashes(args.abm_root) != manifest['source_sha256']:
                    raise RuntimeError('Source changed during seed')
                row['identity_sha256'] = manifest['identity_sha256']
                immutable_json(path, row)
            rows.append(row)
        summary = summarize(rows, config)
        summary['identity_sha256'] = manifest['identity_sha256']
        complete = len(rows) == len(config['seeds'])
        summary['status'] = 'complete' if complete else 'partial'
        if complete:
            immutable_json(args.output / 'summary.json', summary)
        else:
            (args.output / 'progress.json').write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + '\n')
        print(canonical(summary))
        return 0 if complete else 75


if __name__ == '__main__':
    sys.exit(main())
