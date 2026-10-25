#!/usr/bin/env python3
"""Paired response-study runner: daily burn-in, copied branch, passive lineage.

Completed pair files are immutable; incomplete pairs restart from their seed.
No biological checkpoint is saved. Baseline-only jobs compute one trajectory.
"""
import argparse
from copy import deepcopy
import fcntl
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import signal
import sys
import time
HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]
PARENT = ABM / 'research/full_icu_calibration_2026-10-07'
_spec = importlib.util.spec_from_file_location('_response_parent_runner', PARENT / 'run.py')
parent = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(parent)
TICKS_PER_DAY = parent.TICKS_PER_DAY
FIXED_MODEL = dict(parent.FIXED_MODEL)
canonical = parent.canonical
versions = parent.versions
atomic_json = parent.atomic_json

def source_hashes():
    paths = [HERE / 'run.py', HERE / 'instrumented_model.py', *[ABM / 'src' / n for n in ('CPEmodel.py', 'agents.py', 'routes.py')], PARENT / 'run.py']
    if (HERE / 'PROTOCOL.md').exists():
        paths.append(HERE / 'PROTOCOL.md')
    return {str(p.relative_to(ABM)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def normalize_config(raw, stage=None):
    allowed = {'scenarios', 'seeds', 'warmup_days', 'observation_days', 'block_days', 'description', 'model_overrides'}
    if not isinstance(raw, dict) or set(raw) - allowed:
        raise ValueError('Unknown response config keys')
    config = dict(raw)
    config.setdefault('warmup_days', 90.0)
    config.setdefault('observation_days', 28.0)
    config.setdefault('block_days', 1.0)
    if config['block_days'] != 1:
        raise ValueError('Response runner requires daily blocks')
    config['warmup_ticks'] = parent.days_to_ticks(config['warmup_days'], 'warmup_days', True)
    config['observation_ticks'] = parent.days_to_ticks(config['observation_days'], 'observation_days')
    config['block_ticks'] = TICKS_PER_DAY
    seeds = config.get('seeds')
    if not isinstance(seeds, list) or not seeds or any((type(s) is not int or s < 0 for s in seeds)) or (len(set(seeds)) != len(seeds)):
        raise ValueError('Unique nonnegative integer seeds required')
    overrides = dict(config.get('model_overrides', {}))
    if set(overrides) - {'prob_patient_sick', 'initial_environment_contaminated'}:
        raise ValueError('Unsupported model override')
    if 'initial_environment_contaminated' in overrides and type(overrides['initial_environment_contaminated']) is not bool:
        raise ValueError('Environment initialization must be boolean')
    if 'prob_patient_sick' in overrides:
        _probability(overrides['prob_patient_sick'], 'prob_patient_sick')
    config['model_overrides'] = overrides
    ids = set()
    scenarios = []
    for original in config.get('scenarios', []):
        scenario = dict(original)
        if set(scenario) - {'id', 'p', 'h0', 'h1', 'rho', 'e', 'opportunities', 'environment_multiplier', 'initial_environment_contaminated'}:
            raise ValueError('Unsupported scenario field')
        if not isinstance(scenario.get('id'), str) or not scenario['id'] or scenario['id'] in ids:
            raise ValueError('Unique nonempty scenario IDs required')
        ids.add(scenario['id'])
        scenario.setdefault('h0', 0.9)
        scenario.setdefault('h1', scenario['h0'])
        scenario.setdefault('rho', 0.0)
        scenario.setdefault('e', 1.0)
        scenario.setdefault('opportunities', 'zone_separate')
        scenario.setdefault('environment_multiplier', 1.0)
        for name in ('p', 'h0', 'h1', 'rho', 'e', 'environment_multiplier'):
            _probability(scenario.get(name), name)
        if scenario['rho'] == 1:
            raise ValueError('rho must be below1')
        if scenario['opportunities'] not in ('separate', 'shared', 'zone_separate'):
            raise ValueError('Unknown opportunity convention')
        if 'initial_environment_contaminated' in scenario and type(scenario['initial_environment_contaminated']) is not bool:
            raise ValueError('Scenario environment initialization must be boolean')
        scenarios.append(scenario)
    if not scenarios:
        raise ValueError('Nonempty scenarios required')
    config['scenarios'] = scenarios
    return config

def _probability(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or (not math.isfinite(value)) or (not 0 <= value <= 1):
        raise ValueError(name + ' must be finite in[0,1]')

def build_jobs(config):
    return [dict(job_id=i, scenario=scenario, seed=seed) for i, (scenario, seed) in enumerate(((s, seed) for s in config['scenarios'] for seed in config['seeds']))]

def snapshot(model):
    values = parent.snapshot(model)
    values['environment_contaminated_site_days'] = model.response_environment_ticks / TICKS_PER_DAY
    values['environment_contaminated_sites_endpoint'] = model.response_environment_count
    values['response_counts'] = {k: dict(v) for k, v in model.response_counts.items()}
    return values

def difference(start, end):
    values = parent.difference(start, end)
    values['environment_contaminated_site_days'] = end['environment_contaminated_site_days'] - start['environment_contaminated_site_days']
    values['mean_environment_contaminated_fraction'] = values['environment_contaminated_site_days'] / (30 * values['days']) if values['days'] else None
    values['environment_contaminated_sites_endpoint_start'] = start['environment_contaminated_sites_endpoint']
    values['environment_contaminated_sites_endpoint_end'] = end['environment_contaminated_sites_endpoint']
    values['response_counts'] = {kind: {key: count - start['response_counts'].get(kind, {}).get(key, 0) for key, count in counts.items()} for kind, counts in end['response_counts'].items()}
    if any((v < 0 for counts in values['response_counts'].values() for v in counts.values())):
        raise RuntimeError('Response cumulative counter decreased')
    return values

def simulate_pair(job, config, deadline=None, expected_hashes=None, model_factory=None):
    started = time.monotonic()
    result = dict(job, status='error')
    model = None
    try:
        if expected_hashes is not None and source_hashes() != expected_hashes:
            raise RuntimeError('Source changed before pair')
        if model_factory is None:
            spec = importlib.util.spec_from_file_location('_response_instrumented_model', HERE / 'instrumented_model.py')
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
            model_factory = module.InstrumentedModel
        scenario = job['scenario']
        kwargs = dict(FIXED_MODEL, **config.get('model_overrides', {}))
        if 'initial_environment_contaminated' in scenario:
            kwargs['initial_environment_contaminated'] = scenario['initial_environment_contaminated']
        kwargs.update(seed=job['seed'], prob_transmission=scenario['p'], icu_hcw_wash_rate=scenario['h0'], outside_hcw_wash_rate=scenario['h0'], hygiene_opportunities=scenario['opportunities'], hygiene_autocorrelation=scenario['rho'], hygiene_efficacy=scenario['e'], environment_multiplier=scenario['environment_multiplier'])
        model = model_factory(**kwargs)
        initial = snapshot(model)

        def advance(m, ticks):
            for i in range(ticks):
                if i % 128 == 0 and deadline is not None and (time.monotonic() >= deadline):
                    raise TimeoutError('Wall budget: pair incomplete; restart seed')
                m.step()

        def daily(m, ticks, start):
            blocks = []
            previous = start
            origin = start['total_ticks']
            while ticks:
                n = min(ticks, TICKS_PER_DAY)
                advance(m, n)
                current = snapshot(m)
                block = difference(previous, current)
                block['day_start'] = (previous['total_ticks'] - origin) / TICKS_PER_DAY
                block['day_end'] = (current['total_ticks'] - origin) / TICKS_PER_DAY
                blocks.append(block)
                previous = current
                ticks -= n
            return (blocks, previous)
        warmup_blocks, boundary = daily(model, config['warmup_ticks'], initial)
        shared = scenario['h0'] == scenario['h1']
        policy = None if shared else deepcopy(model)
        baseline_blocks, baseline_end = daily(model, config['observation_ticks'], boundary)
        model.assert_response_invariants()
        baseline = dict(observation=difference(boundary, baseline_end), final_snapshot=baseline_end, observation_blocks=baseline_blocks, h=scenario['h0'])
        if shared:
            policy_result = {'shared_with': 'baseline'}
        else:
            parent.set_completion(policy, scenario['h1'])
            policy_blocks, policy_end = daily(policy, config['observation_ticks'], boundary)
            policy.assert_response_invariants()
            policy_result = dict(observation=difference(boundary, policy_end), final_snapshot=policy_end, observation_blocks=policy_blocks, h=scenario['h1'])
        if expected_hashes is not None and source_hashes() != expected_hashes:
            raise RuntimeError('Source changed during pair')
        result.update(status='complete', model_kwargs=kwargs, initial_snapshot=initial, warmup_snapshot=boundary, warmup=difference(initial, boundary), warmup_blocks=warmup_blocks, arms={'baseline': baseline, 'policy': policy_result}, distinct_observation_trajectories=1 if shared else 2, executed_ticks_per_trajectory=config['warmup_ticks'] + config['observation_ticks'], warmup_executions=1)
    except TimeoutError as error:
        result.update(status='timeout', error=str(error))
    except Exception as error:
        import traceback
        result.update(error=str(error), traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic() - started
    return result

def build_manifest(raw, config):
    manifest = dict(schema_version=1, study='full_icu_response_2026-10-08', input_config=raw, resolved_config=config, jobs=build_jobs(config), fixed_model=FIXED_MODEL, source_sha256=source_hashes(), versions=versions(), checkpoint_note='Completed pairs resume; incomplete pairs restart seeded warmup; no biological checkpoint.', lineage_note='Proximate last-pickup attribution and opportunity exposure, not causal decomposition.', environment_note='Tick-start integrated burden; daily endpoints follow synchronized cleaning.')
    manifest['identity_sha256'] = hashlib.sha256(canonical(manifest).encode()).hexdigest()
    return manifest

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--stage', default='evaluate')
    parser.add_argument('--wall-seconds', type=float, default=3300.0)
    args = parser.parse_args(argv)
    raw = json.loads(args.config.read_text())
    config = normalize_config(raw, args.stage)
    manifest = build_manifest(raw, config)
    if not math.isfinite(args.wall_seconds) or args.wall_seconds <= 0:
        parser.error('Positive finite wall budget required')
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'completed').mkdir(exist_ok=True)
    stop = {'requested': False}

    def handle(number, frame):
        stop['requested'] = True
    for name in ('SIGUSR1', 'SIGTERM'):
        signal.signal(getattr(signal, name), handle)
    with (args.output / '.run.lock').open('a+') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        path = args.output / 'manifest.json'
        if path.exists():
            if json.loads(path.read_text()) != manifest:
                raise ValueError('Resume manifest differs; reject mixed source/config/runtime')
        else:
            atomic_json(path, manifest)
        deadline = time.monotonic() + args.wall_seconds
        done = 0
        for job in manifest['jobs']:
            path = args.output / 'completed' / f'{job['job_id']:06d}.json'
            if path.exists():
                row = json.loads(path.read_text())
                if any((row.get(k) != job[k] for k in ('job_id', 'seed', 'scenario'))) or row.get('status') != 'complete' or row.get('identity_sha256') != manifest['identity_sha256']:
                    raise ValueError('Corrupt or mismatched saved pair')
                done += 1
                continue
            if stop['requested'] or time.monotonic() >= deadline:
                break
            row = simulate_pair(job, config, deadline, manifest['source_sha256'])
            row['identity_sha256'] = manifest['identity_sha256']
            if row['status'] != 'complete':
                atomic_json(args.output / 'last_incomplete.json', row)
                break
            atomic_json(path, row)
            done += 1
        complete = done == len(manifest['jobs'])
        atomic_json(args.output / 'progress.json', dict(status='complete' if complete else 'partial', complete=done, planned=len(manifest['jobs'])))
        return 0 if complete else 75
if __name__ == '__main__':
    sys.exit(main())
