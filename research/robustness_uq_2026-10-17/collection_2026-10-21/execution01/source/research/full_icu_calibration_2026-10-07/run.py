#!/usr/bin/env python3
"""Portable synthetic full-ICU batch runner; no private data or auto-fitting.

Only the parent process writes output. Each job has a fresh seeded model.
Warmup always uses h0; h1 is applied at the observation boundary without reset.
"""
import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import hashlib
import importlib.metadata
import json
import math
import multiprocessing
from pathlib import Path
import platform
import statistics
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]
SRC = ABM / 'src'
TICKS_PER_DAY = 2592
COUNTERS = ('acquired_colonizations', 'imported_colonizations',
            'cumul_patients', 'discharged_patients', 'bed_visits',
            'patient_contacts', 'environment_contacts', 'patient_days',
            'susceptible_patient_days', 'hygiene_opportunity_count',
            'hygiene_completion_count', 'hygiene_effective_clear_count',
            'hygiene_decontamination_count')
FIXED_MODEL = dict(num_HCWs=10, num_Patients=30, num_Goo=30,
                   prob_patient_sick=.05, prob_new_patient=.003,
                   isolation_factor=.6, cleaningDay=1., isolate_sick=True,
                   height=11, width=32, stay_extension_days=0,
                   los_mode='synthetic', susceptible_los_mean_days=5.,
                   colonized_los_mean_days=14.,
                   initial_environment_contaminated=True,
                   data_collection_interval=TICKS_PER_DAY)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def source_hashes():
    paths = [HERE / 'run.py', SRC / 'CPEmodel.py', SRC / 'agents.py', SRC / 'routes.py']
    if (HERE / 'protocol.md').exists():
        paths.append(HERE / 'protocol.md')
    return {str(path.relative_to(ABM)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths}


def versions():
    result = {'python': sys.version, 'platform': platform.platform()}
    for package in ('Mesa', 'numpy', 'pandas', 'networkx'):
        try:
            result[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            result[package] = None
    return result


def days_to_ticks(days, name, allow_zero=False):
    if isinstance(days, bool) or not isinstance(days, (int, float)):
        raise ValueError(f'{name} must be numeric')
    if not math.isfinite(days) or days < 0 or (days == 0 and not allow_zero):
        raise ValueError(f'{name} must be finite and positive')
    ticks = round(days * TICKS_PER_DAY)
    if abs(ticks - days * TICKS_PER_DAY) > 1e-7:
        raise ValueError(f'{name} must be an exact number of model ticks')
    return ticks


def normalize_config(raw, stage):
    allowed = {'scenarios', 'seeds', 'warmup_days', 'observation_days',
               'block_days', 'stages', 'description'}
    if set(raw) - allowed:
        raise ValueError(f'Unknown config keys: {sorted(set(raw)-allowed)}')
    config = {k: v for k, v in raw.items() if k != 'stages'}
    if 'stages' in raw:
        stage_config = raw['stages'][stage]
        if set(stage_config) - (allowed - {'stages'}):
            raise ValueError('Unknown stage config keys')
        config.update(stage_config)
    defaults = (1., 1.) if stage == 'pilot' else (7., 14.)
    config.setdefault('warmup_days', defaults[0])
    config.setdefault('observation_days', defaults[1])
    config.setdefault('block_days', 1.)
    config['warmup_ticks'] = days_to_ticks(config['warmup_days'], 'warmup_days', True)
    config['observation_ticks'] = days_to_ticks(config['observation_days'], 'observation_days')
    config['block_ticks'] = days_to_ticks(config['block_days'], 'block_days')
    seeds = config.get('seeds')
    if not isinstance(seeds, list) or not seeds:
        raise ValueError('Provide a nonempty explicit seeds list')
    if any(isinstance(s, bool) or not isinstance(s, int) or s < 0 for s in seeds):
        raise ValueError('Seeds must be nonnegative integers')
    if len(set(seeds)) != len(seeds):
        raise ValueError('Duplicate seeds are not independent replicates')
    scenarios = config.get('scenarios')
    if not isinstance(scenarios, list) or not scenarios:
        raise ValueError('Provide a nonempty scenarios list')
    ids = set()
    normalized = []
    for original in scenarios:
        scenario = dict(original)
        if set(scenario) - {'id', 'p', 'h0', 'h1', 'rho', 'e', 'opportunities'}:
            raise ValueError(f'Unknown scenario keys: {scenario}')
        if not isinstance(scenario.get('id'), str) or not scenario['id']:
            raise ValueError('Every scenario needs a nonempty string id')
        if scenario['id'] in ids:
            raise ValueError('Scenario ids must be unique')
        ids.add(scenario['id'])
        scenario.setdefault('h0', .90)
        scenario.setdefault('h1', scenario['h0'])
        scenario.setdefault('rho', 0.)
        scenario.setdefault('e', 1.)
        scenario.setdefault('opportunities', 'separate')
        for field in ('p', 'h0', 'h1', 'rho', 'e'):
            value = scenario.get(field)
            if (isinstance(value, bool) or not isinstance(value, (float, int))
                    or not math.isfinite(value) or not 0 <= value <= 1):
                raise ValueError(f'{field} must be finite and in [0,1]')
        if scenario['rho'] == 1:
            raise ValueError('rho must be smaller than 1')
        if scenario['opportunities'] not in ('separate', 'shared', 'zone_separate'):
            raise ValueError('opportunities must be separate, shared, or zone_separate')
        normalized.append(scenario)
    config['scenarios'] = normalized
    return config


def snapshot(model):
    values = {key: getattr(model, key) for key in COUNTERS}
    values['hygiene_counts_by_role'] = {
        role: dict(counts) for role, counts in model.hygiene_counts_by_role.items()}
    values['census'] = len(model.current_patients)
    values['colonized_census'] = sum(p.colonized for p in model.current_patients)
    values['total_ticks'] = model.total_ticks
    return values


def difference(start, end):
    result = {key: end[key] - start[key] for key in COUNTERS}
    if any(value < -1e-9 for value in result.values()):
        raise RuntimeError('A cumulative counter decreased')
    result['hygiene_counts_by_role'] = {
        role: {key: value - start['hygiene_counts_by_role'].get(role, {}).get(key, 0)
               for key, value in counts.items()}
        for role, counts in end['hygiene_counts_by_role'].items()}
    result['ticks'] = end['total_ticks'] - start['total_ticks']
    days = result['ticks'] / TICKS_PER_DAY
    result['days'] = days
    result['mean_census'] = result['patient_days'] / days if days else None
    result['mean_susceptible_census'] = result['susceptible_patient_days'] / days if days else None
    result['colonized_patient_days'] = result['patient_days'] - result['susceptible_patient_days']
    result['mean_colonized_census'] = result['colonized_patient_days'] / days if days else None
    for denominator in ('patient_days', 'susceptible_patient_days'):
        result[f'acquired_per_1000_{denominator}'] = (
            1000 * result['acquired_colonizations'] / result[denominator]
            if result[denominator] else None)
    result['census_start'] = start['census']
    result['census_end'] = end['census']
    result['colonized_census_start'] = start['colonized_census']
    result['colonized_census_end'] = end['colonized_census']
    return result


def set_completion(model, h):
    model.icu_hcw_wash_rate = h
    model.outside_hcw_wash_rate = h
    for worker in model.schedule.agents:
        if worker.isHCW:
            worker.hand_wash_rate = h


def simulate_job(job, config, deadline=None, expected_hashes=None, model_factory=None):
    started = time.monotonic()
    result = dict(job, status='error')
    model = None
    try:
        if expected_hashes is not None and source_hashes() != expected_hashes:
            raise RuntimeError('Source files changed after manifest creation')
        if model_factory is None:
            sys.path.insert(0, str(SRC))
            from CPEmodel import CPE_Model
            model_factory = CPE_Model
        scenario = job['scenario']
        kwargs = dict(FIXED_MODEL, seed=job['seed'], prob_transmission=scenario['p'],
                      icu_hcw_wash_rate=scenario['h0'],
                      outside_hcw_wash_rate=scenario['h0'],
                      hygiene_opportunities=scenario['opportunities'],
                      hygiene_autocorrelation=scenario['rho'],
                      hygiene_efficacy=scenario['e'])
        model = model_factory(**kwargs)
        if model.ticks_in_day != TICKS_PER_DAY:
            raise RuntimeError('Unexpected model tick length')
        initial = snapshot(model)

        def advance(ticks):
            for tick in range(ticks):
                if tick % 128 == 0 and deadline is not None and time.monotonic() >= deadline:
                    raise TimeoutError('Batch wall budget reached; replicate is incomplete')
                model.step()

        advance(config['warmup_ticks'])
        boundary = snapshot(model)
        set_completion(model, scenario['h1'])
        blocks = []
        previous = boundary
        remaining = config['observation_ticks']
        while remaining:
            count = min(remaining, config['block_ticks'])
            advance(count)
            current = snapshot(model)
            block = difference(previous, current)
            block['observation_day_start'] = (previous['total_ticks']-boundary['total_ticks'])/TICKS_PER_DAY
            block['observation_day_end'] = (current['total_ticks']-boundary['total_ticks'])/TICKS_PER_DAY
            blocks.append(block)
            previous = current
            remaining -= count
        model.assert_invariants()
        if expected_hashes is not None and source_hashes() != expected_hashes:
            raise RuntimeError('Source files changed during replicate execution')
        result.update(status='complete', model_kwargs=kwargs,
                      warmup=difference(initial, boundary), warmup_snapshot=boundary,
                      observation=difference(boundary, previous),
                      final_snapshot=previous, observation_blocks=blocks)
    except TimeoutError as error:
        result.update(status='timeout', error=str(error))
    except Exception as error:
        result.update(error=str(error), traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic() - started
    if model is not None:
        result['executed_ticks'] = model.total_ticks
    return result


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    temporary.replace(path)


def summarize(records):
    groups = {}
    for row in records:
        if row['status'] == 'complete':
            groups.setdefault(row['scenario']['id'], []).append(row['observation'])
    result = {}
    for scenario_id, rows in groups.items():
        counts = [row['acquired_colonizations'] for row in rows]
        total_days = sum(row['patient_days'] for row in rows)
        total_susceptible = sum(row['susceptible_patient_days'] for row in rows)
        result[scenario_id] = dict(
            n_complete=len(rows), mean_acquired=statistics.mean(counts),
            sd_acquired=statistics.stdev(counts) if len(rows)>1 else None,
            acquired_total=sum(counts), patient_days_total=total_days,
            susceptible_patient_days_total=total_susceptible,
            pooled_acquired_per_1000_patient_days=1000*sum(counts)/total_days if total_days else None,
            pooled_acquired_per_1000_susceptible_patient_days=1000*sum(counts)/total_susceptible if total_susceptible else None)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--stage', choices=('pilot', 'evaluate'), required=True)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--wall-seconds', type=float, default=3600.)
    args = parser.parse_args(argv)
    if not 1 <= args.workers <= 4:
        parser.error('workers must be between 1 and 4')
    if not math.isfinite(args.wall_seconds) or args.wall_seconds <= 0:
        parser.error('wall-seconds must be finite and positive')
    raw = json.loads(args.config.read_text())
    config = normalize_config(raw, args.stage)
    hashes = source_hashes()
    jobs = [dict(job_id=index, scenario=scenario, seed=seed)
            for index, (scenario, seed) in enumerate(
                (scenario, seed) for scenario in config['scenarios'] for seed in config['seeds'])]
    # mkdir is exclusive: a concurrent invocation or any previous output fails.
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = dict(schema_version=1, stage=args.stage, input_config=raw,
                    resolved_config=config, fixed_model=FIXED_MODEL,
                    config_sha256=hashlib.sha256(canonical(raw).encode()).hexdigest(),
                    source_sha256=hashes, versions=versions(), jobs=jobs,
                    execution={'workers':args.workers,'wall_seconds':args.wall_seconds},
                    rng_note='Shared seed labels across mechanisms do not establish eventwise common random numbers. Use fresh evaluate seeds for held-out estimates.',
                    policy_note='Initialize h0; switch every worker and both role-group fields to h1 at observation start. Preserve biological and behavioral states.',
                    time_note='Patient exposures accrue at tick start; observation excludes warmup cumulative counters. Blocks are observation-relative.')
    atomic_json(args.output/'manifest.json', manifest)
    records = []
    started = time.monotonic()
    deadline = started + args.wall_seconds
    next_job = 0
    pending = {}
    status = 'running'

    def progress():
        value = dict(status=status, planned=len(jobs), submitted=next_job,
                     finished=len(records), complete=sum(r['status']=='complete' for r in records),
                     errors=sum(r['status']=='error' for r in records),
                     timeouts=sum(r['status']=='timeout' for r in records),
                     unscheduled=len(jobs)-next_job, elapsed_seconds=time.monotonic()-started,
                     summary=summarize(records))
        atomic_json(args.output/'progress.json', value)
        print(canonical({k:v for k,v in value.items() if k != 'summary'}), flush=True)

    progress()
    with (args.output/'replicates.jsonl').open('x') as stream:
        with ProcessPoolExecutor(max_workers=args.workers,
                                 mp_context=multiprocessing.get_context('spawn')) as executor:
            while pending or next_job < len(jobs):
                while (len(pending)<args.workers and next_job<len(jobs)
                       and time.monotonic()<deadline):
                    job = jobs[next_job]
                    future = executor.submit(simulate_job, job, config, deadline, hashes)
                    pending[future] = job
                    next_job += 1
                if not pending:
                    break
                completed, _ = wait(pending, timeout=10., return_when=FIRST_COMPLETED)
                for future in completed:
                    job = pending.pop(future)
                    try:
                        row = future.result()
                    except Exception as error:
                        row = dict(job, status='error', error=str(error), traceback=traceback.format_exc())
                    stream.write(canonical(row)+'\n')
                    stream.flush()
                    records.append(row)
                progress()
    if len(records)==len(jobs) and all(r['status']=='complete' for r in records):
        status = 'complete'
    elif any(r['status']=='error' for r in records):
        status = 'error'
    else:
        status = 'budget_exhausted'
    atomic_json(args.output/'summary.json', summarize(records))
    progress()
    return 0 if status == 'complete' else 2


if __name__ == '__main__':
    raise SystemExit(main())
