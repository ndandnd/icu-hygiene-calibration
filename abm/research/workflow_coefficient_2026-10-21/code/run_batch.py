#!/usr/bin/env python3
"""Resumable local batch: alternative-workflow p=0 coefficients on the 256 seeds.

One immutable JSON record per seed in ``records/<seed>.json``, written
atomically (temporary file, fsync, rename) by the worker that simulated it.
Completed, valid records are skipped; an interrupted seed restarts from
initialization (no biological checkpoint). A ``records/.run.lock`` prevents two
concurrent runners. At most 10 worker processes.

Science (fixed): archived bridge_all.json config resolved by the frozen
``run.normalize`` (3 primary profiles, h .90->.95, cap 2, 3% carriage, 90+56
days, zone_separate), simulated by the frozen ``run.simulate_seed`` with the
adapter's alternative-workflow class, plus the 19 other Table S5 pairs.

Surface plan (PROTOCOL "Budget and stopping rule"):
  full      all 22 pairs on all 256 seeds
  fallback  22 pairs on bank A (128 seeds); 3 primary pairs on bank B
  auto      the decision recorded in ../BENCHMARK.json (default)

Usage:
  python -B code/run_batch.py --workers 10            # launch / resume
  python -B code/run_batch.py --status                # progress, no simulation
  python -B code/run_batch.py --plan                  # print plan, no simulation
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import fcntl  # noqa: E402
import hashlib  # noqa: E402
import importlib.metadata  # noqa: E402
import json  # noqa: E402
import multiprocessing  # noqa: E402
import os  # noqa: E402
from pathlib import Path  # noqa: E402
import platform  # noqa: E402
import signal  # noqa: E402
import statistics  # noqa: E402
import time  # noqa: E402

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
import workflow_background as wb  # noqa: E402

STUDY = CODE.parent
RECORDS = STUDY / 'records'
BENCHMARK = STUDY / 'BENCHMARK.json'
MAX_WORKERS = 10
SURFACE_KINDS = {'full22': tuple(wb.SURFACE_PAIRS),
                 'primary3': tuple(p for p in wb.SURFACE_PAIRS if p.get('legacy_id'))}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def sha_text(text):
    return hashlib.sha256(text.encode()).hexdigest()


def runtime_identity():
    versions = {name: importlib.metadata.version(name) for name in ('Mesa', 'numpy', 'pandas', 'networkx')}
    return dict(python=sys.version, implementation=platform.python_implementation(),
                system=platform.system(), machine=platform.machine(), **versions)


def science():
    raw, config = wb.archived_scientific_config()
    wb.verify_surface_pairs()
    return dict(workflow=wb.WORKFLOW_NAME, alternative=True, resolved_config=config,
                surface_pairs=[dict(p) for p in wb.SURFACE_PAIRS],
                surface_histogram_cap=wb.SURFACE_HISTOGRAM_CAP,
                surface_kinds={k: [p['id'] for p in v] for k, v in SURFACE_KINDS.items()})


def source_hashes():
    values = dict(wb.verify_frozen())
    values['research/workflow_coefficient_2026-10-21/code/workflow_background.py'] = wb.adapter_sha256()
    values['extended_manifest(part_00/manifest.json)'] = wb.sha256(wb.EXTENDED_MANIFEST)
    return values


def build_manifest():
    sci = science()
    sources = source_hashes()
    manifest = dict(schema_version=1, study='workflow_coefficient_2026-10-21',
                    science=sci, science_sha256=sha_text(canonical(sci)),
                    source_sha256=sources, runtime=runtime_identity(),
                    note='Identity binds science, frozen sources, adapter and runtime. run_batch.py hash is recorded per seed but not part of identity.')
    manifest['identity_sha256'] = sha_text(canonical(manifest))
    return manifest


def plan(surface):
    banks = wb.seed_banks()
    seedbook = json.loads((wb.BRIDGE / 'seedbook.json').read_text())
    a, b = seedbook['seed_blocks']['A'], seedbook['seed_blocks']['B']
    order = [s for pair in zip(a, b) for s in pair]  # A0,B0,A1,B1,... both banks progress together
    if sorted(order) != sorted(seedbook['all_seeds']) or len(order) != 256:
        raise RuntimeError('Seed order does not cover the 256-seed seedbook')
    if surface == 'full':
        return [(seed, banks[seed], 'full22') for seed in order]
    if surface == 'fallback':
        return [(seed, banks[seed], 'full22' if banks[seed] == 'A' else 'primary3') for seed in order]
    raise ValueError('surface must be full or fallback')


def resolve_surface(choice):
    if choice != 'auto':
        return choice
    if not BENCHMARK.exists():
        raise SystemExit('BENCHMARK.json missing: run the benchmark or pass --surface explicitly')
    decision = json.loads(BENCHMARK.read_text())['decision']['surface']
    if decision not in ('full', 'fallback'):
        raise SystemExit('BENCHMARK.json has no valid decision')
    return decision


def atomic_json(path, value):
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w') as handle:
        handle.write(json.dumps(value, indent=1, sort_keys=True, allow_nan=False) + '\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    descriptor = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def seal(row):
    row['payload_sha256'] = sha_text(canonical({k: v for k, v in row.items() if k != 'payload_sha256'}))
    return row


def validate_record(path, manifest, seed, kind):
    """Return the record if complete and bound to this manifest, else raise."""
    row = json.loads(path.read_text())
    payload = {k: v for k, v in row.items() if k != 'payload_sha256'}
    if row.get('payload_sha256') != sha_text(canonical(payload)):
        raise ValueError(f'{path}: payload digest mismatch')
    if row.get('identity_sha256') != manifest['identity_sha256']:
        raise ValueError(f'{path}: bound to a different manifest identity (stale adapter/source/runtime?)')
    if row.get('seed') != seed or row.get('status') != 'complete' or row['workflow']['alternative'] is not True:
        raise ValueError(f'{path}: wrong seed/status/workflow')
    required = {p['id'] for p in SURFACE_KINDS[kind]}
    if not required <= set(row['pairs']):
        raise ValueError(f'{path}: missing required surface pairs')
    for legacy, pair_id in wb.PRIMARY_PAIR_IDS.items():
        if {k: row['pairs'][pair_id][k] for k in ('baseline', 'policy')} != row['profiles'][legacy]:
            raise ValueError(f'{path}: primary pair/profile disagreement')
    return row


# ---------------------------------------------------------------- workers
_FACTORIES = {}
_MANIFEST = None


def _init_worker(manifest):
    global _MANIFEST
    _MANIFEST = manifest
    signal.signal(signal.SIGINT, signal.SIG_IGN)  # parent handles interruption


def run_seed(task):
    seed, bank, kind, records = task
    started_wall, started_cpu = time.monotonic(), time.process_time()
    load_start = os.getloadavg()
    if kind not in _FACTORIES:
        _FACTORIES[kind] = wb.make_factory(alternative=True, surface_pairs=SURFACE_KINDS[kind])
    _, config = wb.archived_scientific_config()
    config = dict(config, seeds=[seed])
    row = wb.simulate_seed(seed, config, _FACTORIES[kind])
    if source_hashes() != _MANIFEST['source_sha256']:
        raise RuntimeError('Frozen source or adapter changed during seed ' + str(seed))
    row.update(schema_version=1, study='workflow_coefficient_2026-10-21', bank=bank, surface=kind,
               cpu_seconds=time.process_time() - started_cpu,
               wall_seconds=time.monotonic() - started_wall,
               load_average_start=list(load_start), load_average_end=list(os.getloadavg()),
               completed_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
               runtime=runtime_identity(), pid=os.getpid(),
               source_sha256=dict(_MANIFEST['source_sha256'],
                                  **{'research/workflow_coefficient_2026-10-21/code/run_batch.py': wb.sha256(Path(__file__))}),
               science_sha256=_MANIFEST['science_sha256'], identity_sha256=_MANIFEST['identity_sha256'])
    seal(row)
    path = Path(records) / f'{seed}.json'
    if path.exists():
        raise RuntimeError(f'{path} appeared during simulation (duplicate runner?)')
    atomic_json(path, row)
    return dict(seed=seed, bank=bank, surface=kind, cpu_seconds=row['cpu_seconds'], wall_seconds=row['wall_seconds'])


# ---------------------------------------------------------------- main
def status(records, tasks):
    done = [(seed, bank, kind) for seed, bank, kind in tasks if (records / f'{seed}.json').exists()]
    rows = [json.loads((records / f'{seed}.json').read_text()) for seed, _, _ in done]
    walls = [r['wall_seconds'] for r in rows]
    cpus = [r['cpu_seconds'] for r in rows]
    info = dict(planned=len(tasks), completed=len(done),
                completed_by_bank={b: sum(1 for _, bank, _ in done if bank == b) for b in 'AB'},
                failed=sorted(p.name for p in records.glob('failed_*.json')),
                temporaries=sorted(p.name for p in records.glob('*.tmp')),
                mean_wall_seconds=statistics.fmean(walls) if walls else None,
                mean_cpu_seconds=statistics.fmean(cpus) if cpus else None,
                cpu_hours_so_far=sum(cpus) / 3600)
    return info


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--workers', type=int, default=MAX_WORKERS)
    parser.add_argument('--surface', choices=('auto', 'full', 'fallback'), default='auto')
    parser.add_argument('--seeds', type=int, nargs='*', help='restrict to these seeds (benchmark)')
    parser.add_argument('--records', type=Path, default=RECORDS)
    parser.add_argument('--stop-after-hours', type=float, default=None,
                        help='submit no new seeds after this many hours (running seeds finish)')
    parser.add_argument('--retry-failed', action='store_true',
                        help='rerun seeds that have records/failed_<seed>.json')
    parser.add_argument('--plan', action='store_true')
    parser.add_argument('--status', action='store_true')
    args = parser.parse_args(argv)
    if not 1 <= args.workers <= MAX_WORKERS:
        parser.error(f'--workers must be in 1..{MAX_WORKERS}')
    surface = resolve_surface(args.surface)
    tasks = plan(surface)
    if args.seeds:
        wanted = set(args.seeds)
        if not wanted <= {seed for seed, _, _ in tasks}:
            parser.error('unknown seed requested')
        tasks = [t for t in tasks if t[0] in wanted]
    records = args.records
    records.mkdir(parents=True, exist_ok=True)
    if args.status:
        print(json.dumps(dict(surface=surface, **status(records, tasks)), indent=1))
        return 0
    if args.plan:
        for seed, bank, kind in tasks:
            print(seed, bank, kind)
        print(f'{len(tasks)} seeds, surface={surface}, workers={args.workers}')
        return 0

    manifest = build_manifest()
    with (records / '.run.lock').open('a+') as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit('Another run_batch.py holds records/.run.lock')
        manifest_path = records / 'manifest.json'
        if manifest_path.exists():
            existing = json.loads(manifest_path.read_text())
            if existing['identity_sha256'] != manifest['identity_sha256']:
                raise SystemExit('records/manifest.json identity differs from current science/sources/runtime; '
                                 'refusing to mix records. Move stale records aside after review.')
        else:
            atomic_json(manifest_path, manifest)
        pending = []
        skipped_failed = []
        for seed, bank, kind in tasks:
            path = records / f'{seed}.json'
            if path.exists():
                validate_record(path, manifest, seed, kind)
            elif (records / f'failed_{seed}.json').exists() and not args.retry_failed:
                skipped_failed.append(seed)
            else:
                pending.append((seed, bank, kind, str(records)))
        if skipped_failed:
            print(f'{len(skipped_failed)} seeds have failure records and are not retried '
                  f'(inspect records/failed_*.json; --retry-failed to rerun): {skipped_failed}', flush=True)
        log = (records / 'batch.log').open('a')
        log.write(f'{time.strftime("%Y-%m-%dT%H:%M:%S")} start surface={surface} pending={len(pending)} '
                  f'of {len(tasks)} workers={args.workers} identity={manifest["identity_sha256"][:12]}\n')
        log.flush()
        print(f'surface={surface}: {len(tasks) - len(pending)} complete, {len(pending)} pending, workers={args.workers}', flush=True)
        if not pending:
            return 0
        deadline = None if args.stop_after_hours is None else time.monotonic() + 3600 * args.stop_after_hours
        context = multiprocessing.get_context('spawn')
        interrupted, failures, done = False, 0, 0

        def on_term(signum, frame):
            raise KeyboardInterrupt
        previous = signal.signal(signal.SIGTERM, on_term)
        workers = min(args.workers, len(pending))
        pool = context.Pool(processes=workers, initializer=_init_worker, initargs=(manifest,))
        queue, inflight = list(pending), {}
        try:
            while queue or inflight:
                # Bounded in-flight dispatch: a seed starts only when a worker is free.
                while queue and len(inflight) < workers and (deadline is None or time.monotonic() < deadline):
                    task = queue.pop(0)
                    inflight[task[0]] = (task, pool.apply_async(run_seed, (task,)))
                if not inflight:
                    break  # deadline reached with nothing running
                finished = [seed for seed, (_, result) in inflight.items() if result.ready()]
                if not finished:
                    time.sleep(2)
                    continue
                for seed in finished:
                    task, result = inflight.pop(seed)
                    stamp = time.strftime('%Y-%m-%dT%H:%M:%S')
                    try:
                        info = result.get()
                    except Exception as error:  # reported, never silently retried
                        failures += 1
                        atomic_json(records / f'failed_{seed}.json',
                                    dict(seed=seed, bank=task[1], surface=task[2], error=repr(error),
                                         identity_sha256=manifest['identity_sha256'], at=stamp))
                        line = f'{stamp} FAILED seed={seed} {error!r}'
                    else:
                        done += 1
                        line = (f'{stamp} done seed={info["seed"]} bank={info["bank"]} surface={info["surface"]} '
                                f'cpu={info["cpu_seconds"]:.0f}s wall={info["wall_seconds"]:.0f}s '
                                f'[{done + failures}/{len(pending)} this run]')
                    log.write(line + '\n')
                    log.flush()
                    print(line, flush=True)
        except KeyboardInterrupt:
            interrupted = True
            pool.terminate()
        else:
            pool.close()
        finally:
            pool.join()
            signal.signal(signal.SIGTERM, previous)
            log.write(f'{time.strftime("%Y-%m-%dT%H:%M:%S")} stop done={done} failed={failures} '
                      f'not_started={len(queue)} interrupted={interrupted}\n')
            log.close()
        if interrupted or queue:
            return 75
        return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
