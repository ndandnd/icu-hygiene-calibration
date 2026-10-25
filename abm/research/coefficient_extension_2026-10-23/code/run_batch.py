#!/usr/bin/env python3
"""Resumable local batch: original-workflow p=0 coefficients for every extension pair.

Study: research/coefficient_extension_2026-10-23 (../PROTOCOL.md). Modelled on
research/workflow_coefficient_2026-10-21/code/run_batch.py (not edited).

Science (fixed): archived bridge_all.json config resolved by the frozen
``run.normalize`` (3 primary profiles, h .90->.95, cap 2, 3% carriage, 90+56
days, zone_separate), simulated by the frozen ``run.simulate_seed`` with the
verified adapter's ``make_factory(alternative=False, surface_pairs=...)``:
the 22 archived pairs first, then the new ``x_`` pairs (pairs.py).

Records: one sealed JSON per seed, ``<records>/<seed>.json``, written once by
the process that simulated it (temporary file, fsync, hard link to the final
name so an existing record is never overwritten, directory fsync). Each record
stores the per-seed gate (gate.py) against the archived extended record. A
gate failure is stored in the record (``gate.passed = false``), not raised.

Processes: each seed runs in its own freshly spawned process (at most
``--workers`` <= 10 at a time); the parent waits on pipes and process sentinels
with a timeout, so it never blocks on a lost worker:
  * a worker that dies without reporting (SIGKILL, OOM) -> failed_<seed>.json;
  * a worker past ``--seed-timeout-hours`` -> cooperative TimeoutError from the
    frozen horizon loop; a hard kill 15 minutes later if it is still alive;
  * a worker whose parent disappears exits within 5 s (orphan watchdog), so a
    killed runner cannot leave writers behind;
  * a worker that finds a valid record (before or after simulating) returns
    'already_complete' instead of failing.
Failures are reported, never silently retried (``--retry-failed`` reruns them).
SIGTERM/SIGINT: stop dispatching, terminate running workers, log, exit 75;
rerunning the same command resumes. ``<records>/.run.lock`` (flock) prevents
two runners on one records folder.

Plans (PROTOCOL "Budget and stopping rule"):
  full      every pair on all 256 seeds                         (kind 'all')
  fallback  bank A: every pair; bank B: no fine grid (50 pairs) ('all' / 'no_fine_grid')
  auto      the decision recorded in ../BENCHMARK.json (default)

Usage (from the study folder; always -B with PYTHONDONTWRITEBYTECODE=1):
  python -B code/run_batch.py --workers 10                 # launch / resume
  python -B code/run_batch.py --status                     # progress, no simulation
  python -B code/run_batch.py --plan                       # print plan, no simulation
  python -B code/run_batch.py --surface full --seeds S1 S2 --records gate_checks/benchmark/records --workers 2
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import fcntl  # noqa: E402
import importlib.metadata  # noqa: E402
import json  # noqa: E402
import multiprocessing  # noqa: E402
from multiprocessing.connection import wait as connection_wait  # noqa: E402
import os  # noqa: E402
from pathlib import Path  # noqa: E402
import platform  # noqa: E402
import signal  # noqa: E402
import statistics  # noqa: E402
import threading  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
import adapter  # noqa: E402
from adapter import wb  # noqa: E402
import gate  # noqa: E402
import pairs  # noqa: E402

STUDY = CODE.parent
STUDY_NAME = 'coefficient_extension_2026-10-23'
RECORDS = STUDY / 'records'
BENCHMARK = STUDY / 'BENCHMARK.json'
MAX_WORKERS = 10
DEFAULT_SEED_TIMEOUT_HOURS = 6.0
HARD_KILL_GRACE_SECONDS = 900
WATCHDOG_SECONDS = 5
POLL_SECONDS = 5
PLAN_RULES = {
    'full': 'every pair (22 archived + 160 new) on all 256 seeds',
    'fallback': 'fine grid on bank A only (128 seeds); everything else on all 256 seeds',
}
CODE_FILES = ('run_batch.py', 'pairs.py', 'gate.py', 'closed_forms.py', 'adapter.py')
# Fields that legitimately differ between two runs of the same seed.
NON_SCIENCE_FIELDS = ('cpu_seconds', 'process_cpu_seconds_total', 'wall_seconds', 'elapsed_seconds', 'load_average_start', 'load_average_end',
                      'completed_at_utc', 'pid', 'payload_sha256', 'runtime', 'code_sha256', 'started_at_utc')

canonical = gate.canonical
sha_text = gate.sha_text


def runtime_identity():
    versions = {name: importlib.metadata.version(name) for name in ('Mesa', 'numpy', 'pandas', 'networkx')}
    return dict(python=sys.version, implementation=platform.python_implementation(),
                system=platform.system(), machine=platform.machine(), **versions)


def science():
    _, config = wb.archived_scientific_config()
    wb.verify_surface_pairs()
    return dict(workflow='original_swap_isolation', alternative=False, resolved_config=config,
                pairs=pairs.pair_table(),
                kinds={kind: [p['id'] for p in members] for kind, members in pairs.KINDS.items()},
                transition_roots=pairs.TRANSITION_ROOTS,
                new_pair_histogram_cap=wb.SURFACE_HISTOGRAM_CAP,
                plans=PLAN_RULES)


def source_hashes():
    values = dict(wb.verify_frozen())
    values[adapter.ADAPTER_RELATIVE] = wb.adapter_sha256()
    values['extended_manifest(part_00/manifest.json)'] = wb.sha256(wb.EXTENDED_MANIFEST)
    return values


def code_hashes():
    return {'research/' + STUDY_NAME + '/code/' + name: wb.sha256(CODE / name) for name in CODE_FILES}


def build_manifest():
    sci = science()
    manifest = dict(schema_version=1, study=STUDY_NAME, science=sci, science_sha256=sha_text(canonical(sci)),
                    source_sha256=source_hashes(), runtime=runtime_identity(),
                    note=('Identity binds science (pairs, kinds, roots, config), frozen sources, the adapter, the '
                          'archived extended manifest and the runtime. Study code hashes are recorded per seed.'))
    manifest['identity_sha256'] = sha_text(canonical(manifest))
    return manifest


# ------------------------------------------------------------------- plan
def seed_order():
    seedbook = json.loads((wb.BRIDGE / 'seedbook.json').read_text())
    a, b = seedbook['seed_blocks']['A'], seedbook['seed_blocks']['B']
    order = [s for pair in zip(a, b) for s in pair]  # A0,B0,A1,B1,...
    if sorted(order) != sorted(seedbook['all_seeds']) or len(order) != 256 or len(set(order)) != 256:
        raise RuntimeError('Seed order does not cover the 256-seed seedbook')
    return order


def plan(surface):
    banks = wb.seed_banks()
    order = seed_order()
    if surface == 'full':
        return [(seed, banks[seed], 'all') for seed in order]
    if surface == 'fallback':
        return [(seed, banks[seed], 'all' if banks[seed] == 'A' else 'no_fine_grid') for seed in order]
    raise ValueError('surface must be full or fallback')


def resolve_surface(choice, benchmark=BENCHMARK):
    if choice != 'auto':
        return choice
    if not benchmark.exists():
        raise SystemExit('BENCHMARK.json missing: run the benchmark or pass --surface explicitly')
    decision = json.loads(benchmark.read_text())['decision']['surface']
    if decision not in ('full', 'fallback'):
        raise SystemExit('BENCHMARK.json has no valid decision')
    return decision


# ---------------------------------------------------------------- records
def fsync_directory(path):
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def write_temporary(path, text):
    temporary = path.with_name(f'.{path.name}.{os.getpid()}.tmp')
    with temporary.open('w') as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    return temporary


def atomic_json(path, value):
    """Atomic replace (manifest, failure notes)."""
    temporary = write_temporary(path, json.dumps(value, indent=1, sort_keys=True, allow_nan=False) + '\n')
    os.replace(temporary, path)
    fsync_directory(path.parent)


def publish_once(path, value):
    """Atomic, no-clobber publication; FileExistsError if ``path`` already exists."""
    temporary = write_temporary(path, json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')
    try:
        os.link(temporary, path)
    finally:
        temporary.unlink()
    fsync_directory(path.parent)


def seal(row):
    row['payload_sha256'] = sha_text(canonical({k: v for k, v in row.items() if k != 'payload_sha256'}))
    return row


def validate_record(path, manifest, seed, kind):
    """Return the record if sealed, bound to ``manifest`` and covering ``kind``; else raise.

    A stored gate failure does not make a record invalid (the seed is then
    excluded and reported by summarize.py).
    """
    row = json.loads(Path(path).read_text())
    payload = {k: v for k, v in row.items() if k != 'payload_sha256'}
    if row.get('payload_sha256') != sha_text(canonical(payload)):
        raise ValueError(f'{path}: payload digest mismatch')
    if row.get('identity_sha256') != manifest['identity_sha256']:
        raise ValueError(f'{path}: bound to a different manifest identity (stale adapter/source/runtime/pairs?)')
    if row.get('seed') != seed or row.get('status') != 'complete' or row.get('study') != STUDY_NAME:
        raise ValueError(f'{path}: wrong seed/status/study')
    if row['workflow']['alternative'] is not False or row['workflow']['name'] != 'original_swap_isolation':
        raise ValueError(f'{path}: not the original workflow')
    required = [p['id'] for p in pairs.KINDS[kind]]
    if row.get('surface') not in pairs.KINDS or not set(required) <= set(row['pairs']):
        raise ValueError(f'{path}: missing required pairs for kind {kind}')
    config = manifest['science']['resolved_config']
    if row.get('executed_ticks') != config['warmup_ticks'] + config['observation_ticks']:
        raise ValueError(f'{path}: wrong horizon')
    for legacy, pair_id in gate.PRIMARY_PAIR_IDS.items():
        if not gate.bitwise_equal({k: row['pairs'][pair_id][k] for k in gate.ARMS}, row['profiles'][legacy]):
            raise ValueError(f'{path}: primary pair/profile disagreement')
    if not isinstance(row.get('gate'), dict) or type(row['gate'].get('passed')) is not bool:
        raise ValueError(f'{path}: no stored gate result')
    return row


def science_view(row):
    return {k: v for k, v in row.items() if k not in NON_SCIENCE_FIELDS}


# ---------------------------------------------------------------- worker
def run_seed(task, manifest):
    """Simulate one seed, run its gate and publish its record (runs in a fresh process)."""
    seed, bank, kind, records, timeout_seconds = task
    records = Path(records)
    path = records / f'{seed}.json'
    if path.exists():
        row = validate_record(path, manifest, seed, kind)
        return dict(seed=seed, bank=bank, surface=kind, status='already_complete', gate_passed=row['gate']['passed'],
                    note='valid record found before simulating')
    started_wall, started_cpu = time.monotonic(), time.process_time()
    started_utc = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    load_start = os.getloadavg()
    factory = wb.make_factory(alternative=False, surface_pairs=pairs.KINDS[kind])
    _, config = wb.archived_scientific_config()
    config = dict(config, seeds=[seed])
    row = wb.simulate_seed(seed, config, factory, time.monotonic() + timeout_seconds)
    if source_hashes() != manifest['source_sha256']:
        raise RuntimeError('Frozen source or adapter changed during seed ' + str(seed))
    extended = gate.load_extended(seed)
    gate_result = gate.run_gate(row, extended['row'])
    gate_result['archived_extended_record'] = dict(path=extended['path'], file_sha256=extended['file_sha256'])
    row.update(schema_version=1, study=STUDY_NAME, bank=bank, surface=kind, gate=gate_result,
               cpu_seconds=time.process_time() - started_cpu, wall_seconds=time.monotonic() - started_wall,
               process_cpu_seconds_total=time.process_time(),
               load_average_start=list(load_start), load_average_end=list(os.getloadavg()),
               started_at_utc=started_utc, completed_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
               runtime=runtime_identity(), pid=os.getpid(), code_sha256=code_hashes(),
               source_sha256=manifest['source_sha256'], science_sha256=manifest['science_sha256'],
               identity_sha256=manifest['identity_sha256'])
    seal(row)
    info = dict(seed=seed, bank=bank, surface=kind, gate_passed=gate_result['passed'],
                cpu_seconds=row['cpu_seconds'], wall_seconds=row['wall_seconds'])
    try:
        publish_once(path, row)
    except FileExistsError:
        existing = validate_record(path, manifest, seed, kind)
        if canonical(science_view(existing)) != canonical(science_view(row)):
            raise RuntimeError(f'{path} appeared during simulation and differs from this run')
        return dict(info, status='already_complete', note='an identical valid record appeared during simulation')
    return dict(info, status='complete')


def _orphan_watchdog(parent_pid):
    while True:
        time.sleep(WATCHDOG_SECONDS)
        if os.getppid() != parent_pid:
            os._exit(70)


def child_main(task, manifest, connection, parent_pid):
    signal.signal(signal.SIGINT, signal.SIG_IGN)  # the runner handles interruption
    threading.Thread(target=_orphan_watchdog, args=(parent_pid,), daemon=True).start()
    try:
        message = ('ok', run_seed(task, manifest))
    except BaseException as error:  # reported to the parent, which writes failed_<seed>.json
        message = ('error', repr(error), traceback.format_exc())
    try:
        connection.send(message)
    finally:
        connection.close()


# ---------------------------------------------------------------- status
def lock_is_held(records):
    lock_path = records / '.run.lock'
    if not lock_path.exists():
        return False
    with lock_path.open('a+') as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        return False


def status(records, tasks):
    if not records.exists():
        return dict(planned=len(tasks), completed=0, records_dir_exists=False)
    manifest_path = records / 'manifest.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else None
    done, invalid, gate_failed, cpus, walls = [], [], [], [], []
    for seed, bank, kind in tasks:
        path = records / f'{seed}.json'
        if not path.exists():
            continue
        try:
            row = validate_record(path, manifest, seed, kind) if manifest else None
        except Exception as error:  # reported, not raised
            invalid.append(dict(seed=seed, error=repr(error)))
            continue
        if row is None:
            invalid.append(dict(seed=seed, error='no manifest'))
            continue
        done.append((seed, bank, kind))
        cpus.append(row['cpu_seconds'])
        walls.append(row['wall_seconds'])
        if not row['gate']['passed']:
            gate_failed.append(seed)
    pending = len(tasks) - len(done) - len(invalid)
    return dict(planned=len(tasks), completed=len(done), pending=pending,
                completed_by_bank={b: sum(1 for _, bank, _ in done if bank == b) for b in 'AB'},
                completed_by_kind={k: sum(1 for _, _, kind in done if kind == k) for k in pairs.KINDS},
                gate_failed=gate_failed, invalid_records=invalid,
                failed=sorted(p.name for p in records.glob('failed_*.json')),
                temporaries=sorted(p.name for p in records.glob('.*.tmp')),
                runner_active=lock_is_held(records),
                mean_wall_seconds=statistics.fmean(walls) if walls else None,
                mean_cpu_seconds=statistics.fmean(cpus) if cpus else None,
                cpu_hours_so_far=sum(cpus) / 3600)


# ---------------------------------------------------------------- main
class Running:
    def __init__(self, task, process, connection):
        self.task, self.process, self.connection = task, process, connection
        self.started = time.monotonic()
        self.killed_for_timeout = False


def stop_process(process, grace=10.0):
    if process.is_alive():
        process.terminate()
        process.join(grace)
    if process.is_alive():
        process.kill()
        process.join(grace)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--workers', type=int, default=MAX_WORKERS)
    parser.add_argument('--surface', choices=('auto', 'full', 'fallback'), default='auto')
    parser.add_argument('--seeds', type=int, nargs='*', help='restrict to these seeds (benchmark, determinism)')
    parser.add_argument('--records', type=Path, default=RECORDS)
    parser.add_argument('--seed-timeout-hours', type=float, default=DEFAULT_SEED_TIMEOUT_HOURS,
                        help='per-seed wall-clock limit (cooperative), hard kill 15 min later')
    parser.add_argument('--stop-after-hours', type=float, default=None,
                        help='submit no new seeds after this many hours (running seeds finish)')
    parser.add_argument('--retry-failed', action='store_true', help='rerun seeds that have failed_<seed>.json')
    parser.add_argument('--plan', action='store_true')
    parser.add_argument('--status', action='store_true')
    args = parser.parse_args(argv)
    if not 1 <= args.workers <= MAX_WORKERS:
        parser.error(f'--workers must be in 1..{MAX_WORKERS}')
    if not args.seed_timeout_hours > 0:
        parser.error('--seed-timeout-hours must be positive')
    surface = resolve_surface(args.surface)
    tasks = plan(surface)
    if args.seeds:
        wanted = set(args.seeds)
        if not wanted <= {seed for seed, _, _ in tasks}:
            parser.error('unknown seed requested')
        tasks = [t for t in tasks if t[0] in wanted]
    records = args.records.resolve()
    if args.status:
        print(json.dumps(dict(surface=surface, records=str(records), **status(records, tasks)), indent=1))
        return 0
    if args.plan:
        for seed, bank, kind in tasks:
            print(seed, bank, kind, len(pairs.KINDS[kind]))
        arms = sum(2 * len(pairs.KINDS[kind]) for _, _, kind in tasks)
        print(f'{len(tasks)} seeds, surface={surface}, workers={args.workers}, coefficient arms in all={arms}')
        return 0

    manifest = build_manifest()
    records.mkdir(parents=True, exist_ok=True)
    with (records / '.run.lock').open('a+') as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit(f'Another run_batch.py holds {records}/.run.lock')
        manifest_path = records / 'manifest.json'
        if manifest_path.exists():
            existing = json.loads(manifest_path.read_text())
            if existing['identity_sha256'] != manifest['identity_sha256']:
                raise SystemExit('manifest.json identity differs from current science/sources/runtime; '
                                 'refusing to mix records. Move stale records aside after review.')
        else:
            atomic_json(manifest_path, manifest)
        for stale in records.glob('.*.tmp'):  # left by a killed worker; no writer can be alive (lock + watchdog)
            stale.unlink()
        pending, skipped_failed = [], []
        timeout_seconds = 3600 * args.seed_timeout_hours
        for seed, bank, kind in tasks:
            path = records / f'{seed}.json'
            if path.exists():
                validate_record(path, manifest, seed, kind)
            elif (records / f'failed_{seed}.json').exists() and not args.retry_failed:
                skipped_failed.append(seed)
            else:
                pending.append((seed, bank, kind, str(records), timeout_seconds))
        if skipped_failed:
            print(f'{len(skipped_failed)} seeds have failure records and are not retried '
                  f'(inspect failed_*.json; --retry-failed to rerun): {skipped_failed}', flush=True)
        log = (records / 'batch.log').open('a')

        def emit(line):
            line = f'{time.strftime("%Y-%m-%dT%H:%M:%S")} {line}'
            log.write(line + '\n')
            log.flush()
            print(line, flush=True)

        emit(f'start surface={surface} pending={len(pending)} of {len(tasks)} workers={args.workers} '
             f'identity={manifest["identity_sha256"][:12]} timeout_h={args.seed_timeout_hours:g}')
        if not pending:
            log.close()
            return 0
        deadline = None if args.stop_after_hours is None else time.monotonic() + 3600 * args.stop_after_hours
        context = multiprocessing.get_context('spawn')
        stopping = {'signal': None}

        def on_signal(signum, frame):
            stopping['signal'] = signum
        previous = {s: signal.signal(s, on_signal) for s in (signal.SIGTERM, signal.SIGINT)}
        queue, running = list(pending), {}
        done = failures = already = gate_failures = 0
        workers = min(args.workers, len(pending))

        def fail(task, reason, detail=None):
            seed = task[0]
            atomic_json(records / f'failed_{seed}.json',
                        dict(seed=seed, bank=task[1], surface=task[2], error=reason, detail=detail,
                             identity_sha256=manifest['identity_sha256'], at=time.strftime('%Y-%m-%dT%H:%M:%S')))
            emit(f'FAILED seed={seed} {reason}')

        def supersede_failure(seed):
            old = records / f'failed_{seed}.json'
            if old.exists():
                os.replace(old, records / f'failed_{seed}.superseded_{time.strftime("%Y%m%dT%H%M%S")}.json')

        try:
            while queue or running:
                while (queue and len(running) < workers and stopping['signal'] is None
                       and (deadline is None or time.monotonic() < deadline)):
                    task = queue.pop(0)
                    receiver, sender = context.Pipe(duplex=False)
                    process = context.Process(target=child_main, args=(task, manifest, sender, os.getpid()),
                                              name=f'seed-{task[0]}')
                    process.start()
                    sender.close()
                    running[task[0]] = Running(task, process, receiver)
                if stopping['signal'] is not None or not running:
                    break
                connection_wait([r.connection for r in running.values()] + [r.process.sentinel for r in running.values()],
                                timeout=POLL_SECONDS)
                now = time.monotonic()
                for seed, entry in list(running.items()):
                    message = None
                    try:
                        if entry.connection.poll():
                            message = entry.connection.recv()
                    except (EOFError, OSError):
                        message = None
                    if message is None and entry.process.is_alive():
                        if now - entry.started > entry.task[4] + HARD_KILL_GRACE_SECONDS:
                            stop_process(entry.process)
                            entry.connection.close()
                            del running[seed]
                            failures += 1
                            fail(entry.task, f'hard timeout after {now - entry.started:.0f}s; worker killed')
                        continue
                    entry.process.join(60)
                    if entry.process.is_alive():
                        stop_process(entry.process)
                    entry.connection.close()
                    del running[seed]
                    if message is None:
                        path = records / f'{seed}.json'
                        try:
                            row = validate_record(path, manifest, seed, entry.task[2]) if path.exists() else None
                        except Exception:
                            row = None
                        if row is not None:
                            already += 1
                            supersede_failure(seed)
                            emit(f'worker for seed={seed} exited (code {entry.process.exitcode}) without reporting, '
                                 f'but its record is valid')
                        else:
                            failures += 1
                            fail(entry.task, f'worker lost: exit code {entry.process.exitcode}, no result')
                    elif message[0] == 'ok':
                        info = message[1]
                        supersede_failure(seed)
                        if info['status'] == 'already_complete':
                            already += 1
                            emit(f'already complete seed={seed} ({info.get("note")})')
                        else:
                            done += 1
                            gate_failures += not info['gate_passed']
                            emit(f'done seed={seed} bank={info["bank"]} surface={info["surface"]} '
                                 f'gate={"pass" if info["gate_passed"] else "FAIL"} cpu={info["cpu_seconds"]:.0f}s '
                                 f'wall={info["wall_seconds"]:.0f}s [{done + already + failures}/{len(pending)} this run]')
                    else:
                        failures += 1
                        fail(entry.task, message[1], message[2])
        finally:
            for seed, entry in running.items():
                stop_process(entry.process)
                entry.connection.close()
            for stale in records.glob('.*.tmp'):
                stale.unlink()
            for s, handler in previous.items():
                signal.signal(s, handler)
            emit(f'stop done={done} already_complete={already} failed={failures} gate_failed={gate_failures} '
                 f'interrupted_running={len(running)} not_started={len(queue)} signal={stopping["signal"]}')
            log.close()
        if stopping['signal'] is not None or queue or running:
            return 75
        return 1 if failures or gate_failures else 0


if __name__ == '__main__':
    sys.exit(main())
