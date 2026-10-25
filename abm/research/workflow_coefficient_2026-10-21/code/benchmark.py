#!/usr/bin/env python3
"""Apply the PROTOCOL.md budget rule to the 3-seed benchmark; write ../BENCHMARK.json.

Rule (PROTOCOL.md "Budget and stopping rule", fixed before outcomes):
  cap 80 local CPU-hours, at most 10 concurrent processes; if the full 22-row
  surface on 256 seeds is projected above the cap, run the 3 primary profiles on
  all 256 seeds and the full surface on bank A only (128 seeds).

Projection: per-seed CPU seconds (process_time of the worker that simulated the
seed, alternative mode, full 22-row surface) x 256 / 3600. The decision uses
the larger of the mean-based and max-based projections. CPU time measured on a
loaded 6P+6E machine includes efficiency-core time, so the projection is
conservative. The decision is computed here, not typed by hand.

Usage: python -B code/benchmark.py [--records gate_checks/benchmark/records]
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
from pathlib import Path  # noqa: E402
import statistics  # noqa: E402
import time  # noqa: E402

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
import run_batch  # noqa: E402
import workflow_background as wb  # noqa: E402

STUDY = CODE.parent
CAP_CPU_HOURS = 80.0
MAX_PROCESSES = 10
N_SEEDS = 256


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--records', type=Path, default=STUDY / 'gate_checks/benchmark/records')
    parser.add_argument('--output', type=Path, default=run_batch.BENCHMARK)
    args = parser.parse_args(argv)
    manifest = json.loads((args.records / 'manifest.json').read_text())
    rows = []
    for path in sorted(args.records.glob('[0-9]*.json')):
        row = run_batch.validate_record(path, manifest, int(path.stem), 'full22')
        if row['surface'] != 'full22' or row['workflow']['alternative'] is not True:
            raise SystemExit('Benchmark records must be alternative-workflow, full 22-row surface')
        rows.append(dict(seed=row['seed'], bank=row['bank'], cpu_seconds=row['cpu_seconds'],
                         wall_seconds=row['wall_seconds'], frozen_simulate_elapsed_seconds=row['elapsed_seconds'],
                         load_average_start=row['load_average_start'], load_average_end=row['load_average_end'],
                         completed_at_utc=row['completed_at_utc'], record_sha256=sha256(path),
                         ordered_bridge_events_sha256=row['ordered_bridge_events_sha256']))
    if len(rows) != 3:
        raise SystemExit(f'Expected 3 benchmark records, found {len(rows)}')
    cpu = [r['cpu_seconds'] for r in rows]
    wall = [r['wall_seconds'] for r in rows]
    projected_mean = statistics.fmean(cpu) * N_SEEDS / 3600
    projected_max = max(cpu) * N_SEEDS / 3600
    projected = max(projected_mean, projected_max)
    surface = 'full' if projected <= CAP_CPU_HOURS else 'fallback'
    # Secondary evidence: original-mode Gate 1 reruns that carried the same 22 pairs.
    gate1 = []
    for path in sorted((STUDY / 'gate_checks/gate1').glob('*_surface.json')):
        g = json.loads(path.read_text())
        gate1.append(dict(seed=g['seed'], cpu_seconds=g['cpu_seconds'], wall_seconds=g['wall_seconds'],
                          mode='original workflow, 22 pairs'))
    result = dict(
        schema_version=1, study='workflow_coefficient_2026-10-21',
        rule=('PROTOCOL.md: cap 80 local CPU-hours, at most 10 concurrent processes. If the full 22-row surface '
              'on 256 seeds is projected above the cap, run the 3 primary profiles on all 256 seeds and the full '
              'surface on bank A only (128 seeds).'),
        cap_cpu_hours=CAP_CPU_HOURS, max_concurrent_processes=MAX_PROCESSES,
        benchmark=dict(mode='alternative workflow, full 22-row surface (3 primary profile arms + 19 extra pairs = 44 coefficient arms)',
                       seeds=rows, concurrent_benchmark_processes=3,
                       note=('Run concurrently with 8 Gate 1 reruns and unrelated jobs from other sessions on a '
                             '12-core (6 performance + 6 efficiency) machine; see load averages. CPU time includes '
                             'efficiency-core time, so the projection is conservative.')),
        projection=dict(mean_cpu_seconds_per_seed=statistics.fmean(cpu), max_cpu_seconds_per_seed=max(cpu),
                        projected_cpu_hours_full_surface_256_mean_based=projected_mean,
                        projected_cpu_hours_full_surface_256_max_based=projected_max,
                        projected_cpu_hours_used_for_decision=projected,
                        mean_wall_seconds_per_seed_under_benchmark_load=statistics.fmean(wall),
                        projected_wall_hours_at_10_workers_if_wall_per_seed_unchanged=statistics.fmean(wall) * N_SEEDS / MAX_PROCESSES / 3600,
                        gate1_original_mode_22_pair_runs=gate1),
        decision=dict(surface=surface,
                      meaning=('full: all 22 pairs on all 256 seeds' if surface == 'full' else
                               'fallback: 3 primary profiles on all 256 seeds; full surface on bank A (128 seeds)'),
                      rationale=f'projected {projected:.1f} CPU-hours {"<=" if surface == "full" else ">"} cap {CAP_CPU_HOURS:g}'),
        fallback_plan_for_record='3 primary profiles on all 256 seeds; full 22-row surface on bank A only (128 seeds).',
        sources=dict(benchmark_py=sha256(Path(__file__)), run_batch_py=sha256(CODE / 'run_batch.py'),
                     adapter=wb.adapter_sha256(), records_manifest_identity=manifest['identity_sha256']),
        written_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    temporary = args.output.with_name(args.output.name + '.tmp')
    temporary.write_text(json.dumps(result, indent=1, sort_keys=True) + '\n')
    os.replace(temporary, args.output)
    print(json.dumps(dict(decision=surface, projected_cpu_hours=round(projected, 2),
                          mean_cpu_seconds=round(statistics.fmean(cpu), 1), output=str(args.output))))
    return 0


if __name__ == '__main__':
    sys.exit(main())
