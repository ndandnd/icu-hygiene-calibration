#!/usr/bin/env python3
"""Apply the PROTOCOL.md budget rule to the 2-seed benchmark; write ../BENCHMARK.json.

Rule (PROTOCOL.md "Budget and stopping rule", fixed before outcomes):
  Benchmark 2 seeds first. Cap: 60 local CPU-hours, at most 10 concurrent
  processes, run with PYTHONDONTWRITEBYTECODE=1. If the projection is above
  the cap, apply this fallback, fixed now: the fine grid on bank A only (128
  seeds), everything else on all 256.

Projection (same convention as the workflow-coefficient study): per-seed CPU
seconds of the full surface (kind 'all', 364 coefficient arms) x 256 / 3600,
using the larger of the mean-based and the max-based projection. The per-seed
CPU is the whole worker process (``process_cpu_seconds_total``: imports,
frozen-hash checks, simulation, gate and record writing), which is at least
the simulation-only ``cpu_seconds``. The decision is computed here, not typed.

The fallback cost is reported only as a labelled estimate (it does not enter
the decision): a linear arm model whose fixed part is taken from the
workflow-coefficient study's measurements on this machine.

Usage: python -B code/benchmark.py [--records gate_checks/benchmark/records]
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
from pathlib import Path  # noqa: E402
import statistics  # noqa: E402
import time  # noqa: E402

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from adapter import wb  # noqa: E402
import pairs  # noqa: E402
import run_batch  # noqa: E402

STUDY = CODE.parent
CAP_CPU_HOURS = 60.0
MAX_PROCESSES = 10
N_SEEDS = 256
N_BENCHMARK = 2
RULE = ('PROTOCOL.md: Benchmark 2 seeds first. Cap: 60 local CPU-hours, at most 10 concurrent processes, run with '
        'PYTHONDONTWRITEBYTECODE=1. If the projection is above the cap, apply this fallback, fixed now: the fine grid on '
        'bank A only (128 seeds), everything else on all 256.')
# Workflow-coefficient study (same machine, loaded): 6 arms 154.2-156.1 CPU-s (6 Gate 1 runs), 44 arms 218.6-224.3.
PRIOR_ARM_MODEL = dict(source='research/workflow_coefficient_2026-10-21/gate_checks/gate1/*.json',
                       six_arm_mean_cpu_seconds=154.65, forty_four_arm_mean_cpu_seconds=221.4)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--records', type=Path, default=STUDY / 'gate_checks/benchmark/records')
    parser.add_argument('--output', type=Path, default=run_batch.BENCHMARK)
    args = parser.parse_args(argv)
    records = args.records.resolve()
    manifest = json.loads((records / 'manifest.json').read_text())
    current = run_batch.build_manifest()
    if manifest['identity_sha256'] != current['identity_sha256']:
        raise SystemExit('Benchmark records were produced under a different identity than the current code')
    rows = []
    for path in sorted(records.glob('[0-9]*.json')):
        row = run_batch.validate_record(path, manifest, int(path.stem), 'all')
        if row['surface'] != 'all':
            raise SystemExit('Benchmark records must carry the full surface (kind all)')
        site = 'local' if list(wb.ABM.glob(f'results/full_icu_confirmation_2026-10-08/bridge_local/part_*/completed/{row["seed"]}.json')) else 'cluster'
        rows.append(dict(seed=row['seed'], bank=row['bank'], archived_original_site=site,
                         process_cpu_seconds_total=row['process_cpu_seconds_total'], cpu_seconds=row['cpu_seconds'],
                         wall_seconds=row['wall_seconds'], frozen_simulate_elapsed_seconds=row['elapsed_seconds'],
                         load_average_start=row['load_average_start'], load_average_end=row['load_average_end'],
                         started_at_utc=row['started_at_utc'], completed_at_utc=row['completed_at_utc'],
                         gate_passed=row['gate']['passed'], gate_counts=row['gate']['counts'],
                         record_sha256=wb.sha256(path), coefficient_arms=2 * len(row['pairs'])))
    if len(rows) != N_BENCHMARK:
        raise SystemExit(f'Expected {N_BENCHMARK} benchmark records, found {len(rows)}')
    cpu = [r['process_cpu_seconds_total'] for r in rows]
    wall = [r['wall_seconds'] for r in rows]
    projected_mean = statistics.fmean(cpu) * N_SEEDS / 3600
    projected_max = max(cpu) * N_SEEDS / 3600
    projected = max(projected_mean, projected_max)
    surface = 'full' if projected <= CAP_CPU_HOURS else 'fallback'
    # Labelled estimate only: C(arms) = a + b * arms, b and a/b from the prior study, scaled to this benchmark.
    b_prior = (PRIOR_ARM_MODEL['forty_four_arm_mean_cpu_seconds'] - PRIOR_ARM_MODEL['six_arm_mean_cpu_seconds']) / 38
    a_over_b = (PRIOR_ARM_MODEL['six_arm_mean_cpu_seconds'] - 6 * b_prior) / b_prior
    arms_all, arms_nofine = 2 * len(pairs.KINDS['all']), 2 * len(pairs.KINDS['no_fine_grid'])
    nofine_fraction = (a_over_b + arms_nofine) / (a_over_b + arms_all)
    fallback_estimate = max(cpu) * (128 + 128 * nofine_fraction) / 3600
    result = dict(
        schema_version=1, study=run_batch.STUDY_NAME, rule=RULE,
        cap_cpu_hours=CAP_CPU_HOURS, max_concurrent_processes=MAX_PROCESSES,
        benchmark=dict(mode=f'original workflow, kind all: {len(pairs.KINDS["all"])} pairs = {arms_all} coefficient arms '
                            '(6 frozen profile arms + 358 extra arms)',
                       seeds=rows, concurrent_benchmark_processes=N_BENCHMARK,
                       seed_choice='first local-site and first cluster-site seed (archived original records) in plan order',
                       note=('Shared 12-core (6 performance + 6 efficiency) machine with jobs from other sessions; see load '
                             'averages. CPU time includes efficiency-core time.')),
        projection=dict(per_seed_cpu_measure='process_cpu_seconds_total (whole worker process)',
                        mean_cpu_seconds_per_seed=statistics.fmean(cpu), max_cpu_seconds_per_seed=max(cpu),
                        projected_cpu_hours_full_plan_256_mean_based=projected_mean,
                        projected_cpu_hours_full_plan_256_max_based=projected_max,
                        projected_cpu_hours_used_for_decision=projected,
                        mean_wall_seconds_per_seed_under_benchmark_load=statistics.fmean(wall),
                        projected_wall_hours_at_10_workers_if_wall_per_seed_unchanged=statistics.fmean(wall) * N_SEEDS / MAX_PROCESSES / 3600,
                        fallback_plan_cpu_hours_estimate=dict(
                            value=fallback_estimate, label='ESTIMATE ONLY, not used for the decision',
                            method=(f'bank A {128} x max per-seed CPU (kind all) + bank B 128 x the same x '
                                    f'{nofine_fraction:.4f}, the ratio (a/b + {arms_nofine})/(a/b + {arms_all}) of a linear arm '
                                    f'model with a/b = {a_over_b:.1f} arm-equivalents from the prior study'),
                            prior_arm_model=PRIOR_ARM_MODEL)),
        decision=dict(surface=surface,
                      meaning=run_batch.PLAN_RULES[surface],
                      rationale=f'projected {projected:.2f} CPU-hours {"<=" if surface == "full" else ">"} cap {CAP_CPU_HOURS:g}'),
        fallback_plan_for_record=run_batch.PLAN_RULES['fallback'],
        sources=dict(benchmark_py=wb.sha256(Path(__file__)), code=run_batch.code_hashes(), adapter=wb.adapter_sha256(),
                     records_manifest_identity=manifest['identity_sha256']),
        written_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    temporary = args.output.with_name('.' + args.output.name + '.tmp')
    temporary.write_text(json.dumps(result, indent=1, sort_keys=True, allow_nan=False) + '\n')
    os.replace(temporary, args.output)
    print(json.dumps(dict(decision=surface, projected_cpu_hours=round(projected, 2), cap=CAP_CPU_HOURS,
                          max_cpu_seconds=round(max(cpu), 1), mean_cpu_seconds=round(statistics.fmean(cpu), 1),
                          fallback_estimate_cpu_hours=round(fallback_estimate, 1), output=str(args.output))))
    return 0


if __name__ == '__main__':
    sys.exit(main())
