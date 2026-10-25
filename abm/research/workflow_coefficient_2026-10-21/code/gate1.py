#!/usr/bin/env python3
"""Gate 1: original-workflow adapter reproduces archived p=0 bridge records.

One seed per invocation (run several in parallel). Uses the archived partition
config of that seed (3 profiles, 3% carriage override, 90+56 days, cap 2),
resolved by the frozen ``run.normalize`` and simulated by the frozen
``run.simulate_seed`` with the adapter's ``alternative=False`` class.

``--surface`` additionally attaches the 19 non-primary Table S5 pairs and
compares all 22 pairs with the archived extended (author-notes) record of the
same seed. Output: gate_checks/gate1/<seed>[_surface].json (atomic).
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402
from pathlib import Path  # noqa: E402
import platform  # noqa: E402
import time  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workflow_background as wb  # noqa: E402

STUDY = wb.STUDY
ARCHIVE = wb.ABM / 'results/full_icu_confirmation_2026-10-08'
EXACT_KEYS = ('seed', 'status', 'model_kwargs', 'event_counts', 'executed_ticks',
              'ordered_bridge_events_sha256', 'observation_patient_days',
              'observation_susceptible_patient_days', 'observation_admissions',
              'observation_imported_colonizations', 'total_admissions', 'total_imported_colonizations')


def archived_record(seed):
    paths = sorted(ARCHIVE.glob(f'bridge_*/part_*/completed/{seed}.json'))
    if len(paths) != 1:
        raise RuntimeError(f'Expected one archived record for seed {seed}, found {len(paths)}')
    return paths[0]


def partition_config(seed):
    seedbook = json.loads((wb.BRIDGE / 'seedbook.json').read_text())
    for site, partitions in seedbook['partitions'].items():
        for partition in partitions:
            if seed in partition['seeds']:
                return site, partition, wb.ABM / partition['config_path']
    raise RuntimeError('Seed not in seedbook partitions')


def relative_error(new, old):
    if new == old:
        return 0.0
    return abs(new - old) / max(abs(old), abs(new))


def compare_arms(new_arms, old_arms):
    out = {}
    for arm in ('baseline', 'policy'):
        new, old = new_arms[arm], old_arms[arm]
        out[arm] = dict(K_new=new['K'], K_archived=old['K'], K_relative_error=relative_error(new['K'], old['K']),
                        K_from_warmup_relative_error=relative_error(new['K_from_warmup_pickups'], old['K_from_warmup_pickups']),
                        arm_result_bitwise_equal=new == old)
    return out


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    with temporary.open('w') as handle:
        handle.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--surface', action='store_true')
    parser.add_argument('--output-dir', type=Path, default=STUDY / 'gate_checks/gate1')
    args = parser.parse_args(argv)

    record_path = archived_record(args.seed)
    record_sha = wb.sha256(record_path)
    archived = json.loads(record_path.read_text())
    site, partition, config_path = partition_config(args.seed)
    raw = json.loads(config_path.read_text())
    config = wb.bridge_run.normalize(raw)
    config = dict(config, seeds=[args.seed])
    surface = wb.SURFACE_PAIRS if args.surface else ()
    factory = wb.make_factory(alternative=False, surface_pairs=surface)

    load_before = os.getloadavg()
    cpu0, wall0 = time.process_time(), time.monotonic()
    row = wb.simulate_seed(args.seed, config, factory)
    cpu, wall = time.process_time() - cpu0, time.monotonic() - wall0

    exact = {key: row[key] == archived[key] for key in EXACT_KEYS}
    profiles = {pid: compare_arms(row['profiles'][pid], archived['profiles'][pid]) for pid in archived['profiles']}
    max_rel = max(arm['K_relative_error'] for p in profiles.values() for arm in p.values())
    result = dict(
        gate='gate1_original_mode_reproduction', seed=args.seed, site=site,
        partition_config=str(config_path.relative_to(wb.ABM)), bank=wb.seed_banks()[args.seed],
        archived_record=str(record_path.relative_to(wb.ABM)), archived_record_sha256=record_sha,
        archived_config_sha256=wb.sha256(config_path), surface_attached=bool(args.surface),
        exact_field_equal=exact,
        ordered_bridge_events_sha256=dict(new=row['ordered_bridge_events_sha256'],
                                          archived=archived['ordered_bridge_events_sha256']),
        profiles=profiles, max_K_relative_error=max_rel,
        profiles_bitwise_equal=row['profiles'] == archived['profiles'],
        workflow=row['workflow'],
        pass_=all(exact.values()) and max_rel <= 1e-9,
        cpu_seconds=cpu, wall_seconds=wall, frozen_simulate_elapsed_seconds=row['elapsed_seconds'],
        load_average_before=load_before, load_average_after=os.getloadavg(),
        adapter_sha256=wb.adapter_sha256(), gate_script_sha256=wb.sha256(Path(__file__)),
        frozen_sha256=wb.verify_frozen(),
        runtime=dict(python=sys.version, platform=platform.platform(), executable=sys.executable,
                     dont_write_bytecode=sys.dont_write_bytecode))
    if args.surface:
        extended_paths = sorted(wb.EXTENDED_RESULTS.glob(f'part_*/completed/{args.seed}.json'))
        if len(extended_paths) != 1:
            raise RuntimeError('Expected one archived extended record')
        extended = json.loads(extended_paths[0].read_text())
        pair_checks = {}
        for pair in wb.SURFACE_PAIRS:
            pair_checks[pair['id']] = compare_arms(row['pairs'][pair['id']], extended['pairs'][pair['id']])
            for arm in ('baseline', 'policy'):
                new, old = row['pairs'][pair['id']][arm], extended['pairs'][pair['id']][arm]
                # The extended record also carries the unweighted raw histogram
                # fields; compare all keys the two arm results share.
                shared = sorted(set(new) & set(old))
                pair_checks[pair['id']][arm]['shared_keys_bitwise_equal'] = all(new[k] == old[k] for k in shared)
                pair_checks[pair['id']][arm]['shared_keys'] = shared
        surface_max = max(arm['K_relative_error'] for p in pair_checks.values() for arm in p.values())
        result['surface'] = dict(
            extended_record=str(extended_paths[0]), extended_record_sha256=wb.sha256(extended_paths[0]),
            ordered_hash_equal=row['ordered_bridge_events_sha256'] == extended['ordered_bridge_events_sha256'],
            pairs=pair_checks, max_K_relative_error=surface_max,
            all_shared_keys_bitwise_equal=all(arm['shared_keys_bitwise_equal'] for p in pair_checks.values() for arm in p.values()))
        result['pass_'] = result['pass_'] and result['surface']['ordered_hash_equal'] and surface_max <= 1e-9
    result['pass'] = result.pop('pass_')
    if not math.isfinite(max_rel):
        result['pass'] = False
    out = args.output_dir / (f'{args.seed}_surface.json' if args.surface else f'{args.seed}.json')
    atomic_json(out, result)
    print(json.dumps(dict(seed=args.seed, surface=args.surface, passed=result['pass'], max_K_relative_error=max_rel,
                          cpu_seconds=round(cpu, 1), wall_seconds=round(wall, 1))))
    return 0 if result['pass'] else 1


if __name__ == '__main__':
    sys.exit(main())
