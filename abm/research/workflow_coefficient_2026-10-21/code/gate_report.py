#!/usr/bin/env python3
"""Collect Gate 1 / Gate 2 evidence into gate_checks/GATES.json (no simulation).

Gate 1: every gate_checks/gate1/<seed>[_surface].json written by gate1.py.
Gate 2: unit-test log(s) and the full-horizon evidence in gate_checks/gate2/,
plus the no-swap counters of the alternative-mode benchmark records.
Also re-lists frozen-tree bytecode files and compares them with
gate_checks/pycache_snapshot_start.txt.
"""
import sys

sys.dont_write_bytecode = True

import hashlib  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
from pathlib import Path  # noqa: E402
import time  # noqa: E402

CODE = Path(__file__).resolve().parent
STUDY = CODE.parent
ABM = STUDY.parents[1]
GATES = STUDY / 'gate_checks'
FROZEN_TREES = ('src', 'research/full_icu_confirmation_2026-10-08', 'research/full_icu_response_2026-10-08',
                'research/robustness_uq_2026-10-17', 'research/author_notes_2026-10-10')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pyc_listing():
    lines = []
    for tree in FROZEN_TREES:
        for path in (ABM / tree).rglob('*.pyc'):
            stat = path.stat()
            lines.append(f'{int(stat.st_mtime)} {stat.st_size} {path.relative_to(ABM)}')
    return sorted(lines, key=lambda line: line.split(' ', 2)[2])


def main():
    gate1 = []
    for path in sorted((GATES / 'gate1').glob('*.json')):
        g = json.loads(path.read_text())
        profiles = {pid: {arm: dict(K_new=v['K_new'], K_archived=v['K_archived'], rel=v['K_relative_error'])
                          for arm, v in arms.items()} for pid, arms in g['profiles'].items()}
        entry = dict(file=path.name, seed=g['seed'], site=g['site'], bank=g['bank'], passed=g['pass'],
                     ordered_hash_equal=g['ordered_bridge_events_sha256']['new'] == g['ordered_bridge_events_sha256']['archived'],
                     ordered_bridge_events_sha256=g['ordered_bridge_events_sha256']['new'],
                     all_exact_fields_equal=all(g['exact_field_equal'].values()),
                     profiles_bitwise_equal=g['profiles_bitwise_equal'], max_K_relative_error=g['max_K_relative_error'],
                     K=profiles, adapter_sha256=g['adapter_sha256'], cpu_seconds=g['cpu_seconds'],
                     original_mode_relocation_audit=g['workflow']['counters_total'])
        if 'surface' in g:
            entry['surface'] = dict(ordered_hash_equal_extended=g['surface']['ordered_hash_equal'],
                                    pairs=len(g['surface']['pairs']),
                                    max_K_relative_error=g['surface']['max_K_relative_error'],
                                    all_shared_keys_bitwise_equal=g['surface']['all_shared_keys_bitwise_equal'])
        gate1.append(entry)
    adapters = {e['adapter_sha256'] for e in gate1}
    gate2_files = sorted(p.name for p in (GATES / 'gate2').glob('*')) if (GATES / 'gate2').exists() else []
    full = [json.loads(p.read_text()) for p in sorted((GATES / 'gate2').glob('full_horizon_*.json'))]
    bench = []
    for path in sorted((GATES / 'benchmark/records').glob('[0-9]*.json')):
        row = json.loads(path.read_text())
        c = row['workflow']['counters_total']
        bench.append(dict(seed=row['seed'], swaps=c['audit.swaps'], moves_into_occupied_bed=c['audit.moves_into_occupied_bed'],
                          moves_from_isolation_bed=c['audit.moves_from_isolation_bed'],
                          relocations=c['workflow_events.isolation_relocation'],
                          deferral_patient_ticks=c['workflow_events.isolation_deferred_no_empty_bed'],
                          discharge_clears=c['environment_events.discharge_clear'],
                          relocation_clears=c['environment_events.relocation_clear']))
    start = (GATES / 'pycache_snapshot_start.txt').read_text().splitlines()
    before = sorted(start[1:], key=lambda line: line.split(' ', 2)[2])
    now = pyc_listing()
    report = dict(
        schema_version=1, written_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        adapter_sha256_now=sha256(CODE / 'workflow_background.py'),
        gate1=dict(seeds=len(gate1), all_passed=bool(gate1) and all(e['passed'] for e in gate1),
                   sites={s: sum(1 for e in gate1 if e['site'] == s) for s in ('local', 'cluster')},
                   banks={b: sum(1 for e in gate1 if e['bank'] == b) for b in 'AB'},
                   max_K_relative_error=max((e['max_K_relative_error'] for e in gate1), default=None),
                   adapter_sha256_used=sorted(adapters),
                   adapter_unchanged_since=adapters == {sha256(CODE / 'workflow_background.py')},
                   entries=gate1),
        gate2=dict(files=gate2_files, full_horizon_evidence=full, benchmark_record_no_swap=bench),
        frozen_bytecode=dict(snapshot_at_start=start[0] if start else None, files_at_start=len(before),
                             files_now=len(now), unchanged=before == now,
                             added=sorted(set(now) - set(before)), removed=sorted(set(before) - set(now))))
    out = GATES / 'GATES.json'
    temporary = out.with_name(out.name + '.tmp')
    temporary.write_text(json.dumps(report, indent=1, sort_keys=True) + '\n')
    os.replace(temporary, out)
    print(json.dumps(dict(gate1_seeds=len(gate1), gate1_all_passed=report['gate1']['all_passed'],
                          gate1_max_rel=report['gate1']['max_K_relative_error'],
                          adapter_unchanged=report['gate1']['adapter_unchanged_since'],
                          full_horizon_files=len(full), benchmark_records=len(bench),
                          frozen_bytecode_unchanged=report['frozen_bytecode']['unchanged'])))
    return 0


if __name__ == '__main__':
    sys.exit(main())
