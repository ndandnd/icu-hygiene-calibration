#!/usr/bin/env python3
"""Existing-record precision planning only. No sampling or model execution."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[2]
DEFAULT_BANK = WORKSPACE / 'research/robustness_uq_2026-10-17/collection_2026-10-21/execution01/banks/policy/completed'

def calculate(bank):
    grouped, members = {}, []
    for path in sorted(bank.glob('job_*.json')):
        raw = path.read_bytes()
        members.append({'name': path.name, 'sha256': hashlib.sha256(raw).hexdigest()})
        record = json.loads(raw)
        assert record['status'] == 'complete', path
        observation = record['arms']['baseline']['observation']
        grouped.setdefault(record['scenario']['id'], []).append({
            'count': observation['acquired_colonizations'],
            'days': observation['patient_days'],
            'seconds': record['elapsed_seconds'],
            'seed': record['seed'],
            'p': record['scenario']['p'],
        })
    assert len(members) == 1536
    assert set(grouped) == {'independent', 'clustered', 'imperfect'}
    # Digest of canonical JSON array containing ordered filenames and raw-byte hashes.
    canonical_members = json.dumps(members, sort_keys=True, separators=(',', ':')).encode()
    output = {
        'schema': 'existing_seed_delta_precision_v1',
        'input_bank': str(bank),
        'input_records': len(members),
        'input_member_digest_algorithm': 'sha256(canonical UTF-8 JSON ordered name/sha256 array; sort_keys=True,separators=(comma,colon))',
        'input_member_digest': hashlib.sha256(canonical_members).hexdigest(),
        'input_members': members,
        'method': {
            'rate': 'R=1000*sum(C_i)/sum(D_i)',
            'influence_sd': 'sample_sd(1000*C_i-R*D_i)/mean(D_i), ddof=1',
            'standard_error': 'influence_sd/sqrt(n)',
            'planned_n': 'ceil((1.96*influence_sd/half_width)^2)',
            'confidence_level': .95,
            'normal_critical_value': 1.96,
            'half_width_units': 'acquisitions per 1000 patient-days',
        },
        'limits': [
            'No new seeds, random draws, scientific simulations, or root searches are performed.',
            'Planning approximates ordinary pointwise 95% precision at current p; it is not guaranteed coverage or a simultaneous family.',
            'Point estimates within a target tolerance do not establish true-rate matching.',
            'Whole CI containment in target +/-0.02 requires half-width less than 0.02 and a center leaving sufficient room; a +/-0.08 CI cannot certify it.',
            'Variances and runtimes can change when p moves toward a rate of 3.',
            'Cost uses recorded elapsed wall seconds for paired records multiplied by assumed one allocated CPU; exact process CPU time was not recorded.',
            'Costs are approximate one-CPU equivalents, not measured CPU usage, requested allocation, or wall-clock forecasts.',
            'Costs exclude root search, calibration, new independent validation beyond the indicated bank, retries, and multiplicity adjustments.',
        ],
        'profiles': {},
    }
    for key, rows in sorted(grouped.items()):
        n = len(rows)
        assert n == 512
        assert len({r['seed'] for r in rows}) == n
        assert len({r['p'] for r in rows}) == 1
        count = sum(r['count'] for r in rows)
        days = sum(r['days'] for r in rows)
        rate = 1000*count/days
        mean_days = days/n
        influence_sd = statistics.stdev(1000*r['count'] - rate*r['days'] for r in rows)/mean_days
        mean_seconds = statistics.mean(r['seconds'] for r in rows)
        planned = {}
        for half_width in (.08, .02):
            target_n = math.ceil((1.96*influence_sd/half_width)**2)
            planned[str(half_width)] = {
                'n_total_fresh_bank_equivalent': target_n,
                'approx_one_cpu_equivalent_hours': target_n*mean_seconds/3600,
                'additional_seeds_if_current_bank_reusable': max(0,target_n-n),
                'additional_approx_one_cpu_equivalent_hours_if_reusable': max(0,target_n-n)*mean_seconds/3600,
            }
        output['profiles'][key] = {
            'n': n, 'p': rows[0]['p'], 'pooled_count': count,
            'pooled_patient_days': days, 'rate_per_1000_patient_days': rate,
            'mean_patient_days': mean_days, 'influence_sd': influence_sd,
            'delta_standard_error': influence_sd/math.sqrt(n),
            'recorded_mean_paired_elapsed_wall_seconds': mean_seconds,
            'recorded_sum_paired_elapsed_wall_seconds': sum(r['seconds'] for r in rows),
            'assumed_allocated_cpus_per_record': 1,
            'planned_half_widths': planned,
        }
    output['totals_by_half_width'] = {
        str(h): {
            'n_total_fresh_bank_equivalent': sum(p['planned_half_widths'][str(h)]['n_total_fresh_bank_equivalent'] for p in output['profiles'].values()),
            'approx_one_cpu_equivalent_hours': sum(p['planned_half_widths'][str(h)]['approx_one_cpu_equivalent_hours'] for p in output['profiles'].values()),
        } for h in (.08,.02)
    }
    return output

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bank', type=Path, default=DEFAULT_BANK)
    parser.add_argument('--output', type=Path, default=HERE/'PRECISION_PLANNING.json')
    args = parser.parse_args()
    result = calculate(args.bank.resolve())
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'input_records': result['input_records'], 'input_member_digest': result['input_member_digest'], 'totals_by_half_width': result['totals_by_half_width']}, indent=2))
