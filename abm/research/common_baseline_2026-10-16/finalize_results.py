#!/usr/bin/env python3
"""Assemble reviewed common-baseline outputs; no simulation or bootstrap runs.

Run from any directory: python3 finalize_results.py --output OUTPUT_DIRECTORY
The two written files are RESULTS.md and results_final.json. Existing output
files are replaced only with deterministic products of the bound review inputs.
The original protocol, implementations and review evidence are never modified.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name):
    return json.loads((ROOT / name).read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    a = read('impl_a/results.json')
    b = read('impl_b/results.json')
    comparison = read('reconcile/comparison.json')['summary']
    stats = read('review_2026-10-17/statistics/review.json')
    provenance = read('review_2026-10-17/provenance_rules/AUDIT.json')
    assert stats['status'] == 'accepted_with_documented_secondary_correction_and_claim_limits'
    assert provenance['status'] == 'passed_with_runtime_limitation'
    assert comparison['n_points_outside_tol'] == comparison['n_intervals_outside_tol'] == 0
    assert stats['primary_A'] == a['primary'] and stats['contrasts_A'] == a['contrasts']
    assert stats['primary_unchanged'] and stats['original_eight_failed_gates_preserved']
    assert stats['raw_records_rehashed'] == 3328
    assert provenance['source_maps_identical_across_four_banks']
    assert provenance['paired_seed_identity_and_config_checks']
    assert provenance['route_review']['D1_original_summary_reproduced_exactly']
    for name, expected in stats['files'].items():
        assert sha(ROOT / name) == expected, f'Changed reviewed input: {name}'
    bound = list(stats['files']) + [
        'review_2026-10-17/statistics/review.json',
        'review_2026-10-17/provenance_rules/audit.py',
        'review_2026-10-17/provenance_rules/AUDIT.json',
        'review_2026-10-17/provenance_rules/REVIEW.md',
        'finalize_results.py',
    ]
    profiles = a['primary']['profiles']
    absolute = {
        name: {
            'unit': 'acquisitions per 1000 patient-days',
            'policy_rate': 3 * (1 - row['F_star']),
            'rate_prevented': 3 * row['F_star'],
            'rate_prevented_ci95': [3 * row['ci95_lo'], 3 * row['ci95_hi']],
            'observed_count_at_p_star': None,
            'interpretation': 'Derived from target 3 and interpolated F; not new observed counts.',
        } for name, row in profiles.items()
    }
    result = {
        'schema': 'reviewed_common_baseline_results_v1',
        'date': '2026-10-17',
        'status': 'accepted_post_hoc_secondary_conditional_interpolation',
        'reporting_implementation': 'A; B and third recomputation used for validation only',
        'simulations_added': 0,
        'cluster_operations': 0,
        'previous_matching_checks': {'failed': 8, 'changed': 0},
        'primary': a['primary'],
        'contrasts': a['contrasts'],
        'derived_absolute_rates': absolute,
        'descriptive_sensitivity': {
            'S1': a['S1'], 'S2': a['S2'],
            'S3_corrected_overlap': stats['S3_corrected_descriptive'],
            'S3_original_intervals': {
                'status': 'retained_in_original_results_but_not_accepted_for_inference',
                'reason': 'Independent upper256 is a subset of holdout lower512; old rule ignored overlap.',
                'amendment': 'review_2026-10-17/statistics/AMENDMENT_S3_2026-10-17.md',
            },
            'S4': a['S4'],
            'Linux_only_post_review': stats['linux_subset_post_review_descriptive'],
        },
        'D1_passive_routes': provenance['route_review'],
        'input_review': {
            'records_verified': 3328,
            'original_implementations_comparison': comparison,
            'third_max_point_difference': stats['third_primary_max_point_error_vs_A'],
            'third_max_interval_difference_pp': stats['third_primary_max_interval_endpoint_error_pp_vs_A'],
            'third_max_contrast_interval_difference_pp': stats['third_contrast_max_interval_endpoint_error_pp_vs_A'],
            'third_was_blind': False,
            'source_maps_identical': True,
            'input_runtime_limitation': provenance['limitations'][:2],
            'reporting_analysis_runtime': a['software'],
        },
        'allowed_claims': stats['allowed_claims'],
        'not_allowed_claims': stats['not_allowed_claims'],
        'files_sha256': {name: sha(ROOT / name) for name in sorted(set(bound))},
    }
    names = {'independent': 'Independent misses', 'clustered': 'Clustered misses', 'imperfect': 'Imperfect clearing'}
    lines = [
        '# Reviewed common-baseline results — 17 October 2026', '',
        '**Post hoc secondary analysis of existing simulations.** Under log-rate',
        'interpolation to a baseline target of 3.0 acquisitions per 1,000 patient-days,',
        'the three profiles retain different estimated responses to increasing hygiene',
        'completion from 90% to 95%. All eight earlier matching checks remain failed.',
        'The target is imposed by the interpolation, not established by fresh simulation.', '',
        '## Primary results', '',
        '| Profile | Baseline bracket | Implied p* | Prevention [95% MC interval] | Absolute rate prevented [95% MC interval] |',
        '|---|---:|---:|---:|---:|',
    ]
    for name, r in profiles.items():
        v = absolute[name]
        lines.append(f"| {names[name]} | {r['B_lo']:.3f}–{r['B_hi']:.3f} | {r['p_star']:.6f} | {100*r['F_star']:.1f}% [{100*r['ci95_lo']:.1f}, {100*r['ci95_hi']:.1f}] | {v['rate_prevented']:.2f} [{v['rate_prevented_ci95'][0]:.2f}, {v['rate_prevented_ci95'][1]:.2f}] |")
    lines += ['', 'Rates are per 1,000 patient-days. The absolute decrease is 3F*, not an',
              'observed count difference at p*. No runs were added at these inferred p values.',
              'The original 512-pair care-zone count benefits remain 3.13, 2.46 and 1.24',
              'mean carriers prevented per 56-day ICU run at their original probabilities.', '',
              '| Contrast | Difference in prevention (percentage points) | 98.333% MC interval |',
              '|---|---:|---:|']
    for c in a['contrasts']['pairs']:
        lines.append(f"| {c['pair']} | {100*c['estimate']:.1f} | [{100*c['ci_lo']:.1f}, {100*c['ci_hi']:.1f}] |")
    lines += ['', 'All three intervals exclude zero in this secondary family of three.',
              'They are Bonferroni-adjusted percentile bootstrap intervals, not clinical',
              'confidence bounds or a test of actual baseline equivalence.', '',
              '## Estimator and uncertainty', '',
              'For each profile, use 256 shared seed pairs at the center and 1.25-times-center',
              'probabilities. Compute each control/policy rate as a ratio of pooled counts to',
              'patient-days. Interpolate log control and log policy rates linearly in log p:', '',
              'w = (log(3) − log(B_lower)) / (log(B_upper) − log(B_lower));',
              'p* = p_lower^(1−w) p_upper^w;',
              'F* = 1 − P_lower^(1−w) P_upper^w / 3.', '',
              'Resample whole seed vectors 50,000 times within A/B banks of 128. Shared',
              'conditions retain their pairing; the independent and sensitivity modules use',
              'disjoint seed lists, and clustered/imperfect conditions share sensitivity draws.',
              'The bootstrap re-estimates w and p* but conditions on the simulated endpoint',
              'probabilities and the interpolation rule. It excludes interpolation error,',
              'original grid/parameter-selection uncertainty, structural uncertainty and clinical',
              'input uncertainty. See the preserved [protocol](PROTOCOL.md).', '',
              '| Profile | Extrapolating draws / 50,000 | Undefined draws |',
              '|---|---:|---:|']
    for name, r in profiles.items():
        lines.append(f"| {names[name]} | {r['n_extrapolated']:,} ({100*r['frac_extrapolated']:.3f}%) | {r['undefined_draws']} |")
    lines += ['', 'No draw was dropped as undefined. Extrapolating draws are retained.', '',
              '## Review findings and corrections', '',
              '- Two independently written implementations agree within the protocol tolerances;',
              '  all 140 reconciled point comparisons pass and the largest interval endpoint',
              '  difference across the reconciled analyses is 0.181 percentage points.',
              '- A third implementation independently reread the raw data and ran 50,000 draws.',
              '  It was written after inspecting A/B, so it is not described as blind.',
              '  Primary point agreement is within 5.56e-16 and primary interval endpoint',
              '  agreement within 0.058 percentage points (contrasts within 0.104).',
              '- Independently checked 3,328 archived records, identical relevant simulator',
              '  source maps, matching scenario settings apart from declared parameters, seed',
              '  identities, and all 3,072 original holdout arms used for route counters.',
              '- S3 wrongly called the independent lower512 and upper256 samples unmatched.',
              '  The upper seeds are a subset of the lower. The original files remain intact;',
              '  a [dated amendment](review_2026-10-17/statistics/AMENDMENT_S3_2026-10-17.md)',
              '  and an overlap-aware four-stratum bootstrap correct this descriptive check.',
              '  The primary analysis is unchanged.',
              '- Report implementation A consistently. Its rounded independent lower95 bound',
              '  is 62.1%; the two independent contrast upper bounds are 19.2 and 46.2 points.',
              '  The handoff mixed rounded endpoints from the two bootstrap realizations.', '',
              '## Descriptive robustness checks', '',
              '- S1 arithmetic interpolation shifts F* by at most 0.03 percentage points.',
              '- S2 three-point curvature shifts clustered and imperfect F* by −0.21 and',
              '  −0.07 points. Independent lacks a third matched condition for this check.',
              '- Corrected S3 estimates: independent67.8% [64.8,70.8], clustered55.5%',
              '  [52.3,58.7], imperfect27.5% [23.0,32.0]. Point estimates are unchanged;',
              '  intervals now preserve the actual overlap. These are descriptive95% intervals.',
              '- S4 at target2.7 extrapolates in72.0% of independent and84.1% of imperfect',
              '  draws; it does not validate interpolation across the whole target band.',
              '- A review-requested Linux-only subset of192 independent pairs (96/bank)',
              '  gives66.64% [61.33,72.04], a −0.055-point change, with18.846% extrapolating',
              '  draws and none undefined. This changes seeds as well as platform mix and',
              '  cannot establish platform equivalence.', '',
              'These checks do not provide a bound on interpolation or model error.', '',
              '## Runtime limitation', '',
              'The primary independent lower bank contains64 macOS/Python3.12.2 and192',
              'Linux/Python3.12.3 records; all upper and sensitivity records use Linux.',
              'Scientific source hashes and recorded package versions agree, but exact',
              'cross-platform replay is not verified. Platform effects could affect marginal',
              'estimates, not merely pairing efficiency. The intervals condition on this',
              'recorded execution mixture; the earlier handoff assertion is withdrawn.', '',
              '## Verified passive pathway counts', '',
              '| Profile | Arm | Site pickup, zero opportunities / all acquisitions | Share |',
              '|---|---|---:|---:|']
    for row in provenance['route_review']['summary']:
        lines.append(f"| {row['profile']} | {row['arm']} | {row['env0']:,}/{row['total']:,} | {100*row['site_zero_opportunity_fraction']:.1f}% |")
    lines += ['', 'These are labels of realized uncleared hand-contamination episodes. They',
              'are not causal shares or a fixed hygiene-insensitive fraction. Counts can',
              'fall while their share rises, as in the independent profile. Earlier-admission',
              'labels do not separately identify contamination inherited through isolation swaps.',
              'The source audit verifies synchronized daily cleaning only, no terminal or',
              'relocation cleaning, fixed bed sites,16 isolation/14 shared beds, occupant',
              'swaps, and the .6 multiplier on isolation patients and sites.', '',
              '## Outputs and reproduction', '',
              '- [Machine-readable final results](results_final.json), including exact values,',
              '  all descriptive checks, claim limits and bound input hashes.',
              '- [Statistical review](review_2026-10-17/statistics/REVIEW.md).',
              '- [Provenance and model-rule review](review_2026-10-17/provenance_rules/REVIEW.md).',
              '- [Third recomputation](review_2026-10-17/statistics/recompute.py),',
              '  [Linux subset](review_2026-10-17/statistics/linux_subset_diagnostic.py), and',
              '  [independent provenance/counter audit](review_2026-10-17/provenance_rules/audit.py).', '',
              'To reproduce this assembly without rerunning statistics:', '',
              '```sh', 'python3 finalize_results.py --output /path/to/fresh/results_directory', '```', '',
              'The numerical scripts require the archived local raw inputs; paths and runtime',
              'details are preserved in their code/results. No new simulation, cluster task,',
              'gate change, commit, publication or external message is implied by this report.', '']
    # Keep ordinary prose spacing around numeric tokens readable.
    text = '\n'.join(lines)
    for old, new in [('independent67','independent 67'),('clustered55','clustered 55'),
                     ('imperfect27','imperfect 27'),('descriptive95','descriptive 95'),
                     ('target2.7','target 2.7'),('in72.0','in 72.0'),('and84.1','and 84.1'),
                     ('of192','of 192'),('gives66.64','gives 66.64'),('with18.846','with 18.846'),
                     ('contains64','contains 64'),('Python3.12','Python 3.12'),('and192','and 192'),
                     ('sites,16','sites, 16'),('isolation/14','isolation / 14')]:
        text = text.replace(old, new)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'results_final.json').write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
    (args.output / 'RESULTS.md').write_text(text)
    print(json.dumps({'status': result['status'], 'outputs': {name: sha(args.output / name) for name in ['RESULTS.md', 'results_final.json']}}))


if __name__ == '__main__':
    main()
