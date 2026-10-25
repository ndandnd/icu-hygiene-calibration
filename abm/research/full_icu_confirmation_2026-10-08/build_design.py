#!/usr/bin/env python3
"""Deterministic confirmation configs and baseline-only selection; never simulate."""
import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]
GRID = (.13496571354322734, .15, .165)
OLD_GRID = (.08097942812593639, .10797257083458187, .13496571354322734)
FIXED_P = {'persistent': .07602158239372622, 'imperfect': .07483995991409355}
PROFILE = {'independent': (0., 1.), 'persistent': (.5, 1.), 'imperfect': (0., .8)}
FAMILIES = {'calibration': 128, 'holdout': 512, 'environment': 512,
            'agent_crn_pilot_optional': 64}


def finite_number(value, label, low=0., high=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f'{label} must be finite numeric')
    if value < low or (high is not None and value > high):
        raise ValueError(f'{label} outside allowed bounds')
    return value


def choose(points, target=3., fallback=.15):
    """Pure bounded selection from exactly the declared three measurements."""
    finite_number(target, 'target')
    if target <= 0: raise ValueError('target must be positive')
    if fallback != .15: raise ValueError('Only the declared .15 fallback is permitted')
    checked = sorted((finite_number(p, 'p', high=1.), finite_number(r, 'rate'))
                     for p, r in points)
    if tuple(p for p, _ in checked) != GRID:
        raise ValueError('Expected exactly the declared three-point grid')
    exact = [p for p, r in checked if r == target]
    brackets = [(a, b) for a, b in zip(checked, checked[1:])
                if a[1] < target < b[1]]
    pair = None; weights = None
    if exact:
        p = exact[0]; method = 'lowest_p_exact_hit'; matched = True
    elif brackets:
        pair = brackets[0]
        (pl, rl), (pu, ru) = pair
        w = (target-rl)/(ru-rl)
        p = pl+(pu-pl)*w; weights = [1-w, w]
        method = 'lowest_p_increasing_adjacent_bracket'; matched = True
    else:
        p = fallback; method = 'fixed_reference_after_calibration_failure'; matched = False
    return {'target_rate_per_1000_patient_days': target,
            'points': [{'p': p, 'rate': r} for p, r in checked],
            'candidate_p': p, 'method': method, 'valid_target_calibration': matched,
            'selected_bracket': pair, 'weights': weights, 'exact_hits': exact,
            'increasing_brackets': brackets,
            'nonmonotone': any(b[1] < a[1] for a, b in zip(checked, checked[1:])),
            'allowed_extension': None,
            'scope': 'Baseline-only candidate selection, not verified baseline matching.'}


def select_baseline_report(report, expected_seeds):
    """Validate a core-report-shaped three-point baseline-only calibration report."""
    if len(expected_seeds) != 128 or len(set(expected_seeds)) != 128:
        raise ValueError('Expected exact 128 distinct calibration seeds')
    rows = report.get('scenarios', [])
    if len(rows) != 3: raise ValueError('Expected exactly three calibration rows')
    points = []
    for row in rows:
        s = row['scenario']
        if (s.get('h0'), s.get('h1'), s.get('rho'), s.get('e'),
            s.get('opportunities'), s.get('environment_multiplier', 1.)) != (.9, .9, 0., 1., 'zone_separate', 1.):
            raise ValueError('Expected independent baseline-only full-environment mechanism')
        if s.get('initial_environment_contaminated', False) is not False:
            raise ValueError('Calibration must start clean')
        if row.get('n') != 128 or sorted(row.get('seeds', [])) != sorted(expected_seeds):
            raise ValueError('Calibration report must contain the exact new 128 seeds')
        a = row['arms']['baseline']
        count = a['acquired_total']; exposure = a['patient_days_total']
        if type(count) is not int or count < 0: raise ValueError('Invalid baseline count')
        finite_number(exposure, 'baseline patient-days')
        if exposure == 0: raise ValueError('Positive baseline exposure required')
        rate = 1000*count/exposure
        stored = finite_number(a['rate_per_1000_patient_days'], 'stored baseline rate')
        if not math.isclose(rate, stored, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError('Stored baseline rate disagrees with count/exposure')
        points.append((s['p'], rate))
    return choose(points)


def seedbook(proposal):
    blocks = deepcopy(proposal['seed_blocks'])
    used = set()
    allocations = {}
    for name, n in FAMILIES.items():
        banks = blocks[name]
        if set(banks) != {'A', 'B'}: raise ValueError('Expected banks A and B')
        for bank in ('A', 'B'):
            values = banks[bank]
            if len(values) != n//2: raise ValueError('Incorrect bank size')
            for seed in values:
                if type(seed) is not int or not 0 < seed < 2**32 or seed in used:
                    raise ValueError('Seeds must be unique positive uint32 values across banks/families')
                used.add(seed)
        a, b = banks['A'], banks['B']
        allocations[name] = {'all': a+b, 'local': a[:n//16]+b[:n//16],
                             'cluster': a[n//16:]+b[n//16:]}
    return {'schema_version': 1, 'study': HERE.name, 'seed_blocks': blocks,
            'allocations': allocations, 'seed_generation': proposal['seed_generation'],
            'bank_site_note': 'One eighth of each bank local, seven eighths cluster. Banks are sample partitions, not site effects.',
            'optional_note': 'Reserved agent CRN seeds are not part of the primary execution budget.'}


def scenario(identifier, p, profile='independent', h1=.95, env=1.):
    finite_number(p, 'p', high=1.)
    rho, e = PROFILE[profile]
    return {'id': identifier, 'p': p, 'h0': .9, 'h1': h1, 'rho': rho,
            'e': e, 'opportunities': 'zone_separate', 'environment_multiplier': env}


def config(scenarios, seeds, description):
    return {'description': description, 'warmup_days': 90, 'observation_days': 56,
            'block_days': 1, 'model_overrides': {'prob_patient_sick': .03,
                                               'initial_environment_contaminated': False},
            'scenarios': scenarios, 'seeds': seeds}


def holdout_candidates():
    return {'schema_version': 1, 'status': 'independent_candidate_awaits_baseline_calibration',
            'study': HERE.name, 'runnable': False,
            'candidates': {'independent': {'p': None, 'rho': 0., 'e': 1.},
                           **{name: {'p': p, 'rho': PROFILE[name][0], 'e': PROFILE[name][1]}
                              for name, p in FIXED_P.items()}},
            'selection_grid': list(GRID), 'fallback': .15,
            'holdout_seeds_per_profile': 512, 'expected_paired_jobs': 1536}


def build_initial(proposal):
    book = seedbook(proposal)
    calibration = [scenario(f'independent_calibration_p{i}', p, h1=.9)
                   for i, p in enumerate(GRID)]
    environment = [scenario(f'independent_low{label}_env{env}', p, env=float(env))
                   for label, p in zip(('075', '100', '125'), OLD_GRID) for env in (1, 0)]
    outputs = {'seedbook.json': book, 'holdout_candidates.json': holdout_candidates()}
    for family, cases in (('calibration', calibration), ('environment', environment)):
        for site, seeds in book['allocations'][family].items():
            outputs[f'configs/{family}_{site}.json'] = config(cases, seeds, f'{HERE.name}: {family}, {site}; fixed finite horizon')
    outputs['output_scopes.json'] = {
        'schema_version': 1, 'study': HERE.name,
        'backend': 'research/full_icu_response_2026-10-08/run.py',
        'results_root': f'results/{HERE.name}',
        'banks': {f'{family}_{site}': f'results/{HERE.name}/{family}_{site}'
                  for family in ('calibration', 'holdout', 'environment') for site in ('local', 'cluster')},
        'analysis_namespaces': {family: f'research/{HERE.name}/{family}_report'
                                for family in ('calibration', 'holdout', 'environment')},
        'backend_metadata_note': 'Backend bank study label is historical; new-study wrapper must bind protocol/config/seedbook/output scope before reporting.',
        'counts': {'calibration_baseline_jobs': 384, 'holdout_paired_jobs': 1536,
                   'environment_paired_jobs': 3072, 'primary_paired_jobs': 4608},
        'no_simulation_in_builder': True}
    return outputs


def build_holdout(proposal, report):
    book = seedbook(proposal)
    selected = select_baseline_report(report, book['allocations']['calibration']['all'])
    ps = {'independent': selected['candidate_p'], **FIXED_P}
    cases = [scenario(f'{name}_target_candidate', p, name) for name, p in ps.items()]
    outputs = {'calibration_selection.json': selected}
    for site, seeds in book['allocations']['holdout'].items():
        outputs[f'configs/holdout_{site}.json'] = config(cases, seeds, f'{HERE.name}: fresh locked holdout, {site}')
    return outputs


def write_outputs(outputs, destination):
    """Preflight every path, then create new artifacts; never overwrite an earlier selection."""
    encoded = {name: json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'
               for name, value in outputs.items()}
    for name in encoded:
        if (destination/name).exists(): raise FileExistsError(destination/name)
    for name, raw in encoded.items():
        p = destination/name; p.parent.mkdir(parents=True, exist_ok=True)
        with p.open('x') as stream: stream.write(raw)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proposal', type=Path, default=HERE/'candidate_design.json')
    parser.add_argument('--output', type=Path, default=HERE)
    parser.add_argument('--calibration-report', type=Path,
                        help='Omit to materialize initial artifacts; supply once to lock holdout candidates.')
    args = parser.parse_args()
    proposal = json.loads(args.proposal.read_text())
    if args.calibration_report is None:
        outputs = build_initial(proposal)
    else:
        raw = args.calibration_report.read_bytes()
        outputs = build_holdout(proposal, json.loads(raw))
        outputs['calibration_selection.json']['report_sha256'] = hashlib.sha256(raw).hexdigest()
        outputs['calibration_selection.json']['proposal_sha256'] = hashlib.sha256(args.proposal.read_bytes()).hexdigest()
    write_outputs(outputs, args.output)
    print(json.dumps({'status': 'artifacts_created_no_simulation', 'files': sorted(outputs)}))


if __name__ == '__main__': main()
