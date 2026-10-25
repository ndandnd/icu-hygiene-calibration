#!/usr/bin/env python
"""Reconcile impl_a and impl_b of common_baseline_2026-10-16 (read-only diagnostics).

Reads both results.json files, aligns their (different) schemas, and compares
- every point estimate (protocol tolerance 1e-10 absolute; relative also reported),
- every bootstrap interval bound (protocol tolerance 0.5 percentage points),
- reproduction-gate recomputations against each other,
- frac_extrapolated (Monte Carlo quantity; difference vs binomial SE, no protocol tolerance).
Also recomputes all point estimates a third way, directly from the archived report
values (no raw records, no code from either implementation), as an anchor.

Writes reconcile/comparison.json only.
"""
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
A = json.loads((ROOT / 'impl_a/results.json').read_text())
B = json.loads((ROOT / 'impl_b/results.json').read_text())
RS_REPORT = (Path(__file__).resolve().parents[4] / 'abm/research/rule_surface_results_2026-10-12/completed_analysis/report.json')
SENS_REPORT = (Path(__file__).resolve().parents[4] / 'research/author_notes_2026-10-10/analysis_2026-10-11/sensitivity_replay_20261011_1015_utc/report.json')
PROFILES = ('independent', 'clustered', 'imperfect')
PAIRS = ('independent-clustered', 'independent-imperfect', 'clustered-imperfect')
DRAWS = 50000
TOL_POINT = 1e-10
TOL_CI = 0.005  # 0.5 percentage points on the F scale

points, intervals, extrap, notes = [], [], [], []


def pt(section, name, a, b):
    d = abs(a - b)
    points.append(dict(section=section, quantity=name, impl_a=a, impl_b=b, abs_diff=d,
                       rel_diff=d / abs(b) if b else d, within_tol=d <= TOL_POINT))


def ci(section, name, a, b, scale_pp=True):
    """a, b: [lo, hi]."""
    d = [abs(a[0] - b[0]), abs(a[1] - b[1])]
    intervals.append(dict(section=section, quantity=name, impl_a=a, impl_b=b,
                          abs_diff_lo=d[0], abs_diff_hi=d[1], max_abs_diff=max(d),
                          max_abs_diff_pp=100 * max(d) if scale_pp else None,
                          on_F_scale=scale_pp,
                          within_tol=(max(d) <= TOL_CI) if scale_pp else None))


def fx(section, name, fa, fb):
    se = math.sqrt(max(fa * (1 - fa) + fb * (1 - fb), 1e-300) / DRAWS)
    extrap.append(dict(section=section, quantity=name, impl_a=fa, impl_b=fb, diff=fa - fb,
                       binomial_se_of_diff=se, z=(fa - fb) / se if se > 0 else 0.0))


# ------------------------------------------------------------- primary
Ap = A['primary']['profiles']
Bp = B['primary']
for h in PROFILES:
    for k in ('B_lo', 'B_hi', 'P_lo', 'P_hi', 'F_lo', 'F_hi', 'p_lo', 'p_hi', 'w', 'p_star', 'F_star'):
        pt('primary', f'{h}.{k}', Ap[h][k], Bp[h][k])
    ci('primary', f'{h}.F_star_ci95', [Ap[h]['ci95_lo'], Ap[h]['ci95_hi']], Bp[h]['ci95'])
    ci('primary', f'{h}.p_star_ci95', Ap[h]['p_star_ci95'], Bp[h]['p_star_ci95'], scale_pp=False)
    ci('primary', f'{h}.w_ci95', Ap[h]['w_ci95'], Bp[h]['w_ci95'], scale_pp=False)
    fx('primary', f'{h}.frac_extrapolated', Ap[h]['frac_extrapolated'], Bp[h]['frac_extrapolated'])
    pt('primary', f'{h}.undefined_draws', Ap[h]['undefined_draws'], Bp[h]['undefined_draws'])
    pt('primary', f'{h}.point_extrapolates', int(Ap[h]['point_extrapolated']), int(Bp[h]['point_extrapolates']))
Ac = {c['pair']: c for c in A['contrasts']['pairs']}
Bc = {c['pair']: c for c in B['contrasts']}
for p in PAIRS:
    pt('primary_contrasts', f'{p}.estimate', Ac[p]['estimate'], Bc[p]['estimate'])
    ci('primary_contrasts', f'{p}.ci98.333', [Ac[p]['ci_lo'], Ac[p]['ci_hi']], Bc[p]['ci'])
    pt('primary_contrasts', f'{p}.level', Ac[p]['level'], Bc[p]['level'])
    pt('primary_contrasts', f'{p}.resolved', int(Ac[p]['resolved']), int(Bc[p]['resolved']))

# ------------------------------------------------------------- S1
for h in PROFILES:
    a = A['S1']['profiles'][h]
    b = B['sensitivity']['S1'][h]
    pt('S1', f'{h}.F_star', a['F_star'], b['F_star'])
    pt('S1', f'{h}.minus_primary', a['minus_primary'], b['minus_primary'])
    ci('S1', f'{h}.F_star_ci95', [a['ci95_lo'], a['ci95_hi']], b['ci95'])
notes.append('S1 contrasts: computed by impl_a only (not required by protocol); no impl_b counterpart.')

# ------------------------------------------------------------- S2
for h in ('clustered', 'imperfect'):
    a = A['S2']['profiles'][h]
    b = B['sensitivity']['S2'][h]
    pt('S2', f'{h}.F_star', a['F_star'], b['F_star'])
    pt('S2', f'{h}.p_star', a['p_star'], b['p_star'])
    pt('S2', f'{h}.minus_primary', a['minus_primary'], b['minus_primary'])
    # impl_a stores (c0,c1,c2) for c0+c1x+c2x^2; impl_b stores np.polyfit order (c2,c1,c0)
    for i, nm in enumerate(('c0', 'c1', 'c2')):
        av = a['lnB_quadratic_coefs_in_lnp'][i]
        bv = b['quadratic_coef_lnB_in_lnp'][2 - i]
        d = abs(av - bv)
        # polynomial coefficients in raw ln p are ill-conditioned (|ln p|~2.6); report relative diff
        points.append(dict(section='S2', quantity=f'{h}.lnB_quadratic_{nm}', impl_a=av, impl_b=bv, abs_diff=d,
                           rel_diff=d / abs(bv), within_tol=d / abs(bv) <= 1e-10,
                           note='raw-ln-p polynomial coefficient; compared on relative scale (ill-conditioned basis)'))
    pt('S2', f'{h}.undefined_draws', a['undefined_draws'], b['undefined_draws'])
    ci('S2', f'{h}.F_star_ci95', [a['ci95_lo'], a['ci95_hi']], b['ci95'])
notes.append('S2 clustered-imperfect contrast: computed by impl_a only (not required by protocol); no impl_b counterpart.')

# ------------------------------------------------------------- S3 (main: A/B 256 banks)
A3 = A['S3']['main_AB256']['profiles']
B3 = B['sensitivity']['S3']['profiles']
for h in PROFILES:
    for k in ('B_lo', 'B_hi', 'F_lo', 'F_hi', 'p_lo', 'p_hi', 'w', 'p_star', 'F_star'):
        pt('S3', f'{h}.{k}', A3[h][k], B3[h][k])
    ci('S3', f'{h}.F_star_ci95', [A3[h]['ci95_lo'], A3[h]['ci95_hi']], B3[h]['ci95'])
    fx('S3', f'{h}.frac_extrapolated', A3[h]['frac_extrapolated'], B3[h]['frac_extrapolated'])
A3c = {c['pair']: c for c in A['S3']['main_AB256']['contrasts']}
B3c = {c['pair']: c for c in B['sensitivity']['S3']['contrasts_descriptive']}
for p in PAIRS:
    pt('S3_contrasts', f'{p}.estimate', A3c[p]['estimate'], B3c[p]['estimate'])
    ci('S3_contrasts', f'{p}.ci98.333', [A3c[p]['ci_lo'], A3c[p]['ci_hi']], B3c[p]['ci'])
    pt('S3_contrasts', f'{p}.excludes_zero', int(A3c[p]['resolved']), int(B3c[p]['resolved']))
# S3 variant: unstratified 512 (both ran it)
A3p = A['S3']['variant_pooled512']['profiles']
B3u = B['sensitivity']['S3']['alternative_conventions_descriptive']['unstratified_512']['profiles']
for h in PROFILES:
    pt('S3_unstratified512', f'{h}.F_star', A3p[h]['F_star'], B3u[h]['F_star'])
    ci('S3_unstratified512', f'{h}.F_star_ci95', [A3p[h]['ci95_lo'], A3p[h]['ci95_hi']], B3u[h]['ci95'])
    fx('S3_unstratified512', f'{h}.frac_extrapolated', A3p[h]['frac_extrapolated'], B3u[h]['frac_extrapolated'])
A3pc = {c['pair']: c for c in A['S3']['variant_pooled512']['contrasts']}
B3uc = {c['pair']: c for c in B['sensitivity']['S3']['alternative_conventions_descriptive']['unstratified_512']['contrasts']}
for p in PAIRS:
    ci('S3_unstratified512_contrasts', f'{p}.ci98.333', [A3pc[p]['ci_lo'], A3pc[p]['ci_hi']], B3uc[p]['ci'])
notes.append('S3 four-strata-128 convention: impl_b only (descriptive); no impl_a counterpart.')

# ------------------------------------------------------------- S4
for T in ('2.7', '3.3'):
    for h in PROFILES:
        a = A['S4']['targets'][T][h]
        b = B['sensitivity']['S4'][T]['profiles'][h]
        for k in ('w', 'p_star', 'F_star'):
            pt(f'S4_{T}', f'{h}.{k}', a[k], b[k])
        pt(f'S4_{T}', f'{h}.point_extrapolates', int(a['point_extrapolated']), int(b['point_extrapolates']))
        ci(f'S4_{T}', f'{h}.F_star_ci95', [a['ci95_lo'], a['ci95_hi']], b['ci95'])
        fx(f'S4_{T}', f'{h}.frac_extrapolated', a['frac_extrapolated'], b['frac_extrapolated'])
    ac = {c['pair']: c for c in A['S4']['contrasts'][T]}
    bc = {c['pair']: c for c in B['sensitivity']['S4'][T]['contrasts']}
    for p in PAIRS:
        pt(f'S4_{T}_contrasts', f'{p}.estimate', ac[p]['estimate'], bc[p]['estimate'])
        ci(f'S4_{T}_contrasts', f'{p}.ci98.333', [ac[p]['ci_lo'], ac[p]['ci_hi']], bc[p]['ci'])
        pt(f'S4_{T}_contrasts', f'{p}.excludes_zero', int(ac[p]['resolved']), int(bc[p]['resolved']))

# ------------------------------------------------------------- reproduction gate: A vs B recomputed values
name_map = {c['condition']: c['condition'] for c in A['reproduction_gate']['conditions']}
gate_cmp = []
for c in A['reproduction_gate']['conditions']:
    nm = c['condition']
    if nm.endswith(' [holdout report]'):
        bname = 'holdout_report_' + nm.replace('old_full_', '').replace(' [holdout report]', '')
    else:
        bname = nm
    bc = B['reproduction_gate']['conditions'][bname]
    row = dict(condition=nm, impl_b_key=bname, a_passed=c['passed'], b_passed=bc['passed'],
               a_max_rel_err_vs_archived=c['max_rel_err'], b_max_rel_err_vs_archived=bc['max_rel_err'])
    for k in ('control_rate', 'policy_rate', 'F'):
        row[f'{k}_a_minus_b_rel'] = abs(c[k] - bc['recomputed'][k]) / abs(bc['recomputed'][k])
        row[f'{k}_archived_same'] = c[f'archived_{k}'] == bc['archived'][k]
    gate_cmp.append(row)
gate_summary = dict(
    impl_a_passed=A['reproduction_gate']['passed'], impl_b_passed=B['reproduction_gate']['passed'],
    n_conditions_a=len(A['reproduction_gate']['conditions']), n_conditions_b=len(B['reproduction_gate']['conditions']),
    max_rel_err_vs_archived_a=max(c['max_rel_err'] for c in A['reproduction_gate']['conditions']),
    max_rel_err_vs_archived_b=B['reproduction_gate']['max_rel_err'],
    max_rel_diff_a_vs_b_recomputed=max(r[f'{k}_a_minus_b_rel'] for r in gate_cmp for k in ('control_rate', 'policy_rate', 'F')),
    archived_targets_identical=all(r[f'{k}_archived_same'] for r in gate_cmp for k in ('control_rate', 'policy_rate', 'F')),
    impl_a_hash_mismatches=len(A['reproduction_gate']['hash_mismatches']),
    impl_a_hash_bulk=[(h['label'], h['n_files'], h['mismatches']) for h in A['reproduction_gate']['hash_checks']],
    impl_a_single_file_hash_checks_all_match=all(x['match'] for x in A['reproduction_gate']['single_file_hash_checks']),
    impl_a_n_single_file_hash_checks=len(A['reproduction_gate']['single_file_hash_checks']),
    impl_b_hash_checks_all_passed=all(v['passed'] for v in B['hash_checks'].values()),
    conditions=gate_cmp)

# ------------------------------------------------------------- third recomputation from archived report values
rs = json.loads(RS_REPORT.read_text())
sr = json.loads(SENS_REPORT.read_text())
arch = {}
for e in rs['scenarios']:
    arch[e['id']] = dict(B=e['control_rate'], P=e['policy_rate'], p=e['scenario']['p'])
for e in sr['scenarios']:
    arch[e['scenario']['id']] = dict(B=e['control_rate'], P=e['policy_rate'], p=e['scenario']['p'])


def est(lo, hi, T=3.0):
    lo, hi = arch[lo], arch[hi]
    w = (math.log(T) - math.log(lo['B'])) / (math.log(hi['B']) - math.log(lo['B']))
    Flo, Fhi = 1 - lo['P'] / lo['B'], 1 - hi['P'] / hi['B']
    Fs = 1 - (1 - Flo) ** (1 - w) * (1 - Fhi) ** w
    ps = math.exp((1 - w) * math.log(lo['p']) + w * math.log(hi['p']))
    return dict(w=w, F_star=Fs, p_star=ps, S1=(1 - w) * Flo + w * Fhi)


brk = dict(independent=('independent_center_matched256', 'independent_upper'),
           clustered=('clustered_center', 'clustered_high'), imperfect=('imperfect_center', 'imperfect_high'))
brk3 = dict(independent=('old_full_independent_target_candidate', 'independent_upper'),
            clustered=('old_full_persistent_target_candidate', 'clustered_high'),
            imperfect=('old_full_imperfect_target_candidate', 'imperfect_high'))
third = []


def t3(section, name, val, a, b):
    third.append(dict(section=section, quantity=name, from_archived=val, impl_a=a, impl_b=b,
                      max_abs_diff=max(abs(val - a), abs(val - b))))


for h in PROFILES:
    e = est(*brk[h])
    for k in ('w', 'F_star', 'p_star'):
        t3('primary', f'{h}.{k}', e[k], Ap[h][k], Bp[h][k])
    t3('S1', f'{h}.F_star', e['S1'], A['S1']['profiles'][h]['F_star'], B['sensitivity']['S1'][h]['F_star'])
    e3 = est(*brk3[h])
    for k in ('w', 'F_star', 'p_star'):
        t3('S3', f'{h}.{k}', e3[k], A3[h][k], B3[h][k])
    for T in (2.7, 3.3):
        e4 = est(*brk[h], T=T)
        for k in ('w', 'F_star', 'p_star'):
            t3(f'S4_{T}', f'{h}.{k}', e4[k], A['S4']['targets'][str(T)][h][k], B['sensitivity']['S4'][str(T)]['profiles'][h][k])
# S2 from archived values: Lagrange quadratic in ln p, root nearest log-linear p*
for h in ('clustered', 'imperfect'):
    ids = [f'{h}_low', f'{h}_center', f'{h}_high']
    x = np.array([math.log(arch[i]['p']) for i in ids])
    yB = np.array([math.log(arch[i]['B']) for i in ids])
    yP = np.array([math.log(arch[i]['P']) for i in ids])
    xc = x[1]
    u = x - xc  # centred basis for conditioning
    cB = np.linalg.solve(np.vander(u, 3), yB)
    cP = np.linalg.solve(np.vander(u, 3), yP)
    r = np.roots(cB - np.array([0, 0, math.log(3.0)]))
    r = [float(v.real) for v in r if abs(v.imag) < 1e-12 and u[0] <= v.real <= u[2]]
    xref = math.log(est(*brk[h])['p_star']) - xc
    ur = min(r, key=lambda v: abs(v - xref))
    F2 = 1 - math.exp(float(np.polyval(cP, ur))) / 3.0
    t3('S2', f'{h}.F_star', F2, A['S2']['profiles'][h]['F_star'], B['sensitivity']['S2'][h]['F_star'])
    t3('S2', f'{h}.p_star', math.exp(ur + xc), A['S2']['profiles'][h]['p_star'], B['sensitivity']['S2'][h]['p_star'])
    t3('S2', f'{h}.curvature_lnB(c2)', float(cB[0]), A['S2']['profiles'][h]['curvature_lnB'],
       B['sensitivity']['S2'][h]['quadratic_coef_lnB_in_lnp'][0])
    t3('S2', f'{h}.n_in_range_roots', len(r), 1, B['sensitivity']['S2'][h]['in_range_roots_at_point'])

# ------------------------------------------------------------- summary
pts_tol = [p for p in points if not p.get('note')]
summary = dict(
    point_tolerance=TOL_POINT, interval_tolerance_pp=100 * TOL_CI,
    n_point_comparisons=len(points), n_points_outside_tol=sum(not p['within_tol'] for p in points),
    max_point_abs_diff=max(p['abs_diff'] for p in pts_tol),
    max_point_abs_diff_item=max(pts_tol, key=lambda p: p['abs_diff'])['quantity'],
    n_F_scale_interval_comparisons=sum(1 for i in intervals if i['on_F_scale']),
    n_intervals_outside_tol=sum(1 for i in intervals if i['on_F_scale'] and not i['within_tol']),
    max_interval_bound_diff_pp=max(i['max_abs_diff_pp'] for i in intervals if i['on_F_scale']),
    max_interval_bound_diff_item=max((i for i in intervals if i['on_F_scale']), key=lambda i: i['max_abs_diff'])['section'] + ':' +
    max((i for i in intervals if i['on_F_scale']), key=lambda i: i['max_abs_diff'])['quantity'],
    max_interval_bound_diff_pp_by_section={s: max(i['max_abs_diff_pp'] for i in intervals if i['section'] == s and i['on_F_scale'])
                                           for s in dict.fromkeys(i['section'] for i in intervals if i['on_F_scale'])},
    mean_interval_bound_diff_pp=float(np.mean([100 * d for i in intervals if i['on_F_scale'] for d in (i['abs_diff_lo'], i['abs_diff_hi'])])),
    max_frac_extrapolated_abs_diff=max(abs(e['diff']) for e in extrap),
    max_frac_extrapolated_abs_z=max(abs(e['z']) for e in extrap),
    third_recomputation_max_abs_diff=max(t['max_abs_diff'] for t in third if 'n_in_range' not in t['quantity']),
    reproduction_gate=dict((k, v) for k, v in gate_summary.items() if k != 'conditions'),
)
out = dict(summary=summary, notes=notes, points=points, intervals=intervals, frac_extrapolated=extrap,
           reproduction_gate_comparison=gate_summary, third_recomputation_from_archived_values=third)
(HERE / 'comparison.json').write_text(json.dumps(out, indent=1))
print(json.dumps(summary, indent=1))
print('\nworst 8 F-scale interval bound diffs (pp):')
for i in sorted((i for i in intervals if i['on_F_scale']), key=lambda i: -i['max_abs_diff'])[:8]:
    print(f"  {i['section']:28s} {i['quantity']:40s} {i['max_abs_diff_pp']:.3f}")
print('\nfrac_extrapolated diffs:')
for e in extrap:
    print(f"  {e['section']:20s} {e['quantity']:32s} A={e['impl_a']:.5f} B={e['impl_b']:.5f} z={e['z']:+.2f}")
print('\nthird recomputation:')
for t in third:
    print(f"  {t['section']:8s} {t['quantity']:30s} {t['from_archived']!r:>24} maxdiff={t['max_abs_diff']:.2e}")
