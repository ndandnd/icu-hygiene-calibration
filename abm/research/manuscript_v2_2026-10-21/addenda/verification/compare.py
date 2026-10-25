"""Compare addenda.json (the other agent's output values) with verify_results.json
(this verifier's independent recomputation) and emit markdown tables.

Status rules
  point values / counts : 'confirmed' if they agree to the addenda's display
                          precision (and here, in fact, to <=1e-9); else DISCREPANT.
  bootstrap endpoints   : 'confirmed (MC)' if |diff| <= 0.03 x interval width + 0.01
                          (two independent 20,000-draw realizations with different
                          streams; addenda seed 20261021, verifier seed 777001).
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
A = json.load(open(os.path.join(HERE, '..', 'addenda.json')))
V = json.load(open(os.path.join(HERE, 'verify_results.json')))
rows = {}
stats = {'point': [0, 0], 'interval': [0, 0], 'count': [0, 0]}
max_rel = {}


def fmt(x, d):
    return f'{x:.{d}f}'


def point(sec, label, a, v, d=2, tol=1e-9):
    ok = abs(a - v) <= tol * max(1.0, abs(a))
    disp_ok = fmt(a, d) == fmt(v, d)
    st = 'confirmed' if ok else ('confirmed (display)' if disp_ok else 'DISCREPANT')
    if 'UNROUNDED' in label and st == 'DISCREPANT':
        st = 'DISCREPANT (minor: input rounding)'
    stats['point'][0 if st.startswith('confirmed') else 1] += 1
    rows.setdefault(sec, []).append((label, fmt(a, d), fmt(v, d), f'{a - v:+.1e}', st))


def count(sec, label, a, v):
    st = 'confirmed' if int(a) == int(v) else 'DISCREPANT'
    stats['count'][0 if st == 'confirmed' else 1] += 1
    rows.setdefault(sec, []).append((label, str(int(a)), str(int(v)), f'{int(a) - int(v):+d}', st))


def interval(sec, label, a, v, d=2, fam='other', small=False):
    a = [a['lo'], a['hi']] if isinstance(a, dict) else list(a)
    w = max(abs(a[1] - a[0]), 1e-12)
    diffs = [abs(a[0] - v[0]), abs(a[1] - v[1])]
    tol = (0.05 if small else 0.03) * w + 0.01
    st = ('confirmed (MC, small-count discreteness)' if small and max(diffs) > 0.03 * w + 0.01 else 'confirmed (MC)') \
        if max(diffs) <= tol else 'DISCREPANT'
    stats['interval'][0 if st != 'DISCREPANT' else 1] += 1
    max_rel[fam] = max(max_rel.get(fam, 0), max(diffs) / w)
    rows.setdefault(sec, []).append((label, f'[{fmt(a[0], d)}, {fmt(a[1], d)}]', f'[{fmt(v[0], d)}, {fmt(v[1], d)}]',
                                     f'{max(diffs):.2f}', st))


# ------------------------------------------------------------------ conditions
cmap = {'sens:clustered_low': 'sens_clustered_low', 'sens:clustered_center': 'sens_clustered_center',
        'sens:clustered_high': 'sens_clustered_high', 'sens:imperfect_low': 'sens_imperfect_low',
        'sens:imperfect_center': 'sens_imperfect_center', 'sens:imperfect_high': 'sens_imperfect_high',
        'ind:center_matched256': 'ind_center256', 'ind:upper': 'ind_upper',
        'env:independent_low075_env1': 'env_independent_low075_env1',
        'env:independent_low100_env1': 'env_independent_low100_env1',
        'env:independent_low125_env1': 'env_independent_low125_env1',
        'hold:independent': 'hold_independent', 'hold:clustered': 'hold_clustered', 'hold:imperfect': 'hold_imperfect',
        'alt:independent': 'alt_independent', 'alt:clustered': 'alt_clustered', 'alt:imperfect': 'alt_imperfect'}
S0 = '1a. Input conditions (ADDENDA.md section 1 input table)'
for ak, vk in cmap.items():
    a, v = A['conditions'][ak], V['conditions'][vk]
    point(S0, f'{ak} p', a['p'], v['p'], 6)
    count(S0, f'{ak} seeds', a['n_seeds'], v['n'])
    point(S0, f'{ak} B', a['control_rate'], v['B'], 3)
    interval(S0, f'{ak} B 95%', a['control_rate_ci95'], v['B95'], 3, 'cond_B')
    point(S0, f'{ak} F %', a['F_pct'], 100 * v['F'], 2)
    interval(S0, f'{ak} F 95%', a['F_pct_ci95'], v['F95'], 2, 'cond_F')

# ------------------------------------------------------------------ slopes
S1 = '1b. Slopes (ADDENDA.md section 1 slope table)'
smap = {('independent', 'primary'): 'independent/primary', ('independent', 'env_low_mid'): 'independent/env_seg1',
        ('independent', 'env_mid_high'): 'independent/env_seg2', ('independent', 'env_ols3'): 'independent/env_ols',
        ('clustered', 'primary'): 'clustered/primary', ('clustered', 'low_center'): 'clustered/low_center',
        ('clustered', 'ols3'): 'clustered/ols', ('imperfect', 'primary'): 'imperfect/primary',
        ('imperfect', 'low_center'): 'imperfect/low_center', ('imperfect', 'ols3'): 'imperfect/ols'}
for (prof, est), vk in smap.items():
    a, v = A['elasticity']['slopes'][prof][est], V['slopes'][vk]
    lab = f'{prof} {est}'
    point(S1, lab + ' dF/dlnB', a['beta_B'], v['dFdlnB'])
    interval(S1, lab + ' dF/dlnB 95%', a['beta_B_ci95'], v['dFdlnB95'], 2, 'slope')
    point(S1, lab + ' per +10%', a['beta_B_per_10pct'], v['per10'])
    interval(S1, lab + ' per +10% 95%', a['beta_B_per_10pct_ci95'], v['per10_95'], 2, 'slope')
    point(S1, lab + ' dF/dlnp', a['beta_p'], v['dFdlnp'])
    interval(S1, lab + ' dF/dlnp 95%', a['beta_p_ci95'], v['dFdlnp95'], 2, 'slope')
    point(S1, lab + ' dlnB/dlnp', a['gamma'], v['dlnBdlnp'])
    interval(S1, lab + ' dlnB/dlnp 95%', a['gamma_ci95'], v['dlnBdlnp95'], 2, 'slope')
S1c = '1c. Consistency with common-baseline F*'
for prof in ('independent', 'clustered', 'imperfect'):
    a, v = A['elasticity']['consistency_with_common_baseline'][prof], V['consistency'][prof]
    point(S1c, prof + ' F(center)+beta*ln(3/Bc) %', a['F_at_3_arithmetic_pct'], v['linear'], 2)
    point(S1c, prof + ' minus common-baseline F* (pp)', a['difference_pp'], v['diff'], 3, tol=1e-6)

# ------------------------------------------------------------------ bias
S2 = '1d. Implied bias (ADDENDA.md bias table)'
for wf in ('original', 'alternative'):
    for prof in ('independent', 'clustered', 'imperfect'):
        a, v = A['bias']['results'][wf][prof], V['bias'][f'{wf}/{prof}']
        lab = f'{wf} {prof}'
        point(S2, lab + ' B_obs', a['control_rate'], v['B_obs'], 3)
        point(S2, lab + ' vs 3.0 (%)', a['shortfall_pct'], v['pct_vs_3'], 1)
        pw, pr = a['by_slope']['piecewise'], a['by_slope']['primary']
        point(S2, lab + ' bias piecewise', pw['bias_pp'], v['piecewise'])
        interval(S2, lab + ' bias piecewise joint 95%', pw['bias_pp_ci95_joint'], v['piecewise95'], 2, 'bias')
        point(S2, lab + ' bias primary', pr['bias_pp'], v['primary'])
        interval(S2, lab + ' bias primary joint 95%', pr['bias_pp_ci95_joint'], v['primary95'], 2, 'bias')
        interval(S2, lab + ' bias primary beta-only 95%', pr['bias_pp_ci95_beta_only'], v['beta_only95'], 2, 'bias')
        rg = a['point_range_across_slope_estimators_pp']
        point(S2, lab + ' range min (all estimators)', rg[0], v['range'][0])
        point(S2, lab + ' range max (all estimators)', rg[1], v['range'][1])
        point(S2, lab + ' F observed %', pw['F_observed_pct'], v['F_obs'], 1)
        point(S2, lab + ' F at 3.0 (piecewise) %', pw['F_adjusted_to_3_pct'], v['F_at3_pw'], 1)
        interval(S2, lab + ' F at 3.0 joint 95%', pw['F_adjusted_to_3_pct_ci95_joint'], v['F_at3_pw95'], 1, 'bias')
S3 = '1e. Contrast bias and adjusted contrasts'
for wf in ('original', 'alternative'):
    for pair in ('independent-clustered', 'independent-imperfect', 'clustered-imperfect'):
        a, v = A['bias']['contrasts'][wf][pair], V['contrast_bias'][f'{wf}/{pair}']
        lab = f'{wf} {pair}'
        point(S3, lab + ' observed (pp)', a['observed_contrast_pp'], v['observed'], 1)
        interval(S3, lab + ' observed 98.333%', a['observed_contrast_ci98333'], v['observed98'], 1, 'contrast')
        point(S3, lab + ' contrast bias piecewise', a['piecewise']['contrast_bias_pp'], v['bias_pw'])
        interval(S3, lab + ' contrast bias piecewise joint 95%', a['piecewise']['contrast_bias_ci95_joint'], v['bias_pw95'], 2, 'contrast')
        point(S3, lab + ' contrast bias primary', a['primary']['contrast_bias_pp'], v['bias_prim'])
        interval(S3, lab + ' contrast bias primary joint 95%', a['primary']['contrast_bias_ci95_joint'], v['bias_prim95'], 2, 'contrast')
        point(S3, lab + ' adjusted (piecewise)', a['piecewise']['adjusted_contrast_pp'], v['adjusted_pw'], 1)
        interval(S3, lab + ' adjusted joint 98.333%', a['piecewise']['adjusted_contrast_ci98333_joint'], v['adjusted_pw98'], 1, 'contrast')

# ------------------------------------------------------------------ routes
S4 = '3. Route labels (ADDENDA.md section 3, both workflows)'
lmap = {'patient': 'patient', 'own_current_0': 'own0', 'own_current_pos': 'ownpos', 'own_current_all': 'own_any',
        'other_zone': 'other', 'carry_over': 'carry_all', 'carry_over_0': 'carry0', 'carry_over_pos': 'carrypos',
        'all_except_carry': 'all_except_carry', 'paper_metric_site_0': 'manuscript_metric', 'all': 'all'}
for wf in ('original', 'alternative'):
    for prof in ('independent', 'clustered', 'imperfect'):
        a, v = A['routes']['results'][wf][prof], V['routes'][f'{wf}/{prof}']
        point(S4, f'{wf} {prof} patient-day rel. diff (%)', 100 * a['patient_days_relative_difference'], 100 * v['pd_rel_diff'], 4)
        for ak, vk in lmap.items():
            x, y = a['labels'][ak], v[vk]
            lab = f'{wf} {prof} {ak}'
            count(S4, lab + ' control n', x['control'], y['ctrl'])
            count(S4, lab + ' policy n', x['policy'], y['pol'])
            point(S4, lab + ' control share %', x['control_share_pct'], y['share_ctrl'], 1)
            point(S4, lab + ' policy share %', x['policy_share_pct'], y['share_pol'], 1)
            if x['reduction_count_pct'] is not None:
                point(S4, lab + ' reduction %', x['reduction_count_pct'], y['reduction'], 1)
                interval(S4, lab + ' reduction 95%', x['reduction_count_pct_ci95'], y['reduction95'], 1, 'route', small=x['control'] < 100)
            point(S4, lab + ' contribution to F (pp)', x['contribution_to_F_pp'], y['contribution'], 1)
            if x['small_count_flag'] != y['dagger']:
                rows[S4].append((lab + ' dagger flag', str(x['small_count_flag']), str(y['dagger']), '', 'DISCREPANT'))

# ------------------------------------------------------------------ selection (arithmetic from addenda inputs + own bootstrap)
S5 = '2. Selection propagation (ADDENDA.md section 2)'
for wf in ('original', 'alternative'):
    for prof in ('independent', 'clustered', 'imperfect'):
        a = A['selection']['results'][wf][prof]
        v = V['selection'][f'{wf}/{prof}']
        lab = f'{wf} {prof}'
        p0, (lo, hi) = a['p_selected'], a['p_ci95']
        al, ah = math.log(lo / p0), math.log(hi / p0)
        b, bl, bh = a['dF_dlnp'], a['dF_dlnp_ci95']['lo'], a['dF_dlnp_ci95']['hi']
        # arithmetic using the addenda's own inputs
        point(S5, lab + ' interval vs p, low (%)', a['dlnp_pct'][0], 100 * (lo / p0 - 1), 1)
        point(S5, lab + ' interval vs p, high (%)', a['dlnp_pct'][1], 100 * (hi / p0 - 1), 1)
        point(S5, lab + ' dF/dlnp (vs own recompute)', b, v['slope'], 2)
        interval(S5, lab + ' dF/dlnp 95% (own bootstrap)', a['dF_dlnp_ci95'], v['slope95'], 1, 'slope')
        sp = sorted([b * al, b * ah])
        point(S5, lab + ' shift point slope, min (arith.)', a['implied_F_shift_pp_point'][0], sp[0])
        point(S5, lab + ' shift point slope, max (arith.)', a['implied_F_shift_pp_point'][1], sp[1])
        interval(S5, lab + ' shift with slope uncertainty (own draws, per-draw envelope)',
                 a['implied_F_shift_pp_with_slope_uncertainty'], v['shift_slopeunc_envelope'], 2, 'selection')
        S = max(abs(x * y) for x in (b, bl, bh) for y in (al, ah))
        point(S5, lab + ' S (arith.)', a['selection_half_width_pp_envelope'], S, 2)
        point(S5, lab + ' MC half-width (archived, rounded-interval)', a['conditional_MC_half_width_pp'], v['mc_half_rounded_iv'], 2)
        point(S5, lab + ' MC half-width vs UNROUNDED archived interval', a['conditional_MC_half_width_pp'], v['mc_half'], 2, tol=1e-4)
        mc = a['conditional_MC_half_width_pp']
        point(S5, lab + ' combined sqrt(MC^2+S^2) (arith., addenda MC input)', a['quadrature_half_width_pp_envelope'], math.sqrt(mc * mc + S * S), 2)
        point(S5, lab + ' combined with UNROUNDED archived MC half-width', a['quadrature_half_width_pp_envelope'], math.sqrt(v['mc_half'] ** 2 + S * S), 2, tol=1e-4)
        if 'sensitivity_beta_B_times_alt_gamma' in a:
            s, g = a['sensitivity_beta_B_times_alt_gamma'], v['alt_grid']
            point(S5, lab + ' alt grid B m1', s['alt_grid_rates'][0], g['B_m1'], 2)
            point(S5, lab + ' alt grid B m2', s['alt_grid_rates'][1], g['B_m2'], 2)
            point(S5, lab + ' alt grid dlnB/dlnp', s['alt_grid_gamma'], g['dlnBdlnp'], 2)
            interval(S5, lab + ' alt grid dlnB/dlnp 95%', s['alt_grid_gamma_ci95'], g['dlnBdlnp95'], 2, 'slope')
            point(S5, lab + ' dF/dlnp with alt gamma', s['dF_dlnp'], g['dFdlnp'], 1)
            interval(S5, lab + ' dF/dlnp with alt gamma 95%', s['dF_dlnp_ci95'], g['dFdlnp95'], 1, 'slope')
            sp2 = sorted(g['shift_point'])
            point(S5, lab + ' shift (alt gamma) min', s['implied_F_shift_pp_point'][0], sp2[0])
            point(S5, lab + ' shift (alt gamma) max', s['implied_F_shift_pp_point'][1], sp2[1])

# ------------------------------------------------------------------ headline
S6 = '0. Headline numbers'
h = A['headline']
allb = [V['bias'][k][e] for k in V['bias'] for e in ('piecewise', 'primary')]
point(S6, 'max abs bias, local estimators', h['max_abs_bias_pp_local_estimators'], max(abs(x) for x in allb))
point(S6, 'max abs bias, any estimator', h['max_abs_bias_pp_any_estimator'],
      max(max(abs(x) for x in V['bias'][k]['range']) for k in V['bias']))
lo = min(min(V['bias'][k]['piecewise95'][0], V['bias'][k]['primary95'][0]) for k in V['bias'])
hi = max(max(V['bias'][k]['piecewise95'][1], V['bias'][k]['primary95'][1]) for k in V['bias'])
interval(S6, 'joint 95% envelope, local estimators', h['joint_ci95_range_pp_local_estimators'], [lo, hi], 2, 'bias')
point(S6, 'max abs contrast bias, local', h['max_abs_contrast_bias_pp_local'],
      max(max(abs(V['contrast_bias'][k]['bias_pw']), abs(V['contrast_bias'][k]['bias_prim'])) for k in V['contrast_bias']))
lo = min(min(V['contrast_bias'][k]['bias_pw95'][0], V['contrast_bias'][k]['bias_prim95'][0]) for k in V['contrast_bias'])
hi = max(max(V['contrast_bias'][k]['bias_pw95'][1], V['contrast_bias'][k]['bias_prim95'][1]) for k in V['contrast_bias'])
interval(S6, 'contrast-bias joint 95% envelope', h['contrast_bias_joint_ci95_range_pp_local'], [lo, hi], 2, 'contrast')
excl = all(V['contrast_bias'][k]['adjusted_pw98'][0] > 0 for k in V['contrast_bias'])
rows[S6].append(('all six adjusted contrasts exclude 0 (98.333%)', str(h['all_adjusted_contrasts_exclude_zero_98333']), str(excl), '',
                 'confirmed' if excl == h['all_adjusted_contrasts_exclude_zero_98333'] else 'DISCREPANT'))
inc = all(V['slopes'][k][m][0] < 0 < V['slopes'][k][m][1] for k in V['slopes'] for m in ('dFdlnB95',))
rows[S6].append(('every dF/dlnB 95% interval includes 0', str(h['all_slope_ci95_include_zero']), str(inc), '',
                 'confirmed' if inc == h['all_slope_ci95_include_zero'] else 'DISCREPANT'))
pcts = [x for k in V['selection'] for x in V['selection'][k]['pct']]
point(S6, 'selected-p interval, min % of p', h['selected_p_interval_pct_range'][0], min(pcts), 1)
point(S6, 'selected-p interval, max % of p', h['selected_p_interval_pct_range'][1], max(pcts), 1)
point(S6, 'max abs selection shift, point slope', h['max_abs_selection_F_shift_pp_point_slope'],
      max(abs(x) for k in V['selection'] for x in V['selection'][k]['shift_point']))
point(S6, 'max abs selection shift, with slope uncertainty (own draws, MC)', h['max_abs_selection_F_shift_pp_with_slope_uncertainty'],
      max(abs(x) for k in V['selection'] for x in V['selection'][k]['shift_slopeunc_envelope']), 2, tol=0.01)

out = []
for sec in sorted(rows):
    out.append(f'\n### {sec}\n')
    out.append('| Item | Addenda | Verifier | abs diff | Status |')
    out.append('| --- | ---: | ---: | ---: | --- |')
    for r in rows[sec]:
        out.append('| ' + ' | '.join(r) + ' |')
open(os.path.join(HERE, 'comparison_tables.md'), 'w').write('\n'.join(out) + '\n')
summary = dict(stats=stats, max_endpoint_diff_over_width=max_rel,
               discrepant=[(sec, r[0], r[1], r[2]) for sec in rows for r in rows[sec] if r[4].startswith('DISCREPANT')],
               display_only=[(sec, r[0], r[1], r[2]) for sec in rows for r in rows[sec] if r[4] == 'confirmed (display)'])
json.dump(summary, open(os.path.join(HERE, 'comparison_summary.json'), 'w'), indent=1)
print(json.dumps(summary, indent=1))
