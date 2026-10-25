#!/usr/bin/env python
"""Render ADDENDA.md from addenda.json (written by addenda_analysis.py). Every number in the
Markdown comes from the JSON; nothing is typed by hand except the archived displays quoted
for comparison (manuscript/verification intervals), which are labelled as such."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
J = json.loads((HERE / 'addenda.json').read_text())
JSON_SHA = hashlib.sha256((HERE / 'addenda.json').read_bytes()).hexdigest()

PROFS = ('independent', 'clustered', 'imperfect')
PNAME = dict(independent='Independent misses', clustered='Clustered misses', imperfect='Imperfect clearing')
WNAME = dict(original='Original', alternative='Alternative')
# archived conditional Monte Carlo 95% intervals for F (manuscript / verification RESULTS.md), display values
ARCHIVED_F = dict(original=dict(independent=(67.6, 64.1, 70.8), clustered=(55.9, 52.1, 59.5),
                                imperfect=(27.0, 21.6, 32.1)),
                  alternative=dict(independent=(76.5, 74.1, 78.7), clustered=(53.7, 49.7, 57.5),
                                   imperfect=(27.3, 22.2, 31.9)))
ARCHIVED_CONTRASTS = dict(original=(11.7, 40.6, 28.9), alternative=(22.8, 49.2, 26.4))


def n(x, d=2, sign=False):
    if x is None:
        return 'n/a'
    if round(x, d) == 0:
        x = 0.0
    return f'{x:+.{d}f}' if sign else f'{x:.{d}f}'


def iv(c, d=2):
    if c is None:
        return 'n/a'
    return f'[{n(c["lo"], d)}, {n(c["hi"], d)}]'


def table(header, rows, align=None):
    align = align or ['l'] + ['r'] * (len(header) - 1)
    sep = ['---:' if a == 'r' else '---' for a in align]
    out = ['| ' + ' | '.join(header) + ' |', '| ' + ' | '.join(sep) + ' |']
    out += ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]
    return '\n'.join(out)


C = J['conditions']
EL = J['elasticity']
SL = EL['slopes']
LAB = EL['slope_labels']
BI = J['bias']['results']
CB = J['bias']['contrasts']
SE = J['selection']['results']
SI = J['selection']['inputs']
RO = J['routes']['results']
CR = J['compliance_range']
md = []
w = md.append

# ------------------------------------------------------------------ derived headline numbers
loc_keys = ('piecewise', 'primary')
pt_loc = [abs(BI[wf][p]['by_slope'][k]['bias_pp']) for wf in BI for p in PROFS for k in loc_keys]
pt_all = [abs(v['bias_pp']) for wf in BI for p in PROFS for v in BI[wf][p]['by_slope'].values()]
jt_loc = [x for wf in BI for p in PROFS for k in loc_keys
          for x in (BI[wf][p]['by_slope'][k]['bias_pp_ci95_joint']['lo'],
                    BI[wf][p]['by_slope'][k]['bias_pp_ci95_joint']['hi'])]
hw = [(ARCHIVED_F[wf][p][2] - ARCHIVED_F[wf][p][1]) / 2 for wf in ARCHIVED_F for p in PROFS]
contr = [c for v in ARCHIVED_CONTRASTS.values() for c in v]
all_slope_ci_cover_zero = all(s['beta_B_ci95']['lo'] < 0 < s['beta_B_ci95']['hi']
                              for p in SL.values() for s in p.values())
sel_pt = max(max(abs(x) for x in SE[wf][p]['implied_F_shift_pp_point']) for wf in SE for p in PROFS)
sel_sl = max(max(abs(x) for x in SE[wf][p]['implied_F_shift_pp_with_slope_uncertainty']) for wf in SE for p in PROFS)
adj_excl0 = all(CB[wf][c][k]['adjusted_contrast_ci98333_joint']['lo'] > 0 for wf in CB for c in CB[wf]
                for k in loc_keys)
max_cb = max(abs(CB[wf][c][k]['contrast_bias_pp']) for wf in CB for c in CB[wf] for k in loc_keys)
cb_jt = [x for wf in CB for c in CB[wf] for k in loc_keys
         for x in (CB[wf][c][k]['contrast_bias_ci95_joint']['lo'], CB[wf][c][k]['contrast_bias_ci95_joint']['hi'])]

# cross-check the render-side derivations against the analysis script's headline block
H = J['headline']
assert abs(H['max_abs_bias_pp_local_estimators'] - max(pt_loc)) < 1e-12
assert abs(H['max_abs_bias_pp_any_estimator'] - max(pt_all)) < 1e-12
assert H['joint_ci95_range_pp_local_estimators'] == [min(jt_loc), max(jt_loc)]
assert all(abs(x - y) < 1e-12 for x, y in zip(H['archived_MC_F_half_width_pp_range'], [min(hw), max(hw)]))
assert H['all_slope_ci95_include_zero'] == all_slope_ci_cover_zero
assert abs(H['max_abs_contrast_bias_pp_local'] - max_cb) < 1e-12
assert H['contrast_bias_joint_ci95_range_pp_local'] == [min(cb_jt), max(cb_jt)]
assert H['all_adjusted_contrasts_exclude_zero_98333'] == adj_excl0
assert abs(H['max_abs_selection_F_shift_pp_point_slope'] - sel_pt) < 1e-12
assert abs(H['max_abs_selection_F_shift_pp_with_slope_uncertainty'] - sel_sl) < 1e-12
assert [list(v) for v in ARCHIVED_CONTRASTS.values()] == list(H['archived_primary_contrasts_pp'].values())
assert all(list(ARCHIVED_F[wf][p][1:]) == H['archived_MC_F_ci95_pct'][wf][p] for wf in ARCHIVED_F for p in PROFS)

# ------------------------------------------------------------------ header
w('# Manuscript v2 addenda: baseline mismatch, selection uncertainty, route labels and study flow')
w('')
w(f'Analysis-only addenda for manuscript v2 (folder dated 21 October 2026; computed '
  f'{J["log"][0].split("analysis ")[-1]}). **No simulation was run.** Every number '
  f'below is recomputed from archived, hash-checked records of a synthetic agent-based ICU model, or read from '
  f'archived tables. Nothing here is clinical evidence. Machine-readable values are in '
  f'[addenda.json](addenda.json) (SHA256 `{JSON_SHA[:12]}…`), produced by '
  f'[addenda_analysis.py](addenda_analysis.py); this page is rendered from that file by '
  f'[render_addenda_md.py](render_addenda_md.py).')
w('')
w(f'Bootstrap: {J["bootstrap"]["draws"]:,} whole-seed draws, seed {J["bootstrap"]["seed"]}, percentile 95% '
  f'intervals. One independent random stream per seed bank ({J["bootstrap"]["streams"]}). Conditions that share '
  f'seeds share draws, so their covariance is kept. Rates are pooled counts over pooled patient-days in every '
  f'draw. No draw was undefined in any reported statistic unless a table says so.')
w('')
w('Results of the alternative-workflow rare-transmission coefficient study (`workflow_coefficient_2026-10-21`) '
  'are not yet available. They appear below only as **[PENDING]** placeholders.')
w('')

# ------------------------------------------------------------------ summary
w('## Summary')
w('')
w(f'1. **The baseline shortfalls imply small biases in F, of uncertain sign.** Using the local slope of F against '
  f'log baseline rate, each observed shortfall from 3.0 per 1,000 patient-days moves F by at most '
  f'{max(pt_loc):.2f} percentage points (pp) under the two local estimators (at most {max(pt_all):.2f} pp under '
  f'any estimator tried). Joint 95% intervals stay within [{min(jt_loc):.2f}, {max(jt_loc):.2f}] pp. For '
  f'comparison, the conditional Monte Carlo half-widths of F are {min(hw):.2f} to {max(hw):.2f} pp and the '
  f'prespecified between-profile contrasts are {min(contr):.1f} to {max(contr):.1f} pp.')
w(f'2. **The direction is not resolved.** '
  + ('Every slope interval, for every profile and every estimator, includes zero. '
     if all_slope_ci_cover_zero else 'Not every slope interval includes zero; see section 1. ')
  + 'The point slopes differ in sign between profiles and, for the independent profile, between the target-region '
    'segment and the lower-baseline environmental segments. These results bound the size of the bias; they do '
    'not establish whether matched baselines would raise or lower F.')
w(f'3. **The contrasts survive the adjustment.** Shifting each profile to 3.0 changes the between-profile '
  f'contrasts by at most {max_cb:.2f} pp (joint 95% intervals within [{min(cb_jt):.2f}, {max(cb_jt):.2f}] pp). '
  + ('All six adjusted contrasts (two workflows x three pairs) still exclude zero at 98.333%. '
     if adj_excl0 else 'Not all adjusted contrasts exclude zero at 98.333%; see section 1. ')
  + 'This is a post hoc, descriptive calculation, not a matched-baseline result.')
w(f'4. **Selection uncertainty in p is small next to Monte Carlo error.** The archived selected-p intervals '
  f'span {min(SE[wf][p]["dlnp_pct"][0] for wf in SE for p in PROFS):.1f}% to '
  f'+{max(SE[wf][p]["dlnp_pct"][1] for wf in SE for p in PROFS):.1f}% of p. Multiplied by dF/dlog p, they imply F shifts of at most {sel_pt:.2f} pp at '
  f'the point slope, or {sel_sl:.2f} pp allowing for slope uncertainty.')
w('5. **Route labels.** In the original workflow, the pre-admission carry-over label fell by far less under the '
  'policy than the other labels in the independent profile. Under the alternative workflow that label became '
  'rare. Current-stay own-zone site pickups occur almost entirely with zero intervening opportunities, yet fell '
  'about as much as the patient-pickup label. The labels are descriptive, not causal.')
w('')
w('The transfer of original-workflow slopes to the alternative workflow is an assumption, not an estimate. '
  'Section 1 lists the assumptions.')
w('')

# ------------------------------------------------------------------ 1. elasticity
w('## 1. Baseline-mismatch elasticity')
w('')
w('### What was estimated')
w('')
w('The question is how much the prevented fraction F (completion 90% to 95%) would change if the control '
  'acquisition rate (the "baseline", B, per 1,000 patient-days) were exactly 3.0 instead of its observed value. '
  'The slope used is')
w('')
w('> beta_B = dF / d ln B, in percentage points of F per unit of natural-log baseline.')
w('')
w('A 10% higher baseline changes F by beta_B x ln(1.1) = 0.095 x beta_B pp. The implied bias of the observed F is')
w('')
w('> bias = F(B_obs) - F(3.0) ~ beta_B x ln(B_obs / 3.0).')
w('')
w('A positive bias means the observed F is higher than F would be at 3.0. F at 3.0 is then about F_obs - bias.')
w('')
w('Slopes come from conditions that share seeds, so the finite differences are paired. They vary p and read '
  'off both B and F. All conditions are original-workflow, care-zone model runs at completion .90 vs .95.')
w('')
rows = []
bank_of = {
    'sens:': ('Author-notes sensitivity bank (A/B 128)', 256),
    'ind:': ('Original holdout seeds A[:128]+B[:128]', 256),
    'env:': ('Environmental comparison bank (A/B 256)', 512),
    'hold:': ('Original holdout (A/B 256)', 512),
    'alt:': ('Alternative policy bank', 512),
}
for k, c in C.items():
    pre = k.split(':')[0] + ':'
    rows.append([c['condition'], bank_of[pre][0], c['n_seeds'], f'{c["p"]:.6f}',
                 f'{n(c["control_rate"], 3)} {iv(c["control_rate_ci95"], 3)}',
                 f'{n(c["F_pct"])} {iv(c["F_pct_ci95"])}'])
w(table(['Condition', 'Seed bank', 'Seeds', 'p', 'Control rate B [95%]', 'F % [95%]'], rows,
        ['l', 'l', 'r', 'r', 'r', 'r']))
w('')
w('Notes on these inputs:')
w('')
w('- **Clustered and imperfect.** The six sensitivity conditions (0.75, 1 and 1.25 x the fixed p; supplement '
  'Table S2) share one bank of 256 seeds.')
w('- **Independent.** The matched center is the original holdout at p = .14328 restricted to the 256 seeds '
  '(first 128 of each holdout bank) used by the independent upper condition at 1.25 x p. This is the same pair '
  'the common-baseline analysis used.')
w('- **Environmental comparison.** Its three full-model (environment enabled) points share 512 seeds among '
  'themselves but not with the holdout. Their baselines (0.40 to 2.29) lie below the target region. They are '
  'reported as secondary slopes only.')
w('- **Holdout stratification.** The original holdout is resampled in four 128-seed strata (A-first, A-rest, '
  'B-first, B-rest) rather than the manuscript\'s two A/B strata. This keeps the overlap with the independent '
  'upper seeds exact. The recomputed holdout F intervals differ from the published 200,000-draw intervals by at '
  f'most {max(abs(C["hold:" + p]["F_pct_ci95"][e] - ARCHIVED_F["original"][p][1 + i]) for p in PROFS for i, e in enumerate(("lo", "hi"))):.2f} pp.')
w('- **Reproduction.** All point values reproduce the archived reports (sensitivity report, common-baseline '
  f'results, alternative `analysis.json`) to 1e-12. There were {J["reproduction_failures"]} reproduction failures '
  f'among {len(J["reproduction_checks"])} checks.')
w('')
w('### Slopes')
w('')
rows = []
order = dict(independent=('primary', 'env_low_mid', 'env_mid_high', 'env_ols3'),
             clustered=('primary', 'low_center', 'ols3'), imperfect=('primary', 'low_center', 'ols3'))
for p in PROFS:
    for k in order[p]:
        s = SL[p][k]
        rates = ' to '.join(n(x, 2) for x in s['control_rates'])
        rows.append([PNAME[p] if k == 'primary' else '', LAB[k], rates,
                     f'{n(s["beta_B"])} {iv(s["beta_B_ci95"])}',
                     f'{n(s["beta_B_per_10pct"])} {iv(s["beta_B_per_10pct_ci95"])}',
                     f'{n(s["beta_p"])} {iv(s["beta_p_ci95"])}',
                     f'{n(s["gamma"])} {iv(s["gamma_ci95"])}'])
w(table(['Profile', 'Estimator', 'Baselines spanned', 'dF/d ln B (pp) [95%]', 'pp per +10% baseline [95%]',
         'dF/d ln p (pp) [95%]', 'd ln B / d ln p [95%]'], rows, ['l', 'l', 'r', 'r', 'r', 'r', 'r']))
w('')
w('How to read the slopes:')
w('')
below = [f'{WNAME[wf].lower()} {p}' for wf in BI for p in PROFS if BI[wf][p]['by_slope']['piecewise']['observed_below_center']]
w('- **Primary estimator.** The center to 1.25 x center segment is the only one whose baselines bracket 3.0, '
  'for every profile. Some observed baselines lie slightly below the sensitivity center ('
  + ', '.join(below) + '). For those, the piecewise estimator in the next table uses the low-to-center slope '
  'for the part of the span below the center.')
w('- **No resolved slope.** '
  + ('All intervals include zero. ' if all_slope_ci_cover_zero else '')
  + f'The largest primary per-10% effect is {max(abs(SL[p]["primary"]["beta_B_per_10pct"]) for p in PROFS):.2f} pp.')
w(f'- **Independent bracket is wide.** Its primary segment spans a factor of '
  f'{SL["independent"]["primary"]["control_rates"][1] / SL["independent"]["primary"]["control_rates"][0]:.1f} in '
  f'baseline (2.78 to 7.29), against {SL["clustered"]["primary"]["control_rates"][1] / SL["clustered"]["primary"]["control_rates"][0]:.1f} '
  f'and {SL["imperfect"]["primary"]["control_rates"][1] / SL["imperfect"]["primary"]["control_rates"][0]:.1f} for '
  f'the others. The environmental full-model segments below the target have the opposite point sign. Local '
  f'linearity is therefore weakest for this profile, although its observed shortfalls are the smallest.')
w('- **Upper-minus-center contrast.** The independent primary slope is the published upper-minus-matched-center '
  'contrast (supplement Table S19, +4.7 pp [-3.7, 13.7] at 99.167%) divided by its log-baseline span.')
w('')
cons = EL['consistency_with_common_baseline']
w('**Consistency check with the common-baseline analysis.** Applying the primary slope from the center to 3.0 '
  'reproduces the common-baseline F* to within '
  f'{max(abs(v["difference_pp"]) for v in cons.values()):.3f} pp '
  f'({", ".join(f"{PNAME[p].split()[0].lower()} {n(cons[p]["F_at_3_arithmetic_pct"])} vs {n(cons[p]["common_baseline_F_star_pct"])}" for p in PROFS)}). '
  'The small residual is the difference between interpolating F linearly and interpolating the control and '
  'policy log rates separately.')
w('')

w('### Implied bias from each observed baseline shortfall')
w('')
w('Two local estimators are shown, with joint 95% intervals (slope and observed baseline both resampled). '
  'Their point estimates agree whenever the observed baseline lies above the sensitivity center; their intervals '
  'can differ slightly because the center itself is resampled.')
w('')
w('- **Piecewise** uses the low-to-center slope for any part of the span below the sensitivity center and the '
  'center-to-high slope above it. It is the most local reading.')
w('- **Primary** uses the center-to-high slope throughout.')
w('')
w('"Beta-only" intervals fix the shortfall at its point value. "Range" is the spread of point estimates across '
  'all estimators in the slope table.')
w('')
rows = []
for wf in ('original', 'alternative'):
    for p in PROFS:
        b = BI[wf][p]
        pw, pr = b['by_slope']['piecewise'], b['by_slope']['primary']
        rows.append([WNAME[wf] if p == 'independent' else '', PNAME[p],
                     f'{n(b["control_rate"], 3)}', f'{n(b["shortfall_pct"], 1)}%',
                     f'{n(pw["bias_pp"], 2, True)} {iv(pw["bias_pp_ci95_joint"])}',
                     f'{n(pr["bias_pp"], 2, True)} {iv(pr["bias_pp_ci95_joint"])}',
                     f'{iv(pr["bias_pp_ci95_beta_only"])}',
                     f'{n(b["point_range_across_slope_estimators_pp"][0], 2, True)} to '
                     f'{n(b["point_range_across_slope_estimators_pp"][1], 2, True)}',
                     f'{n(pw["F_observed_pct"], 1)}', f'{n(pw["F_adjusted_to_3_pct"], 1)} {iv(pw["F_adjusted_to_3_pct_ci95_joint"], 1)}'])
w(table(['Workflow', 'Profile', 'Observed B', 'vs 3.0', 'Bias, piecewise [joint 95%]', 'Bias, primary [joint 95%]',
         'Primary, beta-only 95%', 'Range, all estimators', 'F observed %', 'F at 3.0 (piecewise) % [joint 95%]'],
        rows, ['l', 'l', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']))
w('')
w('In words (piecewise estimator, joint 95% interval):')
w('')
for wf in ('original', 'alternative'):
    for p in PROFS:
        b = BI[wf][p]
        pw = b['by_slope']['piecewise']
        direction = 'higher' if pw['bias_pp'] > 0 else 'lower'
        w(f'- **{WNAME[wf]}, {PNAME[p].lower()}.** The observed baseline {n(b["control_rate"], 2)} is '
          f'{abs(b["shortfall_pct"]):.1f}% below 3.0. The local slope implies the observed F '
          f'({n(pw["F_observed_pct"], 1)}%) is about {abs(pw["bias_pp"]):.2f} pp {direction} than at 3.0 '
          f'(interval {iv(pw["bias_pp_ci95_joint"])} pp'
          + (', which includes zero).' if pw['bias_pp_ci95_joint']['lo'] < 0 < pw['bias_pp_ci95_joint']['hi'] else ').')
          + (' This uses an original-workflow slope (transfer assumption).' if wf == 'alternative' else ''))
w('')
w('**Why "F at 3.0" here differs from the common-baseline F\\*.** This table adjusts the full 512-seed holdout '
  'F. The common-baseline analysis interpolated from the 256-seed sensitivity or matched-center conditions. '
  'Those are different seed samples: center F ' + ', '.join(n(c['F_pct'], 1) for c in (C['ind:center_matched256'],
  C['sens:clustered_center'], C['sens:imperfect_center'])) + '% against holdout F ' + ', '.join(
  n(C['hold:' + p]['F_pct'], 1) for p in PROFS) + '% at the same p. The two adjusted values therefore differ by '
  'that sampling difference, not by the adjustment.')
w('')
w('### Implied bias in the prespecified contrasts')
w('')
w('The contrast bias is bias(first) - bias(second). The adjusted contrast is the observed contrast minus that '
  'bias. Observed-contrast intervals are recomputed here so they can be combined draw by draw; they differ '
  'slightly from the archived intervals (stream and stratification). The archived point contrasts are '
  + ', '.join(f'{WNAME[wf].lower()} {"/".join(str(x) for x in v)}' for wf, v in ARCHIVED_CONTRASTS.items()) + ' pp.')
w('')
rows = []
for wf in ('original', 'alternative'):
    for c, r in CB[wf].items():
        a, b = c.split('-')
        rows.append([WNAME[wf] if c.startswith('independent-clustered') else '', f'{a} - {b}',
                     f'{n(r["observed_contrast_pp"], 1)} {iv(r["observed_contrast_ci98333"], 1)}',
                     f'{n(r["piecewise"]["contrast_bias_pp"], 2, True)} {iv(r["piecewise"]["contrast_bias_ci95_joint"])}',
                     f'{n(r["primary"]["contrast_bias_pp"], 2, True)} {iv(r["primary"]["contrast_bias_ci95_joint"])}',
                     f'{n(r["piecewise"]["adjusted_contrast_pp"], 1)} {iv(r["piecewise"]["adjusted_contrast_ci98333_joint"], 1)}'])
w(table(['Workflow', 'Contrast', 'Observed pp [98.333%]', 'Contrast bias, piecewise [joint 95%]',
         'Contrast bias, primary [joint 95%]', 'Adjusted to 3.0, piecewise [joint 98.333%]'], rows,
        ['l', 'l', 'r', 'r', 'r', 'r']))
w('')
w('### Assumptions and limits of section 1')
w('')
w('1. **Local log-linearity.** F is assumed linear in ln B between the observed baseline and 3.0, within the '
  'segment or segments used. The sensitivity design has only three points per profile, and its curvature is '
  'unresolved.')
w('2. **Moving p is the right way to reach 3.0.** The slopes describe what happens to F when p moves the '
  'baseline. A baseline difference arising some other way (seeds, workflow rules, admission mix) need not '
  'carry the same F change.')
w('3. **Transfer to the alternative workflow is an assumption.** No paired p-sensitivity bank exists for the '
  'alternative workflow, so its bias rows reuse original-workflow slopes. The alternative grid supports a '
  'separate check of the baseline-on-p slope only (section 2). It cannot check dF/d ln B.')
w('4. **Monte Carlo uncertainty only.** The intervals cover Monte Carlo sampling in a fixed synthetic model. '
  'They exclude interpolation-model error, structural uncertainty and clinical-input uncertainty.')
w('5. **Not a matched-baseline result.** The adjusted values are descriptive and post hoc. They do not turn any '
  'failed target-band check into a pass.')
w('')

# ------------------------------------------------------------------ 2. selection propagation
w('## 2. Selection-uncertainty propagation')
w('')
w('Each workflow selected p once from baseline acquisitions on a grid. The archived bootstrap intervals for the selected p '
  'are conditional on valid selections, and every bootstrap draw here was valid. The policy intervals in the '
  'manuscript and the verification report do not include this uncertainty. The calculation here is')
w('')
w('> implied F shift = dF/d ln p x ln(p_bound / p_selected).')
w('')
w('dF/d ln p is the primary (center-to-high) slope from section 1. For the alternative workflow it is the '
  'original-workflow slope (transfer assumption); a sensitivity check follows the table.')
w('')
w('The columns are:')
w('')
w('- **Point slope:** the two shifts at the point slope.')
w('- **With slope uncertainty:** the 2.5th percentile of the lower-end shift and the 97.5th percentile of the '
  'upper-end shift across slope draws.')
w('- **Combined half-width:** sqrt(MC^2 + S^2). MC is the archived Monte Carlo F half-width. S is the largest '
  '|shift| over the slope\'s point value and 95% endpoints times the two p-interval endpoints. It treats the two '
  'sources as independent and symmetric, so it is descriptive only.')
w('')
rows = []
for wf in ('original', 'alternative'):
    for p in PROFS:
        r = SE[wf][p]
        rows.append([WNAME[wf] if p == 'independent' else '', PNAME[p],
                     f'{r["p_selected"]:.5f} [{r["p_ci95"][0]:.5f}, {r["p_ci95"][1]:.5f}]',
                     f'{n(r["dlnp_pct"][0], 1)}% / +{n(r["dlnp_pct"][1], 1)}%',
                     f'{n(r["dF_dlnp"], 1)} {iv(r["dF_dlnp_ci95"], 1)}',
                     f'{n(r["implied_F_shift_pp_point"][0], 2, True)} to {n(r["implied_F_shift_pp_point"][1], 2, True)}',
                     f'{n(r["implied_F_shift_pp_with_slope_uncertainty"][0], 2, True)} to {n(r["implied_F_shift_pp_with_slope_uncertainty"][1], 2, True)}',
                     f'{n(r["conditional_MC_half_width_pp"], 2)}',
                     f'{n(r["quadrature_half_width_pp_envelope"], 2)}'])
w(table(['Workflow', 'Profile', 'Selected p [conditional 95%]', 'Interval vs p', 'dF/d ln p (pp) [95%]',
         'F shift, point slope (pp)', 'F shift, with slope uncertainty (pp)', 'MC F half-width (pp)',
         'Combined half-width (pp)'], rows, ['l', 'l', 'r', 'r', 'r', 'r', 'r', 'r', 'r']))
w('')
w('Sources for the selected-p intervals:')
w('')
w(f'- **Original, independent.** Supplement Table S8, pooled 128-seed row. The fallback-inclusive interval is '
  f'identical, with 0 failed draws.')
w(f'- **Original, clustered and imperfect.** Supplement Table S14, pooled64 row, with 0 failed draws. These are '
  f'retrospective intervals from the older 27-condition exploratory grid.')
w(f'- **Alternative.** Verification RESULTS.md section 3 (`selection.json`), with 0 failed draws.')
w('')
w('Two caveats on these sources:')
w('')
w('- **Fragile imperfect bracket.** Imperfect clearing chose the lower bracket in '
  f'{SI["alternative"]["imperfect"]["lower_bracket_draws"][min(SI["alternative"]["imperfect"]["lower_bracket_draws"], key=lambda k: float(k.split(":")[0]))]:,} '
  'of 20,000 selection draws. The archived conditional interval already includes those draws.')
w('- **Different selection procedures.** The original clustered and imperfect intervals come from an exploratory '
  'design (64 seeds) that differs from the one that fixed p for the holdout.')
w('')
w('**Sensitivity for the alternative workflow.** The transfer assumption can be replaced for the baseline-on-p '
  'part only. Combine the original-workflow dF/d ln B with the alternative grid\'s own log-log slope of B on p, '
  'taken over the m1 to m2 bracket that contains each selected p.')
w('')
rows = []
for p in PROFS:
    s = SE['alternative'][p]['sensitivity_beta_B_times_alt_gamma']
    rows.append([PNAME[p], f'{n(s["alt_grid_rates"][0], 2)} to {n(s["alt_grid_rates"][1], 2)}',
                 f'{n(s["alt_grid_gamma"], 2)} {iv(s["alt_grid_gamma_ci95"])}',
                 f'{n(SL[p]["primary"]["gamma"], 2)}',
                 f'{n(s["dF_dlnp"], 1)} {iv(s["dF_dlnp_ci95"], 1)}',
                 f'{n(s["implied_F_shift_pp_point"][0], 2, True)} to {n(s["implied_F_shift_pp_point"][1], 2, True)}'])
w(table(['Profile', 'Alternative grid B, m1 to m2', 'Alternative d ln B/d ln p [95%]',
         'Original d ln B/d ln p (primary)', 'dF/d ln p (pp) [95%]', 'F shift, point slope (pp)'], rows))
w('')
w('Either way, selection uncertainty moves F by well under 1 pp at the point slope. The slope uncertainty, not '
  'the p interval, dominates the wider bound.')
w('')

# ------------------------------------------------------------------ 3. routes
w('## 3. Route labels in both workflows')
w('')
w('Each acquisition carries one passive label. The label marks the pickup that began the infecting hand\'s '
  'current contamination episode, with the number of hygiene opportunities between that pickup and the '
  'transfer:')
w('')
w('- **Patient pickup:** `patient_pickup`.')
w('- **Own-zone site, current stay:** `own_zone_environment_pickup`, split into 0 and more than 0 intervening '
  'opportunities.')
w('- **Other-zone site:** `other_zone_environment_pickup`.')
w('- **Pre-admission carry-over:** `own_zone_environment_pickup_after_readmission`. The own bedside site\'s '
  'contamination was deposited before the current admission.')
w('')
w('Counts are pooled over 512 paired seeds per profile and arm (control .90, policy .95):')
w('')
w('- **Original:** the 8 October holdout records.')
w('- **Alternative:** the collected policy bank (vacancy cleaning plus empty-bed-only isolation).')
w('')
rr_gap = max(abs(L['reduction_count_pct'] - L['reduction_rate_pct']) for wf in RO for p in PROFS
             for L in RO[wf][p]['labels'].values() if L['reduction_count_pct'] is not None)
w('"Reduction" is 1 - policy/control count, with a paired whole-seed 95% interval. Patient-days differ by less '
  f'than {100 * max(abs(RO[wf][p]["patient_days_relative_difference"]) for wf in RO for p in PROFS):.3f}% between '
  f'arms, so count and rate reductions agree to within {rr_gap:.2f} pp. "Contribution" is (control rate - policy rate of the '
  f'label) / total control rate. It sums to F across the four labels.')
w('')
w('Cells marked † have fewer than 30 control acquisitions. Their reductions are not interpretable, and '
  'intervals at 0 or 100 are degenerate.')
w('')
w('### Side by side (control and policy, both workflows)')
w('')
sbs = (('patient', 'Patient pickup'), ('own_current_0', 'Own-zone, current stay, 0 opp.'),
       ('own_current_pos', 'Own-zone, current stay, >0 opp.'), ('other_zone', 'Other-zone site'),
       ('carry_over', 'Carry-over (pre-admission)'), ('all', 'All'))
rows = []
for p in PROFS:
    for i, (g, name) in enumerate(sbs):
        cells = [PNAME[p] if i == 0 else '', name]
        for wf in ('original', 'alternative'):
            x = RO[wf][p]['labels'][g]
            dag = ' †' if x['small_count_flag'] else ''
            red = 'n/a' if x['reduction_count_pct'] is None else f'{n(x["reduction_count_pct"], 0)}%{dag}'
            cells += [f'{x["control"]:,} ({n(x["control_share_pct"], 1)}%)',
                      f'{x["policy"]:,} ({n(x["policy_share_pct"], 1)}%)', red]
        rows.append(cells)
w(table(['Profile', 'Label', 'Original control', 'Original policy', 'Original reduction',
         'Alternative control', 'Alternative policy', 'Alternative reduction'], rows,
        ['l', 'l', 'r', 'r', 'r', 'r', 'r', 'r']))
w('')
w('Full detail with intervals, the carry-over split by opportunities, and derived rows follows.')
w('')
label_rows = (('patient', 'Patient pickup'), ('own_current_0', 'Own-zone site, current stay, 0 opp.'),
              ('own_current_pos', 'Own-zone site, current stay, >0 opp.'), ('other_zone', 'Other-zone site'),
              ('carry_over', 'Pre-admission carry-over (all)'), ('carry_over_0', '  of which 0 opp.'),
              ('carry_over_pos', '  of which >0 opp.'), ('all', 'All acquisitions'),
              ('all_except_carry', 'All except carry-over'),
              ('paper_metric_site_0', 'Site pickup, 0 opp. (manuscript metric = own 0 + carry-over 0)'))
for wf in ('original', 'alternative'):
    w(f'### {WNAME[wf]} workflow')
    w('')
    rows = []
    for p in PROFS:
        L = RO[wf][p]['labels']
        for i, (g, name) in enumerate(label_rows):
            x = L[g]
            dag = ' †' if x['small_count_flag'] else ''
            red = 'n/a' if x['reduction_count_pct'] is None else f'{n(x["reduction_count_pct"], 1)}%{dag}'
            redci = '' if x['reduction_count_pct_ci95'] is None else ' ' + iv(x['reduction_count_pct_ci95'], 1)
            rows.append([PNAME[p] if i == 0 else '', name,
                         f'{x["control"]:,} ({n(x["control_share_pct"], 1)}%)',
                         f'{x["policy"]:,} ({n(x["policy_share_pct"], 1)}%)',
                         red + redci, n(x['contribution_to_F_pp'], 1)])
    w(table(['Profile', 'Label', 'Control n (share)', 'Policy n (share)', 'Reduction [95%]',
             'Contribution to F (pp)'], rows, ['l', 'l', 'r', 'r', 'r', 'r']))
    w('')
oi, ai = RO['original']['independent']['labels'], RO['alternative']['independent']['labels']
w('### What the table shows (descriptive)')
w('')
w(f'- **Original workflow, independent profile.** Carry-over fell {n(oi["carry_over"]["reduction_count_pct"], 1)}% '
  f'{iv(oi["carry_over"]["reduction_count_pct_ci95"], 1)} under the policy. Patient, current-stay own-zone and '
  f'other-zone labels fell {n(oi["patient"]["reduction_count_pct"], 1)}%, '
  f'{n(oi["own_current_all"]["reduction_count_pct"], 1)}% and {n(oi["other_zone"]["reduction_count_pct"], 1)}%. '
  f'The carry-over share therefore rose from {n(oi["carry_over"]["control_share_pct"], 1)}% to '
  f'{n(oi["carry_over"]["policy_share_pct"], 1)}%. Excluding carry-over, acquisitions fell '
  f'{n(oi["all_except_carry"]["reduction_count_pct"], 1)}%; all labels together fell '
  f'{n(oi["all"]["reduction_count_pct"], 1)}%. This is bookkeeping, not a counterfactual, because each '
  f'acquisition also seeds later ones.')
w(f'- **Alternative workflow.** Carry-over control counts were '
  + ', '.join(f'{RO["alternative"][p]["labels"]["carry_over"]["control"]}' for p in PROFS)
  + ' (independent, clustered, imperfect), against '
  + ', '.join(f'{RO["original"][p]["labels"]["carry_over"]["control"]}' for p in PROFS)
  + f' in the original. The current-stay own-zone share was {n(ai["own_current_all"]["control_share_pct"], 1)}% '
  f'(original {n(oi["own_current_all"]["control_share_pct"], 1)}%) for the independent profile. It fell '
  f'{n(ai["own_current_all"]["reduction_count_pct"], 1)}% under the policy.')
own0_red = [RO[wf][p]['labels']['own_current_0']['reduction_count_pct'] for wf in RO for p in PROFS]
gap = [(abs(RO[wf][p]['labels']['own_current_0']['reduction_count_pct']
            - RO[wf][p]['labels']['patient']['reduction_count_pct']), wf, p) for wf in RO for p in PROFS]
gmax = max(gap)
w('- **Zero-opportunity current-stay site pickups respond to hygiene.** Current-stay own-zone pickups with more '
  'than zero opportunities number 0 to '
  f'{max(max(RO[wf][p]["labels"]["own_current_pos"]["control"], RO[wf][p]["labels"]["own_current_pos"]["policy"]) for wf in RO for p in PROFS)} '
  f'per arm. The zero-opportunity current-stay label fell {min(own0_red):.1f}% to {max(own0_red):.1f}%. In each '
  f'profile and workflow this was within {gmax[0]:.1f} pp of the patient-pickup label (largest gap: '
  f'{gmax[1]}, {gmax[2]}), with overlapping intervals. A zero-opportunity final step does not by itself make a '
  'route insensitive to hygiene.')
w('- **Different conditions.** The original and alternative counts come from different p, seeds and studies. '
  'No between-workflow difference here has a declared estimand or interval. The design cannot attribute any '
  'difference to vacancy cleaning or to the isolation rule separately.')
w('- **Labels are not causal.** They describe realized paths, not causal route fractions.')
w('- **Swap-inherited contamination is mislabelled.** In the original workflow, a patient swapped into a '
  'carrier\'s former bed was not admitted to that bed. Its pickups of the carrier\'s contamination are therefore '
  'labelled own-zone current stay, not carry-over.')
w('')

# ------------------------------------------------------------------ 4. study flow, compliance
w('## 4. Study flow and compliance range')
w('')
w('### Study flow')
w('')
w('Counts are archived records or paired jobs. A paired job holds both branches and is one record.')
w('')
rows = [[r['step'], r['study'], r['question'], r['p_selection'], r['n'],
         r['outcome'].replace('[PENDING: workflow_coefficient_2026-10-21 results not yet available]', '**[PENDING]**'),
         r['source']] for r in J['study_flow']]
w(table(['#', 'Study', 'Question', 'Selection of p', 'n', 'Outcome', 'Source'], rows,
        ['r', 'l', 'l', 'l', 'l', 'l', 'l']))
w('')
w('### Compliance range (rare-transmission limit, supplement Table S5)')
w('')
w('F* = 1 - mean K(h1) / mean K(h0) is the p -> 0 limit of F from 256 transmission-free histories of the '
  'original workflow. Intervals are the archived normal/delta Monte Carlo intervals. q2 is the stationary '
  'two-opportunity comparator, which has no interval. The h-pairs are the three in Table S5 that include all '
  'three profiles.')
w('')
rows = []
for r in CR['by_h_pair']:
    cells = [f'{r["h0"]:.2f} to {r["h1"]:.2f}']
    for p in PROFS:
        x = r[p]
        cells.append(f'{n(x["F_star_pct"])} [{n(x["F_star_ci95"][0])}, {n(x["F_star_ci95"][1])}] (q2 {n(x["q2_response_pct"])})')
    cells += [n(r['independent_minus_clustered_pp'], 1), n(r['independent_minus_imperfect_pp'], 1),
              n(r['clustered_minus_imperfect_pp'], 1), '**[PENDING]**']
    rows.append(cells)
w(table(['Completion h0 to h1', 'Independent F* % [MC 95%] (q2)', 'Clustered F* %', 'Imperfect F* %',
         'Ind - clu (pp)', 'Ind - imp (pp)', 'Clu - imp (pp)', 'Alternative-workflow F*'], rows,
        ['l', 'r', 'r', 'r', 'r', 'r', 'r', 'l']))
w('')
w('K values (mean bridge weights) by h-pair:')
w('')
rows = []
for r in CR['by_h_pair']:
    rows.append([f'{r["h0"]:.2f} to {r["h1"]:.2f}'] + [f'{n(r[p]["K0"])} / {n(r[p]["K1"])}' for p in PROFS])
w(table(['h-pair', 'Independent K(h0) / K(h1)', 'Clustered', 'Imperfect'], rows))
w('')
w('What the compliance-range table shows:')
w('')
ordered = all(r['independent']['F_star_pct'] > r['clustered']['F_star_pct'] > r['imperfect']['F_star_pct']
              for r in CR['by_h_pair'])
w('- **The gaps depend on the completion range.** '
  + ('The ordering independent > clustered > imperfect holds in all three ranges, but the gaps change. '
     if ordered else 'The ordering of the profiles changes across ranges. ')
  + 'Clustered minus imperfect is '
  + ', '.join(n(r['clustered_minus_imperfect_pp'], 1) for r in CR['by_h_pair'])
  + ' pp and independent minus clustered is '
  + ', '.join(n(r['independent_minus_clustered_pp'], 1) for r in CR['by_h_pair'])
  + ' pp, for the three ranges in table order.')
w('- **Not finite-p results.** These are limits, not finite-p holdout results, and they do not predict '
  'absolute rates.')
w('- **Pending.** The alternative-workflow column waits on `workflow_coefficient_2026-10-21`.')
w('')

# ------------------------------------------------------------------ placeholders + provenance
w('## Placeholders')
w('')
w('- **[PENDING] Alternative-workflow p = 0 coefficient.** This covers alternative-workflow F* and its paired '
  'difference from the original on the same backgrounds. Source: '
  '`abm/research/workflow_coefficient_2026-10-21/` (protocol only; no results yet). It fills study-flow step 13 '
  'and the last column of the compliance table.')
w('')
w('## Provenance and reproduction')
w('')
w('Hash checks run by the analysis script:')
w('')
rows = [[h['label'], f'{h["n_files"]:,}', h['mismatches']] for h in J['hash_checks']]
w(table(['Hash check', 'Files', 'Mismatches'], rows))
w('')
w(f'- **Reproduction checks.** {len(J["reproduction_checks"])} checks against archived reports, the common-baseline '
  f'JSON, the alternative `analysis.json`, the supplement displays (Tables S8, S14) and the verification route '
  f'counts. There were {J["reproduction_failures"]} failures.')
w(f'- **Runtime.** Python {J["runtime"]["python"]}, NumPy {J["runtime"]["numpy"]} (local macOS). Analysis '
  f'script SHA256 `{J["script_sha256"][:16]}…`.')
w('- **To reproduce.** Run `.venv/bin/python addenda_analysis.py` (about 1 to 3 minutes; it rewrites '
  'addenda.json), then `.venv/bin/python render_addenda_md.py`. The archived inputs must be present at the '
  'paths recorded in addenda.json `inputs`.')
w('- **No side effects.** No simulation, cluster job, manuscript edit or commit was made. The manuscript and '
  'supplement were read, not edited.')
w('')

(HERE / 'ADDENDA.md').write_text('\n'.join(md) + '\n')
print('wrote', HERE / 'ADDENDA.md', len(md), 'lines')
