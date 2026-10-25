"""Generate LaTeX table fragments and a number-source map for supplement_v2.tex.

Read-only with respect to every source. Writes only into this folder:
  tab_*.tex fragments (\\input by supplement_v2.tex) and NUMBERS_SOURCE_supplement.md.
Asserts headline values against the reviewed RESULTS.md files so that a wrong
field path fails loudly.
"""
import json, math, hashlib, pathlib
from fractions import Fraction
from decimal import Decimal, getcontext

HERE = pathlib.Path(__file__).resolve().parent
ABM = pathlib.Path(__file__).resolve().parents[3]
ROOT = pathlib.Path(__file__).resolve().parents[4]
SRC = {
    'coef': ABM / 'research/workflow_coefficient_2026-10-21/summary.json',
    'coef_gates': ABM / 'research/workflow_coefficient_2026-10-21/gate_checks/GATES.json',
    'coef_bench': ABM / 'research/workflow_coefficient_2026-10-21/BENCHMARK.json',
    'hold': ABM / 'research/full_icu_confirmation_2026-10-08/holdout_report/report.json',
    'alt': ROOT / 'research/robustness_uq_2026-10-17/collection_2026-10-21/execution01/outputs/analysis.json',
    'add': ABM / 'research/manuscript_v2_2026-10-21/addenda/addenda.json',
    'mathv': ABM / 'research/mathematical_extension_2026-10-17/verification.json',
    'rev1': ABM / 'research/manuscript_v2_2026-10-21/reviews_round1/analysis/round1_numbers.json',
}
J = {k: json.load(open(v)) for k, v in SRC.items()}
SHA = {k: hashlib.sha256(open(v, 'rb').read()).hexdigest() for k, v in SRC.items()}
NUM = []  # (where, text, source key, field)


def note(where, text, src, field):
    NUM.append((where, text, src, field))


def r(x, d=1):
    """Round half away from zero; LaTeX minus for negatives."""
    q = 10 ** d
    v = math.copysign(math.floor(abs(x) * q + 0.5) / q, x)
    if v == 0:
        v = 0.0
    s = f"{abs(v):.{d}f}"
    return ("\\mn " + s) if v < 0 else s


def rs(x, d=1):
    """Signed with explicit plus."""
    s = r(x, d)
    return s if s.startswith('\\mn') or float(s) == 0 else '+' + s


def ci(lo, hi, d=1, scale=1.0):
    return f"[{r(lo*scale, d)}, {r(hi*scale, d)}]"


def write(name, body):
    (HERE / name).write_text(body)


# ---------------------------------------------------------------------------
# 1. Workflow coefficient study (summary.json)
# ---------------------------------------------------------------------------
C = J['coef']
P = C['primary']['paired']
S = C['surface']['pairs']
PROF = [('iid_perfect', 'Independent misses', 'independent'),
        ('persistent_perfect', 'Clustered misses', 'clustered'),
        ('iid_imperfect', 'Imperfect clearing', 'imperfect')]

# asserts against RESULTS.md
alt_pct = {k: 100 * P[k]['alternative']['limiting_response'] for k, _, _ in PROF}
org_pct = {k: 100 * P[k]['original']['limiting_response'] for k, _, _ in PROF}
assert [round(alt_pct[k], 2) for k, _, _ in PROF] == [75.18, 53.65, 27.99], alt_pct
assert [round(org_pct[k], 2) for k, _, _ in PROF] == [75.13, 53.65, 27.99], org_pct
assert C['status'] == 'complete' and C['self_checks_pass'] and C['n_alternative_records'] == 256
assert C['failed_records'] == [] and C['missing_planned_seeds'] == []

# 1a. primary coefficient table, both workflows
rows = []
for k, lab, _ in PROF:
    for wf, wlab in (('original', 'Original'), ('alternative', 'Alternative')):
        b = P[k][wf]
        k0, k1 = b['K0'], b['K1']
        fr = b['limiting_response'] * 100
        lo, hi = [100 * x for x in b['limiting_response_normal_95_interval']]
        blo, bhi = [100 * x for x in b['F_bootstrap_95']['interval']]
        rows.append(f"{lab if wf=='original' else ''} & {wlab} & {k0['mean']:.2f} [{k0['normal_approximation_95_interval'][0]:.1f}, {k0['normal_approximation_95_interval'][1]:.1f}] & "
                    f"{k1['mean']:.2f} [{k1['normal_approximation_95_interval'][0]:.1f}, {k1['normal_approximation_95_interval'][1]:.1f}] & "
                    f"{fr:.2f} [{lo:.2f}, {hi:.2f}] & [{blo:.2f}, {bhi:.2f}]\\\\")
        note('Table coefficient-primary', f"{lab} {wlab}: K0 {k0['mean']:.2f}, K1 {k1['mean']:.2f}, F* {fr:.2f}", 'coef', f"primary.paired.{k}.{wf}.(K0,K1,limiting_response,limiting_response_normal_95_interval,F_bootstrap_95)")
write('tab_coef_primary.tex', '\n'.join(rows) + '\n')

# 1b. 22-row surface: coefficients and paired differences
order = []
for key, v in S.items():
    pr = v['parameters']
    order.append((key, pr))
# order like old Table S5: h-pairs .90->.95 first by rho, e; then .80->.90, .95->.99
def sortkey(item):
    key, pr = item
    h0 = pr['h0']
    grp = {0.9: 0, 0.8: 1, 0.95: 2}[h0]
    return (grp, pr['rho'], -pr['e'])
order.sort(key=sortkey)

rows_a, rows_b = [], []
dmax = (None, -1)
n_small = 0
for key, pr in order:
    v = S[key]
    o, a = v['original'], v['alternative']
    fo, fa = 100 * o['limiting_response'], 100 * a['limiting_response']
    d = v['delta_F_alt_minus_orig']
    de = 100 * d['estimate']
    dlo, dhi = [100 * x for x in d['bootstrap_95']['interval']]
    kr = v['delta_K0_alt_minus_orig']
    krel = 100 * kr['relative_change_of_mean']
    klo, khi = [100 * x for x in kr['relative_change_bootstrap_95']['interval']]
    q2 = 100 * v['stationary_q2_response_for_comparison_only']
    # independent check of the q2 comparator
    h0, h1, rho, e = pr['h0'], pr['h1'], pr['rho'], pr['e']
    q2f = lambda h: (1 - h * e) ** 2 + e * e * rho * h * (1 - h)
    assert abs(100 * (1 - q2f(h1) / q2f(h0)) - q2) < 1e-9, key
    assert abs(fa - fo - de) < 1e-9
    if abs(de) > dmax[1]:
        dmax = (key, abs(de))
    if abs(de) < 0.06:
        n_small += 1
    hp = f"{h0:.2f}$\\to${h1:.2f}"
    rows_a.append(f"{hp} & {rho:.2f} & {e:.1f} & {o['K0']['mean']:.2f} & {o['K1']['mean']:.2f} & {fo:.2f} & {a['K0']['mean']:.2f} & {a['K1']['mean']:.2f} & {fa:.2f}\\\\")
    geo = J['rev1']['renewal_comparator']['rows'][key]
    assert abs(geo['icu_original_pct'] - fo) < 1e-9 and abs(geo['F_q2_pct'] - q2) < 1e-9, key
    rows_b.append(f"{hp} & {rho:.2f} & {e:.1f} & {fo:.2f} & {fa:.2f} & {rs(de,3)} [{rs(dlo,3)}, {rs(dhi,3)}] & {r(krel,1)} [{r(klo,1)}, {r(khi,1)}] & {q2:.2f} & {geo['F_geo_pct']:.2f}\\\\")
    note('Table surface22 b (F_geo)', f"{key}: F_geo {geo['F_geo_pct']:.2f}, error vs original {geo['err_geo_original_pp']:+.2f} pp", 'rev1', f"renewal_comparator.rows.{key}")
    note('Tables surface22 a/b', f"{key}: F* orig {fo:.2f}, alt {fa:.2f}, dF {de:+.3f} [{dlo:+.3f},{dhi:+.3f}], K0 rel {krel:.1f} [{klo:.1f},{khi:.1f}], q2 {q2:.2f}", 'coef', f"surface.pairs.{key}.(original|alternative).limiting_response; delta_F_alt_minus_orig.bootstrap_95; delta_K0_alt_minus_orig.relative_change_*; stationary_q2_response_for_comparison_only")
assert len(rows_a) == 22
assert n_small == 21, n_small
assert dmax[0] == 'h0.95_to_0.99_rho0_e1' and round(dmax[1], 3) == 0.272, dmax
write('tab_surface22_coef.tex', '\n'.join(rows_a) + '\n')
write('tab_surface22_diff.tex', '\n'.join(rows_b) + '\n')

# 1c. opportunity-count mix (primary profiles, both arms)
rows = []
for k, lab, _ in PROF:
    mix = P[k]['K_by_opportunities_mean']
    for arm, alab in (('baseline', 'control (90\\%)'), ('policy', 'policy (95\\%)')):
        cells = []
        for wf in ('original', 'alternative'):
            m = mix[wf][arm]
            tot = sum(m.values())
            assert m['1'] == 0.0
            sh = {kk: 100 * vv / tot for kk, vv in m.items()}
            cells.append(f"{sh['0']:.2f} & {sh['2']:.2f} & {sh['>=3']:.2f}")
            note('Table kmix', f"{lab} {arm} {wf}: k0 {sh['0']:.2f}%, k2 {sh['2']:.2f}%, k>=3 {sh['>=3']:.2f}%", 'coef', f"primary.paired.{k}.K_by_opportunities_mean.{wf}.{arm}")
        rows.append(f"{lab if arm=='baseline' else ''} & {alab} & {cells[0]} & {cells[1]}\\\\")
write('tab_kmix.tex', '\n'.join(rows) + '\n')
# check RESULTS.md k-mix figures
m = P['iid_perfect']['K_by_opportunities_mean']
tot = sum(m['original']['baseline'].values())
assert round(100 * m['original']['baseline']['2'] / tot, 2) == 98.95
tot = sum(m['original']['policy'].values())
assert round(100 * m['original']['policy']['0'] / tot, 2) == 0.28

# 1d. prespecified prediction checks
pr1 = C['predictions']['prediction_1_small_change']['by_profile']
pr2 = C['predictions']['prediction_2_sign_reversal']['by_workflow']
pr3 = C['predictions']['prediction_3_crossing']
assert C['predictions']['prediction_1_small_change']['holds_all'] and C['predictions']['prediction_2_sign_reversal']['holds'] and pr3['holds']
assert round(100 * pr2['alternative']['D_e1'], 2) == -21.53 and round(100 * pr2['original_same_seeds']['D_e1'], 2) == -21.48
assert round(100 * pr2['alternative']['D_e0p8'], 2) == 9.02 and round(100 * pr2['original_same_seeds']['D_e0p8'], 2) == 9.02
assert round(pr3['by_workflow']['alternative']['e_star'], 4) == 0.9187 and round(pr3['by_workflow']['original_same_seeds']['e_star'], 4) == 0.9187
rows = []
for k, lab, _ in PROF:
    b = pr1[k]
    lo, hi = b['paired_bootstrap_95_pp']
    rows.append(f"1. Primary limit changes by $<2$ pp: {lab.lower()} & {b['published_original_percent']:.1f} & {b['F_alt_percent']:.2f} & {rs(b['F_alt_minus_exact_original_pp'],3)} [{rs(lo,3)}, {rs(hi,3)}] & Held\\\\")
    note('Table predictions', f"{lab}: alt {b['F_alt_percent']:.2f}, diff vs same-seed {b['F_alt_minus_exact_original_pp']:+.3f}, vs published {b['F_alt_minus_published_pp']:+.3f}", 'coef', f"predictions.prediction_1_small_change.by_profile.{k}")
for wf, wlab in (('original_same_seeds', 'original'), ('alternative', 'alternative')):
    b = pr2[wf]
    d1, d8 = 100 * b['D_e1'], 100 * b['D_e0p8']
    i1 = [100 * x for x in b['D_e1_bootstrap_95']['interval']]
    i8 = [100 * x for x in b['D_e0p8_bootstrap_95']['interval']]
    note('Table predictions', f"{wlab}: D(1) {d1:.2f} [{i1[0]:.2f},{i1[1]:.2f}], D(.8) {d8:.2f} [{i8[0]:.2f},{i8[1]:.2f}]", 'coef', f"predictions.prediction_2_sign_reversal.by_workflow.{wf}")
for wf in ('original_same_seeds', 'alternative'):
    b = pr3['by_workflow'][wf]
    note('Table predictions', f"{wf}: e* {b['e_star']:.4f} [{b['e_star_bootstrap_95']['interval'][0]:.4f},{b['e_star_bootstrap_95']['interval'][1]:.4f}], distance {b['distance_from_reference']:.4f}", 'coef', f"predictions.prediction_3_crossing.by_workflow.{wf}")
b2o, b2a = pr2['original_same_seeds'], pr2['alternative']
rows.append(f"2. Sign reversal, $\\Delta(1)=F^*(\\rho{{=}}.5,1)-F^*(0,1)$ & --- & {r(100*b2o['D_e1'],2)} / {r(100*b2a['D_e1'],2)} & original / alternative; all 20,000 draws $<0$ & Held\\\\")
rows.append(f"2. Sign reversal, $\\Delta(.8)=F^*(\\rho{{=}}.5,.8)-F^*(0,.8)$ & --- & {rs(100*b2o['D_e0p8'],2)} / {rs(100*b2a['D_e0p8'],2)} & original / alternative; all 20,000 draws $>0$ & Held\\\\")
e_o, e_a = pr3['by_workflow']['original_same_seeds'], pr3['by_workflow']['alternative']
rows.append(f"3. Crossing $e^*$ (protocol two-point statistic, $e=.9$ to 1) within .05 of .9235 & .9235 & {e_o['e_star']:.4f} / {e_a['e_star']:.4f} & distance {e_o['distance_from_reference']:.4f} / {e_a['distance_from_reference']:.4f} & Held\\\\")
write('tab_predictions.tex', '\n'.join(rows) + '\n')

# 1e. post hoc finite-p comparison
ph = C['post_hoc_finite_p']['by_profile']
rows = []
for k, lab, _ in PROF:
    b = ph[k]
    flo, fhi = [100 * x for x in b['finite_p_F_95']]
    la = b['leading_order_alternative_K']
    lo_ = b['leading_order_original_K']
    rows.append(f"{lab} & {b['selected_p']:.5f} & {100*b['finite_p_F']:.1f} [{flo:.1f}, {fhi:.1f}] & {100*b['limit_alternative']:.2f} & {rs(b['finite_p_F_minus_limit_alternative_pp'],2)} & {la['control_rate']:.2f} / {b['measured_control_rate']:.2f} ({la['control_ratio_to_measured']:.2f}) & {la['policy_rate']:.2f} / {b['measured_policy_rate']:.2f} ({la['policy_ratio_to_measured']:.2f})\\\\")
    assert b['limit_alternative_inside_finite_p_95']
    assert abs(b['finite_p_F_minus_limit_alternative_pp']) < 1.4
    note('Table posthoc', f"{lab}: finite F {100*b['finite_p_F']:.1f}, alt limit {100*b['limit_alternative']:.2f}, diff {b['finite_p_F_minus_limit_alternative_pp']:+.2f}, LO ctrl {la['control_rate']:.2f}/{b['measured_control_rate']:.2f} ratio {la['control_ratio_to_measured']:.2f}", 'coef', f"post_hoc_finite_p.by_profile.{k}")
write('tab_posthoc.tex', '\n'.join(rows) + '\n')
ratios = [ph[k]['leading_order_alternative_K']['control_ratio_to_measured'] for k, _, _ in PROF]
assert [round(x, 2) for x in ratios] == [0.31, 0.59, 0.60], ratios

# 1f. compliance range, both workflows
hpairs = [(0.8, 0.9), (0.9, 0.95), (0.95, 0.99)]
pk = {('iid_perfect'): (0, 1), ('persistent_perfect'): (0.5, 1), ('iid_imperfect'): (0, 0.8)}
def find(h0, rho, e):
    for key, pr in order:
        if pr['h0'] == h0 and pr['rho'] == rho and pr['e'] == e:
            return S[key]
    raise KeyError
rows = []
for h0, h1 in hpairs:
    for wf, wlab in (('original', 'Original'), ('alternative', 'Alternative')):
        fs = [100 * find(h0, *pk[k])[wf]['limiting_response'] for k, _, _ in PROF]
        q2 = [100 * find(h0, *pk[k])['stationary_q2_response_for_comparison_only'] for k, _, _ in PROF]
        g = (fs[0] - fs[1], fs[0] - fs[2], fs[1] - fs[2])
        q2txt = f" (q$_2$ {q2[0]:.2f} / {q2[1]:.2f} / {q2[2]:.2f})" if wf == 'original' else ''
        rows.append(f"{f'{h0:.2f}$\\to${h1:.2f}' if wf=='original' else ''} & {wlab} & {fs[0]:.2f} & {fs[1]:.2f} & {fs[2]:.2f} & {g[0]:.1f} & {g[1]:.1f} & {g[2]:.1f}\\\\")
        note('Table compliance', f"{h0}->{h1} {wf}: {fs[0]:.2f}/{fs[1]:.2f}/{fs[2]:.2f}; gaps {g[0]:.1f}/{g[1]:.1f}/{g[2]:.1f}", 'coef', 'surface.pairs.<h-pair,profile>.<wf>.limiting_response')
write('tab_compliance.tex', '\n'.join(rows) + '\n')

# 1g. gates / benchmark / diagnostics facts
G = J['coef_gates']['gate1']
assert G['all_passed'] and len(G['entries']) == 8
sites = [(e['seed'], e['bank'], e['site'], 'surface' in e) for e in G['entries']]
assert sum(1 for s in sites if s[2] == 'cluster') == 4 and sum(1 for s in sites if s[1] == 'A') == 4
assert sum(1 for s in sites if s[3]) == 2
B = J['coef_bench']
assert B['decision']['surface'] == 'full'
note('S4 gates', f"Gate 1 seeds {sites}", 'coef_gates', 'gate1.entries[*].(seed,bank,site,surface)')
note('S4 budget', f"projected {B['projection']['projected_cpu_hours_used_for_decision']:.1f} CPU-h vs cap {B['cap_cpu_hours']}", 'coef_bench', 'projection.projected_cpu_hours_used_for_decision; cap_cpu_hours')
wd = C['workflow_diagnostics']
note('S4 diagnostics', f"relocations per obs window {wd['per_record_observation_mean']['workflow_events.isolation_relocation']:.2f}; deferral patient-days {wd['deferral_patient_days_per_observation_window_mean']:.2f}; swaps 0; carrier contacts in shared beds share {wd['observation_carrier_contacts_share_in_shared_beds']:.4f}", 'coef', 'workflow_diagnostics')
ex = C['exposure']['observation_patient_days']
note('S4 exposure', f"patient-days orig {ex['original']['mean']:.1f}, alt {ex['alternative']['mean']:.1f}, paired {ex['paired_alt_minus_orig']['mean']:.2f} {ex['paired_alt_minus_orig']['normal_approximation_95_interval']}", 'coef', 'exposure.observation_patient_days')
note('S4 history', f"seeds with identical ordered event hash: {C['history_identity']['n_seeds_identical_ordered_event_hash']} of {C['history_identity']['n_seeds']}", 'coef', 'history_identity')
banks = C['primary']['banks']
bmax = max(abs(100 * banks[k][wf]['B_minus_A_F']['estimate']) for k, _, _ in PROF for wf in ('alternative', 'original_same_seeds'))
assert bmax < 0.002, bmax
note('S4 banks', f"max |B-A| F* = {bmax:.4f} pp", 'coef', 'primary.banks.*.*.B_minus_A_F.estimate')
FACTS = dict(
    reloc=wd['per_record_observation_mean']['workflow_events.isolation_relocation'],
    defer=wd['deferral_patient_days_per_observation_window_mean'],
    shared_carrier=100 * wd['observation_carrier_contacts_share_in_shared_beds'],
    pd_orig=ex['original']['mean'], pd_alt=ex['alternative']['mean'],
    pd_diff=ex['paired_alt_minus_orig']['mean'], pd_diff_ci=ex['paired_alt_minus_orig']['normal_approximation_95_interval'],
    bmax=bmax, n_identical=C['history_identity']['n_seeds_identical_ordered_event_hash'],
)

# ---------------------------------------------------------------------------
# 2. Baseline-matching checks (original holdout report; alternative analysis)
# ---------------------------------------------------------------------------
H = J['hold']
bm = H['baseline_matching']
assert bm['established'] is False and bm['valid_calibration'] is True
names = {
    'independent_absolute_baseline': ('Original holdout', 'Independent baseline rate', '[2.7, 3.3]'),
    'persistent_absolute_baseline': ('Original holdout', 'Clustered baseline rate', '[2.7, 3.3]'),
    'imperfect_absolute_baseline': ('Original holdout', 'Imperfect baseline rate', '[2.7, 3.3]'),
    'independent_over_persistent_baseline_ratio': ('Original holdout', 'Independent/clustered rate ratio', '[0.909, 1.1]'),
    'independent_over_imperfect_baseline_ratio': ('Original holdout', 'Independent/imperfect rate ratio', '[0.909, 1.1]'),
    'persistent_over_imperfect_baseline_ratio': ('Original holdout', 'Clustered/imperfect rate ratio', '[0.909, 1.1]'),
}
rows = []
for i, g in enumerate(bm['numerical_gates'], 1):
    st, lab, band = names[g['name']]
    lo, hi = g['uncertainty']['intervals']['0.991666666667']
    assert g['passed'] is False
    dd = 4 if 'ratio' in g['name'] else 3
    band = '[0.9091, 1.1000]' if 'ratio' in g['name'] else band
    rows.append(f"{i} & {st} & {lab} & {g['estimate']:.{dd}f} & [{lo:.{dd}f}, {hi:.{dd}f}] & {band} & Not met\\\\")
    note('Table gates', f"{lab}: {g['estimate']:.3f} [{lo:.3f},{hi:.3f}] (99.167%)", 'hold', f"baseline_matching.numerical_gates[name={g['name']}].uncertainty.intervals['0.991666666667']")
write('tab_gates_original.tex', '\n'.join(rows) + '\n')

A = J['alt']
ALTP = {s['scenario']['id']: s for s in A['scenarios']}
rows = []
for pid, lab in (('independent', 'Independent'), ('clustered', 'Clustered misses'), ('imperfect', 'Imperfect clearing')):
    s = ALTP[pid]
    lo, hi = s['control_rate95']['interval']
    inside = s['descriptive_baseline95_inside_target_band']
    rows.append(f"{lab} & {s['control_rate']:.3f} & [{lo:.3f}, {hi:.3f}] & {'inside' if inside else 'not inside'}\\\\")
    note('Table alt band', f"{lab}: {s['control_rate']:.3f} [{lo:.3f},{hi:.3f}] inside={inside}", 'alt', f"scenarios[id={pid}].control_rate95; descriptive_baseline95_inside_target_band")
write('tab_gates_alt.tex', '\n'.join(rows) + '\n')

# holdout results both workflows
HP = {s['scenario']['id']: s for s in H['scenarios']}
rows = []
for (hid, aid, lab) in (('independent_target_candidate', 'independent', 'Independent misses'), ('persistent_target_candidate', 'clustered', 'Clustered misses'), ('imperfect_target_candidate', 'imperfect', 'Imperfect clearing')):
    s = HP[hid]
    c0 = s['arms']['baseline']; c1 = s['arms']['policy']
    F = 100 * s['rate_prevented_fraction']; Fi = [100 * x for x in s['rate_prevented_fraction_uncertainty']['intervals']['0.95']]
    ab = s['absolute_prevented_rate']; abi = s['absolute_prevented_rate_uncertainty']['intervals']['0.95']
    mp = s['mean_prevented']; mpi = s['mean_prevented_uncertainty']['intervals']['0.95']
    r0i = c0['rate_uncertainty']['intervals']['0.95']; r1i = c1['rate_uncertainty']['intervals']['0.95']
    rows.append(f"{lab} & Original & {s['scenario']['p']:.5f} & {c0['rate_per_1000_patient_days']:.2f} [{r0i[0]:.2f}, {r0i[1]:.2f}] & {c1['rate_per_1000_patient_days']:.2f} [{r1i[0]:.2f}, {r1i[1]:.2f}] & {c0['acquired_total']:,} $\\to$ {c1['acquired_total']:,} & {F:.1f} [{Fi[0]:.1f}, {Fi[1]:.1f}] & {ab:.2f} [{abi[0]:.2f}, {abi[1]:.2f}] & {mp:.2f} [{mpi[0]:.2f}, {mpi[1]:.2f}]\\\\")
    note('Table holdout both', f"{lab} original: p {s['scenario']['p']:.5f}, R0 {c0['rate_per_1000_patient_days']:.2f}, R1 {c1['rate_per_1000_patient_days']:.2f}, F {F:.1f}, abs {ab:.2f}, mean prevented {mp:.2f}", 'hold', f"scenarios[id={hid}]")
    t = ALTP[aid]
    pd = J['add']['routes']['results']['alternative'][aid]['labels']['all']
    rows.append(f" & Alternative & {t['scenario']['p']:.5f} & {t['control_rate']:.2f} [{t['control_rate95']['interval'][0]:.2f}, {t['control_rate95']['interval'][1]:.2f}] & {t['policy_rate']:.2f} [{t['policy_rate95']['interval'][0]:.2f}, {t['policy_rate95']['interval'][1]:.2f}] & {pd['control']:,} $\\to$ {pd['policy']:,} & {100*t['F']:.1f} [{100*t['F95']['interval'][0]:.1f}, {100*t['F95']['interval'][1]:.1f}] & {t['absolute_rate_reduction']:.2f} [{t['absolute_rate_reduction95']['interval'][0]:.2f}, {t['absolute_rate_reduction95']['interval'][1]:.2f}] & {t['mean_carriers_prevented']:.2f} [{t['mean_carriers_prevented95']['interval'][0]:.2f}, {t['mean_carriers_prevented95']['interval'][1]:.2f}]\\\\")
    note('Table holdout both', f"{lab} alternative: p {t['scenario']['p']:.5f}, R0 {t['control_rate']:.2f}, R1 {t['policy_rate']:.2f}, F {100*t['F']:.1f}, abs {t['absolute_rate_reduction']:.2f}, mean prevented {t['mean_carriers_prevented']:.2f}", 'alt', f"scenarios[id={aid}]; counts from addenda routes.results.alternative.{aid}.labels.all")
write('tab_holdout_both.tex', '\n'.join(rows) + '\n')
assert [round(100 * ALTP[p]['F'], 1) for p in ('independent', 'clustered', 'imperfect')] == [76.5, 53.7, 27.3]

rows = []
for c in H['primary_profile_F_contrasts']:
    pass
CN = {'independent_minus_clustered_F': 'Independent $-$ clustered', 'independent_minus_imperfect_F': 'Independent $-$ imperfect', 'clustered_minus_imperfect_F': 'Clustered $-$ imperfect'}
alt_c = {c['name']: c for c in A['primary_contrasts']}
orig_c = {}
nm = {'independent': 'independent', 'persistent': 'clustered', 'imperfect': 'imperfect'}
for c in H['primary_profile_F_contrasts']:
    orig_c[f"{nm[c['first']]}_minus_{nm[c['second']]}_F"] = c
for name, lab in CN.items():
    o = orig_c[name]['F_difference']; a = alt_c[name]
    oi = o['uncertainty']['intervals']['0.983333333333']; ai = a['uncertainty']['interval']
    rows.append(f"{lab} & {100*o['estimate']:.1f} [{100*oi[0]:.1f}, {100*oi[1]:.1f}] & {100*a['estimate']:.1f} [{100*ai[0]:.1f}, {100*ai[1]:.1f}]\\\\")
    note('Table contrasts both', f"{lab}: original {100*o['estimate']:.1f} [{100*oi[0]:.1f},{100*oi[1]:.1f}], alternative {100*a['estimate']:.1f} [{100*ai[0]:.1f},{100*ai[1]:.1f}]", 'hold/alt', 'primary_profile_F_contrasts / primary_contrasts')
write('tab_contrasts_both.tex', '\n'.join(rows) + '\n')

# ---------------------------------------------------------------------------
# 3. Addenda: elasticity, bias, contrasts, selection, routes, study flow
# ---------------------------------------------------------------------------
AD = J['add']
EL = AD['elasticity']['slopes']
labs = AD['elasticity']['slope_labels']
rows = []
for prof, plab, keys in (('independent', 'Independent', ['primary', 'env_low_mid', 'env_mid_high', 'env_ols3']),
                         ('clustered', 'Clustered', ['primary', 'low_center', 'ols3']),
                         ('imperfect', 'Imperfect', ['primary', 'low_center', 'ols3'])):
    first = True
    for k in keys:
        s = EL[prof][k]
        cr = s['control_rates']
        span = '--'.join(f"{x:.2f}" for x in cr)
        short = {'primary': 'Center$\\to$1.25$\\times$ (primary)', 'low_center': '0.75$\\times\\to$center',
                 'ols3': 'Three-point OLS', 'env_low_mid': 'Env.\\ bank .081$\\to$.108',
                 'env_mid_high': 'Env.\\ bank .108$\\to$.135', 'env_ols3': 'Env.\\ bank, OLS'}[k]
        bi = s['beta_B_ci95']; pi = s['beta_B_per_10pct_ci95']; bp = s['beta_p_ci95']; gi = s['gamma_ci95']
        rows.append(f"{plab if first else ''} & {short} & {span} & {r(s['beta_B'],2)} [{r(bi['lo'],2)}, {r(bi['hi'],2)}] & {r(s['beta_B_per_10pct'],2)} [{r(pi['lo'],2)}, {r(pi['hi'],2)}] & {r(s['beta_p'],2)} [{r(bp['lo'],2)}, {r(bp['hi'],2)}] & {s['gamma']:.2f} [{gi['lo']:.2f}, {gi['hi']:.2f}]\\\\")
        note('Table slopes', f"{prof} {k}: beta_B {s['beta_B']:.2f} [{bi['lo']:.2f},{bi['hi']:.2f}]", 'add', f"elasticity.slopes.{prof}.{k}")
        assert bi['lo'] < 0 < bi['hi']
        first = False
write('tab_slopes.tex', '\n'.join(rows) + '\n')

BI = AD['bias']['results']
rows = []
for wf, wlab in (('original', 'Original'), ('alternative', 'Alternative')):
    first = True
    for prof, plab in (('independent', 'Independent misses'), ('clustered', 'Clustered misses'), ('imperfect', 'Imperfect clearing')):
        b = BI[wf][prof]
        pw = b['by_slope']['piecewise']; pm = b['by_slope']['primary']
        allpts = [v['bias_pp'] for v in b['by_slope'].values()]
        rows.append(f"{wlab if first else ''} & {plab} & {b['control_rate']:.3f} & {r(b['shortfall_pct'],1)}\\% & {rs(pw['bias_pp'],2)} [{r(pw['bias_pp_ci95_joint']['lo'],2)}, {r(pw['bias_pp_ci95_joint']['hi'],2)}] & {rs(pm['bias_pp'],2)} [{r(pm['bias_pp_ci95_joint']['lo'],2)}, {r(pm['bias_pp_ci95_joint']['hi'],2)}] & {rs(min(allpts),2)} to {rs(max(allpts),2)} & {pw['F_observed_pct']:.1f} & {pw['F_adjusted_to_3_pct']:.1f} [{pw['F_adjusted_to_3_pct_ci95_joint']['lo']:.1f}, {pw['F_adjusted_to_3_pct_ci95_joint']['hi']:.1f}]\\\\")
        note('Table bias', f"{wf} {prof}: B {b['control_rate']:.3f}, piecewise {pw['bias_pp']:+.2f}, primary {pm['bias_pp']:+.2f}", 'add', f"bias.results.{wf}.{prof}.by_slope.(piecewise|primary)")
        first = False
write('tab_bias.tex', '\n'.join(rows) + '\n')

BC = AD['bias']['contrasts']
rows = []
CL = {'independent-clustered': 'Independent $-$ clustered', 'independent-imperfect': 'Independent $-$ imperfect', 'clustered-imperfect': 'Clustered $-$ imperfect'}
for wf, wlab in (('original', 'Original'), ('alternative', 'Alternative')):
    first = True
    for ck, clab in CL.items():
        c = BC[wf][ck]
        oc = c['observed_contrast_ci98333']; pw = c['piecewise']; pm = c['primary']
        rows.append(f"{wlab if first else ''} & {clab} & {c['observed_contrast_pp']:.1f} [{oc['lo']:.1f}, {oc['hi']:.1f}] & {rs(pw['contrast_bias_pp'],2)} [{r(pw['contrast_bias_ci95_joint']['lo'],2)}, {r(pw['contrast_bias_ci95_joint']['hi'],2)}] & {rs(pm['contrast_bias_pp'],2)} [{r(pm['contrast_bias_ci95_joint']['lo'],2)}, {r(pm['contrast_bias_ci95_joint']['hi'],2)}] & {pw['adjusted_contrast_pp']:.1f} [{pw['adjusted_contrast_ci98333_joint']['lo']:.1f}, {pw['adjusted_contrast_ci98333_joint']['hi']:.1f}]\\\\")
        assert pw['adjusted_contrast_ci98333_joint']['lo'] > 0
        note('Table contrast bias', f"{wf} {ck}: adjusted {pw['adjusted_contrast_pp']:.1f}", 'add', f"bias.contrasts.{wf}.{ck}")
        first = False
write('tab_contrast_bias.tex', '\n'.join(rows) + '\n')

# selection propagation with half-widths from UNROUNDED archived F intervals
SE = AD['selection']['results']
unrounded = {
    ('original', 'independent'): HP['independent_target_candidate']['rate_prevented_fraction_uncertainty']['intervals']['0.95'],
    ('original', 'clustered'): HP['persistent_target_candidate']['rate_prevented_fraction_uncertainty']['intervals']['0.95'],
    ('original', 'imperfect'): HP['imperfect_target_candidate']['rate_prevented_fraction_uncertainty']['intervals']['0.95'],
    ('alternative', 'independent'): ALTP['independent']['F95']['interval'],
    ('alternative', 'clustered'): ALTP['clustered']['F95']['interval'],
    ('alternative', 'imperfect'): ALTP['imperfect']['F95']['interval'],
}
rows = []
HW = {}
for wf, wlab in (('original', 'Original'), ('alternative', 'Alternative')):
    first = True
    for prof, plab in (('independent', 'Independent misses'), ('clustered', 'Clustered misses'), ('imperfect', 'Imperfect clearing')):
        s = SE[wf][prof]
        lo, hi = unrounded[(wf, prof)]
        mc = 100 * (hi - lo) / 2
        Ssel = s['selection_half_width_pp_envelope']
        comb = math.sqrt(mc ** 2 + Ssel ** 2)
        HW[(wf, prof)] = (mc, comb)
        dl = s['dlnp_pct']; sp = s['implied_F_shift_pp_point']; su = s['implied_F_shift_pp_with_slope_uncertainty']; dfi = s['dF_dlnp_ci95']
        rows.append(f"{wlab if first else ''} & {plab.split()[0]} & {s['p_selected']:.5f} \\newline [{s['p_ci95'][0]:.5f}, {s['p_ci95'][1]:.5f}] & {r(dl[0],1)}\\%/+{dl[1]:.1f}\\% & {r(s['dF_dlnp'],1)} \\newline [{r(dfi['lo'],1)}, {r(dfi['hi'],1)}] & {r(sp[0],2)} to +{sp[1]:.2f} & {r(su[0],2)} to +{su[1]:.2f} & {mc:.2f} & {comb:.2f}\\\\")
        note('Table selection propagation', f"{wf} {prof}: MC half-width (unrounded) {mc:.3f}, S {Ssel:.3f}, combined {comb:.3f}; point shift {sp}; slope-uncertainty envelope {su}", 'add+hold/alt', f"selection.results.{wf}.{prof}; MC from unrounded archived F 95% interval")
        first = False
write('tab_selection_prop.tex', '\n'.join(rows) + '\n')
assert round(HW[('alternative', 'clustered')][0], 2) == 3.88 and round(HW[('alternative', 'clustered')][1], 2) == 4.10
assert round(HW[('original', 'imperfect')][0], 2) == 5.24
assert round(HW[('original', 'independent')][1], 2) == 3.73 and round(HW[('original', 'clustered')][1], 2) == 3.95
mcs = [v[0] for v in HW.values()]
assert round(min(mcs), 2) == 2.30 and round(max(mcs), 2) == 5.24

rows = []
for prof, plab in (('independent', 'Independent misses'), ('clustered', 'Clustered misses'), ('imperfect', 'Imperfect clearing')):
    s = SE['alternative'][prof]['sensitivity_beta_B_times_alt_gamma']
    g = AD['selection']['alt_grid_gamma'][prof]
    orig_gamma = EL[prof]['primary']['gamma']
    di = s['dF_dlnp_ci95']; sp = s['implied_F_shift_pp_point']
    rows.append(f"{plab} & {g['B_m1']:.2f} to {g['B_m2']:.2f} & {g['gamma']:.2f} [{g['gamma_ci95']['lo']:.2f}, {g['gamma_ci95']['hi']:.2f}] & {orig_gamma:.2f} & {r(s['dF_dlnp'],1)} [{r(di['lo'],1)}, {r(di['hi'],1)}] & {r(sp[0],2)} to +{sp[1]:.2f}\\\\")
    note('Table selection alt-gamma', f"{prof}: alt gamma {g['gamma']:.2f}", 'add', f"selection.alt_grid_gamma.{prof}; selection.results.alternative.{prof}.sensitivity_beta_B_times_alt_gamma")
write('tab_selection_altgamma.tex', '\n'.join(rows) + '\n')

# route tables
RT = AD['routes']['results']
LBL = [('patient', 'Patient pickup'), ('own_current_0', 'Own-zone site, current stay, 0 opp.'),
       ('own_current_pos', 'Own-zone site, current stay, $>0$ opp.'), ('other_zone', 'Other-zone site'),
       ('carry_over', 'Pre-admission carry-over'), ('all', 'All acquisitions'),
       ('all_except_carry', 'All except carry-over'), ('paper_metric_site_0', 'Site pickup, 0 opp.\\ (own 0 + carry-over 0)')]
UNDEF = []
for wf in ('original', 'alternative'):
    rows = []
    for prof, plab in (('independent', 'Independent misses'), ('clustered', 'Clustered misses'), ('imperfect', 'Imperfect clearing')):
        L = RT[wf][prof]['labels']
        first = True
        for key, lab in LBL:
            x = L[key]
            dag = '$^\\dagger$' if (x['small_count_flag'] and key not in ('all', 'all_except_carry', 'paper_metric_site_0')) else ''
            if x['reduction_count_pct'] is None:
                red = 'n/a'
            else:
                ci_ = x['reduction_count_pct_ci95']
                red = f"{x['reduction_count_pct']:.1f}{dag} [{r(ci_['lo'],1)}, {ci_['hi']:.1f}]"
                if ci_['undefined_draws']:
                    red += f"$^\\ddagger$"
                    UNDEF.append((wf, prof, key, ci_['undefined_draws']))
            rows.append(f"{plab if first else ''} & {lab} & {x['control']:,} ({x['control_share_pct']:.1f}) & {x['policy']:,} ({x['policy_share_pct']:.1f}) & {red} & {r(x['contribution_to_F_pp'],1)}\\\\")
            note(f'Table routes {wf}', f"{prof} {key}: {x['control']} -> {x['policy']}, red {x['reduction_count_pct']}", 'add', f"routes.results.{wf}.{prof}.labels.{key}")
            first = False
        rows.append('\\addlinespace')
    write(f'tab_routes_{wf}.tex', '\n'.join(rows[:-1]) + '\n')
# carry-over >0 opportunities for original clustered (1 -> 0) is not a displayed row; record its undefined draws
co = RT['original']['clustered']['labels']['carry_over_pos']
if co['reduction_count_pct_ci95'] and co['reduction_count_pct_ci95']['undefined_draws']:
    UNDEF.append(('original', 'clustered', 'carry_over_pos', co['reduction_count_pct_ci95']['undefined_draws']))
UNDEF_SET = {(a, b, c): d for a, b, c, d in UNDEF}
assert UNDEF_SET.get(('original', 'clustered', 'own_current_pos')) == 53
assert UNDEF_SET.get(('original', 'clustered', 'carry_over_pos')) == 7417
assert UNDEF_SET.get(('alternative', 'clustered', 'own_current_pos')) == 56
note('Routes undefined draws', str(UNDEF), 'add', 'routes.results.*.*.labels.*.reduction_count_pct_ci95.undefined_draws')

# study flow
SF = AD['study_flow']
rows = []
for st in SF:
    out = st['outcome']
    n = st['n']
    if st['step'] == 13:
        out = (f"$F^*$ {alt_pct['iid_perfect']:.2f}/{alt_pct['persistent_perfect']:.2f}/{alt_pct['iid_imperfect']:.2f}\\% versus "
               f"{org_pct['iid_perfect']:.2f}/{org_pct['persistent_perfect']:.2f}/{org_pct['iid_imperfect']:.2f}\\% on the same seeds; all 22 rows within 0.3 pp; all three prespecified reproduction checks held")
        n = '256 backgrounds, full 22-row surface (budget rule resolved to the full surface)'
    rows.append((st['step'], st['study'], st['question'], st['p_selection'], n, out))
FLOW = rows
assert not any('PENDING' in str(x) for row in FLOW for x in row)

# slope input conditions (original workflow, used for the elasticity slopes)
CONDLAB = [('sens:clustered_low', 'Clustered, 0.75$\\times$', 'Sensitivity bank (A/B 128)'),
           ('sens:clustered_center', 'Clustered, center', 'Sensitivity bank (A/B 128)'),
           ('sens:clustered_high', 'Clustered, 1.25$\\times$', 'Sensitivity bank (A/B 128)'),
           ('sens:imperfect_low', 'Imperfect, 0.75$\\times$', 'Sensitivity bank (A/B 128)'),
           ('sens:imperfect_center', 'Imperfect, center', 'Sensitivity bank (A/B 128)'),
           ('sens:imperfect_high', 'Imperfect, 1.25$\\times$', 'Sensitivity bank (A/B 128)'),
           ('ind:center_matched256', 'Independent, center (matched)', 'Holdout seeds A[:128]+B[:128]'),
           ('ind:upper', 'Independent, 1.25$\\times$ (matched)', 'Holdout seeds A[:128]+B[:128]'),
           ('env:independent_low075_env1', 'Independent, environmental bank', 'Environmental bank (A/B 256)'),
           ('env:independent_low100_env1', 'Independent, environmental bank', 'Environmental bank (A/B 256)'),
           ('env:independent_low125_env1', 'Independent, environmental bank', 'Environmental bank (A/B 256)')]
rows = []
for key, lab, bank in CONDLAB:
    c = AD['conditions'][key]
    ri = c['control_rate_ci95']; fi = c['F_pct_ci95']
    rows.append(f"{lab} & {bank} & {c['n_seeds']} & {c['p']:.6f} & {c['control_rate']:.3f} [{ri['lo']:.3f}, {ri['hi']:.3f}] & {c['F_pct']:.2f} [{fi['lo']:.2f}, {fi['hi']:.2f}]\\\\")
    note('Table slope inputs', f"{key}: B {c['control_rate']:.3f}, F {c['F_pct']:.2f}", 'add', f"conditions.{key}")
write('tab_slope_inputs.tex', '\n'.join(rows) + '\n')

# ---------------------------------------------------------------------------
# 3b. Round-1 revision tables (reviews_round1/analysis/round1_numbers.json)
# ---------------------------------------------------------------------------
R1 = J['rev1']
LAB3 = [('independent', 'Independent misses'), ('clustered', 'Clustered misses'), ('imperfect', 'Imperfect clearing')]
# source status (imported vs acquired source carriers, patient-pickup acquisitions)
rows = []
for wf, wlab in (('original', 'Original'), ('alternative', 'Alternative')):
    for p_, lab in LAB3:
        im = R1['source_status'][wf][p_]['imported_source']
        ac = R1['source_status'][wf][p_]['acquired_source']
        rows.append(f"{wlab if p_=='independent' else ''} & {lab} & {im['control']:,} $\\to$ {im['policy']:,} & {r(im['estimate_pct'])} {ci(*im['ci95_pct'])} & {im['limit_pct']:.2f} & "
                    f"{ac['control']:,} $\\to$ {ac['policy']:,} & {r(ac['estimate_pct'])} {ci(*ac['ci95_pct'])} & {r(ac['two_gate_reference_pct'])}\\\\")
        assert im['limit_inside_ci'] and ac['reference_inside_ci']
        note('Table source status', f"{wf} {p_}: imported {im['control']}->{im['policy']} fall {im['estimate_pct']:.2f} {im['ci95_pct']}; acquired {ac['control']}->{ac['policy']} fall {ac['estimate_pct']:.2f} {ac['ci95_pct']}", 'rev1', f"source_status.{wf}.{p_}")
    if wf == 'original':
        rows.append('\\addlinespace')
write('tab_source_status.tex', '\n'.join(rows) + '\n')
# gap decomposition, original independent
GD = R1['gap_decomposition_original_independent']
rows = []
GLAB = [('patient', 'Patient pickup'), ('own_current_all', 'Own site, current stay'), ('other_zone', 'Other-zone site'), ('carry_over', 'Carry-over')]
for k, lab in GLAB:
    d = GD['labels'][k]
    rows.append(f"{lab} & {r(100*d['control_share'])} & {r(100*d['reduction'])} & {rs(d['gap_contribution_pp'],2)}\\\\")
    note('Table gap decomposition', f"{k}: share {100*d['control_share']:.2f}, reduction {100*d['reduction']:.2f}, contribution {d['gap_contribution_pp']:+.3f}", 'rev1', f"gap_decomposition_original_independent.labels.{k}")
F_count = org_pct['iid_perfect'] - GD['limit_minus_F_count_pp']
assert abs(F_count - 100 * (1 - 770 / 2375)) < 1e-9
rows.append(f"\\midrule All labels & 100.0 & {r(F_count)} & {rs(GD['sum_pp'],2)}\\\\")
assert abs(GD['sum_pp'] - GD['limit_minus_F_count_pp']) < 1e-9
write('tab_gap.tex', '\n'.join(rows) + '\n')
# face validity and resource metrics, alternative control arms
FV = R1['face_validity_alternative_control']
fv_rows = [
    ('Colonized share of patient-days (\\%)', 'colonized_share_of_patient_days_pct', 1),
    ('Mean fraction of bedside sites contaminated (\\%)', 'mean_contaminated_site_fraction_pct', 1),
    ('Mean census (occupancy \\%)', None, None),
    ('Hygiene opportunities per occupied patient-day', 'opportunities_per_occupied_patient_day', 1),
    ('Care-zone visits per occupied patient-day', 'bed_visits_per_patient_day', 1),
    ('Patient contacts per occupied patient-day', 'patient_contacts_per_patient_day', 1),
    ('Site contacts per occupied patient-day', 'site_contacts_per_patient_day', 1),
    ('Opportunities met with a contaminated hand (\\%)', 'opportunities_met_with_contaminated_hand_pct', 2),
    ('Hand pickups from bedside sites (\\%)', 'share_of_hand_pickups_from_sites_pct', 1),
    ('Extra completed actions per run (policy $-$ control)', 'extra_completed_actions_per_run', 0),
    ('Acquisitions prevented per run', 'acquisitions_prevented_per_run', 2),
    ('Extra completed actions per acquisition prevented', 'extra_completed_actions_per_acquisition_prevented', -2),
    ('Natural-episode $q_1$ (stationary)', 'natural_episode_q1', 4),
    ('Natural-episode $q_2$ (stationary)', 'natural_episode_q2', 4),
]
rows = []
for lab, key, d in fv_rows:
    cells = []
    for p_, _ in LAB3:
        v = FV[p_]
        if key is None:
            cells.append(f"{v['mean_census']:.1f} ({v['occupancy_pct']:.1f})")
        elif key.startswith('natural_episode'):
            st = v['stationary_q1' if key.endswith('q1') else 'stationary_q2']
            cells.append(f"{v[key]:.4f} ({st:.4f})")
        elif d == -2:
            cells.append(f"{round(v[key], -2):,.0f}")
        elif d == 0:
            cells.append(f"{v[key]:,.0f}")
        else:
            cells.append(f"{v[key]:.{d}f}")
    rows.append(f"{lab} & " + ' & '.join(cells) + '\\\\')
    note('Table face validity', f"{lab}: " + ' / '.join(cells), 'rev1', f"face_validity_alternative_control.*.{key}")
write('tab_face_validity.tex', '\n'.join(rows) + '\n')
# alternative observation halves
rows = []
for p_, lab in LAB3:
    h = R1['alternative_observation_halves'][p_]
    c, q_ = h['baseline'], h['policy']
    rows.append(f"{lab} & ${c['first']:.2f}\\to {c['last']:.2f}$ & {ci(*c['change_ci95'], d=2)} & ${q_['first']:.2f}\\to {q_['last']:.2f}$ & {ci(*q_['change_ci95'], d=2)}\\\\")
    note('Table alt halves', f"{p_}: control {c['first']:.3f}->{c['last']:.3f} {c['change_ci95']}; policy {q_['first']:.3f}->{q_['last']:.3f} {q_['change_ci95']}", 'rev1', f"alternative_observation_halves.{p_}")
write('tab_alt_halves.tex', '\n'.join(rows) + '\n')
# e* by completion range and isolated F2
rows = []
for row in R1['e_star_by_range']:
    es = f"{row['e_star']:.3f}" if row['e_star'] is not None else 'none'
    rows.append(f"{row['h0']:.2f}$\\to${row['h1']:.2f} & {row['h0_plus_h1']:.2f} & {es} & {row['F2_independent_pct']:.1f} & {row['F2_clustered_pct']:.1f} & {row['F2_imperfect_pct']:.1f} & {row['F2_combined_pct']:.1f}\\\\")
    note('Table e* range', f"{row['h0']}->{row['h1']}: e* {es}; F2 {row['F2_independent_pct']:.2f}/{row['F2_clustered_pct']:.2f}/{row['F2_imperfect_pct']:.2f}/{row['F2_combined_pct']:.2f}", 'rev1', 'e_star_by_range')
assert abs(R1['e_star_by_range'][5]['e_star'] - 0.92352808870981) < 1e-12
write('tab_estar_range.tex', '\n'.join(rows) + '\n')
# WHO-style transition compliance (closed form)
rows = []
for t, lab, (rho_, e_) in (('independent', 'Independent misses', (0, 1)), ('clustered', 'Clustered misses', (0.5, 1)),
                           ('imperfect', 'Imperfect clearing', (0, 0.8)), ('combined', 'Combined', (0.5, 0.8))):
    tc = R1['transition_compliance_pct'][t]
    for h_ in ('0.9', '0.95'):
        hh = float(h_)
        miss = (1 - hh) * (rho_ + (1 - rho_) * (1 - hh))
        assert abs(100 * (1 - miss) - tc[h_]) < 1e-9
        q2e = f"{miss:.4f}" if e_ == 1 else '---'
        rows.append(f"{lab + f' ({rho_:g}, {e_:g})' if h_ == '0.9' else ''} & {hh:.2f} & {tc[h_]:.2f} & {q2e}\\\\")
    note('Table WHO transition', f"{t}: {tc}", 'rev1', f"transition_compliance_pct.{t}")
write('tab_who.tex', '\n'.join(rows) + '\n')
# Isolated-bridge projection check (n = 8,000). Round 2 (23 Oct 2026): the round-1
# approximate ICU surrogate columns were removed; the ICU projection through the
# computed coefficient grid is tab_extension_projection.tex (make_extension_tables.py).
IP = R1['icu_projection_approximate']['results']
CW = R1['continuous_widths_n8000_pct']
rows = []
for t, lab in (('independent', 'Independent'), ('clustered', 'Clustered'), ('imperfect', 'Imperfect'), ('combined', 'Combined')):
    v = IP[t]['8000']
    rows.append(f"{lab} & {CW[t]['q1_only']:.1f} & {CW[t]['q1_and_q2']:.2f} & {v['median_width_pp_isolated_k2_same_grid_check']:.2f}\\\\")
    note('Table ICU projection (isolated check)', f"{t}: iso q1 {CW[t]['q1_only']:.2f}, iso joint {CW[t]['q1_and_q2']:.3f}, regrid {v['median_width_pp_isolated_k2_same_grid_check']:.3f}", 'rev1', f"continuous_widths_n8000_pct.{t}; icu_projection_approximate.results.{t}.8000.median_width_pp_isolated_k2_same_grid_check")
    assert abs(v['median_width_pp_isolated_k2_same_grid_check'] - CW[t]['q1_and_q2']) < 0.02, t
write('tab_icu_projection.tex', '\n'.join(rows) + '\n')

# ---------------------------------------------------------------------------
# 4. Exact checks of the isolated-bridge formulas
# ---------------------------------------------------------------------------
def q(k, h, rho, e):
    # stationary q_k = pi D (P D)^(k-1) 1 with exact Fractions
    h, rho, e = Fraction(h), Fraction(rho), Fraction(e)
    Pm = [[rho + (1 - rho) * (1 - h), (1 - rho) * h], [(1 - rho) * (1 - h), rho + (1 - rho) * h]]
    D = [Fraction(1), 1 - e]
    v = [(1 - h) * D[0], h * D[1]]
    for _ in range(k - 1):
        v = [(v[0] * Pm[0][0] + v[1] * Pm[1][0]) * D[0], (v[0] * Pm[0][1] + v[1] * Pm[1][1]) * D[1]]
    return v[0] + v[1]
F = lambda k, rho, e, h0='9/10', h1='19/20': 1 - q(k, h1, rho, e) / q(k, h0, rho, e)
EX = {
    'F2_ind': F(2, 0, 1), 'F2_clu': F(2, Fraction(1, 2), 1), 'F2_imp': F(2, 0, Fraction(4, 5)),
    'F2_rho1_e1': F(2, 1, 1), 'F2_clu_imp': F(2, Fraction(1, 2), Fraction(4, 5)), 'F2_rho1_e08': F(2, 1, Fraction(4, 5)),
    'F6_ind': F(6, 0, 1), 'F6_clu': F(6, Fraction(1, 2), 1), 'F6_imp': F(6, 0, Fraction(4, 5)),
}
assert EX['F2_ind'] == Fraction(3, 4) and EX['F2_clu'] == Fraction(23, 44) and EX['F2_rho1_e1'] == Fraction(1, 2)
assert EX['F2_imp'] == Fraction(13, 49) and EX['F2_clu_imp'] == Fraction(43, 134) and EX['F2_rho1_e08'] == Fraction(6, 17)
q1 = {n: q(1, '9/10', *pp) for n, pp in (('ind', (0, 1)), ('clu', (Fraction(1, 2), 1)), ('imp', (0, Fraction(4, 5))))}
q2v = {n: q(2, '9/10', *pp) for n, pp in (('ind', (0, 1)), ('clu', (Fraction(1, 2), 1)), ('imp', (0, Fraction(4, 5))))}
assert q1 == {'ind': Fraction(1, 10), 'clu': Fraction(1, 10), 'imp': Fraction(7, 25)}
assert q2v == {'ind': Fraction(1, 100), 'clu': Fraction(11, 200), 'imp': Fraction(49, 625)}
getcontext().prec = 40
h0, h1 = Decimal('0.9'), Decimal('0.95')
a, b = (h0 * (1 - h0)).sqrt(), (h1 * (1 - h1)).sqrt()
estar = (a - b) / (h1 * a - h0 * b)
assert str(estar).startswith('0.92352808870981')
f6 = [float(100 * EX[k]) for k in ('F6_ind', 'F6_clu', 'F6_imp')]
assert [round(x, 1) for x in f6] == [98.4, 60.4, 60.3], f6
note('S2/S3 exact', f"F2 = 3/4, 23/44, 13/49; rho=1: 1/2, 6/17; (0.5,0.8): 43/134; F6 = {f6}; e* = {estar}", 'computed (Fraction/Decimal) in this script', 'cross-check: mathematical_extension_2026-10-17/verification.json')
MV = J['mathv']
note('S3 checker', f"verification.json keys {list(MV.keys())[:8]}", 'mathv', 'top-level')

# ---------------------------------------------------------------------------
# 5. Write facts JSON and the number-source map
# ---------------------------------------------------------------------------
FACTS.update(dict(estar=str(estar)[:12], F6=f6, gate1=sites, undefined=UNDEF,
                  mc_halfwidths={f'{a}/{b}': round(v[0], 3) for (a, b), v in HW.items()},
                  combined={f'{a}/{b}': round(v[1], 3) for (a, b), v in HW.items()}))
(HERE / 'facts.json').write_text(json.dumps(FACTS, indent=1, default=str))
(HERE / 'study_flow.json').write_text(json.dumps(FLOW, indent=1))

lines = ['# Number sources for supplement_v2.tex', '',
         'Generated by `supplement_numbers/make_supplement_tables.py` (read-only on sources).',
         'Every generated table fragment (`supplement_numbers/tab_*.tex`) is listed below with its source file and field.',
         'Tables reused verbatim from the 17 October supplement/manuscript are listed in the second part.', '',
         '## Source files (sha256 at generation time)', '']
for k, v in SRC.items():
    lines.append(f'- `{k}`: `{v.relative_to(ROOT)}` (sha256 `{SHA[k][:16]}...`)')
lines += ['', '## Generated numbers', '', '| Where | Value(s) | Source | Field |', '|---|---|---|---|']
for w, t, s, f in NUM:
    lines.append(f'| {w} | {t.replace("|", "/")} | {s} | `{f}` |')
(HERE / 'NUMBERS_SOURCE_supplement_generated.md').write_text('\n'.join(lines) + '\n')
print('ok', len(NUM), 'numbers;', FACTS)
