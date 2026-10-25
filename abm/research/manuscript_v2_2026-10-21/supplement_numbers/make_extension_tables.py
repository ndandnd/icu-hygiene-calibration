#!/usr/bin/env python3
"""Generate extension-only tables from the accepted, preserved coefficient summary.

Reads the original summary; performs no inference or simulation. Legacy JSON names
containing 'exact' are retained in provenance only: ICU projection is interpolated.
"""
from pathlib import Path
from decimal import Decimal, ROUND_HALF_UP
import hashlib
import json

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent.parent / 'coefficient_extension_2026-10-23/summary.json'
S = json.loads(SOURCE.read_text())
assert S['status'] == 'complete' and S['self_checks_pass']
assert S['n_records'] == S['n_passing'] == S['n_fine_grid_passing'] == 256
assert S['planned_surface'] == 'full'
assert not S['excluded_gate_failures'] and not S['record_problems']
assert len(S['pairs']) == 182 and S['bootstrap']['all']['draws'] == 20000
P = S['predictions']
files = {}
provenance = []

def num(x, dp=2, pct=False):
    if x is None:
        return 'undefined'
    value = Decimal(str(x)) * (100 if pct else 1)
    return format(value.quantize(Decimal(1).scaleb(-dp), rounding=ROUND_HALF_UP), f'.{dp}f')

def ci(block, dp=2, pct=True):
    interval = block.get('interval')
    if interval is None:
        return 'undefined'
    return '[' + ', '.join(num(v, dp, pct) for v in interval) + ']'

def write(name, rows, source_fields):
    target = HERE / name
    text = '\n'.join(rows) + '\n'
    target.write_text(text)
    files[name] = dict(sha256=hashlib.sha256(text.encode()).hexdigest(), rows=len(rows))
    provenance.append(dict(output=name, source_fields=source_fields, rows=rows))

p1 = P['prediction_1_crossing']
p2 = P['prediction_2_completion_ranges']
p3 = P['prediction_3_transition_policy']
p4 = P['prediction_4_renewal_comparator']
w = p4['rows'][p4['worst_row']]
# Round 2 (23 Oct 2026): criteria in the protocol's wording; the failing component of
# prediction 2 is named with its bootstrap sign fraction. Decisions are copied from summary.json.
fails2 = [k for k, v in p2['components'].items() if not v['holds']]
assert fails2 == ['e0.8_h0.7_to_0.8_e0.8'], fails2
f2 = p2['D_rows']['h0.7_to_0.8_e0.8']
assert f2['fraction_draws_positive'] == 1.0 and f2['bootstrap_95']['defined_draws'] == 20000
bs1 = p1['bootstrap_sign_change_counts']
assert bs1['draws_one_sign_change'] == bs1['draws'] == 20000
rows = [
    r'1 & ICU crossing $e^*$ (root of $D(e)=F^*(\rho{=}.5)-F^*(\rho{=}0)$ at $h$ .90$\to$.95) lies in $[.91,.935]$ & '
    + num(p1['e_star'], 5) + ' ' + ci(p1['e_star_bootstrap_95'], 5, False) + '; one sign change in all ' + f"{bs1['draws']:,}" + ' draws & ' + ('Held' if p1['holds'] else 'Failed') + r'\\',
    r'2 & At $e=.8$: $D<0$ for .60$\to$.70 and .70$\to$.80, $D>0$ for .80$\to$.90; at $e=1$, $D<0$ in every range & '
    + r'.70$\to$.80, $e=.8$: $D=+' + num(f2['D'], 2, True) + r'$ pp, predicted $<0$ (all ' + f"{f2['bootstrap_95']['defined_draws']:,}" + r' draws $>0$); other ' + str(sum(v['holds'] for v in p2['components'].values())) + ' of ' + str(len(p2['components'])) + ' components as predicted & ' + ('Held' if p2['holds'] else 'Failed') + r'\\',
    r'3 & Transition-defined policy: $|D(e=1)|<5$ pp (against ' + num(p3['per_boundary_abs_D_e1_pp'], 1) + r' pp per boundary) & $+' + num(p3['D']['e1']['D'], 2, True) + '$ pp ' + ci(p3['D']['e1']['bootstrap_95'], 2) + ' & ' + ('Held' if p3['holds'] else 'Failed') + r'\\',
    r'4 & Renewal comparator: $\max|F_{\rm geo}-F^*|\le3.5$ pp over the ' + str(p4['n_rows']) + '-row grid & ' + num(p4['max_abs_error_pp']) + r' pp ($\rho=' + num(w['rho']) + '$, $e=' + num(w['e']) + '$) & ' + ('Held' if p4['holds'] else 'Failed') + r'\\',
]
write('tab_extension_predictions.tex', rows, ['predictions'])
rows = []
for e, b in p1['D_by_e'].items():
    rows.append(' & '.join([num(float(e)), num(b['F_b'], 2, True), num(b['F_a'], 2, True), num(b['D'], 3, True), ci(b['bootstrap_95'], 3)]) + r'\\')
write('tab_extension_crossing.tex', rows, ['predictions.prediction_1_crossing.D_by_e'])
rows = []
for key, b in sorted(p2['D_rows'].items(), key=lambda kv: (kv[1]['h0'], kv[1]['h1'], kv[1]['e'])):
    rows.append(' & '.join([num(b['h0']) + r'$\to$' + num(b['h1']), num(b['e']), num(b['F_b'], 2, True), num(b['F_a'], 2, True), num(b['D'], 3, True), ci(b['bootstrap_95'], 3), num(b['isolated_k2_D'], 3, True)]) + r'\\')
write('tab_extension_ranges.tex', rows, ['predictions.prediction_2_completion_ranges.D_rows'])
rows = []
for key, b in p3['D'].items():
    rows.append(' & '.join([num(float(key[1:])), num(b['F_b'], 2, True), num(b['F_a'], 2, True), num(b['D'], 3, True), ci(b['bootstrap_95'], 3), num(b['isolated_k2_D'], 3, True)]) + r'\\')
write('tab_extension_transition.tex', rows, ['predictions.prediction_3_transition_policy.D'])
rows = []
for key, b in sorted(p4['rows'].items(), key=lambda kv: (kv[1]['rho'], kv[1]['e'])):
    rows.append(' & '.join([num(b['rho']), num(b['e']), num(b['F_star'], 2, True), '[' + ', '.join(num(v, 2, True) for v in b['F_star_bootstrap_95']) + ']', num(b['F_geo'], 2, True), num(b['F_geo_minus_F_star'], 3, True)]) + r'\\')
assert len(rows) == 132
write('tab_extension_fine_grid.tex', rows, ['predictions.prediction_4_renewal_comparator.rows'])
rows = []
Q = S['secondary_exact_icu_projection']
assert Q['status'] == 'computed' and Q['code_check']['passed']
for truth in ('independent', 'clustered', 'imperfect', 'combined'):
    b = Q['by_truth'][truth]
    assert b['defined'] + b['empty'] + b['out_of_domain']['total'] == b['experiments'] == 1000
    rows.append(' & '.join([truth.capitalize(), num(b['median_width_pp']), num(b['approximate_surface_median_pp_on_same_in_domain_sets']), str(b['defined']), str(b['empty']), str(b['out_of_domain']['total'])]) + r'\\')
write('tab_extension_projection.tex', rows, ['secondary_exact_icu_projection.by_truth'])

assert len(files) == 6
# Bind the post hoc interpolation diagnostic and prose values to the audit.
AUDIT = SOURCE.parent / 'review_2026-10-23/inference/independent_results.json'
A = json.loads(AUDIT.read_text())
checks = A['interpolation_holdout_checks']
assert len(checks) == 10
max_interp_error = max(abs(v['interpolation_error_pp']) for v in checks)
assert num(max_interp_error, 3) == '0.656'
for truth in ('independent', 'clustered', 'imperfect', 'combined'):
    assert A['projection'][truth]['domain_classification_disagreements'] == 0
    assert abs(A['projection'][truth]['median_width_pp'] - Q['by_truth'][truth]['median_width_pp']) < 1e-10
prose_values = {
    'cpu_hours_records': dict(value=S['compute']['cpu_hours_records'], display=num(S['compute']['cpu_hours_records']), field='compute.cpu_hours_records'),
    'process_cpu_hours_records': dict(value=S['compute']['process_cpu_hours_records'], display=num(S['compute']['process_cpu_hours_records']), field='compute.process_cpu_hours_records'),
    'crossing': dict(value=p1['e_star'], display=num(p1['e_star'], 5), field='predictions.prediction_1_crossing.e_star'),
    'crossing_interval': dict(value=p1['e_star_bootstrap_95']['interval'], display=ci(p1['e_star_bootstrap_95'], 5, False), field='predictions.prediction_1_crossing.e_star_bootstrap_95'),
    'renewal_max_error_pp': dict(value=p4['max_abs_error_pp'], display=num(p4['max_abs_error_pp']), field='predictions.prediction_4_renewal_comparator.max_abs_error_pp'),
    'interpolation_holdout_max_pp': dict(value=max_interp_error, display=num(max_interp_error, 3), source='independent_audit', field='max(abs(interpolation_holdout_checks[*].interpolation_error_pp))'),
}
receipt = dict(source=str(SOURCE.relative_to(HERE.parents[3])), source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
               script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               rounding='decimal half away from zero from repr(value); no intermediate rounding',
               description='Extension-only descriptive tables, no new inference. ICU projection remains interpolated and discretized.',
               files=files, provenance=provenance, prose_values=prose_values, independent_audit=dict(path=str(AUDIT.relative_to(HERE.parents[3])), sha256=hashlib.sha256(AUDIT.read_bytes()).hexdigest()))
(HERE / 'extension_number_provenance.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(dict(source_sha256=receipt['source_sha256'], tables=files), indent=2))
