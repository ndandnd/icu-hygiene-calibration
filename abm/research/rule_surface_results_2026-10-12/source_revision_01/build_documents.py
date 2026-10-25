#!/usr/bin/env python3
"""Reconstruct a provisional source revision from exact preserved inputs; no statistics."""
import argparse
import difflib
import hashlib
import json
from pathlib import Path
import re


def sha(raw): return hashlib.sha256(raw).hexdigest()


def require(ok, message):
    if not ok: raise ValueError(message)


def identity(d):
    return sha(json.dumps({k:v for k,v in d.items() if k != 'identity_sha256'}, sort_keys=True, separators=(',', ':'), allow_nan=False).encode())


def number(value, digits=2):
    if value is None: return 'undefined'
    if digits == 2 and 0 < abs(value) < .1: digits = 3
    return f'{value:.{digits}f}'


def interval(values, digits=2, factor=1):
    if values is None: return 'unavailable'
    return '[' + ', '.join(number(v*factor, digits) for v in values) + ']'


def estimate(value, uncertainty, digits=2, factor=1):
    return number(None if value is None else value*factor, digits) + r'\,' + interval(uncertainty['interval'], digits, factor)


def table(caption, label, columns, header, rows):
    return (r'\begin{longtable}{' + columns + '}\n' + r'\caption{' + caption + r'}\label{' + label + r'}\\' + '\n' +
            r'\toprule' + '\n' + header + r'\\\midrule\endfirsthead' + '\n' + header + r'\\\midrule\endhead' + '\n' +
            '\n'.join(row + r'\\' for row in rows) + '\n' + r'\bottomrule\end{longtable}' + '\n\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path(__file__).resolve().parent)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); root = a.root.resolve(); out = a.output.resolve()
    require(not out.exists() and not out.is_relative_to(root/'input'), 'Fresh output outside preserved inputs required')
    manifest_raw = (root/'INPUTS.json').read_bytes(); manifest = json.loads(manifest_raw)
    require(sha(manifest_raw) == '9077ad20a077b83e3900ec874bc1fb668b9e42e3453ffdd89ff66aa8fe6ce4b0' and
            manifest['identity_sha256'] == identity(manifest) and manifest['canonical_install_or_publication_authorized_by_this_snapshot'] is False,
            'Preserved provisional input manifest')
    inputs = {}
    for name, entry in manifest['inputs'].items():
        source = root/name
        require(not any(q.is_symlink() for q in (source,*source.parents)), 'Symlink input refused')
        raw = source.read_bytes(); require(sha(raw) == entry['sha256'] and len(raw) == entry['size'], 'Changed preserved input: '+name)
        inputs[name] = raw
    require(sha(inputs['input/report.json']) == '9189e84179ca027e54ec03d2214137cfe55889ebaf11b9bd6276e0c56510f30f' and
            sha(inputs['input/selection.json']) == '13ecd228596a82994e1ae859d4c6bead6d0c275be54c23d02e41cb5ceddf7705', 'Exact accepted report/selection inputs')
    report = json.loads(inputs['input/report.json']); selection = json.loads(inputs['input/selection.json'])
    require(report['primary_family_size'] == 6 and report['primary_interval_level'] == 1-.05/6 and
            len(report['scenarios']) == 11 and len(report['primary_contrasts']) == 6, 'Fixed report scope')
    scenarios = {s['id']:s for s in report['scenarios']}
    manuscript = inputs['input/base/manuscript.tex'].decode(); supplement = inputs['input/base/supplement.tex'].decode()
    ledger = []

    def edit(text, old, new, label, target='main'):
        require(text.count(old) == 1, 'Unique source anchor: '+label)
        ledger.append(dict(target=target, label=label, before=old, after=new, before_sha256=sha(old.encode()), after_sha256=sha(new.encode())))
        return text.replace(old, new, 1)

    manuscript = edit(manuscript, r'''In follow-up runs, the inherited per-contact washing clock produced no
sampled baseline match, and varying transmission probability by 25\% in
either direction produced no resolved change in prevention for the clustered misses
and imperfect clearing profiles, although the intervals are too wide to show
insensitivity. The results show what baseline fitting leaves undetermined;
they are not a clinical validation.''', r'''An initial inherited-clock sweep found no sampled baseline match. A later
baseline-only selection and fresh policy comparison gave prevention of
72.1\% and 64.0\% for clustered misses and imperfect clearing, but both new
target-band checks failed. Disabling surface exchange increased independent
relative prevention by 9.7 percentage points while lowering the absolute
number prevented. Varying transmission probability by 25\% in the earlier
local sensitivity study produced no resolved change for the other two
profiles; a new independent upper-probability contrast was also unresolved.
The results concern fixed synthetic configurations, not forecasts at a
demonstrated common baseline or a clinical validation.''', 'abstract_new_followup')
    manuscript = edit(manuscript, '''clearing after a completed action. We apply the same baseline-only calibration strategy to each hygiene
mechanism: select its transmission probability using baseline acquisitions,
then hold that probability fixed while testing a common improvement in completion. A small calculation''', '''clearing after a completed action. We also compare rules for when washing
opportunities occur. The calibration principle is the same for each
mechanism and rule: choose the transmission probability using baseline
acquisitions, then hold it fixed while testing a common improvement in
completion. A small calculation''', 'introduction_profile_and_rule')
    methods = r'''
\subsection{Separate opportunity-rule and surface study}
We then designed a separate opportunity-rule and surface study after those
results and the research comments were known. Under the inherited
before/after-object clock, a baseline-only grid contained nine probabilities
for clustered misses and ten for imperfect clearing, with 128 fresh shared
seeds each (2,432 records). We selected the lowest exact target hit or the
lowest increasing adjacent bracket, interpolated within that bracket, and
froze the result before 512 fresh policy/control pairs per selected profile.
There was no fallback, extrapolation, grid expansion or reselection. An
independent $p=1$ comparison used 512 fresh pairs and was not fitted.
Surface exchange was disabled during warmup and observation at each original
probability using the same 512 holdout seeds per profile. The independent
upper probability used 256 pairs, matched to the first 128 seeds in each old
bank. The fixed study added 5,760 records; old configurations were not refitted.

This separate analysis uses 50,000 whole-seed draws, seed 20261011, and a
six-contrast family at 99.167\%: two clock-minus-care-zone comparisons,
three surface-disabled-minus-full comparisons and the independent
upper-minus-matched-center comparison. The old vectors use four fixed
128-seed strata (bank A/B core and remainder); fresh clock vectors use
independent A/B banks of 256. Other intervals are descriptive 95\%.
Each selected-clock target check requires its entire 99.167\% baseline
interval inside 2.7--3.3; it does not establish equivalence with the old
configurations. Policy intervals condition on the once-selected probability.
Undefined intervals, selection failures and old matching failures are retained.
New four-stratum intervals do not replace the archived original intervals.
'''
    manuscript = edit(manuscript, '\\section{What hand survival predicts}', methods+'\n\\section{What hand survival predicts}', 'separate_study_methods')
    manuscript = edit(manuscript, 'Under the inherited per-contact clock, the highest sampled baseline rates were',
                      'In the initial coarse inherited-clock sweep, the highest sampled baseline rates were', 'initial_clock_scope')
    manuscript = edit(manuscript, 'We made no policy forecast with this clock.',
                      'That initial sweep contained no policy forecast; the separately selected-clock comparisons follow below.', 'initial_no_policy_scope')
    names = {
        'clock_minus_zone_clustered_F': 'Inherited clock minus care-zone: clustered misses',
        'clock_minus_zone_imperfect_F': 'Inherited clock minus care-zone: imperfect clearing',
        'surfaces_off_minus_full_independent_F': 'Surface-disabled minus full: independent',
        'surfaces_off_minus_full_persistent_F': 'Surface-disabled minus full: clustered misses',
        'surfaces_off_minus_full_imperfect_F': 'Surface-disabled minus full: imperfect clearing',
        'upper_minus_center_independent_F': 'Upper minus matched center: independent',
    }
    primary_rows = [names[c['name']]+' & '+number(c['estimate']*100,1)+' & '+interval(c['uncertainty']['interval'],1,100)
                    for c in report['primary_contrasts']]
    primary = (r'''\begin{table}[htbp]\centering\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{>{\raggedright\arraybackslash}p{.52\linewidth}rr}\toprule
Comparison & Change in $F$ (pp) & 99.167\% interval\\\midrule
'''+ '\n'.join(row+r'\\' for row in primary_rows) + r'''
\bottomrule\end{tabular}
\caption{The separate study's six prespecified prevention contrasts. Positive values mean greater relative prevention in the first configuration. Intervals containing zero are unresolved. Clock comparisons change both the opportunity rule and its selected transmission probability; baseline matching is not established.}\label{tab:rule-surface-primary}
\end{table}
''')
    results = r'''
\subsection{Opportunity rules and surfaces at fixed configurations}
The fine grid selected $p=.283$ for clustered misses and $.345$ for
imperfect clearing (exact values are archived). Fresh 512-pair comparisons
gave prevention of 72.1\% [69.4,74.6] and 64.0\% [60.2,67.5], respectively,
versus 55.9\% and 27.0\% in the original care-zone configurations.
Both prespecified target checks failed: baseline rates were 3.09 and 2.85,
with 99.167\% intervals [2.85,3.33] and [2.62,3.09]. The first crosses
the upper 3.3 edge; the second crosses the lower 2.7 edge. Selecting a point
near the target therefore did not demonstrate baseline matching.
Table~\ref{tab:rule-surface-primary} gives all six primary contrasts.
'''+primary+r'''
At the original independent probability, disabling surface exchange lowered
the baseline from 2.83 to 0.70 per 1,000 patient-days while increasing
relative prevention from 67.6\% to 77.3\%. Mean carriers prevented fell
from 3.13 to 0.88 per 56-day run: a larger fraction of a smaller burden
can protect fewer patients. The independent prevention contrast was resolved;
the contrasts for clustered misses and imperfect clearing were not.
With surfaces disabled, baseline rates were 0.70, 1.41 and 1.66 across the
three profiles, compared with leading predictions of 0.69, 1.44 and 1.59.
This is a descriptive comparison of point estimates, not an error bound or
equivalence test.
These matched-seed worlds retain downstream feedback and are not identical
contact histories or causal route-attribution experiments.

The independent upper probability, $p=.179$, raised the baseline from 2.78
to 7.29 in the matched 256-seed comparison. Prevention point estimates rose
from 66.3\% to 71.0\%, but the 4.7-point contrast remained unresolved.
This upward point direction is opposite a simple dilution explanation;
its mechanism is not established. The original full 512-seed center remains
2.83 and 67.6\%, and is not replaced by the matched subset.

Under the inherited clock, the independent $p=1$ boundary comparison had a
baseline of 1.70 and prevention of 98.3\% [97.5,99.0]. It is unfitted and
does not establish a global maximum or impossibility of fitting. Its control
rate rose from 1.58 to 1.82 between observation halves, with a descriptive
95\% change interval [0.01,0.47]. This unadjusted diagnostic does not
establish equilibrium or sufficient warmup. Supplementary Tables~S15--S22
retain the grid, selection failures, condition results, target gates and
bank/half-window diagnostics.
'''
    manuscript = edit(manuscript, '\\subsection{Retrospective diagnostics}', results+'\n\\subsection{Retrospective diagnostics}', 'new_results_and_primary_table')
    manuscript = edit(manuscript, '''Direct simulation at the selected probabilities preserved the differences
between the three profiles, but the baseline-matching check failed. The
experiment therefore does not yet show different forecasts at an established
common baseline.''', '''The original confirmation preserved differences between the three profiles,
but its baseline-matching check failed. Both selected-clock target checks
also failed. The combined experiments therefore still do not show different
forecasts at an established common baseline.''', 'discussion_matching_scope')
    manuscript = edit(manuscript, '''The sensitivity contrasts were not resolved, the
independent profile was varied only below the target, and the prespecified
matching decision stands.''', '''Those earlier sensitivity contrasts were not resolved, and that earlier
independent series remained below the target. The new upper-probability
comparison extends the baseline above the target, but its prevention contrast
also remains unresolved. The original matching decision stands.''', 'discussion_upper_coverage')
    manuscript = edit(manuscript, '''reach: under the inherited clock, independent completion remains below the target
at every sampled probability, including $p=1$, while the other two profiles
bracket it locally.''', '''reach: the initial inherited-clock sweep kept independent completion below
the target at every sampled probability, including $p=1$, while the other
two profiles bracketed it locally. The finer baseline-only grid permitted
new selected-clock policy comparisons, but their fresh target checks failed.
The clock comparisons change both the opportunity convention and its
calibrated $p$, so their resolved differences do not isolate a pure clock
effect. A common completion percentage also refers to different modeled
opportunity denominators.''', 'discussion_clock_p_and_denominator')
    manuscript = edit(manuscript, '''the tested probabilities. The practical use of this work is to identify which''', r'''the tested probabilities. The separate follow-up reinforces the importance
of the opportunity rule and environmental feedback, while leaving common-baseline
equivalence and several prevention contrasts unresolved. Stationary hand
survival and the fixed-horizon $p\to0$ limit do not become exact predictions
for these finite-$p$ policies. The practical use of this work is to identify which''', 'conclusion_new_scope')
    manuscript = edit(manuscript, '''Frozen configurations, seed lists, source code, analysis scripts and all
8,576 simulation records from the confirmation and follow-up studies are
archived with checksums; the primary reports replay byte for byte in the
pinned analysis environment.''', '''Frozen configurations, seed lists, source code, analysis scripts and
14,336 archived records are preserved with checksums: the earlier 8,576
records plus 5,760 new records in this separate study. Its analysis reuses
1,536 original holdout records, which are not counted again; paired branches
are not separate records. The coefficient replay also reuses its original
background sample. Primary reports replay byte for byte in the pinned
analysis environment.''', 'availability_record_scope')

    labels = {
        'old_full_independent_target_candidate':'Care-zone: independent (old 512)',
        'old_full_persistent_target_candidate':'Care-zone: clustered misses (old 512)',
        'old_full_imperfect_target_candidate':'Care-zone: imperfect clearing (old 512)',
        'independent_surfaces_off':'Surface-disabled: independent',
        'persistent_surfaces_off':'Surface-disabled: clustered misses',
        'imperfect_surfaces_off':'Surface-disabled: imperfect clearing',
        'independent_upper':'Independent: upper (matched 256)',
        'independent_center_matched256':'Independent: center (matched 256)',
        'clustered_selected_clock':'Inherited clock: clustered misses',
        'imperfect_selected_clock':'Inherited clock: imperfect clearing',
        'independent_boundary_p1':'Inherited clock: independent boundary',
    }
    appendix = r'''
\clearpage
\section{Opportunity-rule and matched-surface follow-up}\label{app:rule-surface-followup}
This separate study was frozen after the earlier results and research
comments were known. Its 5,760 new records comprise a 19-condition baseline
grid (128 shared seeds per condition), two selected-clock comparisons
(512 fresh pairs each), an unfitted independent $p=1$ boundary (512 fresh
pairs), three surface-disabled conditions at the original fixed probabilities
(512 matched old seeds each), and an independent upper condition (256 matched
old seeds). All other declared model inputs, the 90-day warmup and 56-day
observation horizon remain unchanged. Removing surface exchanges includes
warmup and downstream feedback. No old configuration or matching decision
was replaced.

Baseline selection uses the lowest exact hit, otherwise the lowest increasing
adjacent bracket and linear interpolation; no fallback, extrapolation or
post-policy reselection is allowed. Policy intervals condition on the one
selected $p$. Selection uncertainty resamples the complete baseline grid
within A/B64 and remains a separate diagnostic; a failed draw makes its
unconditional interval unavailable. The effect analysis uses 50,000 whole-seed
draws, seed 20261011 and chunks of 500. Shared old vectors use four fixed
128-seed strata (A-core, A-remainder, B-core, B-remainder); matched upper
comparisons use the cores only. Fresh clock A/B256 vectors are independent
of old vectors. All six primary contrasts retain their 99.167\% family level;
other intervals are descriptive 95\%. A target check requires the entire
99.167\% baseline interval inside [2.7,3.3], and does not establish equivalence
to another configuration. Old matching failures stand.

Tables~S15--S22 copy point summaries and intervals from the exact frozen
selection and complete report. Displayed $p$ is rounded to three significant
figures; exact probabilities, denominators and resampling streams are archived.
No new grid-rate intervals are inferred from the selection output. The old
full-512 rows below use newly recomputed four-stratum intervals; they do not
replace the older two-bank intervals elsewhere in this paper.
{\footnotesize\setlength{\tabcolsep}{3pt}
'''
    grid_rows = []
    for profile,label in [('clustered','Clustered misses'),('imperfect','Imperfect clearing')]:
        for x in selection['profiles'][profile]['grid_points']:
            grid_rows.append(label+' & '+f"{x['p']:.3g}"+' & '+str(x['acquired_total'])+' & '+f"{x['patient_days_total']:,.0f}"+' & '+number(x['rate'],2))
    require(len(grid_rows)==19,'Exact nineteen grid point rows')
    appendix += table('Baseline-only fine grid under the inherited clock: 128 seeds per condition, A/B64. Rates per 1,000 patient-days; patient-days are rounded for display. These point summaries select a probability; they do not certify target-band matching.',
        'tab:new-rule-grid','lrrrr','Profile & $p$ & Acquisitions & Patient-days & Rate',grid_rows)
    selection_rows=[]
    for profile,label in [('clustered','Clustered misses'),('imperfect','Imperfect clearing')]:
        x=selection['profiles'][profile];u=x['selection_uncertainty']
        selection_rows.append(label+' & '+f"{x['selected_p']:.3g}"+' & '+interval(u['conditional_on_valid_selection95'],3)+' & '+str(u['failed_draws'])+'/50,000 & '+interval(u['unconditional95']['interval'],3))
    appendix += table('Once-selected probability and selection uncertainty. Conditional intervals retain valid bracket/hit selections; any failed draw makes the unconditional interval unavailable. No selection variability is propagated into unexecuted policy outcomes.',
        'tab:new-rule-selection','lrrrr','Profile & Selected $p$ & Conditional 95\\% & Failed draws & Unconditional 95\\%',selection_rows)
    condition_rows=[];absolute_rows=[];bank_rows=[];half_rows=[]
    for x in report['scenarios']:
        lab=labels[x['id']]
        condition_rows.append(lab+r'\newline $p='+f"{x['scenario']['p']:.3g}"+r'$, $n='+str(x['n'])+'$ & $'+estimate(x['control_rate'],x['control_rate95'])+'$ & $'+estimate(x['policy_rate'],x['policy_rate95'])+'$ & $'+estimate(x['F'],x['F95'],1,100)+'$')
        absolute_rows.append(lab+' & '+str(x['n'])+' & $'+estimate(x['absolute_rate_reduction'],x['absolute_rate_reduction95'])+'$ & $'+estimate(x['mean_carriers_prevented'],x['mean_carriers_prevented95'])+'$')
        ba,bb=x['banks']['A'],x['banks']['B']
        bank_rows.append(lab+' & '+str(ba['n'])+'/'+str(bb['n'])+' & '+number(ba['control_rate'])+'/'+number(bb['control_rate'])+' & '+number(ba['F']*100,1)+'/'+number(bb['F']*100,1))
        strings=[]
        for arm in ('baseline','policy'):
            y=x['observation_halves'][arm]
            strings.append('$'+number(y['first28_rate'])+r'\to '+number(y['last28_rate'])+r'$\newline $'+interval(y['last_minus_first95']['interval'])+'$')
        half_rows.append(lab+' & '+strings[0]+' & '+strings[1])
    appendix += table(r'All condition results: control/policy acquisition rates per 1,000 patient-days and prevention fraction in percent, each with descriptive 95\% intervals. All old full-512 points are unchanged; their four-stratum intervals are new diagnostics. The p=1 condition is unfitted.',
        'tab:new-rule-conditions',r'>{\raggedright\arraybackslash}p{.32\linewidth}rrr','Condition & Control rate & Policy rate & $F$ (\\%)',condition_rows)
    appendix += table(r'Absolute pooled rate reduction per 1,000 patient-days and mean paired count reduction per 56-day ICU run, with descriptive 95\% intervals. Branch exposures remain separate; paired branches are not independent records.',
        'tab:new-rule-absolute',r'>{\raggedright\arraybackslash}p{.38\linewidth}rrr','Condition & $n$ & Rate reduction & Mean count reduction',absolute_rows)
    appendix += table(r'The fixed six-member primary family. Estimates and 99.167\% intervals are percentage points of prevention. An interval containing zero is unresolved, not equivalence or insensitivity.',
        'tab:new-rule-contrasts',r'>{\raggedright\arraybackslash}p{.56\linewidth}rr','Contrast & Estimate (pp) & 99.167\\% interval',primary_rows)
    target_rows=[]
    for n,label in [('clustered_selected_clock','Clustered misses'),('imperfect_selected_clock','Imperfect clearing')]:
        x=scenarios[n];require(x['sampled_baseline_target_gate'] is False,'Preserve both failed target gates')
        target_rows.append(label+' & '+number(x['control_rate'])+' & '+interval(x['baseline_target_interval']['interval'])+' & Failed')
    appendix += table(r'Fresh selected-clock target checks: the entire family-level interval must lie inside [2.7,3.3]. Both checks fail. The clustered descriptive 95\% interval is inside the band, but it cannot replace this wider prespecified gate. No cross-configuration equivalence follows.',
        'tab:new-rule-targets','lrrl','Profile & Baseline rate & 99.167\\% interval & Target gate',target_rows)
    appendix += table('A/B bank point diagnostics: control rates per 1,000 patient-days and prevention fractions in percent. These point summaries are not formal between-bank equivalence tests. Matched upper/center banks contain 128 seeds each; full/fresh conditions contain 256 each.',
        'tab:new-rule-banks',r'>{\raggedright\arraybackslash}p{.42\linewidth}rrr','Condition & A/B seeds & Control A/B & $F$ A/B (\\%)',bank_rows)
    appendix += table(r'Observation-half diagnostics. Each cell gives the first-to-last 28-day rate and a descriptive 95\% interval for last minus first (per 1,000 patient-days). These multiple unadjusted checks do not establish equilibrium or sufficient warmup.',
        'tab:new-rule-halves',r'>{\raggedright\arraybackslash}p{.34\linewidth}>{\raggedright\arraybackslash}p{.29\linewidth}>{\raggedright\arraybackslash}p{.29\linewidth}',
        'Condition & Control first/last; change interval & Policy first/last; change interval',half_rows)
    appendix += r'''\endgroup
'''.replace(r'\endgroup','}')
    appendix += r'''
Every reported condition/primary/half-window bootstrap interval has zero
undefined draws. This differs from the single failed clustered selection draw
in Table~S16; that unconditional selection interval remains unavailable.
The boundary control half-window rise and the older imperfect-control decline
are descriptive finite-horizon diagnostics, not proof of stationarity.
The independent upper point response rises while its primary interval includes
zero; a universal simple dilution trend is not established. Stationary
$q_2/q_6$, the fixed-horizon $p\to0$ coefficient and these finite-$p$ ICU
responses remain separate quantities. None of these tables certifies clinical
inputs or changes the original failed-matching conclusion.
'''
    supplement = edit(supplement, '\\end{document}', appendix+'\n\\end{document}', 'new_supplement_s15_to_s22', 'supplement')
    original_main = inputs['input/base/manuscript.tex'].decode()
    original_supplement = inputs['input/base/supplement.tex'].decode()
    require(manuscript.split('\\appendix',1)[1] == original_main.split('\\appendix',1)[1], 'Entire original proof/appendix/bibliography suffix preserved')
    for kind in ('figure','table'):
        blocks = re.findall(r'\\begin\{'+kind+r'\}.*?\\end\{'+kind+r'\}', original_main, re.S)
        require(all(block in manuscript for block in blocks), 'Original '+kind+' blocks unchanged')
    require(supplement.startswith(original_supplement.split('\\end{document}',1)[0]), 'Entire original supplement source prefix preserved')
    for name, entry in manifest['inputs'].items():
        require(sha((root/name).read_bytes()) == entry['sha256'], 'Input changed during source construction')
    out.mkdir(parents=True); (out/'figures').mkdir()
    output_bytes={'manuscript.tex':manuscript.encode(),'supplement.tex':supplement.encode(),
                  'figures/weekly_control.png':inputs['input/assets/weekly_control.png']}
    for name,raw in output_bytes.items():
        with (out/name).open('xb') as f:f.write(raw)
    for name,new in [('manuscript.tex',manuscript),('supplement.tex',supplement)]:
        old=inputs['input/base/'+name].decode()
        diff=''.join(difflib.unified_diff(old.splitlines(keepends=True),new.splitlines(keepends=True),fromfile='preserved/base/'+name,tofile='candidate/'+name))
        with (out/(name+'.diff')).open('x') as f:f.write(diff)
    result=dict(schema_version=1,kind='reconstructible_provisional_rule_surface_document_revision',
        inputs_manifest_sha256=sha(manifest_raw),inputs_identity=manifest['identity_sha256'],
        builder_sha256=sha(Path(__file__).read_bytes()),inputs_sha256={n:sha(v) for n,v in inputs.items()},
        outputs_sha256={n:sha(v) for n,v in output_bytes.items()},edits=ledger,
        old_supplement_source_preserved_as_exact_prefix=supplement.startswith(inputs['input/base/supplement.tex'].decode().removesuffix('\\end{document}\n')),
        original_proofs_figures_bibliography_and_result_tables_changed=False,
        old_numerical_results_not_replaced=True,new_old_full_intervals_explicitly_separate=True,
        numerical_analysis_or_simulation_run=False,canonical_install_or_publication=False,
        status='provisional_source_pending_exact_replay_final_source_and_pdf_review')
    result['identity_sha256']=identity(result)
    with (out/'SOURCE_MANIFEST.json').open('x') as f:f.write(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'],source_manifest_sha256=sha((out/'SOURCE_MANIFEST.json').read_bytes()),identity_sha256=result['identity_sha256'],outputs=result['outputs_sha256'])))


if __name__=='__main__': main()
