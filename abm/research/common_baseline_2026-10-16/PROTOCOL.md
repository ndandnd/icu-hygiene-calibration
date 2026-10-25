# Prevention at a common baseline: post hoc secondary analysis

Written 16 October 2026, before any common-baseline estimate below was computed.
No new simulation. All inputs are archived, hash-listed records.

## Status and what was already seen

This analysis was declared after all baseline-matching gates had failed: six
original checks plus two fresh selected-clock checks. It is post hoc and
secondary. It does not replace or soften any failed gate, which remain reported
unchanged.

Seen before writing this protocol:

- The published point estimates of control rate and prevention fraction for
  every input condition, from the manuscript and the archived reports.
- A rough log-interpolated point calculation from those published numbers:
  about 68%, 55% and 30%.

No bootstrap, interval or contrast for the estimand below had been computed.

## Question

At the study's own prespecified target baseline of 3.0 acquisitions per 1,000
patient-days, what fraction of acquisitions does raising completion from 90% to
95% prevent in each care-zone profile? Do the three profiles differ?

## Estimand

For profile h and transmission probability p:

- B_h(p) and P_h(p) are the expected control and policy acquisition rates per
  1,000 patient-days over the 56-day observation window. Each is estimated as a
  ratio of sums (total acquisitions divided by total patient-days, times 1,000),
  matching the archived reports.
- F_h(p) = 1 − P_h(p)/B_h(p).

The common-baseline prevention is F*_h = F_h(p*_h), where B_h(p*_h) = 3.0.

## Inputs

Each profile uses two bracketing conditions run on identical seed lists (common
random numbers) within one study module. Every condition holds its frozen
profile settings fixed except p.

| Profile | Lower bracket (p) | Upper bracket (p) | Seeds | Module |
|---|---|---|---|---|
| Independent (ρ=0, e=1) | independent_center_matched256 (0.14327559) | independent_upper (0.17909449) | 256 shared | rule_surface_followup_2026-10-11 |
| Clustered (ρ=.5, e=1) | clustered_center (0.07602158) | clustered_high (0.09502698) | 256 shared | author_notes sensitivity |
| Imperfect (ρ=0, e=.8) | imperfect_center (0.07484) | imperfect_high (0.09355) | 256 shared | author_notes sensitivity |

The independent lower bracket consists of the original holdout records that the
rule-surface report identifies as its matched 256-seed center. The sensitivity
module uses one 256-seed list for all six of its conditions.

**Reproduction gate.** Before anything else, each implementation must
reproduce, from the raw records, the archived control rate, policy rate and F of
every input condition to a relative error of 1e-12. Archived values come from
`abm/research/rule_surface_results_2026-10-12/completed_analysis/report.json`
and
`research/author_notes_2026-10-10/analysis_2026-10-11/sensitivity_replay_20261011_1015_utc/report.json`.
Each implementation must also verify input file SHA256 hashes wherever an
archived manifest lists them. If this gate fails, stop and report.

## Primary estimator

Interpolate log rate linearly in log p between the two brackets, separately for
control and policy (a local power law, consistent with the observed slopes of
2–4). With

    w = (ln 3 − ln B_lo) / (ln B_hi − ln B_lo),

the interpolated control rate equals 3 exactly, and

    F* = 1 − (1 − F_lo)^(1−w) · (1 − F_hi)^w,

with implied p* = exp((1−w) ln p_lo + w ln p_hi).

## Uncertainty

- Use a whole-seed paired bootstrap with 50,000 draws and bootstrap seed
  20261016.
- In each draw, resample seeds with replacement within each seed bank (A and B,
  128 seeds each, as the archived reports define them). Apply the same
  resampled seed vector to every condition that shares that seed list.
- Recompute B, P, F and w from ratios of sums in each draw.
- The independent module and the sensitivity module have different seed lists,
  so they are resampled independently of each other. The clustered and
  imperfect conditions share the sensitivity seed list, so they are resampled
  synchronously.
- Report 95% percentile intervals for each F*_h.
- Report the three pairwise contrasts F*_i − F*_j with 98.333% percentile
  intervals, giving a Bonferroni family of 3 that matches the paper's
  between-profile family.
- Report the fraction of draws in which 3.0 falls outside the bracket (B_lo ≥ 3
  or B_hi ≤ 3), where the formula extrapolates.
- Report the number of undefined draws (zero policy count, or B_hi = B_lo).

## Reporting and decision rules (fixed in advance)

- Report every F*_h and every contrast, whatever the result.
- A contrast is "resolved" only if its 98.333% interval excludes zero.
- Unresolved pairs are stated as unresolved and are not described as differing.
- The original and fresh matching gates stay in the manuscript unchanged.
- This analysis is labelled post hoc and secondary wherever it appears.

## Descriptive sensitivity analyses (none can change the primary result)

- **S1:** arithmetic interpolation of F in ln B between the same brackets.
- **S2:** clustered and imperfect only. Fit a quadratic in ln p through the low,
  center and high conditions for ln B and ln P, then evaluate at B = 3. This is
  a curvature check.
- **S3:** use the original 512-seed holdout at the center p as the lower bracket
  for each profile. Its seeds are not matched to the upper bracket, so the two
  are resampled independently. This shows how much the result depends on the
  seed set. For example, the imperfect center is 27.0% in the 512-seed holdout
  but 29.8% in the 256-seed sensitivity center.
- **S4:** the primary estimator evaluated at B* = 2.7 and B* = 3.3, the edges of
  the target band.

## Descriptive tabulation (no inference)

- **D1:** from the original 512-seed holdout records, by profile and arm, the
  share of acquisitions in each `acquisitions_by_route_and_opportunities` cell,
  especially environment pickup with 0 opportunities.
- Before tabulating, verify the counters' semantics in the code.
- These are passive labels, not causal route shares.

## Replication

- Two implementations are written independently, each blind to the other's
  code.
- Point estimates must agree to 1e-10.
- Bootstrap intervals must agree within 0.5 percentage points (Monte Carlo
  tolerance). Any larger difference is diagnosed before results are used.

## Compute and limits

- Local computation only. No cluster jobs.
- Frozen sources, configurations, seedbooks and archived outputs are read-only.
- Outputs go only under `abm/research/common_baseline_2026-10-16/`.
- The interpolation assumes the rates are locally log-linear in p between the
  brackets. S2 checks that assumption where three points exist.
- Like the archived policy intervals, these intervals condition on the
  simulated p values. They describe Monte Carlo uncertainty only, not clinical
  or structural uncertainty.
