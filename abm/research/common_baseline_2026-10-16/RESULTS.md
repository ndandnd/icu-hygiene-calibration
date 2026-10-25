# Reviewed common-baseline results — 17 October 2026

**Post hoc secondary analysis of existing simulations.** Under log-rate
interpolation to a baseline target of 3.0 acquisitions per 1,000 patient-days,
the three profiles retain different estimated responses to increasing hygiene
completion from 90% to 95%. All eight earlier matching checks remain failed.
The target is imposed by the interpolation, not established by fresh simulation.

## Primary results

| Profile | Baseline bracket | Implied p* | Prevention [95% MC interval] | Absolute rate prevented [95% MC interval] |
|---|---:|---:|---:|---:|
| Independent misses | 2.777–7.289 | 0.145853 | 66.7% [62.1, 71.2] | 2.00 [1.86, 2.14] |
| Clustered misses | 2.630–5.029 | 0.079545 | 55.0% [50.7, 59.4] | 1.65 [1.52, 1.78] |
| Imperfect clearing | 2.809–5.928 | 0.076328 | 30.0% [23.6, 36.6] | 0.90 [0.71, 1.10] |

Rates are per 1,000 patient-days. The absolute decrease is 3F*, not an
observed count difference at p*. No runs were added at these inferred p values.
The original 512-pair care-zone count benefits remain 3.13, 2.46 and 1.24
mean carriers prevented per 56-day ICU run at their original probabilities.

| Contrast | Difference in prevention (percentage points) | 98.333% MC interval |
|---|---:|---:|
| independent-clustered | 11.7 | [4.0, 19.2] |
| independent-imperfect | 36.7 | [26.8, 46.2] |
| clustered-imperfect | 25.0 | [15.3, 34.5] |

All three intervals exclude zero in this secondary family of three.
They are Bonferroni-adjusted percentile bootstrap intervals, not clinical
confidence bounds or a test of actual baseline equivalence.

## Estimator and uncertainty

For each profile, use 256 shared seed pairs at the center and 1.25-times-center
probabilities. Compute each control/policy rate as a ratio of pooled counts to
patient-days. Interpolate log control and log policy rates linearly in log p:

w = (log(3) − log(B_lower)) / (log(B_upper) − log(B_lower));
p* = p_lower^(1−w) p_upper^w;
F* = 1 − P_lower^(1−w) P_upper^w / 3.

Resample whole seed vectors 50,000 times within A/B banks of 128. Shared
conditions retain their pairing; the independent and sensitivity modules use
disjoint seed lists, and clustered/imperfect conditions share sensitivity draws.
The bootstrap re-estimates w and p* but conditions on the simulated endpoint
probabilities and the interpolation rule. It excludes interpolation error,
original grid/parameter-selection uncertainty, structural uncertainty and clinical
input uncertainty. See the preserved [protocol](PROTOCOL.md).

| Profile | Extrapolating draws / 50,000 | Undefined draws |
|---|---:|---:|
| Independent misses | 2,267 (4.534%) | 0 |
| Clustered misses | 21 (0.042%) | 0 |
| Imperfect clearing | 1,922 (3.844%) | 0 |

No draw was dropped as undefined. Extrapolating draws are retained.

## Review findings and corrections

- Two independently written implementations agree within the protocol tolerances;
  all 140 reconciled point comparisons pass and the largest interval endpoint
  difference across the reconciled analyses is 0.181 percentage points.
- A third implementation independently reread the raw data and ran 50,000 draws.
  It was written after inspecting A/B, so it is not described as blind.
  Primary point agreement is within 5.56e-16 and primary interval endpoint
  agreement within 0.058 percentage points (contrasts within 0.104).
- Independently checked 3,328 archived records, identical relevant simulator
  source maps, matching scenario settings apart from declared parameters, seed
  identities, and all 3,072 original holdout arms used for route counters.
- S3 wrongly called the independent lower512 and upper256 samples unmatched.
  The upper seeds are a subset of the lower. The original files remain intact;
  a [dated amendment](review_2026-10-17/statistics/AMENDMENT_S3_2026-10-17.md)
  and an overlap-aware four-stratum bootstrap correct this descriptive check.
  The primary analysis is unchanged.
- Report implementation A consistently. Its rounded independent lower95 bound
  is 62.1%; the two independent contrast upper bounds are 19.2 and 46.2 points.
  The handoff mixed rounded endpoints from the two bootstrap realizations.

## Descriptive robustness checks

- S1 arithmetic interpolation shifts F* by at most 0.03 percentage points.
- S2 three-point curvature shifts clustered and imperfect F* by −0.21 and
  −0.07 points. Independent lacks a third matched condition for this check.
- Corrected S3 estimates: independent 67.8% [64.8,70.8], clustered 55.5%
  [52.3,58.7], imperfect 27.5% [23.0,32.0]. Point estimates are unchanged;
  intervals now preserve the actual overlap. These are descriptive 95% intervals.
- S4 at target 2.7 extrapolates in 72.0% of independent and 84.1% of imperfect
  draws; it does not validate interpolation across the whole target band.
- A review-requested Linux-only subset of 192 independent pairs (96/bank)
  gives 66.64% [61.33,72.04], a −0.055-point change, with 18.846% extrapolating
  draws and none undefined. This changes seeds as well as platform mix and
  cannot establish platform equivalence.

These checks do not provide a bound on interpolation or model error.

## Runtime limitation

The primary independent lower bank contains 64 macOS/Python 3.12.2 and 192
Linux/Python 3.12.3 records; all upper and sensitivity records use Linux.
Scientific source hashes and recorded package versions agree, but exact
cross-platform replay is not verified. Platform effects could affect marginal
estimates, not merely pairing efficiency. The intervals condition on this
recorded execution mixture; the earlier handoff assertion is withdrawn.

## Verified passive pathway counts

| Profile | Arm | Site pickup, zero opportunities / all acquisitions | Share |
|---|---|---:|---:|
| Independent | baseline | 1,373/2,375 | 57.8% |
| Independent | policy | 547/770 | 71.0% |
| Clustered misses | baseline | 563/2,255 | 25.0% |
| Clustered misses | policy | 275/995 | 27.6% |
| Imperfect clearing | baseline | 554/2,352 | 23.6% |
| Imperfect clearing | policy | 440/1,718 | 25.6% |

These are labels of realized uncleared hand-contamination episodes. They
are not causal shares or a fixed hygiene-insensitive fraction. Counts can
fall while their share rises, as in the independent profile. Earlier-admission
labels do not separately identify contamination inherited through isolation swaps.
The source audit verifies synchronized daily cleaning only, no terminal or
relocation cleaning, fixed bed sites, 16 isolation / 14 shared beds, occupant
swaps, and the .6 multiplier on isolation patients and sites.

## Outputs and reproduction

- [Machine-readable final results](results_final.json), including exact values,
  all descriptive checks, claim limits and bound input hashes.
- [Statistical review](review_2026-10-17/statistics/REVIEW.md).
- [Provenance and model-rule review](review_2026-10-17/provenance_rules/REVIEW.md).
- [Third recomputation](review_2026-10-17/statistics/recompute.py),
  [Linux subset](review_2026-10-17/statistics/linux_subset_diagnostic.py), and
  [independent provenance/counter audit](review_2026-10-17/provenance_rules/audit.py).

To reproduce this assembly without rerunning statistics:

```sh
python3 finalize_results.py --output /path/to/fresh/results_directory
```

The numerical scripts require the archived local raw inputs; paths and runtime
details are preserved in their code/results. No new simulation, cluster task,
gate change, commit, publication or external message is implied by this report.
