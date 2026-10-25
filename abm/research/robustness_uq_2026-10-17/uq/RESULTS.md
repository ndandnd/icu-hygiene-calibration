# Executed measurement-to-prediction uncertainty benchmark

17 October 2026. **Synthetic controlled-assay experiment, not ICU data.** The
prospective [protocol](PROTOCOL.md) was saved before assay outcomes were drawn.
All 12,000 repeated experiments completed under the fixed design. No cluster,
simulator, old study or manuscript was changed by this module.

**Updated scope: the finite grid understates continuous interval widths for k=2.**
A separate [post hoc continuous-k=2 analysis](continuous_k2_2026-10-17/RESULTS.md)
now removes the finite-grid restriction for the 90% to 95% two-opportunity target.
At 8,000 measurements per assay, the joint median widths become 6.27, 1.13, 5.71
and 6.69 percentage points across the four truths, versus the grid's 5.10, 0.91,
3.54 and 4.91. The information gain persists: the corresponding q1-only continuous
widths are 25.04, 25.03, 11.45 and 11.45 points. Use this continuous sensitivity
when interpreting that target's precision. The original numerical tables,
findings, raw run and figure below remain preserved; k=6 remains grid-conditional.

[View the separate continuous-bounds figure](continuous_k2_2026-10-17/continuous_uncertainty.png)
and [independent review](continuous_k2_2026-10-17/REVIEW.md). This sensitivity uses
the same synthetic controlled assays with known h, adds no measurements and
retains the doubled assay-budget caveat.

## Main finding

Observing survival through two hygiene opportunities can substantially reduce
uncertainty in predicted hygiene benefit beyond observing one-opportunity
survival. The gain depends on the mechanism and measurement sample size. The
result has a genuine finite-sample confidence interpretation **within the declared
4,141-point stationary Markov family**, with known mean completion and controlled
assay assumptions. It is not a clinical posterior or a full-ICU predictive interval.

[View the uncertainty-contraction figure](uncertainty_contraction.png).

Each assay estimates a Bernoulli survival probability. For each of four known
truths and sample sizes 500, 2,000 and 8,000 per assay, we repeated the experiment
1,000 times. The two 97.5% exact Clopper–Pearson intervals define a joint parameter
set with at least 95% coverage. Enumerating every retained model gives sharp
prediction intervals over that finite set. Each interval includes uncertainty in
both clearing effectiveness and completion correlation; no prior is assigned.

At 8,000 measurements **per assay**, for completion increasing from 90% to 95%:

| Known truth | k=2 interval width: q1 only → q1+q2 | k=6 interval width: q1 only → q1+q2 | All-six-target coverage | Empty joint sets |
|---|---:|---:|---:|---:|
| Independent misses | 25.00 → 5.10 pp | 48.44 → 2.39 pp | 95.2% | 26/1,000 |
| Clustered misses | 25.00 → 0.91 pp | 48.44 → 3.68 pp | 96.1% | 15/1,000 |
| Imperfect clearing | 10.83 → 3.54 pp | 16.21 → 5.84 pp | 96.2% | 9/1,000 |
| Combined clustering and imperfect clearing | 10.83 → 4.91 pp | 16.21 → 7.57 pp | 97.3% | 0/1,000 |

Widths are medians among defined intervals, in percentage points of prevention.
Empty sets remain in the coverage denominator as failures; they are not replaced
or discarded. Across all 12,000 experiments there were 110 empty joint sets.
The guarantee covers all six predictions simultaneously because each is projected
from the same confidence set: policy completion .925, .950 and .975, each with
bridges of k=2 and k=6 opportunities. No policy-level measurement was fitted.

The empirical all-target coverage has Monte Carlo uncertainty. At n=8,000 its
exact 95% intervals are [93.69%,96.44%], [94.71%,97.21%], [94.82%,97.30%] and
[96.10%,98.21%], in table order. These intervals describe precision of the repeated
synthetic experiment. Parameter-set coverage in these cases was 95.2%, 94.7%,
95.5% and 96.0%; a realized fraction below 95% is compatible with its sampling
error and is not a reason to change the design. The mathematical coverage
statement comes from the confidence construction, not an observed pass/fail gate.

The second assay is not equally informative in every case. With n=500, the
combined mechanism's median k=2 interval width did not decrease, although other
targets did improve. The protocol's nested sets ensure that adding q2 never
widens any defined interval; they do not guarantee a large gain in every sample.
The joint experiment uses twice the total assay budget of q1 alone. This is a
comparison of adding a different observable, not an equal-cost allocation study
or proof that this sampling allocation is optimal. The q1-only set retains the
same 97.5% marginal interval; the joint confidence guarantee is at least 95%.

## A limitation demonstrated, not just stated

The prespecified non-Markov alternative has the same h, q1 and q2 as the combined
Markov mechanism, including at policy h values. It mixes episode-level completion
propensities, so its longer-run dependence differs. At h=.95 its k=6 prevention is
56.42%, versus 58.73% in the matched Markov model. Its k=2 prevention is identical.

Reusing the identical-distribution q1/q2 measurements, the Markov-based k=6
intervals include that alternative truth in only 86.8% of n=8,000 experiments
(Monte Carlo 95% interval [84.54%,88.84%]). There is no nominal coverage promise
for this alternative; it lies outside the assumed Markov family. At lower sample
sizes the wider intervals cover it more often. Better measurement can therefore
produce narrower, confidently model-dependent predictions unless the dependence
assumption is also assessed. The alternative and all sample sizes were specified
before drawing the measurements; no case was added after inspecting coverage.

## What this adds to the research

This is a concrete value-of-information result: the q2 observation can identify
dependence that average completion and one-opportunity survival cannot. It
complements the exact inverse formulas and sharp bridge-weight bounds in the
mathematical extension. It supplies a reproducible example of uncertainty
contraction with finite noisy measurements and simultaneous coverage, rather than
merely drawing a sensitivity surface.

It does not establish that these assays are available or practicable clinically,
or that their independent-initialization assumption holds for naturally selected
contamination episodes. Mean compliance was treated as known; uncertainty in h,
measurement error, policy-induced changes in rho/e and non-Markov mechanisms are
outside the primary confidence claim. The discrete grid matters: exact extrema
over it are not certified extrema or confidence bounds over continuous rho/e.

For the new ICU robustness study, keep this UQ demonstration separate from its
paired policy Monte Carlo intervals. Baseline selection bootstraps quantify
selection variability, conditional on the candidate grid and selector. Policy
outcomes at one selected p cannot propagate uncertainty across unobserved p
values without extra simulations or an explicitly justified response model.

## Verification and reproducibility

- [Machine-readable results](run_001/results.json), [raw counts, intervals and bounds](run_001/raw_benchmark.npz),
  and [execution receipt](run_001/receipt.json).
- [Runner](run_benchmark.py) binds the saved protocol and source hashes, actual
  Python/NumPy/SciPy versions, and output hashes. It refuses to overwrite a run.
- A [separate algorithmic verifier](verify_benchmark.py) reproduced every random
  count from the declared seed streams, explicitly enumerated completion words
  instead of using the runner's recursion, and recomputed every defined bound.
- [Verification receipt](run_001/independent_verification.json): all 12,000
  experiments reproduced; 285,924 projected endpoint comparisons passed.
  Maximum prediction discrepancy was 2.45e-15, below the declared 2e-12 tolerance.
  Exact-binomial tail identities for the assay intervals agreed within 4.91e-15.
  This verifier was written after results were seen, so it is not called blind.
- A [separate scientific reviewer](../REVIEW_UQ.md) checked the coverage argument,
  model-class limits and heterogeneous counterexample, then recomputed every
  summary width, contraction and alternative coverage from the raw archive:
  834 additional assertions passed, with no mismatch. This review was not blind.

The recorded runtime is Python 3.12.2, NumPy 1.26.4 and SciPy 1.13.1 on macOS.
This small synthetic benchmark has its own receipt and does not alter the pinned
runtime requirements of any ICU experiment.

Reproduce into a new directory:

```sh
python3 run_benchmark.py --output /absolute/path/to/fresh/run_directory
python3 verify_benchmark.py /absolute/path/to/fresh/run_directory
```

The verifier resolves the protocol/runner next to the run directory, as in the
preserved `uq/run_001` layout. To use an unrelated destination, preserve that
layout or adapt only the input-file location resolution before rerunning and
record the changed verifier hash. The scientific computation itself needs no
private raw ICU data.
