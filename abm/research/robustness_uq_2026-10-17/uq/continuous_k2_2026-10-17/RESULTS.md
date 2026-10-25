# Continuous bounds preserve the information gain, with wider intervals

17 October 2026. **Post hoc sensitivity using the same saved measurements.** No
new random draws, cases or sample sizes were added. The original finite-grid
protocol, raw archive and numerical results remain unchanged.

For two-opportunity bridges at completion 90% to 95%, the projection now covers
every rho,e in [0,1]. Adding the two-opportunity assay still reduces uncertainty,
but the finite grid had narrowed some intervals materially. The continuous
results should be used for this target when discussing numerical precision.

[View the continuous versus grid comparison figure](continuous_uncertainty.png).

At 8,000 measurements per assay:

| Known truth | q1-only continuous width | q1+q2 grid width | q1+q2 continuous width | Continuous k=2 coverage | Empty joint sets |
|---|---:|---:|---:|---:|---:|
| Independent misses | 25.04 pp | 5.10 pp | 6.27 pp | 95.2% | 22/1,000 |
| Clustered misses | 25.03 pp | 0.91 pp | 1.13 pp | 97.0% | 15/1,000 |
| Imperfect clearing | 11.45 pp | 3.54 pp | 5.71 pp | 97.5% | 1/1,000 |
| Combined clustering and imperfect clearing | 11.45 pp | 4.91 pp | 6.69 pp | 99.8% | 0/1,000 |

Widths are medians among defined intervals, in percentage points of prevention.
Empty sets count as coverage failures. Coverage here concerns **one k=2 target**;
do not compare it with the original six-target simultaneous coverage as though
they were the same outcome. The original grid k=2 coverages were 95.2%, 96.4%,
96.2% and 97.5%, respectively. These empirical fractions have sampling variation;
the confidence argument comes from the construction, not observed performance.

The q1-only set uses a 97.5% marginal assay interval; the q1+q2 set has at least
95% joint parameter coverage under the observation model. The joint experiment
uses twice as many total measurements. This is not an equal-cost comparison or
an optimized allocation study. The continuous result does not add measurements
or statistical information; it removes a restriction on the candidate parameters.

## All fixed sample sizes

| Truth | Measurements per assay | q1-only continuous width | q1+q2 continuous width | Empty joint sets |
|---|---:|---:|---:|---:|
| Independent | 500 | 25.71 pp | 19.28 pp | 15 |
| Independent | 2,000 | 25.16 pp | 11.50 pp | 16 |
| Independent | 8,000 | 25.04 pp | 6.27 pp | 22 |
| Clustered | 500 | 25.71 pp | 5.68 pp | 14 |
| Clustered | 2,000 | 25.16 pp | 2.45 pp | 10 |
| Clustered | 8,000 | 25.03 pp | 1.13 pp | 15 |
| Imperfect | 500 | 19.02 pp | 18.35 pp | 1 |
| Imperfect | 2,000 | 14.06 pp | 10.51 pp | 1 |
| Imperfect | 8,000 | 11.45 pp | 5.71 pp | 1 |
| Combined | 500 | 19.01 pp | 18.97 pp | 0 |
| Combined | 2,000 | 14.07 pp | 12.98 pp | 0 |
| Combined | 8,000 | 11.45 pp | 6.69 pp | 0 |

Across all 12,000 experiments there are 95 empty joint sets, versus 110 on the
finite grid. The 15 newly nonempty sets were retained; no replacement measurements
were drawn. Every defined finite-grid interval is contained in its continuous
counterpart, within the separately checked floating-point comparison tolerance.

## Why these are continuous extrema

This is an analytic optimization, not a denser parameter sweep. At fixed e the
survival ratio is a ratio of affine functions of rho, so its extremes occur at
the feasible rho endpoints. The assay constraints give a single interval of
feasible e values. On each clipped rho boundary the prevention function is
nondecreasing in e: each piece fixes rho at zero or one, or fixes baseline q2 at
one of its assay limits. The derivatives in the [dated protocol and proof](POSTHOC_PROTOCOL.md)
are nonnegative. Consequently the global minimum uses the smallest feasible e,
and the global maximum the largest; only the two rho endpoints need evaluation.

The proof holds for 0<h0<h1<1 and the declared stationary Markov/constant-clearing
family. It provides exact sharp projection formulas in real arithmetic over the
whole rectangle. With exact input confidence limits, the usual set-projection
argument supplies the corresponding model-conditional coverage guarantee for
continuous parameters. There is no remaining finite-grid restriction for k=2.

## Numerical certification and independent checks

- The solver uses 60-digit Decimal interval arithmetic with directed rounding.
  Every square-root enclosure is checked by exact rational squaring.
- It evaluated 24,000 confidence sets: 23,842 feasible and 158 empty, counting
  q1-only and joint sets separately. There were **zero unresolved decisions**;
  no precision escalation beyond 60 digits was required.
- Each sharp extremum has a saved decimal enclosure, with maximum width
  **1e-58**. The final floating-point endpoints round outward; the largest
  outward padding is approximately **1.67e-16** in prevention fraction.
- An independent 100-digit scalar implementation, using a different rational
  expression for prevention, verified all 47,684 feasible extrema and all 158
  empty-set decisions. It also checked coverage against exact decimal truth
  parameters and verified the floating outputs enclose the scalar extrema.
- Every available grid endpoint was checked for inclusion, and q1+q2 intervals
  were checked for nesting inside q1-only intervals.
- A [separate scientific and exact-rational review](REVIEW.md) independently
  checked the proof, all 47,684 Decimal and floating extrema, every feasible/empty
  decision and all reported summaries. It found zero grid-bound noncontainment
  and no numerical or scientific mismatch. The review binds the input hashes.

The numerical certificate is for optimization **given the supplied assay limits**.
Those limits are the original floating-point Clopper–Pearson beta quantiles,
widened outward by 3e-12. Their exact-binomial tail identities were checked in
the original benchmark, but their evaluation was not itself formally certified
by interval beta inversion. Do not describe the whole statistical pipeline as
formally machine-certified. The widening covers the original grid's 2e-12
comparison tolerance and its floating rounding; it supplies no new information.

The independent scalar verifier initially treated roundoff of about 3.57e-101
around an exact point-valued certificate [.75,.75] as a failure. A documented
[verification-only repair](VERIFICATION_REPAIR.md) allows 1e-90 scalar comparison
error. The directed-rounding solver, confidence bounds and all scientific
outputs were unchanged. This repair was recorded before the verification rerun.

## Limits that remain

This result covers k=2 only. The saved k=6 intervals remain finite-grid results;
no continuous k=6 claim has been made. The non-Markov counterexample with matching
q1/q2 and a different k=6 response remains a valid limitation. Known h,
independently initialized assays, stationary first-order Markov completion,
homogeneous clearing and unchanged rho/e across policy arms are still required.
Nothing here estimates full-ICU finite-p, clinical, calibration or workflow
uncertainty. The original model-class limitation and doubled measurement budget
must accompany any presentation of the contraction result.

## Files and reproduction

- [Continuous results](run_001/results.json), [floating bounds](run_001/continuous_bounds.npz),
  [decimal certificates](run_001/decimal_certificates.json.gz), [execution receipt](run_001/receipt.json).
- [Solver](run_continuous.py), [independent verifier](verify_continuous.py),
  [verification receipt](run_001/independent_verification.json).

The solver binds the saved input and protocol hashes and refuses to overwrite
its `run_001` destination. To reproduce, copy this module without its output
directory into a fresh location preserving the relative layout alongside the
original `uq/run_001` measurement archive, then run:

```sh
python3 run_continuous.py
python3 verify_continuous.py
```

Both require Python and NumPy; verification additionally uses mpmath. No private
ICU inputs, cluster access, new random measurements or manuscript edits are needed.
