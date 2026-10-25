# Overnight author-note follow-up

10 October 2026. This is a new, bounded follow-up designed after the completed
8 October confirmation outcomes and review of manuscript_revised_3. It does
not modify, refit, or extend that frozen study. All outcomes will be retained.
The user authorized cluster execution and automatic requeue on interruption.

## Positive-transmission modules

The unchanged response backend uses 30 beds, 19 staff role slots, 3% admission
carriage, initially clean hands and sites, synthetic geometric stays with
means 5/14 days selected at admission, immediate carrier detection/isolation,
isolation multiplier .6, daily complete site cleaning, 90 warmup days and 56
observed days. A common p governs hand/patient/site exchanges. This is a
synthetic colonization model, without clinical infection or mortality.

1. **Original opportunity clock:** independent completion (rho=0,e=1),
   clustered misses (rho=.5,e=1), and imperfect clearing (rho=0,e=.8), each at
   p=.15,.30,.60,1.00 and completion .90. Use the original `separate` before/
   after-object contact convention. Run 128 fresh seeds, shared across the 12
   conditions, in two disjoint banks of 64. These are baseline-only runs:
   h0=h1=.90, with one trajectory per job. Total 1,536 jobs.
2. **Probability sensitivity:** clustered misses and imperfect clearing at
   .75,1.00,1.25 times their fixed probabilities .07602158239372622 and
   .07483995991409355. Use zone entry/exit opportunities and completion
   .90 -> .95. Run 256 fresh paired seeds, shared across all six conditions,
   in two disjoint banks of 128. Total 1,536 paired jobs. The new center runs
   provide same-bank comparisons; old centers are not substituted.

Fresh positive-p seeds are deterministically generated from the study/module
namespace and audited against explicit seeds in accessible public research
JSON files. They are disjoint between modules and from the audited old seeds.
Shared seeds support paired contrasts but do not promise identical contacts:
conditional RNG consumption and endogenous isolation can change trajectories.

Primary clock estimates are pooled acquisitions per 1,000 patient-days, with
whole-seed synchronous bootstrap intervals at 99.583333333% for the 12-rate
family. A sampled point establishes approximate numerical target matching
only when its interval is wholly within [2.7,3.3]. Below/above-target points
and adjacent local brackets are reported separately. No global infeasibility,
monotonicity, or p=1 upper-bound claim is allowed. No new p is selected.

Primary sensitivity estimates are F=1-(pooled policy rate)/(pooled control
rate), and four within-profile low-minus-center/high-minus-center contrasts.
Their family uses 98.75% intervals. Also report 95% descriptive rates, F,
paired count reductions, and rates in both observation halves. An interval
containing zero does not establish insensitivity. A practical +/-5 percentage
point band may be shown descriptively, not substituted for an equivalence test.

All new bootstrap analyses use 50,000 draws and fixed seed 20261010. Resample
whole seeds, independently between statistical banks, preserving all condition
and branch pairing. Undefined draws are counted; interval certification is
unavailable if any required bootstrap statistic is undefined. These are
conditional Monte Carlo intervals, excluding input/selection uncertainty.

## Transmission-free coefficient replay

Replay the original 256 p=0 background seeds, preserving the original .90
completion/rho=0/e=1 epidemiological event sequence and background laws.
Integrate 22 coefficient pairs: rho=[0,.25,.5,.75] by e=[1,.9,.8,.6] for
.90 -> .95, and the three original profiles for .80 -> .90 and .95 -> .99.
Each coefficient schedule has its own warmup/control completion and policy
transition. The background simulator itself retains its original settings.
Original K values and ordered event hashes must reconcile with old records.

Record unweighted bridge bins 0..12 and >=13, observed patient/site contacts
and washing opportunities by staff role, and exposure. Preserve warmup versus
observation distinctions and shift-reset truncation. Do not store huge event
traces. Compute 1-mean(K1)/mean(K0), not the mean of per-seed ratios. The 22
limits are descriptive 95% calculations, without a familywise direction claim.
This extends exact conditional integration on existing background seeds; it
is not fresh independent epidemiological validation or a finite-p certificate.

## Existing-data author notes

Saved-data analyses are retrospective diagnostics, clearly separated from new
simulations: weekly pooled control rates over all 146 days; existing coarse
bridge counts/opportunity totals; bootstrap uncertainty of three-point rate
slopes with the six environmental conditions resampled together; and
uncertainty of the independent calibration's declared bracketing/interpolation
rule, including no-bracket frequency. No equilibrium claim follows from a
weekly plot; no asymptotic scaling law follows from three finite probabilities.
All hand/site ablations include downstream consequences, not isolated pickup
mechanisms or attributable route counts. Literature, parameter provenance and
human author declarations require documentary/author work, not compute.

## Execution, failure handling and preservation

Freeze source, configurations, seedbook and analysis plans before full-horizon
execution. Isolated output banks and atomic completed files reject changed
inputs, runtime or partition on resume. Interrupted incomplete trajectories
restart from their seed; completed ones are validated and skipped. Retryable
timeouts/interruption return 75; deterministic model/source errors return 2
and stop. Slurm requeues are bounded to 12 restarts. No manual duplicate jobs,
sample extensions, outcome-driven grid changes or source changes are permitted.

At most 48 simulation workers run concurrently: clock 24, sensitivity 16 and
coefficient replay 8. One CPU per worker; serial BLAS. A short implementation
probe precedes each production module. Resource-only changes, if a probe
requires them, are recorded before execution. Private cluster addresses and
execution receipts stay outside public commits. Old releases, advisor archives,
the canonical manuscript, paused goal and paused automation remain unchanged.
