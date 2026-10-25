# Direct confirmation of ICU hygiene forecasts

8 October 2026. This new study was designed after reading the completed response
study and its critique. It does not change that study or retrospectively register
its results. The following design is accepted before any new full-horizon run.

## Fixed scientific plan

`DESIGN_REVIEW.md` specifies the estimands, selection rule, seed construction,
uncertainty calculation, failure handling and scientific qualifications. Its
proposed scientific choices are adopted here, including the six baseline gates.
The executable configurations and seedbook govern the exact values.

1. Run 128 fresh baseline seeds at each of three independent-completion contact
   probabilities: .13496571354322734, .15 and .165. Select the lowest adjacent
   increasing bracket of 3 acquisitions per 1,000 patient-days (or lowest exact
   hit), and interpolate a candidate probability. If none exists, retain .15 as
   an explicitly unmatched reference. No further grid extension is allowed.
2. Lock that candidate before generating holdout outcomes. Run 512 fresh paired
   seeds for each of the three profiles. The other probabilities are fixed at
   .07602158239372622 (persistent completion) and .07483995991409355 (imperfect
   clearing). These directly test the earlier interpolations; they are not refit.
3. At each of the three original independent-profile probabilities
   (.08097942812593639, .10797257083458187, .13496571354322734), run 512 fresh
   paired seeds with and without environmental transmission. Retain all six
   conditions regardless of their results.

The budget is 384 baseline jobs and 4,608 control/policy pairs. Every condition
has 90 warmup days and 56 observation days, 3% admission carriage, clean initial
hands and bed environments, and the original frozen model and instrumented
runner. Control completion remains .90; policy completion changes to .95.
The calibration jobs have .90 in both arms and execute one trajectory.

Primary estimates pool counts and exposure before taking ratios. Whole-seed
bootstrap resampling preserves all arm and condition pairs. Use 200,000 draws.
The six approximate-baseline gates have Bonferroni 99.1666667% intervals: each
rate must lie within [2.7,3.3] and each pairwise ratio within [1/1.1,1.1]. These
are declared numerical tolerances, not clinical standards. All must pass, and
calibration must have found a bracket, before claiming approximate matching.
Undefined draws fail a gate. Bootstrap coverage is approximate.

The three between-profile prevented-fraction contrasts and the three environment
contrasts are separate primary families, each using 98.3333333% intervals for
joint directional claims. Also show descriptive 95% intervals, rates, absolute
prevented counts, baseline residuals, and attained interval widths. A failed
matching gate or imprecise result is an outcome, not a reason to retune or add
seeds. Holdout uncertainty is conditional on the selected probability;
calibration-selection uncertainty is reported separately.

Each family has two disjoint seed banks. Report both and bank B minus bank A;
resample independently between banks and synchronously across conditions within
each bank. Preserve first/last 28-day control-window contrasts. Neither a small
bank difference nor a small window difference establishes equivalence or
stationarity. No longer-warmup or six-opportunity experiment is added this round.

## Computation and immutable records

The scientific backend is unchanged. The new wrapper binds this protocol,
seedbook, configurations, source hashes and output scope to each bank. Completed
pairs are atomic and immutable; interrupted pairs restart from their seed.
Changes to source, configuration, partition or runtime refuse resume. Earlier
archives remain untouched. No unfinished or failed record enters an analysis.

Before execution, the technical site allocation is revised to one eighth of
each seed bank locally and seven eighths on the cluster, with at most six local
and 48 cluster workers. Both statistical banks occur at both sites. This change
addresses observed throughput and preserves every scientific seed and condition.
Cluster jobs use one CPU each and scheduler-managed memory. A runtime-only probe
and exact short replay check precede cluster production. Different site runtime
builds are recorded honestly; no old manifest is reused under a new build.

Separate RNG streams remain a separately tested adapter. They are not introduced
into this primary replication. Endogenous isolation and patient relocation mean
that no coupling can promise identical care contacts after disease paths differ.

## Interpretation and writing decisions

The proof is a fixed-horizon small-transmission expansion. Its remainder grows
with opportunities for transmission; small contact probability alone supplies
no accuracy guarantee over this long horizon. A large pL makes the elementary
gate-sparsity bound uninformative; it is not a necessary boundary of accuracy.
Environmental acquisition logs warn that omitted paths matter, but do not measure
the approximation's error or identify a causal fraction.

A separate p=0 coefficient module will be frozen after independent review. It
must match this study's 3% admissions and 90+56-day horizon. It will evaluate all
three profiles and both policies on the same backgrounds, including warmup
pickups, longer bridges, isolation weights and worker resets. No completed K
result is assumed here. Its separately frozen budget will not depend on policy
outcomes.

Environmental bypass and amplification through newly colonized patients are
testable explanations, not two universally signed corrections. Depletion,
isolation, turnover and competing routes can alter their effects. A measured
total environment ablation is not a decomposition of those mechanisms.

The ICU manuscript will retain its ICU proof and findings. The hand-only
certificate, algorithm and record illustrations move to a separate note. Both
the paper and reader guide will distinguish interpolation from fresh measurement,
simulation uncertainty from uncertain clinical inputs, and WHO-moment audits
from the model's entry/exit opportunities. Plain-language editing must preserve
these distinctions and every numerical qualification.
