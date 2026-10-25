# Workflow robustness and uncertainty study — 17 October 2026

The user authorized this new study after reviewing the prior studies and their
limitations. It is prospective for the new outputs, not for the known old results.
The old freezes, raw records, failed matching checks and current manuscript are
preserved. This protocol does not turn a synthetic workflow into clinical truth.

## Question and workflow

Does disagreement among the three hygiene profiles persist when every bedside
site is cleared as its patient leaves by discharge or relocation, and detected
patients move only to empty isolation beds? If none is empty, the patient stays
in the original bed, with the original location-dependent isolation multiplier.
There are no occupant swaps and no added contact-precaution effect. Site cleaning
clears contamination and its origin label without resetting the daily cleaning
clock. Initial admission to never-occupied beds adds no separate cleaning event.
All other source, admission, staffing, transfer, hygiene and timing rules stay as
in the preserved model. The changes apply during warmup and observation.

This is one combined alternative world. Its comparison cannot attribute a change
separately to vacancy cleaning versus isolation. Original-versus-alternative
results are descriptive structural sensitivity, not an isolated causal effect.

## Fixed workload and selection

`design.json`, `seedbook.json` and `configs/grid.json` are generated without reading
new outcomes. Use the three original care-zone probabilities multiplied by
0.5, 1, 2 and 4, with 128 shared fresh seeds: 1,536 baseline records. Each record
has 90 warmup days and 56 observation days at completion 0.90. Retain the
original independent (rho=0,e=1), clustered (rho=.5,e=1) and imperfect
(rho=0,e=.8) profiles. Sites start clean; admission carriage is .03.

Select each probability once using pooled baseline acquisitions per 1,000
patient-days. Target 3. Select the lowest-p exact target hit; otherwise use the
lowest-p adjacent increasing bracket with both rates positive. Interpolate log
rate linearly against log p. A zero lower rate is not a positive logarithmic
bracket. No extrapolation, pseudocount, alternate interpolator, fallback, grid
expansion or selection using policy outcomes. Preserve any failure and omit that
profile's selected-policy simulations.

After complete grid validation, hash-bind the selection and derived policy
configuration to the source/design freeze and every baseline input record.
Each selected profile receives 512 fresh paired control/policy seeds, shared
across profiles and disjoint from baseline seeds. Clone state after warmup;
completion is .90 versus .95, preserving memory. Maximum 1,536 pairs. Total
scientific ceiling is 3,072 records/pairs, with up to nine distinct short
one-day/one-day replay fixtures outside this scientific bank. Replay of the
same fixtures checks implementation and adds no inferential observations.

## Inference fixed before production

Use pooled rates R=1000 sum(acquisitions)/sum(patient-days), prevented fraction
F=1-R_policy/R_control, absolute rate change and mean acquisitions prevented.
Resample whole seed vectors across all available conditions/profiles, preserving
branch pairing. Grid and policy banks are independent. Use 20,000 bootstrap
draws, fixed seed 20261017. Resample the complete 128-seed grid bank or 512-seed
policy bank, with shared weights across conditions. Use independent streams
constructed by SeedSequence([20261017,0]) for grid and [20261017,1] for policy.
The three primary within-alternative contrasts are
independent-minus-clustered, independent-minus-imperfect and clustered-minus-
imperfect. Report 98.333333% intervals for this fixed Bonferroni family; missing
profiles do not redistribute alpha. Other intervals are descriptive 95%.

Report baseline 95% intervals and whether each lies wholly inside [2.7,3.3].
These are descriptive target-band checks, not profile-equivalence tests. Do not
refit or add seeds to pass. Reproduce selection under grid-seed bootstrap and
report selected-p uncertainty, invalid/no-bracket counts, and conditional-on-valid
intervals explicitly. If any bootstrap statistic is undefined, its unconditional
interval is unavailable; retain conditional summaries with failure counts.

Policy intervals are conditional on the once-selected p. They do not propagate
calibration uncertainty through an unobserved policy response surface and do not
include interpolation-model error. The broad factor-two grid can mispredict the
target; the independent holdout measures that discrepancy. No equality claim is
created by interpolation, selection or a nonsignificant baseline difference.

## Separate measurement uncertainty analysis

`uq/PROTOCOL.md` defines an independently seeded local synthetic measurement
benchmark. It studies uncertainty in the stationary hygiene law, with exact
simultaneous binomial confidence sets projected onto a declared mechanism family
and interventions withheld from measurement. It does not reuse full-ICU outcomes,
invent a clinical prior, or merge its envelopes into ICU Monte Carlo intervals.
Report feasible-family discretization, coverage, empty sets and misspecification
limits. No favorable case selection or adaptive increase in assay sample size.

## Implementation and execution gates

Use an isolated subclass of the hash-checked instrumented model; leave maintained
and earlier frozen source unchanged. Check vacancy cleaning, isolation occupancy,
blocked moves, site/lineage counters and reference replay before production.
Validate source/configuration/seed identities, full ordered Cartesian products,
runtime manifests, daily/cumulative counters and route histograms before analysis.
All production records use one pinned Linux runtime; local implementation fixtures
are separate. Source/config hashes are checked before and after each record.

The execution plan records runtime measurements, concurrency, wall-clock and retry
ceilings before launch. Respect OPTE priority and existing running/held work; do
not release or reprioritize another project. Current OPTE activity defers fresh
ABM submission. No completed-study monitor is revived. Unknown submission outcomes
fail closed and require receipt/accounting diagnosis. Completed records are atomic
and immutable; incomplete records restart from their seed within a bounded budget.
No scientific repair after outcomes without a separately preserved amendment.

Independent numerical/source review precedes manuscript claims. Preserve null,
reversed and undefined findings. No collaborator contact, Google Doc writes,
authorship invention, submission or payment is part of this authorization.
