# Hygiene response, acquisition scale, and environmental routes

8 October 2026. This is a new, exploratory follow-up prompted by a critique of
the completed study. All earlier outcomes were available when this protocol
was written. Earlier releases remain unchanged. It is not a new preregistration
of those results, and no post-hoc result will be presented as prospective.

## Questions and endpoints

1. How well does the probability that contaminated hands survive a two-opportunity
   care-zone gap predict the prevented fraction in the full ICU?
2. Does the prevented fraction stay stable as the common contact probability
   changes, including acquisition rates on the scale reported in Korean ICUs?
3. Does the answer depend on initially contaminated bed environments, persistent
   environmental transmission, or the convention used to count hygiene opportunities?

The principal outcome is the ratio of pooled acquisition rates, expressed as
`F = 1 - (policy acquisitions / policy patient-days) /
(control acquisitions / control patient-days)`. Counts, count-based prevented
fractions, susceptible-patient-day rates, census and importation are also retained.
Ratios of averages are used; per-run fractions are undefined when a control run
has zero acquisitions and must not be averaged or silently discarded.
All results are conditional synthetic colonization forecasts. Literature sets
a scale to investigate; it does not validate the simulator or identify parameters.

The hand-survival comparator is explicitly an approximation: stationary
`q_k = pi D (P D)^(k-1) 1`. Completion-state conditioning, environmental pickup
after entry hygiene, onward spread, turnover and isolation feedback can invalidate
its use as a whole-ICU forecast. The finite hand-only certificate has its original
scope and is not a certificate for this environmental model.

## Initial pilot, frozen before its outcomes

The initial pilot uses the three previously studied mechanisms:
independent completion/complete clearing `(rho,e)=(0,1)`, persistent completion
`(.5,1)`, and independent completion/imperfect clearing `(0,.8)`.
Each uses care-zone hygiene, baseline completion .90, common contact probabilities
`.025, .05, .10`, and eight newly allocated seeds. All 72 runs have 90 warmup
days followed by 28 baseline observation days. There is no policy intervention
in this pilot. Admission carriage is .03, initial bed environments are clean,
and the other previously declared synthetic ICU settings remain fixed.

Daily warmup and observation records include acquisition, importation, patient
exposure, colonized census, and tick-integrated environmental contamination.
End-of-day environmental snapshots alone are unsuitable because daily cleaning
is synchronized. The pilot will locate a contact-probability range around
3 acquisitions per 1,000 patient-days, motivated by (not fitted to individual
records from) Korean CRE surveillance. CRE and CPE endpoints are not interchangeable.

If a mechanism brackets this target, use piecewise linear interpolation of the
pooled baseline rate for a candidate center. If it does not, do not extrapolate
and call it fitted. At most one additional bracket extension is allowed, with
probabilities and seeds recorded before execution. Pilot policy outcomes are
never generated or used for selecting a center. A pilot center is only a design
choice; the next study must measure the actual rate at and around it.

## Main bank: freeze after the pilot and before new policy outcomes

The complete scenario table, seed list, sample size, site allocation and file
hashes will be recorded in `main_design.json` before its first simulation.
The intended design is a small mechanistic response study, not an exhaustive
factorial: several contact probabilities around each pilot center, the prior
high-transmission setting as a stress test, and selected initialization,
environmental-route and hygiene-clock controls. A lower-completion comparison
may halve missed opportunities (.70 to .85), using the same relative change as
.90 to .95. Its baseline and policy probabilities must be explicit.

The main bank uses 90 warmup days and longer observation than the pilot.
Control and policy branches start from identical saved warmup state, including
all RNG and hygiene-memory state. Daily warmup windows 31–60 and 61–90 and
early/late control windows will be compared descriptively. Lack of a significant
difference is not proof of equilibrium. Clean/dirty initialization is a separate
sensitivity check, not a reason to discard inconvenient trajectories.

The design records its exact finite computational budget. No sample extension
or parameter retuning is allowed based on a desired sign or significance.
Any later genuinely required follow-up is a new labeled experiment.

## Analysis

Show acquisition rates against contact probability, and prevented fractions
against the hand-survival prediction. Report every scenario, including undefined
ratios, wide intervals and poor agreement. Bootstrap whole seed blocks, preserving
all arms for a seed, for 95% pointwise Monte Carlo intervals. These are descriptive
exploratory intervals, not simultaneous confidence over a selected surface.
If reporting an interpolated forecast at the target, reselect the bracket and
recalculate interpolation in each bootstrap draw. Report the fraction of draws
without a bracket; no silent extrapolation or deletion.

Route logging is passive at its default settings: it consumes no random draws
and must reproduce original seeded trajectories. It records the last successful
hand pickup and intervening hygiene opportunities, with bed-environment history
where available. This is proximate source attribution, not a causal proportion
of infections prevented by deleting a route. Explicit environmental ablation is
a distinct model intervention and may alter downstream contacts and exposure.

## Verification and execution

Preserve frozen core files and earlier data. The new instrumented adapter is
versioned separately. Require seeded trajectory/RNG agreement, clone-versus-full
run agreement, counter conservation and restart identity checks before production.
Completed paired jobs are atomic checkpoints; incomplete jobs restart from their
seed after interruption. Source/configuration changes must refuse resume.
Use at most six local workers and twelve cluster workers, one CPU per simulation,
subject to scheduler availability. Never cancel, inspect results of, or modify
unrelated jobs. Infrastructure paths and credentials stay outside public files.

Retrospective reanalysis of the old data is labeled separately, including the
prevented-fraction calculation and the anomalous 32-seed training bank audit.
Old failed fitting remains documented in the repository and a concise manuscript
statement; it does not need a second full narrative in the main paper.
