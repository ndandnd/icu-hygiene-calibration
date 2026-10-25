# Workflow robustness study: verified results

Study `robustness_uq_2026-10-17`, execution 01. Report written 21 October 2026 by
an independent scientific reviewer. It draws on the two verifications in this
folder (`frozen_replay/` and `blind/`) and on this review's own re-tallies from the
raw records. These are results from a synthetic agent-based model. They are not
clinical evidence and the alternative workflow is not a validated clinical
practice.

## Bottom line

1. **The disagreement between the hygiene profiles persists.** Under the
   alternative workflow, all three prespecified contrasts exclude zero at
   98.333%. The ordering is unchanged: independent misses > clustered misses >
   imperfect clearing.
2. **Prevented fractions** (completion 90% to 95%) were 76.5%, 53.7% and
   27.3%. These are conditional on the once-selected transmission probability p.
3. **Rare-transmission limits.** The three fractions lie within 1.4 percentage
   points (pp) of the limits computed for the *original* workflow (75.1%,
   53.6%, 28.0%), and each limit lies inside its 95% interval. This agreement is
   descriptive only. The limit was not recomputed for this workflow, the
   comparison was not prespecified, and the intervals cannot tell the limit
   apart from the simpler two-opportunity benchmark. As before, the limit does
   not predict absolute rates.
4. **Target band.** Only the independent profile's baseline passed the
   descriptive target-band check. The clustered and imperfect holdout baselines
   were 2.62 and 2.63 per 1,000 patient-days, below the band's lower edge of
   2.7.
5. **Route labels (descriptive).** In the original workflow, one label hardly
   responded to hand hygiene: contamination at a bed that predates the current
   admission ("carry-over"). Under the alternative it almost disappeared.
   Site pickups contaminated during the current stay kept the same share, and
   they fell under the policy about as much as hand-borne routes did.
6. **Attribution.** The design changes vacancy cleaning and the isolation
   rule together. No difference from the original workflow can be attributed
   to either rule.
7. **Verification.** The cluster's numbers were reproduced locally, bit for
   bit, by a cross-platform replay and confirmed by a blind reimplementation.
   This was not an exact-runtime replay. The frozen analyzer refused the local
   runtime, and that gate was not bypassed.

## 1. Verification status

The headline numbers in this report come from the cluster analyzer's output,
`outputs/analysis.json` (SHA256 `063c0d95f376…579f`). It ran on the pinned
runtime (Linux, Python 3.12.3, NumPy 1.26.4) as Slurm job 366402. The
selection, `selection.json` (`1dcb9eac…6ecf`), ran as job 363325. After
collection, the following checks were made.

| Check | Result |
|---|---|
| Transport | All 3,583 collected files match their SHA256 receipts. |
| Frozen inputs | 71 of 71 hash checks pass. Freeze identity is `33fddd6c…8da6`. The derived policy freeze (`13d1aa1c…566b`) binds `selection.json`, `configs/policy.json` and all 1,536 grid record hashes. |
| Records | 1,536 grid records and 1,536 policy pairs were checked twice: by the unmodified frozen `load_bank` and by an independent pass. Both found complete scenario-by-seed products and 0 failures. Every record shows one runtime: Linux, Python 3.12.3 (GCC 13.3.0), NumPy 1.26.4, Mesa 1.2.1. Counter conservation, route histograms and the site mass balance all hold. No task restarted. |
| Selection replay | The frozen code, run cross-platform, reproduces all 636 numeric fields bitwise. The blind reimplementation gives identical counts and brackets, and p values within 0 to 3 ulp; the difference comes from summation order only. |
| Analysis replay | The frozen code, run cross-platform, reproduces all 139 numeric fields bitwise. The file's SHA256 matches the cluster's once the runtime descriptor is restored. The blind point estimates agree within 2e-15. The blind multinomial variant reproduces the cluster's bootstrap weight-stream hash (`8aad664b…`), with intervals within 4e-15. |
| Routes and workflow counters | This review re-tallied them from the raw policy records. They are identical to the blind tallies. |
| Frozen analyzer CLI, local runtime | **Refused** in both stages, exit 1, nothing written: `ValueError: Pinned analysis Python3.12.3/NumPy1.26.4 required` (`analyze.py` line 344). The gate was not bypassed. |

What these checks do and do not establish:

- **Cross-platform, not exact-runtime.** Both local replays ran on macOS arm64
  with Python 3.12.2 and NumPy 1.26.4.
  - The frozen-code replay calls the unmodified analyzer functions in
    production order.
  - It evaluated the two runtime-pin checks but did not enforce them. Both are
    false on this host.
  - No second run in the pinned runtime was made after collection.
  - `REVIEW_HANDOFF.md` lists a "matching-runtime analysis replay" in the
    closeout. Whether that item must be closed before manuscript use is the
    user's decision.
- **Records were not re-simulated.** No full-horizon record was re-simulated.
  Their contents rest on three things: the cluster run, its nine-fixture
  replay probe, and the integrity and invariant checks above.
- **Not peer review.** These are automated and delegated checks, not external
  human peer review.

## 2. What was run

The alternative workflow changes two rules together:

- **Vacancy cleaning.** Each bedside site is cleared when its patient leaves,
  by discharge or by relocation.
- **Empty-bed-only isolation.** Detected carriers move only into an empty
  isolation bed. If none is free, the carrier stays in its shared bed, with
  that bed's multiplier. There are no occupant swaps.

Everything else is inherited from the original model:

- 30 beds and 3% imported carriage.
- 90 warmup days and 56 observation days.
- Sites start clean, with daily scheduled cleaning.
- Zone-separate hygiene opportunities.
- Isolation multiplier .6.

The design ran in two stages:

- **Grid.** 4 probabilities × 3 profiles × 128 seeds, giving 1,536 baseline
  records.
- **Holdout.** One p was selected per profile from the grid. Each profile then
  received 512 fresh paired control/policy seeds (completion .90 vs .95), for
  1,536 pairs.

The scientific ceiling of 3,072 records was used exactly. Recorded compute was
66 CPU-hours for the grid and 101 for the policy stage, against a ceiling of
768.

## 3. Selection of p (grid stage, 128 seeds per point)

Baseline acquisitions per 1,000 patient-days at 0.5, 1, 2 and 4 times each
profile's original p. In every profile the selection used the lowest positive,
increasing adjacent bracket, interpolated log-log. There were no exact hits, no
nonincreasing segments and no zero lower rates.

| Profile | m0.5 | m1 | m2 | m4 | Bracket, log weight | Selected p (× m1) | Conditional 95% for p | Bootstrap draws choosing lower bracket |
|---|---:|---:|---:|---:|---|---|---|---:|
| Independent | 0.205 (43) | 1.919 (403) | 17.03 (3,574) | 96.73 (20,303) | m1–m2, w=.205 | 0.165111 (1.152) | [0.15999, 0.17033] | 0 / 20,000 |
| Clustered misses | 0.433 (91) | 2.463 (517) | 25.15 (5,280) | 161.72 (33,952) | m1–m2, w=.085 | 0.080633 (1.061) | [0.07801, 0.08332] | 7 / 20,000 |
| Imperfect clearing | 0.538 (113) | 2.828 (594) | 30.07 (6,312) | 168.58 (35,387) | m1–m2, w=.025 | 0.076144 (1.017) | [0.07338, 0.07862] | 2,861 / 20,000 (14.3%) |

Counts of acquisitions are in parentheses. The selection bootstrap had no failed
or undefined draws. The p intervals are conditional on valid draws, which here
were all draws. The policy intervals below do not propagate this p uncertainty,
and nothing in the study measures interpolation-model error.

## 4. Holdout results (512 paired seeds per profile)

| Profile | Control rate [95%] | In [2.7, 3.3]? | Policy rate [95%] | Prevented fraction F [95%] | Rate reduction [95%] | Carriers prevented per 56-day run [95%] |
|---|---|---|---|---|---|---|
| Independent | 2.904 [2.756, 3.052] | yes | 0.682 [0.622, 0.746] | **76.5% [74.1, 78.7]** | 2.22 [2.07, 2.37] | 3.64 [3.39, 3.89] |
| Clustered misses | 2.621 [2.480, 2.765] | **no** | 1.213 [1.122, 1.306] | **53.7% [49.7, 57.5]** | 1.41 [1.25, 1.56] | 2.31 [2.06, 2.56] |
| Imperfect clearing | 2.626 [2.484, 2.769] | **no** | 1.910 [1.796, 2.028] | **27.3% [22.2, 31.9]** | 0.72 [0.56, 0.87] | 1.17 [0.93, 1.42] |

Additional details:

- **Units.** Rates are acquisitions per 1,000 patient-days in the 56-day
  observation window.
- **Acquisition totals** (control → policy): 2,438 → 573, 2,201 → 1,019 and
  2,205 → 1,604.
- **Bootstrap.** 20,000 synchronous whole-seed draws, with 0 undefined draws for
  any statistic.
- **Target-band check.** It is descriptive, not a test of equivalence.

## 5. Primary contrasts (fixed Bonferroni family, 98.333% intervals)

| Contrast in F | Estimate | 98.333% interval |
|---|---:|---|
| Independent minus clustered | +22.8 pp | [17.1, 28.5] |
| Independent minus imperfect | +49.2 pp | [42.7, 56.0] |
| Clustered minus imperfect | +26.4 pp | [18.9, 34.3] |

All three intervals exclude zero. The contrasts compare three fixed
configurations whose baselines differ (2.90 vs 2.62 and 2.63). They are not
contrasts at a common baseline.

## 6. Comparison with the rare-transmission limits

### How the limits were computed

The manuscript's limits (Appendix "Full ICU model coefficient") come from 256
transmission-free (p = 0) histories of the **original** model. The relevant
code is `full_icu_confirmation_2026-10-08/bridge/background.py`, a
`BackgroundModel` subclass of the hash-pinned original `CPE_Model` with
`isolate_sick=True` and occupant swaps. Every other setting is the same as in
this study:

- 3% imports (config override), isolation multiplier .6.
- Lengths of stay of 5 and 14 days.
- Zone-separate opportunities.
- 90 + 56 days, completion .90 → .95.

These histories keep admissions, turnover, screening and isolation of imported
carriers.

K is a sum over every "bridge": a pickup from an imported carrier followed by
contact with a susceptible patient by the same worker. Each bridge is weighted
by the probability that the hand was not cleared in between. The limit is
F → 1 − K(.95)/K(.90).

### Could the alternative workflow change K?

- **Vacancy cleaning cannot.** With clean starting sites and p = 0, no site is
  ever contaminated, so clearing does nothing and consumes no random numbers.
  Site routes enter only beyond leading order.
- **Empty-bed-only isolation can.**
  - In the original, a detected carrier moves to isolation at once, swapping
    with a non-carrier if needed.
  - In the alternative, the carrier waits in its shared bed until an isolation
    bed is empty. While waiting, its contacts carry weight 1 rather than .6,
    and they fall at different points in staff routes.
  - No non-carrier is displaced.
  - The transmission-free history therefore differs.
- **How often it bites.** In this bank, a detected patient waited for a free
  isolation bed 1.7 patient-days per 56-day control record. That is about 1%
  of colonized patient-days and 0.1% of all patient-days.
  - The figure is the per-tick counter divided by 2,592 ticks per day.
  - It bounds how often the background differs, not how much K differs.
  - It suggests, but does not establish, that this workflow's limit is close
    to the original values.

### Comparison

| Profile | Original-workflow limit [MC 95%] | Stationary 2-opportunity comparator | Alternative F [95%] | F minus limit | Limit inside interval? | Original-workflow holdout F [95%] |
|---|---|---:|---|---:|---|---|
| Independent | 75.13 [75.10, 75.17] | 75.0 | 76.5 [74.1, 78.7] | +1.4 pp | yes | 67.6 [64.1, 70.8], limit outside |
| Clustered misses | 53.65 [53.64, 53.65] | 52.3 | 53.7 [49.7, 57.5] | +0.1 pp | yes | 55.9 [52.1, 59.5] |
| Imperfect clearing | 27.99 [27.99, 28.00] | 26.5 | 27.3 [22.2, 31.9] | −0.7 pp | yes | 27.0 [21.6, 32.1] |

**Absolute rates are still not predicted.** The leading-order rate is
1000 p² K / D⁰. Here D⁰ = 1,640.4 patient-days, and K is the original
workflow's value taken at the new p. The resulting rates are:

| Profile | Leading-order control rate | Measured control rate | Ratio | Leading-order policy rate | Measured policy rate | Ratio |
|---|---:|---:|---:|---:|---:|---:|
| Independent | 0.92 | 2.90 | 0.32 | 0.23 | 0.68 | 0.34 |
| Clustered misses | 1.62 | 2.62 | 0.62 | 0.75 | 1.21 | 0.62 |
| Imperfect clearing | 1.64 | 2.63 | 0.63 | 1.18 | 1.91 | 0.62 |

The original workflow's ratios were 0.25, 0.54 and 0.57.

### Judgment

"The analytic limit predicts the full-ICU response under the alternative
workflow" is **not defensible as worded**, for four reasons:

1. The limit was computed for a different workflow. Its isolation rule can
   change K.
2. The comparison was not prespecified. `PROTOCOL.md` does not mention it.
3. Intervals of ±2.3 to ±4.9 pp that contain the limit do not establish
   equality.
4. The stationary two-opportunity values (75.0/52.3/26.5%) also lie inside all
   three intervals. So the result does not specifically support the full
   coefficient calculation over the simpler benchmark.

Suggested wording:

> Under the alternative workflow, the full-ICU prevented fractions (76.5%,
> 53.7% and 27.3%) were within 1.4 percentage points of the rare-transmission
> limits computed from the original workflow's transmission-free histories
> (75.1%, 53.6% and 28.0%), and each limit lay inside the corresponding 95%
> Monte Carlo interval. This comparison is descriptive and was not
> prespecified. The coefficient was not recomputed for the alternative
> workflow: vacancy cleaning cannot change it, but the empty-bed-only isolation
> rule can. In the original workflow the independent profile's interval
> excluded its limit; under the alternative it contains it. As before, the
> leading-order calculation underestimated absolute rates, by factors of about
> 1.6 to 3.1.

## 7. Route labels (descriptive)

Each acquisition carries a passive label. The label marks the pickup that began
the infecting hand's current contamination episode, and how many hygiene
opportunities separated that pickup from the transfer. There are four labels:

- **Patient**: pickup from a patient.
- **Other-zone site**: pickup from a site in another care zone.
- **Own-zone site**: pickup from the patient's own bedside site, contaminated
  during the current stay.
- **Carry-over**: pickup from the own bedside site, where the contamination
  was deposited before this patient's admission.

All own-zone and carry-over pickups in the independent profile occurred with
zero opportunities before transfer. In every profile, patient pickups and
other-zone site pickups always had at least two.

### Shares of control acquisitions

The paper's metric is "site pickup with no intervening opportunity". It equals
own-zone|0 plus carry-over|0.

| Profile, arm | Paper's metric | of which own-zone, current stay | of which carry-over |
|---|---|---|---|
| Independent, original, control | 57.8% (1,373/2,375) | 38.9% | 18.9% |
| Independent, original, policy | 71.0% (547/770) | 29.6% | 41.4% |
| Independent, alternative, control | 41.7% (1,017/2,438) | 39.0% | 2.7% |
| Independent, alternative, policy | 42.1% (241/573) | 37.2% | 4.9% |
| Clustered, original, control / policy | 25.0% / 27.6% | 21.7% / 22.3% | 3.3% / 5.3% |
| Clustered, alternative, control / policy | 21.9% / 23.2% | 20.9% / 22.5% | 1.0% / 0.7% |
| Imperfect, original, control / policy | 23.6% / 25.6% | 20.9% / 22.2% | 2.6% / 3.4% |
| Imperfect, alternative, control / policy | 20.6% / 20.9% | 19.9% / 20.5% | 0.7% / 0.4% |

### Fall in acquisitions under the policy, by label

The table shows control → policy counts. Patient-days differ by less than
0.01% between arms, so count falls are close to rate falls.

| Profile, workflow | Patient | Other-zone site | Own-zone site | Carry-over | All | All except carry-over |
|---|---|---|---|---|---|---|
| Independent, original | 658→143 (−78%) | 344→80 (−77%) | 923→228 (−75%) | **450→319 (−29%)** | −67.6% | −76.6% |
| Independent, alternative | 824→195 (−76%) | 597→137 (−77%) | 952→213 (−78%) | 65→28 (−57%) | −76.5% | −77.0% |
| Clustered, original | 1,361→590 (−57%) | 324→130 (−60%) | 495→222 (−55%) | 75→53 (−29%) | −55.9% | −56.8% |
| Clustered, alternative | 1,342→640 (−52%) | 371→142 (−62%) | 466→230 (−51%) | 22→7 | −53.7% | −53.6% |
| Imperfect, original | 1,467→1,047 (−29%) | 331→231 (−30%) | 492→382 (−22%) | 62→58 (−6%) | −27.0% | −27.5% |
| Imperfect, alternative | 1,404→1,029 (−27%) | 347→238 (−31%) | 438→330 (−25%) | 16→7 | −27.3% | −27.0% |

In this table, "Own-zone site" counts own-zone pickups at any number of
opportunities, so the clustered counts differ slightly from the zero-opportunity
shares above (for example 495 against 489 in the original control arm). The
original-workflow counts come from
`full_icu_confirmation_2026-10-08/holdout_report/report.json`. The
alternative-workflow counts were re-tallied from the 1,536 policy records.

### What this shows, descriptively

- **The original metric mixed two labels that behaved very differently.**
  - Own-zone pickups contaminated during the current stay fell by 75% under
    the policy for the independent profile, about the size of the limit.
  - Carry-over fell by only 29%.
  - The whole rise in the paper's metric, from 57.8% to 71.0%, came from
    carry-over (18.9% → 41.4%). The current-stay share fell.
  - Excluding carry-over-labelled acquisitions, the original independent count
    fell 76.6%, close to the 75.1% limit. This is bookkeeping, not a
    counterfactual: each acquisition also seeds later ones.
- **A zero-opportunity final step does not by itself make a route insensitive
  to hygiene.** In both workflows, current-stay own-zone pickups responded like
  the hand routes.
  - These sites can only have been contaminated during the stay by a hand
    arriving from elsewhere. Such a hand crossed hygiene opportunities on the
    way.
  - The label that hardly responded was carry-over: contamination already at
    the bed before admission.
  - This is consistent with a chain containing no hygiene opportunity at all.
    For example, a carrier's own care could contaminate its bedside site, and
    the next occupant could pick it up. The labels do not record the donor, so
    this is not shown directly.
- **Under the alternative, carry-over almost disappeared.** Control counts fell
  from 450 to 65 (independent), 75 to 22 (clustered) and 62 to 16 (imperfect).
  - The current-stay own-zone share stayed the same (38.9% → 39.0%,
    independent).
  - The paper's combined metric fell (57.8% → 41.7%), but it no longer rose
    under the policy (41.7% → 42.1%).
  - By the model's rules, any remaining carry-over must have been deposited
    while the bed stood empty, because a move-in does not trigger cleaning.
- **Why carry-over mattered most for independent misses (hypothesis, not
  tested).**
  - That profile needs a larger p to reach the target (.143 vs about .076).
  - A route through a clean site needs four transmission events, not two.
  - Its weight relative to hand bridges therefore scales roughly with p².
    That is about 3.5 times larger for independent than for the other
    profiles at their original p.

### Limits on these labels

- **Not causal.** The labels are realized-path descriptions, not causal route
  fractions.
- **Swap-inherited contamination is mislabelled.** In the original workflow, a
  patient swapped into a carrier's former bed was not admitted to that bed. Any
  pickup of the carrier's contamination is therefore labelled own-zone, not
  carry-over. This is the limitation the manuscript already states.
- **Different conditions.** The original and alternative label counts come
  from different p, different seeds and different studies.
- **No declared comparison.** None of these between-workflow differences has a
  declared estimand or interval.

### Workflow counters (mean per 56-day record, control / policy)

| Counter | Independent | Clustered | Imperfect |
|---|---|---|---|
| Isolation relocations (equal to relocation clears) | 6.74 / 5.29 | 6.62 / 5.80 | 6.70 / 6.21 |
| Detected-patient-days waiting in a shared bed (no isolation bed empty) | 1.68 / 1.25 | 1.68 / 1.42 | 1.81 / 1.58 |
| Discharge clears that removed contamination | 5.90 / 4.31 | 2.31 / 1.75 | 2.22 / 1.88 |
| Relocation clears that removed contamination | 2.52 / 1.63 | 1.08 / 0.70 | 1.03 / 0.81 |
| Contaminated vacancy clears as % of depositions | 7.6 / 6.2 | 7.9 / 6.7 | 8.1 / 7.4 |

Notes on the counters:

- **Waiting time.** The counter `isolation_deferred_no_empty_bed` increments
  once per tick for each detected patient in a shared bed while no isolation
  bed is empty. It counts patient-ticks, not patients. The table divides it by
  2,592 ticks per day.
- **Other counters.** Discharge clears equal discharges, about 311 per record.
  Rejections of an occupied destination were 0 in every arm, and scheduled
  clears were 1,680 in every record.

## 8. What this design cannot attribute

- **Every original-versus-alternative difference.** This includes:
  - the independent F, which moved from 67.6% to 76.5%;
  - the independent-minus-clustered contrast, from 11.7 to 22.8 pp;
  - the selected p, which was higher by 15%, 6% and 2%;
  - the lower baseline at the original p (independent grid m1 1.92
    [1.70, 2.15] against the original holdout's 2.83 [2.66, 3.01]);
  - the collapse of carry-over labels.

  Each mixes two rule changes, a different selected p, fresh seeds and a
  separate study. None was declared as an estimand, and none has an interval.
  Do not say that a gap "grew", "shrank" or "was removed".
- **Vacancy cleaning versus isolation rule.** The protocol says this
  explicitly. The two rules are also entangled. With swaps, a relocation never
  leaves a bed vacant, because the displaced patient moves straight into the
  carrier's former bed. A "vacancy cleaning with swaps" cell would therefore
  need an extra rule, such as cleaning at swap, before any factorial design
  could separate them.
- **What can be said follows from the model's rules, not from estimated
  effects:**
  - Vacancy cleaning removes pre-admission contamination, except deposits made
    while a bed is empty.
  - Vacancy cleaning cannot change the p = 0 coefficient.
  - The isolation rule can change the p = 0 coefficient.

## 9. The missed target band for clustered and imperfect

**What happened.** The clustered and imperfect holdout baselines were 2.621
[2.480, 2.765] and 2.626 [2.484, 2.769]. Both point estimates are below 2.7,
so the descriptive check failed for both, and the flags in `analysis.json` are
correct. These are two more failed target checks and must be preserved. The
protocol forbids refitting or adding seeds to pass.

**Why it probably happened.** The 128-seed grid estimates at m1 are imprecise,
and acquisitions are overdispersed (variance/mean about 1.9). This review's
descriptive bootstrap gives 2.46 [2.18, 2.77] for clustered and 2.83
[2.51, 3.15] for imperfect.

- **Imperfect.** The selected p is only 1.7% above m1, yet the holdout (2.63)
  came in below the grid's m1 estimate (2.83). The grid point was probably
  high by chance.
- **Clustered.** Grid error and log-log curvature cannot be separated. The
  log-log slope is 2.5 below m1 and 3.35 above it.
- The protocol anticipated both kinds of error. Neither is propagated into the
  policy intervals.

**Does it matter for F?** Probably by about 1 pp or less, but this was not
measured for this workflow.

- The only evidence on flatness is the original workflow's p-sensitivity table
  (Supplementary Table S3, 256 pairs per point):
  - clustered: F of 60.5, 55.3 and 53.6% at baselines 1.30, 2.63 and 5.03;
  - imperfect: F of 28.9, 29.8 and 32.7% at baselines 1.25, 2.81 and 5.93.
- The steepest point-estimate slope there is about 7.4 pp per unit of log
  baseline. A 13% shortfall is 0.135 log units, or about 1 pp.
- Those sensitivity contrasts were unresolved. Even the largest change their
  intervals allowed, about 17 pp across about 0.7 log units, scales to about
  3 pp.
- Both figures are well below the F interval half-widths (4 to 5 pp) and the
  contrasts (22 to 49 pp).
- A lower p also moves F toward, not away from, the p → 0 limit.

**Where it does matter:**

1. **No common baseline.** The profiles were not compared at a common
   baseline, and nothing here establishes baseline matching.
2. **Absolute measures.** Absolute reductions and carriers prevented scale
   roughly with the baseline. Comparing them across profiles carries the 11%
   baseline gap between the independent profile and the other two.
3. **Clustered vs imperfect.** Their baselines happen to be within 0.005 of
   each other. That is a coincidence of these runs, not matching.

## 10. Context: the original workflow (different estimand)

| | Original workflow holdout (October 8 study) | Alternative workflow (this study) |
|---|---|---|
| p (independent / clustered / imperfect) | .14328 / .07602 / .07484 | .16511 / .08063 / .07614 |
| Control rate | 2.83 / 2.69 / 2.80 | 2.90 / 2.62 / 2.63 |
| F | 67.6 / 55.9 / 27.0% | 76.5 / 53.7 / 27.3% |
| Independent − clustered (98.333%) | 11.7 [5.5, 18.0] | 22.8 [17.1, 28.5] |
| Independent − imperfect | 40.6 [33.1, 48.4] | 49.2 [42.7, 56.0] |
| Clustered − imperfect | 28.9 [21.1, 37.0] | 26.4 [18.9, 34.3] |
| Carriers prevented per run | 3.13 / 2.46 / 1.24 | 3.64 / 2.31 / 1.17 |

These are separate studies with different p and different seeds. Side-by-side
display is context only and supports no claim of change.

## 11. Caveats to carry with any use

- **One combined workflow.** This is one illustrative alternative workflow, not
  a clinical practice and not a family of workflows. "Robust to workflow" would
  overgeneralize from a single alternative.
- **Conditional on one p.** The policy intervals condition on the once-selected
  p. They exclude selection uncertainty and interpolation-model error.
- **Fragile bracket for imperfect.** Imperfect clearing chose the lower bracket
  in 14.3% of selection-bootstrap draws. Its conditional p interval is still
  narrow, at about ±3.5%.
- **Precision of p.** Selected p is defined only to about 1 to 4 ulp, which
  depends on the summation order of patient-days. This has no practical
  consequence. All 1,536 policy records use the cluster's p exactly.
- **Monte Carlo only.** Intervals describe Monte Carlo sampling in a fixed
  synthetic model. They do not cover clinical-input or structural uncertainty.
- **Counter semantics.** The isolation-deferral counter counts patient-ticks,
  not patients.
- **Manuscript text would need updating.** The current sentence "We did not
  test terminal cleaning or alternatives to the isolation-swap rule" (main
  text, discussion of environmental persistence) would need revision if these
  results are integrated. Any revision must still state that the two rules were
  changed together.

## 12. Claims

**Supported**, with the qualifiers shown:

- "Under a combined alternative workflow, with vacated bedside sites cleared
  and isolation only into empty beds, all three prespecified between-profile
  contrasts in prevented fraction excluded zero at 98.333%. The ordering was
  independent > clustered > imperfect."
- "Prevented fractions were 76.5% [74.1, 78.7], 53.7% [49.7, 57.5] and 27.3%
  [22.2, 31.9]. They are conditional on the once-selected transmission
  probability."
- "Only the independent profile's holdout baseline lay inside the descriptive
  target band. The clustered and imperfect baselines (2.62, 2.63) fell below
  it, so the profiles were compared as fixed configurations, not at a matched
  baseline."
- "These fractions were within 1.4 percentage points of the rare-transmission
  limits computed for the original workflow, and each limit lay inside the
  corresponding 95% interval. This comparison was descriptive and not
  prespecified, and the coefficient was not recomputed for the alternative
  workflow."
- "In the original workflow, acquisitions labelled as pickups of contamination
  deposited before the current admission fell by only 29% under the policy,
  against 75 to 78% for other labels in the independent profile. Under the
  alternative workflow this label became rare. The labels are descriptive, not
  causal."
- "The design cannot attribute any difference from the original workflow to
  vacancy cleaning or to the isolation rule separately."

**Not supported:**

- "The analytic limit predicts the full-ICU response under the alternative
  workflow" or "the limit is confirmed or validated".
- "Vacancy cleaning reduced the surface route", "terminal cleaning explains the
  independent profile's shortfall", or any attribution to either rule.
- "The independent-versus-clustered gap grew (or doubled)", "the gap to the
  limit closed", or any significance claim about original-versus-alternative
  differences.
- "The hygiene-bypassing surface route was eliminated". Current-stay
  zero-opportunity site pickups kept a 39% share.
- "The conclusions are robust to workflow assumptions" or "robust across
  workflows".
- "Baselines were matched at 3 per 1,000 patient-days" or "the profiles had
  equal baselines".
- "The results were replayed in the pinned runtime" or "reproduced byte for
  byte in the pinned analysis environment". The local replays were
  cross-platform.
- "The policy intervals include calibration or selection uncertainty".

## 13. Provenance

| Item | Location |
|---|---|
| Cluster outputs | `research/robustness_uq_2026-10-17/collection_2026-10-21/execution01/outputs/analysis.json` and `final_status.json` (in the project root, not under `abm/`) |
| Selection | `…/execution01/source/research/robustness_uq_2026-10-17/selection.json` |
| Frozen-code cross-platform replay | `abm/research/robustness_uq_2026-10-17/verification_2026-10-21/frozen_replay/` (`VERIFICATION_SUMMARY.json`, `comparison.json`, `cli_attempt/`) |
| Blind reimplementation | `abm/research/robustness_uq_2026-10-17/verification_2026-10-21/blind/` (`PRE_COMPARISON_RECEIPT.json`, `blind_results.json`) |
| Original-workflow limits | `abm/research/full_icu_confirmation_2026-10-08/bridge/aggregate_report/summary.json`: K0 55.51 / 409.13 / 464.91, K1 13.80 / 189.65 / 334.77, D⁰ 1,640.36 |
| Original-workflow route counts | `abm/research/full_icu_confirmation_2026-10-08/holdout_report/report.json` |
| Analyzer and freeze hashes | `analyze.py` `2fda7ea7…b6ab5ac`; parent freeze identity `33fddd6c…8da6` |

The re-tallies of route labels and counters, the grid-rate bootstrap and the
leading-order rate calculations in sections 6, 7 and 9 were done by this review
for interpretation only. They are not analyzer output.
