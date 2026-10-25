# Leading-order coefficient extension: protocol

Written 23 October 2026, before any outcome of this study existed.

- **Scope.** Local only. No cluster jobs, and no edit to any frozen source,
  configuration, seedbook, record or accepted result.
- **Origin.** It answers round-1 referee items O-2, O-6 and O-10
  (`../manuscript_v2_2026-10-21/reviews_round1/OPEN_ITEMS.md`). It also adds
  one new question: does a policy defined on transition (WHO-style) compliance
  remove the effect of clustering?

## Design

- **Backgrounds.** The same 256 seeds as Table S5 and the workflow-coefficient
  study (`../full_icu_confirmation_2026-10-08/bridge/seedbook.json`). Settings
  are unchanged: 3% admission carriage, p = 0, 90 + 56 days, `zone_separate`,
  the frozen 30-bed model and histogram cap 2.
- **Workflow and engine.** Original workflow, using the verified adapter
  `../workflow_coefficient_2026-10-21/code/workflow_background.py` (sha256
  939aabd5…) unchanged, via `make_factory(alternative=False,
  surface_pairs=...)`.
- **Arms.** The 22 archived Table S5 pairs are evaluated first. Every new pair
  is an extra coefficient arm fed the identical event stream. Arms consume no
  random numbers.
- **Gate (every seed).** The ordered event hash, the 3 primary profiles and all
  22 archived pairs must equal the archived extended records **bitwise**. Any
  mismatch is a failure, and that seed is excluded and reported. One seed
  rerun in a fresh process must reproduce every new arm bitwise.

### New pairs (ids prefixed `x_`)

1. **Crossing location (O-2).** h .90→.95, ρ ∈ {0, .5} × e ∈ {.91, .92, .93,
   .94, .95, .96, .98}.
2. **Completion ranges (O-10).**
   - h .60→.70 and .70→.80, with (ρ, e) ∈ {(0, 1), (.5, 1), (0, .8),
     (.5, .8)};
   - h .80→.90 and .95→.99, with (ρ, e) = (.5, .8).
3. **Transition-defined policy (new).**
   - Here h is per boundary (zone exit or entry). Each profile's h0 and h1
     give WHO-style transition compliance 1 − (1−h)(ρ + (1−ρ)(1−h)) of 90%
     and then 95%.
   - ρ = 0: h .683772 → .776393. ρ = .5: h .829180 → .908392. Exact roots are
     computed in code and recorded.
   - Each is evaluated at e = 1 and at e = .8 (4 rows).
4. **Fine grid (O-6).** h .90→.95, with
   - ρ ∈ {0, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95};
   - e ∈ {.5, .6, .7, .75, .78, .8, .82, .85, .9, .95, .98, 1}.

   That is 132 rows. A row that duplicates an archived pair is still computed
   and must equal it bitwise, as an additional check.

## Predictions written before any outcome

The isolated two-opportunity values come from closed forms computed on 23
October.

1. **Crossing.** The ICU crossing e\* is the root of D(e) = F\*(ρ=.5) − F\*(ρ=0)
   at h .90→.95, found by linear interpolation between adjacent grid e values.
   It lies in [.91, .935]. (Isolated k=2 root: .9235.)
2. **Completion ranges.**
   - At e = .8, D < 0 for .60→.70 and .70→.80 (isolated −3.5 and −2.4 pp), and
     D > 0 for .80→.90 (isolated +1.2 pp).
   - At e = 1, D < 0 in every range.
3. **Transition-defined policy.** At e = 1, |F\*(ρ=.5) − F\*(ρ=0)| < 5 pp,
   against 21.5 pp under the per-boundary policy. (Isolated k=2: exactly 0.)
4. **Renewal comparator.** Over the 132-row fine grid, the largest
   |F_geo − F\*| is at most 3.5 pp. F_geo is the closed form defined in the
   manuscript (supplement S4.2); its largest deviation on the 22 archived rows
   was 2.78 pp.

All four are reported whatever the outcome.

## Secondary, specified now: exact ICU projection of the measurement study

- Reuse the saved assay confidence limits at 8,000 per assay. Use the same
  feasible-set construction as supplement S11 (e from the q1 interval; at each
  e, ρ from the q2 interval), with 81 e values × 61 ρ values.
- Replace the approximate surface with bilinear interpolation of the computed
  ICU F\* over the fine grid. The domain is limited to ρ ≤ .95 and e ≥ .5;
  sets beyond it are counted and reported, never extrapolated.
- Report median widths for the four truths, and compare them with the
  approximate values 6.2, 0.9, 6.7 and 11.7, labelled as a comparison.
- This is a projection of saved measurements. It adds no measurement or
  coverage claim.

## Budget and stopping rule

- Benchmark 2 seeds first.
- Cap: 60 local CPU-hours, at most 10 concurrent processes, run with
  `PYTHONDONTWRITEBYTECODE=1`.
- If the projection is above the cap, apply this fallback, fixed now: the fine
  grid on bank A only (128 seeds), everything else on all 256.
- No seed replacement and no change to settings after outcomes. Failures are
  reported.

## Inference

- F\* = 1 − mean(K1)/mean(K0) over seeds.
- Whole-seed bootstrap: 20,000 draws, `numpy.random.default_rng(20261023)`.
- Delta-method normal intervals as in the earlier studies.
- Crossings: report the bootstrap interval of the interpolated root.

## Limitations stated in advance

- These are fixed-horizon p → 0 limits. They give no finite-p error bound, and
  the finite-p reversal remains untested (O-3).
- Transition-defined compliance assumes the stationary Markov law and the
  model's exit/entry opportunity rule.
