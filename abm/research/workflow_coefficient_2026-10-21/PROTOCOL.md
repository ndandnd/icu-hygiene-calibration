# Rare-transmission coefficient under the alternative workflow: protocol

Written 21 October 2026, before any outcome of this study existed. The study is
local only: no cluster jobs, and no change to any frozen source, configuration,
seedbook, record or accepted result.

## Status of what is already known

- The verified alternative-workflow full-ICU results are already seen (see
  `../robustness_uq_2026-10-17/verification_2026-10-21/RESULTS.md`): prevented
  fractions 76.5%, 53.7% and 27.3%.
- The original-workflow coefficient limits are also known: 75.1%, 53.6% and
  28.0%, with the extended surface in supplement Table S5.
- So any comparison between this study's alternative-workflow limits and those
  finite-p results is **post hoc**, and will be labelled that way.
- What is prospective here is the coefficient values themselves, plus the
  predictions written below.

## Question

The alternative workflow clears a bedside site when its patient leaves, and
relocates a detected carrier only into an empty isolation bed (no occupant
swap). Under that workflow, what are the leading-order coefficients K_h and the
limiting prevention F* = 1 − mean(K1)/mean(K0)? Compute them:

- for the three profiles: (ρ, e) = (0, 1), (.5, 1), (0, .8);
- over the extended surface of supplement Table S5 (the same 22 rows).

How do they differ from the original-workflow coefficients on the **same 256
background seeds**?

## Design

- **Backgrounds:** the 256 seeds in
  `../full_icu_confirmation_2026-10-08/bridge/seedbook.json` (banks A and B,
  128 each). Every scientific setting matches that study: 3% admission carriage,
  clean hands and sites, p = 0, 90 + 56 days, `zone_separate`, the frozen
  30-bed model, histogram cap 2, and the same profile set and h-pairs as the
  archived extended run behind Table S5.
- **The only change** is the alternative workflow rules, taken from the frozen
  robustness `workflow_variant.py` semantics:
  - vacated sites are cleared at discharge and relocation;
  - a detected patient moves only into an empty isolation bed;
  - if no isolation bed is empty, the patient stays in place.
- At p = 0, site clearing cannot change hand bridges, because clean sites
  inject no environmental mass. It is kept only for fidelity. The isolation
  rule can change contact histories and isolation weights, so it can change K.
- **Implementation:** a new adapter that composes the frozen bridge
  `background.py` hooks with the frozen workflow-variant relocation and clearing
  rules, without editing either file.
  - **Gate 1:** with the alternative rules switched off, the adapter must
    reproduce the archived original-workflow records for at least 6 seeds:
    identical `ordered_bridge_events_sha256`, and K within 1e-9 relative.
  - **Gate 2:** with the rules switched on, it must pass unit tests showing the
    rule semantics:
    - no swaps;
    - moves only into empty isolation beds;
    - deferral when none is empty;
    - clears at discharge and relocation.
  - **Gate 3:** an independent code review before the batch run.
- **Paired comparison:** the same seeds give
  ΔK = K_alt − K_orig and ΔF* per profile. Seed-level paired bootstrap (whole
  seeds, 20,000 draws, numpy seed 20261021) gives 95% intervals. Nominal
  normal/delta intervals are reported as in the original study.

## Predictions written before any outcome

1. **Small change.** The three primary limits change by less than 2 percentage
   points from 75.1%, 53.6% and 28.0%. The only mechanism is the isolation rule
   changing contact histories, which affects K0 and K1 nearly proportionally.
2. **Sign reversal.** For 90%→95% completion, F*(ρ=.5, e) − F*(ρ=0, e) is
   negative at e = 1 and positive at e = .8, as it is on the original
   backgrounds (Table S5).
3. The crossing e* (by linear interpolation in e between .9 and 1, at ρ = .5)
   lies within .05 of the isolated two-opportunity value .9235.

All three are reported, whatever the outcome.

## Budget and stopping rule

- Benchmark 3 seeds first.
- Cap: 80 local CPU-hours, at most 10 concurrent processes.
- If the full 22-row surface on 256 seeds is projected above the cap, then (a
  rule fixed now, before outcomes) run:
  - the 3 primary profiles on all 256 seeds;
  - the full surface on bank A only (128 seeds).
- No seed replacement, no extension, no change to settings after outcomes.
- Failures are reported as failures.

## Outputs (all under this folder)

- adapter code, tests and review;
- per-seed records;
- `summary.json`;
- `RESULTS.md`, containing:
  - primary limits, both workflows, paired differences;
  - the full surface for both workflows;
  - the prediction checks;
  - the post hoc comparison with the verified finite-p results
    (76.5/53.7/27.3), labelled post hoc;
  - leading-order absolute rates versus the finite-p control rates (as in the
    original study).
- Limitation: the coefficient is a fixed-horizon p → 0 limit. It gives no
  finite-p error bound.
