# Leading-order coefficient extension: results

Computed 23 October 2026. The batch finished at 09:32 EDT and `summary.json`
was written at about 18:05 EDT. Protocol: `PROTOCOL.md`, written before any
outcome. Launch policies: `LAUNCH_NOTES.md`. Code review: `CODE_REVIEW.md`
(approve with minor notes). Everything ran locally; nothing frozen was edited.

## Execution and checks

- **Completion.** 256 of 256 seeds ran the full plan of 182 pairs (364
  coefficient arms) per seed, under the original workflow.
  - No failures and no exclusions. Every seed passed the bitwise gate: event
    hash, 3 primary profiles, all 22 archived Table S5 pairs, and 12 duplicated
    rows, all identical to the archived extended records.
- **Compute.** 54.9 CPU-hours, about 6 h 25 min wall at 10 workers.
  - The benchmark had projected 48.9 CPU-hours, so actual cost was about 12%
    above it. The difference reflects machine load. It was still under the
    60-hour cap, which governed only the plan decision.
- **Pre-summary checks** (`LAUNCH_NOTES.md`): all passed.
  - All 256 records were written by the same code, and the code hashes equal
    those in the review.
  - Every record has `surface == 'all'`.
  - The 2 benchmark seeds recomputed in the batch, and the determinism rerun,
    are bitwise identical in every science field.
- **Self-checks in `summarize.py`:** pass.
- **Independent recomputation (parent session, from raw records).** The
  transition-policy differences, the crossing and the 70→80% e=.8 difference
  agree with `summary.json` to the printed digits.
- **Disclosure.** Before launch, the implementer saw outcomes on 2 seeds when
  testing `summarize.py`. The plan and protocol did not change afterwards.

## Prediction checks (written before outcomes)

| # | Prediction | Outcome | Verdict |
|---|---|---|---|
| 1 | ICU crossing e\* at h .90→.95 lies in [.91, .935] | e\* = .9336 [.9336, .9337]; unique sign change; one sign change in all 20,000 draws | **Holds** (near the upper edge) |
| 2 | At e = .8: D < 0 for 60→70% and 70→80%, D > 0 for 80→90%; at e = 1, D < 0 in every range | 60→70%: −0.94 pp. **70→80%: +0.25 pp.** 80→90%: +4.38 pp. D(e=1) < 0 in all five ranges | **Fails** (one component, by 0.25 pp) |
| 3 | Transition-defined policy: \|D(e=1)\| < 5 pp (against 21.5 pp per boundary) | +0.12 pp [0.12, 0.12] | **Holds** |
| 4 | Renewal comparator: max \|F_geo − F\*\| ≤ 3.5 pp over the 132-row fine grid | 4.32 pp (ρ=.95, e=.5) | **Fails** |

Notes:

- **Prediction 1.** The earlier linear interpolation between e=.9 and e=1
  (.919) understated the crossing, because D(e) is strongly nonlinear in e.
  - D(e) by e: .5: +3.9; .6: +5.6; .7: +7.7; .8: +9.0 (maximum); .85: +8.3;
    .9: +5.0; .91: +3.8; .92: +2.4; .93: +0.7; .94: −1.3; .95: −3.5;
    .96: −6.2; .98: −12.8; 1: −21.5 pp.
  - The ICU crossing (.934) is **above** the isolated two-opportunity root
    (.9235), even though 26% of the ρ=.5 bridge mass has 3 or more
    opportunities, and longer isolated bridges cross lower.
  - The renewal comparator's crossing (.936, round 1) is close to the ICU
    value.
- **Prediction 2.** At every completion range the ICU difference sits 0.6 to
  5.3 pp above the isolated two-opportunity value. So in the full ICU model the
  reversal extends to lower completion than the isolated bridge suggests: the
  e = .8 crossing lies between the 60→70% and 70→80% ranges, not between
  70→80% and 80→90%. The failing component is small (+0.25 pp) but its sign is
  unambiguous: all 20,000 draws are positive.
- **Prediction 4.**
  - Within the region of the published surface (ρ ≤ .75, all 12 e values from
    .5 to 1), the maximum error is 2.48 pp.
  - Errors grow toward high ρ and low e: 3.9 pp at ρ=.9, e=.5, and 4.3 pp at
    ρ=.95, e=.5.
  - F_geo exceeds F\* in all 132 rows. 49 rows are within 0.5 pp and 81
    within 1 pp.

## Completion ranges (F\*, %, 256 seeds)

D = F\*(ρ=.5) − F\*(ρ=0) in percentage points, ICU against the isolated
two-opportunity bridge.

| Completion | ρ=0, e=1 | ρ=.5, e=1 | ρ=0, e=.8 | ρ=.5, e=.8 | D(e=1) ICU / isolated | D(e=.8) ICU / isolated |
|---|---|---|---|---|---|---|
| 60%→70% | 47.4 | 36.8 | 33.9 | 32.9 | -10.6 / -13.4 | -0.9 / -3.5 |
| 70%→80% | 57.7 | 43.2 | 37.1 | 37.4 | -14.5 / -17.1 | +0.3 / -2.4 |
| 80%→90% | 75.7 | 57.1 | 42.4 | 46.8 | -18.6 / -20.8 | +4.4 / +1.2 |
| 90%→95% | 75.1 | 53.6 | 28.0 | 37.0 | -21.5 / -22.7 | +9.0 / +5.6 |
| 95%→99% | 95.7 | 81.2 | 25.9 | 42.5 | -14.6 / -15.2 | +16.6 / +11.3 |

## Transition-defined (WHO-style) policy

- **Construction.** Each profile's per-boundary h0 and h1 are set so that
  compliance per patient-to-patient transition goes from 90% to 95%. That is,
  1 − (1−h)(ρ + (1−ρ)(1−h)) equals .90 and then .95.
  - ρ=0: h .6838 → .7764.
  - ρ=.5: h .8292 → .9084.
  - Exact values are in `summary.json`.

| e | F\*(ρ=0) | F\*(ρ=.5) | D (pp) [95%] | Isolated k=2 D |
|---|---|---|---|---|
| 1 | 52.36% | 52.48% | +0.12 [0.12, 0.12] | 0 (exact) |
| .8 | 34.18% | 41.62% | +7.45 [7.44, 7.45] | +6.11 |

- **Interpretation, descriptive.** Under perfect clearing, the 21.5 pp clustering
  effect at equal per-boundary completion vanishes almost entirely when the
  profiles are matched on transition compliance. This is expected: at e = 1 the
  transition miss probability equals q2. Under 80% clearing a clustering effect
  of +7.4 pp remains.
- So "same compliance" is not one condition. The forecast gap between clustered
  and independent misses depends on which compliance is matched. Matching
  transition compliance removes the gap only when clearing is perfect.

## Secondary: exact ICU projection of the measurement study

- Saved assay limits at 8,000 per assay, projected through bilinear
  interpolation of the computed fine-grid F\*.
- No set fell outside the domain, and none was extrapolated.

| Truth | Exact median width (pp) [10th–90th] | Approximate (round 1, same sets) | Defined / empty |
|---|---|---|---|
| Independent | 5.26 [2.69, 7.26] | 6.17 | 978 / 22 |
| Clustered | 0.98 [0.86, 1.11] | 0.92 | 985 / 15 |
| Imperfect | 6.66 [3.99, 8.22] | 6.66 | 999 / 1 |
| Combined (ρ=.5, e=.8) | 11.62 [11.15, 12.28] | 11.67 | 1000 / 0 |

- The code check reproduced the certified isolated-bridge widths (to 2e-9 pp)
  and the round-1 approximate values exactly.
- The exact projection replaces the approximation. It is a projection of saved
  measurements and makes no new coverage claim.

## Limitations

- These are fixed-horizon p → 0 limits with no finite-p error bound. The
  finite-p reversal remains untested (O-3).
- The transition-defined policy relies on the stationary Markov law and the
  model's exit/entry opportunity rule (one exit and one entry draw per move).
- The study ran on macOS, under the original workflow only. At p = 0 the
  original and alternative workflows differ by less than 0.3 pp on all 22
  archived rows (workflow-coefficient study).

## Files

`summary.json`, `records/` (256 sealed records), `code/`, `gate_checks/`,
`BENCHMARK.json`, `CODE_REVIEW.md`, `review_scratch/`.
