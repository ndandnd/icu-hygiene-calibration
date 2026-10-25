# Rare-transmission coefficient under the alternative workflow: results

Computed 22 October 2026 (batch finished 05:41 EDT; `summary.json` written 12:11 UTC).
Protocol: `PROTOCOL.md`, written 21 October before any outcome. This is a local
study. Nothing frozen was edited.

## Status of gates and execution

- **Gate 1 (reproduction): passed.**
  - With the alternative rules off, the adapter reproduced the archived
    original-workflow records **bitwise** on 8 seeds: 4 local and 4 cluster,
    4 in bank A and 4 in bank B. Two of them also matched all 22 surface pairs.
  - The Gate 3 reviewer reproduced 2 further seeds with an independent script.
  - Evidence: `gate_checks/GATES.json` and `review_scratch/`.
- **Gate 2 (rule semantics): passed.**
  - 14 unit tests passed; the opt-in full-horizon test was also run and passed.
  - The reviewer ran a lockstep test against the frozen robustness model for
    the full 378,432 ticks on 2 seeds. There was no state divergence.
- **Gate 3 (independent code review): approve with minor notes**
  (`CODE_REVIEW.md`).
  - All notes are operational or reporting notes. None affects a coefficient.
- **Budget rule.** The benchmark projected 15.7 CPU-hours against an 80-hour
  cap, so the full surface ran: 256 seeds × 22 pairs.
  - Actual cost: 14.3 CPU-hours, 10 workers, 1 h 46 min wall.
  - 256 of 256 records completed, with 0 failures and no seed replacement.
- **Runtime.** macOS arm64, CPython 3.12.2, NumPy 1.26.4, Mesa 1.2.1. The
  original records were produced on the cluster (Linux, 3.12.3).
  - In Gate 1 and the reviewer's rerun, 10 seeds were reproduced bitwise in
    original mode on macOS. Five of them were originally produced on the
    cluster: 4018965857, 2050894151, 154101125, 2312492055 and 576264056.
    Only those five are cross-platform evidence.
- **Determinism check.** The 3 benchmark seeds were recomputed in the batch.
  Their science fields are identical to the benchmark records.
- **Disclosure.** Before launch, the implementation agent saw 3-seed outcomes
  while testing `summarize.py` on the benchmark records. Nothing in the
  protocol, settings or code changed afterwards, and the budget decision was
  mechanical.
- **Independent recomputation (parent session).** F* for all 22 pairs was
  recomputed directly from the raw records of both workflows. It agrees with
  `summary.json` to every printed digit.

## Primary profiles (h .90 → .95, 256 shared seeds)

| Profile | Original K0 → K1 | Alternative K0 → K1 | F* original | F* alternative | ΔF* (pp) [95% paired bootstrap] |
|---|---|---|---|---|---|
| Independent (ρ=0, e=1) | 55.51 → 13.80 | 53.37 → 13.25 | 75.13% | 75.18% | +0.051 [+0.019, +0.090] |
| Clustered (ρ=.5, e=1) | 409.13 → 189.65 | 393.16 → 182.22 | 53.65% | 53.65% | +0.004 [+0.000, +0.009] |
| Imperfect (ρ=0, e=.8) | 464.91 → 334.77 | 447.05 → 321.94 | 27.99% | 27.99% | −0.008 [−0.012, −0.004] |

- Under the alternative workflow, K falls by about 4% in both arms and in every
  row. For example, independent K0 changes by −3.9% [−6.8, −0.8].
- Because control and policy fall in proportion, the ratio, and therefore F*,
  is essentially unchanged.
- The tiny paired intervals arise because F* depends on a background mainly
  through how many hygiene opportunities its bridges contain (next section).
  That mix barely varies across seeds.
- Bank A and bank B agree: the B−A difference in F* is at most 0.002 pp in
  either workflow.

## Why the limit is insensitive to the workflow

At leading order, K_h is a sum over pickup-to-transfer bridges, each weighted by
its hygiene-survival probability over the opportunities it contains. The
pooled share of K by opportunity count (k) is:

| Profile, arm | Original k=0 / 2 / ≥3 | Alternative k=0 / 2 / ≥3 |
|---|---|---|
| Independent, control | 0.07% / 98.95% / 0.98% | 0 / 99.02% / 0.98% |
| Independent, policy | 0.28% / 99.47% / 0.25% | 0 / 99.75% / 0.25% |
| Clustered, control | 0.01% / 73.83% / 26.16% | 0 / 73.93% / 26.07% |
| Imperfect, control | 0.01% / 92.62% / 7.37% | 0 / 92.68% / 7.32% |

(k=1 mass is 0 throughout.)

- **The mix is almost the same in both workflows.** Leading-order bridges are
  dominated by k=2 bridges, so F* sits close to the stationary two-opportunity
  value (q2 column below).
- **What the alternative workflow removes.** It has no k=0 (hygiene-insensitive)
  leading-order bridges. Those occur only in the original workflow, and are
  rare even there.
- **Per-seed variation.** The independent profile's per-seed K1/K0 has SD
  0.00003 in the alternative workflow, against 0.0034 in the original, whose
  range is .2481–.2786. Some original seeds carry a k=0 component; no
  alternative seed does.
- **The largest surface change.** It is at h .95→.99, ρ=0, e=1: 95.74% →
  96.01%, which moves onto the q2 value of 96.0%. Here K1 is smallest, so the
  removed k=0 mass matters most.
- **This is descriptive.** The workflow comparison changes vacancy cleaning and
  the isolation rule together. At p=0, site clearing cannot inject mass, so
  any difference arises through the isolation rule's effect on contact
  histories, but this study does not isolate it causally.

## Full surface (Table S5 rows; 256 shared seeds)

The q2 column is the stationary two-opportunity comparator, shown for
comparison only.

| h0→h1 | ρ | e | F* original (%) | F* alternative (%) | ΔF* (pp) [95% paired bootstrap] | K0 relative change (%) [95%] | q2 comparator (%) |
|---|---|---|---|---|---|---|---|
| 0.90→0.95 | 0.0 | 1.0 | 75.13 | 75.18 (iid_perfect) | +0.051 [+0.019, +0.090] | -3.9 [-6.8, -0.8] | 75.0 |
| 0.90→0.95 | 0.0 | 0.9 | 42.60 | 42.61 | +0.001 [-0.005, +0.007] | -3.8 [-6.8, -0.8] | 41.8 |
| 0.90→0.95 | 0.0 | 0.8 | 27.99 | 27.99 (iid_imperfect) | -0.008 [-0.012, -0.004] | -3.8 [-6.8, -0.8] | 26.5 |
| 0.90→0.95 | 0.0 | 0.6 | 14.93 | 14.92 | -0.005 [-0.010, -0.001] | -3.9 [-6.8, -0.9] | 12.6 |
| 0.90→0.95 | 0.25 | 1.0 | 56.73 | 56.73 | +0.005 [-0.002, +0.014] | -3.9 [-6.8, -0.8] | 55.8 |
| 0.90→0.95 | 0.25 | 0.9 | 45.51 | 45.50 | -0.006 [-0.012, -0.000] | -3.9 [-6.8, -0.8] | 43.6 |
| 0.90→0.95 | 0.25 | 0.8 | 32.49 | 32.48 | -0.012 [-0.019, -0.005] | -3.9 [-6.8, -0.8] | 29.7 |
| 0.90→0.95 | 0.25 | 0.6 | 17.19 | 17.18 | -0.006 [-0.013, +0.001] | -3.9 [-6.8, -0.9] | 13.9 |
| 0.90→0.95 | 0.5 | 1.0 | 53.65 | 53.65 (persistent_perfect) | +0.004 [+0.000, +0.009] | -3.9 [-6.9, -0.9] | 52.3 |
| 0.90→0.95 | 0.5 | 0.9 | 47.55 | 47.55 | -0.003 [-0.009, +0.004] | -3.9 [-6.8, -0.9] | 44.5 |
| 0.90→0.95 | 0.5 | 0.8 | 37.01 | 37.00 | -0.009 [-0.019, +0.001] | -3.9 [-6.8, -0.9] | 32.1 |
| 0.90→0.95 | 0.5 | 0.6 | 20.54 | 20.53 | -0.002 [-0.012, +0.008] | -3.9 [-6.8, -0.9] | 15.1 |
| 0.90→0.95 | 0.75 | 1.0 | 52.64 | 52.64 | +0.006 [+0.003, +0.009] | -3.8 [-6.8, -0.8] | 50.8 |
| 0.90→0.95 | 0.75 | 0.9 | 49.60 | 49.61 | +0.005 [-0.001, +0.010] | -3.8 [-6.8, -0.8] | 45.0 |
| 0.90→0.95 | 0.75 | 0.8 | 42.64 | 42.64 | +0.004 [-0.008, +0.015] | -3.8 [-6.8, -0.8] | 33.9 |
| 0.90→0.95 | 0.75 | 0.6 | 27.06 | 27.07 | +0.013 [-0.003, +0.029] | -3.8 [-6.8, -0.8] | 16.2 |
| 0.80→0.90 | 0.0 | 1.0 | 75.71 | 75.72 | +0.007 [-0.001, +0.017] | -3.8 [-6.8, -0.8] | 75.0 |
| 0.80→0.90 | 0.0 | 0.8 | 42.43 | 42.42 | -0.014 [-0.021, -0.007] | -3.9 [-6.8, -0.8] | 39.5 |
| 0.80→0.90 | 0.5 | 1.0 | 57.11 | 57.11 | +0.003 [-0.002, +0.009] | -3.9 [-6.8, -0.9] | 54.2 |
| 0.95→0.99 | 0.0 | 1.0 | 95.74 | 96.01 | +0.272 [+0.108, +0.471] | -4.0 [-7.0, -1.0] | 96.0 |
| 0.95→0.99 | 0.0 | 0.8 | 25.93 | 25.92 | -0.005 [-0.008, -0.001] | -3.8 [-6.8, -0.8] | 24.9 |
| 0.95→0.99 | 0.5 | 1.0 | 81.16 | 81.18 | +0.016 [+0.006, +0.028] | -3.9 [-6.9, -0.9] | 80.8 |

- Every difference is less than 0.3 pp, and 21 of 22 are less than 0.06 pp.
- Several of the paired intervals exclude zero. With intervals this narrow,
  that is not of practical importance.

## Prediction checks (written before outcomes)

1. **The primary limits change by less than 2 pp: HOLDS.**
   - The changes are +0.05, +0.004 and −0.008 pp.
   - Against the published rounded values (75.1/53.6/28.0), the changes are
     +0.08, +0.05 and −0.01 pp.
   - Every paired 95% interval lies inside ±2 pp.
2. **The sign reversal persists: HOLDS.**
   - D(e) = F*(ρ=.5, e) − F*(ρ=0, e), at h .90→.95.
   - Alternative workflow: D(1) = −21.53 pp and D(.8) = +9.02 pp.
   - Original workflow, same seeds: −21.48 and +9.02 pp.
   - In all 20,000 bootstrap draws, D(1) < 0 and D(.8) > 0, in both workflows.
3. **The crossing e\* is within .05 of .9235: HOLDS.**
   - The crossing is found by linear interpolation of D between e=.9 and e=1.
   - Alternative workflow: e\* = .9187 [.9187, .9187].
   - Original workflow, same seeds: e\* = .9187 [.9187, .9188].
   - Distance from .9235: .0048.
   - The protocol value .9235 is the exact root of the stationary q2 crossing
     (.92353). Applying the same linear interpolation to the q2 comparator
     gives .9108. Reported for context only.

## Post hoc: comparison with the verified finite-p alternative-workflow results

These finite-p results were seen before this protocol was written, so the
comparison is post hoc and not prespecified.

| Profile | Selected p | Finite-p F [95%] (512 pairs) | Alternative limit | Finite-p − limit | Leading-order control rate / measured |
|---|---|---|---|---|---|
| Independent | .16511 | 76.5% [74.1, 78.7] | 75.18% | +1.31 pp | 0.89 / 2.90 (0.31) |
| Clustered | .08063 | 53.7% [49.7, 57.5] | 53.65% | +0.05 pp | 1.56 / 2.62 (0.59) |
| Imperfect | .07614 | 27.3% [22.2, 31.9] | 27.99% | −0.73 pp | 1.58 / 2.63 (0.60) |

- Each alternative limit lies inside the finite-p 95% interval, and all are
  within 1.4 pp.
- The leading-order term 1000·p²·mean(K_h)/D0 accounts for only 31–60% of the
  measured absolute control rate. Higher-order terms, including environmental
  routes, are therefore large in absolute terms. In this workflow they shrink
  the control and policy rates nearly in proportion. That is an observation,
  not a guaranteed property.
- Wording allowed: "descriptively, each finite-p prevented fraction was within
  1.4 pp of the alternative-workflow leading-order limit, and inside its 95%
  interval".
- Wording not allowed: that the limit "predicts", "validates" or "confirms" the
  finite-p response; any causal attribution to cleaning or to isolation; any
  error bound.

## Limitations

- The coefficient is a fixed-horizon p → 0 limit. It gives no finite-p error
  bound.
- The alternative workflow is one combined synthetic world. It does not
  represent clinical practice.
- Records ran on macOS. Bitwise cross-platform identity was shown only for
  original mode, on 5 cluster-produced seeds; 5 more seeds were reproduced
  macOS to macOS. Alternative mode was not tested across platforms.
- The coefficient intervals cover background Monte Carlo variation only.

## Files

- `summary.json`: all numbers; sha256 in `RESULTS_SHA256.txt`.
- `records/`: 256 sealed per-seed records.
- `code/`, `gate_checks/`, `CODE_REVIEW.md`, `review_scratch/`, `BENCHMARK.json`.
