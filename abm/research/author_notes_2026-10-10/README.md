# ICU author-note follow-up

This separate study addresses useful author notes in an externally rewritten
draft. It preserves the completed 8 October confirmation study, including its
failed baseline-matching outcome. No manuscript source or old result changed.

| Module | Fixed computation | Purpose |
| --- | --- | --- |
| Original clock | 12 conditions x 128 fresh baseline seeds | Test sampled acquisition rates under before/after-contact washing at p=.15,.30,.60,1 |
| Probability sensitivity | 6 conditions x 256 fresh control/policy pairs | Measure response changes at +/-25% of the two existing fitted probabilities |
| Coefficient sensitivity | Same 256 old p=0 backgrounds, 22 coefficient pairs | Generalize hand-survival calculations and describe contacts/bridges by role |
| Saved-data diagnostics | Complete prior synthetic records | Weekly trajectory, calibration-selection uncertainty, rate-slope uncertainty and coarse bridge counts |

Read `PROTOCOL.md`, `design.json` and `DESIGN_REVIEW.md` before interpreting
results. `freeze.json` binds the new executable/configuration closure. Runtime
and job identifiers are held in the local research journal, outside public
source. The remote arrays share at most 48 workers (24/16/8 respectively), use
atomic completed records, and requeue only retryable interruption/timeouts.
They stop after 12 restarts or three consecutive retries without progress.

The new positive-p banks are written to `results/author_notes_2026-10-10/clock`
and `sensitivity`; coefficient replays go to `bridge/part_*`. Each production
module has a fixed budget. No new fitting, condition selection, automatic sample
extension, clinical recommendation or equilibrium claim is permitted.

## Saved-data results already completed

`saved_analysis/outputs` contains separately hashed JSON/CSV results and a
weekly trajectory PNG. These are retrospective descriptive diagnostics:

- Independent surface-disabled log slope: 2.0275, 95% interval [1.6697,2.4180].
- Surface-enabled minus surface-disabled log slope: 4.3311 [3.7469,5.1028].
  Both are three-point finite-probability fits, not established asymptotic laws.
- Independent calibration selection: p=.143276, conditional 95% interval
  [.140045,.147920]. Pooled resampling had no bracket failures; bank A had
  1.48%, explicitly retained. This does not repair holdout baseline matching.
- The old backgrounds contain 107,763,665 unweighted, reset-surviving bridges.
  Their coarse histogram is 0:10, 1:0, 2:3,158,292, >=3:104,605,363. Their
  hygiene-weighted contribution is a different quantity.
- Observation contains 87,088,900 washing opportunities over 419,933.0876
  patient-days, or 207.3876 opportunities per patient-day. This is a modeled
  clock intensity, not comparable compliance data without harmonized counting.

The 146-day weekly plot labels the week crossing the warmup boundary and the
final six-day week. It does not establish stationarity. The new bridge replay
adds fine histograms and role-specific counters absent from old summaries.

## Reproduction and continuation

Use the pinned Mesa 1.2.1 / NumPy 1.26.4 research environment. Verify all frozen
bytes before resuming. A fresh runtime is checked against five short positive-p
numerical reference cases. Bridge replays reconcile each original seed's event
hash, exposure and three original coefficient pairs before publishing it.
`analyze_positive.py` refuses incomplete banks and resamples whole paired seed
vectors independently within statistical banks. `bridge/VALIDATION.md` records
independent tiny tests and bounded performance measurements.

Parameter provenance, literature comparisons, repository DOI and confirmed
human declarations require documentary or author work. Compute cannot supply
those facts. The advisor discussion draft remains the prior canonical document
until the new results have been reviewed and their qualifications preserved.
