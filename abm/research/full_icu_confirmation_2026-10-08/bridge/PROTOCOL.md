# Fixed full-ICU p=0 coefficient study

8 October 2026. This separate coefficient study follows the completed response study and its critique. Its budget, inputs and partitions are fixed before the full-horizon p=0 batch. It uses the independently reviewed bridge implementation and the frozen scientific model. It adds no positive-p simulation and does not alter the response or confirmation study.

## Estimand and fixed budget

Run **256 fresh epidemiological backgrounds**, statistical bank A with 128 seeds and bank B with 128 disjoint seeds. The exact seed lists are in `bridge/seedbook.json` and the already frozen parent `seed_audit.json`. They are disjoint from the prior seed audit and the 1,216 declared primary/optional family seeds. The audit covers accessible JSON seed fields in the ABM and parent research trees; external absent archives are not certified. Every background evaluates all three `(rho,e)=(0,1),(.5,1),(0,.8)` profiles and both schedules together by exact conditional integration.

All backgrounds have **3% admission carriage, initially clean hands/sites, p=0, 90×2592 warmup ticks and 56×2592 observation ticks**. Completion is h0=.90 throughout warmup and control; the policy changes to h1=.95 after the final warmup tick. Use `zone_separate`, the frozen 30-bed/staffing/topology/demographic model, and histogram cap 2: exact contributions at 0,1,2 opportunities and an exact aggregate >=3 tail. Other fixed model settings are those in the independently reviewed `background.py`. Do not use its 5%/90+28 defaults: every executable full-batch config contains the explicit overrides above.

The estimand is the expectation over zero-transmission demographic/contact backgrounds of the imported-carrier-to-susceptible same-worker bridge coefficient K_h. The recurrence retains initialized/persistent completion memory, pre-screening isolation weights, warmup source contacts, role resets, admission identities, contact episodes and actual opportunity order. It integrates hygiene histories conditional on a sampled background; floating-point arithmetic and sampling over backgrounds remain. Clean sites inject no environmental two-success mass, although site contacts influence opportunities. The source and policy coefficients share each background.

The primary limiting response is **1−mean(K1)/mean(K0)**, undefined if mean(K0)=0. Report K0, K1, their paired difference, the ratio of means, limiting response, warmup-source contributions, opportunity contributions, exposure and event counts. Display the stationary q2 response as a labelled comparison. It is not automatically the ICU coefficient. Do not substitute the mean of per-seed ratios or exclude zero-coefficient seeds.

The reviewed runner reports seed-level Monte Carlo SE and nominal 95% normal/delta intervals. These are descriptive asymptotic intervals, not exact coverage claims, clinical-input uncertainty, or the primary confirmation study's baseline-equivalence gates. Pool all 256 backgrounds for the primary estimate after completeness/source checks. Show both statistical banks and their B−A estimates as Monte Carlo diagnostics, preserving within-background profile/arm covariance and independence between banks. A bank difference is not a compute-site effect; absence of a detectable difference does not prove bank or RNG equivalence. Any later bootstrap must resample whole background vectors, retain zero/undefined ratios and declare its analysis source/seed; it may not change this simulation budget based on outcomes.

No sample extension, replacement of valid seeds, longer horizon, positive-p convergence series, or mechanism refit is authorized by this protocol. Wide intervals and undefined responses are study outcomes.

## Portable fixed partitions and restart

Allocate one eighth of each statistical bank locally: **32 local seeds =16A+16B**; **224 cluster seeds =112A+112B**. Local seed order is A[:16]+B[:16], cluster order A[16:]+B[16:]. Partition each site list deterministically by `site_seeds[index::partition_count]`. There are **6 local partitions** (sizes 6,6,5,5,5,5) and **32 cluster partitions** (7 seeds each). Both statistical banks occur at both sites and in every partition. Statistical bank membership is independent of the technical site/partition labels.

A partition has its own config and output scope:

- `bridge/configs/local_part_00.json` through `local_part_05.json` → `results/full_icu_confirmation_2026-10-08/bridge_local/part_00` through `part_05`.
- `bridge/configs/cluster_part_00.json` through `cluster_part_31.json` → `results/full_icu_confirmation_2026-10-08/bridge_cluster/part_00` through `part_31`.

`bridge_all.json`, `bridge_local.json` and `bridge_cluster.json` document the complete scientific/site seed lists; the partition wrapper executes only its declared partition config. The exact config/path/output membership is in `bridge/seedbook.json`. Public receipts contain only ABM-root-relative paths.

Before launch, the parent creates **freeze_bridge.json** with the existing `freeze.py --phase bridge`, binding this protocol, seedbook, configs, wrapper, reviewed bridge sources/tests/review, initial freeze and existing science chain. The outer wrapper calls the parent `array_worker.verify_freeze`, requires a bridge-phase freeze with all needed binding/scientific files, validates exact science/partition membership and primary seed non-overlap, then writes immutable `bridge_binding.json` in that partition scope. It rechecks the binding before invoking the backend and after return. Inner runtime/source/config manifests belong to the unchanged `bridge/run.py`; existing output directories cannot be resumed under a changed runtime/platform, source, config or seed partition.

Call a partition only after the parent freeze exists, for example:

```
.venv/bin/python research/full_icu_confirmation_2026-10-08/bridge/array_worker.py --site local --task-index 0 --freeze research/full_icu_confirmation_2026-10-08/freeze_bridge.json --wall-seconds 3300
```

The wrapper enforces a finite wall budget no larger than **3300 seconds**. Slurm SIGUSR1 and SIGTERM raise a caught interruption and return **75**. A seed interrupted before immutable publication restarts from initialization; saved complete seeds validate and skip on resume. An interruption during an atomic temporary write can leave a temporary file, which the backend replaces on retry; it does not create a complete seed record. There is no biological checkpoint. Each partition's existing backend `.run.lock` prevents duplicate active writers. Exit 75 requires resubmission of the same partition/config/output, not a new seed or changed budget. Hard kills cannot guarantee an exit code but retain the same immutable-output/restart rule.

Use one CPU per cluster partition. Runtime/memory requests and retry scheduling are technical choices checked before production; they do not alter seeds, scientific settings or n. Completed shard summaries must be combined from individual seed rows, not by averaging unequal-sized shard responses or confidence intervals. Verify complete exact 256-seed membership, common science/source hashes, both outer/inner bindings and legitimate site runtime identities before pooled analysis. Runtime differences are recorded rather than erased.

## Interpretation and limits

The result computes K for this fixed 90+56-day, clean-reservoir stochastic ICU and its ideal independent-draw law. The legacy RNG conditionally consumes draws; common distributional background laws do not imply alternative-profile literal seeded replay or an eventwise positive-p coupling. Source-conditional behavior at positive p, saturation, isolation, competing routes and reservoir feedback remain outside the coefficient.

The expansion is A_h(p)=p²K_h+O(p³) with horizon fixed. Its remainder depends on potential transmission slots and reward bounds. Neither this computation nor small p alone supplies a useful finite-p relative-error bound or establishes stationarity after warmup. Environmental last-pickup proportions are proximal descriptions, not a lower bound on approximation error or a causal attributable fraction. A later comparison to measured positive-p forecasts must retain their baseline residuals and finite-horizon interpretation.

The independent implementation signoff and bounded nine-test suite are in `INDEPENDENT_REVIEW.md`. Preparing this protocol, configurations and wrapper does not freeze or launch the full batch. No completed K estimate is assumed by the manuscript until the separately frozen batch is finished and checked.
