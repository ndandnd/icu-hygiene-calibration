# Clustered hand-hygiene lapses and decontamination efficacy in calibrated ICU transmission models

Nathan Cho · Cornell Tech · nc437@cornell.edu

Code and synthetic data release **v1.0.2**: research software, plots and **17,920 completed synthetic simulation records** supporting the study named above. These are simulator-generated histories and numerical outputs, not patient records. Funding: none. Conflicts of interest: none.

The study tests how clustered hand-hygiene lapses and incomplete clearing affect transmission under two ICU workflow rules. All failed checks and predictions remain in the deposited results: the original finite-transmission holdout failed all six baseline-matching interval checks; two of the extension's four prespecified predictions failed; its finite-transmission reversal remains untested. The conclusions are conditional model results, not clinical validation. Full publication documents are preserved outside this repository.

## Technical report and detailed tables

The [technical reproducibility report](TECHNICAL_REPRODUCIBILITY_REPORT.pdf) retains the complete technical record, including the detailed calibration checks, historical sweeps, coefficient grids and uncertainty qualifications. It is supporting documentation, separate from the main paper and shorter journal supplement. The main paper and core supplement PDFs are reserved for a preprint/journal deposit.

- [Full 132-point coefficient grid, CSV](technical_report/companion_data/coefficient_grid132.csv)
- [Original row fields and full precision, JSON](technical_report/companion_data/coefficient_grid132_original_fields.json)
- [Grid provenance and units](technical_report/companion_data/PROVENANCE.json)
- [Rebuild the technical report](technical_report/source/README.txt)

All old synthetic simulation banks are unchanged. The report states the assumption that hand contamination alone leaves movement, screening and exposure unchanged before acquisition, and explains why irrelevant-event cancellations do not supply a horizon-independent error bound.

## Start here

From this repository root, use Python 3.12. The following checks require only the standard library and do not run simulations:

```sh
python -B tools/exact_checks.py
python -B abm/research/referee_perspective_2026-10-24/worker_robustness/replay_pooled.py --output worker_replay.json
python -B abm/research/manuscript_v2_2026-10-21/supplement_numbers/make_supplement_tables.py
python -B abm/research/manuscript_v2_2026-10-21/supplement_numbers/make_extension_tables.py
```

The first command verifies 1,659 exact rational-arithmetic identities and counterexamples. The second independently replays eight worker profiles, 16 coefficient boxes and 24 conservative contrasts from the saved pooled coefficients; it does not revalidate their raw-record provenance. The table builders read the deposited machine-readable numerical sources and assert the reported values. They write generated outputs beside their scripts. Existing saved scientific outputs should be preserved if byte comparisons are needed.

The clean local candidate passed the exact checker, pooled replay and 237-number supplement table generation on 24 October 2026. `REPRODUCTION_CHECKS.json` records the specific validation scope. Full simulations and all long bootstrap analyses were not rerun for this deposit.

## Synthetic data

Raw banks are compressed into independently extractable archives, each below 90 MB. `DATA_INDEX.json` lists every study, extraction path, archive hash and member hash. There are no duplicate transport archives, remote execution logs or old Git histories. Verify all compressed data without writing extracted banks:

```sh
python -B tools/unpack_data.py --verify-only
```

Extract all banks, or choose one or more study names from the table:

```sh
python -B tools/unpack_data.py --study coefficient_extension_2026-10-23
python -B tools/unpack_data.py
```

Extraction verifies compressed hashes and every member hash, rejects unsafe paths and links, and refuses to overwrite different existing files. Extraction expands several GB of JSON.

| Study name | Completed records | Scientific use |
|---|---:|---|
| original_confirmation | 5,248 | Calibration 384; holdout 1,536; environment 3,072; bridge histories 256 |
| original_clock | 1,536 | Earlier clock sensitivity |
| sensitivity | 1,536 | Earlier sensitivity analysis |
| extended_bridge | 256 | Original-workflow 22-pair coefficient surface |
| rule_clock | 3,456 | Complete canonical grid and policy banks from the latest policy collection |
| rule_surfaces_off | 1,536 | Environmental pathway sensitivity |
| rule_boundary | 512 | Clock boundary sensitivity |
| rule_upper | 256 | Independent upper-bound study |
| alternative_finite_p | 3,072 | Alternative-workflow finite-transmission grid and policy banks |
| workflow_coefficient_2026-10-21 | 256 | Same-background comparison between workflows |
| coefficient_extension_2026-10-23 | 256 | Expanded 182-pair coefficient study and measurement projection |

Older study files keep their scientific directory names. Some are called `author_notes` or `review_...` because those were their historical analysis locations; only scientific source/numerical files are included there, not private correspondence or review conversations. The `review_.../inference` implementation is an independent numerical implementation, not a prose review.

## Numerical reanalysis and figures

Create an environment and install the direct scientific requirements:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

After extracting the extension bank, its independent numerical implementation recomputes the 182 pooled estimates, 20,000 common-seed bootstrap draws, four predictions and assay projection:

```sh
python -B abm/research/coefficient_extension_2026-10-23/review_2026-10-23/inference/independent_inference.py
```

This is a longer numerical computation, not an ICU simulation. It writes results beside the script. Use `python -B abm/research/manuscript_v2_2026-10-21/figures/make_figures.py` to rebuild the publication plots from the deposited numerical sources. This repository contains code, synthetic data, generated plots/tables and [technical execution extracts](execution_records/INDEX.md). Complete paper and Supplementary Information PDFs and full manuscript source snapshots are excluded. The author plans separate preprint dissemination; a preprint URL will be added when available.

## Full simulation source and scope

`abm/src/` contains the maintained scientific modules and the exact three frozen model files used by the studies. The original confirmation, response, bridge, workflow variant, expanded coefficient, mathematical, UQ, common-baseline and follow-up analyses are included under `abm/research/`. Seedbooks, resolved configurations, manifests, source hashes and saved results accompany the relevant studies. `abm/requirements-lock.txt` preserves the historical environment snapshot; `requirements.txt` lists the minimal direct packages for scientific reanalysis. Historical lock entries are preserved for provenance, not a promise that every future platform can install every historical version.

For example, the expanded coefficient runner can inspect its plan with `python -B abm/research/coefficient_extension_2026-10-23/code/run_batch.py --plan`. A fresh simulation requires a **new output directory**, chosen seeds/workers and the matching frozen source hashes. The original positive-transmission runner is `abm/research/full_icu_confirmation_2026-10-08/array_worker.py`; its `--help` documents the configuration, output and partition arguments. Runs use the deposited study-specific configurations and seedbooks. Positive-transmission confirmation runs use 90 warmup days and 56 observation days. Runtime depends strongly on workload and hardware; consult each deposited `BENCHMARK.json`/runtime report and protocol before scheduling the complete banks. The expanded runner documents per-seed timeout controls. No short runtime or full end-to-end rerun is claimed.

Some immutable records contain historical absolute provenance paths. They are not needed for the portable proof/table/pooled replay commands. Frozen source checks deliberately preserve byte identities. Historical full-provenance analyzers may need those references mapped to the corresponding extracted `DATA_INDEX.json` member paths; they should not be used by silently editing the immutable records. Current figure/table scripts use repository-relative discovery. `REPRODUCTION_SCOPE.md` distinguishes the portable checks from these archival workflows.

See `RIGHTS.md` for the source notice policy; no new blanket license is asserted for upstream or third-party code. `SHA256SUMS` binds every delivered regular file except itself and subsequently added Git metadata. Regenerating outputs changes their hashes.
