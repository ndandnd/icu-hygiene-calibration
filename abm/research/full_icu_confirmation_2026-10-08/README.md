# Full ICU confirmation study

This completed synthetic study supports an advisor discussion draft. It does
not establish a clinical forecast or readiness for submission. Human authorship,
declarations and clinical assessment remain unconfirmed.

## What was run

| Module | Complete records | Purpose |
|---|---:|---|
| Calibration | 384 baseline runs | Extend the independent profile's baseline grid and select its probability before policy results |
| Holdout | 1,536 paired runs | Directly simulate three fixed profiles, 512 pairs each |
| Environment | 3,072 paired runs | Compare hand/site transmission enabled and disabled at three probabilities, 512 pairs per condition |
| Hand-survival coefficient | 256 histories | Integrate washing survival over transmission-free ICU histories |

All positive-transmission runs use 90 warmup days and 56 observation days.
The policy raises completion from 90% to 95%. The frozen protocol, analysis,
seed lists and execution-only resource amendment remain unchanged.

## Main findings

The holdout baseline rates are 2.829, 2.686 and 2.801 acquisitions per 1,000
patient-days for independent completion, persistent completion and imperfect
clearing. **All six declared baseline-matching interval checks fail.** A valid
calibration bracket did not establish the required match in fresh simulations.
These configurations predict prevention of 67.58%, 55.87% and 26.96%, or mean
count reductions of 3.135, 2.461 and 1.238 acquired carriers per 56-day run.
The profile contrasts cannot be interpreted as effects at an established equal
baseline.

For independent completion, full-minus-environment-disabled differences in the
fraction prevented are -10.20, -11.74 and -14.31 percentage points. The adjusted
98.333% intervals include zero at the lowest probability and exclude it at the
other two. The full model nevertheless has greater absolute prevention because
its starting burden is higher. The complete point estimates and intervals are
in the environmental report and figure.

The computed rare-transmission limits are 75.13%, 53.65% and 27.99%. The last
two approximate the measured prevention fractions, but the leading absolute
rates underestimate all three full-model holdouts. Agreement of a percentage
does not validate its underlying rates. The clean-start theorem supplies no
useful finite-probability error guarantee for these long runs.

The imperfect holdout control rate declines between observation halves at
descriptive 95% confidence. Two descriptive precision goals are narrowly missed.
These limitations remain reported; no probability, sample size or threshold was
changed after reading the outcomes.

## Read the work

- Main paper: `manuscript/confirmation/final_02/main.pdf` and `main.tex`.
- Plain-English guide: `guide/confirmation/final_02/plain_english_guide.pdf`.
- Separate hand-only mathematics: `manuscript/preconfirmation_pdf/certificate_note.pdf`.
- Earlier full tables: `manuscript/confirmation/final_02/antecedent_results.tex`, an uncompiled supporting fragment also preserved in the earlier full manuscript.
- Fresh figures and nine-condition rate comparison: `figures/confirmation_final_02/`.
- Independent computational reviews: `holdout_report/INDEPENDENT_REVIEW.md`, `environment_report/INDEPENDENT_REVIEW.md` and `bridge/RESULT_REVIEW.md`.

The final paper keeps four figures and separates the certificate theory and old
experiment tables. Earlier releases and preparation previews remain preserved.
The guide explains the same outcomes with patient counts, a worked washing
example and figure-by-figure instructions.

## Reproduce the results

See `REPRODUCIBILITY.md` for the complete release and replay instructions.
`public_release/release_manifest.json` binds the supplied raw records and sources.
The scientific record contains no newly recovered clinical or survey inputs.
The confidence intervals describe Monte Carlo uncertainty conditional on the
model and selected probabilities; they exclude uncertainty about a real ICU.
