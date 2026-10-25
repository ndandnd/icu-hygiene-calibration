"""Prospective finite-grid controlled-assay benchmark; see PROTOCOL.md.

Run with Python 3, NumPy and SciPy:
    python3 run_benchmark.py --output /absolute/path/to/fresh/run_directory
No ICU inputs or simulator are loaded. Refuses to overwrite an output directory.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.special import betaincinv


BASE_SEED = 2026101703
N_EXPERIMENTS = 1000
SAMPLE_SIZES = (500, 2000, 8000)
TRUTHS = (
    ("independent", 0.0, 1.0),
    ("clustered", 0.5, 1.0),
    ("imperfect", 0.0, 0.8),
    ("combined", 0.5, 0.8),
)
H0 = 0.90
TARGETS = [(h, k) for h in (0.925, 0.95, 0.975) for k in (2, 6)]
TOL = 2e-12


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def survival(h: float, rho, e, k: int):
    """Stationary vector recursion, using all Markov paths implicitly."""
    rho, e = np.asarray(rho), np.asarray(e)
    if k == 0:
        return np.ones(np.broadcast_shapes(rho.shape, e.shape))
    miss, complete = np.broadcast_arrays(1-h+np.zeros_like(rho), h*(1-e))
    for _ in range(k-1):
        miss, complete = (
            miss*(1-h+rho*h) + complete*(1-h)*(1-rho),
            (miss*h*(1-rho)+complete*(h+rho*(1-h)))*(1-e),
        )
    return miss+complete


def cp_interval(counts, n: int, alpha: float):
    """Equal-tailed exact Clopper-Pearson, including endpoint cases."""
    counts = np.asarray(counts)
    lo = np.zeros(counts.shape, dtype=float)
    hi = np.ones(counts.shape, dtype=float)
    positive = counts > 0
    not_all = counts < n
    lo[positive] = betaincinv(counts[positive], n-counts[positive]+1, alpha/2)
    hi[not_all] = betaincinv(counts[not_all]+1, n-counts[not_all], 1-alpha/2)
    return lo, hi


def coverage_summary(covered):
    count = int(np.sum(covered))
    lo, hi = cp_interval(np.array([count]), len(covered), .05)
    return {"covered": count, "total": len(covered), "fraction": count/len(covered),
            "monte_carlo_95_interval": [float(lo[0]), float(hi[0])]}


def width_summary(width):
    finite = width[np.isfinite(width)]
    return {"defined": len(finite), "undefined": len(width)-len(finite),
            "quantiles_10_50_90": np.quantile(finite, [.1, .5, .9]).tolist()
            if len(finite) else [None, None, None]}


def projected_bounds(mask, prediction):
    has_candidate = mask.any(axis=1)
    lower = np.empty((len(mask), prediction.shape[1]))
    upper = np.empty_like(lower)
    for j in range(prediction.shape[1]):
        lower[:, j] = np.where(mask, prediction[:, j], np.inf).min(axis=1)
        upper[:, j] = np.where(mask, prediction[:, j], -np.inf).max(axis=1)
        lower[~has_candidate, j] = np.nan
        upper[~has_candidate, j] = np.nan
        # Extrema must be witnessed by a retained candidate; this checks every row.
        assert np.all(np.any(mask[has_candidate] &
                             (np.abs(prediction[:, j]-lower[has_candidate, j, None]) <= TOL), axis=1))
        assert np.all(np.any(mask[has_candidate] &
                             (np.abs(prediction[:, j]-upper[has_candidate, j, None]) <= TOL), axis=1))
    return lower, upper


def heterogeneous_survival(h: float, k: int):
    rho, e = .5, .8
    weight = (1-h)/(1-h+rho*h)
    return weight*(1-e*(1-rho)*h)**k+(1-weight)*(1-e)**k


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    here = Path(__file__).resolve().parent
    started = datetime.now(timezone.utc).isoformat()
    rho, e = np.meshgrid(np.arange(41)/40, np.arange(101)/100, indexing="ij")
    rho, e = rho.ravel(), e.ravel()
    q1 = survival(H0, rho, e, 1)
    q2 = survival(H0, rho, e, 2)
    assert np.max(np.abs(q1-(1-H0*e))) < TOL
    assert np.max(np.abs(q2-(q1*q1+e*e*rho*H0*(1-H0)))) < TOL
    predictions = np.column_stack([
        1-survival(h, rho, e, k)/survival(H0, rho, e, k)
        for h, k in TARGETS
    ])
    assert np.all(np.isfinite(predictions))
    assert np.min(predictions) >= -TOL and np.max(predictions) <= 1+TOL
    alternative = np.array([
        1-heterogeneous_survival(h, k)/heterogeneous_survival(H0, k)
        for h, k in TARGETS
    ])
    for h in [H0]+[target[0] for target in TARGETS]:
        for k in (1, 2):
            assert abs(heterogeneous_survival(h, k)-float(survival(h, .5, .8, k))) < TOL
    target_labels = [f"h={h:.3f},k={k}" for h, k in TARGETS]
    archive = {"rho_grid": rho, "e_grid": e, "q1_grid": q1, "q2_grid": q2,
               "prediction_grid": predictions, "alternative_predictions": alternative}
    rows = []
    for truth_index, (name, true_rho, true_e) in enumerate(TRUTHS):
        indices = np.flatnonzero((abs(rho-true_rho) <= TOL) & (abs(e-true_e) <= TOL))
        assert len(indices) == 1
        index = indices[0]
        truth_prediction = predictions[index]
        for n in SAMPLE_SIZES:
            rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([BASE_SEED, truth_index, n])))
            count1 = rng.binomial(n, q1[index], N_EXPERIMENTS)
            count2 = rng.binomial(n, q2[index], N_EXPERIMENTS)
            lo1, hi1 = cp_interval(count1, n, .025)
            lo2, hi2 = cp_interval(count2, n, .025)
            single = (q1[None, :] >= lo1[:, None]-TOL) & (q1[None, :] <= hi1[:, None]+TOL)
            joint = single & (q2[None, :] >= lo2[:, None]-TOL) & (q2[None, :] <= hi2[:, None]+TOL)
            assert np.all(~joint | single)
            prefix = f"{name}_n{n}"
            archive.update({f"{prefix}_counts": np.column_stack([count1, count2]),
                            f"{prefix}_confidence_limits": np.column_stack([lo1, hi1, lo2, hi2])})
            row = {"truth": name, "rho": true_rho, "e": true_e, "n_per_assay": n,
                   "true_q1": float(q1[index]), "true_q2": float(q2[index]),
                   "true_predictions": truth_prediction.tolist(), "information_sets": {}}
            bound_sets = []
            for label, mask in (("q1_only", single), ("q1_and_q2", joint)):
                low, high = projected_bounds(mask, predictions)
                defined = mask.any(axis=1)
                covered = defined[:, None] & (low-TOL <= truth_prediction) & (truth_prediction <= high+TOL)
                information = {"empty_sets": int(np.sum(~defined)),
                               "candidate_count_quantiles_10_50_90": np.quantile(mask.sum(axis=1), [.1,.5,.9]).tolist(),
                               "parameter_coverage": coverage_summary(mask[:, index]),
                               "simultaneous_prediction_coverage": coverage_summary(covered.all(axis=1)),
                               "targets": []}
                for j, target in enumerate(target_labels):
                    information["targets"].append({"target": target,
                        "coverage": coverage_summary(covered[:, j]),
                        "width": width_summary(high[:, j]-low[:, j])})
                if name == "combined":
                    alternative_covered = defined[:, None] & (low-TOL <= alternative) & (alternative <= high+TOL)
                    information["non_markov_alternative_no_coverage_guarantee"] = {
                        "targets": [{"target": target, "coverage": coverage_summary(alternative_covered[:, j])}
                                    for j, target in enumerate(target_labels)],
                        "simultaneous_prediction_coverage": coverage_summary(alternative_covered.all(axis=1))}
                row["information_sets"][label] = information
                archive[f"{prefix}_{label}_bounds"] = np.stack([low, high], axis=-1)
                archive[f"{prefix}_{label}_candidate_counts"] = mask.sum(axis=1)
                archive[f"{prefix}_{label}_parameter_covered"] = mask[:, index]
                archive[f"{prefix}_{label}_prediction_covered"] = covered
                bound_sets.append((low, high))
            (low1, high1), (low2, high2) = bound_sets
            valid = np.isfinite(low1) & np.isfinite(low2)
            assert np.all(low2[valid] >= low1[valid]-TOL)
            assert np.all(high2[valid] <= high1[valid]+TOL)
            width1, width2 = high1-low1, high2-low2
            row["contraction"] = []
            for j, target in enumerate(target_labels):
                ok = valid[:, j]
                delta = width1[ok,j]-width2[ok,j]
                row["contraction"].append({"target": target, "both_defined": int(np.sum(ok)),
                    "strictly_narrower_count": int(np.sum(delta > TOL)),
                    "median_absolute_width_reduction": float(np.median(delta)) if len(delta) else None})
            rows.append(row)
            print(json.dumps({"finished": prefix, "joint_empty": row["information_sets"]["q1_and_q2"]["empty_sets"]}), flush=True)
    archive_path = out / "raw_benchmark.npz"
    np.savez_compressed(archive_path, **archive)
    result = {"scope": "Finite-grid stationary controlled-assay confidence-set projection; no ICU or clinical validation",
              "design": {"baseline_h": H0, "grid_size": len(rho), "seed": BASE_SEED,
                         "experiments_per_case": N_EXPERIMENTS, "sample_sizes_per_assay": SAMPLE_SIZES,
                         "targets": target_labels, "assay_marginal_coverage_lower_bound": .975,
                         "joint_parameter_and_all_target_coverage_lower_bound": .95,
                         "comparison_tolerance": TOL},
              "unconstrained_grid_prediction_bounds": np.stack([predictions.min(axis=0), predictions.max(axis=0)], axis=-1).tolist(),
              "non_markov_alternative_predictions": alternative.tolist(), "cases": rows}
    result_path = out / "results.json"
    result_path.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    receipt = {"started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
               "python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__,
               "inputs": {"PROTOCOL.md": sha256(here/"PROTOCOL.md"), "run_benchmark.py": sha256(Path(__file__))},
               "outputs": {archive_path.name: sha256(archive_path), result_path.name: sha256(result_path)},
               "status": "completed", "replicate_experiments": len(rows)*N_EXPERIMENTS,
               "all_grid_predictions_finite": True, "nestedness_checked_every_experiment_target": True,
               "extrema_witness_checked_every_defined_experiment_target": True}
    (out/"receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")


if __name__ == "__main__":
    main()
