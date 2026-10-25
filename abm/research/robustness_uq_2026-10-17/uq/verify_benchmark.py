"""Independent algorithmic replay of the archived finite-grid UQ benchmark.

Uses explicit completion-word enumeration rather than the runner's vector
recursion; verifies Clopper-Pearson limits via binomial tail identities. Written
after viewing benchmark results, so this is not a blind verification.
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np
import scipy
from scipy.stats import binom


TOL = 2e-12


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def explicit_survival(h, rho, efficacy, k):
    answer = np.zeros_like(rho)
    transition = ((1-h+rho*h, h*(1-rho)),
                  ((1-h)*(1-rho), h+rho*(1-h)))
    for word in product((0, 1), repeat=k):
        probability = np.full_like(rho, (1-h, h)[word[0]])
        for previous, current in zip(word, word[1:]):
            probability *= transition[previous][current]
        answer += probability*(1-efficacy)**sum(word)
    return answer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory", type=Path)
    args = parser.parse_args()
    root = args.run_directory.resolve()
    verification_path = root/"independent_verification.json"
    if verification_path.exists():
        raise FileExistsError(verification_path)
    result = json.loads((root/"results.json").read_text())
    receipt = json.loads((root/"receipt.json").read_text())
    for file, expected in receipt["inputs"].items():
        assert digest(root.parent/file) == expected, file
    for file, expected in receipt["outputs"].items():
        assert digest(root/file) == expected, file
    data = np.load(root/"raw_benchmark.npz", allow_pickle=False)
    rho, e = data["rho_grid"], data["e_grid"]
    assert len(rho) == 4141
    q1 = 1-.9*e
    q2 = (1-.9*e)**2 + e*e*rho*.9*.1
    assert np.max(abs(q1-data["q1_grid"])) < TOL
    assert np.max(abs(q2-data["q2_grid"])) < TOL
    predictions = np.column_stack([
        1-explicit_survival(h,rho,e,k)/explicit_survival(.9,rho,e,k)
        for h in (.925,.95,.975) for k in (2,6)
    ])
    max_prediction_difference = float(np.max(abs(predictions-data["prediction_grid"])))
    assert max_prediction_difference < TOL
    # Projection membership uses recorded runner q to preserve its specified
    # exact floating boundary rule; their independent formulas agree above.
    q1, q2 = data["q1_grid"], data["q2_grid"]
    max_tail_error = 0.0
    comparisons = 0
    truth_names = ("independent", "clustered", "imperfect", "combined")
    for case in result["cases"]:
        name, n = case["truth"], case["n_per_assay"]
        prefix = f"{name}_n{n}"
        counts, limits = data[f"{prefix}_counts"], data[f"{prefix}_confidence_limits"]
        idx = np.flatnonzero((rho==case["rho"]) & (e==case["e"]))
        assert len(idx) == 1
        idx = idx[0]
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([2026101703, truth_names.index(name), n])))
        assert np.array_equal(counts[:,0], rng.binomial(n,q1[idx],1000))
        assert np.array_equal(counts[:,1], rng.binomial(n,q2[idx],1000))
        for assay in (0,1):
            successes = counts[:,assay]
            lower, upper = limits[:,2*assay], limits[:,2*assay+1]
            nonzero, nonfull = successes>0, successes<n
            errors = np.concatenate([
                abs(binom.sf(successes[nonzero]-1,n,lower[nonzero])-.0125),
                abs(binom.cdf(successes[nonfull],n,upper[nonfull])-.0125)])
            max_tail_error = max(max_tail_error,float(errors.max()))
            assert np.all(lower[~nonzero]==0) and np.all(upper[~nonfull]==1)
        for label in ("q1_only","q1_and_q2"):
            saved_bounds = data[f"{prefix}_{label}_bounds"]
            saved_counts = data[f"{prefix}_{label}_candidate_counts"]
            all_covered = []
            parameter_covered = []
            for i,(lo1,hi1,lo2,hi2) in enumerate(limits):
                valid = (lo1-TOL<=q1) & (q1<=hi1+TOL)
                if label == "q1_and_q2":
                    valid &= (lo2-TOL<=q2) & (q2<=hi2+TOL)
                assert int(valid.sum()) == saved_counts[i]
                parameter_covered.append(bool(valid[idx]))
                if not valid.any():
                    assert np.all(np.isnan(saved_bounds[i]))
                    all_covered.append(False)
                    continue
                low = predictions[valid].min(axis=0)
                high = predictions[valid].max(axis=0)
                assert np.max(abs(low-saved_bounds[i,:,0])) < TOL
                assert np.max(abs(high-saved_bounds[i,:,1])) < TOL
                covered = (low-TOL<=predictions[idx]) & (predictions[idx]<=high+TOL)
                assert np.array_equal(covered,data[f"{prefix}_{label}_prediction_covered"][i])
                all_covered.append(bool(covered.all()))
                comparisons += 12
            info = case["information_sets"][label]
            assert sum(all_covered) == info["simultaneous_prediction_coverage"]["covered"]
            assert sum(parameter_covered) == info["parameter_coverage"]["covered"]
            assert sum(saved_counts==0) == info["empty_sets"]
    assert max_tail_error < 2e-10
    verification = {"status":"passed", "not_blind":True,
        "method":"Explicit completion-word probabilities, direct per-experiment extrema, exact binomial CI tail identities, and fixed-seed redraw",
        "all_12000_experiment_seeds_and_counts_reproduced":True,
        "projected_endpoint_comparisons":comparisons,
        "max_prediction_discrepancy":max_prediction_difference,
        "max_binomial_tail_identity_error":max_tail_error,
        "prediction_tolerance":TOL,"binomial_tail_tolerance":2e-10,
        "python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__,
        "verifier_sha256":digest(Path(__file__)),"input_receipt_sha256":digest(root/"receipt.json")}
    verification_path.write_text(json.dumps(verification,indent=2)+"\n")
    print(json.dumps(verification))


if __name__ == "__main__":
    main()
