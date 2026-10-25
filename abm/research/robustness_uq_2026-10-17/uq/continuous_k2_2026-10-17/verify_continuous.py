"""Independent 100-digit scalar verification of continuous-k2 certificates.

Uses the equivalent direct numerator formula for F and scalar boundary formulas,
not the interval implementation. Written after the analytic derivation, without
using new assay draws. This numerical verification is not a new theorem about
the inherited Clopper-Pearson beta quantiles.
"""
import gzip
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np


mp.mp.dps = 100
H0,H1 = mp.mpf(9)/10,mp.mpf(19)/20
GUARD = mp.mpf(3)/10**12
SCALAR_CHECK_TOLERANCE = mp.mpf("1e-90")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scalar_projection(limits,information):
    a,b,c,d = (mp.mpf(float(x)) for x in limits)
    l1,u1 = max(0,a-GUARD),min(1,b+GUARD)
    l2,u2 = (mp.mpf(0),mp.mpf(1)) if information=="q1_only" else (max(0,c-GUARD),min(1,d+GUARD))
    lo = max(0,(1-u1)/H0,(1-mp.sqrt(u2))/H0)
    hi = min(1,(1-l1)/H0)
    if l2>1-H0:hi=min(hi,1-mp.sqrt((l2-(1-H0))/H0))
    if lo>hi:return None
    def candidates(e):
        if e==0:return [mp.mpf(0),mp.mpf(0)]
        a0=(1-H0*e)**2;b0=e**2*H0*(1-H0)
        rlo=max(0,min(1,(l2-a0)/b0));rhi=min(1,max(0,(u2-a0)/b0))
        return [(H1-H0)*e*(2-(rho+(1-rho)*(H0+H1))*e)/(a0+rho*b0)
                for rho in (rlo,rhi)]
    return min(candidates(lo)),max(candidates(hi))


def main():
    here=Path(__file__).resolve().parent
    out=here/"run_001"
    target=out/"independent_verification.json"
    if target.exists():raise FileExistsError(target)
    source=json.loads((here.parent/"run_001/results.json").read_text())
    data=np.load(here.parent/"run_001/raw_benchmark.npz",allow_pickle=False)
    current=np.load(out/"continuous_bounds.npz",allow_pickle=False)
    receipt=json.loads((out/"receipt.json").read_text())
    for path,sha in receipt["inputs"].items():assert digest(here.parent/path)==sha
    for path,sha in receipt["outputs"].items():assert digest(out/path)==sha
    with gzip.open(out/"decimal_certificates.json.gz","rt") as stream:certificates=json.load(stream)
    defined,empty=0,0
    max_difference=mp.mpf(0)
    for case in source["cases"]:
        prefix=f"{case['truth']}_n{case['n_per_assay']}"
        limits=data[f"{prefix}_confidence_limits"]
        rho,e=mp.mpf(str(case['rho'])),mp.mpf(str(case['e']))
        truth_f=(H1-H0)*e*(2-(rho+(1-rho)*(H0+H1))*e)/((1-H0*e)**2+rho*e**2*H0*(1-H0))
        for information in ("q1_only","q1_and_q2"):
            key=f"{prefix}_{information}"
            for i,lim in enumerate(limits):
                scalar=scalar_projection(lim,information)
                cert=certificates[key][i]
                if scalar is None:
                    assert cert["status"]=="empty" and current[key+"_status"][i]==1
                    empty+=1
                    continue
                assert cert["status"]=="feasible" and current[key+"_status"][i]==0
                defined+=1
                bounds=current[key+"_bounds"][i]
                for name,value in zip(("minimum","maximum"),scalar):
                    lo,hi=map(mp.mpf,cert[name])
                    assert lo-SCALAR_CHECK_TOLERANCE<=value<=hi+SCALAR_CHECK_TOLERANCE,(key,i,name,str(lo-value),str(value-hi))
                    max_difference=max(max_difference,abs(value-(lo+hi)/2))
                assert mp.mpf(float(bounds[0]))<=scalar[0]
                assert mp.mpf(float(bounds[1]))>=scalar[1]
                exact_covered=scalar[0]<=truth_f<=scalar[1]
                assert bool(exact_covered)==bool(current[key+"_covered"][i])
    report={"status":"passed","method":"100-digit independent scalar formulas and alternative rational numerator for prevention",
        "decimal_working_precision":100,"defined_sets_checked":defined,"empty_sets_checked":empty,
        "scalar_verification_tolerance":str(SCALAR_CHECK_TOLERANCE),
        "extremum_certificate_checks":2*defined,"max_difference_from_decimal_certificate_midpoint":str(max_difference),
        "floating_outputs_enclose_scalar_extrema":True,"coverage_checked_against_exact_decimal_truth_parameters":True,
        "python":sys.version,"mpmath":mp.__version__,"verifier_sha256":digest(Path(__file__)),
        "source_receipt_sha256":digest(out/"receipt.json"),
        "limitation":"Does not certify inherited beta quantiles; it verifies continuous optimization for supplied widened numeric limits."}
    target.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report))


if __name__=="__main__":main()
