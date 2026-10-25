"""Closed-form continuous k=2 confidence projection with directed rounding.

Separate post-hoc sensitivity of saved measurements. No random draws. See
POSTHOC_PROTOCOL.md for the analytic endpoint proof and inferential limits.
"""
from __future__ import annotations
from dataclasses import dataclass
from decimal import Decimal as D, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np


PRECISION = 60
GUARD = "0.000000000003"


def rounded(operation, mode):
    with localcontext() as context:
        context.prec = PRECISION
        context.rounding = mode
        return operation()


@dataclass(frozen=True)
class Interval:
    lo: D
    hi: D

    @staticmethod
    def point(value):
        x = value if isinstance(value, D) else D(str(value))
        return Interval(x, x)

    def __add__(self, other):
        if not isinstance(other, Interval): other = Interval.point(other)
        return Interval(rounded(lambda: self.lo+other.lo, ROUND_FLOOR),
                        rounded(lambda: self.hi+other.hi, ROUND_CEILING))

    __radd__ = __add__

    def __sub__(self, other):
        if not isinstance(other, Interval): other = Interval.point(other)
        return Interval(rounded(lambda: self.lo-other.hi, ROUND_FLOOR),
                        rounded(lambda: self.hi-other.lo, ROUND_CEILING))

    def __rsub__(self, other):
        return Interval.point(other)-self

    def __mul__(self, other):
        if not isinstance(other, Interval): other = Interval.point(other)
        pairs = [(a,b) for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return Interval(min(rounded(lambda a=a,b=b:a*b,ROUND_FLOOR) for a,b in pairs),
                        max(rounded(lambda a=a,b=b:a*b,ROUND_CEILING) for a,b in pairs))

    __rmul__ = __mul__

    def __truediv__(self, other):
        if not isinstance(other, Interval): other = Interval.point(other)
        if other.lo <= 0 <= other.hi: raise ArithmeticError("Interval denominator includes zero")
        pairs = [(a,b) for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return Interval(min(rounded(lambda a=a,b=b:a/b,ROUND_FLOOR) for a,b in pairs),
                        max(rounded(lambda a=a,b=b:a/b,ROUND_CEILING) for a,b in pairs))

    def sqrt(self):
        if self.lo < 0: raise ArithmeticError("Negative square-root argument")
        with localcontext() as context:
            context.prec = PRECISION
            low = self.lo.sqrt().next_minus() if self.lo else D(0)
            high = self.hi.sqrt().next_plus() if self.hi else D(0)
        # Exact rational checks certify the square-root enclosure independently
        # of the decimal square-root rounding convention.
        assert Fraction(low)**2 <= Fraction(self.lo)
        assert Fraction(high)**2 >= Fraction(self.hi)
        return Interval(low,high)

    def width(self):
        return rounded(lambda:self.hi-self.lo,ROUND_CEILING)

    def strings(self):
        return [str(self.lo),str(self.hi)]


ZERO, ONE = Interval.point(0), Interval.point(1)
H0, H1 = Interval.point("0.9"), Interval.point("0.95")


def imax(*values):
    return Interval(max(v.lo for v in values),max(v.hi for v in values))


def imin(*values):
    return Interval(min(v.lo for v in values),min(v.hi for v in values))


def clip(value):
    return imax(ZERO,imin(ONE,value))


def survival(h,e,rho):
    x = 1-h*e
    return x*x+rho*e*e*h*(1-h)


def prevention(e,rho):
    return 1-survival(H1,e,rho)/survival(H0,e,rho)


def rho_endpoints(e,l2,u2):
    a = (1-H0*e)*(1-H0*e)
    b = e*e*H0*(1-H0)
    return clip((l2-a)/b),clip((u2-a)/b)


def at_endpoint(e,l2,u2):
    if e.lo == 0 and e.hi == 0:
        return ZERO,ZERO
    if e.lo <= 0: raise ArithmeticError("Cannot certify positive efficacy endpoint")
    low_rho,high_rho = rho_endpoints(e,l2,u2)
    return prevention(e,low_rho),prevention(e,high_rho)


@lru_cache(maxsize=None)
def solve(limits,information):
    global PRECISION
    for precision in (60,120,240):
        PRECISION = precision
        numbers = [Interval.point(D.from_float(float(x))) for x in limits]
        l1,u1 = clip(numbers[0]-GUARD),clip(numbers[1]+GUARD)
        l2,u2 = (ZERO,ONE) if information=="q1_only" else (clip(numbers[2]-GUARD),clip(numbers[3]+GUARD))
        low_e = imax(ZERO,(1-u1)/H0,(1-u2.sqrt())/H0)
        high_e = imin(ONE,(1-l1)/H0)
        boundary = 1-H0
        if l2.lo > boundary.hi:
            high_e = imin(high_e,1-((l2-boundary)/H0).sqrt())
        elif l2.hi > boundary.lo:
            continue
        if low_e.lo > high_e.hi:
            return {"status":"empty","precision":precision,"e_low":low_e.strings(),"e_high":high_e.strings()}
        if low_e.hi > high_e.lo:
            continue
        try:
            minimum = imin(*at_endpoint(low_e,l2,u2))
            maximum = imax(*at_endpoint(high_e,l2,u2))
        except ArithmeticError:
            continue
        assert minimum.width() < D("1e-40") and maximum.width() < D("1e-40")
        if minimum.lo < 0:
            minimum = imax(ZERO,minimum)
        if maximum.hi > 1:
            maximum = imin(ONE,maximum)
        return {"status":"feasible","precision":precision,"e_low":low_e.strings(),"e_high":high_e.strings(),
                "minimum":minimum.strings(),"maximum":maximum.strings()}
    return {"status":"unresolved","reason":"Could not certify feasibility/endpoint computation within240 digits"}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def float_outer(decimal_value,direction):
    value = float(D(decimal_value))
    return float(np.nextafter(value,direction))


def width_stats(values):
    finite = values[np.isfinite(values)]
    return {"defined":len(finite),"undefined":len(values)-len(finite),
            "quantiles_10_50_90":np.quantile(finite,[.1,.5,.9]).tolist() if len(finite) else None}


def main():
    here = Path(__file__).resolve().parent
    out = here/"run_001"
    out.mkdir(exist_ok=False)
    prior = here.parent/"run_001"
    source = json.loads((prior/"results.json").read_text())
    archive = np.load(prior/"raw_benchmark.npz",allow_pickle=False)
    all_certificates = {}; raw = {}; summaries = []
    largest_enclosure = D(0); largest_float_padding = D(0); max_precision = 0
    n_grid_enclosure_checks = 0; unresolved_total = 0
    for case in source["cases"]:
        prefix = f"{case['truth']}_n{case['n_per_assay']}"
        limits = archive[f"{prefix}_confidence_limits"]
        truth_f = case["true_predictions"][2]
        row = {"truth":case["truth"],"n_per_assay":case["n_per_assay"],"information_sets":{}}
        bounds_by_information = []
        for information in ("q1_only","q1_and_q2"):
            certificates = [solve(tuple(lims),information) for lims in limits]
            all_certificates[f"{prefix}_{information}"] = certificates
            bounds = np.full((len(limits),2),np.nan)
            status = np.zeros(len(limits),dtype=np.int8)
            for i,certificate in enumerate(certificates):
                if certificate["status"]=="unresolved":
                    status[i]=2;unresolved_total+=1;continue
                max_precision=max(max_precision,certificate["precision"])
                if certificate["status"]=="empty":
                    status[i]=1;continue
                bounds[i]=[float_outer(certificate["minimum"][0],-np.inf),
                           float_outer(certificate["maximum"][1],np.inf)]
                for name in ("minimum","maximum"):
                    lo,hi=map(D,certificate[name])
                    largest_enclosure=max(largest_enclosure,hi-lo)
                largest_float_padding=max(largest_float_padding,
                    D(certificate["minimum"][0])-D.from_float(bounds[i,0]),
                    D.from_float(bounds[i,1])-D(certificate["maximum"][1]))
            old_bounds=archive[f"{prefix}_{information}_bounds"][:,2,:]
            old_defined=np.isfinite(old_bounds[:,0])
            assert np.all(status[old_defined]==0)
            assert np.all(bounds[old_defined,0]<=old_bounds[old_defined,0]+2e-12)
            assert np.all(bounds[old_defined,1]>=old_bounds[old_defined,1]-2e-12)
            n_grid_enclosure_checks+=int(old_defined.sum())*2
            covered=(status==0)&(bounds[:,0]<=truth_f)&(truth_f<=bounds[:,1])
            width=bounds[:,1]-bounds[:,0];old_width=old_bounds[:,1]-old_bounds[:,0]
            jointly_defined=np.isfinite(width)&np.isfinite(old_width)
            summary={"certified_empty":int(np.sum(status==1)),"unresolved":int(np.sum(status==2)),
                "covered":int(covered.sum()),"total":len(limits),"coverage":float(covered.mean()),
                "continuous_width":width_stats(width),"grid_width":width_stats(old_width),
                "median_width_increase_among_both_defined":float(np.median(width[jointly_defined]-old_width[jointly_defined])),
                "maximum_width_increase_among_both_defined":float(np.max(width[jointly_defined]-old_width[jointly_defined])),
                "grid_empty_continuous_nonempty":int(np.sum(~old_defined&(status==0)))}
            row["information_sets"][information]=summary
            raw[f"{prefix}_{information}_bounds"]=bounds
            raw[f"{prefix}_{information}_status"]=status
            raw[f"{prefix}_{information}_covered"]=covered
            bounds_by_information.append(bounds)
        single,joint=bounds_by_information
        ok=np.isfinite(joint[:,0]);assert np.all(np.isfinite(single[ok]))
        assert np.all(joint[ok,0]>=single[ok,0]-2e-12) and np.all(joint[ok,1]<=single[ok,1]+2e-12)
        summaries.append(row)
        print(json.dumps({"completed":prefix,"joint_empty":row['information_sets']['q1_and_q2']['certified_empty']}),flush=True)
    np.savez_compressed(out/"continuous_bounds.npz",**raw)
    with gzip.open(out/"decimal_certificates.json.gz","wt",encoding="utf-8") as stream:
        json.dump(all_certificates,stream,separators=(",",":"))
    result={"scope":"Post-hoc continuous-domain k=2 sensitivity; no new draws; h0=.90,h1=.95",
        "conditional_coverage":"Exact confidence-set projection theorem over continuous stationary Markov rho/e; numeric optimization is outward-certified for supplied widened assay limits. CP quantiles themselves are inherited numeric inputs.",
        "guard_on_assay_limits":GUARD,"largest_decimal_extremum_enclosure_width":str(largest_enclosure),
        "largest_outward_float_padding":str(largest_float_padding),"maximum_decimal_precision":max_precision,
        "unresolved_total":unresolved_total,"grid_endpoint_enclosure_checks":n_grid_enclosure_checks,"cases":summaries}
    (out/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    receipt={"status":"completed","python":sys.version,"numpy":np.__version__,
        "inputs":{str(p.relative_to(here.parent)):digest(p) for p in
            (here/"POSTHOC_PROTOCOL.md",Path(__file__),prior/"results.json",prior/"raw_benchmark.npz")},
        "outputs":{p.name:digest(p) for p in out.iterdir() if p.is_file()},
        "new_random_draws":0,"source_run_modified":False}
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")


if __name__=="__main__":main()
