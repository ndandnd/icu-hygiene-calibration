"""Independent exact checks of proposed results; no ICU data or simulator edits.

Run from any directory with Python 3 and SymPy. Outputs a dated research receipt
alongside this file. All numerical examples are constructed, not observations.
"""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib
import json
import math
import platform
import sys

import sympy as sp


def survival(h, rho, e, k):
    """Enumerate every stationary completion word, independent of closed forms."""
    if k == 0:
        return Q(1)
    pi = (1-h, h)
    P = ((1-h+rho*h, h*(1-rho)),
         ((1-h)*(1-rho), h+rho*(1-h)))
    ans = Q(0)
    for word in product((0, 1), repeat=k):
        probability = pi[word[0]]
        for i in range(1, k):
            probability *= P[word[i-1]][word[i]]
        ans += probability * (1-e)**sum(word)
    return ans


def encode(x):
    if isinstance(x, Q):
        return {"exact": str(x), "decimal": float(x)}
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [encode(v) for v in x]
    return x


def main():
    checks = 0
    for h, rho, e in product((Q(1,10), Q(1,2), Q(9,10), Q(19,20)),
                              (Q(0), Q(1,5), Q(1,2), Q(1)),
                              (Q(1,5), Q(1,2), Q(4,5), Q(1))):
        q1, q2 = (survival(h, rho, e, k) for k in (1,2))
        assert q1 == 1-h*e
        assert q2 == q1*q1 + e*e*rho*h*(1-h)
        recovered_e = (1-q1)/h
        recovered_rho = (q2-q1*q1)/(recovered_e**2*h*(1-h))
        assert (recovered_e, recovered_rho) == (e, rho)
        checks += 3

    h0, h1 = Q(9,10), Q(19,20)
    interaction = []
    for e in (Q(4,5), Q(1)):
        Fs = [1-survival(h1,rho,e,2)/survival(h0,rho,e,2)
              for rho in (Q(0), Q(1,2), Q(1))]
        a0, a1 = (1-h0*e)**2, (1-h1*e)**2
        b0, b1 = e*e*h0*(1-h0), e*e*h1*(1-h1)
        derivative_numerator = b1*a0-b0*a1
        for left, right in zip(Fs,Fs[1:]):
            assert (right-left)*derivative_numerator < 0
            checks += 1
        interaction.append({"efficacy":e,"rho":[0,Q(1,2),1],"F":Fs,
                            "ratio_derivative_numerator":derivative_numerator})
    t0, t1 = math.sqrt(float(h0*(1-h0))), math.sqrt(float(h1*(1-h1)))
    threshold = (t0-t1)/(float(h1)*t0-float(h0)*t1)

    # Same raw weighted count, mean bridge length AND baseline coefficient.
    ks = (1,2,3,4)
    spectra = ((Q(11,42),Q(3,28),Q(1,2),Q(11,84)),
               (Q(5,21),Q(11,28),Q(0),Q(31,84)))
    examples = []
    for A in spectra:
        K0 = sum(a*(1-h0)**k for a,k in zip(A,ks))
        K1 = sum(a*(1-h1)**k for a,k in zip(A,ks))
        assert sum(A) == 1
        assert sum(a*k for a,k in zip(A,ks)) == Q(5,2)
        assert K0 == Q(1111,40000)
        examples.append({"A":A,"K0":K0,"K1":K1,"F":1-K1/K0})
        checks += 3
    assert examples[0]["F"] != examples[1]["F"]
    checks += 1
    M = sp.Matrix([[1 for k in ks],list(ks),[sp.Rational(1,10)**k for k in ks]])
    qpolicy = sp.Matrix([[sp.Rational(1,20)**k for k in ks]])
    assert M.rank() == 3 and M.col_join(qpolicy).rank() == 4
    for delta in M.nullspace():
        assert (qpolicy*delta)[0] != 0
    checks += 2

    # Enumerate all vertices of the one-moment probability polytope. Its
    # extreme distributions use at most two support points. Independent
    # rational enumeration checks the adjacent-point/end-point sharp bounds.
    moment_examples = []
    for lo,hi in ((0,4),(1,7),(2,6)):
        for mu in (Q(lo),Q(hi),Q(lo+hi,2),Q(3*lo+hi,4)):
            for ratio in (Q(1,2), Q(6,7)):
                values=[]
                for a in range(lo,hi+1):
                    if Q(a)==mu:
                        values.append(ratio**a)
                    for b in range(a+1,hi+1):
                        if a<=mu<=b:
                            mass_a=(b-mu)/(b-a)
                            values.append(mass_a*ratio**a+(1-mass_a)*ratio**b)
                m=math.floor(mu)
                frac=mu-m
                lower=(1-frac)*ratio**m+frac*ratio**(m+1)
                upper=(hi-mu)/(hi-lo)*ratio**lo+(mu-lo)/(hi-lo)*ratio**hi
                assert min(values)==lower and max(values)==upper
                checks += 1
                if (lo,hi,mu,ratio)==(2,6,Q(4),Q(1,2)):
                    moment_examples.append({"support":[lo,hi],"weighted_mean":mu,
                                            "F_bounds":[1-upper,1-lower]})

    report={"status":"all_checks_passed","assertion_groups":checks,
            "scope":"Exact constructed stationary-bridge examples; no clinical data, no ICU validation, no novelty claim.",
            "python":sys.version,"platform":platform.platform(),"sympy":sp.__version__,
            "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "interaction":interaction,"clearing_threshold_e":threshold,
            "same_baseline_raw_mean_counterexample":examples,
            "sharp_weighted_mean_examples":moment_examples}
    output=Path(__file__).with_name("verification.json")
    output.write_text(json.dumps(encode(report),indent=2)+"\n")
    print(json.dumps({"status":report["status"],"assertion_groups":checks,
                      "threshold_e":threshold,"output":str(output)}))


if __name__ == "__main__":
    main()
