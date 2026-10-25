#!/usr/bin/env python3
"""Independent exact-arithmetic O12 checks; no production-code imports.

Enumerate binary completion histories, rather than evaluating the proposed
matrix formula. All comparisons use fractions.Fraction equality.
"""
from fractions import Fraction as Q
from itertools import product
from math import comb
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent


def history(k, h, rho):
    """Probability and rho derivative for each complete binary history."""
    if k == 0:
        yield (), Q(1), Q(0)
        return
    pi = [1 - h, h]
    for xs in product((0, 1), repeat=k):
        prob, deriv = pi[xs[0]], Q(0)
        for x, y in zip(xs, xs[1:]):
            dp = Q(x == y) - pi[y]
            p = pi[y] + rho * dp
            prob, deriv = prob * p, deriv * p + prob * dp
        yield xs, prob, deriv


def enumerated(k, h, rho, e):
    q, dq = Q(0), Q(0)
    coefficients = [Q(0)] * (k + 1)
    for xs, p, dp in history(k, h, rho):
        n = sum(xs)
        survival = (1 - e) ** n
        q += p * survival
        dq += dp * survival
        for j in range(n + 1):
            coefficients[j] += p * (-1) ** j * comb(n, j)
    return q, dq, coefficients + [Q(0)] * max(0, 3 - len(coefficients))


def q2ratio_coefficient(n, d):
    assert n[0] == d[0] == 1
    return n[2] - n[1] * d[1] + d[1] ** 2 - d[2]


def at_one(k, h, rho):
    return Q(1) if k == 0 else (1-h) * (rho + (1-rho)*(1-h)) ** (k-1)


def serialize(x):
    if isinstance(x, Q):
        return {"exact": str(x), "decimal": float(x)}
    raise TypeError(type(x))


def main():
    counts = {"local_derivative": 0, "small_e_coefficient": 0,
              "perfect_clearing_strict_sign": 0, "spectrum_identity": 0,
              "degenerate_zero": 0}
    hs = [Q(1,10), Q(2,5), Q(3,5), Q(4,5), Q(9,10), Q(19,20)]
    es = [Q(1,5), Q(4,5), Q(1)]
    rhos = [Q(1,100), Q(1,4), Q(1,2), Q(9,10), Q(99,100)]
    for h in hs:
        for k in range(8):
            for e in es:
                q, dq, _ = enumerated(k,h,Q(0),e)
                a, b = 1-h*e, e*e*h*(1-h)
                expected = Q(0) if k < 2 else (k-1)*b*a**(k-2)
                assert q == a**k and dq == expected
                counts["local_derivative"] += 1
    pairs = [(h0,h1) for h0 in hs for h1 in hs if h0 < h1]
    for h0,h1 in pairs:
        for k in range(8):
            for rho in rhos:
                c0 = enumerated(k,h0,rho,Q(1))[2]
                c1 = enumerated(k,h1,rho,Q(1))[2]
                independent0 = [Q(1), -k*h0, comb(k,2)*h0*h0]
                independent1 = [Q(1), -k*h1, comb(k,2)*h1*h1]
                got = -q2ratio_coefficient(c1,c0) + q2ratio_coefficient(independent1,independent0)
                S = sum((k-d)*rho**d for d in range(1,k))
                assert got == (h1-h0)*(h0+h1-1)*S
                counts["small_e_coefficient"] += 1
                delta = (1-h1)**k/(1-h0)**k - at_one(k,h1,rho)/at_one(k,h0,rho)
                if k >= 2:
                    assert delta < 0
                    counts["perfect_clearing_strict_sign"] += 1
                else:
                    assert delta == 0
                    counts["degenerate_zero"] += 1
    spectra = [{0:Q(1)}, {1:Q(1)}, {2:Q(1)}, {0:Q(1),2:Q(1)},
               {1:Q(1),6:Q(10000)}, {2:Q(1),6:Q(10000)},
               {0:Q(2),1:Q(3),2:Q(1,5),4:Q(11),7:Q(8)}]
    for h0,h1 in pairs:
        for e in es:
            for A in spectra:
                Ks, Ds = [], []
                for h in [h0,h1]:
                    terms = {k:enumerated(k,h,Q(0),e)[:2] for k in A}
                    Ks.append(sum(A[k]*terms[k][0] for k in A))
                    Ds.append(sum(A[k]*terms[k][1] for k in A))
                actual = (Ks[1]*Ds[0]-Ks[0]*Ds[1])/Ks[0]**2
                a0,a1 = 1-h0*e,1-h1*e
                t0,t1 = e*e*h0*(1-h0)/a0**2,e*e*h1*(1-h1)/a1**2
                E0 = sum(A[k]*a0**k*max(k-1,0) for k in A)/Ks[0]
                E1 = sum(A[k]*a1**k*max(k-1,0) for k in A)/Ks[1]
                assert actual == Ks[1]/Ks[0]*(t0*E0-t1*E1)
                assert E1 <= E0
                if t0 >= t1:
                    assert actual >= 0
                if t0 > t1 and any(k>=2 for k in A):
                    assert actual > 0
                counts["spectrum_identity"] += 1
    examples = []
    for A,rho in [({0:1,2:1},Q(1,2)),({1:1,6:10000},Q(1,2)),({2:1,6:10000},Q(1,100))]:
        def K(h,r): return sum(a*at_one(k,h,r) for k,a in A.items())
        independent = 1-K(Q(19,20),Q(0))/K(Q(9,10),Q(0))
        clustered = 1-K(Q(19,20),rho)/K(Q(9,10),rho)
        assert clustered > independent
        examples.append({"spectrum":A,"rho":rho,"F_independent":independent,
                         "F_clustered":clustered,"D_at_e_1":clustered-independent})
    # With A0=A2=1, the ratio-derivative numerator divided by e^2 is
    # c0+c1*e+c2*e^2. Negative discriminant and c2>0 prove D>0 at ALL e>0.
    h0,h1=Q(9,10),Q(19,20)
    b0,b1=h0*(1-h0),h1*(1-h1)
    coeff=[2*(b0-b1),-2*b0*h1+2*b1*h0,b0*h1*h1-b1*h0*h0]
    discriminant=coeff[1]**2-4*coeff[2]*coeff[0]
    assert coeff[2]>0 and discriminant<0
    assert coeff == [Q(340,4000),Q(-342,4000),Q(171,4000)]
    sources = []  # Private editorial provenance is excluded from the public package.
    result={"passed":True,"arithmetic":"exact fractions; complete binary-history enumeration",
            "counts":counts,"total":sum(counts.values()),"mixture_counterexamples":examples,
            "all_e_counterexample_numerator_coefficients":coeff,
            "all_e_counterexample_discriminant":discriminant,
            "source_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            "check_script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'checks.json').write_text(json.dumps(result,default=serialize,indent=2)+'\n')
    print(json.dumps({"passed":True,"counts":counts,"total":sum(counts.values())}))


if __name__=='__main__':
    main()
