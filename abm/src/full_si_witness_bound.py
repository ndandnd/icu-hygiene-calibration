"""Finite genealogical extra-acquisition union bound for permanent full SI.

Only valid for initial clean hands, one uniform source, permanent patient
colonization, independent perfect before/after hygiene, identical alpha/beta,
transfer retaining contamination and newly positive sources at FUTURE visits.
No environmental sources, turnover or correlated hygiene. A mathematical
witness bound, not a clinical estimate. Scalar floating-point arithmetic.
"""
from collections import Counter
from itertools import combinations
import math
from numbers import Real
from order_motifs import _inputs
from first_generation import mean_source_expected_acquisitions
from permutation_exact import exact_permuted_expected_acquisitions


def _arguments(sequence,h,alpha,beta,pool):
    x,alpha,beta,pool=_inputs(sequence,alpha,beta,pool)
    if not isinstance(h,Real) or not math.isfinite(h) or not 0<=h<=1:
        raise ValueError('h must be a finite probability')
    if len(set(x))>15:
        raise ValueError('bounded implementation supports at most15 labels')
    return x,(1-h)**2,alpha*beta,pool


def genealogical_extra_bound(sequence,h,alpha,beta,source_pool_size=15,*,max_states=1000000):
    """V_ge2 over all distinct-bed, strictly nonoverlapping ancestral bridges.

    Forward state is (ancestral bed subset, final infection-entry index).
    Multiple hand acquisition/transfer witnesses are summed by a union bound.
    Source mixing weight1/pool is included once at the initial bridge.
    """
    x,q,product,pool=_arguments(sequence,h,alpha,beta,source_pool_size)
    if type(max_states) is not int or max_states<1:
        raise ValueError('positive integer state budget required')
    labels={label:i for i,label in enumerate(dict.fromkeys(x))}
    bits=[1<<labels[label] for label in x];m=len(labels);n=len(x)
    if m<3 or q==0 or product==0 or max(Counter(x).values())<2:
        return 0.
    at=[{} for _ in x];states=0
    def put(j,mask,value):
        nonlocal states
        if not value:
            return
        if mask not in at[j]:
            states+=1
            if states>max_states:
                raise RuntimeError('genealogical witness state budget exhausted')
            at[j][mask]=value
        else:
            at[j][mask]+=value
    for j in range(n):
        for i in range(j):
            if x[i]!=x[j]:
                put(j,bits[i]|bits[j],product*q**(j-i)/pool)
    contributions=[]
    for j in range(n):
        for mask,value in at[j].items():
            if mask.bit_count()>=3:
                contributions.append(value)
        # Future clean-hand acquisition must occur strictly after infection j.
        # Sum its candidate source visits, discounting only each hand bridge.
        survival_sum=0.
        for end in range(j+1,n):
            survival_sum=q*(survival_sum+(1. if end-1>j and x[end-1]==x[j] else 0.))
            edge=product*survival_sum
            if not edge:
                continue
            for mask,value in at[j].items():
                if not mask&bits[end]:
                    put(end,mask|bits[end],value*edge)
    return math.fsum(contributions)


def positional_bridge_sums(n,q,max_generation):
    """T_g=sum_(i1<j1<...<ig<jg) q**sum(j_l-i_l), all positions distinct."""
    if type(n) is not int or n<1 or type(max_generation) is not int or max_generation<0 or not 0<=q<=1:
        raise ValueError('invalid positional bridge arguments')
    closed=[1.]+[0.]*max_generation;opened=[0.]*max_generation
    for _ in range(n):
        # Old closed values ensure a target and next source cannot coincide.
        new_closed=[closed[0]]+[closed[g]+q*opened[g-1] for g in range(1,max_generation+1)]
        opened=[q*opened[g]+closed[g] for g in range(max_generation)]
        closed=new_closed
    return closed


def ancestral_label_weight(counts,generation):
    """C_g=sum distinct orderedbedpaths n_first*n_last*prod (n_mid)_2."""
    if generation<1 or generation>=len(counts):
        return 0
    total=0
    for middle in combinations(range(len(counts)),generation-1):
        factor=math.prod(counts[i]*(counts[i]-1) for i in middle)
        outside=[c for i,c in enumerate(counts) if i not in middle]
        total+=factor*(sum(outside)**2-sum(c*c for c in outside))
    return math.factorial(generation-1)*total


def exact_null_genealogical_extra_bound(sequence,h,alpha,beta,source_pool_size=15):
    """Exact uniform-multiset mean of V_ge2, not exact mean acquisition.

    Label probabilities C_g/(N)_(2g) and position sums T_g factor because
    every selection of2g distinct positions has the same multiset law.
    """
    x,q,product,pool=_arguments(sequence,h,alpha,beta,source_pool_size)
    counts=list(Counter(x).values());n=len(x);maximum=min(len(counts)-1,n//2)
    if maximum<2 or q==0 or product==0:
        return 0.
    times=positional_bridge_sums(n,q,maximum)
    return math.fsum(product**g*ancestral_label_weight(counts,g)*times[g]/math.prod(range(n-2*g+1,n+1))/pool
        for g in range(2,maximum+1))


def full_si_interval(sequence,h,alpha,beta,source_pool_size=15,*,random_order=False):
    """Return verified first-generation F <= full SI <= min(countcap,F+V)."""
    x,q,product,pool=_arguments(sequence,h,alpha,beta,source_pool_size)
    if random_order:
        first=exact_permuted_expected_acquisitions(x,h,alpha,beta,source_pool_size=pool)
        extra=exact_null_genealogical_extra_bound(x,h,alpha,beta,pool)
    else:
        first=mean_source_expected_acquisitions(x,h,alpha,beta)*len(set(x))/pool
        extra=genealogical_extra_bound(x,h,alpha,beta,pool)
    m=len(set(x));cap=m*(m-1)/pool
    return dict(lower=first,upper=min(cap,first+extra),extra_union_bound=extra,source_pool_size=pool)
