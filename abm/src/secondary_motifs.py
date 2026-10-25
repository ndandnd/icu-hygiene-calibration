"""Onward-transmission q² coefficient in the finite permanent-SI assay.

Independent before/after perfect hygiene, q=(1-h)**2, uniform fixed initial
source; newly acquired beds contaminate clean hands on FUTURE entries only.
No environment, time-based clearance, multiple workers, or patient turnover.
This is an exact coefficient statement, not finite-q error control.
"""
from collections import Counter
from order_motifs import _inputs, _history, hygiene_coefficients


def relay_count(sequence):
    """Count ordered nonoverlapping edges s->t then t->u, all labels distinct.

    Edge starts i,j require j>=i+2. A touching triple s,t,u is excluded: its
    retained-hand path already belongs to the first-generation coefficient.
    Intermediate-bed colonization persists between the two edges.
    """
    x=_history(sequence)
    earlier=Counter()
    total=0
    for j in range(len(x)-1):
        if j>=2:
            s,t=x[j-2],x[j-1]
            if s!=t:
                earlier[(s,t)]+=1
        t,u=x[j],x[j+1]
        if t!=u:
            total+=sum(count for (s,end),count in earlier.items() if end==t and s not in (t,u))
    return total


def exact_null_relay_mean(sequence):
    """Uniform-multiset mean of relay_count; N<4 gives zero."""
    x=_history(sequence);n=len(x)
    if n<4:
        return 0.
    counts=list(Counter(x).values())
    sum_squares=sum(c*c for c in counts)
    numerator=sum(c*(c-1)*((n-c)**2-(sum_squares-c*c)) for c in counts)
    return numerator/(2*n*(n-1))


def full_si_hygiene_coefficients(sequence,alpha,beta,source_pool_size=15,*,random_order=False):
    """Return coefficients a,b for full-SI R(q)=a*q+b*q²+O(q³).

    a equals first-generation a. Additional b is (alpha*beta)² times the
    ordered two-acquisition relay count divided by the fixed source pool.
    """
    x,alpha,beta,pool=_inputs(sequence,alpha,beta,source_pool_size)
    a,b=hygiene_coefficients(x,alpha,beta,pool,random_order=random_order)
    relay=exact_null_relay_mean(x) if random_order else relay_count(x)
    return a,b+(alpha*beta)**2*relay/pool
