"""Second-order hygiene expansion for the restricted first-generation assay.

All outputs are conditional mathematical quantities, not clinical estimates.
One fixed uniform source, clean initial hands, independent perfect hygiene,
constant source/transfer probabilities, no secondary/environmental transmission.
The expectation R(q)=a*q+b*q**2+O(q**3), q=(1-h)**2, is polynomial for
a finite sequence; this module returns its first two coefficients, not a
validated second-order approximation outside the stated finite certificates.
"""
from collections import Counter
import math
from numbers import Integral, Real


def _history(sequence):
    if isinstance(sequence, (str, bytes)):
        raise ValueError('sequence must be an iterable of nonempty string labels')
    try:
        history = tuple(sequence)
    except TypeError as exc:
        raise ValueError('sequence must be an iterable of labels') from exc
    if not history or any(not isinstance(x, str) or not x.strip() for x in history):
        raise ValueError('sequence must contain nonempty string labels')
    return history


def _inputs(sequence, alpha, beta, source_pool_size):
    history = _history(sequence)
    for value in (alpha, beta):
        if not isinstance(value, Real) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError('alpha and beta must be finite probabilities')
    if (not isinstance(source_pool_size, Integral) or isinstance(source_pool_size, bool)
            or source_pool_size < len(set(history))):
        raise ValueError('source pool must be an integer at least the visited label count')
    return history, float(alpha), float(beta), int(source_pool_size)


def motif_counts(sequence):
    """Counts of short different-bed motifs and repeated ordered adjacencies."""
    history = _history(sequence)
    adjacent = Counter(zip(history, history[1:]))
    triples = zip(history, history[1:], history[2:])
    lag2 = repeated_source = repeated_target = 0
    for x, y, z in triples:
        lag2 += x != z
        repeated_source += x == y and x != z
        repeated_target += x != y and y == z
    return {
        'adjacent_different': sum(v for (x, y), v in adjacent.items() if x != y),
        'lag2_different': lag2,
        'source_source_target': repeated_source,
        'source_target_target': repeated_target,
        'repeated_adjacent_pairs': sum(v*(v-1)//2 for (x, y), v in adjacent.items() if x != y),
    }


def exact_null_motif_means(sequence):
    """Exact uniform-multiset expectations; degenerate histories give zeros."""
    history = _history(sequence)
    n = len(history)
    if n < 2:
        return {k: 0. for k in motif_counts(history)}
    counts = tuple(Counter(history).values())
    falling = [c*(c-1) for c in counts]
    d = 1-sum(falling)/(n*(n-1))
    repeated = sum(c*(c-1)*(n-c) for c in counts)/(n*(n-1)) if n >= 3 else 0.
    repeated_pairs = (sum(falling)**2-sum(v*v for v in falling))/(2*n*(n-1))
    return {
        'adjacent_different': (n-1)*d,
        'lag2_different': (n-2)*d if n >= 3 else 0.,
        'source_source_target': repeated,
        'source_target_target': repeated,
        'repeated_adjacent_pairs': repeated_pairs,
    }


def hygiene_coefficients(sequence, alpha, beta, source_pool_size=15, *, random_order=False):
    """Return (a,b) in R(q)=a*q+b*q**2+O(q**3), raw acquisitions.

    For the fixed order, b decomposes into lag-two hand bridges, source
    saturation, target absorption on s,t,t, and repeated adjacent-target
    absorption. Random-order coefficients use exact multiset motif moments.
    Unvisited source cases contribute zero through the declared denominator.
    """
    history, alpha, beta, pool = _inputs(sequence, alpha, beta, source_pool_size)
    counts = exact_null_motif_means(history) if random_order else motif_counts(history)
    product = alpha*beta
    a = product*counts['adjacent_different']/pool
    b = (product*(counts['lag2_different']-alpha*counts['source_source_target']
                  -beta*counts['source_target_target'])
         -product**2*counts['repeated_adjacent_pairs'])/pool
    return a, b


def second_order_matched_gap(recorded_second, null_second, q0, q1, baseline_difference=0.):
    """Finite-q second-order prediction, not a guaranteed error-controlled gap.

    The exact baseline difference is recorded minus null. Parameters are held
    fixed across policy; remainders contribute an unquantified O(q**3) term
    unless a separate finite certificate is evaluated. q0 must be positive.
    """
    values = (recorded_second, null_second, q0, q1, baseline_difference)
    if any(not isinstance(v, Real) or not math.isfinite(v) for v in values):
        raise ValueError('coefficients and probabilities must be finite')
    if not 0 < q0 <= 1 or not 0 <= q1 <= q0:
        raise ValueError('require 0 <= q1 <= q0 <= 1 with q0 positive')
    lam = q1/q0
    return ((null_second-recorded_second)*(q1*q1-lam*q0*q0)
            +(1-lam)*baseline_difference)
