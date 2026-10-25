"""Elementary lower/upper certificates for the restricted first-generation replay.

The source is uniform over distinct visited labels. Independent perfect hygiene
occurs before and after every entry; contamination comes only from one fixed
source, and acquired patients never transmit. See research/first_generation_bounds.md.
These are model-specific probability bounds, not clinical predictions or a novel
general temporal-network theorem. The original replay module remains unchanged.
"""
from collections import Counter
from collections.abc import Iterable
import math

from first_generation import _history, _probability, first_generation_union_bound


def lower_source_average(sequence: Iterable[str], h: float,
                         alpha: float, beta: float) -> float:
    """Lower bound on expected distinct acquisitions with a uniform visited source.

    For each ordered pair of different labels, count adjacent source-to-target
    visits c. Its sufficient-event probability is 1-(1-alpha*beta*q)**c, where
    q=(1-h)**2. Those events are independent for a fixed ordered pair because
    their entries and latent Bernoulli draws cannot overlap. Other target events
    need not be independent: summing marginal lower bounds is valid.

    An unvisited-source universe must be accounted for externally by multiplying
    by visited_count/universe_count. Complexity is O(N) time and O(N) space.
    """
    history = _history(sequence)
    h, alpha, beta = (_probability(h, 'h'), _probability(alpha, 'alpha'),
                      _probability(beta, 'beta'))
    u = alpha*beta*(1-h)**2
    if u == 0.:
        return 0.
    counts = Counter((source, target) for source, target in zip(history, history[1:])
                     if source != target)
    if u == 1.:
        return len(counts)/len(set(history))
    log_failure = math.log1p(-u)
    return math.fsum(-math.expm1(count*log_failure) for count in counts.values()) / len(set(history))


def first_generation_bounds(sequence: Iterable[str], h: float,
                            alpha: float, beta: float) -> tuple[float, float]:
    """Return (lower, raw upper) bounds for the uniform visited-source expectation.

    The raw upper is the existing alpha*beta*B/M union bound and can exceed the
    target-count cap. The iterable is consumed once, so generator inputs work.
    """
    history = _history(sequence)
    return (lower_source_average(history, h, alpha, beta),
            first_generation_union_bound(history, h, alpha, beta))
