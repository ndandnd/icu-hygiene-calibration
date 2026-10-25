"""Exact uniform-multiset permutation expectation for first-generation carriage.

This module shares the restricted mechanism of first_generation.py: a clean
initial worker, one fixed colonized source, initially susceptible targets, no
secondary transmission, no environment, and independent perfect hygiene before
and after every entry. Transfer does not clear worker contamination. It evaluates
uniform random order conditional on each record's bed frequencies, not clinical
risk or a biological calibration. Original order affects the recorded assay but
is intentionally irrelevant to this counterfactual.
"""
from collections import Counter
from collections.abc import Iterable
from functools import lru_cache
import math
from numbers import Integral, Real


def _probability(value, name):
    if not isinstance(value, Real) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f'{name} must be a finite probability in [0,1]')
    return float(value)


@lru_cache(maxsize=8)
def _category_evaluator(h: float, alpha: float, beta: float):
    """Memoize clean/dirty future-risk pairs for this parameter family.

    At most eight parameter families are retained. Each family's state cache
    contains only the remaining-count triples reached by callers; count-class
    aggregation lets records and target pairs reuse those states.
    """
    q = (1-h)**2
    source_dirty_from_clean = alpha*(1-h)
    source_dirty_from_dirty = (1-h)*((1-h)+h*alpha)
    target_immediate_from_dirty = beta*(1-h)
    target_remaining_dirty = q*(1-beta)
    target_remaining_clean = 1-target_immediate_from_dirty-target_remaining_dirty

    @lru_cache(maxsize=None)
    def risk(a: int, b: int, o: int) -> tuple[float, float]:
        # The target is still susceptible on entry to this state. If no target
        # visits remain, ever-acquisition risk is zero for either hand state.
        if b == 0:
            return 0., 0.
        total = a+b+o
        clean_terms, dirty_terms = [], []
        if a:
            child_clean, child_dirty = risk(a-1, b, o)
            clean_terms.append(a/total*((1-source_dirty_from_clean)*child_clean
                                       + source_dirty_from_clean*child_dirty))
            dirty_terms.append(a/total*((1-source_dirty_from_dirty)*child_clean
                                       + source_dirty_from_dirty*child_dirty))
        if b:
            child_clean, child_dirty = risk(a, b-1, o)
            clean_terms.append(b/total*child_clean)
            dirty_terms.append(b/total*(target_immediate_from_dirty
                                       + target_remaining_clean*child_clean
                                       + target_remaining_dirty*child_dirty))
        if o:
            child_clean, child_dirty = risk(a, b, o-1)
            clean_terms.append(o/total*child_clean)
            dirty_terms.append(o/total*((1-q)*child_clean+q*child_dirty))
        # Sum immediate acquisition, not one minus near-unit survival mass.
        clean_risk, dirty_risk = math.fsum(clean_terms), math.fsum(dirty_terms)
        return min(1., max(0., clean_risk)), min(1., max(0., dirty_risk))
    return risk


def exact_permuted_expected_acquisitions(sequence: Iterable[str], h: float,
                                        alpha: float, beta: float,
                                        source_pool_size: int | None = 15) -> float:
    """Exact expected distinct acquisitions under a uniform multiset shuffle.

    `source_pool_size` gives the equally likely source-bed pool. Its size must
    be at least the distinct visited-bed count; an unvisited source causes zero
    acquisition. None uses the distinct visited beds as the source pool, matching
    first_generation.mean_source_expected_acquisitions. The default 15 is the
    illustrative wing pool used by the recovered-record pilot.

    For fixed source/target counts (a,b), all other visits are equivalent under
    the no-secondary mechanism. Category sequences with counts (a,b,N-a-b) are
    uniformly sampled by the original full-label multiset shuffle. A memoized
    dynamic program returns target-ever-acquired risk for clean/dirty initial
    hands. The public result uses clean hands and sums target marginals divided
    by the source pool, without assuming target independence. Equal frequency
    classes have multiplicities m_a*(m_b - 1[a=b]). Complexity per target-count
    triple is at most O((a+1)*(b+1)*(o+1)); states are shared across count classes
    and records. The implementation is recursive and intended for finite shift
    histories rather than unrestricted very long sequences.
    """
    if isinstance(sequence, (str, bytes)):
        raise ValueError('sequence must contain bed labels, not a single string')
    try:
        history = tuple(sequence)
    except TypeError as exc:
        raise ValueError('sequence must be an iterable of bed labels') from exc
    if not history or any(not isinstance(label, str) or not label.strip() for label in history):
        raise ValueError('sequence must contain at least one nonempty string bed label')
    h, alpha, beta = (_probability(h, 'h'), _probability(alpha, 'alpha'),
                      _probability(beta, 'beta'))
    counts = Counter(history)
    pool = len(counts) if source_pool_size is None else source_pool_size
    if not isinstance(pool, Integral) or isinstance(pool, bool) or pool < len(counts):
        raise ValueError('source_pool_size must be an integer at least the distinct visited-bed count, or None')
    if len(counts) == 1 or h == 1 or alpha == 0 or beta == 0:
        return 0.
    frequencies = Counter(counts.values())
    evaluate = _category_evaluator(h, alpha, beta)
    contributions = []
    for a, source_multiplicity in frequencies.items():
        for b, target_multiplicity in frequencies.items():
            pairs = source_multiplicity*(target_multiplicity-int(a == b))
            if pairs:
                target_risk, _ = evaluate(a, b, len(history)-a-b)
                contributions.append(pairs*target_risk)
    return math.fsum(contributions)/int(pool)
