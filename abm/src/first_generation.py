"""Exact first-generation acquisition in a single-HCW contact replay.

This deliberately small model is separate from the spatial ICU simulation. One
fixed source patient is colonized throughout; every other patient starts
susceptible. The worker starts clean, can acquire contamination only from that
source, retains contamination during transfer, and has independent perfect
hygiene opportunities before and after *every* sequence entry. Newly acquired
patients never transmit: all secondary transmission is disabled. No environment,
clinical disease, admission, discharge, screening, or multi-worker process exists.

The expectation is the number of distinct nonsource patients newly colonized,
not the number of transmission attempts or a clinical infection estimate. Target
outcomes can be correlated. Linearity of expectation permits their marginal
probabilities to be computed separately without assuming their independence.
"""
from collections.abc import Iterable
import math
from numbers import Real


def _probability(value: float, name: str) -> float:
    if not isinstance(value, Real) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f'{name} must be a finite probability in [0, 1]')
    return float(value)


def _history(sequence: Iterable[str]) -> tuple[str, ...]:
    if isinstance(sequence, (str, bytes)):
        raise ValueError('sequence must contain bed labels, not a single string')
    try:
        history = tuple(sequence)
    except TypeError as exc:
        raise ValueError('sequence must be an iterable of bed labels') from exc
    if not history or any(not isinstance(label, str) or not label.strip() for label in history):
        raise ValueError('sequence must contain at least one nonempty string bed label')
    return history


def _visited_label(label: str, history: tuple[str, ...], name: str) -> None:
    if not isinstance(label, str) or not label.strip() or label not in history:
        raise ValueError(f'{name} must be a nonempty bed label present in sequence')


def _target_probability(history: tuple[str, ...], source: str, target: str,
                        h: float, alpha: float, beta: float) -> float:
    # Substochastic masses: clean/dirty worker AND this target still susceptible.
    clean, dirty = 1., 0.
    removed_masses = []
    for label in history:
        clean, dirty = clean + h*dirty, (1-h)*dirty  # before-contact hygiene
        if label == source:
            acquisition = alpha*clean
            clean, dirty = clean-acquisition, dirty+acquisition
        elif label == target:
            removed_masses.append(beta*dirty)
            dirty *= 1-beta  # remove paths at this target's first acquisition
        # Other patient contacts do not change contamination or this target.
        clean, dirty = clean + h*dirty, (1-h)*dirty  # after-contact hygiene
    # In exact arithmetic this equals 1 - (clean + dirty). Accumulating the
    # removed mass avoids cancellation and spurious nonzero boundary results.
    return max(0., min(1., math.fsum(removed_masses)))


def exact_target_acquisition_probability(sequence: Iterable[str], source: str,
                                       target: str, h: float, alpha: float,
                                       beta: float) -> float:
    """Probability that a particular initially susceptible target is acquired.

    `source` and `target` must be distinct visited labels. `h` is the adherence
    probability of each independent before/after perfect-hygiene opportunity;
    `alpha` is source-to-clean-worker acquisition probability and `beta` is
    dirty-worker-to-susceptible-patient transfer probability per contact.
    Repeated visits to the same target count only its first acquisition.
    """
    history = _history(sequence)
    _visited_label(source, history, 'source')
    _visited_label(target, history, 'target')
    if source == target:
        raise ValueError('source and susceptible target must be different beds')
    h, alpha, beta = (_probability(h, 'h'), _probability(alpha, 'alpha'),
                      _probability(beta, 'beta'))
    return _target_probability(history, source, target, h, alpha, beta)


def exact_expected_acquisitions(sequence: Iterable[str], source: str,
                                h: float, alpha: float, beta: float) -> float:
    """Exact expected distinct acquired patients for one fixed visited source.

    The fixed source is never counted as newly acquired. An absent source is an
    error rather than an implied outside source. Complexity is O(M*N) for M
    distinct beds and N contacts. Patient outcomes need not be independent.
    """
    history = _history(sequence)
    _visited_label(source, history, 'source')
    h, alpha, beta = (_probability(h, 'h'), _probability(alpha, 'alpha'),
                      _probability(beta, 'beta'))
    targets = (label for label in dict.fromkeys(history) if label != source)
    return math.fsum(_target_probability(history, source, target, h, alpha, beta)
                     for target in targets)


def mean_source_expected_acquisitions(sequence: Iterable[str], h: float,
                                     alpha: float, beta: float) -> float:
    """Exact expectation under a uniformly sampled distinct visited source bed.

    Uniformity over distinct beds is an explicit pilot assumption; a bed visited
    more often is not more likely to be the source. This is not fitted source
    prevalence. Complexity is O(M**2 * N).
    """
    history = _history(sequence)
    h, alpha, beta = (_probability(h, 'h'), _probability(alpha, 'alpha'),
                      _probability(beta, 'beta'))
    labels = tuple(dict.fromkeys(history))
    return math.fsum(_target_probability(history, source, target, h, alpha, beta)
                     for source in labels for target in labels if source != target) / len(labels)


def discounted_switch_score(sequence: Iterable[str], h: float) -> float:
    """B = sum_{i<j, bed_i != bed_j} q**(j-i), with q=(1-h)**2.

    B is an elementary geometrically weighted contact-opportunity score, not a
    novel theorem or a clinical risk model. A linear-time recurrence computes
    it without enumerating contact pairs; order matters when 0 < q < 1.
    """
    history = _history(sequence)
    h = _probability(h, 'h')
    q = (1-h)**2
    total_weight = 0.
    by_label: dict[str, tuple[float, int]] = {}
    contributions = []
    for j, label in enumerate(history):
        total_weight *= q
        previous_weight, previous_position = by_label.get(label, (0., j))
        same_weight = previous_weight * q**(j-previous_position)
        contributions.append(max(0., total_weight-same_weight))
        total_weight += 1.
        by_label[label] = (same_weight+1., j)
    return math.fsum(contributions)


def first_generation_union_bound(sequence: Iterable[str], h: float,
                                 alpha: float, beta: float) -> float:
    """Return the raw bound alpha*beta*B/M for a uniformly sampled source.

    Proof under this module's assumptions: any acquisition at target visit j
    requires an earlier source-to-worker acquisition at a source visit i, followed
    by no successful wash after i, before j, or at either opportunity between
    them. These are 2*(j-i) independent opportunities, with survival q**(j-i).
    The witness event has probability at most alpha*beta*q**(j-i), since the
    worker's chance of being clean at i is at most one. A union bound over source
    and target contact pairs bounds first acquisition of each distinct target.
    Summing target marginal bounds and averaging over the M equally likely
    distinct source beds gives alpha*beta*B/M.

    Initial contamination, secondary transmission, environmental acquisition,
    unequal source weights, nonconstant probabilities, or dependent hygiene
    can invalidate this particular expression. It is an upper bound, not the
    expectation. It can exceed M-1; no artificial clipping is applied.
    """
    history = _history(sequence)
    h, alpha, beta = (_probability(h, 'h'), _probability(alpha, 'alpha'),
                      _probability(beta, 'beta'))
    return alpha*beta*discounted_switch_score(history, h) / len(set(history))
