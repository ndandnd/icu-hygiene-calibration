"""Exact first-generation assay with stationary Markov hygiene failures.

The chain runs across every before/after opportunity. Its marginal failure
probability is 1-h and its lag-one autocorrelation is rho. Correlation is an
assumed mechanism, not an empirically estimated parameter. Four substochastic
states track worker contamination and the previous wash outcome while a target
remains susceptible. Patients acquired during the assay never transmit.
"""
from collections import Counter
from collections.abc import Iterable
from functools import lru_cache
import math
from numbers import Integral

from first_generation import _history, _probability


def _validated(sequence, h, alpha, beta, rho, source_pool_size):
    history = _history(sequence)
    values = tuple(_probability(value, name) for value, name in
                   [(h, 'h'), (alpha, 'alpha'), (beta, 'beta'), (rho, 'rho')])
    labels = tuple(dict.fromkeys(history))
    pool = len(labels) if source_pool_size is None else source_pool_size
    if not isinstance(pool, Integral) or isinstance(pool, bool) or pool < len(labels):
        raise ValueError('source pool must be an integer at least the visited-label count')
    return history, labels, int(pool), values


def failed_run_probability(h: float, rho: float, opportunities: int) -> float:
    """Stationary probability that a specified run contains only failures."""
    h, rho = _probability(h, 'h'), _probability(rho, 'rho')
    if not isinstance(opportunities, Integral) or isinstance(opportunities, bool) or opportunities < 1:
        raise ValueError('opportunities must be a positive integer')
    epsilon = 1-h
    return epsilon*(rho+(1-rho)*epsilon)**(opportunities-1)


@lru_cache(maxsize=16)
def _kernels(h, alpha, beta, rho):
    # Index = 2*worker_dirty + previous_wash_failed. A successful wash sends
    # either hand state to index0. Failed washes retain the current hand state.
    pf = ((1-rho)*(1-h), rho+(1-rho)*(1-h))

    def wash(mass):
        result = [0.]*4
        for state, value in enumerate(mass):
            fail = pf[state % 2]
            result[0] += value*(1-fail)
            result[2*(state//2)+1] += value*fail
        return result

    kernels = []
    for kind in range(3):  # source, target, other
        rows, losses = [], []
        for state in range(4):
            mass = [float(i == state) for i in range(4)]
            mass = wash(mass)
            loss = 0.
            if kind == 0:
                for last in (0, 1):
                    acquisition = alpha*mass[last]
                    mass[last] -= acquisition
                    mass[2+last] += acquisition
            elif kind == 1:
                loss = beta*(mass[2]+mass[3])
                mass[2] *= 1-beta
                mass[3] *= 1-beta
            rows.append(tuple(wash(mass)))
            losses.append(loss)
        kernels.append((tuple(rows), tuple(losses)))
    return tuple(kernels)


def correlated_expected_acquisitions(sequence: Iterable[str], h: float,
                                     alpha: float, beta: float, rho: float,
                                     source_pool_size: int | None = 15) -> float:
    """Expected distinct acquisitions, uniform source pool, fixed visit order."""
    history, labels, pool, values = _validated(sequence, h, alpha, beta, rho, source_pool_size)
    h, alpha, beta, rho = values
    if len(labels) == 1 or h == 1 or alpha == 0 or beta == 0:
        return 0.
    kernels = _kernels(h, alpha, beta, rho)
    risks = []
    for source in labels:
        for target in labels:
            if source == target:
                continue
            mass = [h, 1-h, 0., 0.]
            losses = []
            for label in history:
                matrix, immediate = kernels[0 if label == source else 1 if label == target else 2]
                losses.append(math.fsum(mass[i]*immediate[i] for i in range(4)))
                mass = [math.fsum(mass[i]*matrix[i][j] for i in range(4)) for j in range(4)]
            risks.append(math.fsum(losses))
    return math.fsum(risks)/pool


@lru_cache(maxsize=8)
def _category_evaluator(h, alpha, beta, rho):
    kernels = _kernels(h, alpha, beta, rho)

    @lru_cache(maxsize=None)
    def risk(a, b, o):
        if b == 0:
            return (0.,)*4
        total = a+b+o
        terms = [[] for _ in range(4)]
        counts = (a, b, o)
        for kind, count in enumerate(counts):
            if not count:
                continue
            remaining = list(counts)
            remaining[kind] -= 1
            child = risk(*remaining)
            matrix, immediate = kernels[kind]
            for state in range(4):
                terms[state].append(count/total*(immediate[state] +
                    math.fsum(matrix[state][j]*child[j] for j in range(4))))
        return tuple(math.fsum(values) for values in terms)
    return risk


def correlated_permuted_expected_acquisitions(sequence: Iterable[str], h: float,
                                              alpha: float, beta: float, rho: float,
                                              source_pool_size: int | None = 15) -> float:
    """Exact uniform multiset-order mean using source/target/other counts.

The hygiene chain is independent of the random order. Source/target identities
with equal visit counts share the same four-state remaining-count recursion.
The state reduction requires first-generation transmission; it does not apply
unchanged when newly acquired patients can transmit.
"""
    history, labels, pool, values = _validated(sequence, h, alpha, beta, rho, source_pool_size)
    h, alpha, beta, rho = values
    if len(labels) == 1 or h == 1 or alpha == 0 or beta == 0:
        return 0.
    evaluate = _category_evaluator(h, alpha, beta, rho)
    frequencies = Counter(Counter(history).values())
    terms = []
    for a, source_multiplicity in frequencies.items():
        for b, target_multiplicity in frequencies.items():
            pairs = source_multiplicity*(target_multiplicity-int(a == b))
            if pairs:
                vector = evaluate(a, b, len(history)-a-b)
                terms.append(pairs*(h*vector[0]+(1-h)*vector[1]))
    return math.fsum(terms)/pool


def correlated_union_bound(sequence: Iterable[str], h: float, alpha: float,
                           beta: float, rho: float, source_pool_size: int | None = 15) -> float:
    """Unclipped source-acquisition/transfer witness union bound."""
    history, _, pool, values = _validated(sequence, h, alpha, beta, rho, source_pool_size)
    h, alpha, beta, rho = values
    witnesses = [failed_run_probability(h, rho, 2*(j-i))
                 for i in range(len(history)) for j in range(i+1, len(history))
                 if history[i] != history[j]]
    return alpha*beta*math.fsum(witnesses)/pool
