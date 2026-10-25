"""Exact finite-state full-SI extension of the one-worker contact assay.

Patient colonization is permanent and newly colonized beds may contaminate the
worker at future visits. A mixture over all possible initial sources is propagated
once. The mechanism has independent perfect hygiene before/after every entry,
contamination retained on transfer, and no clinical infection, clearance, turnover,
environment or concurrent staffing. It is separate from the spatial ICU model.
"""
from collections.abc import Iterable
import math
from numbers import Integral, Real
import time

import numpy as np


def full_si_expected_acquisitions(sequence: Iterable[str], h: float,
                                  alpha: float, beta: float,
                                  source_pool_size: int | None = 15, *,
                                  deadline: float | None = None) -> float:
    """Exact expected nonsource first acquisitions, uniform initial-source pool.

    There are 2**M masks for M distinct visited beds, with worker-clean/dirty
    masses c,d. Initial mass is 1/pool on each single-source visited mask; absent
    sources contribute (pool-M)/pool at the all-susceptible mask. Sum transfer
    mass into colonized masks counts each patient acquisition once. Complexity
    O(N*2**M), supported M<=15; numerical arithmetic uses float64. None selects
    the distinct visited beds as the uniform source pool.
    """
    if isinstance(sequence, (str, bytes)):
        raise ValueError('sequence must contain bed labels, not a single string')
    try:
        history = tuple(sequence)
    except TypeError as exc:
        raise ValueError('sequence must be an iterable of bed labels') from exc
    if not history or any(not isinstance(label, str) or not label.strip() for label in history):
        raise ValueError('history requires nonempty string bed labels')
    for name, value in [('h', h), ('alpha', alpha), ('beta', beta)]:
        if not isinstance(value, Real) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f'{name} must be a finite probability in [0,1]')
    labels = tuple(dict.fromkeys(history))
    visited = len(labels)
    if visited > 15:
        raise ValueError('This bounded exact evaluator supports at most15 visited beds')
    pool = visited if source_pool_size is None else source_pool_size
    if not isinstance(pool, Integral) or isinstance(pool, bool) or pool < visited:
        raise ValueError('source pool must be an integer at least the visited-bed count, or None')
    if deadline is not None and time.perf_counter() >= deadline:
        raise TimeoutError('Exact structural-test wall budget exhausted')
    if visited == 1 or h == 1 or alpha == 0 or beta == 0:
        return 0.
    count = 1 << visited
    masks = np.arange(count, dtype=np.uint32)
    clean = np.zeros(count, dtype=np.float64)
    dirty = np.zeros_like(clean)
    clean[np.left_shift(1, np.arange(visited))] = 1/int(pool)
    clean[0] = (int(pool)-visited)/int(pool)
    contact_indices = {}
    for bed, label in enumerate(labels):
        bit = 1 << bed
        susceptible = np.flatnonzero((masks & bit) == 0)
        contact_indices[label] = (susceptible, susceptible | bit)
    transfers = []
    for label in history:
        if deadline is not None and time.perf_counter() >= deadline:
            raise TimeoutError('Exact structural-test wall budget exhausted')
        clean += h*dirty
        dirty *= 1-h
        susceptible, colonized = contact_indices[label]
        # Both event masses are snapshots of the pre-contact distribution.
        # Disjoint masks distinguish clean-worker acquisition from dirty-worker
        # transfer. Newly transferred mass cannot transmit again this entry.
        acquisition = alpha*clean[colonized]
        transfer = beta*dirty[susceptible]
        clean[colonized] -= acquisition
        dirty[colonized] += acquisition
        dirty[susceptible] -= transfer
        dirty[colonized] += transfer
        transfers.append(float(transfer.sum()))
        clean += h*dirty
        dirty *= 1-h
    total_mass = float(clean.sum()+dirty.sum())
    if not math.isclose(total_mass, 1., rel_tol=1e-11, abs_tol=1e-12):
        raise ArithmeticError(f'Full-SI probability mass was not conserved: {total_mass}')
    if np.any(clean < -1e-14) or np.any(dirty < -1e-14):
        raise ArithmeticError('Full-SI probability mass became negative')
    return math.fsum(transfers)
