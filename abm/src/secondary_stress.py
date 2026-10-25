"""Seeded paired secondary-transmission stress test of a restricted assay.

Full SI permits acquired beds to contaminate a clean worker. The paired restricted
arm permits only the fixed initial source to do so. Both count first acquisition
of each nonsource bed, retain contamination on transfer, and use the same source
and four independent entry uniforms. Every record is a separate hypothetical
assay; no joint roster, clinical infection, turnover, clearance, or environment.
"""
from collections.abc import Iterable
import math
from numbers import Integral, Real
import time

import numpy as np


def paired_record_counts(sequence: Iterable[str], h: float, alpha: float,
                         beta: float, *, seed: int, replicates: int,
                         source_pool_size: int = 15, chunk_size: int = 8192,
                         deadline: float | None = None):
    """Return (full, restricted) integer counts for coupled replicate assays.

    Sources are uniform over the full pool, including unvisited beds. Labels map
    to arbitrary visited indices; unused pool indices represent absent sources.
    Newly acquired patients can transmit from subsequent visits in full SI only.
    A monotone coupling guarantees full counts are at least restricted counts.
    `deadline` optionally sets an absolute time.perf_counter budget boundary.
    """
    if isinstance(sequence, (str, bytes)):
        raise ValueError('sequence must be an iterable of bed labels')
    try:
        history = tuple(sequence)
    except TypeError as exc:
        raise ValueError('sequence must be an iterable of bed labels') from exc
    if not history or any(not isinstance(x, str) or not x.strip() for x in history):
        raise ValueError('history requires nonempty string bed labels')
    for name, value in [('h', h), ('alpha', alpha), ('beta', beta)]:
        if not isinstance(value, Real) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f'{name} must be a finite probability in [0,1]')
    for name, value in [('seed', seed), ('replicates', replicates),
                        ('source_pool_size', source_pool_size), ('chunk_size', chunk_size)]:
        if not isinstance(value, Integral) or isinstance(value, bool) or value < (0 if name == 'seed' else 1):
            raise ValueError(f'{name} must be an integer in its supported range')
    labels = tuple(dict.fromkeys(history))
    if source_pool_size < len(labels):
        raise ValueError('source pool must include every visited bed')
    positions = {label: i for i, label in enumerate(labels)}
    path = [positions[label] for label in history]
    rng = np.random.default_rng(int(seed))
    full_counts = np.empty(replicates, dtype=np.int32)
    restricted_counts = np.empty(replicates, dtype=np.int32)
    for start in range(0, replicates, chunk_size):
        if deadline is not None and time.perf_counter() >= deadline:
            raise TimeoutError('Predeclared stress-test wall budget exhausted')
        n = min(chunk_size, replicates-start)
        source = rng.integers(source_pool_size, size=n)
        present = source < len(labels)
        full = np.zeros((n, len(labels)), dtype=bool)
        rows = np.flatnonzero(present)
        full[rows, source[present]] = True
        restricted = full.copy()
        full_hand = np.zeros(n, dtype=bool)
        restricted_hand = np.zeros(n, dtype=bool)
        for bed in path:
            before_survives = rng.random(n) >= h
            full_hand &= before_survives
            restricted_hand &= before_survives
            acquisition_draw = rng.random(n) < alpha
            transfer_draw = rng.random(n) < beta
            # Compute both directions from the pre-contact state. Only an
            # initially dirty worker transfers; transfer never clears it.
            full_acquisition = ~full_hand & full[:, bed] & acquisition_draw
            full_transfer = full_hand & ~full[:, bed] & transfer_draw
            restricted_acquisition = ~restricted_hand & (source == bed) & acquisition_draw
            restricted_transfer = restricted_hand & ~restricted[:, bed] & transfer_draw
            full[:, bed] |= full_transfer
            restricted[:, bed] |= restricted_transfer
            full_hand |= full_acquisition
            restricted_hand |= restricted_acquisition
            after_survives = rng.random(n) >= h
            full_hand &= after_survives
            restricted_hand &= after_survives
        a = full.sum(axis=1, dtype=np.int32)-present
        b = restricted.sum(axis=1, dtype=np.int32)-present
        if np.any(a < b):
            raise AssertionError('Monotone full/restricted coupling was violated')
        full_counts[start:start+n] = a
        restricted_counts[start:start+n] = b
    return full_counts, restricted_counts
