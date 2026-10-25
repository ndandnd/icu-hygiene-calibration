"""Bounded exact full-SI multiset null using exchangeable occupancy classes.

Each active bed is represented by (remaining visits, colonized status). Beds
with zero remaining visits are omitted: first acquisitions already earned their
reward, and these beds cannot affect future contacts. The worker state is kept
exactly as a clean/dirty pair. This quotient requires identical alpha/beta across
beds, permanent colonization and independent perfect before/after hygiene.
It is conditional finite-state mathematics, not a clinical inference method.
"""
from collections import Counter
import math
from numbers import Integral, Real
import time


class StateLimitExceeded(RuntimeError):
    """The predeclared maximum number of distinct occupancy states was reached."""


def exact_full_si_permuted_expected_acquisitions(
        sequence, h, alpha, beta, source_pool_size=15, *, max_states=250_000,
        deadline=None, diagnostics=None):
    """Uniform multiset-order and uniform-source expected first acquisitions.

    `None` chooses the visited-source pool; otherwise absent initial sources
    contribute zero. `max_states` bounds all started memo states, including the
    active recursion stack. `deadline` is an absolute time.perf_counter deadline.
    Optional mutable `diagnostics` receives counts/timing even on budget failure.
    No persistent global cache or simulation randomness is used. This recursive
    implementation supports at most 200 entries; longer histories fail explicitly.
    """
    if isinstance(sequence, (str, bytes)):
        raise ValueError('sequence must be an iterable of bed labels, not a string')
    try:
        history = tuple(sequence)
    except TypeError as exc:
        raise ValueError('sequence must be an iterable of bed labels') from exc
    if not history or any(not isinstance(x, str) or not x.strip() for x in history):
        raise ValueError('history requires nonempty string bed labels')
    if len(history) > 200:
        raise ValueError('bounded recursive evaluator supports at most 200 entries')
    for name, value in [('h', h), ('alpha', alpha), ('beta', beta)]:
        if not isinstance(value, Real) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f'{name} must be a finite probability in [0,1]')
    frequencies = Counter(history)
    visited = len(frequencies)
    pool = visited if source_pool_size is None else source_pool_size
    if not isinstance(pool, Integral) or isinstance(pool, bool) or pool < visited:
        raise ValueError('source pool must be an integer at least the visited-bed count, or None')
    if not isinstance(max_states, Integral) or isinstance(max_states, bool) or max_states < 1:
        raise ValueError('max_states must be a positive integer')
    if deadline is not None and (not isinstance(deadline, Real) or not math.isfinite(deadline)):
        raise ValueError('deadline must be a finite absolute perf_counter time or None')
    if diagnostics is not None and not isinstance(diagnostics, dict):
        raise ValueError('diagnostics must be a dict or None')
    start = time.perf_counter()
    memo = {}
    started = 0
    status = 'running'
    max_r = max(frequencies.values())
    multiplicities = Counter(frequencies.values())
    wash_failure = float(1-h)
    q = wash_failure*wash_failure
    source_dirty_from_clean = float(alpha)*wash_failure
    source_dirty_from_dirty = wash_failure*(wash_failure+float(h)*float(alpha))
    immediate = float(beta)*wash_failure
    susceptible_dirty_no_transfer = q*(1-float(beta))
    susceptible_clean_no_transfer = float(h)+wash_failure*(1-float(beta))*float(h)

    def child(state, r, colonized, becomes_colonized):
        updated = list(state)
        updated[(max_r if colonized else 0)+r-1] -= 1
        if r > 1:
            updated[(max_r if becomes_colonized else 0)+r-2] += 1
        return tuple(updated)

    def future(state):
        nonlocal started
        if deadline is not None and time.perf_counter() >= deadline:
            raise TimeoutError('Exact occupancy wall budget exhausted')
        if state in memo:
            return memo[state]
        if started >= max_states:
            raise StateLimitExceeded('Exact occupancy state budget exhausted')
        started += 1
        remaining = sum((i % max_r+1)*count for i, count in enumerate(state))
        if not remaining:
            memo[state] = (0., 0.)
            return memo[state]
        clean_terms, dirty_terms = [], []
        for i, count in enumerate(state):
            if not count:
                continue
            r, colonized = i % max_r+1, i >= max_r
            weight = r*count/remaining
            if colonized:
                f0, f1 = future(child(state, r, True, True))
                clean = (1-source_dirty_from_clean)*f0+source_dirty_from_clean*f1
                dirty = (1-source_dirty_from_dirty)*f0+source_dirty_from_dirty*f1
            else:
                f0, f1 = future(child(state, r, False, False))
                g0, g1 = future(child(state, r, False, True))
                clean = f0
                # Immediate reward and future contributions use the original
                # before-wash/contact/after-wash kernel. Transfer retains dirty
                # hands; the newly positive bed becomes a source only later.
                dirty = (immediate+susceptible_clean_no_transfer*f0
                         +susceptible_dirty_no_transfer*f1
                         +immediate*float(h)*g0+immediate*wash_failure*g1)
            clean_terms.append(weight*clean)
            dirty_terms.append(weight*dirty)
        memo[state] = (math.fsum(clean_terms), math.fsum(dirty_terms))
        return memo[state]

    try:
        if deadline is not None and start >= deadline:
            raise TimeoutError('Exact occupancy wall budget exhausted')
        if visited == 1 or h == 1 or alpha == 0 or beta == 0:
            answer = 0.
        else:
            all_susceptible = [multiplicities.get(r, 0) for r in range(1, max_r+1)]+[0]*max_r
            terms = []
            for r, count in sorted(multiplicities.items()):
                initial = all_susceptible.copy()
                initial[r-1] -= 1
                initial[max_r+r-1] += 1
                terms.append(count/int(pool)*future(tuple(initial))[0])
            answer = math.fsum(terms)
        status = 'complete'
        return answer
    except StateLimitExceeded:
        status = 'state_limit'
        raise
    except TimeoutError:
        status = 'time_limit'
        raise
    finally:
        if diagnostics is not None:
            diagnostics.update(status=status, states_started=started, memoized_states=len(memo),
                               max_states=int(max_states), runtime_seconds=time.perf_counter()-start,
                               visited_beds=visited, entries=len(history), max_frequency=max_r,
                               source_pool_size=int(pool), absent_source_probability=(int(pool)-visited)/int(pool))
