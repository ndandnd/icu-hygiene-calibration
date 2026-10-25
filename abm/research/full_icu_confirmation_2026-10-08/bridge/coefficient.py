"""Exact conditional p² coefficient; pure Python, no simulation dependencies.

Rows are absent, previous miss, previous complete. ν is the unselected
completion-state law; b contains summed pickup weights times surviving state
mass. Injection b += cν is linear, so no pair enumeration is required.
"""
from collections import Counter
import math


def probability(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(name + " must be finite in [0,1]")
    return float(value)


def advance(vector, h, rho, efficacy=0.0):
    """Row-vector QD; efficacy=0 gives the ordinary state transition Q."""
    absent, miss, complete = vector
    p_miss = h * (1 - rho)
    p_complete = h + rho * (1 - h)
    yes = absent * h + miss * p_miss + complete * p_complete
    no = absent * (1 - h) + miss * (1 - p_miss) + complete * (1 - p_complete)
    return [0.0, no, yes * (1 - efficacy)]


class WorkerMass:
    def __init__(self, cap):
        self.nu = [1.0, 0.0, 0.0]
        self.b = [0.0, 0.0, 0.0]
        self.bins = None if cap is None else [[0.0, 0.0, 0.0] for _ in range(cap + 2)]
        self.raw = None if cap is None else [0 for _ in range(cap + 2)]
        self.pickups_since_reset = 0


class BridgeCoefficient:
    """One hygiene schedule evaluated on a common epidemiological background.

    event order defines whether a wash lies in a bridge. ``observed`` is true
    only after warmup. Patient contacts must refer to individual admissions.
    ``carrier`` means imported in the clean p=0 background. No site injections
    exist under clean initial sites. Fixed cap gives O(events) time and memory.
    """
    def __init__(self, h0, h1, rho=0.0, efficacy=1.0, histogram_cap=8):
        self.h0 = probability(h0, 'h0')
        self.h1 = probability(h1, 'h1')
        self.rho = probability(rho, 'rho')
        if self.rho == 1:
            raise ValueError('rho must be below 1')
        self.efficacy = probability(efficacy, 'efficacy')
        if histogram_cap is not None and (type(histogram_cap) is not int or histogram_cap < 0):
            raise ValueError('histogram_cap must be a nonnegative integer or None')
        self.cap = histogram_cap
        self.workers = {}
        self.k = 0.0
        self.k_warmup_pickups = 0.0
        self.warmup_mass = {}
        self.k_by_worker = Counter()
        self.k_by_opportunities = None if self.cap is None else [0.0] * (self.cap + 2)
        self.raw_pairs_by_opportunities = None if self.cap is None else [0] * (self.cap + 2)
        self.raw_pairs_without_reset = 0
        self.counts = Counter()

    def worker(self, uid):
        if uid not in self.workers:
            self.workers[uid] = WorkerMass(self.cap)
            self.warmup_mass[uid] = [0.0, 0.0, 0.0]
        return self.workers[uid]

    def opportunity(self, uid, observed):
        state = self.worker(uid)
        h = self.h1 if observed else self.h0
        state.nu = advance(state.nu, h, self.rho)
        state.b = advance(state.b, h, self.rho, self.efficacy)
        self.warmup_mass[uid] = advance(self.warmup_mass[uid], h, self.rho, self.efficacy)
        if state.bins is not None:
            updated = [[0.0, 0.0, 0.0] for _ in state.bins]
            raw = [0 for _ in state.raw]
            for i, mass in enumerate(state.bins):
                j = min(i + 1, self.cap + 1)
                transformed = advance(mass, h, self.rho, self.efficacy)
                updated[j] = [a + b for a, b in zip(updated[j], transformed)]
                raw[j] += state.raw[i]
            state.bins, state.raw = updated, raw
        self.counts['opportunities_observed' if observed else 'opportunities_warmup'] += 1

    def reset(self, uid):
        state = self.worker(uid)
        state.b = [0.0, 0.0, 0.0]
        self.warmup_mass[uid] = [0.0, 0.0, 0.0]
        state.pickups_since_reset = 0
        if state.bins is not None:
            state.bins = [[0.0, 0.0, 0.0] for _ in state.bins]
            state.raw = [0 for _ in state.raw]
        self.counts['implicit_resets'] += 1

    def patient_contact(self, uid, patient_id, carrier, weight, observed):
        if type(carrier) is not bool:
            raise ValueError('carrier must be boolean')
        weight = probability(weight, 'isolation weight')
        state = self.worker(uid)
        self.counts['imported_carrier_contacts' if carrier else 'susceptible_contacts'] += 1
        if carrier:
            injected = [weight * v for v in state.nu]
            state.b = [a + b for a, b in zip(state.b, injected)]
            if not observed:
                self.warmup_mass[uid] = [a + b for a, b in zip(self.warmup_mass[uid], injected)]
            state.pickups_since_reset += 1
            if state.bins is not None:
                state.bins[0] = [a + b for a, b in zip(state.bins[0], injected)]
                state.raw[0] += 1
        elif observed:
            contribution = weight * sum(state.b)
            self.k += contribution
            self.k_by_worker[str(uid)] += contribution
            self.k_warmup_pickups += weight * sum(self.warmup_mass[uid])
            self.raw_pairs_without_reset += state.pickups_since_reset
            if state.bins is not None:
                for i, mass in enumerate(state.bins):
                    self.k_by_opportunities[i] += weight * sum(mass)
                    self.raw_pairs_by_opportunities[i] += state.raw[i]

    def result(self):
        labels = None if self.cap is None else [str(i) for i in range(self.cap + 1)] + [f'>={self.cap + 1}']
        histogram = None if labels is None else dict(zip(labels, self.k_by_opportunities))
        if histogram is not None and not math.isclose(sum(histogram.values()), self.k, rel_tol=2e-11, abs_tol=1e-10):
            raise RuntimeError('Bridge histogram fails coefficient reconciliation')
        return {'K': self.k, 'K_from_warmup_pickups': self.k_warmup_pickups,
                'K_by_worker': dict(self.k_by_worker), 'K_by_opportunities': histogram,
                'raw_pairs_without_reset': self.raw_pairs_without_reset,
                'raw_pairs_without_reset_by_opportunities': None if labels is None else dict(zip(labels, self.raw_pairs_by_opportunities)),
                'event_counts': dict(self.counts),
                'completion_law_endpoint': {str(uid): state.nu for uid, state in self.workers.items()}}


def stationary_q2(h, rho=0.0, efficacy=1.0):
    law = [0.0, 1 - h, h]
    for _ in range(2):
        law = advance(law, h, rho, efficacy)
    return sum(law)
