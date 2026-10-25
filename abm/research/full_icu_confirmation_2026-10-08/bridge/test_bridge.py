"""Exact tiny survival enumeration and bounded frozen-engine replay checks."""
import importlib.util
import json
import math
import unittest

from coefficient import BridgeCoefficient, stationary_q2


def recurrence(events, h0=.4, h1=.8, rho=.5, e=.7, cap=2):
    result = BridgeCoefficient(h0, h1, rho, e, cap)
    for kind, uid, observed, details in events:
        if kind == 'opportunity':
            result.opportunity(uid, observed)
        elif kind == 'reset':
            result.reset(uid)
        else:
            result.patient_contact(uid, *details, observed)
    return result.result()


def enumerate_histories(events, h0=.4, h1=.8, rho=.5, e=.7):
    """Explicit completion/efficacy paths; no matrices or aggregated state mass."""
    branches = [(1., {}, {}, 0., {})]
    for kind, uid, observed, details in events:
        updated = []
        for chance, previous, alive, reward, histogram in branches:
            if kind == 'opportunity':
                h = h1 if observed else h0
                old = previous.get(uid)
                probability = h if old is None else (h + rho*(1-h) if old else h*(1-rho))
                for complete, pcomplete in ((False, 1-probability), (True, probability)):
                    for clear, pefficacy in ((False, 1-e), (True, e)) if complete else ((False, 1.),):
                        prob = chance * pcomplete * pefficacy
                        if not prob:
                            continue
                        memory = dict(previous, **{str(uid): complete})
                        # uid are strings in this independent enumerator.
                        sources = dict(alive)
                        sources[uid] = [] if clear else [(weight, n+1) for weight, n in alive.get(uid, [])]
                        updated.append((prob, memory, sources, reward, dict(histogram)))
            elif kind == 'reset':
                sources = dict(alive)
                sources[uid] = []
                updated.append((chance, previous, sources, reward, histogram))
            else:
                patient_id, carrier, weight = details
                sources = dict(alive)
                hist = dict(histogram)
                if carrier:
                    sources[uid] = alive.get(uid, []) + [(weight, 0)]
                elif observed:
                    for pickup, n in alive.get(uid, []):
                        reward += weight * pickup
                        hist[n] = hist.get(n, 0.) + weight * pickup
                updated.append((chance, previous, sources, reward, hist))
        branches = updated
    assert math.isclose(sum(row[0] for row in branches), 1., abs_tol=1e-12)
    histogram = {}
    for chance, _, _, _, counts in branches:
        for n, value in counts.items():
            histogram[n] = histogram.get(n, 0.) + chance * value
    return sum(chance*reward for chance, _, _, reward, _ in branches), histogram


class ExactBridgeTests(unittest.TestCase):
    def setUp(self):
        self.events = [
            ('patient', 'A', False, ('import0', True, .6)),  # before first-ever opportunity
            ('opportunity', 'A', False, ()),
            ('patient', 'B', False, ('import1', True, 1.)),
            ('opportunity', 'A', False, ()),
            ('patient', 'A', False, ('import2', True, 1.)),
            ('opportunity', 'B', False, ()),
            ('patient', 'A', True, ('sus0', False, .6)),
            ('opportunity', 'A', True, ()),
            ('patient', 'A', True, ('sus1', False, 1.)),
            ('reset', 'B', True, ()),
            ('patient', 'B', True, ('sus2', False, 1.)),
            ('reset', 'A', True, ()),
            ('patient', 'A', True, ('import3', True, .6)),
            ('opportunity', 'A', True, ()),
            ('patient', 'A', True, ('sus3', False, 1.)),
        ]

    def test_matches_exact_history_enumeration_all_profiles(self):
        for rho, efficacy in ((0., 1.), (.5, 1.), (0., .8), (.75, .4)):
            with self.subTest(rho=rho, efficacy=efficacy):
                answer, hist = enumerate_histories(self.events, .9, .95, rho, efficacy)
                estimated = recurrence(self.events, .9, .95, rho, efficacy)
                self.assertAlmostEqual(answer, estimated['K'], places=12)
                for label, mass in estimated['K_by_opportunities'].items():
                    expected = sum(v for k, v in hist.items() if k >= 3) if label == '>=3' else hist.get(int(label), 0.)
                    self.assertAlmostEqual(mass, expected, places=12)

    def test_tail_conservation_and_histogram_disabled(self):
        events = [('patient', 'A', False, ('source', True, .6))]
        events += [('opportunity', 'A', True, ())] * 6
        events += [('patient', 'A', True, ('recipient', False, .6))]
        exact, _ = enumerate_histories(events)
        for cap in (0, 1, 2, 8, None):
            result = recurrence(events, cap=cap)
            self.assertAlmostEqual(result['K'], exact, places=12)
            if cap is not None:
                self.assertEqual(sum(result['raw_pairs_without_reset_by_opportunities'].values()), 1)
                self.assertAlmostEqual(sum(result['K_by_opportunities'].values()), exact, places=12)

    def test_reset_retains_memory_and_switch_retains_warmup_mass(self):
        state = BridgeCoefficient(.9, .1, .5, 1., 2)
        state.opportunity('A', False)
        state.patient_contact('A', 'imported', True, .6, False)
        state.patient_contact('A', 'susceptible', False, .6, True)
        self.assertAlmostEqual(state.k_warmup_pickups, .36)
        old_nu = list(state.worker('A').nu)
        state.reset('A')
        self.assertEqual(old_nu, state.worker('A').nu)
        state.opportunity('A', True)
        self.assertAlmostEqual(state.worker('A').nu[2], .5)
        self.assertEqual(state.worker('A').b, [0., 0., 0.])

    def test_nonstationary_q2_is_not_coefficient(self):
        events = [('opportunity', 'A', False, ()), ('patient', 'A', False, ('S', True, 1.)),
                  ('opportunity', 'A', True, ()), ('opportunity', 'A', True, ()),
                  ('patient', 'A', True, ('R', False, 1.))]
        result = recurrence(events, .9, .1, .5, 1.)['K']
        self.assertNotAlmostEqual(result, stationary_q2(.1, .5, 1.))

    def test_cross_worker_and_outside_observation(self):
        events = [('patient', 'A', False, ('S', True, 1.)),
                  ('patient', 'A', False, ('R', False, 1.)),
                  ('patient', 'B', True, ('R2', False, 1.))]
        self.assertEqual(recurrence(events)['K'], 0.)

    def test_exact_gate_enumeration_third_order_repeated_transfer(self):
        # S pickup, A transfer, B transfer. One dirty hand persists and A/B each
        # acquire with p²; exactly two gates contributes 2p²(1-p), three 2p³.
        for p in (.01, .1, .3):
            expected = 0.
            two_gate = 0.
            for mask in range(8):
                gates = [(mask >> j) & 1 for j in range(3)]
                chance = math.prod(p if gate else 1-p for gate in gates)
                acquired = gates[0] * (gates[1] + gates[2])
                expected += chance * acquired
                if sum(gates) == 2:
                    two_gate += chance * acquired
            self.assertAlmostEqual(expected, 2*p*p)
            self.assertAlmostEqual(two_gate, 2*p*p*(1-p))

    def test_clean_environment_requires_four_gates(self):
        # After deposition, original worker clears; a distinct worker picks up
        # and transfers. No patient-source bridge joins the two workers.
        for p in (.01, .1, .3):
            expected = 0.
            for mask in range(16):
                gates = [(mask >> j) & 1 for j in range(4)]
                expected += math.prod(p if gate else 1-p for gate in gates) * math.prod(gates)
            self.assertAlmostEqual(expected, p**4)


@unittest.skipUnless(importlib.util.find_spec('mesa'), 'Mesa runtime unavailable')
class FrozenBackgroundTests(unittest.TestCase):
    def test_role_reset_predicates_clear_mass_only(self):
        from background import FIXED_MODEL, make_background_model
        Background = make_background_model()
        kwargs = dict(FIXED_MODEL, seed=142, hygiene_opportunities='zone_separate')
        for uid, tick, finishing in ((-1, 864, False), (-12, 108, False),
                                     (-17, 1296, True), (-18, 1, True)):
            with self.subTest(worker=uid, tick=tick):
                model = Background(warmup_ticks=0, h0=.9, h1=.95, **kwargs)
                worker = next(a for a in model.schedule.agents if a.unique_id == uid)
                # Suppress movement/contact: this test isolates the frozen
                # step's role-boundary assignments from hygiene opportunities.
                worker.move = lambda: None
                worker.spread = lambda: None
                if hasattr(worker, 'leaveICU'):
                    worker.leaveICU = lambda: None
                if finishing:
                    worker.activated = True
                    worker._route_entry_completed = True
                    worker.moveTick = len(worker.path)-1
                model.tick = tick
                for arm in model.bridge_arms:
                    arm.opportunity(uid, False)
                    arm.patient_contact(uid, 'imported', True, .6, False)
                laws = [list(arm.worker(uid).nu) for arm in model.bridge_arms]
                worker.step()
                self.assertEqual(model.bridge_counts['reset'], 1)
                for arm, law in zip(model.bridge_arms, laws):
                    self.assertEqual(arm.worker(uid).b, [0., 0., 0.])
                    self.assertEqual(arm.worker(uid).nu, law)

    def test_same_seed_replays_unmodified_p0_engine(self):
        from background import FIXED_MODEL, load_core, make_background_model
        _, plain_core = load_core(tag='_bridge_plain_test')
        Background = make_background_model()
        for opportunities, rho, efficacy in (('separate', 0., 1.), ('shared', .5, 1.), ('zone_separate', 0., .8)):
            with self.subTest(opportunities=opportunities):
                kwargs = dict(FIXED_MODEL, seed=142, prob_patient_sick=.2,
                              prob_transmission=0., icu_hcw_wash_rate=.9, outside_hcw_wash_rate=.9,
                              hygiene_opportunities=opportunities,
                              hygiene_autocorrelation=rho, hygiene_efficacy=efficacy)
                plain = plain_core.CPE_Model(**kwargs)
                observed = Background(warmup_ticks=128, h0=.9, h1=.95, histogram_cap=2, **kwargs)
                for _ in range(2700):
                    plain.step()
                    observed.step()
                attrs = ('unique_id', 'pos', 'colonized', 'stay', 'positive', 'isolated',
                         '_hygiene_previous_completion', '_hygiene_zone', 'moveTick', 'activated',
                         'patient_contacts', 'environment_contacts', 'hygiene_opportunity_count')
                def state(model):
                    return [[getattr(agent, attr, None) for attr in attrs] for agent in model.schedule.agents]
                self.assertEqual(state(plain), state(observed))
                self.assertEqual(json.dumps(plain.rng.bit_generator.state), json.dumps(observed.rng.bit_generator.state))
                self.assertEqual(plain.random.getstate(), observed.random.getstate())
                self.assertEqual(plain.patient_days, observed.patient_days)
                self.assertEqual(plain.susceptible_patient_days, observed.susceptible_patient_days)
                result = observed.bridge_result()
                self.assertGreater(result['event_counts']['reset'], 0)
                self.assertGreater(result['profiles']['configured']['baseline']['K'], 0)


if __name__ == '__main__':
    unittest.main()
