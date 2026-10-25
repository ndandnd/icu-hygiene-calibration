"""Gate 2: alternative-workflow semantics of the adapter (plus original-mode identity).

Run from the ABM root (fast tests, ~1-4 min depending on load):
    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -B -m unittest -v \
        research/workflow_coefficient_2026-10-21/code/test_workflow_background.py
Full-horizon no-swap test (one seed, 90+56 days, several minutes):
    WFC_FULL_HORIZON=1 ... -k full_horizon
Evidence JSON from the full-horizon test goes to gate_checks/gate2/.
BenchmarkRecordTests re-reads every alternative record in records/ and
gate_checks/benchmark/records/ (skipped when there are none).
"""
import sys

sys.dont_write_bytecode = True

import json  # noqa: E402
import os  # noqa: E402
from pathlib import Path  # noqa: E402
import time  # noqa: E402
import unittest  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workflow_background as wb  # noqa: E402

STUDY = wb.STUDY
GATE2 = STUDY / 'gate_checks/gate2'
PRIMARY = [dict(id='iid_perfect', rho=0., e=1.), dict(id='persistent_perfect', rho=.5, e=1.),
           dict(id='iid_imperfect', rho=0., e=.8)]
_CLASSES = {}


def model_class(alternative):
    if alternative not in _CLASSES:
        _CLASSES[alternative] = wb.make_model_class(alternative=alternative)
    return _CLASSES[alternative]


def fixture_model(alternative, seed=3080160514, **changes):
    """Tick-0 fixture: no carriers, no admissions; all 30 beds occupied at init."""
    kwargs = dict(wb.bridge_background.FIXED_MODEL, seed=seed, prob_patient_sick=0., prob_new_patient=0.,
                  hygiene_opportunities='zone_separate', hygiene_autocorrelation=0., hygiene_efficacy=1.)
    kwargs.update(changes)
    return model_class(alternative)(warmup_ticks=10, h0=.9, h1=.95, histogram_cap=2,
                                    profiles=PRIMARY, **kwargs)


def patient_at(model, bed):
    return next(p for p in model.current_patients if p.pos == bed.pos)


def positions(model):
    return {p.unique_id: p.pos for p in model.current_patients}


def discharge(model, patient):
    model.discharged.append(patient)
    model._remove_discharged()
    model.refresh_beds()


def dirty(model, bed_position):
    site = model.response_site_by_bed_position[bed_position]
    if not site.colonized:
        model.response_environment_count += 1
    site.colonized = True
    return site


class FrozenCompositionTests(unittest.TestCase):
    def test_frozen_hashes_paths_and_no_bytecode(self):
        self.assertEqual(wb.verify_frozen(), wb.FROZEN_SHA256)
        self.assertTrue(sys.dont_write_bytecode)
        self.assertEqual(Path(wb.workflow_variant.__file__).resolve(),
                         (wb.ABM / 'research/robustness_uq_2026-10-17/workflow_variant.py').resolve())
        self.assertEqual(Path(wb.bridge_background.__file__).resolve(), (wb.BRIDGE / 'background.py').resolve())
        self.assertIs(wb.bridge_background.BridgeCoefficient, wb.bridge_coefficient.BridgeCoefficient)

    def test_surface_pairs_equal_archived_extended_run(self):
        archived = wb.verify_surface_pairs()
        self.assertEqual(len(archived), 22)
        legacy = {p['legacy_id']: p['id'] for p in archived if p.get('legacy_id')}
        self.assertEqual(legacy, wb.PRIMARY_PAIR_IDS)
        _, config = wb.archived_scientific_config()
        self.assertEqual(config['profiles'], PRIMARY)
        self.assertEqual((config['h0'], config['h1'], config['histogram_cap']), (.9, .95, 2))
        self.assertEqual((config['warmup_ticks'], config['observation_ticks']), (90*2592, 56*2592))
        self.assertEqual(config['model_overrides'], {'prob_patient_sick': .03})

    def test_alternative_rules_are_frozen_code_objects(self):
        Alt = model_class(True)
        for name in wb.VARIANT_METHODS:
            frozen = wb.FrozenVariant.__dict__[name]
            used = Alt.__dict__[name]
            self.assertIs(used.__code__, frozen.__code__, name)
            self.assertIs(used.__globals__, frozen.__globals__, name)
            if frozen.__code__.co_freevars:
                self.assertEqual(frozen.__code__.co_freevars, ('__class__',))
                self.assertIs(used.__closure__[0].cell_contents, Alt)
                self.assertIs(frozen.__closure__[0].cell_contents, wb.FrozenVariant)

    def test_original_class_keeps_core_rules(self):
        Orig = model_class(False)
        for name in ('_relocate_detected_patients', '_remove_discharged'):
            function = getattr(Orig, name)
            self.assertEqual(function.__qualname__, 'CPE_Model.' + name)
            self.assertEqual(function.__module__, '_bridge_core')
        model = fixture_model(False)
        self.assertFalse(hasattr(model, 'response_counts'))
        self.assertEqual(set(model.workflow_counters()), {'audit.' + k for k in wb.AUDIT_KEYS})


class OriginalEquivalenceTests(unittest.TestCase):
    def test_short_horizon_identical_to_frozen_background(self):
        _, config = wb.archived_scientific_config()
        config = dict(config, warmup_ticks=2592, observation_ticks=2592)
        seed = 3308159835
        frozen = wb.simulate_seed(seed, config, wb.bridge_background.make_background_model())
        duplicate = dict(id='dup_persistent', h0=.9, h1=.95, rho=.5, e=1.)
        adapted = wb.simulate_seed(seed, config, wb.make_factory(alternative=False,
                                                                 surface_pairs=wb.SURFACE_PAIRS + (duplicate,)))
        for key in ('ordered_bridge_events_sha256', 'event_counts', 'profiles', 'executed_ticks',
                    'observation_patient_days', 'observation_admissions', 'total_admissions', 'model_kwargs'):
            self.assertEqual(adapted[key], frozen[key], key)
        mix = adapted['workflow']['patient_contacts_by_period_status_bed']
        self.assertEqual(sum(mix.values()), frozen['event_counts']['patient'])
        self.assertEqual(sum(v for k, v in mix.items() if '|carrier|' in k),
                         frozen['profiles']['iid_perfect']['baseline']['event_counts'].get('imported_carrier_contacts', 0))
        # An extra arm fed the same event stream equals the frozen profile arm (cap aside).
        for arm in ('baseline', 'policy'):
            a, b = adapted['pairs']['dup_persistent'][arm], frozen['profiles']['persistent_perfect'][arm]
            for key in ('K', 'K_from_warmup_pickups', 'K_by_worker', 'raw_pairs_without_reset',
                        'event_counts', 'completion_law_endpoint'):
                self.assertEqual(a[key], b[key], key)
        for legacy, pair_id in wb.PRIMARY_PAIR_IDS.items():
            self.assertEqual({k: adapted['pairs'][pair_id][k] for k in ('baseline', 'policy')},
                             frozen['profiles'][legacy])


class AlternativeRuleTests(unittest.TestCase):
    def test_all_beds_full_at_init(self):
        model = fixture_model(True)
        self.assertEqual(len(model.current_patients), 30)
        self.assertEqual(model.empty_beds, set())
        self.assertEqual(model.response_environment_count, 0)

    def test_deferral_when_no_isolation_bed_is_empty_and_original_swaps(self):
        model = fixture_model(True)
        patient = patient_at(model, model.shared_beds[0])
        patient.positive = True
        before = positions(model)
        model.step()
        self.assertEqual(positions(model), before)
        self.assertFalse(patient.isolated)
        events = model.response_counts['workflow_events']
        self.assertEqual(events['isolation_deferred_no_empty_bed'], 1)
        self.assertEqual(events['isolation_relocation'], 0)
        self.assertEqual(model.response_counts['environment_events']['relocation_clear'], 0)
        self.assertEqual(model.relocation_audit['moves'], 0)
        model.step()  # deferral is counted once per eligible patient per tick
        self.assertEqual(events['isolation_deferred_no_empty_bed'], 2)
        model.check_workflow_invariants()
        model.assert_invariants()

        # Positive control: the frozen original rule swaps with an isolation occupant.
        original = fixture_model(False)
        sick = patient_at(original, original.shared_beds[0])
        target = original.isol_beds[0]
        displaced = patient_at(original, target)
        sick.positive = True
        original.step()
        self.assertEqual(sick.pos, target.pos)
        self.assertEqual(displaced.pos, original.shared_beds[0].pos)
        self.assertEqual(original.relocation_audit['swaps'], 1)
        self.assertEqual(original.relocation_audit['moves_into_occupied_bed'], 1)
        self.assertEqual(original.relocation_audit['moves_from_isolation_bed'], 1)

    def test_move_into_bed_emptied_by_discharge_same_tick(self):
        model = fixture_model(True)
        iso = model.isol_beds[3]
        leaving = patient_at(model, iso)
        leaving.stay = 1  # discharged in the next step, before relocation runs
        shared = model.shared_beds[5]
        patient = patient_at(model, shared)
        patient.positive = True
        others = {uid: pos for uid, pos in positions(model).items() if uid not in (patient.unique_id, leaving.unique_id)}
        model.step()
        environment = model.response_counts['environment_events']
        events = model.response_counts['workflow_events']
        self.assertEqual(model.discharged_patients, 1)
        self.assertEqual(environment['discharge_clear'], 1)
        self.assertEqual(events['isolation_relocation'], 1)
        self.assertEqual(environment['relocation_clear'], 1)
        self.assertEqual(events['isolation_deferred_no_empty_bed'], 0)
        self.assertEqual(patient.pos, iso.pos)
        self.assertTrue(patient.isolated)
        self.assertIn(shared.pos, model.empty_beds)
        self.assertEqual({uid: pos for uid, pos in positions(model).items() if uid != patient.unique_id}, others)
        self.assertEqual(model.relocation_audit['swaps'], 0)
        self.assertEqual(model.relocation_audit['moves_into_occupied_bed'], 0)
        self.assertEqual(model.relocation_audit['moves_shared_to_isolation'], 1)
        model.check_workflow_invariants()
        model.assert_invariants()

    def test_first_empty_isolation_bed_in_frozen_order_then_deferral(self):
        model = fixture_model(True)
        rng_before = model.rng.bit_generator.state
        for index in (5, 2):
            discharge(model, patient_at(model, model.isol_beds[index]))
        movers = [patient_at(model, model.shared_beds[i]) for i in (0, 4, 7)]
        for patient in movers:
            patient.positive = True
        model._relocate_detected_patients()
        model.refresh_beds()
        # Workflow rules (discharge clear, relocation, deferral) draw no random numbers.
        self.assertEqual(model.rng.bit_generator.state, rng_before)
        self.assertEqual(movers[0].pos, model.isol_beds[2].pos)
        self.assertEqual(movers[1].pos, model.isol_beds[5].pos)
        self.assertEqual(movers[2].pos, model.shared_beds[7].pos)
        events = model.response_counts['workflow_events']
        environment = model.response_counts['environment_events']
        self.assertEqual(events['isolation_relocation'], 2)
        self.assertEqual(events['isolation_deferred_no_empty_bed'], 1)
        self.assertEqual(environment['relocation_clear'], 2)
        self.assertEqual(environment['discharge_clear'], 2)
        self.assertEqual(model.relocation_audit['swaps'], 0)
        model.check_workflow_invariants()
        model.assert_invariants()

    def test_vacated_sites_cleared_destination_not_cleaned_at_move_in(self):
        model = fixture_model(True)
        iso, shared = model.isol_beds[1], model.shared_beds[2]
        dirty(model, iso.pos)
        discharge(model, patient_at(model, iso))
        environment = model.response_counts['environment_events']
        self.assertFalse(model.response_site_by_bed_position[iso.pos].colonized)
        self.assertEqual((environment['discharge_clear'], environment['dirty_discharge_clear']), (1, 1))
        destination_site = dirty(model, iso.pos)  # contaminated while vacant
        origin_site = dirty(model, shared.pos)
        patient = patient_at(model, shared)
        patient.positive = True
        model._relocate_detected_patients()
        self.assertEqual(patient.pos, iso.pos)
        self.assertFalse(origin_site.colonized)
        self.assertTrue(destination_site.colonized)
        self.assertEqual((environment['relocation_clear'], environment['dirty_relocation_clear']), (1, 1))
        self.assertEqual(model.response_environment_count, 1)
        # Clean-reservoir check rejects the artificial dirt, as it should at p=0.
        with self.assertRaises(RuntimeError):
            model.check_workflow_invariants()

    def test_occupied_or_non_bed_destination_rejected(self):
        model = fixture_model(True)
        patient = patient_at(model, model.shared_beds[0])
        origin = patient.pos
        with self.assertRaises(ValueError):
            model._move_patient(patient, model.isol_beds[0].pos)
        self.assertEqual(model.response_counts['workflow_events']['occupied_destination_rejected'], 1)
        self.assertEqual(patient.pos, origin)
        with self.assertRaises(ValueError):
            model._move_patient(patient, (0, 5))
        self.assertEqual(model.relocation_audit['moves'], 0)
        with self.assertRaises(RuntimeError):
            model._clear_vacated_site(origin, 'discharge')

    def test_invariant_detects_counter_inconsistency(self):
        model = fixture_model(True)
        model.check_workflow_invariants()
        model.response_counts['environment_events']['discharge_clear'] += 1
        with self.assertRaises(RuntimeError):
            model.check_workflow_invariants()


class AlternativeStressHorizonTests(unittest.TestCase):
    def test_saturated_isolation_short_horizon(self):
        """High carriage saturates isolation: relocations and deferrals both occur, never a swap."""
        days, warmup_days = 10, 4
        results = {}
        for alternative in (True, False):
            kwargs = dict(wb.bridge_background.FIXED_MODEL, seed=1820087184, prob_patient_sick=.3,
                          hygiene_opportunities='zone_separate', hygiene_autocorrelation=0., hygiene_efficacy=1.)
            model = model_class(alternative)(warmup_ticks=warmup_days*2592, h0=.9, h1=.95,
                                             histogram_cap=2, profiles=PRIMARY, **kwargs)
            isolation = model.bed_isolation_positions
            for day in range(days):
                for _ in range(2592):
                    deferred_before = model.response_counts['workflow_events']['isolation_deferred_no_empty_bed'] if alternative else 0
                    model.step()
                    if not alternative:
                        continue
                    # Independent per-tick check of the rule: every detected patient
                    # left in a shared bed was deferred this tick, and a deferral
                    # implies that no isolation bed was empty.
                    waiting = sum(1 for p in model.current_patients if p.positive and p.pos not in isolation)
                    deferred = model.response_counts['workflow_events']['isolation_deferred_no_empty_bed'] - deferred_before
                    self.assertEqual(deferred, waiting)
                    if deferred:
                        self.assertTrue(isolation <= {p.pos for p in model.current_patients})
                    for uid, origin, destination in model._audit_step_moves:
                        mover = next(p for p in model.current_patients if p.unique_id == uid)
                        self.assertTrue(mover.positive)
                        self.assertNotIn(origin, isolation)
                        self.assertIn(destination, isolation)
                model.assert_invariants()
                model.check_workflow_invariants()
                for patient in model.current_patients:
                    self.assertEqual(patient.isolated, patient.pos in model.bed_isolation_positions)
            results[alternative] = model.bridge_result()
        alt = results[True]['workflow']['counters_total']
        self.assertGreater(alt['workflow_events.isolation_relocation'], 0)
        self.assertGreater(alt['workflow_events.isolation_deferred_no_empty_bed'], 0)
        self.assertGreater(alt['environment_events.discharge_clear'], 0)
        self.assertEqual(alt['environment_events.relocation_clear'], alt['workflow_events.isolation_relocation'])
        self.assertEqual(alt['audit.swaps'], 0)
        self.assertEqual(alt['audit.moves_into_occupied_bed'], 0)
        self.assertEqual(alt['audit.moves_from_isolation_bed'], 0)
        self.assertGreater(results[False]['workflow']['counters_total']['audit.swaps'], 0)


@unittest.skipUnless(os.environ.get('WFC_FULL_HORIZON') == '1', 'set WFC_FULL_HORIZON=1 (several minutes)')
class FullHorizonTests(unittest.TestCase):
    SEED = 3308159835  # also a Gate 1 seed: its original-mode audit is in gate_checks/gate1

    def test_full_horizon_no_swap(self):
        _, config = wb.archived_scientific_config()
        config = dict(config, seeds=[self.SEED])
        factory = wb.make_factory(alternative=True)
        captured = []

        def capture(**kwargs):
            captured.append(factory(**kwargs))
            return captured[-1]
        started = time.monotonic()
        cpu = time.process_time()
        row = wb.simulate_seed(self.SEED, config, capture)  # bridge_result enforces invariants
        model = captured[0]
        counters = row['workflow']['counters_total']
        checks = dict(
            swaps_zero=counters['audit.swaps'] == 0,
            no_move_into_occupied=counters['audit.moves_into_occupied_bed'] == 0,
            no_move_out_of_isolation=counters['audit.moves_from_isolation_bed'] == 0,
            all_moves_shared_to_isolation=counters['audit.moves_shared_to_isolation'] == counters['audit.moves'],
            moves_equal_relocations=counters['audit.moves'] == counters['workflow_events.isolation_relocation'],
            relocation_clears_equal_relocations=counters['environment_events.relocation_clear'] == counters['workflow_events.isolation_relocation'],
            discharge_clears_equal_discharges=counters['environment_events.discharge_clear'] == model.discharged_patients,
            no_dirty_clears=counters['environment_events.dirty_discharge_clear'] == counters['environment_events.dirty_relocation_clear'] == 0,
            no_occupied_destination_attempt=counters['workflow_events.occupied_destination_rejected'] == 0,
            relocations_occurred=counters['workflow_events.isolation_relocation'] > 0,
            executed_full_horizon=row['executed_ticks'] == 146*2592)
        gate1 = STUDY / f'gate_checks/gate1/{self.SEED}.json'
        original_audit = json.loads(gate1.read_text())['workflow'] if gate1.exists() else None
        evidence = dict(test='full_horizon_no_swap', seed=self.SEED, checks=checks, counters_total=counters,
                        counters_observation=row['workflow']['counters_observation'],
                        discharged_patients=model.discharged_patients,
                        deferral_patient_days_total=counters['workflow_events.isolation_deferred_no_empty_bed']/2592,
                        primary_K={k: (v['baseline']['K'], v['policy']['K']) for k, v in row['profiles'].items()},
                        ordered_bridge_events_sha256=row['ordered_bridge_events_sha256'],
                        original_mode_audit_same_seed_from_gate1=original_audit,
                        wall_seconds=time.monotonic() - started, cpu_seconds=time.process_time() - cpu,
                        adapter_sha256=wb.adapter_sha256())
        GATE2.mkdir(parents=True, exist_ok=True)
        (GATE2 / f'full_horizon_{self.SEED}.json').write_text(json.dumps(evidence, indent=2, sort_keys=True) + '\n')
        for name, ok in checks.items():
            self.assertTrue(ok, name)


class BenchmarkRecordTests(unittest.TestCase):
    def test_alternative_records_never_swap(self):
        paths = sorted((STUDY / 'records').glob('[0-9]*.json')) + \
            sorted((STUDY / 'gate_checks/benchmark/records').glob('[0-9]*.json'))
        if not paths:
            self.skipTest('no alternative-workflow records yet')
        for path in paths:
            row = json.loads(path.read_text())
            counters = row['workflow']['counters_total']
            self.assertTrue(row['workflow']['alternative'])
            self.assertEqual(counters['audit.swaps'], 0, path.name)
            self.assertEqual(counters['audit.moves_into_occupied_bed'], 0, path.name)
            self.assertEqual(counters['audit.moves_from_isolation_bed'], 0, path.name)
            self.assertEqual(counters['workflow_events.occupied_destination_rejected'], 0, path.name)
            self.assertEqual(counters['environment_events.relocation_clear'],
                             counters['workflow_events.isolation_relocation'], path.name)
            self.assertEqual(row['executed_ticks'], 146*2592)


if __name__ == '__main__':
    unittest.main()
