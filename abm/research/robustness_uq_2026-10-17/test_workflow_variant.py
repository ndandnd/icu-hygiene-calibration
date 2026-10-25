"""Bounded implementation checks, not study outcomes or calibration runs."""

from copy import deepcopy
import hashlib
import importlib.util
from pathlib import Path
import sys
import unittest


HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


variant = load("_test_robustness_workflow", HERE / "workflow_variant.py")
runner = load("_test_robustness_pair_runner", ABM / "research/full_icu_response_2026-10-08/run.py")


class WorkflowTests(unittest.TestCase):
    def kwargs(self, **changes):
        kwargs = dict(
            runner.FIXED_MODEL, seed=77, prob_patient_sick=0,
            initial_environment_contaminated=False, prob_transmission=0.3,
            icu_hcw_wash_rate=0.9, outside_hcw_wash_rate=0.9,
            hygiene_opportunities="zone_separate",
        )
        kwargs.update(changes)
        return kwargs

    def model(self, **changes):
        return variant.WorkflowRobustnessModel(**self.kwargs(**changes))

    def patient_at(self, model, bed):
        return next(p for p in model.current_patients if p.pos == bed.pos)

    def dirty_site(self, model, bed):
        site = model.response_site_by_bed_position[bed.pos]
        if not site.colonized:
            model.response_environment_count += 1
        site.colonized = True
        site._response_environment_origin = {"kind": "fixture", "zone": bed.unique_id}
        return site

    def discharge(self, model, patient):
        model.discharged.append(patient)
        model._remove_discharged()
        model.refresh_beds()

    def mark_positive(self, patient):
        # Detection flags are sufficient for the relocation mechanism under test.
        patient.positive = True

    def test_discharge_cleans_site_immediately_and_keeps_cleaning_clock(self):
        model = self.model()
        bed = model.shared_beds[0]
        patient = self.patient_at(model, bed)
        site = self.dirty_site(model, bed)
        clock = site.clean_tick
        numpy_state = deepcopy(model.rng.bit_generator.state)
        python_state = model.random.getstate()
        self.discharge(model, patient)
        self.assertFalse(site.colonized)
        self.assertIsNone(site._response_environment_origin)
        self.assertEqual(site.clean_tick, clock)
        self.assertEqual(model.response_environment_count, 0)
        self.assertEqual(model.response_counts["environment_events"]["discharge_clear"], 1)
        self.assertEqual(model.response_counts["environment_events"]["dirty_discharge_clear"], 1)
        self.assertEqual(model.response_counts["environment_events"]["scheduled_clear"], 0)
        self.assertEqual(model.rng.bit_generator.state, numpy_state)
        self.assertEqual(model.random.getstate(), python_state)
        self.assertEqual(len(model.current_patients), 29)
        self.assertNotIn(patient, model.schedule.agents)
        self.assertIn(bed.pos, model.empty_beds)
        model.assert_response_invariants()

    def test_clean_discharge_count_does_not_decrement_environment_below_zero(self):
        model = self.model()
        self.discharge(model, self.patient_at(model, model.shared_beds[0]))
        self.assertEqual(model.response_environment_count, 0)
        self.assertEqual(model.response_counts["environment_events"]["discharge_clear"], 1)
        self.assertEqual(model.response_counts["environment_events"]["dirty_discharge_clear"], 0)
        model.assert_response_invariants()

    def test_relocation_uses_empty_bed_cleans_origin_and_preserves_occupants(self):
        model = self.model()
        source, destination = model.shared_beds[0], model.isol_beds[2]
        departing = self.patient_at(model, destination)
        self.discharge(model, departing)
        patient = self.patient_at(model, source)
        self.mark_positive(patient)
        origin_site = self.dirty_site(model, source)
        destination_site = self.dirty_site(model, destination)
        origin_clock, destination_clock = origin_site.clean_tick, destination_site.clean_tick
        positions = {p.unique_id: p.pos for p in model.current_patients}
        numpy_state = deepcopy(model.rng.bit_generator.state)
        model._relocate_detected_patients()
        model.refresh_beds()
        self.assertEqual(patient.pos, destination.pos)
        self.assertTrue(patient.isolated)
        self.assertEqual((patient.x, patient.y), destination.pos)
        for other in model.current_patients:
            if other is not patient:
                self.assertEqual(other.pos, positions[other.unique_id])
        self.assertFalse(origin_site.colonized)
        self.assertIsNone(origin_site._response_environment_origin)
        self.assertTrue(destination_site.colonized)
        self.assertEqual((origin_site.clean_tick, destination_site.clean_tick), (origin_clock, destination_clock))
        self.assertEqual(model.response_counts["environment_events"]["relocation_clear"], 1)
        self.assertEqual(model.response_counts["environment_events"]["dirty_relocation_clear"], 1)
        self.assertEqual(model.response_counts["workflow_events"]["isolation_relocation"], 1)
        self.assertEqual(model.rng.bit_generator.state, numpy_state)
        self.assertEqual(len(model.current_patients), 29)
        self.assertEqual(model.discharged_patients, 1)
        self.assertEqual(model.empty_beds, {source.pos})
        model.assert_response_invariants()

    def test_full_isolation_capacity_leaves_patient_in_place_without_swapping(self):
        model = self.model()
        patient = self.patient_at(model, model.shared_beds[0])
        self.mark_positive(patient)
        positions = {p.unique_id: p.pos for p in model.current_patients}
        model._relocate_detected_patients()
        model.refresh_beds()
        self.assertEqual({p.unique_id: p.pos for p in model.current_patients}, positions)
        self.assertEqual(model.response_counts["workflow_events"]["isolation_deferred_no_empty_bed"], 1)
        self.assertEqual(model.response_counts["workflow_events"]["isolation_relocation"], 0)
        self.assertEqual(model.response_counts["environment_events"]["relocation_clear"], 0)
        model.assert_response_invariants()

    def test_occupied_destination_is_rejected_before_any_move_or_clear(self):
        model = self.model()
        source, destination = model.shared_beds[0], model.isol_beds[0]
        patient = self.patient_at(model, source)
        site = self.dirty_site(model, source)
        positions = {p.unique_id: p.pos for p in model.current_patients}
        with self.assertRaisesRegex(ValueError, "empty destination"):
            model._move_patient(patient, destination.pos)
        self.assertEqual({p.unique_id: p.pos for p in model.current_patients}, positions)
        self.assertTrue(site.colonized)
        self.assertEqual(model.response_counts["workflow_events"]["occupied_destination_rejected"], 1)
        self.assertEqual(model.response_counts["environment_events"]["relocation_clear"], 0)
        model.assert_response_invariants()

    def test_one_empty_isolation_bed_cannot_be_allocated_twice(self):
        model = self.model()
        destination = model.isol_beds[0]
        self.discharge(model, self.patient_at(model, destination))
        first = self.patient_at(model, model.shared_beds[0])
        second = self.patient_at(model, model.shared_beds[1])
        second_origin = second.pos
        self.mark_positive(first)
        self.mark_positive(second)
        model._relocate_detected_patients()
        model.refresh_beds()
        self.assertEqual(first.pos, destination.pos)
        self.assertEqual(second.pos, second_origin)
        self.assertEqual(len({p.pos for p in model.current_patients}), 29)
        self.assertEqual(model.response_counts["workflow_events"]["isolation_relocation"], 1)
        self.assertEqual(model.response_counts["workflow_events"]["isolation_deferred_no_empty_bed"], 1)
        model.assert_response_invariants()

    def test_readmission_uses_already_cleared_site(self):
        model = self.model(prob_new_patient=1)
        bed = model.shared_beds[0]
        old = self.patient_at(model, bed)
        site = self.dirty_site(model, bed)
        self.discharge(model, old)
        self.assertFalse(site.colonized)
        bed.step()
        model.refresh_beds()
        self.assertTrue(bed.filled)
        self.assertNotEqual(self.patient_at(model, bed).unique_id, old.unique_id)
        self.assertFalse(site.colonized)
        self.assertEqual(len(model.current_patients), 30)
        self.assertEqual(model.cumul_patients, 31)
        model.assert_response_invariants()

    def test_full_step_orders_discharge_before_relocation_and_reconciles_census(self):
        model = self.model(prob_transmission=0, prob_new_patient=0)
        source, destination = model.shared_beds[0], model.isol_beds[0]
        patient = self.patient_at(model, source)
        self.mark_positive(patient)
        departing = self.patient_at(model, destination)
        departing.stay = 1
        origin_site = self.dirty_site(model, source)
        destination_site = self.dirty_site(model, destination)
        model.step()
        self.assertEqual(model.discharged_patients, 1)
        self.assertEqual(len(model.current_patients), 29)
        self.assertEqual(patient.pos, destination.pos)
        self.assertEqual(model.empty_beds, {source.pos})
        self.assertFalse(origin_site.colonized)
        self.assertFalse(destination_site.colonized)
        self.assertEqual(model.response_environment_count, 0)
        self.assertEqual(model.response_environment_ticks, 2)  # Tick-start burden.
        self.assertEqual(model.response_counts["environment_events"]["discharge_clear"], 1)
        self.assertEqual(model.response_counts["environment_events"]["relocation_clear"], 1)
        self.assertEqual(model.response_counts["workflow_events"]["isolation_relocation"], 1)
        model.assert_response_invariants()

    def test_scheduled_cleaning_remains_separate(self):
        model = self.model()
        bed = model.shared_beds[0]
        site = self.dirty_site(model, bed)
        site.clean_tick = 1
        site.step()
        self.assertFalse(site.colonized)
        self.assertEqual(site.clean_tick, model.cleaningDay * model.ticks_in_day)
        self.assertEqual(model.response_counts["environment_events"]["scheduled_clear"], 1)
        self.assertEqual(model.response_counts["environment_events"]["dirty_scheduled_clear"], 1)
        self.assertEqual(model.response_counts["environment_events"]["discharge_clear"], 0)
        model.assert_response_invariants()

    def test_reference_kernel_and_frozen_files_are_unchanged(self):
        for name, expected in variant.instrumentation.CORE_SHA256.items():
            self.assertEqual(hashlib.sha256((ABM / "src" / name).read_bytes()).hexdigest(), expected)
        self.assertIs(variant.ReferenceModel._remove_discharged, variant.instrumentation.core.CPE_Model._remove_discharged)
        self.assertIs(variant.ReferenceModel._relocate_detected_patients, variant.instrumentation.core.CPE_Model._relocate_detected_patients)
        reference = variant.ReferenceModel(**self.kwargs(initial_environment_contaminated=True))
        source, destination = reference.shared_beds[0], reference.isol_beds[0]
        positive = self.patient_at(reference, source)
        other = self.patient_at(reference, destination)
        positive.positive = True
        reference._relocate_detected_patients()
        reference.refresh_beds()
        self.assertEqual(positive.pos, destination.pos)
        self.assertEqual(other.pos, source.pos)  # Original occupied-bed swap retained.
        site = next(s for s in reference.response_sites if reference.care_zone_by_position[s.pos] == source.unique_id)
        self.discharge(reference, other)
        self.assertTrue(site.colonized)  # Original vacancy carryover retained.
        reference.assert_response_invariants()

    def test_without_vacancy_or_detection_variant_preserves_all_draws_and_states(self):
        reference = variant.ReferenceModel(**self.kwargs(prob_transmission=0))
        model = self.model(prob_transmission=0)
        for _ in range(300):
            reference.step()
            model.step()
        self.assertEqual(model.rng.bit_generator.state, reference.rng.bit_generator.state)
        self.assertEqual(model.random.getstate(), reference.random.getstate())
        for old, new in zip(reference.schedule.agents, model.schedule.agents):
            self.assertEqual((old.unique_id, old.pos, old.colonized), (new.unique_id, new.pos, new.colonized))
            if old.isHCW:
                self.assertEqual(old.hygiene_rng.bit_generator.state, new.hygiene_rng.bit_generator.state)
        self.assertEqual(runner.parent.snapshot(model), runner.parent.snapshot(reference))
        model.assert_response_invariants()

    def test_deepcopy_and_existing_pair_runner_preserve_variant(self):
        model = self.model()
        copied = deepcopy(model)
        bed = copied.shared_beds[0]
        site = self.dirty_site(copied, bed)
        self.discharge(copied, self.patient_at(copied, bed))
        self.assertFalse(site.colonized)
        self.assertEqual(len(model.current_patients), 30)
        self.assertEqual(model.discharged_patients, 0)
        copied.assert_response_invariants()
        config = runner.normalize_config(dict(
            seeds=[77], warmup_days=0, observation_days=1 / 24,
            model_overrides={"prob_patient_sick": 0, "initial_environment_contaminated": False},
            scenarios=[dict(id="implementation_check", p=0.2, h0=0.9, h1=0.95)],
        ))
        result = runner.simulate_pair(runner.build_jobs(config)[0], config, model_factory=variant.WorkflowRobustnessModel)
        self.assertEqual(result["status"], "complete", result.get("traceback"))
        for arm in ("baseline", "policy"):
            self.assertIn("workflow_events", result["arms"][arm]["observation"]["response_counts"])
        self.assertEqual(result["distinct_observation_trajectories"], 2)


if __name__ == "__main__":
    unittest.main()
