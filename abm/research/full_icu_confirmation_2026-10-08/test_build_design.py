"""Bounded protocol behavior and sample partition checks; no biological simulation."""
from copy import deepcopy
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('confirmation_design', HERE/'build_design.py')
design = importlib.util.module_from_spec(spec)
spec.loader.exec_module(design)


class SelectionTests(unittest.TestCase):
    def test_increasing_bracket_interpolates_without_extrapolation(self):
        result = design.choose(zip(design.GRID, (2., 4., 6.)))
        self.assertEqual(result['candidate_p'], (design.GRID[0]+design.GRID[1])/2)
        self.assertEqual(result['weights'], [0.5, 0.5])
        self.assertTrue(result['valid_target_calibration'])
        self.assertIsNone(result['allowed_extension'])

    def test_lowest_exact_hit_precedes_bracket_and_nonmonotonicity_is_retained(self):
        result = design.choose(zip(design.GRID, (3., 3., 2.)))
        self.assertEqual(result['candidate_p'], design.GRID[0])
        self.assertEqual(result['exact_hits'], list(design.GRID[:2]))
        self.assertTrue(result['nonmonotone'])

    def test_all_below_all_above_and_decreasing_crossing_have_fixed_fallback(self):
        for rates in ((1., 2., 2.9), (3.1, 4., 5.), (4., 3.5, 2.)):
            with self.subTest(rates=rates):
                result = design.choose(zip(design.GRID, rates))
                self.assertEqual(result['candidate_p'], .15)
                self.assertFalse(result['valid_target_calibration'])
                self.assertIsNone(result['selected_bracket'])
                self.assertIsNone(result['allowed_extension'])

    def test_nonmonotonic_curve_retains_valid_first_bracket(self):
        result = design.choose(zip(design.GRID, (2., 4., 2.)))
        self.assertTrue(result['nonmonotone'])
        self.assertTrue(result['valid_target_calibration'])

    def test_rejects_duplicate_extra_grid_negative_and_nonfinite_input(self):
        bad = [list(zip(design.GRID, (2., 4., math.nan))),
               list(zip(design.GRID, (2., 4., -1.))),
               [(design.GRID[0], 2.)]*3,
               list(zip(design.GRID, (2., 4., 5.)))+[(.2, 6.)]]
        for points in bad:
            with self.subTest(points=points):
                with self.assertRaises(ValueError): design.choose(points)

    def report(self):
        seeds = list(range(1001, 1129))
        report = {'scenarios': []}
        for i, (p, rate) in enumerate(zip(design.GRID, (2., 4., 6.))):
            report['scenarios'].append({'scenario': design.scenario(f'cal_{i}', p, h1=.9),
                'n': 128, 'seeds': seeds,
                'arms': {'baseline': {'acquired_total': int(rate*128),
                                     'patient_days_total': 128000.,
                                     'rate_per_1000_patient_days': rate}}})
        return report, seeds

    def test_selection_recomputes_baseline_rate_and_ignores_nonessential_policy_fields(self):
        report, seeds = self.report()
        report['scenarios'][0]['arms']['policy'] = {'unused': 'shared baseline'}
        selected = design.select_baseline_report(report, seeds)
        self.assertTrue(selected['valid_target_calibration'])

    def test_report_rejects_policy_stage_reused_seeds_corrupt_rate_or_missing_exposure(self):
        original, seeds = self.report()
        mutations = [lambda r: r['scenarios'][0]['scenario'].update(h1=.95),
                     lambda r: r['scenarios'][0].update(seeds=list(range(811001, 811129))),
                     lambda r: r['scenarios'][0]['arms']['baseline'].update(rate_per_1000_patient_days=2.1),
                     lambda r: r['scenarios'][0]['arms']['baseline'].update(patient_days_total=0.),
                     lambda r: r['scenarios'][0]['scenario'].update(environment_multiplier=0.)]
        for mutate in mutations:
            report = deepcopy(original); mutate(report)
            with self.assertRaises(ValueError): design.select_baseline_report(report, seeds)


class ArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.proposal = json.loads((HERE/'candidate_design.json').read_text())

    def test_banks_families_and_sites_partition_exact_sample_sizes(self):
        book = design.seedbook(self.proposal)
        all_used = []
        for family, n in design.FAMILIES.items():
            banks = book['seed_blocks'][family]
            alloc = book['allocations'][family]
            self.assertEqual(len(alloc['all']), n)
            self.assertFalse(set(alloc['local']) & set(alloc['cluster']))
            self.assertEqual(set(alloc['all']), set(alloc['local']) | set(alloc['cluster']))
            for bank in ('A', 'B'):
                self.assertEqual(len(set(banks[bank]) & set(alloc['local'])), n//16)
                self.assertEqual(len(set(banks[bank]) & set(alloc['cluster'])), 7*n//16)
            all_used.extend(alloc['all'])
        self.assertEqual(len(all_used), len(set(all_used)))

    def test_primary_counts_and_absent_independent_holdout(self):
        outputs = design.build_initial(self.proposal)
        cal = outputs['configs/calibration_all.json']
        env = outputs['configs/environment_all.json']
        self.assertEqual(len(cal['scenarios'])*len(cal['seeds']), 384)
        self.assertEqual(len(env['scenarios'])*len(env['seeds']), 3072)
        self.assertTrue(all(s['h0'] == s['h1'] == .9 for s in cal['scenarios']))
        self.assertEqual(outputs['holdout_candidates.json']['candidates']['independent']['p'], None)
        self.assertNotIn('configs/holdout_all.json', outputs)
        self.assertEqual(outputs['holdout_candidates.json']['candidates']['persistent']['p'], design.FIXED_P['persistent'])
        self.assertEqual(outputs, design.build_initial(self.proposal))

    def test_locked_holdout_uses_fresh_bank_and_retains_calibration_failure(self):
        book = design.seedbook(self.proposal)
        expected = book['allocations']['calibration']['all']
        report = {'scenarios': []}
        for i, p in enumerate(design.GRID):
            report['scenarios'].append({'scenario': design.scenario(f'cal_{i}', p, h1=.9),
                 'n': 128, 'seeds': expected,
                 'arms': {'baseline': {'acquired_total': 128, 'patient_days_total': 128000.,
                                      'rate_per_1000_patient_days': 1.}}})
        outputs = design.build_holdout(self.proposal, report)
        self.assertFalse(outputs['calibration_selection.json']['valid_target_calibration'])
        hold = outputs['configs/holdout_all.json']
        self.assertEqual(hold['scenarios'][0]['p'], .15)
        self.assertEqual(len(hold['scenarios'])*len(hold['seeds']), 1536)
        self.assertFalse(set(hold['seeds']) & set(expected))
        self.assertEqual(hold['scenarios'][1]['p'], design.FIXED_P['persistent'])
        self.assertEqual(hold['scenarios'][2]['p'], design.FIXED_P['imperfect'])

    def test_duplicate_seed_rejected_before_artifacts(self):
        changed = deepcopy(self.proposal)
        changed['seed_blocks']['holdout']['B'][0] = changed['seed_blocks']['calibration']['A'][0]
        with self.assertRaises(ValueError): design.build_initial(changed)

    def test_no_overwrite_and_preflight_prevents_partial_new_selection(self):
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp)
            design.write_outputs({'locked.json': {'p': .15}}, destination)
            with self.assertRaises(FileExistsError):
                design.write_outputs({'new.json': {}, 'locked.json': {'p': .16}}, destination)
            self.assertFalse((destination/'new.json').exists())
            self.assertEqual(json.loads((destination/'locked.json').read_text()), {'p': .15})


if __name__ == '__main__': unittest.main()
