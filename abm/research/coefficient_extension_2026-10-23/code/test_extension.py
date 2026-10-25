"""Tests of the coefficient-extension study code (no file outside the study is written).

Run from the study folder:
    PYTHONDONTWRITEBYTECODE=1 ../../.venv/bin/python -B code/test_extension.py -v

Fast tests: transition-defined roots, F_geo against closed forms and the
round-1 values, the pair table, the gate (negative controls on archived data),
the projection code check and the runner's plan/seal logic.

Record tests (skipped when the records are absent):
  * gate test on 2 seeds: gate_checks/benchmark/records holds one seed whose
    archived original record is from the local site and one from the cluster;
    the gate is recomputed independently and every archived comparison is
    redone here with plain dict equality;
  * determinism: gate_checks/determinism/records/<seed>.json, a second run of a
    benchmark seed in a fresh process, must equal the benchmark record bitwise
    in every new (x_) arm and in every other science field.
"""
import sys

sys.dont_write_bytecode = True

import copy  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
from decimal import Decimal  # noqa: E402
import decimal  # noqa: E402
from pathlib import Path  # noqa: E402
import unittest  # noqa: E402

import numpy as np  # noqa: E402

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from adapter import wb  # noqa: E402
import closed_forms as cf  # noqa: E402
import gate  # noqa: E402
import pairs  # noqa: E402
import projection  # noqa: E402
import run_batch  # noqa: E402

STUDY = CODE.parent
BENCHMARK_RECORDS = STUDY / 'gate_checks/benchmark/records'
DETERMINISM_RECORDS = STUDY / 'gate_checks/determinism/records'
ROUND1 = projection.checked_json(projection.ROUND1, projection.ROUND1_SHA256)
ULP = 2.0 ** -52


class TransitionRootTests(unittest.TestCase):
    def test_rho0_is_one_minus_sqrt_miss(self):
        with decimal.localcontext() as context:
            context.prec = 60
            for target, miss in (('0.90', Decimal('0.10')), ('0.95', Decimal('0.05'))):
                self.assertLess(abs(cf.transition_h_decimal(0., target) - (1 - miss.sqrt())), Decimal('1e-50'))
                self.assertEqual(cf.transition_h(0., target), float(1 - miss.sqrt()))
        self.assertLessEqual(abs(cf.transition_h(0., '0.90') - (1 - math.sqrt(.1))), 2 * ULP)
        self.assertLessEqual(abs(cf.transition_h(0., '0.95') - (1 - math.sqrt(.05))), 2 * ULP)

    def test_rho_half_is_the_quadratic_root(self):
        # m = 1-h solves .5 m^2 + .5 m - miss = 0, i.e. m = (-1 + sqrt(1 + 8 miss)) / 2.
        with decimal.localcontext() as context:
            context.prec = 60
            for target, miss in (('0.90', Decimal('0.10')), ('0.95', Decimal('0.05'))):
                m = (-1 + (1 + 8 * miss).sqrt()) / 2
                self.assertLess(abs(cf.transition_h_decimal(.5, target) - (1 - m)), Decimal('1e-50'))

    def test_roots_give_the_transition_compliance_and_match_protocol_text(self):
        expected = {0.: (.683772, .776393), .5: (.829180, .908392)}
        for rho, (h0, h1) in expected.items():
            for target, h6 in (('0.90', h0), ('0.95', h1)):
                h = cf.transition_h(rho, target)
                self.assertEqual(round(h, 6), h6)
                self.assertLess(abs(cf.transition_compliance(h, rho) - float(target)), 4 * ULP)
                record = pairs.TRANSITION_ROOTS[f'rho{rho:g}']['h0' if target == '0.90' else 'h1']
                self.assertEqual(record['h_double'], h)
                self.assertLess(abs(record['transition_miss_error_at_double_exact_arithmetic']), 1e-16)

    def test_transition_pairs_use_the_recorded_roots(self):
        self.assertEqual(len(pairs.TRANSITION), 4)
        for p in pairs.TRANSITION:
            roots = pairs.TRANSITION_ROOTS[f'rho{p["rho"]:g}']
            self.assertEqual((p['h0'], p['h1']), (roots['h0']['h_double'], roots['h1']['h_double']))
            self.assertIn(p['e'], (1., .8))

    def test_isolated_two_opportunity_transition_difference_is_zero_at_e1(self):
        # q2 at e=1 is the transition miss probability, so F2 = 1 - .05/.10 = .5 for every rho.
        values = []
        for rho in (0., .5):
            roots = pairs.TRANSITION_ROOTS[f'rho{rho:g}']
            values.append(cf.Fk(2, roots['h0']['h_double'], roots['h1']['h_double'], rho, 1.))
        for value in values:
            self.assertAlmostEqual(value, .5, places=14)
        self.assertLess(abs(values[1] - values[0]), 1e-14)
        self.assertAlmostEqual(ROUND1['transition_policy_F2_e1_pct'], 50.0, places=12)
        self.assertEqual(ROUND1['per_boundary_h_for_transition_compliance_independent']['0.90'], 1 - math.sqrt(.1))


class FgeoTests(unittest.TestCase):
    def test_e1_closed_form(self):
        # e = 1: q_k = (1-h) p^(k-1), p = rho + (1-rho)(1-h)  =>  S(h) = (1-h) p / (1 - p^2).
        def S(h, rho):
            p = rho + (1 - rho) * (1 - h)
            return (1 - h) * p / (1 - p * p)
        for rho in (0., .1, .25, .5, .75, .95):
            for h0, h1 in ((.9, .95), (.6, .7), (.8, .9), (.95, .99)):
                self.assertAlmostEqual(cf.Fgeo(h0, h1, rho, 1.), 1 - S(h1, rho) / S(h0, rho), places=13)
                self.assertAlmostEqual(cf.sum_even(h0, rho, 1.), S(h0, rho), places=13)
        # Independent misses, perfect clearing: 1 - [.0025/.9975]/[.01/.99].
        self.assertAlmostEqual(cf.Fgeo(.9, .95, 0., 1.), 1 - (.0025 / .9975) / (.01 / .99), places=14)

    def test_truncated_direct_sum_over_fine_grid(self):
        for rho in (0., .3, .5, .95):
            for e in (.5, .78, .9, 1.):
                direct = 1 - (sum(cf.qk(2 * j, .95, rho, e) for j in range(1, 3000))
                              / sum(cf.qk(2 * j, .9, rho, e) for j in range(1, 3000)))
                self.assertLess(abs(direct - cf.Fgeo(.9, .95, rho, e)), 1e-10)

    def test_round1_values(self):
        renewal = ROUND1['renewal_comparator']
        for pid, row in renewal['rows'].items():
            p = next(q for q in pairs.ARCHIVED if q['id'] == pid)
            self.assertEqual(100 * cf.Fgeo(p['h0'], p['h1'], p['rho'], p['e']), row['F_geo_pct'])
            self.assertEqual(100 * cf.Fk(2, p['h0'], p['h1'], p['rho'], p['e']), row['F_q2_pct'])
        self.assertEqual(len(renewal['rows']), 22)
        for name, (rho, e) in (('independent', (0., 1.)), ('clustered', (.5, 1.)), ('imperfect', (0., .8))):
            self.assertEqual(100 * cf.Fgeo(.9, .95, rho, e), renewal['primary_pct'][name])
        self.assertEqual(cf.estar(.9, .95), ROUND1['crossings']['isolated_k2_closed_form'])

    def test_vector_form_equals_matrix_form_on_fine_grid(self):
        for rho in pairs.FINE_RHO:
            for e in pairs.FINE_E:
                self.assertLess(abs(float(cf.Fgeo_vec(.9, .95, rho, e)) - cf.Fgeo(.9, .95, rho, e)), 1e-12)

    def test_isolated_completion_range_values_in_protocol(self):
        for (h0, h1), pp in (((.6, .7), -3.5), ((.7, .8), -2.4), ((.8, .9), 1.2)):
            self.assertEqual(round(100 * (cf.Fk(2, h0, h1, .5, .8) - cf.Fk(2, h0, h1, 0., .8)), 1), pp)


class PairTests(unittest.TestCase):
    def test_structure(self):
        self.assertEqual(pairs.ALL[:22], tuple(dict(p) for p in wb.SURFACE_PAIRS))
        self.assertEqual(len(pairs.ALL), 182)
        self.assertEqual(len({p['id'] for p in pairs.ALL}), 182)
        self.assertTrue(all(p['id'].startswith('x_') for p in pairs.NEW))
        self.assertEqual([len(g) for _, g in pairs.GROUPS], [14, 10, 4, 132])
        self.assertEqual(sorted({(p['h0'], p['h1']) for p in pairs.COMPLETION}),
                         [(.6, .7), (.7, .8), (.8, .9), (.95, .99)])
        self.assertEqual({(p['rho'], p['e']) for p in pairs.COMPLETION if (p['h0'], p['h1']) in ((.8, .9), (.95, .99))},
                         {(.5, .8)})
        self.assertEqual({(p['rho'], p['e']) for p in pairs.FINE_GRID},
                         {(r, e) for r in (0., .1, .2, .3, .4, .5, .6, .7, .8, .9, .95)
                          for e in (.5, .6, .7, .75, .78, .8, .82, .85, .9, .95, .98, 1.)})
        self.assertEqual({(p['rho'], p['e']) for p in pairs.CROSSING},
                         {(r, e) for r in (0., .5) for e in (.91, .92, .93, .94, .95, .96, .98)})

    def test_duplicates(self):
        archived = dict(pairs.ARCHIVED_DUPLICATES)
        self.assertEqual(len(archived), 8)
        for new_id, old_id in archived.items():
            new = next(p for p in pairs.NEW if p['id'] == new_id)
            old = next(p for p in pairs.ARCHIVED if p['id'] == old_id)
            self.assertEqual((new['h0'], new['h1'], new['rho'], new['e']), (old['h0'], old['h1'], old['rho'], old['e']))
        self.assertEqual(sum(1 for _, old in pairs.ARCHIVED_DUPLICATES if old in pairs.LEGACY), 3)
        self.assertEqual(len(pairs.NEW_DUPLICATES), 4)

    def test_adapter_accepts_the_surface_with_default_cap(self):
        _, config = wb.archived_scientific_config()
        extra, legacy = wb._validate_surface(pairs.ALL, config['profiles'], config['h0'], config['h1'])
        self.assertEqual(len(extra), 179)
        self.assertEqual(len(legacy), 3)
        self.assertTrue(all(p['histogram_cap'] is None for p in extra))
        self.assertEqual(sum(2 * len(pairs.KINDS[k]) for k in ('all',)), 364)


def _row_from_extended(seed):
    """A record-shaped row assembled from the archived extended record (what a correct run produces)."""
    ext = gate.load_extended(seed)['row']
    row = dict(ordered_bridge_events_sha256=ext['ordered_bridge_events_sha256'],
               profiles={legacy: copy.deepcopy({k: ext['pairs'][pid][k] for k in gate.ARMS})
                         for legacy, pid in gate.PRIMARY_PAIR_IDS.items()},
               pairs={pid: copy.deepcopy(ext['pairs'][pid]) for pid in pairs.ARCHIVED_IDS})
    for new_id, old_id in pairs.ARCHIVED_DUPLICATES:
        arms = copy.deepcopy(ext['pairs'][old_id])
        if old_id in pairs.LEGACY:
            for arm in gate.ARMS:
                for key in pairs.HISTOGRAM_SENSITIVE_KEYS:
                    arms[arm][key] = None
        row['pairs'][new_id] = arms
    for cross_id, fine_id in pairs.NEW_DUPLICATES:
        row['pairs'][fine_id] = copy.deepcopy(ext['pairs']['h0.9_to_0.95_rho0.5_e0.9'])
        row['pairs'][cross_id] = copy.deepcopy(row['pairs'][fine_id])
    for key in gate.EXPOSURE_FIELDS:
        row[key] = ext[key]
    return row, ext


class GateUnitTests(unittest.TestCase):
    seed = 3080160514

    def test_correct_row_passes(self):
        row, ext = _row_from_extended(self.seed)
        result = gate.run_gate(row, ext)
        self.assertTrue(result['passed'])
        self.assertEqual(result['counts'], dict(primary_profiles_equal=3, archived_pairs_equal=22,
                                                archived_duplicates_checked=8, archived_duplicates_equal=8,
                                                internal_duplicates_checked=4, internal_duplicates_equal=4))

    def test_one_ulp_in_one_archived_pair_fails(self):
        row, ext = _row_from_extended(self.seed)
        arm = row['pairs']['h0.9_to_0.95_rho0.75_e0.6']['policy']
        arm['K'] = math.nextafter(arm['K'], math.inf)
        result = gate.run_gate(row, ext)
        self.assertFalse(result['passed'])
        self.assertEqual(result['counts']['archived_pairs_equal'], 21)

    def test_negative_zero_and_hash_and_profile_fail(self):
        row, ext = _row_from_extended(self.seed)
        self.assertEqual(row['profiles']['iid_perfect']['baseline']['K_from_warmup_pickups'], 0.0)
        row['profiles']['iid_perfect']['baseline']['K_from_warmup_pickups'] = -0.0
        result = gate.run_gate(row, ext)
        self.assertFalse(result['passed'])
        self.assertEqual(result['counts']['primary_profiles_equal'], 2)
        row, ext = _row_from_extended(self.seed)
        row['ordered_bridge_events_sha256'] = '0' * 64
        self.assertFalse(gate.run_gate(row, ext)['passed'])

    def test_duplicate_mismatches_fail(self):
        row, ext = _row_from_extended(self.seed)
        row['pairs']['x_fine_h0.9_to_0.95_rho0.5_e0.6']['baseline']['K_by_worker']['-5'] += 1e-9
        result = gate.run_gate(row, ext)
        self.assertFalse(result['passed'])
        self.assertEqual(result['counts']['archived_duplicates_equal'], 7)
        row, ext = _row_from_extended(self.seed)
        # A primary duplicate is compared without the histogram keys only because they are None.
        row['pairs']['x_fine_h0.9_to_0.95_rho0_e1']['policy']['K_by_opportunities'] = {'0': 0.0}
        self.assertFalse(gate.run_gate(row, ext)['passed'])
        row, ext = _row_from_extended(self.seed)
        row['pairs']['x_cross_h0.9_to_0.95_rho0.5_e0.98']['policy']['K'] *= 1 + ULP
        result = gate.run_gate(row, ext)
        self.assertFalse(result['passed'])
        self.assertEqual(result['counts']['internal_duplicates_equal'], 3)

    def test_extended_record_digest_is_verified(self):
        loaded = gate.load_extended(self.seed)
        self.assertEqual(loaded['file_sha256'], wb.sha256(loaded['path']))
        self.assertEqual(loaded['row']['seed'], self.seed)


class ProjectionTests(unittest.TestCase):
    def test_code_check_reproduces_round1(self):
        check = projection.code_check()
        self.assertTrue(check['passed'], json.dumps(check['by_truth'], indent=1))

    def test_grid_surface_interpolates_and_refuses_extrapolation(self):
        R, E = np.array(pairs.FINE_RHO), np.array(pairs.FINE_E)
        Z = np.add.outer(2 * R, 3 * E)  # bilinear interpolation is exact for a bilinear function
        grid = projection.GridSurface(R, E, Z)
        for i, r in enumerate(R):
            for j, e in enumerate(E):
                self.assertEqual(float(grid(np.array([r]), e)[0]), Z[i, j])
        rho = np.linspace(0, .95, 37)
        for e in (.5, .63, .779, .95, 1.):
            np.testing.assert_allclose(grid(rho, e), 2 * rho + 3 * e, rtol=0, atol=1e-14)
        with self.assertRaises(ValueError):
            grid(np.array([.951]), .8)
        with self.assertRaises(ValueError):
            grid(np.array([.5]), .499)

    def test_out_of_domain_sets_are_counted_not_evaluated(self):
        R, E = np.array(pairs.FINE_RHO), np.array(pairs.FINE_E)
        grid = projection.GridSurface(R, E, np.zeros((len(R), len(E))))
        # At e = .8: b = (1 - .72)^2 = .0784, v = .64 x .09 = .0576, so q2 = b + rho v.
        b, v = .0784, .0576
        limits = np.array([[.60, .62, .0, 1.],                    # q1 forces e in [.42, .44] < .5
                           [.279, .281, b + .9 * v, b + 1. * v],  # rho up to about 1 > .95
                           [.279, .281, b + .3 * v, b + .5 * v],  # inside the domain
                           [.279, .281, b + 3. * v, b + 4. * v]])  # rho interval above 1: every slice empty
        result = projection.project(limits, grid, domain=grid)
        self.assertEqual(result['out_of_domain'], dict(total=2, rho_above_max=2, e_below_min=1))
        self.assertEqual((result['defined'], result['empty']), (1, 1))


class RunnerTests(unittest.TestCase):
    def test_plans(self):
        full, fallback = run_batch.plan('full'), run_batch.plan('fallback')
        self.assertEqual(len(full), 256)
        self.assertEqual({k for _, _, k in full}, {'all'})
        self.assertEqual({(b, k) for _, b, k in fallback}, {('A', 'all'), ('B', 'no_fine_grid')})
        self.assertEqual([b for _, b, _ in full[:4]], ['A', 'B', 'A', 'B'])
        self.assertEqual(len({s for s, _, _ in full}), 256)

    def test_seal_detects_tampering(self):
        paths = sorted(BENCHMARK_RECORDS.glob('[0-9]*.json'))
        if not paths:
            self.skipTest('no benchmark record yet')
        manifest = json.loads((BENCHMARK_RECORDS / 'manifest.json').read_text())
        row = run_batch.validate_record(paths[0], manifest, int(paths[0].stem), 'all')
        self.assertEqual(row['study'], run_batch.STUDY_NAME)
        tampered = copy.deepcopy(row)
        tampered['pairs']['x_fine_h0.9_to_0.95_rho0.95_e0.5']['policy']['K'] += 1.0
        temporary = STUDY / 'gate_checks' / f'.tamper_test_{paths[0].name}'
        try:
            temporary.write_text(json.dumps(tampered))
            with self.assertRaises(ValueError):
                run_batch.validate_record(temporary, manifest, row['seed'], 'all')
        finally:
            temporary.unlink(missing_ok=True)


class BytecodeTests(unittest.TestCase):
    def test_no_new_bytecode_in_frozen_trees_or_study(self):
        snapshot = STUDY / 'gate_checks/pycache_snapshot_start.txt'
        lines = snapshot.read_text().splitlines()
        before = sorted(lines[1:])
        root = wb.ABM.parent  # disease_modeling/
        # The trees listed when the snapshot was taken (before any code of this study ran).
        trees = ['abm/src', 'abm/research/full_icu_confirmation_2026-10-08', 'abm/research/full_icu_response_2026-10-08',
                 'abm/research/robustness_uq_2026-10-17', 'abm/research/workflow_coefficient_2026-10-21',
                 'abm/research/manuscript_v2_2026-10-21', 'research/author_notes_2026-10-10',
                 'research/robustness_uq_2026-10-17']
        now = []
        for tree in trees:
            for path in (root / tree).rglob('*.pyc'):
                stat = path.stat()
                now.append(f'{int(stat.st_mtime)} {stat.st_size} {path.relative_to(root)}')
        self.assertEqual(sorted(now), before)
        self.assertEqual(list(STUDY.rglob('__pycache__')), [])
        self.assertEqual(list(STUDY.rglob('*.pyc')), [])


def _records(folder):
    if not (folder / 'manifest.json').exists():
        return None, {}
    manifest = json.loads((folder / 'manifest.json').read_text())
    rows = {}
    for path in sorted(folder.glob('[0-9]*.json')):
        rows[int(path.stem)] = run_batch.validate_record(path, manifest, int(path.stem), 'all')
    return manifest, rows


class GateRecordTests(unittest.TestCase):
    """Gate test on 2 seeds: one archived from the local site, one from the cluster."""

    @classmethod
    def setUpClass(cls):
        cls.manifest, cls.rows = _records(BENCHMARK_RECORDS)
        if len(cls.rows) != 2:
            raise unittest.SkipTest('benchmark records not complete')

    def test_one_local_one_cluster(self):
        sites = set()
        for seed in self.rows:
            local = list(wb.ABM.glob(f'results/full_icu_confirmation_2026-10-08/bridge_local/part_*/completed/{seed}.json'))
            cluster = list(wb.ABM.glob(f'results/full_icu_confirmation_2026-10-08/bridge_cluster/part_*/completed/{seed}.json'))
            self.assertEqual(len(local) + len(cluster), 1)
            sites.add('local' if local else 'cluster')
        self.assertEqual(sites, {'local', 'cluster'})

    def test_stored_gate_passes_and_equals_recomputation(self):
        for seed, row in self.rows.items():
            ext = gate.load_extended(seed)
            recomputed = gate.run_gate(row, ext['row'])
            stored = {k: v for k, v in row['gate'].items() if k != 'archived_extended_record'}
            self.assertTrue(row['gate']['passed'], seed)
            self.assertEqual(gate.canonical(recomputed), gate.canonical(stored))
            self.assertEqual(row['gate']['archived_extended_record']['file_sha256'], ext['file_sha256'])
            self.assertEqual(row['gate']['counts'], dict(primary_profiles_equal=3, archived_pairs_equal=22,
                                                         archived_duplicates_checked=8, archived_duplicates_equal=8,
                                                         internal_duplicates_checked=4, internal_duplicates_equal=4))

    def test_independent_plain_comparison_with_archive(self):
        # Independent of gate.py: raw JSON text of the archived record, plain == on parsed values plus repr of K.
        for seed, row in self.rows.items():
            ext = json.loads(Path(gate.extended_record_path(seed)).read_text())
            self.assertEqual(row['ordered_bridge_events_sha256'], ext['ordered_bridge_events_sha256'])
            for legacy, pid in wb.PRIMARY_PAIR_IDS.items():
                self.assertEqual(row['profiles'][legacy], ext['pairs'][pid])
            for pid in pairs.ARCHIVED_IDS:
                for arm in ('baseline', 'policy'):
                    self.assertEqual(row['pairs'][pid][arm], ext['pairs'][pid][arm])
                    self.assertEqual(repr(row['pairs'][pid][arm]['K']), repr(ext['pairs'][pid][arm]['K']))
            for new_id, old_id in pairs.ARCHIVED_DUPLICATES:
                for arm in ('baseline', 'policy'):
                    new, old = row['pairs'][new_id][arm], ext['pairs'][old_id][arm]
                    self.assertEqual(repr(new['K']), repr(old['K']))
                    self.assertEqual({k: v for k, v in new.items() if k not in pairs.HISTOGRAM_SENSITIVE_KEYS or old_id not in pairs.LEGACY},
                                     {k: v for k, v in old.items() if k not in pairs.HISTOGRAM_SENSITIVE_KEYS or old_id not in pairs.LEGACY})
            for key in gate.EXPOSURE_FIELDS:
                self.assertEqual(row[key], ext[key], key)
            self.assertEqual(len(row['pairs']), 182)


class DeterminismTests(unittest.TestCase):
    """One seed run twice in fresh processes: every new arm bitwise identical."""

    @classmethod
    def setUpClass(cls):
        cls.manifest_a, cls.first = _records(BENCHMARK_RECORDS)
        cls.manifest_b, cls.second = _records(DETERMINISM_RECORDS)
        cls.common = sorted(set(cls.first) & set(cls.second))
        if not cls.common:
            raise unittest.SkipTest('no determinism rerun yet')

    def test_same_identity_different_processes(self):
        self.assertEqual(self.manifest_a['identity_sha256'], self.manifest_b['identity_sha256'])
        for seed in self.common:
            self.assertNotEqual((self.first[seed]['pid'], self.first[seed]['started_at_utc']),
                                (self.second[seed]['pid'], self.second[seed]['started_at_utc']))

    def test_new_arms_bitwise_identical(self):
        for seed in self.common:
            a, b = self.first[seed], self.second[seed]
            new_ids = [p['id'] for p in pairs.NEW]
            self.assertEqual(len(new_ids), 160)
            for pid in new_ids:
                self.assertEqual(gate.canonical(a['pairs'][pid]), gate.canonical(b['pairs'][pid]), pid)

    def test_every_science_field_identical(self):
        for seed in self.common:
            a, b = (run_batch.science_view(r[seed]) for r in (self.first, self.second))
            self.assertEqual(gate.canonical(a), gate.canonical(b))


if __name__ == '__main__':
    unittest.main()
