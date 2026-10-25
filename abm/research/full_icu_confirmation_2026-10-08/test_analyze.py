"""Synthetic analysis/identity checks; no scientific simulator/outcomes used."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('_confirmation_analyzer_tests',HERE/'analyze.py')
analyzer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analyzer)
b = analyzer.builder


def totals(count=3,pd=1000.,days=56):
    return dict(acquired_colonizations=count,imported_colonizations=0,patient_days=pd,
        susceptible_patient_days=.9*pd,environment_contaminated_site_days=days,days=days,
        response_counts=dict(acquisitions_by_route={'patient_pickup':count},
            acquisitions_by_opportunities={'2':count},acquisitions_by_route_and_opportunities={'patient_pickup|2':count}))


def blocks(summary):
    result = []
    count = summary['acquired_colonizations']
    for day in range(summary['days']):
        acquired = count//2 if day == 0 else count-count//2 if day == summary['days']//2 else 0
        row = totals(acquired,summary['patient_days']/summary['days'],1)
        row.update(day_start=day,day_end=day+1)
        result.append(row)
    return result


def fake_row(job,config,baseline=3,policy=1,pd0=1000.,pd1=1000.):
    arms = {}
    for name,count,pd,h in (('baseline',baseline,pd0,job['scenario']['h0']),('policy',policy,pd1,job['scenario']['h1'])):
        observation = totals(count,pd)
        arms[name] = dict(observation=observation,observation_blocks=blocks(observation),h=h,
                          final_snapshot=dict(total_ticks=config['warmup_ticks']+config['observation_ticks']))
    shared = job['scenario']['h0'] == job['scenario']['h1']
    if shared:
        arms['policy'] = dict(shared_with='baseline')
    warmup = totals(0,2700.,90)
    return dict(job,status='complete',model_kwargs=analyzer.expected_model_kwargs(config,job),arms=arms,
        warmup_snapshot=dict(total_ticks=config['warmup_ticks']),warmup=warmup,warmup_blocks=blocks(warmup),
        executed_ticks_per_trajectory=config['warmup_ticks']+config['observation_ticks'],warmup_executions=1,
        distinct_observation_trajectories=1 if shared else 2,execution=dict(python='synthetic',platform='fixture'))


def small_book(module,seeds=(1,2,3,4)):
    a,banks_b = list(seeds[:2]),list(seeds[2:])
    return dict(seed_blocks={module:dict(A=a,B=banks_b)},
        allocations={module:dict(all=list(seeds),local=[a[0],banks_b[0]],cluster=[a[1],banks_b[1]])})


def cases(module):
    if module == 'calibration':
        return [b.scenario(f'independent_calibration_p{i}',p,h1=.9) for i,p in enumerate(b.GRID)]
    if module == 'holdout':
        return [b.scenario(name+'_target_candidate',.15 if name == 'independent' else b.FIXED_P[name],name)
                for name in ('independent','persistent','imperfect')]
    return [b.scenario(f'independent_low{label}_env{env}',p,env=float(env))
            for label,p in zip(('075','100','125'),b.OLD_GRID) for env in (1,0)]


def records_for(module,seeds=(1,2,3,4)):
    raw = b.config(cases(module),list(seeds),'synthetic only')
    config = analyzer.normalized(raw)
    records = {}
    for j,scenario in enumerate(config['scenarios']):
        for i,seed in enumerate(seeds):
            policy = i%2+j if module == 'holdout' else i%2 if j%2 else 1+i%2
            row = fake_row(dict(job_id=j*len(seeds)+i,seed=seed,scenario=scenario),config,3,policy)
            records[(scenario['id'],seed)] = row
    return raw,records


def write_bank(directory,raw,sources):
    config = analyzer.normalized(raw)
    jobs = analyzer.backend.build_jobs(config)
    manifest = dict(schema_version=1,config=raw,resolved_config=config,jobs=jobs,sources=sources,
                    software=dict(python='synthetic'),task_count=1)
    manifest['identity'] = analyzer.hashlib.sha256(analyzer.canonical(manifest).encode()).hexdigest()
    directory.mkdir()
    (directory/'completed').mkdir()
    (directory/'manifest.json').write_text(json.dumps(manifest))
    for job in jobs:
        row = fake_row(job,config)
        row['identity'] = manifest['identity']
        (directory/'completed'/f"job_{job['job_id']:08d}.json").write_text(json.dumps(row))
    return manifest


class BootstrapTests(unittest.TestCase):
    def test_ratio_of_sums_differs_from_mean_individual_rates(self):
        values = np.zeros((1,2,10))
        values[0,0] = [1,10,9,1,20,18,0,5,1,5]
        values[0,1] = [9,90,81,3,180,162,4,45,5,45]
        point,_ = analyzer.bootstrap(values,20,42,3)
        self.assertAlmostEqual(point[0,0],100.)
        self.assertAlmostEqual(point[0,1],20.)
        self.assertAlmostEqual(point[0,2],.8)
        self.assertAlmostEqual(point[0,3],.6)
        self.assertAlmostEqual(point[0,4],3.)

    def test_common_seed_vectors_cancel_profile_noise_and_chunk_size(self):
        raw,records = records_for('holdout')
        ids = [s['id'] for s in raw['scenarios']]
        values = np.asarray([analyzer.condition_values([records[(identifier,seed)] for seed in raw['seeds']]) for identifier in ids])
        point,draws = analyzer.bootstrap(values,200,123,1)
        _,chunked = analyzer.bootstrap(values,200,123,31)
        np.testing.assert_array_equal(draws,chunked)
        np.testing.assert_allclose(draws[:,0,2]-draws[:,1,2],np.full(200,1/3),atol=1e-14)
        self.assertGreater(np.std(draws[:,0,2]),0)
        result = analyzer.holdout_results(ids,point,draws,True)
        self.assertTrue(result['baseline_matching']['established'])
        self.assertEqual(len(result['baseline_matching']['numerical_gates']),6)
        interval = result['primary_profile_F_contrasts'][0]['F_difference']['uncertainty']['intervals'][format(analyzer.PRIMARY_LEVEL,'.12g')]
        self.assertAlmostEqual(interval[0],interval[1])

    def test_undefined_bootstrap_draw_fails_gate_without_filtering(self):
        point = np.array([[3.,2.,.3],[3.,2.,.3],[3.,2.,.3]])
        samples = np.broadcast_to(point,(20,3,3)).copy()
        samples[0,0,0] = np.nan
        ids = [name+'_target_candidate' for name in ('independent','persistent','imperfect')]
        result = analyzer.holdout_results(ids,point,samples,True)
        gate = result['baseline_matching']['numerical_gates'][0]
        self.assertFalse(gate['passed'])
        self.assertFalse(result['baseline_matching']['established'])
        self.assertEqual(gate['uncertainty']['undefined_bootstrap_draws'],1)
        self.assertIsNone(gate['uncertainty']['interval95'])
        self.assertFalse(gate['uncertainty']['conditional_on_defined_draws'])
        self.assertFalse(analyzer.holdout_results(ids,point,np.broadcast_to(point,(20,3,3)),False)['baseline_matching']['established'])

    def test_environment_sign_interactions_and_two_banks(self):
        raw,records = records_for('environment')
        book = small_book('environment')
        report = analyzer.analyze_records(records,raw,book,'environment',200,44,17)
        effects = report['primary_environment_effects']
        self.assertEqual(len(effects),3)
        for effect in effects:
            self.assertLess(effect['full_minus_disabled']['F']['estimate'],0)
        self.assertEqual(report['high_minus_low_F_effect']['estimate'],0)
        self.assertEqual(report['center_minus_mean_endpoints_F_effect']['estimate'],0)
        self.assertEqual(len(report['two_bank_diagnostics']['B_minus_A_environment_F_effects']),3)
        json.dumps(report,allow_nan=False)

    def test_bank_resampling_is_independent_not_paired_between_banks(self):
        raw,records = records_for('holdout')
        report = analyzer.analyze_records(records,raw,small_book('holdout'),'holdout',200,44,17,
                                         dict(valid_target_calibration=True))
        difference = report['two_bank_diagnostics']['B_minus_A']['independent_target_candidate']['F']
        self.assertAlmostEqual(difference['estimate'],0)
        low,high = difference['uncertainty']['interval95']
        self.assertLess(low,0)
        self.assertGreater(high,0)

    def test_first_last_windows_and_zero_early_denominator(self):
        raw,records = records_for('holdout')
        row = records[('independent_target_candidate',1)]
        for block in row['arms']['baseline']['observation_blocks'][:28]:
            block['acquired_colonizations'] = 0
        values = analyzer.condition_values([row])
        point = analyzer.metrics(values[None,:,:].sum(axis=1),1)[0]
        self.assertEqual(point[8],0)
        self.assertTrue(np.isnan(point[11]))

    def test_calibration_schema_feeds_builder_and_reports_failure_draws(self):
        proposal = json.loads((HERE/'candidate_design.json').read_text())
        book = b.seedbook(proposal)
        seeds = book['allocations']['calibration']['all']
        raw = b.config(cases('calibration'),seeds,'synthetic calibration')
        config = analyzer.normalized(raw)
        records = {}
        for job in analyzer.backend.build_jobs(config):
            index = next(i for i,s in enumerate(config['scenarios']) if s['id'] == job['scenario']['id'])
            records[(job['scenario']['id'],job['seed'])] = fake_row(job,config,4+index,4+index,1600.,1600.)
        report = analyzer.analyze_records(records,raw,book,'calibration',100,91,20)
        selection = b.select_baseline_report(report,seeds)
        self.assertTrue(selection['valid_target_calibration'])
        self.assertLess(b.GRID[0],selection['candidate_p'])
        self.assertLess(selection['candidate_p'],b.GRID[1])
        self.assertEqual(report['calibration']['selection'],selection)
        point = np.zeros((3,12)); point[:,0] = [1.,1.5,2.]
        samples = np.broadcast_to(point,(10,3,12)).copy()
        failure = analyzer.calibration_results(point,samples)
        self.assertEqual(failure['selection']['candidate_p'],.15)
        self.assertEqual(failure['draws_without_exact_hit_or_bracket'],10)
        self.assertIsNone(failure['candidate_p_uncertainty_conditional_valid']['interval95'])
        self.assertEqual(failure['fallback_inclusive_candidate_p_uncertainty']['interval95'],[.15,.15])


class StrictInputTests(unittest.TestCase):
    def test_confirmation_receipt_binds_exact_freeze_config_scope_and_partition(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            directory = root/'results/calibration_local'; directory.mkdir(parents=True)
            config_path = root/'configs/calibration_local.json'; config_path.parent.mkdir()
            config = dict(seeds=[1,3],description='synthetic')
            config_path.write_text(json.dumps(config))
            seedbook_path = root/'seedbook.json'; seedbook_path.write_text('{}')
            freeze = dict(identity_sha256='f'*64,files={'configs/calibration_local.json':analyzer.digest(config_path)})
            freeze_path = root/'freeze.json'; freeze_path.write_text(json.dumps(freeze))
            manifest = dict(config=config,task_count=6)
            context = dict(root=root,freeze=freeze,freeze_path=freeze_path,seedbook_path=seedbook_path,
                           module='calibration',output_scopes=dict(banks={'calibration_local':'results/calibration_local'}))
            receipt = dict(schema_version=1,study=HERE.name,freeze_sha256=analyzer.digest(freeze_path),
                freeze_identity='f'*64,config_path='configs/calibration_local.json',config_sha256=analyzer.digest(config_path),
                seedbook_sha256=analyzer.digest(seedbook_path),output_scope='results/calibration_local',task_count=6,
                backend='research/full_icu_response_2026-10-08/array_worker.py')
            receipt_path = directory/'confirmation_binding.json'
            receipt_path.write_text(json.dumps(receipt))
            self.assertEqual(analyzer.verify_binding(directory,manifest,context),'local')
            for key in ('freeze_sha256','freeze_identity','seedbook_sha256','output_scope','task_count','config_sha256'):
                with self.subTest(key=key):
                    wrong = dict(receipt); wrong[key] = 'incorrect'
                    receipt_path.write_text(json.dumps(wrong))
                    with self.assertRaises(ValueError): analyzer.verify_binding(directory,manifest,context)
            receipt_path.write_text(json.dumps(receipt))
            wrong_manifest = dict(manifest,config=dict(config,seeds=[2,4]))
            with self.assertRaisesRegex(ValueError,'runtime manifest'):
                analyzer.verify_binding(directory,wrong_manifest,context)

    def fixture(self,root):
        book = small_book('calibration')
        full = b.config(cases('calibration'),book['allocations']['calibration']['all'],'full synthetic')
        sources = {'fixture.py':'a'*64}
        paths = []
        for site in ('local','cluster'):
            raw = b.config(cases('calibration'),book['allocations']['calibration'][site],site+' synthetic')
            directory = root/site
            write_bank(directory,raw,sources)
            paths.append(directory)
        return paths,full,sources,book

    def test_complete_two_site_banks_and_reject_bad_source_identity(self):
        with tempfile.TemporaryDirectory() as temp:
            args = self.fixture(Path(temp))
            records,receipts = analyzer.load_banks(*args,'calibration')
            self.assertEqual(len(records),12)
            self.assertEqual(len(receipts),2)
            with self.assertRaisesRegex(ValueError,'Scientific source'):
                analyzer.load_banks(args[0],args[1],{'fixture.py':'b'*64},args[3],'calibration')
            path = args[0][0]/'manifest.json'
            manifest = json.loads(path.read_text())
            manifest['task_count'] = 2
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,'Manifest identity'):
                analyzer.load_banks(*args,'calibration')

    def test_missing_job_wrong_model_runtime_and_daily_totals_rejected(self):
        mutations = (
            ('model kwargs',lambda row:row['model_kwargs'].update(prob_patient_sick=.05)),
            ('execution software',lambda row:row['execution'].update(python='wrong')),
            ('Executed trajectory',lambda row:row.update(executed_ticks_per_trajectory=1)),
            ('Daily totals',lambda row:row['arms']['baseline']['observation'].update(patient_days=999.)),
            ('source counters',lambda row:row['arms']['baseline']['observation']['response_counts']['acquisitions_by_route'].update(patient_pickup=5)),
        )
        for message,mutation in mutations:
            with self.subTest(message=message), tempfile.TemporaryDirectory() as temp:
                args = self.fixture(Path(temp))
                path = sorted((args[0][0]/'completed').glob('*.json'))[0]
                row = json.loads(path.read_text()); mutation(row); path.write_text(json.dumps(row))
                with self.assertRaisesRegex(ValueError,message):
                    analyzer.load_banks(*args,'calibration')
        with tempfile.TemporaryDirectory() as temp:
            args = self.fixture(Path(temp))
            next((args[0][0]/'completed').glob('*.json')).unlink()
            with self.assertRaisesRegex(ValueError,'Complete expected job'):
                analyzer.load_banks(*args,'calibration')

    def test_missing_condition_seed_and_duplicate_bank_rejected(self):
        raw,records = records_for('holdout')
        records.pop(next(iter(records)))
        with self.assertRaisesRegex(ValueError,'common condition seed'):
            analyzer.analyze_records(records,raw,small_book('holdout'),'holdout',10,selection=dict(valid_target_calibration=True))
        with tempfile.TemporaryDirectory() as temp:
            args = self.fixture(Path(temp))
            with self.assertRaisesRegex(ValueError,'Duplicate declared'):
                analyzer.load_banks([args[0][0],args[0][0]],args[1],args[2],args[3],'calibration')

    def test_proposal_seedbook_and_scenario_declaration_checked(self):
        proposal = json.loads((HERE/'candidate_design.json').read_text())
        book = b.seedbook(proposal)
        raw = b.config(cases('calibration'),book['allocations']['calibration']['all'],'synthetic')
        analyzer.validate_design('calibration',raw,book,proposal)
        wrong = deepcopy(raw); wrong['scenarios'][0]['p'] = .16
        with self.assertRaisesRegex(ValueError,'prespecified family'):
            analyzer.validate_design('calibration',wrong,book,proposal)
        wrong_book = deepcopy(book); wrong_book['allocations']['calibration']['all'][0] = 8
        with self.assertRaisesRegex(ValueError,'Seedbook differs'):
            analyzer.validate_design('calibration',raw,wrong_book,proposal)


if __name__ == '__main__':
    unittest.main()
