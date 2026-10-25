"""Prospective synthetic accounting and inference checks; no study outcomes."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

spec=importlib.util.spec_from_file_location('_new_workflow_analysis',Path(__file__).with_name('analyze.py'))
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
v=a.accounting
SOFTWARE=dict(v.PACKAGES,python='3.12.3 (artificial fixture)')


def start():
    s={k:(0. if k in v.EXPOSURES else 0) for k in v.COUNTERS}
    s.update(total_ticks=0,census=30,colonized_census=2,cumul_patients=30,imported_colonizations=2,
        environment_contaminated_site_days=0.,environment_contaminated_sites_endpoint=0,
        hygiene_counts_by_role={r:dict.fromkeys(v.HYGIENE,0) for r in ('Nurse','Dr','XrayDr')},
        response_counts={k:{} for k in v.RESPONSE_KINDS})
    return s


def period(before,n,count):
    end=copy.deepcopy(before);blocks=[]
    for day in range(n):
        response={k:{} for k in v.RESPONSE_KINDS}
        response.update(pickups={'patient':2},pickups_by_completion_state={'completed':2},
            acquisitions_by_patient_source_status={'imported':count},acquisitions_by_route={'patient_pickup':count},
            acquisitions_by_opportunities={'0':count},acquisitions_by_route_and_opportunities={'patient_pickup|0':count},
            hand_at_risk_opportunities={'1':2},hand_cleared_opportunities={'1':2},
            environment_events={'scheduled_clear':30,'discharge_clear':1})
        b=dict.fromkeys(v.COUNTERS,0)
        roles={r:dict.fromkeys(v.HYGIENE,0) for r in ('Nurse','Dr','XrayDr')}
        roles['Nurse']=dict(zip(v.HYGIENE,(40,36,2,2)))
        b.update(acquired_colonizations=count,cumul_patients=1,discharged_patients=1,
            bed_visits=40,patient_contacts=40,environment_contacts=20,patient_days=30.,susceptible_patient_days=28.,
            hygiene_opportunity_count=40,hygiene_completion_count=36,hygiene_effective_clear_count=2,
            hygiene_decontamination_count=2,environment_contaminated_site_days=0.,ticks=2592,days=1.,
            day_start=float(day),day_end=float(day+1),hygiene_counts_by_role=roles,response_counts=response,
            mean_census=30.,mean_susceptible_census=28.,colonized_patient_days=2.,mean_colonized_census=2.,
            mean_environment_contaminated_fraction=0.,acquired_per_1000_patient_days=1000*count/30,
            acquired_per_1000_susceptible_patient_days=1000*count/28)
        for key in ('census','colonized_census','environment_contaminated_sites_endpoint'):
            b[key+'_start']=b[key+'_end']=before[key]
        for key in v.COUNTERS+('environment_contaminated_site_days',):end[key]+=b[key]
        end['total_ticks']+=2592
        for group in ('hygiene_counts_by_role','response_counts'):
            for category,counts in b[group].items():
                for key,value in counts.items():end[group][category][key]=end[group][category].get(key,0)+value
        blocks.append(b)
    total={key:end[key]-before[key] for key in v.COUNTERS+('environment_contaminated_site_days',)}
    total.update(ticks=2592*n,days=float(n),mean_census=30.,mean_susceptible_census=28.,
        colonized_patient_days=2.*n,mean_colonized_census=2.,mean_environment_contaminated_fraction=0.,
        acquired_per_1000_patient_days=1000*count/30,acquired_per_1000_susceptible_patient_days=1000*count/28,
        hygiene_counts_by_role={r:{k:end['hygiene_counts_by_role'][r][k]-before['hygiene_counts_by_role'][r][k]
            for k in v.HYGIENE} for r in ('Nurse','Dr','XrayDr')},
        response_counts={kind:{k:end['response_counts'][kind][k]-before['response_counts'][kind].get(k,0)
            for k in end['response_counts'][kind]} for kind in v.RESPONSE_KINDS})
    for key in ('census','colonized_census','environment_contaminated_sites_endpoint'):
        total[key+'_start']=before[key];total[key+'_end']=end[key]
    return total,blocks,end


def fixture(shared=False):
    sc=dict(id='clustered',p=.2,h0=.9,h1=.9 if shared else .95,rho=.5,e=1.,
        opportunities='zone_separate',environment_multiplier=1.)
    cfg=dict(scenarios=[sc],seeds=[123],warmup_days=90,observation_days=56,block_days=1,
             model_overrides=dict(prob_patient_sick=.03,initial_environment_contaminated=False))
    initial=start();warm,blocks,boundary=period(initial,90,1);arms={}
    for name,h,count in (('baseline',sc['h0'],1),('policy',sc['h1'],0)):
        obs,daily,end=period(boundary,56,count)
        arms[name]=dict(h=h,observation=obs,observation_blocks=daily,final_snapshot=end)
    if shared:arms['policy']={'shared_with':'baseline'}
    row=dict(job_id=0,scenario=sc,seed=123,status='complete',identity='synthetic',execution=dict(SOFTWARE,platform='Linux-fixture'),
        distinct_observation_trajectories=1 if shared else 2,warmup_executions=1,executed_ticks_per_trajectory=378432,
        model_kwargs=dict(v.MODEL,seed=123,prob_transmission=.2,icu_hcw_wash_rate=.9,outside_hcw_wash_rate=.9,
            hygiene_opportunities='zone_separate',hygiene_autocorrelation=.5,hygiene_efficacy=1.,environment_multiplier=1.),
        initial_snapshot=initial,warmup_snapshot=boundary,warmup=warm,warmup_blocks=blocks,arms=arms)
    return row,cfg


class Tests(unittest.TestCase):
    def test_valid_accounting_pair_and_baseline_only(self):
        row,cfg=fixture();self.assertEqual(a.validate_records([row],cfg).tolist(),[[[56.,1680.,0.,1680.]]])
        row,cfg=fixture(True);self.assertEqual(a.validate_records([row],cfg).tolist(),[[[56.,1680.,56.,1680.]]])

    def rejected(self,mutation,pattern):
        row,cfg=fixture();mutation(row)
        with self.assertRaisesRegex(ValueError,pattern):a.validate_records([row],cfg)

    def test_reject_omitted_extra_and_wrong_identity_records(self):
        row,cfg=fixture()
        for records in ([],[row,row]):
            with self.assertRaisesRegex(ValueError,'Complete'):a.validate_records(records,cfg)
        self.rejected(lambda r:r.update(seed=124),'declaration')
        self.rejected(lambda r:r.update(job_id=True),'identity')

    def test_runtime_and_daily_mutations_rejected(self):
        self.rejected(lambda r:r['execution'].update(numpy='2.0'),'runtime')
        self.rejected(lambda r:r['warmup_blocks'][0].update(patient_days=29.),'Exposure')
        self.rejected(lambda r:r['warmup_blocks'].pop(),'period')

    def test_workflow_reconciliation_rejects_impossible_actions(self):
        self.rejected(lambda r:r['warmup_blocks'][0]['response_counts']['environment_events'].update(discharge_clear=0),'Discharge')
        self.rejected(lambda r:r['warmup_blocks'][0]['response_counts']['workflow_events'].update(isolation_relocation=1),'Relocation')
        self.rejected(lambda r:r['warmup_blocks'][0]['response_counts']['workflow_events'].update(occupied_destination_rejected=1),'Occupied')
        self.rejected(lambda r:r['warmup_blocks'][0]['response_counts']['environment_events'].update(deposition=1),'mass balance')

    def test_select_exact_positive_log_bracket_and_failures(self):
        chosen=a.choose([(.1,1.),(.2,4.)]);self.assertAlmostEqual(chosen['selected_p'],math_sqrt3_over10())
        self.assertEqual(a.choose([(.1,3.),(.2,3.)])['selected_p'],.1)
        self.assertIsNone(a.choose([(.1,0.),(.2,4.)])['selected_p'])
        self.assertIsNone(a.choose([(.1,4.),(.2,2.)])['selected_p'])
        x=a.choose([(.1,1.),(.2,4.),(.3,2.),(.4,5.)])
        self.assertEqual(len(x['increasing_positive_adjacent_brackets']),2)
        self.assertLess(x['selected_p'],.2)
        self.assertIsNone(a.choose([(.1,1.),(.2,2.)])['selected_p'])

    def test_selection_input_invalidity(self):
        for points in ([('x',1)],[(0,1)],[(.2,1),(.1,4)],[(.1,-1)],[(.1,float('nan'))]):
            with self.assertRaises(ValueError):a.choose(points)

    def test_pooled_ratio_and_undefined(self):
        values=np.array([[2.,10.,1.,20.],[3.,40.,1.,10.]])
        got=a.estimates(values.sum(axis=0),2)
        np.testing.assert_allclose(got,[100.,200/3,1/3,100/3,1.5])
        zero=a.estimates(np.array([0.,10.,0.,10.]),1)
        self.assertTrue(np.isnan(zero[2]));uncertainty=a.interval([.2,np.nan])
        self.assertIsNone(uncertainty['interval'])
        self.assertEqual(uncertainty['undefined_draws'],1)
        self.assertEqual(uncertainty['conditional_on_defined_interval'],[.2,.2])
        self.assertIsNone(a.interval([np.nan])['conditional_on_defined_interval'])

    def test_bootstrap_preserves_joint_seed_dependence_and_streams(self):
        values=np.array([[[1,10,1,10],[2,10,2,10]],[[2,10,1,10],[4,10,2,10]],
                         [[3,10,2,10],[6,10,4,10]],[[4,10,3,10],[8,10,6,10]]],dtype=float)
        draws,h=a.bootstrap_totals(values,0,200)
        np.testing.assert_array_equal(draws[:,1,[0,2]],2*draws[:,0,[0,2]])
        again,h2=a.bootstrap_totals(values,0,200);np.testing.assert_array_equal(draws,again);self.assertEqual(h,h2)
        other,h3=a.bootstrap_totals(values,1,200);self.assertNotEqual(h,h3)

    def test_policy_report_keeps_three_contrasts_on_selection_failure(self):
        values=np.tile(np.array([[[5.,1680.,2.,1680.]]]),(512,1,1))
        report=a.summarize_values(values,dict(scenarios=[dict(id='independent',h1=.95)]))
        self.assertEqual(len(report['primary_contrasts']),3)
        self.assertTrue(all(c['estimate'] is None for c in report['primary_contrasts']))
        self.assertAlmostEqual(report['scenarios'][0]['F'],.6)
        np.testing.assert_allclose(report['scenarios'][0]['F95']['interval'],[.6,.6])
        self.assertEqual(report['primary_contrasts'][0]['uncertainty']['level'],1-.05/3)

    def test_full_frozen_selector_preserves_failed_profile_and_conditional_uncertainty(self):
        cfg=a.read(a.HERE/'configs/grid.json');book=a.read(a.HERE/'seedbook.json')
        values=np.empty((128,12,4))
        rates=[1.,2.,4.,8.]*2+[0.,0.,1.,2.]
        for j,rate in enumerate(rates):values[:,j,:]=[rate,1000.,rate,1000.]
        report=a.select_values(values,cfg,book['policy'])
        self.assertEqual([s['id'] for s in report['policy_config']['scenarios']],['independent','clustered'])
        bad=report['profiles']['imperfect']['selection_uncertainty']
        self.assertEqual(bad['failed_draws'],a.DRAWS);self.assertIsNone(bad['conditional_on_success95'])
        good=report['profiles']['independent']['selection_uncertainty']
        self.assertEqual(good['failed_draws'],0);self.assertTrue(good['conditional_interval'])
        self.assertAlmostEqual(report['profiles']['independent']['selected_p'],
                               cfg['scenarios'][1]['p']*1.5)
        wrong=list(book['policy']);wrong[0]=next(s for s in range(1000) if s not in book['policy']+book['grid'])
        with self.assertRaisesRegex(ValueError,'frozen policy seeds'):a.select_values(values,cfg,wrong)

    def test_immutable_outputs(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td).resolve()/'receipt.json';a.write_new(path,{'a':1});a.write_new(path,{'a':1})
            with self.assertRaisesRegex(ValueError,'overwrite'):a.write_new(path,{'a':2})

    def test_streamed_manifest_binding_and_record_inventory(self):
        row,cfg=fixture()
        paths=set(v.SOURCE_PATHS)|{f'research/{a.HERE.name}/workflow_variant.py',f'research/{a.HERE.name}/worker.py'}
        sources={name:a.digest(a.ABM/name) for name in paths}
        manifest=dict(schema_version=1,study=a.HERE.name,model_factory='WorkflowRobustnessModel',task_count=1,
            config=cfg,resolved_config=dict(cfg,warmup_ticks=233280,observation_ticks=145152,block_ticks=2592),
            jobs=a.jobs_for(cfg),software=SOFTWARE,sources=sources,freeze_identity='selected',parent_freeze_identity='parent')
        manifest['identity']=a.identity(manifest);row['identity']=manifest['identity']
        with tempfile.TemporaryDirectory() as td:
            root=Path(td).resolve();(root/'completed').mkdir()
            (root/'manifest.json').write_text(json.dumps(manifest))
            record_path=root/'completed/job_00000000.json';record_path.write_text(json.dumps(row))
            values,receipt=a.load_bank(root,cfg,sources,'selected','parent')
            self.assertEqual(values.tolist(),[[[56.,1680.,0.,1680.]]])
            self.assertEqual(receipt['complete_input_hashes'],{'completed/job_00000000.json':a.digest(record_path)})
            with self.assertRaisesRegex(ValueError,'execution freeze'):a.load_bank(root,cfg,sources,'other','parent')
            (root/'completed/extra.json').write_text('{}')
            with self.assertRaisesRegex(ValueError,'inventory'):a.load_bank(root,cfg,sources,'selected','parent')

    def test_manifest_receipt_binds_bytes_read_before_record_streaming(self):
        row,cfg=fixture()
        paths=set(v.SOURCE_PATHS)|{f'research/{a.HERE.name}/workflow_variant.py',f'research/{a.HERE.name}/worker.py'}
        sources={name:a.digest(a.ABM/name) for name in paths}
        manifest=dict(schema_version=1,study=a.HERE.name,model_factory='WorkflowRobustnessModel',task_count=1,
            config=cfg,resolved_config=dict(cfg,warmup_ticks=233280,observation_ticks=145152,block_ticks=2592),
            jobs=a.jobs_for(cfg),software=SOFTWARE,sources=sources,freeze_identity='selected',parent_freeze_identity='parent')
        manifest['identity']=a.identity(manifest);row['identity']=manifest['identity']
        with tempfile.TemporaryDirectory() as td:
            root=Path(td).resolve();(root/'completed').mkdir()
            manifest_path=root/'manifest.json';manifest_path.write_text(json.dumps(manifest))
            original_hash=a.digest(manifest_path)
            record_path=root/'completed/job_00000000.json';record_path.write_text(json.dumps(row))
            original_bound_read=v.bound_read

            def swap_during_record_read(path):
                result=original_bound_read(path)
                if path==record_path:
                    replacement=dict(manifest,task_count=2)
                    replacement['identity']=a.identity(replacement)
                    manifest_path.write_text(json.dumps(replacement))
                return result

            with patch.object(v,'bound_read',side_effect=swap_during_record_read):
                _,receipt=a.load_bank(root,cfg,sources,'selected','parent')
            # The receipt must bind the manifest actually interpreted, enabling
            # main's final input-hash recheck to reject the changed file.
            self.assertEqual(receipt['manifest_identity'],manifest['identity'])
            self.assertEqual(receipt['manifest_sha256'],original_hash)
            self.assertNotEqual(receipt['manifest_sha256'],a.digest(manifest_path))


def math_sqrt3_over10():
    return float(np.sqrt(3)/10)


if __name__=='__main__':unittest.main()
