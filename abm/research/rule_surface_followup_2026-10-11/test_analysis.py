"""Prospective artificial accounting and statistical fixtures; no new outcomes."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock
import numpy as np
import analyze as a
import select_clock as s

SOFTWARE=dict(a.PACKAGES,python='3.12.3 (artificial fixture build)')


def ledger_start():
    row={k:(0. if k in a.EXPOSURES else 0) for k in a.COUNTERS}
    row.update(total_ticks=0,census=30,colonized_census=2,cumul_patients=30,imported_colonizations=2,
        environment_contaminated_site_days=0.,environment_contaminated_sites_endpoint=1,
        hygiene_counts_by_role={r:dict.fromkeys(a.HYGIENE,0) for r in ('Nurse','Dr','XrayDr')},
        response_counts={k:{} for k in a.RESPONSE_KINDS})
    return row


def ledger_period(start,n,acquired):
    end=copy.deepcopy(start);blocks=[]
    for day in range(n):
        response={k:{} for k in a.RESPONSE_KINDS}
        response.update(pickups={'patient':2},pickups_by_completion_state={'completed':2},
            acquisitions_by_patient_source_status={'imported':acquired},acquisitions_by_route={'patient_pickup':acquired},
            acquisitions_by_opportunities={'0':acquired},acquisitions_by_route_and_opportunities={'patient_pickup|0':acquired},
            hand_at_risk_opportunities={'1':2},hand_cleared_opportunities={'1':2},
            environment_events={'scheduled_clear':1,'dirty_scheduled_clear':1})
        b=dict.fromkeys(a.COUNTERS,0)
        roles={r:dict.fromkeys(a.HYGIENE,0) for r in ('Nurse','Dr','XrayDr')}
        roles['Nurse']=dict(zip(a.HYGIENE,(40,36,2,2)))
        b.update(acquired_colonizations=acquired,cumul_patients=1,discharged_patients=1,
            bed_visits=40,patient_contacts=40,environment_contacts=20,patient_days=30.,susceptible_patient_days=28.,
            hygiene_opportunity_count=40,hygiene_completion_count=36,hygiene_effective_clear_count=2,
            hygiene_decontamination_count=2,environment_contaminated_site_days=1.,ticks=2592,days=1.,
            day_start=float(day),day_end=float(day+1),hygiene_counts_by_role=roles,response_counts=response,
            mean_census=30.,mean_susceptible_census=28.,colonized_patient_days=2.,mean_colonized_census=2.,
            mean_environment_contaminated_fraction=1/30,acquired_per_1000_patient_days=1000*acquired/30,
            acquired_per_1000_susceptible_patient_days=1000*acquired/28)
        for key in ('census','colonized_census','environment_contaminated_sites_endpoint'):
            b[key+'_start']=b[key+'_end']=start[key]
        for key in a.COUNTERS+('environment_contaminated_site_days',):end[key]+=b[key]
        end['total_ticks']+=2592
        for group in ('hygiene_counts_by_role','response_counts'):
            for category,counts in b[group].items():
                for key,value in counts.items():end[group][category][key]=end[group][category].get(key,0)+value
        blocks.append(b)
    total={key:end[key]-start[key] for key in a.COUNTERS+('environment_contaminated_site_days',)}
    total.update(ticks=2592*n,days=float(n),mean_census=30.,mean_susceptible_census=28.,
        colonized_patient_days=2.*n,mean_colonized_census=2.,mean_environment_contaminated_fraction=1/30,
        acquired_per_1000_patient_days=1000*acquired/30,acquired_per_1000_susceptible_patient_days=1000*acquired/28,
        hygiene_counts_by_role={r:{k:end['hygiene_counts_by_role'][r][k]-start['hygiene_counts_by_role'][r][k]
            for k in a.HYGIENE} for r in ('Nurse','Dr','XrayDr')},
        response_counts={kind:{k:end['response_counts'][kind][k]-start['response_counts'][kind].get(k,0)
            for k in end['response_counts'][kind]} for kind in a.RESPONSE_KINDS})
    for key in ('census','colonized_census','environment_contaminated_sites_endpoint'):
        total[key+'_start']=start[key];total[key+'_end']=end[key]
    return total,blocks,end


def record(shared=False):
    sc=dict(id='toy',p=.2,h0=.9,h1=.9 if shared else .95,rho=.5,e=1.,
        opportunities='separate',environment_multiplier=1.)
    job=dict(job_id=0,scenario=sc,seed=123)
    initial=ledger_start();warm,blocks,boundary=ledger_period(initial,90,1);arms={}
    for name,h,count in (('baseline',sc['h0'],1),('policy',sc['h1'],0)):
        obs,daily,end=ledger_period(boundary,56,count)
        arms[name]=dict(h=h,observation=obs,observation_blocks=daily,final_snapshot=end)
    if shared:arms['policy']={'shared_with':'baseline'}
    row=dict(job,status='complete',identity='fixture',execution=dict(SOFTWARE,platform='artificial'),
        distinct_observation_trajectories=1 if shared else 2,warmup_executions=1,executed_ticks_per_trajectory=378432,
        model_kwargs=dict(a.MODEL,seed=123,prob_transmission=.2,icu_hcw_wash_rate=.9,outside_hcw_wash_rate=.9,
            hygiene_opportunities='separate',hygiene_autocorrelation=.5,hygiene_efficacy=1.,environment_multiplier=1.),
        initial_snapshot=initial,warmup_snapshot=boundary,warmup=warm,warmup_blocks=blocks,arms=arms)
    return row,job


class AnalysisTests(unittest.TestCase):
    def rejected(self,mutate,pattern='.'):
        row,job=record();mutate(row)
        with self.assertRaisesRegex(ValueError,pattern):a.validate_record(row,job,SOFTWARE)

    def test_accounting_pair_and_shared_baseline_have_complete_halves(self):
        row,job=record();values,windows=a.validate_record(row,job,SOFTWARE)
        self.assertEqual(values,[56,1680.,0,1680.]);self.assertEqual(windows,[28,840.,28,840.,0,840.,0,840.])
        row,job=record(True);values,windows=a.validate_record(row,job,SOFTWARE)
        self.assertEqual(values,[56,1680.,56,1680.]);self.assertEqual(windows[:4],windows[4:])

    def test_nonbaseline_shared_branch_is_rejected(self):
        self.rejected(lambda r:r['arms'].update(policy={'shared_with':'baseline'}))
        row,job=record(True);row['arms']['policy']={'shared_with':'policy'}
        with self.assertRaises(ValueError):a.validate_record(row,job,SOFTWARE)

    def test_runtime_and_model_declarations_are_exact(self):
        self.rejected(lambda r:r['execution'].update(python='3.12.3 (other build)'),'runtime')
        self.rejected(lambda r:r['model_kwargs'].update(prob_transmission=.2+1e-12),'declaration')
        self.rejected(lambda r:r.update(distinct_observation_trajectories=True),'Execution')

    def test_period_and_integer_endpoint_schema_reject_mutations(self):
        self.rejected(lambda r:r['warmup_blocks'].pop(),'period')
        self.rejected(lambda r:r['warmup'].update(census_start=30.,census_end=30.),'Period snapshot')
        self.rejected(lambda r:r['arms']['policy']['final_snapshot'].update(total_ticks=378433),'horizon')

    def test_counter_and_exposure_conservation(self):
        self.rejected(lambda r:r['warmup_blocks'][0].update(bed_visits=41),'Integer')
        self.rejected(lambda r:r['warmup_blocks'][0].update(patient_days=29.),'Exposure')
        self.rejected(lambda r:r['initial_snapshot'].update(cumul_patients=31),'Snapshot patient')

    def test_role_histogram_reconciliation_is_not_just_temporal(self):
        self.rejected(lambda r:r['warmup_blocks'][0]['hygiene_counts_by_role']['Nurse'].update(hygiene_completion_count=35),'reconciliation')
        self.rejected(lambda r:r['warmup_blocks'][0]['response_counts']['acquisitions_by_route'].update(patient_pickup=2),'Acquisition')
        self.rejected(lambda r:r['warmup_blocks'][0]['response_counts']['hand_survived_opportunities'].update({'1':1}),'Hand-at-risk')
        self.rejected(lambda r:r['initial_snapshot']['hygiene_counts_by_role'].update(Other={}), 'worker roles')

    def test_zero_control_event_fraction_is_undefined(self):
        point=a.estimates(np.array([0.,10.,1.,20.]),2)
        self.assertEqual(point[0],0.);self.assertEqual(point[1],50.);self.assertTrue(np.isnan(point[2]))
        self.assertEqual(point[4],-.5)
        self.assertIsNone(a.interval(np.array([.2,np.nan]),.95)['interval'])
        self.assertEqual(a.interval(np.array([.2,np.nan]),.95)['undefined_draws'],1)
        self.assertIsNone(a.difference_or_none(None,.3))

    def test_pooled_unequal_exposure_fraction_and_paired_count(self):
        values=np.array([[2.,10.,1.,20.],[3.,40.,1.,10.]])
        point=a.estimates(values.sum(axis=0),2)
        self.assertAlmostEqual(point[0],100.);self.assertAlmostEqual(point[1],200/3)
        self.assertAlmostEqual(point[2],1/3);self.assertAlmostEqual(point[3],100/3);self.assertEqual(point[4],1.5)
        per_seed=a.estimates(values,1)[:,2].mean();self.assertNotAlmostEqual(point[2],per_seed)

    def test_four_strata_fixed_core_manual_oracle_and_weights(self):
        values=np.zeros((8,2,4));windows=np.zeros((8,2,8));core={0,1,4,5}
        for i in range(8):
            v=np.array([i+1.,10.+i,i%3+1.,20.+2*i]);values[i,0]=v
            if i in core:values[i,1]=v
        for j in range(2):
            windows[:,j,:4]=values[:,j]/2;windows[:,j,4:]=values[:,j]/2
            # half columns require baseline (count,exposure,count,exposure), then policy.
            windows[:,j]=np.stack((values[:,j,0]/3,values[:,j,1]/2,2*values[:,j,0]/3,values[:,j,1]/2,
                values[:,j,2]/4,values[:,j,3]/2,3*values[:,j,2]/4,values[:,j,3]/2),axis=-1)
        groups=[[0,1],[2,3],[4,5],[6,7]];draws=1003;seed=79;n=[8,4]
        actual,actualwin,sha=a.boot(values,windows,groups,np.array(n),np.random.default_rng(seed),draws)
        rng=np.random.default_rng(seed);stream=bytearray();oracle=[];winoracle=[]
        for begin in range(0,draws,500):
            size=min(500,draws-begin);weights=[[0]*8 for _ in range(size)]
            for group in groups:
                bank=rng.multinomial(2,[.5,.5],size=size)
                for b in range(size):
                    for ii,i in enumerate(group):weights[b][i]=int(bank[b,ii])
            for weight in weights:
                for w in weight:stream.extend(w.to_bytes(8,'little',signed=True))
                self.assertEqual(sum(weight[i] for i in core),4)
                point=[];win=[]
                for c in range(2):
                    cn,cd,pn,pd=[sum(weight[i]*values[i,c,k] for i in range(8)) for k in range(4)]
                    r0=1000*cn/cd;r1=1000*pn/pd
                    point.append([r0,r1,1-r1/r0,r0-r1,(cn-pn)/n[c]])
                    ww=[sum(weight[i]*windows[i,c,k] for i in range(8)) for k in range(8)]
                    win.append([1000*(ww[2]/ww[3]-ww[0]/ww[1]),1000*(ww[6]/ww[7]-ww[4]/ww[5])])
                oracle.append(point);winoracle.append(win)
        np.testing.assert_allclose(actual,oracle,rtol=1e-14,atol=1e-14)
        np.testing.assert_allclose(actualwin,winoracle,rtol=1e-14,atol=1e-14)
        self.assertEqual(sha,hashlib.sha256(stream).hexdigest())

    def test_bootstrap_rejects_duplicated_or_missing_strata(self):
        values=np.ones((4,1,4));windows=np.ones((4,1,8))
        for groups in ([[0,1],[1,2,3]],[[0,1],[3]]):
            with self.assertRaises(ValueError):a.boot(values,windows,groups,4,np.random.default_rng(3),2)

    def test_incomplete_inventory_precedes_outcome_read_and_symlink_refusal(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp).resolve();(root/'completed').mkdir()
            p=root/'completed/job_00000000.json';p.write_text('not valid JSON')
            with mock.patch.object(a,'bound_read',side_effect=AssertionError('Outcome read')):
                with self.assertRaisesRegex(ValueError,'outcomes were not read'):a.load_bank(root,{'seeds':[1,2],'scenarios':[{}]},'x','fixture')
            link=root/'linked.json';link.symlink_to(p)
            with self.assertRaisesRegex(ValueError,'Symlink'):a.read(link)

    def test_exact_hit_takes_priority_and_lowest_hit_wins(self):
        chosen=s.choose([(0.,1.),(.1,4.),(.2,3.),(.3,3.)])
        self.assertEqual(chosen['method'],'lowest_exact_hit');self.assertEqual(chosen['selected_p'],.2)

    def test_lowest_increasing_adjacent_bracket_and_failed_decrease(self):
        chosen=s.choose([(0.,2.),(.1,4.),(.2,2.),(.3,4.)])
        self.assertAlmostEqual(chosen['selected_p'],.05);self.assertEqual(len(chosen['increasing_brackets']),2)
        self.assertEqual(len(chosen['nonincreasing_adjacent_segments']),1)
        self.assertEqual(s.choose([(0.,4.),(.2,2.)])['status'],'failed')
        self.assertEqual(s.choose([(0.,1.),(.2,2.)])['status'],'failed')

    def test_selection_refuses_boolean_nonfinite_duplicate_unordered_points(self):
        for points in ([(True,2.)],[(.2,float('nan'))],[(.1,2.),(.1,4.)],[(.2,2.),(.1,4.)]):
            with self.assertRaises(ValueError):s.choose(points)

    def test_selection_failure_retained_in_unconditional_interval(self):
        values=np.array([[[0,1000,0,1000],[4,1000,4,1000]],[[0,1000,0,1000],[0,1000,0,1000]]],dtype=float)
        sc=[dict(id='clustered_a',p=.1),dict(id='clustered_b',p=.2)]
        # Supply same two-point profile names to test both specified selectors.
        values=np.concatenate((values,values),axis=1);sc+= [dict(id='imperfect_a',p=.1),dict(id='imperfect_b',p=.2)]
        result,_=s.selection_draws(values,sc,[[0,1]],draws=101,seed=7)
        for row in result.values():
            self.assertGreater(row['failed_draws'],0);self.assertLess(row['failed_draws'],101)
            self.assertIsNone(row['unconditional95']['interval']);self.assertIsNotNone(row['conditional_on_valid_selection95'])

    def test_write_is_idempotent_and_never_reselects(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp).resolve()/'selection.json';s.write_new(path,{'p':.2});original=path.read_bytes()
            s.write_new(path,{'p':.2});self.assertEqual(original,path.read_bytes())
            with self.assertRaisesRegex(ValueError,'no reselection'):s.write_new(path,{'p':.3})

    def test_policy_config_preserves_fixed_seeds_and_only_selected_profiles(self):
        cfg=dict(warmup_days=90,observation_days=56,block_days=1,model_overrides={})
        book=dict(clock_policy=dict(A=list(range(256)),B=list(range(256,512))))
        result=s.policy_config(cfg,book,dict(clustered={'selected_p':.23},imperfect={'selected_p':None}))
        self.assertEqual(result['seeds'],list(range(512)));self.assertEqual(len(result['scenarios']),1)
        self.assertEqual(result['scenarios'][0],dict(id='clustered_selected_clock',p=.23,rho=.5,e=1.,h0=.9,h1=.95,
            opportunities='separate',environment_multiplier=1.))

    def test_resealed_selection_and_config_must_recompute_from_baseline(self):
        cfg=a.read(a.HERE/'configs/clock_grid.json');book=a.read(a.HERE/'seedbook.json')
        values=np.zeros((128,19,4));profiles={}
        for i in range(19):
            count=1.+.45*i if i<9 else 6.+i
            values[:,i]=[count,1000.,count,1000.]
        totals=values.sum(axis=0);rates=1000*totals[:,0]/totals[:,1]
        for profile in ('clustered','imperfect'):
            ii=[i for i,sc in enumerate(cfg['scenarios']) if sc['id'].startswith(profile+'_')]
            profiles[profile]=s.choose([(cfg['scenarios'][i]['p'],float(rates[i])) for i in ii])
            profiles[profile]['grid_points']=[dict(scenario_id=cfg['scenarios'][i]['id'],p=cfg['scenarios'][i]['p'],
                acquired_total=int(totals[i,0]),patient_days_total=float(totals[i,1]),rate=float(rates[i])) for i in ii]
            profiles[profile]['banks']={}
            for label,rows in (('A',range(64)),('B',range(64,128))):
                tt=values[list(rows)].sum(axis=0);rr=1000*tt[:,0]/tt[:,1]
                profiles[profile]['banks'][label]=s.choose([(cfg['scenarios'][i]['p'],float(rr[i])) for i in ii])
        grid=dict(config=cfg,values=values,manifest_sha256='m',binding_sha256='b',input_sha256={'completed/toy':'h'})
        policy=s.policy_config(cfg,book,profiles)
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp).resolve();(root/'configs').mkdir();(root/'select_clock.py').write_bytes((s.HERE/'select_clock.py').read_bytes())
            original=dict(schema_version=1,study=root.name,status='complete_baseline_selection',parent_freeze_identity='parent',
                source_sha256=a.digest(root/'select_clock.py'),policy_outcomes_used=False,bootstrap_draws=50000,bootstrap_seed=20261011,
                seed_banks=dict(grid=book['clock_grid'],policy=book['clock_policy']),profiles=profiles,selected_profiles=['clustered'],
                grid_manifest_sha256='m',grid_binding_sha256='b',grid_input_hashes=grid['input_sha256'])
            def check(selection,derived=policy,extra_file=False):
                (root/'selection.json').write_text(json.dumps(selection));(root/'configs/clock_policy.json').write_text(json.dumps(derived))
                files={f'research/{root.name}/selection.json':a.digest(root/'selection.json'),
                    f'research/{root.name}/configs/clock_policy.json':a.digest(root/'configs/clock_policy.json')}
                if extra_file:files['extra']='h'
                secondary=dict(schema_version=1,study=root.name,parent_freeze_identity='parent',files=files,
                    grid_manifest_sha256='m',grid_input_hashes=grid['input_sha256'])
                secondary['identity_sha256']=a.identity(secondary)
                with mock.patch.object(a,'HERE',root):a.verify_selection(selection,secondary,grid,book,'parent',derived)
            check(original)
            for mutate in (lambda r:r['profiles']['clustered'].update(selected_p=.3),
                           lambda r:r.update(policy_outcomes_used=True),
                           lambda r:r.update(selected_profiles=['imperfect','clustered']),
                           lambda r:r.update(grid_binding_sha256='other'),
                           lambda r:r.update(schema_version=True),
                           lambda r:r['profiles']['clustered']['grid_points'][0].update(acquired_total=9)):
                changed=copy.deepcopy(original);mutate(changed)
                with self.assertRaises(ValueError):check(changed)
            wrong=copy.deepcopy(policy);wrong['scenarios'][0]['p']+=.001
            with self.assertRaises(ValueError):check(original,wrong)
            with self.assertRaises(ValueError):check(original,extra_file=True)

    def test_parse_rejects_duplicate_or_nonfinite_and_exact_bool(self):
        for raw in ('{"x":1,"x":2}','{"x":NaN}'):
            with self.assertRaises(ValueError):a.parse(raw)
        with self.assertRaises(ValueError):a.exact(True,1)


if __name__=='__main__':unittest.main()
