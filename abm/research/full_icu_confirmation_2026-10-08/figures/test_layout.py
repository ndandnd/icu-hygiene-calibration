#!/usr/bin/env python3
"""Bounded tests and visibly labeled synthetic layout PNG/TeX fixtures only."""
import argparse
import copy
import json
import math
from pathlib import Path
import tempfile
import unittest
import build_figures as b


def uncertainty(value, width=.05, undefined=False):
    def ci(factor):return None if undefined else [value-factor*width,value+factor*width]
    return dict(interval95=ci(1),intervals={format(level,'.12g'):ci(factor)
        for level,factor in ((.95,1),(b.PRIMARY,1.3),(b.GATE,1.5))},bootstrap_draws=200000,
        undefined_bootstrap_draws=int(undefined),interval_unavailable_if_any_draw_undefined=True,conditional_on_defined_draws=False)


def est(value,width=.05):return dict(estimate=value,uncertainty=uncertainty(value,width))


def row(identifier,p,rho=0.,e=1.,environment=1.,f=.25):
    return dict(scenario=dict(id=identifier,p=p,rho=rho,e=e,h0=.9,h1=.95,environment_multiplier=environment,opportunities='zone_separate'),
        n=512,seeds=list(range(1000,1512)),arms={name:dict(rate_per_1000_patient_days=value,
            rate_uncertainty=uncertainty(value,.45)) for name,value in (('baseline',3.2),('policy',2.1))},
        rate_prevented_fraction=f,rate_prevented_fraction_uncertainty=uncertainty(f,.14),
        mean_prevented=1.3,mean_prevented_uncertainty=uncertainty(1.3,1.6),
        absolute_prevented_rate=.7,absolute_prevented_rate_uncertainty=uncertainty(.7,.9))


def layouts():
    hr={p+'_target_candidate':row(p+'_target_candidate',prob,*b.PARAMS[p],f=f)
        for p,prob,f in zip(b.PROFILES,(.15,.07,.06),(.2,.3,.4))}
    gates=[]
    for p in b.PROFILES:gates.append(dict(name=p+'_absolute_baseline',passed=False,band=[2.7,3.3],**est(3.2,.45)))
    for a,c in (('independent','persistent'),('independent','imperfect'),('persistent','imperfect')):
        gates.append(dict(name=a+'_over_'+c+'_baseline_ratio',passed=False,band=[1/1.1,1.1],**est(1.,.16)))
    holdout=dict(is_synthetic_fixture=True,baseline_matching=dict(numerical_gates=gates,
        numerical_interval_level=b.GATE,valid_calibration=True,established=False))
    er={f'independent_low{label}_env{env}':row(f'independent_low{label}_env{env}',prob,environment=env,f=.2+.1*i+.25*(1-env))
        for i,(label,prob) in enumerate(zip(('075','100','125'),(.08,.11,.135))) for env in (1,0)}
    environment=dict(is_synthetic_fixture=True,primary_environment_effects=[dict(p=prob,
        full_minus_disabled={name:est(value,width) for name,value,width in (('F',-.25,.16),('mean_prevented',.5,.4),('prevented_rate',.2,.3))})
        for prob in (.08,.11,.135)])
    bridge=dict(is_synthetic_fixture=True,profiles={b.BRIDGE_KEYS[p]:dict(limiting_response=f,
        limiting_response_normal_95_interval=[f-.1,f+.1],stationary_q2_response_for_comparison_only=q,
        K0=dict(mean=100),K1=dict(mean=100*(1-f)))
        for p,f,q in zip(b.PROFILES,(.15,.5,.28),(.75,.65,.35))})
    bridge['observation_patient_days']=dict(mean=1600)
    return holdout,hr,environment,er,bridge


class Guards(unittest.TestCase):
    def test_all_zero_nonpercent_axes_keep_zero_target_and_band(self):
        class Recorder:
            def __init__(self):self.commands=[]
            def line(self,*args,**kwargs):self.commands.append(('line',args))
            def rectangle(self,*args,**kwargs):self.commands.append(('rectangle',args))
            def text(self,*args,**kwargs):self.commands.append(('text',args))
        for refs,band,nonnegative in (((0,),None,False),((0,),None,True),
            ((0,3),None,True),((0,3),(2.7,3.3),True)):
            with self.subTest(refs=refs,band=band,nonnegative=nonnegative):
                scene=Recorder()
                xx=b.axis(scene,10,210,50,100,[(0.,[0.,0.])],'Zero estimand',refs,band=band,nonnegative=nonnegative)
                for value in refs+(tuple(band) if band else ()):
                    self.assertTrue(math.isfinite(xx(value)))
                    self.assertGreaterEqual(xx(value),10)
                    self.assertLessEqual(xx(value),210)
                self.assertTrue(scene.commands)
                self.assertTrue(all(math.isfinite(a) for _,args in scene.commands for a in args if type(a) in (int,float)))

    def test_nine_point_only_leading_rate_formulas(self):
        h,hr,e,er,k=layouts()
        rows=b.leading_rate_rows(hr,er,k)
        self.assertEqual(len(rows),9)
        self.assertAlmostEqual(rows[0]['R0_leading'],1000*.15**2*100/1600)
        self.assertAlmostEqual(rows[0]['R1_leading'],1000*.15**2*85/1600)
        for i in (3,5,7):
            self.assertEqual(rows[i]['R0_leading'],rows[i+1]['R0_leading'])
            self.assertEqual(rows[i]['R1_leading'],rows[i+1]['R1_leading'])
        self.assertEqual(rows[0]['F_measured'],hr['independent_target_candidate']['rate_prevented_fraction'])
        self.assertEqual(rows[0]['limiting_F_normal_CI95_low'],k['profiles']['iid_perfect']['limiting_response_normal_95_interval'][0])
    def test_undefined_interval_never_filtered(self):
        u=uncertainty(.25,undefined=True)
        self.assertEqual(b.metric(.25,u),(.25,None))
        u['intervals']['0.95']=[.1,.4]
        with self.assertRaisesRegex(ValueError,'undefined draw'):b.metric(.25,u)

    def test_negative_fraction_axis_not_clipped(self):
        lo,hi=b.bounds([(-.2,[-.7,.1])],(0,))
        self.assertLess(lo,-.7);self.assertGreater(hi,.1)

    def test_reject_partial_wrong_n_and_synthetic(self):
        # No actual outcome files are read. The valid analyzer/source freeze
        # and an ephemeral artificial report exercise the production guards.
        with tempfile.TemporaryDirectory(dir=b.HERE) as tmp:
            tmp=Path(tmp)
            h,rows,*_=layouts()
            design=dict(warmup_days=90,observation_days=56,model_overrides={'initial_environment_contaminated':False,'prob_patient_sick':.03},
                seeds=list(range(1000,1512)),scenarios=[r['scenario'] for r in rows.values()])
            dp=tmp/'design.json';dp.write_text(json.dumps(design))
            analyzer='research/full_icu_confirmation_2026-10-08/analyze.py'
            source='src/agents.py'
            files={analyzer:b.digest(b.ABM/analyzer),source:b.digest(b.ABM/source),
                str(dp.relative_to(b.ABM)):b.digest(dp),
                'research/full_icu_confirmation_2026-10-08/seedbook.json':b.digest(b.NEW/'seedbook.json'),
                'research/full_icu_confirmation_2026-10-08/candidate_design.json':b.digest(b.NEW/'candidate_design.json')}
            freeze=dict(files=files)
            import hashlib
            freeze['identity_sha256']=hashlib.sha256(b.canonical(freeze).encode()).hexdigest()
            fp=tmp/'freeze.json';fp.write_text(json.dumps(freeze))
            report=dict(h, schema_version=1,status='complete',module='holdout',n_paired_jobs=1536,bootstrap_draws=200000,
                scenarios=list(rows.values()),design_sha256=b.digest(dp),freeze_sha256=b.digest(fp),freeze_identity=freeze['identity_sha256'],
                analysis_freeze_sha256=b.digest(fp),analysis_freeze_identity=freeze['identity_sha256'],source_sha256=files[analyzer],
                seedbook_sha256=files['research/full_icu_confirmation_2026-10-08/seedbook.json'],
                proposal_sha256=files['research/full_icu_confirmation_2026-10-08/candidate_design.json'],
                calibration_selection_sha256=files[source],calibration_report_sha256=files[source],
                primary_family_interval_level=b.PRIMARY,banks=[dict(jobs=1536,manifest_sha256='0'*64,
                    sources={source:files[source]},completed_records_sha256={f'synthetic_{i}.json':'0'*64 for i in range(1536)})])
            rp=tmp/'report.json'
            def validate(value):
                rp.write_text(json.dumps(value))
                return b.validate_report(rp,dp,fp,fp,'holdout')
            with self.assertRaisesRegex(ValueError,'production report'):validate(report)
            # This copy stays ephemeral and is never rendered or published.
            report['is_synthetic_fixture']=False
            self.assertEqual(len(validate(report)[1]),3)
            for mutation,pattern in (({'status':'partial'},'complete fresh'),({'n_paired_jobs':1535},'Exact planned'),
                ({'design_sha256':'0'*64},'binding mismatch')):
                with self.assertRaisesRegex(ValueError,pattern):validate(dict(report,**mutation))
            missing=copy.deepcopy(report);missing['scenarios'][0]['seeds']=missing['scenarios'][0]['seeds'][:-1]
            with self.assertRaisesRegex(ValueError,'ordered seed'):validate(missing)
            duplicates=copy.deepcopy(report);duplicates['scenarios'][1]=duplicates['scenarios'][0]
            with self.assertRaisesRegex(ValueError,'duplicate'):validate(duplicates)

    def test_native_picture_newlines_comment_closed(self):
        with tempfile.TemporaryDirectory(dir=b.HERE) as tmp:
            s=b.load_scene()(40,30);s.text(5,10,'Test',size=9)
            b.save_scene(s,Path(tmp),'native',False)
            text=(Path(tmp)/'native.tex').read_text()
            interior=text.split('\\resizebox{\\linewidth}{!}{%\n',1)[1].rsplit('}\n',1)[0]
            self.assertTrue(all(line.endswith('%') for line in interior.splitlines()))

    def test_bridge_rejects_default_horizon_wrong_n_and_missing_seeds(self):
        import hashlib
        with tempfile.TemporaryDirectory(dir=b.HERE) as tmp:
            tmp=Path(tmp)
            config=b.read(b.NEW/'bridge/configs/bridge_all.json')
            cp=tmp/'config.json';cp.write_text(json.dumps(config))
            sources={f'src/{n}.py':b.digest(b.ABM/f'src/{n}.py') for n in ('CPEmodel','agents','routes')}
            sources.update({f'research/full_icu_confirmation_2026-10-08/bridge/{n}.py':b.digest(b.NEW/f'bridge/{n}.py') for n in ('run','coefficient','background')})
            freeze=dict(files=sources);freeze['identity_sha256']=hashlib.sha256(b.canonical(freeze).encode()).hexdigest()
            fp=tmp/'freeze.json';fp.write_text(json.dumps(freeze))
            *_,r=layouts();r.update(module='bridge',status='complete',n_backgrounds=256,seeds_completed=256,seeds_planned=256,is_synthetic_fixture=False)
            r.update(scientific_sources=sources,source_sha256=sources,analysis_code_sha256='1'*64)
            for row in r['profiles'].values():
                for name in ('K0','K1'):row[name]['normal_approximation_95_interval']=[80,120]
            rp=tmp/'summary.json';bp=tmp/'binding.json'
            def validate(summary=r,conf=config,seed_list=None):
                cp.write_text(json.dumps(conf))
                summary=dict(summary,config_sha256=b.digest(cp),seeds=config['seeds'],freeze_sha256=b.digest(fp),freeze_identity=freeze['identity_sha256'])
                rp.write_text(json.dumps(summary))
                binding=dict(status='complete',summary_sha256=b.digest(rp),config_sha256=b.digest(cp),
                    seeds=config['seeds'] if seed_list is None else seed_list,freeze_sha256=b.digest(fp),
                    freeze_identity=freeze['identity_sha256'],scientific_sources=sources,source_sha256=sources,analysis_code_sha256='1'*64)
                bp.write_text(json.dumps(binding))
                return b.validate_bridge(rp,cp,bp,fp)
            self.assertEqual(validate()['seeds_completed'],256)
            with self.assertRaisesRegex(ValueError,'256-background'):validate(dict(r,seeds_completed=255))
            with self.assertRaisesRegex(ValueError,'90\\+56-day'):validate(conf=dict(config,observation_ticks=28*2592))
            with self.assertRaisesRegex(ValueError,'exact seeds'):validate(seed_list=config['seeds'][:-1])


def render_fixture(output):
    if output.exists():raise FileExistsError('Use a fresh fixture output directory')
    if not output.resolve().is_relative_to(b.HERE):raise ValueError('Fixture must stay in new figures directory')
    output.mkdir(parents=True)
    h,hr,e,er,k=layouts();Scene=b.load_scene()
    for name,scene in (('holdout_rates_and_prevention',b.draw_holdout(Scene,h,hr,True)),
        ('environment_estimands',b.draw_environment(Scene,e,er,True)),
        ('bridge_limit_comparison',b.draw_bridge(Scene,h,hr,k,True))):
        b.save_scene(scene,output,name,False)
    manifest=dict(status='synthetic_layout_fixture_not_for_paper',is_synthetic_fixture=True,
        warning='All values are invented layout inputs. No confirmation outcome files read. Do not use in paper.',
        builder_sha256=b.digest(Path(b.__file__)),test_source_sha256=b.digest(Path(__file__)),
        outputs_sha256={p.name:b.digest(p) for p in sorted(output.iterdir())})
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--render-fixture',type=Path)
    args,rest=p.parse_known_args()
    if args.render_fixture:render_fixture(args.render_fixture)
    else:unittest.main(argv=[__file__]+rest)
