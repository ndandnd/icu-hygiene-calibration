#!/usr/bin/env python3
"""Draw complete confirmation reports; no simulation, interpolation or analysis.

All production inputs are explicit. See DESIGN.md for the bridge aggregation
binding contract. Test layouts use the separate test_layout.py entry point.
"""
import argparse
import csv
from datetime import datetime,timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
NEW = HERE.parent
ABM = HERE.parents[2]
RENDERER = ABM/'research/full_icu_calibration_2026-10-07/figures/build_baseline_response.py'
RENDERER_SHA = 'b647c99116fd5cb9ed9d4411b4dcfe3d95654a537ad013d5769f54fc55a676a8'
PROFILES = ('independent', 'persistent', 'imperfect')
LABELS = {'independent':'Independent', 'persistent':'Persistent', 'imperfect':'80% clearing'}
BRIDGE_KEYS = dict(zip(PROFILES, ('iid_perfect','persistent_perfect','iid_imperfect')))
PARAMS = {'independent':(0.,1.), 'persistent':(.5,1.), 'imperfect':(0.,.8)}
PRIMARY = 1-.05/3
GATE = 1-.05/6
FULL = '0072B2'
DISABLED = 'D55E00'
GRAY = '67747D'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)


def sha(value):
    return isinstance(value,str) and len(value)==64 and all(c in '0123456789abcdef' for c in value)


def number(value):
    if type(value) not in (int,float) or not math.isfinite(value):
        raise ValueError('A finite numeric report field is required')
    return value


def interval(value):
    if value is None:
        return None
    if not isinstance(value,list) or len(value)!=2 or number(value[0])>number(value[1]):
        raise ValueError('Invalid report interval')
    return value


def metric(value, uncertainty, level=.95):
    if value is not None:
        number(value)
    ci = interval(uncertainty['intervals'][format(level,'.12g')])
    undefined = uncertainty['undefined_bootstrap_draws']
    if type(undefined) is not int or undefined < 0 or uncertainty['bootstrap_draws']!=200000:
        raise ValueError('Unexpected bootstrap count')
    if uncertainty['conditional_on_defined_draws'] is not False or uncertainty['interval_unavailable_if_any_draw_undefined'] is not True:
        raise ValueError('Conditional-on-defined intervals are not allowed')
    if undefined and ci is not None:
        raise ValueError('An undefined draw must make the interval unavailable')
    return value, ci


def estimate(field, level=.95):
    return metric(field['estimate'],field['uncertainty'],level)


def fraction(row, level=.95):
    return metric(row['rate_prevented_fraction'],row['rate_prevented_fraction_uncertainty'],level)


def arm_rate(row, name):
    field = row['arms'][name]
    return metric(field['rate_per_1000_patient_days'],field['rate_uncertainty'])


def normalized(s):
    return {
        **dict(h0=.9,h1=.95,rho=0.,e=1.,opportunities='zone_separate',environment_multiplier=1.),**s}


def validate_freeze(path):
    value = read(path)
    ident = value.get('identity_sha256')
    base = {key:item for key,item in value.items() if key!='identity_sha256'}
    if ident!=hashlib.sha256(canonical(base).encode()).hexdigest():
        raise ValueError('Freeze identity mismatch')
    if not value.get('files') or not all(sha(v) for v in value['files'].values()):
        raise ValueError('Freeze requires file SHA256 bindings')
    for name,expected in value['files'].items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts or digest(ABM/relative)!=expected:
            raise ValueError('Frozen source/input changed: '+name)
    return value


def validate_report(path, design_path, freeze_path, analysis_freeze_path, module):
    r,d = read(path),read(design_path)
    freeze,analysis = validate_freeze(freeze_path),validate_freeze(analysis_freeze_path)
    n_conditions = {'holdout':3,'environment':6}[module]
    if r.get('schema_version')!=1 or r.get('status')!='complete' or r.get('module')!=module or r.get('is_synthetic_fixture'):
        raise ValueError('Only a complete fresh production report is allowed')
    if r.get('n_paired_jobs')!=512*n_conditions or r.get('bootstrap_draws')!=200000:
        raise ValueError('Exact planned sample and bootstrap counts are required')
    bindings = {'design_sha256':digest(design_path),'freeze_sha256':digest(freeze_path),
        'freeze_identity':freeze['identity_sha256'],'analysis_freeze_sha256':digest(analysis_freeze_path),
        'analysis_freeze_identity':analysis['identity_sha256']}
    if any(r.get(k)!=v for k,v in bindings.items()):
        raise ValueError('Report input hash/identity binding mismatch')
    if freeze['files'].get(str(Path(design_path).resolve().relative_to(ABM)))!=digest(design_path):
        raise ValueError('Design config is absent from study freeze')
    analyzer = 'research/full_icu_confirmation_2026-10-08/analyze.py'
    if r.get('source_sha256')!=analysis['files'].get(analyzer):
        raise ValueError('Report analyzer does not match the separate analysis freeze')
    for key in ('seedbook_sha256','proposal_sha256'):
        if not sha(r.get(key)) or r[key] not in freeze['files'].values():
            raise ValueError('Report declaration input is absent from freeze')
    if module=='holdout':
        for key in ('calibration_selection_sha256','calibration_report_sha256'):
            if not sha(r.get(key)) or r[key] not in freeze['files'].values():
                raise ValueError('Locked selection evidence is absent from holdout freeze')
    if (d['warmup_days'],d['observation_days'],d['model_overrides'])!=(90,56,{'initial_environment_contaminated':False,'prob_patient_sick':.03}):
        raise ValueError('Actual confirmation 90+56-day, 3% import, clean-start design required')
    seeds = d['seeds']
    if len(seeds)!=512 or len(set(seeds))!=512 or len(d['scenarios'])!=n_conditions:
        raise ValueError('Exact full design coverage required')
    rows = {row['scenario']['id']:row for row in r['scenarios']}
    declared = {s['id']:normalized(s) for s in d['scenarios']}
    expected = ({name+'_target_candidate' for name in PROFILES} if module=='holdout' else
        {f'independent_low{label}_env{env}' for label in ('075','100','125') for env in (0,1)})
    if len(r['scenarios'])!=n_conditions or set(rows)!=expected or set(declared)!=expected:
        raise ValueError('Missing, duplicate or unexpected condition')
    for ident,row in rows.items():
        if row['n']!=512 or row['seeds']!=seeds or normalized(row['scenario'])!=declared[ident]:
            raise ValueError('Condition settings or ordered seed coverage differ from declared design')
        s=row['scenario']; profile=ident.split('_')[0] if module=='holdout' else 'independent'
        if (s['rho'],s['e'],s['h0'],s['h1'],s['opportunities'])!=(*PARAMS[profile],.9,.95,'zone_separate'):
            raise ValueError('Profile mechanism differs from the figure label')
        fraction(row)
        for name in ('baseline','policy'):
            arm_rate(row,name)
        for name in ('mean_prevented','absolute_prevented_rate'):
            metric(row[name],row[name+'_uncertainty'])
    if not r.get('banks') or sum(bank['jobs'] for bank in r['banks'])!=512*n_conditions:
        raise ValueError('Bank receipts do not conserve the exact job count')
    for bank in r['banks']:
        if not sha(bank['manifest_sha256']) or len(bank['completed_records_sha256'])!=bank['jobs'] or not all(sha(v) for v in bank['completed_records_sha256'].values()):
            raise ValueError('Bank receipt lacks completed-record hash bindings')
        if not bank['sources'] or any(freeze['files'].get('research/full_icu_response_2026-10-08/array_worker.py' if k=='array_worker.py' else k)!=v for k,v in bank['sources'].items()):
            raise ValueError('Bank backend sources are not bound by the study freeze')
    if not math.isclose(r['primary_family_interval_level'],PRIMARY):
        raise ValueError('Primary interval level changed')
    if module=='holdout':
        matching = r['baseline_matching']; gates=matching['numerical_gates']
        expected_gates = {p+'_absolute_baseline' for p in PROFILES}|{
            'independent_over_persistent_baseline_ratio','independent_over_imperfect_baseline_ratio','persistent_over_imperfect_baseline_ratio'}
        if len(gates)!=6 or {g['name'] for g in gates}!=expected_gates or not math.isclose(matching['numerical_interval_level'],GATE):
            raise ValueError('All six matching gates are required')
        for gate in gates:
            estimate(gate,GATE)
            if type(gate['passed']) is not bool:
                raise ValueError('Explicit gate status required')
        if type(matching['valid_calibration']) is not bool or matching['established'] is not (matching['valid_calibration'] and all(g['passed'] for g in gates)):
            raise ValueError('Overall matching status contradicts reported gates')
    else:
        effects=r['primary_environment_effects']
        if len(effects)!=3 or [e['p'] for e in effects]!=[rows[f'independent_low{l}_env1']['scenario']['p'] for l in ('075','100','125')]:
            raise ValueError('Environmental contrasts do not match the three p values')
        for effect in effects:
            for name in ('F','mean_prevented','prevented_rate'):
                estimate(effect['full_minus_disabled'][name],PRIMARY)
                estimate(effect['full_minus_disabled'][name])
    return r,rows


def validate_bridge(path, config_path, binding_path, freeze_path):
    r,c,b = read(path),read(config_path),read(binding_path)
    freeze=validate_freeze(freeze_path)
    if r.get('module')!='bridge' or r.get('status')!='complete' or r.get('n_backgrounds')!=256 or r.get('seeds_completed')!=256 or r.get('seeds_planned')!=256 or r.get('is_synthetic_fixture'):
        raise ValueError('Complete actual 256-background bridge summary required')
    science=dict(description='Fixed 256 fresh p=0 backgrounds for exact finite-horizon ICU coefficient; no positive-p runs.',
        warmup_ticks=90*2592,observation_ticks=56*2592,h0=.9,h1=.95,opportunities='zone_separate',histogram_cap=2,
        profiles=[{'id':BRIDGE_KEYS[p],'rho':PARAMS[p][0],'e':PARAMS[p][1]} for p in PROFILES],model_overrides={'prob_patient_sick':.03})
    if {k:v for k,v in c.items() if k!='seeds'}!=science or len(c['seeds'])!=256 or len(set(c['seeds']))!=256:
        raise ValueError('Bridge input must be the actual 3% admission, 90+56-day clean p=0 design')
    required={'summary_sha256':digest(path),'config_sha256':digest(config_path),'seeds':c['seeds'],
        'freeze_sha256':digest(freeze_path),'freeze_identity':freeze['identity_sha256']}
    if b.get('status')!='complete' or any(b.get(k)!=v for k,v in required.items()):
        raise ValueError('Bridge summary/config/exact seeds/freeze binding mismatch')
    if any(r.get(k)!=v for k,v in required.items() if k!='summary_sha256'):
        raise ValueError('Bridge summary metadata disagree with exact config/seeds/freeze binding')
    expected_sources={f'src/{name}.py' for name in ('CPEmodel','agents','routes')}|{
        'research/full_icu_confirmation_2026-10-08/bridge/'+name+'.py' for name in ('run','coefficient','background')}
    if set(b.get('scientific_sources',{}))!=expected_sources or any(freeze['files'].get(k)!=v for k,v in b['scientific_sources'].items()):
        raise ValueError('Bridge scientific sources do not match their freeze')
    if r.get('scientific_sources')!=b['scientific_sources'] or r.get('source_sha256')!=b['scientific_sources'] or b.get('source_sha256')!=b['scientific_sources']:
        raise ValueError('Bridge summary and binding source maps differ')
    analysis_sha=r.get('analysis_code_sha256')
    if not sha(analysis_sha) or b.get('analysis_code_sha256')!=analysis_sha:
        raise ValueError('Bridge aggregation analysis-code hash binding required')
    if set(r['profiles'])!=set(BRIDGE_KEYS.values()):
        raise ValueError('All three bridge profiles required')
    for row in r['profiles'].values():
        if row['limiting_response'] is not None:number(row['limiting_response'])
        interval(row['limiting_response_normal_95_interval'])
        if row['stationary_q2_response_for_comparison_only'] is not None:number(row['stationary_q2_response_for_comparison_only'])
        for k in ('K0','K1'):
            number(row[k]['mean']);interval(row[k]['normal_approximation_95_interval'])
    return r


def load_scene():
    if digest(RENDERER)!=RENDERER_SHA:
        raise ValueError('Shared native renderer hash changed')
    spec=importlib.util.spec_from_file_location('_confirmation_native_scene',RENDERER)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    class PrintScene(m.Scene):
        def text(self,x,y,value,size=8,**kwargs):
            # At a 6.5-inch manuscript width, a 640-point scene scales by
            # 0.731. An 11-point native minimum therefore stays above 8 pt.
            return super().text(x,y,value,size=max(size,11),**kwargs)
    return PrintScene


def scaled(pair, factor=1):
    value,ci=pair
    return (None if value is None else value*factor,None if ci is None else [x*factor for x in ci])


def bounds(pairs, references=()):
    values=list(references)+[x for value,ci in pairs for x in ([value] if ci is None else ci+[value]) if x is not None]
    if not values:return (-1,1)
    lo,hi=min(values),max(values)
    padding=max((hi-lo)*.12,abs(hi)*.03,.001)
    return lo-padding,hi+padding


def axis(scene,x0,x1,y0,y1,pairs,label,references=(),band=None,nonnegative=False,percent=False):
    reference_values=list(references)+(list(band) if band is not None else [])
    lo,hi=bounds(pairs,reference_values)
    actual=[x for value,ci in pairs for x in ([value] if ci is None else ci+[value]) if x is not None]
    anchors=actual+reference_values
    if 0 in references and anchors:
        if min(anchors)>=0:lo=0.
        if max(anchors)<=0:hi=0.
    if percent:
        lo=0 if not actual or min(actual)>=0 else min(lo,0)
        hi=100 if not actual or max(actual)<=100 else max(hi,100)
    if hi<=lo:
        # An exactly zero estimated effect and zero CI is valid. Give it a
        # finite plotting span while keeping the zero point/reference visible.
        lo,hi=(0.,1.) if nonnegative else (-1.,1.)
    divisions=5 if percent else 3
    step=10**math.floor(math.log10((hi-lo)/divisions))
    step*=next(x for x in (1,2,2.5,5,10) if x*step>=(hi-lo)/divisions)
    lo=math.floor(lo/step)*step;hi=math.ceil(hi/step)*step
    if nonnegative:lo=0.
    xx=lambda v:x0+(v-lo)/(hi-lo)*(x1-x0)
    if band:scene.rectangle(xx(band[0]),y0,xx(band[1])-xx(band[0]),y1-y0,'EDF3F7')
    for index in range(int(round((hi-lo)/step))+1):
        v=lo+index*step
        x=xx(v);scene.line(x,y0,x,y1,'E5EAED',.5)
        scene.text(x,y0-12,f'{v:.3g}',size=9,align='center',color=GRAY)
    scene.line(x0,y0,x1,y0,GRAY,.7)
    scene.text((x0+x1)/2,y0-31,label,size=9.5,align='center')
    return xx


def vdash(scene,x,y0,y1,color=GRAY):
    for y in range(int(y0),int(y1),8):scene.line(x,y,x,min(y+4,y1),color,.75)


def point(scene,xx,pair,y,color,marker='o',width=1.3):
    value,ci=pair
    if ci is not None:
        scene.line(xx(ci[0]),y,xx(ci[1]),y,color,width)
        for x in ci:scene.line(xx(x),y-3,xx(x),y+3,color,.8)
    if value is not None:scene.point(xx(value),y,color,5,marker)
    if value is None or ci is None:
        scene.text(xx(value) if value is not None else xx(0),y-12,'estimate unavailable' if value is None else 'CI unavailable',size=9,color=GRAY)


def header(scene,title,subtitle,fixture=False):
    scene.text(16,scene.height-22,title,size=14,weight='bold')
    scene.text(16,scene.height-43,subtitle,size=9.5,color=GRAY)
    if fixture:
        scene.rectangle(10,scene.height-78,scene.width-20,21,'FFF1E8')
        scene.text(scene.width/2,scene.height-67,'SYNTHETIC LAYOUT FIXTURE - NOT STUDY RESULTS',size=11,align='center',weight='bold',color='9C4300')


def draw_holdout(Scene,r,rows,fixture=False):
    s=Scene(640,500)
    header(s,'Fresh holdout: measured ICU rates and prevention','90 warmup + 56 observation days; 512 paired seeds/profile; completion 90% to 95%',fixture)
    s.text(16,395,'A  Acquired colonization rates',size=11.5,weight='bold')
    s.text(428,395,'B  Rate reduction F',size=11.5,weight='bold')
    s.text(16,377,'95% bootstrap CIs; band = 2.7 to 3.3; dashed target = 3',size=9,color=GRAY)
    rates=[arm_rate(rows[p+'_target_candidate'],a) for p in PROFILES for a in ('baseline','policy')]
    fs=[scaled(fraction(rows[p+'_target_candidate']),100) for p in PROFILES]
    xr=axis(s,160,397,205,360,rates,'Acquisitions / 1000 patient-days',(0,3),band=(2.7,3.3),nonnegative=True)
    xf=axis(s,462,621,205,360,fs,'Prevention (%)',(0,),percent=True)
    vdash(s,xr(3),205,360);vdash(s,xf(0),205,360)
    for p,y in zip(PROFILES,(335,280,225)):
        row=rows[p+'_target_candidate']
        s.text(149,y+4,LABELS[p],size=10,align='right')
        s.text(149,y-12,f"p = {row['scenario']['p']:.5f}",size=9,align='right',color=GRAY)
        point(s,xr,arm_rate(row,'baseline'),y+5,FULL,'o')
        point(s,xr,arm_rate(row,'policy'),y-5,DISABLED,'s')
        point(s,xf,scaled(fraction(row),100),y,FULL,'o')
    s.point(170,161,FULL,5,'o');s.text(180,161,'Baseline 90%',size=9.5)
    s.point(293,161,DISABLED,5,'s');s.text(303,161,'Policy 95%',size=9.5)
    m=r['baseline_matching']
    status='ESTABLISHED' if m['established'] else 'NOT ESTABLISHED'
    s.text(16,134,'Baseline matching: '+status,size=11,weight='bold')
    s.text(16,115,'Valid calibration: '+('PASS' if m['valid_calibration'] else 'NOT PASSED')+'; six matching CIs use 99.167%',size=9.5)
    gates={g['name']:g for g in m['numerical_gates']}
    pairs=(('independent','persistent'),('independent','imperfect'),('persistent','imperfect'))
    for p,(a,b),y in zip(PROFILES,pairs,(95,78,61)):
        astatus='PASS' if gates[p+'_absolute_baseline']['passed'] else 'NOT PASSED'
        pstatus='PASS' if gates[a+'_over_'+b+'_baseline_ratio']['passed'] else 'NOT PASSED'
        s.text(16,y,LABELS[p]+' target: '+astatus,size=9)
        s.text(295,y,LABELS[a]+' / '+LABELS[b]+' ratio: '+pstatus,size=9)
    s.text(16,34,'R0/R1 and F are measured at the locked probabilities; points are not forced onto the target.',size=9,color=GRAY)
    s.text(16,17,'Monte Carlo CIs conditional on p; calibration and clinical-input uncertainty excluded.',size=9,color=GRAY)
    return s


def draw_environment(Scene,r,rows,fixture=False):
    s=Scene(640,660)
    header(s,'Environmental intervention: relative and absolute prevention','Independent; 512 paired seeds/condition; 90 + 56 days; full vs environment disabled',fixture)
    s.text(16,553,'A  Rate reduction F',size=11.5,weight='bold')
    s.text(385,553,'B  Full minus disabled F',size=11.5,weight='bold')
    s.text(16,534,'Condition CIs: 95% bootstrap',size=9.5,color=GRAY)
    s.text(385,534,'Outer 98.333%; inner 95% contrast CIs',size=9.5,color=GRAY)
    ordered=[rows[f'independent_low{l}_env{e}'] for l in ('075','100','125') for e in (1,0)]
    fractions=[scaled(fraction(row),100) for row in ordered]
    effects=r['primary_environment_effects']
    contrasts=[scaled(estimate(e['full_minus_disabled']['F'],PRIMARY),100) for e in effects]
    xa=axis(s,134,344,350,510,fractions,'Prevention (%)',(0,),percent=True)
    xb=axis(s,429,621,350,510,contrasts,'Difference (percentage points)',(0,))
    vdash(s,xa(0),350,510);vdash(s,xb(0),350,510)
    for label,effect,y in zip(('075','100','125'),effects,(485,429,373)):
        s.text(121,y,f"p = {effect['p']:.5f}",size=9.5,align='right')
        for env,offset,color,marker in ((1,5,FULL,'o'),(0,-5,DISABLED,'open')):
            point(s,xa,scaled(fraction(rows[f'independent_low{label}_env{env}']),100),y+offset,color,marker)
        point(s,xb,scaled(estimate(effect['full_minus_disabled']['F'],PRIMARY),100),y,GRAY,'o',.8)
        point(s,xb,scaled(estimate(effect['full_minus_disabled']['F']),100),y,FULL,'o',2.5)
    s.point(144,303,FULL,5,'o');s.text(154,303,'Full',size=9.5)
    s.point(230,303,DISABLED,5,'open');s.text(240,303,'Environment disabled',size=9.5)
    s.text(16,274,'C  Absolute prevention rate',size=11.5,weight='bold')
    s.text(363,274,'D  Mean colonizations prevented',size=11.5,weight='bold')
    s.text(16,255,'Baseline minus policy; 95% bootstrap CIs',size=9.5,color=GRAY)
    s.text(363,255,'Paired count difference / 56-day ICU run',size=9.5,color=GRAY)
    def absolute(row,key):return metric(row[key],row[key+'_uncertainty'])
    xc=axis(s,134,344,101,237,[absolute(row,'absolute_prevented_rate') for row in ordered],'Prevented / 1000 patient-days',(0,))
    xd=axis(s,429,621,101,237,[absolute(row,'mean_prevented') for row in ordered],'Mean prevented / 56 days',(0,))
    vdash(s,xc(0),101,237);vdash(s,xd(0),101,237)
    for label,effect,y in zip(('075','100','125'),effects,(214,165,116)):
        s.text(121,y,f"p = {effect['p']:.5f}",size=9.5,align='right')
        for env,offset,color,marker in ((1,5,FULL,'o'),(0,-5,DISABLED,'open')):
            row=rows[f'independent_low{label}_env{env}']
            point(s,xc,absolute(row,'absolute_prevented_rate'),y+offset,color,marker)
            point(s,xd,absolute(row,'mean_prevented'),y+offset,color,marker)
    s.text(16,53,'F = 1 - R1/R0; absolute rate = R0 - R1; count = mean paired C0 - C1.',size=9.5)
    s.text(16,34,'Hand-site transmission is disabled; later patient spread can also change.',size=9,color=GRAY)
    s.text(16,17,'This does not isolate one pickup pathway. No monotonic fit; CIs measure Monte Carlo precision.',size=9,color=GRAY)
    return s


def draw_bridge(Scene,holdout,rows,bridge,fixture=False):
    s=Scene(640,490)
    header(s,'Full ICU model: rare-transmission limit and direct simulations','Clean reservoirs; imports 3%; 90 + 56 days; p = 0 n = 256; holdout n = 512/profile',fixture)
    pairs=[]
    for p in PROFILES:
        b=bridge['profiles'][BRIDGE_KEYS[p]]
        pairs.extend([scaled(fraction(rows[p+'_target_candidate']),100),
            scaled((b['limiting_response'],b['limiting_response_normal_95_interval']),100),
            scaled((b['stationary_q2_response_for_comparison_only'],None),100)])
    x=axis(s,211,621,163,367,pairs,'Rate reduction F (%)',(0,),percent=True)
    vdash(s,x(0),163,367)
    for p,y in zip(PROFILES,(330,260,190)):
        b=bridge['profiles'][BRIDGE_KEYS[p]];row=rows[p+'_target_candidate']
        s.text(198,y+4,LABELS[p],size=10.5,align='right')
        s.text(198,y-13,f"Holdout p = {row['scenario']['p']:.5f}",size=9,align='right',color=GRAY)
        point(s,x,scaled(fraction(row),100),y+13,FULL,'o')
        point(s,x,scaled((b['limiting_response'],b['limiting_response_normal_95_interval']),100),y,DISABLED,'s')
        q=b['stationary_q2_response_for_comparison_only']
        if q is not None:s.point(x(100*q),y-13,GRAY,5,'d')
        else:s.text(215,y-13,'q2 unavailable',size=9,color=GRAY)
    for y,color,marker,label in ((108,FULL,'o','Measured holdout F: 95% paired bootstrap CI'),
        (87,DISABLED,'s','p = 0 coefficient limit: 1 - K1/K0; 95% normal/delta Monte Carlo CI'),
        (66,GRAY,'d','Stationary q2 benchmark: fixed comparator; no Monte Carlo CI')):
        s.point(22,y,color,5,marker);s.text(33,y,label,size=9.5)
    s.text(16,43,'K integrates survival on imported-source bridges in the actual finite-horizon p = 0 ICU.',size=9,color=GRAY)
    s.text(16,25,'No finite-p accuracy guarantee. q2 assumes stationary two-opportunity bridges.',size=9,color=GRAY)
    s.text(16,9,'Normal/delta CI is asymptotic; intervals measure synthetic-model Monte Carlo precision.',size=9,color=GRAY)
    return s


def save_scene(scene, output, name, export_pdf=True):
    # Comment every newline inside the resized picture. Spaces after picture
    # and endgroup otherwise add width and shift the native figure off-center.
    lines=scene.tex+[r'\end{picture}',r'\endgroup']
    inline='\n'.join(line.rstrip()+'%' for line in lines)+'\n'
    portable='% Requires xcolor and graphicx; report-bound primitives only.\n'+r'\resizebox{\linewidth}{!}{%'+'\n'+inline+'}\n'
    (output/(name+'.tex')).write_text(portable)
    standalone=(r'\documentclass{article}'+'\n'+r'\usepackage{xcolor,graphicx}'+'\n'+
        '\\usepackage[paperwidth=%gpt,paperheight=%gpt,margin=20pt]{geometry}\n'%(scene.width+40,scene.height+40)+
        r'\pagestyle{empty}\setlength{\parindent}{0pt}'+'\n'+r'\begin{document}\noindent'+'\n'+portable+r'\end{document}'+'\n')
    (output/(name+'_standalone.tex')).write_text(standalone)
    scene.figure.savefig(output/(name+'.png'),dpi=240,facecolor='white')
    # PDF is vector output from the same primitives. Native TeX compilation is
    # an independent QA step, not simulated by this Matplotlib export.
    if export_pdf:
        fixed_date=datetime(2026,10,8,tzinfo=timezone.utc)
        scene.figure.savefig(output/(name+'.pdf'),facecolor='white',metadata={'Title':name,'Creator':'Confirmation report figure builder',
            'CreationDate':fixed_date,'ModDate':fixed_date})
    import matplotlib.pyplot as plt
    plt.close(scene.figure)


def public_path(path):
    path=Path(path).resolve()
    try:return str(path.relative_to(ABM))
    except ValueError:return path.name


def leading_rate_rows(holdout_rows,environment_rows,bridge):
    """Presentation formula explicitly requested; not a refit or reanalysis.

    The clean p=0 background is shared by full/disabled environment, whose
    environmental acquisition terms enter beyond second order. These are
    point calculations from reported K and exposure means, without CIs.
    """
    days=number(bridge['observation_patient_days']['mean'])
    if days<=0:raise ValueError('A positive reported p=0 exposure mean is required')
    result=[]
    entries=[('holdout',p,holdout_rows[p+'_target_candidate']) for p in PROFILES]
    entries += [('environment','independent',environment_rows[f'independent_low{label}_env{env}'])
        for label in ('075','100','125') for env in (1,0)]
    for module,profile,row in entries:
        s=row['scenario'];k=bridge['profiles'][BRIDGE_KEYS[profile]]
        p=number(s['p']);k0=number(k['K0']['mean']);k1=number(k['K1']['mean'])
        item=dict(module=module,condition=s['id'],profile=profile,p=p,environment_multiplier=s['environment_multiplier'],
            K0_mean=k0,K1_mean=k1,p0_observation_patient_days_mean=days,
            R0_leading=1000*p*p*k0/days,R1_leading=1000*p*p*k1/days,
            limiting_F=k['limiting_response'],stationary_q2_F=k['stationary_q2_response_for_comparison_only'])
        for name,field in (('R0',arm_rate(row,'baseline')),('R1',arm_rate(row,'policy')),('F',fraction(row))):
            value,ci=field;item[name+'_measured']=value
            item[name+'_measured_CI95_low']=None if ci is None else ci[0]
            item[name+'_measured_CI95_high']=None if ci is None else ci[1]
        ci=k['limiting_response_normal_95_interval']
        item['limiting_F_normal_CI95_low']=None if ci is None else ci[0]
        item['limiting_F_normal_CI95_high']=None if ci is None else ci[1]
        result.append(item)
    return result


def save_leading_table(output,rows):
    with (output/'leading_rate_comparison.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    def f(value):return '--' if value is None else f'{value:.3f}'
    def measured(row,key):
        a,b=row[key+'_measured_CI95_low'],row[key+'_measured_CI95_high']
        return f(row[key+'_measured'])+(' [CI unavailable]' if a is None else f' [{f(a)}, {f(b)}]')
    def label(row):
        if row['module']=='holdout':return 'Holdout '+LABELS[row['profile']].replace('%',r'\%')
        return 'Full environment' if row['environment_multiplier']==1 else 'Environment disabled'
    lines=[r'% All numbers use fresh reports; leading rates are deterministic presentation formulas.']
    for key in ('R0','R1'):
        h=key[-1]
        lines += [r'\begin{tabular}{lrrr}\hline',
            r'Condition & $p$ & Measured $R_'+h+r'$ [95\% CI] & Leading $R_'+h+r'$\\\hline']
        for row in rows:lines.append(' & '.join((label(row),f(row['p']),measured(row,key),f(row[key+'_leading'])))+r'\\')
        lines += [r'\hline\end{tabular}',r'\par\smallskip']
    lines += [
        r'\begin{tabular}{lrrrr}\hline',r'Condition & $p$ & Measured $F$ [95\% CI] & Limiting $F$ & Stationary $q_2$ $F$\\\hline']
    for row in rows:lines.append(' & '.join((label(row),f(row['p']),measured(row,'F'),f(row['limiting_F']),f(row['stationary_q2_F'])))+r'\\')
    lines += [r'\hline\end{tabular}',r'\par\smallskip',
        r'Leading rates: $1000p^2\overline K_h/\overline D^{(0)}$. Units: acquisitions per 1000 patient-days.',
        r'Here $\overline K_h$ is the reported mean coefficient and $\overline D^{(0)}$ is the reported mean observation patient-days at $p=0$.',
        r'Leading values are point calculations; measured CIs are paired bootstrap intervals.',
        r'The CSV supplies the limiting-response normal/delta CIs and all source means.',
        r'No finite-$p$ error certification or inference from CI overlap/nonoverlap is implied.']
    (output/'leading_rate_comparison.tex').write_text('\n'.join(lines)+'\n')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('holdout-report','holdout-design','holdout-freeze','environment-report','environment-design','environment-freeze',
                 'analysis-freeze','bridge-report','bridge-config','bridge-binding','bridge-freeze','output'):
        p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args(argv)
    output=args.output.resolve()
    if not output.is_relative_to(HERE.resolve()):raise ValueError('Outputs must stay beneath the new figures directory')
    if output.exists():raise FileExistsError('Use a fresh output directory')
    sources={key:{'path':public_path(value),'sha256':digest(value)} for key,value in vars(args).items() if key!='output'}
    holdout,hr=validate_report(args.holdout_report,args.holdout_design,args.holdout_freeze,args.analysis_freeze,'holdout')
    environment,er=validate_report(args.environment_report,args.environment_design,args.environment_freeze,args.analysis_freeze,'environment')
    bridge=validate_bridge(args.bridge_report,args.bridge_config,args.bridge_binding,args.bridge_freeze)
    table=leading_rate_rows(hr,er,bridge)
    if set(next(iter(hr.values()))['seeds']).intersection(read(args.bridge_config)['seeds']):
        raise ValueError('Bridge and holdout seeds must be disjoint')
    Scene=load_scene()
    scenes=[('holdout_rates_and_prevention',draw_holdout(Scene,holdout,hr)),
        ('environment_estimands',draw_environment(Scene,environment,er)),
        ('bridge_limit_comparison',draw_bridge(Scene,holdout,hr,bridge))]
    if any(digest(value)!=sources[key]['sha256'] for key,value in vars(args).items() if key!='output'):
        raise ValueError('An explicit input changed during figure generation')
    # All gates precede output creation, so a partial report produces no figures.
    output.mkdir(parents=True)
    for name,scene in scenes:save_scene(scene,output,name)
    save_leading_table(output,table)
    manifest=dict(schema_version=1,status='complete_render',is_synthetic_fixture=False,sources=sources,
        builder_sha256=digest(Path(__file__)),native_renderer={'path':public_path(RENDERER),'sha256':RENDERER_SHA},
        estimates_recomputed=False,interpolation=False,simulations_run=False,
        deterministic_presentation_formula='R_h_leading = 1000*p^2*reported_mean_K_h/reported_mean_p0_observation_patient_days; point only',
        derived_table_rows=9,finite_p_error_certified=False,
        counts={'holdout_seeds_per_profile':512,'environment_seeds_per_condition':512,'p0_backgrounds':256},
        confidence={'descriptive':.95,'primary_family':PRIMARY,'matching_gate':GATE,'bridge':'95% asymptotic normal/delta'},
        estimands={'holdout':'Measured pooled R0/R1 and 1-R1/R0 at locked p',
            'environment':'Fixed-p full minus disabled relative response; absolute paired rate/count prevention',
            'bridge':'1-K1/K0 at p=0; stationary q2 comparator; positive-p measured holdout F'},
        native_compilation='Not performed by builder; standalone TeX provided for independent compilation/QA',
        native_minimum_font_pt=11,minimum_font_pt_at_6p5_inch_width=11*468/640,
        outputs_sha256={file.name:digest(file) for file in sorted(output.iterdir()) if file.is_file()})
    (output/'manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2,allow_nan=False)+'\n')
    print(canonical({'status':'complete_render','output':public_path(output),'figures':len(scenes)}))


if __name__=='__main__':main()
