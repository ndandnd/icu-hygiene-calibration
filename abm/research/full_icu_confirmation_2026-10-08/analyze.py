#!/usr/bin/env python3
"""Prespecified confirmation analysis; complete declared banks only.

No simulator calls. Bootstrap multiplicities resample whole seed vectors with
common weights across every condition/arm. Undefined draws remain visible and
make confidence intervals unavailable, including matching gates.
"""
import argparse
from collections import Counter
import csv
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]
RESPONSE = ABM / 'research/full_icu_response_2026-10-08'
BOOTSTRAP_DRAWS = 200000
BOOTSTRAP_SEED = 2026100817
GATE_LEVEL = 1-.05/6
PRIMARY_LEVEL = 1-.05/3
METRICS = ('R0', 'R1', 'F', 'count_F', 'mean_prevented', 'prevented_rate',
           'R0_susceptible', 'R1_susceptible', 'early_R0', 'late_R0',
           'late_minus_early_R0', 'late_over_early_R0')
TOTAL_KEYS = ('acquired_colonizations', 'imported_colonizations', 'patient_days',
              'susceptible_patient_days', 'environment_contaminated_site_days')


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


builder = load_module('_confirmation_design_builder', HERE/'build_design.py')
backend = load_module('_confirmation_frozen_response_runner', RESPONSE/'run.py')


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite(value, label, integer=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError('Invalid nonnegative finite '+label)
    if integer and type(value) is not int:
        raise ValueError('Expected integer '+label)
    return value


def divide(numerator, denominator):
    a, b = np.broadcast_arrays(np.asarray(numerator, dtype=float), np.asarray(denominator, dtype=float))
    return np.divide(a, b, out=np.full(a.shape, np.nan), where=np.isfinite(b) & (b != 0))


def scalar(value):
    return float(value) if np.isfinite(value) else None


def uncertainty(values, levels=(.95, PRIMARY_LEVEL, GATE_LEVEL)):
    values = np.asarray(values, dtype=float)
    undefined = int(np.count_nonzero(~np.isfinite(values)))
    intervals = {}
    for level in levels:
        intervals[format(level, '.12g')] = None if undefined or not len(values) else np.quantile(values, [(1-level)/2, (1+level)/2]).tolist()
    return dict(interval95=intervals.get('0.95'), intervals=intervals,
                attained_interval_widths={key:None if value is None else value[1]-value[0] for key,value in intervals.items()},
                bootstrap_draws=len(values), undefined_bootstrap_draws=undefined,
                interval_unavailable_if_any_draw_undefined=True,
                conditional_on_defined_draws=False)


def estimate(point, samples, levels=(.95, PRIMARY_LEVEL, GATE_LEVEL)):
    return dict(estimate=scalar(point), uncertainty=uncertainty(samples, levels))


def arm(row, name):
    value = row['arms'][name]
    if 'shared_with' in value:
        if name != 'policy' or value != {'shared_with': 'baseline'}:
            raise ValueError('Invalid shared-arm reference')
        return row['arms']['baseline']
    return value


def validate_summary(value, label):
    for key in TOTAL_KEYS:
        finite(value[key], label+' '+key, integer=key in ('acquired_colonizations', 'imported_colonizations'))
    finite(value['days'], label+' days')
    if value['susceptible_patient_days'] > value['patient_days']+1e-7:
        raise ValueError('Susceptible exposure exceeds patient exposure')
    if value['patient_days'] > 30*value['days']+1e-7 or value['environment_contaminated_site_days'] > 30*value['days']+1e-7:
        raise ValueError('Exposure exceeds fixed 30-bed/site topology')
    for key in ('acquisitions_by_route', 'acquisitions_by_opportunities', 'acquisitions_by_route_and_opportunities'):
        counts = value['response_counts'][key]
        if not isinstance(counts, dict):
            raise ValueError('Acquisition histogram must be a mapping')
        for count in counts.values():
            finite(count, label+' source counter', integer=True)
        if sum(counts.values()) != value['acquired_colonizations']:
            raise ValueError('Acquisition source counters fail conservation')


def validate_blocks(blocks, summary, expected_days):
    validate_summary(summary, 'summary')
    boundary = 0.
    for block in blocks:
        validate_summary(block, 'daily block')
        if block['day_start'] != boundary or block['day_end'] != boundary+1 or block['days'] != 1:
            raise ValueError('Daily blocks overlap, omit days, or span multiple days')
        boundary = block['day_end']
    if boundary != expected_days or summary['days'] != expected_days:
        raise ValueError('Daily blocks do not cover declared horizon')
    for key in TOTAL_KEYS:
        if not math.isclose(math.fsum(block[key] for block in blocks), summary[key], rel_tol=1e-10, abs_tol=1e-7):
            raise ValueError('Daily totals disagree with summary: '+key)
    for key in ('acquisitions_by_route', 'acquisitions_by_opportunities', 'acquisitions_by_route_and_opportunities'):
        total = Counter()
        for block in blocks:
            total.update(block['response_counts'][key])
        if +total != +Counter(summary['response_counts'][key]):
            raise ValueError('Daily acquisition histograms disagree with summary')


def normalized(config):
    return backend.normalize_config(config)


def validate_design(module, raw_config, book, proposal, selection=None, calibration_report=None,
                    selection_path=None, calibration_report_path=None, proposal_path=None):
    expected_book = builder.seedbook(proposal)
    if book != expected_book:
        raise ValueError('Seedbook differs from audited proposal allocations')
    seeds = expected_book['allocations'][module]['all']
    if module == 'calibration':
        cases = [builder.scenario(f'independent_calibration_p{i}', p, h1=.9) for i,p in enumerate(builder.GRID)]
    elif module == 'environment':
        cases = [builder.scenario(f'independent_low{label}_env{env}', p, env=float(env))
                 for label,p in zip(('075','100','125'), builder.OLD_GRID) for env in (1,0)]
    else:
        if selection is None or calibration_report is None:
            raise ValueError('Holdout requires locked selection and complete calibration report')
        recomputed = builder.select_baseline_report(calibration_report, expected_book['allocations']['calibration']['all'])
        # Canonical comparison accepts JSON tuple/list roundtrips while
        # distinguishing a boolean gate from an integer/string lookalike.
        if any(canonical(selection.get(key)) != canonical(value) for key,value in recomputed.items()):
            raise ValueError('Locked selection differs from baseline-only calibration')
        if calibration_report.get('module') != 'calibration' or calibration_report.get('status') != 'complete':
            raise ValueError('Holdout requires a validated complete calibration report')
        if calibration_report_path is not None and selection.get('report_sha256') != digest(calibration_report_path):
            raise ValueError('Locked calibration report hash differs')
        if proposal_path is not None and selection.get('proposal_sha256') != digest(proposal_path):
            raise ValueError('Locked proposal hash differs')
        ps = {'independent': recomputed['candidate_p'], **builder.FIXED_P}
        cases = [builder.scenario(f'{name}_target_candidate', p, name) for name,p in ps.items()]
    actual = normalized(raw_config)
    expected = normalized(builder.config(cases, seeds, raw_config.get('description', '')))
    if actual != expected:
        raise ValueError('Analysis config differs from complete prespecified family design')
    return actual


def checked_manifest(manifest):
    identity_key = 'identity' if 'identity' in manifest else 'identity_sha256'
    identity = manifest.get(identity_key)
    unsigned = {key:value for key,value in manifest.items() if key != identity_key}
    if hashlib.sha256(canonical(unsigned).encode()).hexdigest() != identity:
        raise ValueError('Manifest identity does not hash its contents')
    return identity_key, identity


def verify_analysis_freeze(path):
    freeze = json.loads(path.read_text())
    checked_manifest(freeze)
    if freeze.get('study') != HERE.name or freeze.get('scope') != 'confirmation_analyzer_before_outcome_analysis':
        raise ValueError('Wrong analysis freeze scope')
    required = {str(path.relative_to(ABM)) for path in (HERE/'analyze.py',HERE/'test_analyze.py',HERE/'build_design.py')}
    if not required <= set(freeze['files']):
        raise ValueError('Analysis freeze omits analyzer/test/selection dependency')
    for name,sha in freeze['files'].items():
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts or (ABM/relative).is_symlink() or digest(ABM/relative) != sha:
            raise ValueError('Frozen analysis source changed: '+name)
    return freeze


def expected_model_kwargs(config, job):
    s = job['scenario']
    kwargs = dict(backend.FIXED_MODEL, **config['model_overrides'])
    if 'initial_environment_contaminated' in s:
        kwargs['initial_environment_contaminated'] = s['initial_environment_contaminated']
    kwargs.update(seed=job['seed'], prob_transmission=s['p'], icu_hcw_wash_rate=s['h0'],
                  outside_hcw_wash_rate=s['h0'], hygiene_opportunities=s['opportunities'],
                  hygiene_autocorrelation=s['rho'], hygiene_efficacy=s['e'],
                  environment_multiplier=s['environment_multiplier'])
    return kwargs


def verify_binding(directory, manifest, context):
    receipt = json.loads((directory/'confirmation_binding.json').read_text())
    root = Path(context.get('root',ABM)).resolve()
    freeze = context['freeze']
    required = dict(schema_version=1,study=HERE.name,freeze_sha256=digest(context['freeze_path']),
        freeze_identity=freeze['identity_sha256'],seedbook_sha256=digest(context['seedbook_path']),
        output_scope=str(directory.resolve().relative_to(root)),task_count=manifest['task_count'],
        backend='research/full_icu_response_2026-10-08/array_worker.py')
    if any(receipt.get(key) != value for key,value in required.items()):
        raise ValueError('Confirmation receipt differs from frozen study/bank binding')
    name = receipt['config_path']
    relative = Path(name)
    if relative.is_absolute() or '..' in relative.parts:
        raise ValueError('Unsafe confirmation config path')
    path = root/relative
    if path.is_symlink() or receipt['config_sha256'] != digest(path) or freeze['files'].get(name) != digest(path):
        raise ValueError('Confirmation config is not the exact frozen input')
    if json.loads(path.read_text()) != manifest['config']:
        raise ValueError('Confirmation config differs from runtime manifest')
    scopes = context['output_scopes']
    if scopes['banks'].get(path.stem) != receipt['output_scope']:
        raise ValueError('Confirmation output scope differs from frozen allocation')
    family,site = path.stem.rsplit('_',1)
    if family != context['module'] or site not in ('local','cluster'):
        raise ValueError('Unexpected confirmation family/site config')
    return site


def load_banks(paths, expected_config, expected_sources, book, module, confirmation_context=None):
    """Validate headers/completeness before opening any numerical outcome rows."""
    expected = normalized(expected_config)
    jobs_union = {(s['id'], seed): s for s in expected['scenarios'] for seed in expected['seeds']}
    declared, headers, records, receipts = {}, [], {}, []
    for directory in paths:
        directory = Path(directory)
        manifest_path = directory/'manifest.json'
        manifest = json.loads(manifest_path.read_text())
        identity_key, identity = checked_manifest(manifest)
        sources = manifest.get('sources', manifest.get('source_sha256'))
        if sources != expected_sources:
            raise ValueError('Scientific source identities differ from frozen expected sources')
        site = None if confirmation_context is None else verify_binding(directory,manifest,confirmation_context)
        configured = manifest['resolved_config']
        if normalized(manifest.get('config', manifest.get('input_config'))) != configured:
            raise ValueError('Manifest raw/resolved configs disagree')
        if {k:v for k,v in configured.items() if k not in ('seeds','description')} != {k:v for k,v in expected.items() if k not in ('seeds','description')}:
            raise ValueError('Bank differs from family scientific config')
        bank_seeds = configured['seeds']
        allocations = book['allocations'][module]
        if site in ('local', 'cluster'):
            if bank_seeds != allocations[site]:
                raise ValueError('Site seeds differ from preallocated cohort membership')
        elif bank_seeds not in (allocations['local'], allocations['cluster'], allocations['all']):
            raise ValueError('Bank seeds are not a declared cohort/site allocation')
        jobs = {job['job_id']:job for job in manifest['jobs']}
        if len(jobs) != len(manifest['jobs']) or manifest['jobs'] != backend.build_jobs(configured):
            raise ValueError('Manifest jobs differ from full configured Cartesian design')
        for job in manifest['jobs']:
            key = job['scenario']['id'], job['seed']
            if key in declared:
                raise ValueError('Duplicate declared condition/seed across sites')
            declared[key] = job['scenario']
        files = sorted((directory/'completed').glob('*.json'))
        names = {f'job_{job_id:08d}.json' for job_id in jobs} if identity_key == 'identity' else {f'{job_id:06d}.json' for job_id in jobs}
        if {path.name for path in files} != names:
            raise ValueError('Complete expected job files required before analysis')
        headers.append((manifest_path, manifest, identity_key, identity, jobs, files))
    if declared != jobs_union:
        raise ValueError('Banks do not exhaust complete family condition/seed design')
    for manifest_path, manifest, identity_key, identity, jobs, files in headers:
        configured = manifest['resolved_config']
        seen = set()
        completed_hashes = {}
        for path in files:
            payload = path.read_bytes()
            row = json.loads(payload)
            completed_hashes[path.name] = hashlib.sha256(payload).hexdigest()
            job_id = row.get('job_id')
            if job_id not in jobs or job_id in seen:
                raise ValueError('Unexpected or duplicate completed job ID')
            job = jobs[job_id]
            expected_name = f'job_{job_id:08d}.json' if identity_key == 'identity' else f'{job_id:06d}.json'
            if path.name != expected_name:
                raise ValueError('Completed filename differs from job ID')
            if row.get('status') != 'complete' or row.get(identity_key) != identity or any(row.get(key) != value for key,value in job.items()):
                raise ValueError('Completed row differs from its manifest job/identity')
            if row['model_kwargs'] != expected_model_kwargs(configured, job):
                raise ValueError('Actual model kwargs differ from fixed declared model')
            if identity_key == 'identity':
                software = {key:value for key,value in row['execution'].items() if key != 'platform'}
                if software != manifest['software']:
                    raise ValueError('Record execution software differs from manifest')
            if row['warmup_snapshot']['total_ticks'] != expected['warmup_ticks']:
                raise ValueError('Warmup tick differs from declared horizon')
            if row.get('executed_ticks_per_trajectory') != expected['warmup_ticks']+expected['observation_ticks']:
                raise ValueError('Executed trajectory ticks differ from declared horizon')
            if row.get('warmup_executions') != 1:
                raise ValueError('Expected one shared warmup execution')
            shared = job['scenario']['h0'] == job['scenario']['h1']
            if row.get('distinct_observation_trajectories') != (1 if shared else 2):
                raise ValueError('Unexpected observed trajectory count')
            if shared != ('shared_with' in row['arms']['policy']):
                raise ValueError('Shared policy representation differs from declared h')
            validate_blocks(row['warmup_blocks'], row['warmup'], expected['warmup_days'])
            for name in ('baseline','policy'):
                observed = arm(row, name)
                h = job['scenario']['h0'] if name == 'baseline' else job['scenario']['h1']
                if observed['h'] != h:
                    raise ValueError('Observed arm completion differs from declared policy')
                validate_blocks(observed['observation_blocks'], observed['observation'], expected['observation_days'])
                if observed['final_snapshot']['total_ticks'] != expected['warmup_ticks']+expected['observation_ticks']:
                    raise ValueError('Final tick differs from declared horizon')
            seen.add(job_id)
            # Warmup blocks/final snapshots have been validated; retain only
            # observation data used by analysis rather than all parsed records.
            observations = {name:({'shared_with':'baseline'} if 'shared_with' in row['arms'][name]
                                  else {key:row['arms'][name][key] for key in ('observation','observation_blocks','h')})
                            for name in ('baseline','policy')}
            records[(job['scenario']['id'], job['seed'])] = dict(seed=job['seed'],scenario=job['scenario'],arms=observations)
        receipts.append(dict(manifest_sha256=digest(manifest_path), identity=identity, jobs=len(jobs),
                             software=manifest.get('software', manifest.get('versions')), sources=expected_sources,
                             completed_records_sha256=completed_hashes))
    return records, receipts


def condition_values(rows):
    """Seed rows: both-arm count/exposure plus first/last 28 control days."""
    values = []
    for row in rows:
        line = [arm(row,name)['observation'][key] for name in ('baseline','policy')
                for key in ('acquired_colonizations','patient_days','susceptible_patient_days')]
        for start,end in ((0,28),(28,56)):
            blocks = [block for block in arm(row,'baseline')['observation_blocks']
                      if start <= block['day_start'] < end]
            if len(blocks) != 28:
                raise ValueError('First/last 28-day diagnostic requires all daily blocks')
            line += [math.fsum(block[key] for block in blocks) for key in ('acquired_colonizations','patient_days')]
        values.append(line)
    return np.asarray(values, dtype=float)


def metrics(totals, n):
    # totals shape (..., conditions, 10), results (..., conditions, 12).
    c0,p0,s0,c1,p1,s1,ce,pe,cl,pl = np.moveaxis(totals,-1,0)
    r0,r1 = 1000*divide(c0,p0), 1000*divide(c1,p1)
    early,late = 1000*divide(ce,pe), 1000*divide(cl,pl)
    return np.stack((r0,r1,1-divide(r1,r0),1-divide(c1,c0),(c0-c1)/n,r0-r1,
                     1000*divide(c0,s0),1000*divide(c1,s1),early,late,late-early,divide(late,early)),axis=-1)


def bootstrap(values, draws=BOOTSTRAP_DRAWS, seed=BOOTSTRAP_SEED, chunk_size=256):
    """values shape (conditions,seeds,10); no draw×seed×condition tensor."""
    if type(draws) is not int or draws < 1 or type(chunk_size) is not int or chunk_size < 1:
        raise ValueError('Positive integer draws/chunk_size required')
    n_conditions,n,n_values = values.shape
    if n == 0 or n_values != 10 or not np.isfinite(values).all():
        raise ValueError('Bootstrap requires finite nonempty complete seed vectors')
    rng = np.random.default_rng(seed)
    flat = np.transpose(values,(1,0,2)).reshape(n,-1)
    samples = np.empty((draws,n_conditions,len(METRICS)))
    probabilities = np.full(n,1/n)
    for start in range(0,draws,chunk_size):
        stop = min(start+chunk_size,draws)
        multiplicities = rng.multinomial(n,probabilities,size=stop-start)
        totals = (multiplicities @ flat).reshape(stop-start,n_conditions,n_values)
        samples[start:stop] = metrics(totals,n)
    return metrics(values.sum(axis=1),n), samples


def q2(h,rho,e):
    nu = np.array([1-h,h])
    p = np.array([[1-h*(1-rho),h*(1-rho)],[(1-h)*(1-rho),h+rho*(1-h)]])
    d = np.diag([1.,1-e])
    return float(nu @ p @ d @ p @ d @ np.ones(2))


def summarize_condition(rows, point, samples):
    scenario = rows[0]['scenario']
    totals = condition_values(rows).sum(axis=0)
    report = dict(scenario=scenario,n=len(rows),seeds=[row['seed'] for row in rows],arms={})
    for offset,name,j in ((0,'baseline',0),(3,'policy',1)):
        counters = {}
        for kind in ('acquisitions_by_route','acquisitions_by_opportunities','acquisitions_by_route_and_opportunities'):
            counter = Counter()
            for row in rows:
                counter.update(arm(row,name)['observation']['response_counts'][kind])
            counters[kind] = dict(counter)
        report['arms'][name] = dict(acquired_total=int(totals[offset]),mean_acquired=float(totals[offset]/len(rows)),
            patient_days_total=float(totals[offset+1]),susceptible_patient_days_total=float(totals[offset+2]),
            rate_per_1000_patient_days=scalar(point[j]),rate_uncertainty=uncertainty(samples[:,j]),
            rate_per_1000_susceptible_days=scalar(point[6+j]),
            susceptible_rate_uncertainty=uncertainty(samples[:,6+j]),
            route_counts=counters['acquisitions_by_route'],opportunity_counts=counters['acquisitions_by_opportunities'],
            route_opportunity_counts=counters['acquisitions_by_route_and_opportunities'])
    for name,index in (('rate_prevented_fraction',2),('count_prevented_fraction',3),('mean_prevented',4),('absolute_prevented_rate',5)):
        report[name] = scalar(point[index])
        report[name+'_uncertainty'] = uncertainty(samples[:,index])
    q0,q1 = q2(scenario['h0'],scenario['rho'],scenario['e']),q2(scenario['h1'],scenario['rho'],scenario['e'])
    prediction = None if q0 == 0 else 1-q1/q0
    report['stationary_bridge_prediction'] = dict(opportunities=2,q0=q0,q1=q1,prevented_fraction=prediction,
        scope='Fixed stationary q2 benchmark; full ICU bridges and policy transients need not satisfy its assumptions.')
    report['F_minus_stationary_q2'] = estimate(np.nan if prediction is None else point[2]-prediction,
                                              np.full(len(samples),np.nan) if prediction is None else samples[:,2]-prediction)
    report['control_window_diagnostic'] = dict(first=dict(days=[0,28],rate=estimate(point[8],samples[:,8])),
        last=dict(days=[28,56],rate=estimate(point[9],samples[:,9])),
        last_minus_first=estimate(point[10],samples[:,10]),last_over_first=estimate(point[11],samples[:,11]),
        interpretation='Finite-window diagnostic; nonsignificance does not establish stationarity.')
    return report


def contrast(point,samples,index_a,index_b,metric):
    return estimate(point[index_a,metric]-point[index_b,metric],samples[:,index_a,metric]-samples[:,index_b,metric])


def holdout_results(ids,point,samples,valid_calibration):
    index = {identifier:i for i,identifier in enumerate(ids)}
    pairs = list(itertools.combinations(('independent','persistent','imperfect'),2))
    gates,contrasts = [],[]
    for name in ('independent','persistent','imperfect'):
        i = index[name+'_target_candidate']
        result = estimate(point[i,0],samples[:,i,0])
        interval = result['uncertainty']['intervals'][format(GATE_LEVEL,'.12g')]
        gates.append(dict(name=name+'_absolute_baseline',band=[2.7,3.3],**result,
                          passed=interval is not None and result['estimate'] is not None and 2.7 <= interval[0] <= interval[1] <= 3.3,
                          residual_from_target=scalar(point[i,0]-3.)))
    for a,b in pairs:
        i,j = index[a+'_target_candidate'],index[b+'_target_candidate']
        result = estimate(divide(point[i,0],point[j,0]),divide(samples[:,i,0],samples[:,j,0]))
        interval = result['uncertainty']['intervals'][format(GATE_LEVEL,'.12g')]
        gates.append(dict(name=a+'_over_'+b+'_baseline_ratio',band=[1/1.1,1.1],**result,
                          passed=interval is not None and result['estimate'] is not None and 1/1.1 <= interval[0] <= interval[1] <= 1.1))
        contrasts.append(dict(first=a,second=b,F_difference=contrast(point,samples,i,j,2),
                              R0_difference=contrast(point,samples,i,j,0),R1_difference=contrast(point,samples,i,j,1)))
    return dict(baseline_matching=dict(numerical_gates=gates,numerical_interval_level=GATE_LEVEL,
            valid_calibration=valid_calibration,established=valid_calibration and all(g['passed'] for g in gates),
            interpretation='Matching established only if all six interval gates and calibration gate pass; no refitting or sample extension.'),
        primary_profile_F_contrasts=contrasts,primary_family_interval_level=PRIMARY_LEVEL,
        reporting_precision_goal=dict(target_F_descriptive95_halfwidth=.05,fixed_n=True,
                                      no_sample_extension_for_unmet_precision=True),
        uncertainty_scope='Fresh holdout intervals conditional on locked p; calibration selection/interpolation uncertainty and clinical-input uncertainty excluded.')


def environment_results(ids,point,samples):
    index = {identifier:i for i,identifier in enumerate(ids)}
    effects,draw_effects,point_effects = [],[],[]
    for label,p in zip(('075','100','125'),builder.OLD_GRID):
        i,j = index[f'independent_low{label}_env1'],index[f'independent_low{label}_env0']
        effects.append(dict(p=p,full_minus_disabled={METRICS[k]:contrast(point,samples,i,j,k) for k in (0,1,2,4,5)}))
        draw_effects.append(samples[:,i,2]-samples[:,j,2])
        point_effects.append(point[i,2]-point[j,2])
    return dict(primary_environment_effects=effects,primary_family_interval_level=PRIMARY_LEVEL,
        reporting_precision_goal=dict(environment_F_contrast_descriptive95_halfwidth=.10,fixed_n=True,
                                      no_sample_extension_for_unmet_precision=True),
        high_minus_low_F_effect=estimate(point_effects[2]-point_effects[0],draw_effects[2]-draw_effects[0]),
        center_minus_mean_endpoints_F_effect=estimate(point_effects[1]-(point_effects[0]+point_effects[2])/2,
                                                     draw_effects[1]-(draw_effects[0]+draw_effects[2])/2),
        interpretation='Fixed-p total environmental intervention including reservoir/patient feedback; not a causal fraction of a specific pickup pathway.')


def calibration_results(point,samples):
    if not np.isfinite(point[:,0]).all():
        return dict(status='failed_undefined_baseline_exposure',candidate_p=None)
    selected = builder.choose(zip(builder.GRID,point[:,0].tolist()))
    rates = samples[:,:,0]
    finite_rows = np.isfinite(rates).all(axis=1)
    exact = (rates == 3.) & finite_rows[:,None]
    brackets = (rates[:,:-1] < 3.) & (rates[:,1:] > 3.) & finite_rows[:,None]
    any_exact,any_bracket = exact.any(axis=1),brackets.any(axis=1)
    valid = any_exact|any_bracket
    candidates = np.full(len(samples),.15)
    for j,p in reversed(list(enumerate(builder.GRID))):
        candidates[exact[:,j]] = p
    for j in (1,0):
        use = brackets[:,j] & ~any_exact
        fraction = divide(3.-rates[use,j],rates[use,j+1]-rates[use,j])
        candidates[use] = builder.GRID[j]+(builder.GRID[j+1]-builder.GRID[j])*fraction
    conditional = uncertainty(candidates[valid],(.95,))
    conditional['conditional_on_valid_exact_hit_or_bracket'] = True
    return dict(selection=selected,candidate_p_uncertainty_conditional_valid=conditional,
        fallback_inclusive_candidate_p_uncertainty=uncertainty(candidates,(.95,)),
        draws_without_exact_hit_or_bracket=int(np.count_nonzero(~valid)),
        draws_without_increasing_bracket=int(np.count_nonzero(~any_bracket)),
        draws_with_undefined_baseline_rate=int(np.count_nonzero(~finite_rows)),
        nonmonotone_bootstrap_draws=int(np.count_nonzero(np.any(np.diff(rates,axis=1)<0,axis=1)&finite_rows)),
        scope='Candidate-p selection uncertainty published separately; no policy results enter calibration.')


def bank_diagnostics(values,seeds,ids,book,module,draws,seed,chunk):
    banks = {}
    sampled = {}
    for ordinal,name in enumerate(('A','B'),start=1):
        selected = book['seed_blocks'][module][name]
        mask = [i for i,s in enumerate(seeds) if s in set(selected)]
        if len(mask) != len(selected):
            raise ValueError('Missing bank cohort members')
        # Separate addressable RNG streams ensure independent bank resampling.
        point,samples = bootstrap(values[:,mask,:],draws,np.random.SeedSequence([seed,ordinal]),chunk)
        sampled[name] = (point,samples[:,:,:3].copy())
        banks[name] = dict(n=len(mask),seeds=[seeds[i] for i in mask],
            conditions={identifier:{METRICS[k]:estimate(point[i,k],samples[:,i,k],(.95,)) for k in (0,1,2)} for i,identifier in enumerate(ids)})
        if module == 'holdout':
            banks[name]['profile_F_contrasts'] = holdout_results(ids,point,samples,False)['primary_profile_F_contrasts']
        elif module == 'environment':
            banks[name]['environment_effects'] = environment_results(ids,point,samples)['primary_environment_effects']
        del samples
    pa,sa = sampled['A']; pb,sb = sampled['B']
    differences = {identifier:{METRICS[k]:estimate(pb[i,k]-pa[i,k],sb[:,i,k]-sa[:,i,k],(.95,)) for k in (0,1,2)} for i,identifier in enumerate(ids)}
    report = dict(banks=banks,B_minus_A=differences,
                  method='Independent whole-seed resampling within A and B; synchronized condition weights within each bank. Primary estimates pool both banks regardless of disagreement.')
    if module == 'environment':
        index = {identifier:i for i,identifier in enumerate(ids)}
        interactions = []
        for label,p in zip(('075','100','125'),builder.OLD_GRID):
            i,j = index[f'independent_low{label}_env1'],index[f'independent_low{label}_env0']
            interactions.append(dict(p=p,B_minus_A_F_ablation=estimate((pb[i,2]-pb[j,2])-(pa[i,2]-pa[j,2]),
                (sb[:,i,2]-sb[:,j,2])-(sa[:,i,2]-sa[:,j,2]),(.95,))))
        report['B_minus_A_environment_F_effects'] = interactions
    return report


def analyze_records(records,expected_config,book,module,draws=BOOTSTRAP_DRAWS,seed=BOOTSTRAP_SEED,chunk_size=256,selection=None):
    config = normalized(expected_config)
    seeds = config['seeds']
    ids = [s['id'] for s in config['scenarios']]
    expected_keys = {(identifier,s) for identifier in ids for s in seeds}
    if set(records) != expected_keys:
        raise ValueError('Analysis requires complete common condition seed coverage')
    groups = [[records[(identifier,s)] for s in seeds] for identifier in ids]
    values = np.asarray([condition_values(rows) for rows in groups])
    point,samples = bootstrap(values,draws,seed,chunk_size)
    report = dict(schema_version=1,module=module,status='complete',n_paired_jobs=len(records),
        bootstrap_draws=draws,bootstrap_seed=seed,bootstrap_chunk_size=chunk_size,
        bootstrap_method='Whole-seed multinomial multiplicities, common across all condition arms; pooled rates are ratios of summed counts/exposure.',
        uncertainty_scope='Monte Carlo uncertainty in this synthetic model; bootstrap coverage approximate.',
        scenarios=[summarize_condition(rows,point[i],samples[:,i]) for i,rows in enumerate(groups)],
        descriptive_interval_level=.95)
    if module == 'calibration':
        report['calibration'] = calibration_results(point,samples)
    elif module == 'holdout':
        if selection is None:
            raise ValueError('Holdout analysis requires calibration selection status')
        report.update(holdout_results(ids,point,samples,bool(selection['valid_target_calibration'])))
    else:
        report.update(environment_results(ids,point,samples))
    del samples
    report['two_bank_diagnostics'] = bank_diagnostics(values,seeds,ids,book,module,draws,seed,chunk_size)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--module',choices=('calibration','holdout','environment'),required=True)
    parser.add_argument('--bank',type=Path,action='append',required=True)
    parser.add_argument('--design-config',type=Path,required=True)
    parser.add_argument('--seedbook',type=Path,default=HERE/'seedbook.json')
    parser.add_argument('--proposal',type=Path,default=HERE/'candidate_design.json')
    parser.add_argument('--selection',type=Path)
    parser.add_argument('--calibration-report',type=Path)
    parser.add_argument('--freeze',type=Path,required=True,help='New-study freeze bound by each confirmation_binding.json.')
    parser.add_argument('--analysis-freeze',type=Path,default=HERE/'analyze_freeze.json',help='Separate pre-outcome analyzer source freeze.')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--chunk-size',type=int,default=256)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise FileExistsError(args.output)
    analysis_freeze = verify_analysis_freeze(args.analysis_freeze)
    config = json.loads(args.design_config.read_text())
    book = json.loads(args.seedbook.read_text())
    proposal = json.loads(args.proposal.read_text())
    selection = None if args.selection is None else json.loads(args.selection.read_text())
    calibration = None if args.calibration_report is None else json.loads(args.calibration_report.read_text())
    validate_design(args.module,config,book,proposal,selection,calibration,args.selection,args.calibration_report,args.proposal)
    wrapper = load_module('_confirmation_wrapper_for_validation',HERE/'array_worker.py')
    freeze = wrapper.verify_freeze(args.freeze)
    # Require the analyzer and all declaration inputs themselves in this freeze.
    for path in (args.design_config,args.seedbook,args.proposal,HERE/'output_scopes.json'):
        name = str(path.resolve().relative_to(ABM))
        if freeze['files'].get(name) != digest(path):
            raise ValueError('Analysis input is absent from freeze: '+name)
    if args.module == 'holdout':
        for path in (args.selection,args.calibration_report):
            name = str(path.resolve().relative_to(ABM))
            if freeze['files'].get(name) != digest(path):
                raise ValueError('Holdout selection evidence is absent from freeze: '+name)
    sources = backend.source_hashes()
    sources['array_worker.py'] = digest(RESPONSE/'array_worker.py')
    for name,sha in sources.items():
        relative = str((RESPONSE/'array_worker.py').relative_to(ABM)) if name == 'array_worker.py' else name
        if freeze['files'].get(relative) != sha:
            raise ValueError('Scientific backend source is absent from freeze: '+relative)
    context = dict(freeze_path=args.freeze,freeze=freeze,seedbook_path=args.seedbook,module=args.module,
                   output_scopes=json.loads((HERE/'output_scopes.json').read_text()))
    records,receipts = load_banks(args.bank,config,sources,book,args.module,context)
    report = analyze_records(records,config,book,args.module,chunk_size=args.chunk_size,selection=selection)
    if verify_analysis_freeze(args.analysis_freeze) != analysis_freeze:
        raise ValueError('Analysis freeze changed during report generation')
    report.update(banks=receipts,source_sha256=digest(Path(__file__)),design_sha256=digest(args.design_config),
                  seedbook_sha256=digest(args.seedbook),proposal_sha256=digest(args.proposal),
                  freeze_sha256=digest(args.freeze),freeze_identity=freeze['identity_sha256'])
    report.update(analysis_freeze_sha256=digest(args.analysis_freeze),analysis_freeze_identity=analysis_freeze['identity_sha256'],
                  analysis_software=dict(python=sys.version,numpy=np.__version__))
    if args.selection is not None:
        report['calibration_selection_sha256'] = digest(args.selection)
        report['calibration_report_sha256'] = digest(args.calibration_report)
    args.output.mkdir(parents=True)
    (args.output/'report.json').write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
    columns = ['scenario','p','rho','e','h0','h1','n','baseline_rate','policy_rate',
               'count_prevented_fraction','rate_prevented_fraction','mean_prevented','absolute_prevented_rate','prediction']
    with (args.output/'summary.csv').open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=columns)
        writer.writeheader()
        for row in report['scenarios']:
            scenario = row['scenario']
            writer.writerow(dict(scenario=scenario['id'],**{key:scenario[key] for key in ('p','rho','e','h0','h1')},n=row['n'],
                baseline_rate=row['arms']['baseline']['rate_per_1000_patient_days'],policy_rate=row['arms']['policy']['rate_per_1000_patient_days'],
                **{key:row[key] for key in ('count_prevented_fraction','rate_prevented_fraction','mean_prevented','absolute_prevented_rate')},
                prediction=row['stationary_bridge_prediction']['prevented_fraction']))
    print(canonical(dict(module=args.module,complete_jobs=len(records),bootstrap_draws=BOOTSTRAP_DRAWS,
                         baseline_matching=report.get('baseline_matching',{}).get('established'))))
    return 0


if __name__ == '__main__':
    sys.exit(main())
