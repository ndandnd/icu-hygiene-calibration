#!/usr/bin/env python3
"""Frozen descriptive follow-ups from saved confirmation records; no simulation."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
STUDY = 'research/full_icu_confirmation_2026-10-08'
RESULTS = 'results/full_icu_confirmation_2026-10-08'
MODULES = ('holdout', 'environment', 'calibration')


def read(path):
    return json.loads(Path(path).read_text())


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def write(path, obj):
    with Path(path).open('x') as stream:
        stream.write(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False)+'\n')


def uncertainty(samples):
    samples = np.asarray(samples)
    undefined = int(np.count_nonzero(~np.isfinite(samples)))
    return dict(interval95=None if undefined or not len(samples) else np.quantile(samples, [.025, .975]).tolist(),
                bootstrap_draws=len(samples), undefined_bootstrap_draws=undefined,
                interval_unavailable_if_any_draw_undefined=True, conditional_on_defined_draws=False)


def estimate(point, samples):
    return dict(estimate=float(point) if np.isfinite(point) else None, uncertainty=uncertainty(samples))


def pooled_bootstrap(values, config, task, cohort=0):
    """Seed first; synchronized weights across all other dimensions."""
    values = np.asarray(values, dtype=float)
    if len(values) == 0 or not np.isfinite(values).all():
        raise ValueError('Finite complete seed vectors required')
    n = len(values)
    flat = values.reshape(n, -1)
    rng = np.random.default_rng(np.random.SeedSequence([config['bootstrap_seed'], config['stream_ordinals'][task], cohort]))
    samples = np.empty((config['bootstrap_draws'], *values.shape[1:]))
    for start in range(0, config['bootstrap_draws'], config['chunk_size']):
        end = min(start+config['chunk_size'], config['bootstrap_draws'])
        weights = rng.multinomial(n, np.full(n, 1/n), size=end-start)
        samples[start:end] = (weights @ flat).reshape(end-start, *values.shape[1:])
    return values.sum(axis=0), samples


def rate(totals):
    counts, days = totals[..., 0], totals[..., 1]
    return np.divide(1000*counts, days, out=np.full(counts.shape, np.nan), where=days > 0)


def log_slope(outcomes, p):
    y = np.asarray(outcomes, dtype=float)
    valid = np.isfinite(y).all(axis=-1) & (y > 0).all(axis=-1)
    x = np.log(p); centered = x-x.mean()
    out = np.full(y.shape[:-1], np.nan)
    out[valid] = np.log(y[valid]) @ centered / (centered @ centered)
    return out


def select(p, rates, target=3., fallback=.15):
    rates = np.atleast_2d(rates)
    finite = np.isfinite(rates).all(axis=1)
    exact = (rates == target) & finite[:, None]
    bracket = (rates[:, :-1] < target) & (rates[:, 1:] > target) & finite[:, None]
    valid = exact.any(axis=1) | bracket.any(axis=1)
    candidate = np.full(len(rates), fallback)
    # Reverse iteration makes lowest-p selected when more than one is eligible.
    for j in reversed(range(len(p)-1)):
        use = bracket[:, j] & ~exact.any(axis=1)
        candidate[use] = p[j]+(target-rates[use, j])*(p[j+1]-p[j])/(rates[use, j+1]-rates[use, j])
    for j in reversed(range(len(p))):
        candidate[exact[:, j]] = p[j]
    return candidate, valid, exact.any(axis=1), bracket.any(axis=1), finite


def daily_values(row, config):
    blocks = row['warmup_blocks']+row['arms']['baseline']['observation_blocks']
    horizons = ((row['warmup_blocks'], row['warmup'], config['warmup_days']),
                (row['arms']['baseline']['observation_blocks'], row['arms']['baseline']['observation'], config['observation_days']))
    for phase, summary, days in horizons:
        if len(phase) != days:
            raise ValueError('Incomplete daily horizon')
        for i, block in enumerate(phase):
            if (block['day_start'], block['day_end'], block['days']) != (i, i+1, 1):
                raise ValueError('Noncontiguous daily coverage')
        for key in ('acquired_colonizations', 'patient_days'):
            if not math.isclose(math.fsum(b[key] for b in phase), summary[key], abs_tol=1e-7, rel_tol=1e-10):
                raise ValueError('Daily totals fail conservation')
    pairs = np.asarray([[b['acquired_colonizations'], b['patient_days']] for b in blocks])
    if not np.isfinite(pairs).all() or (pairs < 0).any():
        raise ValueError('Invalid count/exposure')
    return np.asarray([pairs[start:start+config['week_length_days']].sum(axis=0)
                       for start in range(0, len(blocks), config['week_length_days'])])


def closure():
    """Use prior independently reviewed report receipts to bind every record."""
    hashes = {}; modules = {}
    def add(relative, expected=None):
        relative = str(relative)
        actual = sha(ROOT/relative)
        if expected is not None and actual != expected:
            raise ValueError('Frozen evidence hash mismatch: '+relative)
        if relative in hashes and hashes[relative] != actual:
            raise ValueError('Changing evidence')
        hashes[relative] = actual
    scopes = read(ROOT/STUDY/'output_scopes.json')
    book = read(ROOT/STUDY/'seedbook.json')
    add(STUDY+'/output_scopes.json'); add(STUDY+'/seedbook.json')
    for module in MODULES:
        report_path = f'{STUDY}/{module}_report/report.json'
        add(report_path)
        report = read(ROOT/report_path)
        manifests = {}
        for site in ('local', 'cluster'):
            directory = scopes['banks'][f'{module}_{site}']
            for name in ('manifest.json', 'confirmation_binding.json'):
                add(directory+'/'+name)
            manifests[sha(ROOT/directory/'manifest.json')] = directory
        paths = []
        for bank in report['banks']:
            directory = manifests[bank['manifest_sha256']]
            for name, expected in sorted(bank['completed_records_sha256'].items()):
                relative = directory+'/completed/'+name
                add(relative, expected); paths.append(relative)
        expected_n = {'holdout': 1536, 'environment': 3072, 'calibration': 384}[module]
        if len(paths) != expected_n or len(paths) != len(set(paths)):
            raise ValueError('Incorrect primary record closure')
        modules[module] = dict(records=paths, report=report_path, seeds=report['scenarios'][0]['seeds'],
                               conditions=[s['scenario'] for s in report['scenarios']], banks=book['seed_blocks'][module])
        for name in ('all', 'local', 'cluster'):
            add(f'{STUDY}/configs/{module}_{name}.json')
    bridge_report = f'{STUDY}/bridge/aggregate_report/summary.json'
    add(bridge_report); add(f'{STUDY}/bridge/aggregate_report/binding.json')
    report = read(ROOT/bridge_report); paths = []
    for shard in report['banks']:
        directory = shard['output_scope']
        for name, key in (('manifest.json','manifest_sha256'), ('summary.json','summary_sha256'), ('bridge_binding.json','bridge_binding_sha256')):
            add(directory+'/'+name, shard[key])
        add(shard['config_path'], shard['config_sha256'])
        for name, expected in sorted(shard['completed_records_sha256'].items()):
            relative = directory+'/completed/'+name; add(relative, expected); paths.append(relative)
    if len(paths) != 256 or len(set(paths)) != 256:
        raise ValueError('Incorrect bridge record closure')
    modules['bridge'] = dict(records=paths, report=bridge_report, seeds=report['seeds'])
    for relative in ('freeze.json', 'freeze_holdout.json', 'freeze_bridge.json', 'analyze_freeze.json',
                     'build_design.py', 'analyze.py', 'PROTOCOL.md', 'calibration_selection.json', 'holdout_candidates.json',
                     'bridge/seedbook.json', 'bridge/PROTOCOL.md', 'bridge/aggregate.py', 'bridge/coefficient.py', 'bridge/background.py'):
        add(STUDY+'/'+relative)
    for name in ('baseline_report.py', 'freeze_validation.py'):
        add('research/full_icu_calibration_2026-10-07/'+name)
    return dict(schema_version=1, files=hashes, modules=modules)


def prepare():
    manifest = closure()
    write(HERE/'input_manifest.json', manifest)
    files = {name: sha(HERE/name) for name in ('METHODS.md', 'config.json', 'analyze_saved.py', 'test_analyze_saved.py', 'input_manifest.json')}
    payload = dict(schema_version=1, scope='prospective_derived_saved_data_analysis', files=files,
                   created_utc=datetime.now(timezone.utc).isoformat(), old_outcomes_already_available=True)
    payload['identity_sha256'] = hashlib.sha256(canonical(payload).encode()).hexdigest()
    write(HERE/'freeze.json', payload)
    print(canonical(dict(status='frozen', files=len(manifest['files']), records={k:len(v['records']) for k,v in manifest['modules'].items()}, identity=payload['identity_sha256'])))


def verify_freeze():
    frozen = read(HERE/'freeze.json'); payload = {k:v for k,v in frozen.items() if k != 'identity_sha256'}
    if hashlib.sha256(canonical(payload).encode()).hexdigest() != frozen['identity_sha256']:
        raise ValueError('Freeze identity mismatch')
    for name, expected in frozen['files'].items():
        if sha(HERE/name) != expected:
            raise ValueError('Changed analysis source/method/config: '+name)
    return frozen, read(HERE/'input_manifest.json')


def load_records(module, manifest):
    spec = manifest['modules'][module]
    for relative in spec['records']:
        if sha(ROOT/relative) != manifest['files'][relative]:
            raise ValueError('Changed completed record: '+relative)
        row = read(ROOT/relative)
        if row['status'] != 'complete' or row['seed'] not in spec['seeds']:
            raise ValueError('Invalid completion identity')
        yield row


def indexed_records(module, manifest):
    spec = manifest['modules'][module]; rows = {}
    for row in load_records(module, manifest):
        key = (row['scenario']['id'], row['seed'])
        if key in rows:
            raise ValueError('Duplicate condition/seed')
        rows[key] = row
    ids = [s['id'] for s in spec['conditions']]
    if set(rows) != {(identifier, seed) for identifier in ids for seed in spec['seeds']}:
        raise ValueError('Incomplete synchronized seed closure')
    # Reconcile with frozen pooled baseline totals in prior report.
    report = read(ROOT/spec['report'])
    for scenario in report['scenarios']:
        values = [rows[(scenario['scenario']['id'], seed)]['arms']['baseline']['observation'] for seed in spec['seeds']]
        for key, prior in (('acquired_colonizations','acquired_total'), ('patient_days','patient_days_total')):
            if not math.isclose(math.fsum(v[key] for v in values), scenario['arms']['baseline'][prior], rel_tol=1e-10, abs_tol=1e-7):
                raise ValueError('Saved records fail reviewed report conservation')
    return spec, rows


def cohorts(spec):
    yield 'pooled', 0, spec['seeds']
    for ordinal, name in enumerate(('A','B'), 1):
        yield name, ordinal, spec['banks'][name]


def weekly(manifest, config):
    spec, rows = indexed_records('holdout', manifest); ids = [s['id'] for s in spec['conditions']]
    values = np.asarray([[daily_values(rows[(identifier, seed)], config) for identifier in ids] for seed in spec['seeds']])
    output = []; lookup = {s:i for i,s in enumerate(spec['seeds'])}
    for name, ordinal, seeds in cohorts(spec):
        point, samples = pooled_bootstrap(values[[lookup[s] for s in seeds]], config, 'weekly', ordinal)
        pr, sr = rate(point), rate(samples)
        for i, identifier in enumerate(ids):
            for j, start in enumerate(range(0, 146, 7)):
                end = min(start+7, 146)
                phase = 'warmup' if end <= 90 else 'observation' if start >= 90 else 'mixed'
                output.append(dict(cohort=name, n_seeds=len(seeds), profile=identifier, week=j+1, day_start=start, day_end=end,
                                   days=end-start, phase=phase, partial_week=end-start < 7,
                                   acquired_total=int(point[i,j,0]), patient_days_total=float(point[i,j,1]),
                                   rate_per_1000_patient_days=estimate(pr[i,j], sr[:,i,j])))
    return dict(rows=output, interpretation='Finite-horizon control trajectory; does not establish equilibrium or stationarity.')


def slopes(manifest, config):
    spec, rows = indexed_records('environment', manifest)
    scenarios = spec['conditions']; ids = [s['id'] for s in scenarios]
    p = sorted({s['p'] for s in scenarios})
    full = [next(i for i,s in enumerate(scenarios) if s['p']==v and s['environment_multiplier']==1) for v in p]
    off = [next(i for i,s in enumerate(scenarios) if s['p']==v and s['environment_multiplier']==0) for v in p]
    values = np.asarray([[[rows[(identifier,seed)]['arms']['baseline']['observation'][key] for key in ('acquired_colonizations','patient_days')] for identifier in ids] for seed in spec['seeds']])
    output = {}; lookup = {s:i for i,s in enumerate(spec['seeds'])}
    for name, ordinal, seeds in cohorts(spec):
        point, samples = pooled_bootstrap(values[[lookup[s] for s in seeds]], config, 'slopes', ordinal)
        pr, sr = rate(point), rate(samples)
        difference, draw_difference = pr[full]-pr[off], sr[:,full]-sr[:,off]
        output[name] = dict(n_seeds=len(seeds), p=p, env_off_rates=pr[off].tolist(), full_minus_off_rates=difference.tolist(),
                            env_off_log_slope=estimate(log_slope(pr[off],p), log_slope(sr[:,off],p)),
                            full_minus_off_log_slope=estimate(log_slope(difference,p), log_slope(draw_difference,p)))
    return dict(cohorts=output, interpretation='Exploratory finite-p log slopes; not evidence of an asymptotic mechanism. Environmental ablation includes feedback and is not route counts.')


def calibration(manifest, config):
    spec, rows = indexed_records('calibration', manifest)
    ids = [s['id'] for s in spec['conditions']]; p = [s['p'] for s in spec['conditions']]
    values = np.asarray([[[rows[(identifier,seed)]['arms']['baseline']['observation'][key] for key in ('acquired_colonizations','patient_days')] for identifier in ids] for seed in spec['seeds']])
    output = {}; lookup = {s:i for i,s in enumerate(spec['seeds'])}
    for name, ordinal, seeds in cohorts(spec):
        point, samples = pooled_bootstrap(values[[lookup[s] for s in seeds]], config, 'calibration', ordinal)
        pr, sr = rate(point), rate(samples)
        selected, valid, exact, bracket, finite = select(p,sr,config['calibration_target_rate'], config['calibration_fallback_p'])
        point_p, point_valid, *_ = select(p,pr)
        conditional = uncertainty(selected[valid]); conditional['conditional_on_valid_exact_hit_or_bracket']=True
        output[name] = dict(n_seeds=len(seeds), p_grid=p, point_baseline_rates=pr.tolist(), selected_p=float(point_p[0]),
                            valid_point_bracket_or_hit=bool(point_valid[0]), selection_method='lowest_p_exact_hit_else_lowest_p_increasing_adjacent_bracket_else_fixed_fallback',
                            conditional_selected_p_uncertainty=conditional, fallback_inclusive_selected_p_uncertainty=uncertainty(selected),
                            exact_hit_draws=int(exact.sum()), no_increasing_bracket_draws=int((~bracket).sum()),
                            no_valid_hit_or_bracket_draws=int((~valid).sum()), no_valid_hit_or_bracket_frequency=float((~valid).mean()),
                            undefined_baseline_rate_draws=int((~finite).sum()), nonmonotone_draws=int((np.any(np.diff(sr,axis=1)<0,axis=1)&finite).sum()))
    return dict(cohorts=output, scope='Selection uncertainty for iid p from128 seeds per3p, separate from conditional holdout validation. Other2 profile p choices originate in older study; their uncertainty is not estimated here.',
                interpolation_verification=dict(confirmation='linear p versus baseline rate', prior_validation='linear p versus acquisition count', prior_refinement='linear p versus acquisition count', logarithmic=False))


def bridge(manifest, config):
    spec = manifest['modules']['bridge']; rows = list(load_records('bridge',manifest))
    if {r['seed'] for r in rows} != set(spec['seeds']) or len(rows)!=256:
        raise ValueError('Bridge seed membership differs')
    hist = {label:0 for label in ('0','1','2','>=3')}; raw = 0; observed = 0; warmup = 0; days = 0.
    for row in rows:
        arms = [arm for profile in row['profiles'].values() for arm in profile.values()]
        first = arms[0]
        for arm in arms:
            if (arm['raw_pairs_without_reset'], arm['raw_pairs_without_reset_by_opportunities'], arm['event_counts']) != (first['raw_pairs_without_reset'], first['raw_pairs_without_reset_by_opportunities'], first['event_counts']):
                raise ValueError('Raw counts differ across profiles/arms')
            if sum(arm['raw_pairs_without_reset_by_opportunities'].values()) != arm['raw_pairs_without_reset']:
                raise ValueError('Raw histogram fails conservation')
        for label in hist: hist[label] += first['raw_pairs_without_reset_by_opportunities'][label]
        raw += first['raw_pairs_without_reset']; observed += first['event_counts']['opportunities_observed']; warmup += first['event_counts']['opportunities_warmup']; days += row['observation_patient_days']
        if row['event_counts']['opportunity'] != first['event_counts']['opportunities_observed']+first['event_counts']['opportunities_warmup']:
            raise ValueError('Opportunity periods fail conservation')
    return dict(n_backgrounds=256, raw_unweighted_pairs_total=raw, raw_unweighted_pairs_histogram=hist,
                raw_histogram_fractions={k:v/raw for k,v in hist.items()}, observation_days=56, warmup_days=90,
                opportunities_observed_total=observed, opportunities_warmup_total=warmup, observation_patient_days_total=days,
                observed_opportunities_per_patient_day=observed/days, warmup_opportunities_per_patient_day=None,
                opportunity_period_scope='Observation numerator and denominator both span56 days. Warmup90-day numerator has no saved patient-day denominator.',
                pair_period_scope='Observed susceptible depositions define pair endpoints; eligible pickups and intervening opportunities can originate in warmup.',
                patient_contact_limitation='Saved whole146-day patient event counts and weighted K_by_worker do not identify per-role patient contact counts; bridge replay needed.',
                interpretation='Fixed-sample observed p=0 raw totals, unweighted by hygiene; no inference about positive-p route shares.')


def flatten_csv(path, rows):
    flat=[]
    for row in rows:
        item = dict(row); result = item.pop('rate_per_1000_patient_days')
        interval=result['uncertainty']['interval95'];item.update(rate_per_1000_patient_days=result['estimate'], interval95_low=None if interval is None else interval[0], interval95_high=None if interval is None else interval[1], undefined_bootstrap_draws=result['uncertainty']['undefined_bootstrap_draws'])
        flat.append(item)
    with path.open('x', newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)


def run(task):
    frozen, manifest = verify_freeze(); config = read(HERE/'config.json')
    output=HERE/'outputs'; output.mkdir(exist_ok=True)
    target = output/(task+'.json')
    if target.exists(): raise FileExistsError(target)
    module={'weekly':'holdout','slopes':'environment','calibration':'calibration','bridge':'bridge'}[task]
    # Validate all non-record inputs once before computation; selected record hashes are checked during loading.
    for relative, expected in manifest['files'].items():
        if '/completed/' not in relative and sha(ROOT/relative)!=expected:
            raise ValueError('Changed input declaration/report/source: '+relative)
    result=globals()[task](manifest,config)
    result.update(schema_version=1, task=task, bootstrap_draws=None if task=='bridge' else config['bootstrap_draws'], bootstrap_seed=config['bootstrap_seed'],
                  descriptive_interval_level=.95, uncertainty_scope='Synthetic-model Monte Carlo uncertainty; excludes clinical-input uncertainty.', analysis_freeze_identity=frozen['identity_sha256'])
    verify_freeze();write(target,result)
    if task=='weekly':flatten_csv(output/'weekly.csv',result['rows'])
    if task=='slopes':
        with (output/'slopes.csv').open('x',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['cohort','n_seeds','outcome','slope','interval95_low','interval95_high','undefined_draws'])
            for name, cohort in result['cohorts'].items():
                for outcome in ('env_off_log_slope','full_minus_off_log_slope'):
                    est=cohort[outcome]; interval=est['uncertainty']['interval95']
                    writer.writerow([name,cohort['n_seeds'],outcome,est['estimate'],None if interval is None else interval[0],None if interval is None else interval[1],est['uncertainty']['undefined_bootstrap_draws']])
    if task=='calibration':
        with (output/'calibration.csv').open('x',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['cohort','n_seeds','selected_p','no_valid_hit_or_bracket_draws','no_valid_hit_or_bracket_frequency','conditional95_low','conditional95_high','fallback_inclusive95_low','fallback_inclusive95_high'])
            for name,cohort in result['cohorts'].items():
                conditional=cohort['conditional_selected_p_uncertainty']['interval95']; fallback=cohort['fallback_inclusive_selected_p_uncertainty']['interval95']
                writer.writerow([name,cohort['n_seeds'],cohort['selected_p'],cohort['no_valid_hit_or_bracket_draws'],cohort['no_valid_hit_or_bracket_frequency'],None if conditional is None else conditional[0],None if conditional is None else conditional[1],None if fallback is None else fallback[0],None if fallback is None else fallback[1]])
    if task=='bridge':
        with (output/'bridge_histogram.csv').open('x',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['opportunity_bin','raw_unweighted_pairs','fraction']);writer.writerows((k,v,result['raw_histogram_fractions'][k]) for k,v in result['raw_unweighted_pairs_histogram'].items())
    inputs={k:v for k,v in manifest['files'].items() if '/completed/' not in k or k in set(manifest['modules'][module]['records'])}
    receipt=dict(schema_version=1,task=task,status='complete',analysis_freeze_identity=frozen['identity_sha256'], analysis_freeze_sha256=sha(HERE/'freeze.json'),
                 config_sha256=sha(HERE/'config.json'),source_sha256=sha(Path(__file__)),input_hashes=inputs,
                 output_hashes={p.name:sha(p) for p in output.glob(task+'*')}, software=dict(python=sys.version,numpy=np.__version__), created_utc=datetime.now(timezone.utc).isoformat())
    write(output/(task+'_receipt.json'),receipt)
    print(canonical(dict(task=task,status='complete',output=str(target.relative_to(ROOT)),analysis_freeze_identity=frozen['identity_sha256'])))


def main():
    global ROOT
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--abm-root',type=Path,default=ROOT);parser.add_argument('action',choices=('prepare','run'));parser.add_argument('--task',choices=('weekly','slopes','calibration','bridge','all'),default='all');args=parser.parse_args()
    ROOT=args.abm_root.resolve()
    if args.action=='prepare':prepare()
    else:
        for task in ('weekly','slopes','calibration','bridge') if args.task=='all' else (args.task,):run(task)


if __name__ == '__main__':main()
