#!/usr/bin/env python3
"""Protocol outputs of workflow_coefficient_2026-10-21 (summary JSON).

Read-only inputs:
  alternative workflow   <records>/<seed>.json + <records>/manifest.json, each
                         record validated by run_batch.validate_record
                         (payload digest, manifest identity, seed, surface).
  original workflow,     results/full_icu_confirmation_2026-10-08/bridge_{local,
    3 primary profiles   cluster}/part_*/completed/<seed>.json (256 archived).
  original workflow,     the archived extended run behind supplement Table S5
    22-row surface       (256 records, payload digests verified).
  post hoc finite-p      robustness execution01 analysis.json/selection.json
                         (hash-pinned; used only in the labelled post hoc block).

Statistics (PROTOCOL.md "Paired comparison"):
  * F* = 1 - mean(K1)/mean(K0); never a mean of seed ratios; zero seeds kept.
  * Nominal intervals exactly as the original study (frozen run.summarize):
    mean +- 1.96 SE; delta SE of F* = sd((K1_i - r K0_i)/mean K0)/sqrt(n).
    Paired workflow difference dF* uses the analogous influence function
    psi_i = -(y_a - r_a x_a)/mean x_a + (y_o - r_o x_o)/mean x_o.
  * Seed-level paired bootstrap: whole seeds, 20,000 draws, one index matrix
    numpy.random.default_rng(20261021).integers(0, n, (20000, n)) per seed set
    (seeds in seedbook order), used synchronously for every statistic on that
    set, so both workflows always see the same resampled seeds. Percentile 95%
    intervals (numpy linear interpolation); undefined draws (mean K0 = 0, or no
    crossing) are counted and excluded.
  * Pairing is on the seeds present in both sources. If the fallback plan ran,
    the surface is summarized on the bank-A seeds that carry it.

Self-checks against the archive run every time (original F*, D0, aggregate
summary, Table S5 parsed from supplement.tex, frozen-summarizer agreement).
A failed self-check is reported in the output and gives exit code 1.

Usage:
  python -B code/summarize.py                      # records/ -> summary.json (or progress_summary.json if incomplete)
  python -B code/summarize.py --records gate_checks/benchmark/records --output gate_checks/benchmark/summary_test.json
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402
from pathlib import Path  # noqa: E402
import re  # noqa: E402
import statistics  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
import workflow_background as wb  # noqa: E402
import run_batch  # noqa: E402

STUDY = wb.STUDY
ABM = wb.ABM
ORIGINAL_RESULTS = ABM / 'results/full_icu_confirmation_2026-10-08'
ORIGINAL_AGGREGATE = wb.BRIDGE / 'aggregate_report/summary.json'
ORIGINAL_AGGREGATE_SHA256 = '46318f3aa337c489cb9e977b02f7d9e582a76069e24506d75a68677ad7abac4a'
SUPPLEMENT = ABM / 'manuscript/supplement.tex'
ROBUSTNESS = ABM.parent / 'research/robustness_uq_2026-10-17/collection_2026-10-21/execution01'
FINITE_P_ANALYSIS = ROBUSTNESS / 'outputs/analysis.json'
FINITE_P_ANALYSIS_SHA256 = '063c0d95f37661de931d0b433168cdb0e178affe41b5017e7c5eb49aaec5579f'
FINITE_P_SELECTION = ROBUSTNESS / 'source/research/robustness_uq_2026-10-17/selection.json'
FINITE_P_SELECTION_SHA256 = '1dcb9eac1f84ac24ab5fbee814fac4711e64eeb8d12beb1b653255d08ed96ecf'

BOOTSTRAP_DRAWS = 20000
BOOTSTRAP_SEED = 20261021
TICKS_PER_DAY = 2592
PRIMARY = (('iid_perfect', 'h0.9_to_0.95_rho0_e1'),
           ('persistent_perfect', 'h0.9_to_0.95_rho0.5_e1'),
           ('iid_imperfect', 'h0.9_to_0.95_rho0_e0.8'))
# PROTOCOL.md "Predictions written before any outcome".
PUBLISHED_PRIMARY_PERCENT = {'iid_perfect': 75.1, 'persistent_perfect': 53.6, 'iid_imperfect': 28.0}
PREDICTION1_TOLERANCE_PP = 2.0
PREDICTION3_REFERENCE = .9235
PREDICTION3_TOLERANCE = .05
# Original-study reference values (RESULT_REVIEW.md, verified to the stated digits).
ORIGINAL_REFERENCE_F = {'iid_perfect': 0.75131890, 'persistent_perfect': 0.53647154, 'iid_imperfect': 0.27993874}
ORIGINAL_REFERENCE_D0 = 1640.3636233482832
# Post hoc reference (robustness verification RESULTS.md section 6): leading-order
# control rates with ORIGINAL-workflow K at the alternative workflow's selected p.
RESULTS_MD_ORIGINAL_K_LEADING_CONTROL = {'iid_perfect': 0.92, 'persistent_perfect': 1.62, 'iid_imperfect': 1.64}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def finite(value):
    if value is None:
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def close(a, b, rel=1e-12, abs_=0.0):
    return a is not None and b is not None and math.isclose(a, b, rel_tol=rel, abs_tol=abs_)


# ----------------------------------------------------------------- loading
def seedbook():
    book = json.loads((wb.BRIDGE / 'seedbook.json').read_text())
    return book['all_seeds'], wb.seed_banks()


def load_original():
    rows = {}
    for path in sorted(ORIGINAL_RESULTS.glob('bridge_*/part_*/completed/*.json')):
        row = json.loads(path.read_text())
        if row.get('status') != 'complete' or row['seed'] in rows or path.stem != str(row['seed']):
            raise RuntimeError('Invalid or duplicate archived original record: ' + str(path))
        if row['executed_ticks'] != 146 * TICKS_PER_DAY:
            raise RuntimeError('Archived original record has a wrong horizon: ' + str(path))
        rows[row['seed']] = row
    return rows


def load_extended():
    if sha256(wb.EXTENDED_MANIFEST) != wb.EXTENDED_MANIFEST_SHA256:
        raise RuntimeError('Archived extended manifest changed')
    rows = {}
    for path in sorted(wb.EXTENDED_RESULTS.glob('part_*/completed/*.json')):
        row = json.loads(path.read_text())
        payload = {k: v for k, v in row.items() if k != 'payload_sha256'}
        if hashlib.sha256(canonical(payload).encode()).hexdigest() != row.get('payload_sha256'):
            raise RuntimeError('Extended record payload digest mismatch: ' + str(path))
        if row.get('status') != 'complete' or row['seed'] in rows:
            raise RuntimeError('Invalid or duplicate extended record: ' + str(path))
        rows[row['seed']] = row
    return rows


def load_alternative(records):
    manifest_path = records / 'manifest.json'
    if not manifest_path.exists():
        raise SystemExit(f'{manifest_path} missing: no alternative-workflow run in {records}')
    manifest = json.loads(manifest_path.read_text())
    rows = {}
    for path in sorted(records.glob('[0-9]*.json')):
        seed = int(path.stem)
        kind = json.loads(path.read_text())['surface']
        rows[seed] = run_batch.validate_record(path, manifest, seed, kind)
    failed = sorted(p.name for p in records.glob('failed_*.json'))
    current = run_batch.build_manifest()
    return rows, manifest, failed, dict(
        manifest_identity_sha256=manifest['identity_sha256'],
        manifest_sha256=sha256(manifest_path),
        current_code_identity_sha256=current['identity_sha256'],
        current_code_matches_manifest=current['identity_sha256'] == manifest['identity_sha256'],
        record_sha256={str(seed): sha256(records / f'{seed}.json') for seed in rows})


def planned_tasks():
    if not run_batch.BENCHMARK.exists():
        return None, None
    surface = run_batch.resolve_surface('auto')
    return surface, run_batch.plan(surface)


# --------------------------------------------------------------- statistics
def mean_uncertainty(values):
    """Frozen run.mean_uncertainty (statistics.fmean / stdev / 1.96)."""
    return wb.bridge_run.mean_uncertainty(list(values))


def ratio_block(x, y):
    """Frozen run.summarize formulas for one (K0, K1) vector pair."""
    x, y = [float(v) for v in x], [float(v) for v in y]
    n = len(x)
    xmean, ymean = statistics.fmean(x), statistics.fmean(y)
    ratio = None if xmean == 0 else ymean / xmean
    response = None if ratio is None else 1 - ratio
    residual = None if ratio is None else [(b - ratio * a) / xmean for a, b in zip(x, y)]
    se = None if residual is None or n < 2 else statistics.stdev(residual) / math.sqrt(n)
    return dict(n=n, K0=mean_uncertainty(x), K1=mean_uncertainty(y),
                paired_K1_minus_K0=mean_uncertainty([b - a for a, b in zip(x, y)]),
                K1_over_K0=ratio, limiting_response=response,
                limiting_response_delta_monte_carlo_se=se,
                limiting_response_normal_95_interval=None if se is None else [response - 1.96 * se, response + 1.96 * se])


class Bootstrap:
    """Synchronous whole-seed resampling for one ordered seed set."""

    def __init__(self, seeds):
        self.seeds = list(seeds)
        n = len(self.seeds)
        rng = np.random.default_rng(BOOTSTRAP_SEED)
        index = rng.integers(0, n, size=(BOOTSTRAP_DRAWS, n))
        self.index_sha256 = hashlib.sha256(index.astype('<i8').tobytes()).hexdigest()
        flat = (np.arange(BOOTSTRAP_DRAWS)[:, None] * n + index).ravel()
        self.weights = np.bincount(flat, minlength=BOOTSTRAP_DRAWS * n).reshape(BOOTSTRAP_DRAWS, n).astype(float)
        self.n = n

    def means(self, values):
        return self.weights @ np.asarray(values, dtype=float) / self.n

    def describe(self):
        return dict(draws=BOOTSTRAP_DRAWS, numpy_seed=BOOTSTRAP_SEED, n_seeds=self.n,
                    method='numpy.random.default_rng(seed).integers(0, n, (draws, n)); whole seeds in seedbook order; '
                           'percentile 95% (numpy linear interpolation); one index matrix for every statistic on this seed set',
                    index_matrix_sha256_int64LE=self.index_sha256)


def percentile_interval(samples, level=.95):
    samples = np.asarray(samples, dtype=float)
    defined = samples[np.isfinite(samples)]
    undefined = int(samples.size - defined.size)
    if defined.size == 0:
        return dict(level=level, interval=None, defined_draws=0, undefined_draws=undefined)
    tail = (1 - level) / 2 * 100
    low, high = np.percentile(defined, [tail, 100 - tail])
    return dict(level=level, interval=[float(low), float(high)], defined_draws=int(defined.size),
                undefined_draws=undefined)


def boot_response(boot, x, y):
    mx, my = boot.means(x), boot.means(y)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(mx == 0, np.nan, 1 - my / mx)


def paired_block(seeds, orig, alt, boot):
    """orig/alt: dict(K0=[...], K1=[...]) aligned with ``seeds``; boot over the same seeds."""
    xo, yo, xa, ya = (np.asarray(v, dtype=float) for v in (orig['K0'], orig['K1'], alt['K0'], alt['K1']))
    n = len(seeds)
    out = dict(n=n, original=ratio_block(xo, yo), alternative=ratio_block(xa, ya))
    fo, fa = out['original']['limiting_response'], out['alternative']['limiting_response']
    for arm, a, o in (('K0', xa, xo), ('K1', ya, yo)):
        diff = a - o
        block = mean_uncertainty(diff)
        block['bootstrap_95'] = percentile_interval(boot.means(diff))
        block['relative_change_of_mean'] = None if o.mean() == 0 else float(a.mean() / o.mean() - 1)
        with np.errstate(divide='ignore', invalid='ignore'):
            rel = boot.means(a) / boot.means(o) - 1
        block['relative_change_bootstrap_95'] = percentile_interval(np.where(np.isfinite(rel), rel, np.nan))
        out['delta_' + arm + '_alt_minus_orig'] = block
    delta = None if fo is None or fa is None else fa - fo
    se = None
    if delta is not None and n > 1:
        ro, ra = out['original']['K1_over_K0'], out['alternative']['K1_over_K0']
        psi = -(ya - ra * xa) / statistics.fmean(xa.tolist()) + (yo - ro * xo) / statistics.fmean(xo.tolist())
        se = statistics.stdev(psi.tolist()) / math.sqrt(n)
    bo, ba = boot_response(boot, xo, yo), boot_response(boot, xa, ya)
    out['delta_F_alt_minus_orig'] = dict(
        estimate=delta, delta_monte_carlo_se=se,
        normal_95_interval=None if se is None else [delta - 1.96 * se, delta + 1.96 * se],
        bootstrap_95=percentile_interval(ba - bo))
    out['original']['F_bootstrap_95'] = percentile_interval(bo)
    out['alternative']['F_bootstrap_95'] = percentile_interval(ba)
    return out


def q2_response(h0, h1, rho, e):
    q0 = wb.bridge_coefficient.stationary_q2(h0, rho, e)
    q1 = wb.bridge_coefficient.stationary_q2(h1, rho, e)
    return None if q0 == 0 else 1 - q1 / q0


def crossing(d_low, d_high, e_low=.9, e_high=1.):
    """Linear-interpolation root of D(e) on [e_low, e_high]; None without a sign change."""
    if d_low is None or d_high is None or d_low == d_high or (d_low > 0) == (d_high > 0):
        return None
    return e_low + (e_high - e_low) * d_low / (d_low - d_high)


def q2_crossing():
    """Exact root in e of q2-response(rho=.5,e) - q2-response(rho=0,e), h .90->.95."""
    def g(e):
        return q2_response(.9, .95, .5, e) - q2_response(.9, .95, 0., e)
    low, high = .9, 1.
    if (g(low) > 0) == (g(high) > 0):
        return None
    for _ in range(200):
        mid = (low + high) / 2
        if (g(mid) > 0) == (g(low) > 0):
            low = mid
        else:
            high = mid
    linear = crossing(g(.9), g(1.))
    return dict(exact_root=(low + high) / 2, linear_interpolation_between_e_0p9_and_1=linear,
                g_at_e0p9=g(.9), g_at_e1=g(1.),
                note='g(e) = q2 response(rho=.5, e) - q2 response(rho=0, e) at h .90->.95 (stationary two-opportunity comparator).')


# ------------------------------------------------------------- self-checks
def parse_table_s5():
    text = SUPPLEMENT.read_text()
    start = text.index(r'\label{tab:author-coefficients}')
    body = text[start:text.index(r'\bottomrule', start)]
    body = body.split(r'\midrule\endhead', 1)[1]
    rows = []
    pattern = re.compile(r'([\d.]+) & ([\d.]+) & ([\d.]+) & ([\d.]+) & ([\d.]+) & ([\d.]+) & '
                         r'([\d.]+) \[([\d.]+), ([\d.]+)\] & ([\d.]+)\\\\')
    for line in body.strip().splitlines():
        match = pattern.fullmatch(line.strip())
        if not match:
            raise RuntimeError('Unparsed Table S5 row: ' + line)
        v = [float(g) for g in match.groups()]
        rows.append(dict(h0=v[0], h1=v[1], rho=v[2], e=v[3], K0=v[4], K1=v[5], F=v[6], F_low=v[7], F_high=v[8], q2=v[9]))
    return rows


def pair_id(h0, h1, rho, e):
    return f'h{h0:g}_to_{h1:g}_rho{rho:g}_e{e:g}'


def self_checks(original, extended, all_seeds):
    checks = {}
    rows = [original[s] for s in all_seeds]
    _, config = wb.archived_scientific_config()
    frozen = wb.bridge_run.summarize(rows, config)
    aggregate_ok = sha256(ORIGINAL_AGGREGATE) == ORIGINAL_AGGREGATE_SHA256
    aggregate = json.loads(ORIGINAL_AGGREGATE.read_text())
    for legacy, _ in PRIMARY:
        mine = ratio_block([r['profiles'][legacy]['baseline']['K'] for r in rows],
                           [r['profiles'][legacy]['policy']['K'] for r in rows])
        fz, ag = frozen['profiles'][legacy], aggregate['profiles'][legacy]
        checks['original_' + legacy] = dict(
            F=mine['limiting_response'], reference_8_digits=ORIGINAL_REFERENCE_F[legacy],
            matches_reference=abs(mine['limiting_response'] - ORIGINAL_REFERENCE_F[legacy]) < 5e-9,
            equals_frozen_run_summarize=all(mine[k] == fz[k] for k in (
                'K0', 'K1', 'paired_K1_minus_K0', 'K1_over_K0', 'limiting_response',
                'limiting_response_delta_monte_carlo_se', 'limiting_response_normal_95_interval')),
            matches_archived_aggregate_summary=aggregate_ok and all(close(mine[a], ag[a]) for a in (
                'limiting_response', 'limiting_response_delta_monte_carlo_se', 'K1_over_K0'))
            and close(mine['K0']['mean'], ag['K0']['mean']) and close(mine['K1']['mean'], ag['K1']['mean']))
    d0 = frozen['observation_patient_days']['mean']
    checks['original_D0'] = dict(value=d0, reference=ORIGINAL_REFERENCE_D0, matches=close(d0, ORIGINAL_REFERENCE_D0))
    # Extended archive: legacy pairs equal original profiles; Table S5 reproduced.
    legacy_equal = all(extended[s]['pairs'][pid][arm] == original[s]['profiles'][legacy][arm]
                       for s in all_seeds for legacy, pid in PRIMARY for arm in ('baseline', 'policy'))
    hash_equal = all(extended[s]['ordered_bridge_events_sha256'] == original[s]['ordered_bridge_events_sha256'] for s in all_seeds)
    table = parse_table_s5()
    mismatches = []
    for entry in table:
        pid = pair_id(entry['h0'], entry['h1'], entry['rho'], entry['e'])
        block = ratio_block([extended[s]['pairs'][pid]['baseline']['K'] for s in all_seeds],
                            [extended[s]['pairs'][pid]['policy']['K'] for s in all_seeds])
        low, high = block['limiting_response_normal_95_interval']
        q2 = q2_response(entry['h0'], entry['h1'], entry['rho'], entry['e'])
        got = dict(K0=round(block['K0']['mean'], 2), K1=round(block['K1']['mean'], 2),
                   F=round(100 * block['limiting_response'], 2), F_low=round(100 * low, 2),
                   F_high=round(100 * high, 2), q2=round(100 * q2, 2))
        bad = {k: (got[k], entry[k]) for k in got if abs(got[k] - entry[k]) > 0.0051}
        if bad:
            mismatches.append(dict(pair=pid, differences=bad))
    checks['extended_archive'] = dict(
        legacy_pairs_bitwise_equal_original_profiles=legacy_equal, ordered_event_hashes_equal_original=hash_equal,
        table_s5_rows_parsed=len(table), table_s5_pairs_equal_surface=sorted(pair_id(r['h0'], r['h1'], r['rho'], r['e']) for r in table) == sorted(p['id'] for p in wb.SURFACE_PAIRS),
        table_s5_reproduced_to_2_decimals=not mismatches, table_s5_mismatches=mismatches,
        supplement_sha256=sha256(SUPPLEMENT))
    passed = (all(c['matches_reference'] and c['equals_frozen_run_summarize'] and c['matches_archived_aggregate_summary']
                  for k, c in checks.items() if k.startswith('original_i') or k.startswith('original_p'))
              and checks['original_D0']['matches'] and legacy_equal and hash_equal and not mismatches
              and len(table) == 22 and checks['extended_archive']['table_s5_pairs_equal_surface'])
    return passed, checks


# --------------------------------------------------------------- sections
def vectors(rows, seeds, getter):
    return dict(K0=[getter(rows[s])['baseline']['K'] for s in seeds],
                K1=[getter(rows[s])['policy']['K'] for s in seeds])


def histogram_means(rows, seeds, legacy):
    out = {}
    for arm in ('baseline', 'policy'):
        hists = [rows[s]['profiles'][legacy][arm]['K_by_opportunities'] for s in seeds]
        out[arm] = {label: statistics.fmean(h[label] for h in hists) for label in hists[0]}
    return out


def exposure(rows, seeds):
    keys = ('observation_patient_days', 'observation_susceptible_patient_days', 'observation_admissions',
            'observation_imported_colonizations', 'total_admissions')
    return {k: [float(rows[s][k]) for s in seeds] for k in keys}


def bank_split(rows, seeds, banks, getter):
    out = {}
    blocks = {}
    for bank in ('A', 'B'):
        members = [s for s in seeds if banks[s] == bank]
        if len(members) >= 2:
            v = vectors(rows, members, getter)
            blocks[bank] = ratio_block(v['K0'], v['K1'])
            out[bank] = dict(n=len(members), F=blocks[bank]['limiting_response'],
                             F_delta_se=blocks[bank]['limiting_response_delta_monte_carlo_se'],
                             K0_mean=blocks[bank]['K0']['mean'], K1_mean=blocks[bank]['K1']['mean'])
    if set(blocks) == {'A', 'B'}:
        a, b = blocks['A'], blocks['B']
        est = b['limiting_response'] - a['limiting_response']
        se = math.sqrt(a['limiting_response_delta_monte_carlo_se'] ** 2 + b['limiting_response_delta_monte_carlo_se'] ** 2)
        out['B_minus_A_F'] = dict(estimate=est, monte_carlo_se=se, normal_approximation_95_interval=[est - 1.96 * se, est + 1.96 * se])
    return out


def predictions(primary, surface_boot, surface_seeds, alt, extended):
    """The three prespecified checks; all are reported whatever the outcome."""
    result = {}
    p1 = {}
    for legacy, _ in PRIMARY:
        block = primary['paired'][legacy]
        f_alt = block['alternative']['limiting_response']
        f_orig = block['original']['limiting_response']
        diff_published = 100 * f_alt - PUBLISHED_PRIMARY_PERCENT[legacy]
        ci = block['delta_F_alt_minus_orig']['bootstrap_95']['interval']
        p1[legacy] = dict(F_alt_percent=100 * f_alt, published_original_percent=PUBLISHED_PRIMARY_PERCENT[legacy],
                          F_alt_minus_published_pp=diff_published,
                          F_alt_minus_exact_original_pp=100 * (f_alt - f_orig),
                          paired_bootstrap_95_pp=None if ci is None else [100 * ci[0], 100 * ci[1]],
                          holds=abs(diff_published) < PREDICTION1_TOLERANCE_PP,
                          paired_95_interval_inside_tolerance=None if ci is None else (
                              -PREDICTION1_TOLERANCE_PP < 100 * ci[0] and 100 * ci[1] < PREDICTION1_TOLERANCE_PP))
    result['prediction_1_small_change'] = dict(
        statement='The three primary limits change by less than 2 pp from 75.1%, 53.6%, 28.0% (point estimates).',
        by_profile=p1, holds_all=all(v['holds'] for v in p1.values()))
    if not surface_seeds:
        result['prediction_2_sign_reversal'] = result['prediction_3_crossing'] = dict(
            status='not_evaluable', reason='no seed carries the 22-row surface')
        return result

    def d_values(source_rows, seeds, getter):
        values = {}
        boots = {}
        for e in (1., .9, .8):
            pid5, pid0 = pair_id(.9, .95, .5, e), pair_id(.9, .95, 0., e)
            v5, v0 = vectors(source_rows, seeds, lambda r, p=pid5: getter(r, p)), vectors(source_rows, seeds, lambda r, p=pid0: getter(r, p))
            f5 = ratio_block(v5['K0'], v5['K1'])['limiting_response']
            f0 = ratio_block(v0['K0'], v0['K1'])['limiting_response']
            values[e] = f5 - f0
            boots[e] = boot_response(surface_boot, v5['K0'], v5['K1']) - boot_response(surface_boot, v0['K0'], v0['K1'])
        return values, boots

    alt_d, alt_boot = d_values(alt, surface_seeds, lambda r, p: r['pairs'][p])
    orig_d, orig_boot = d_values(extended, surface_seeds, lambda r, p: r['pairs'][p])
    p2 = {}
    for name, d, b in (('alternative', alt_d, alt_boot), ('original_same_seeds', orig_d, orig_boot)):
        p2[name] = dict(
            D_e1=d[1.], D_e1_bootstrap_95=percentile_interval(b[1.]), fraction_draws_D_e1_negative=float(np.mean(b[1.] < 0)),
            D_e0p8=d[.8], D_e0p8_bootstrap_95=percentile_interval(b[.8]), fraction_draws_D_e0p8_positive=float(np.mean(b[.8] > 0)),
            holds=d[1.] < 0 and d[.8] > 0)
    result['prediction_2_sign_reversal'] = dict(
        statement='For .90->.95, D(e) = F*(rho=.5,e) - F*(rho=0,e) is negative at e=1 and positive at e=.8.',
        n_seeds=len(surface_seeds), by_workflow=p2, holds=p2['alternative']['holds'])
    p3 = {}
    for name, d, b in (('alternative', alt_d, alt_boot), ('original_same_seeds', orig_d, orig_boot)):
        estar = crossing(d[.9], d[1.])
        low, high = b[.9], b[1.]
        with np.errstate(divide='ignore', invalid='ignore'):
            draws = np.where((low > 0) != (high > 0), .9 + .1 * low / (low - high), np.nan)
        p3[name] = dict(D_e0p9=d[.9], D_e1=d[1.], e_star=estar, e_star_bootstrap_95=percentile_interval(draws),
                        distance_from_reference=None if estar is None else abs(estar - PREDICTION3_REFERENCE),
                        holds=estar is not None and abs(estar - PREDICTION3_REFERENCE) <= PREDICTION3_TOLERANCE)
    result['prediction_3_crossing'] = dict(
        statement='e* (linear interpolation of D(e) between e=.9 and 1, rho=.5 vs 0, h .90->.95) lies within .05 of .9235.',
        reference_value=PREDICTION3_REFERENCE, tolerance=PREDICTION3_TOLERANCE, n_seeds=len(surface_seeds),
        comparator_recomputed=q2_crossing(), by_workflow=p3, holds=p3['alternative']['holds'])
    return result


def post_hoc(primary_alt, primary_orig_all, d0_alt, d0_orig):
    if sha256(FINITE_P_ANALYSIS) != FINITE_P_ANALYSIS_SHA256 or sha256(FINITE_P_SELECTION) != FINITE_P_SELECTION_SHA256:
        return dict(status='unavailable', reason='finite-p inputs missing or changed')
    analysis = json.loads(FINITE_P_ANALYSIS.read_text())
    selection = json.loads(FINITE_P_SELECTION.read_text())
    selected = {s['id']: s for s in selection['policy_config']['scenarios']}
    out = {}
    for entry in analysis['scenarios']:
        scenario = entry['scenario']
        if selected[scenario['id']]['p'] != scenario['p']:
            raise RuntimeError('Finite-p selection/analysis disagree')
        legacy = next(name for name, pid in PRIMARY
                      if wb.PRIMARY_PAIR_IDS[name] == pair_id(scenario['h0'], scenario['h1'], scenario['rho'], scenario['e']))
        p = scenario['p']
        alt, orig = primary_alt[legacy], primary_orig_all[legacy]
        row = dict(finite_p_profile=scenario['id'], selected_p=p, finite_p_n_pairs=entry['n'],
                   finite_p_F=entry['F'], finite_p_F_95=entry['F95']['interval'],
                   limit_alternative=alt['limiting_response'], limit_original=orig['limiting_response'],
                   finite_p_F_minus_limit_alternative_pp=100 * (entry['F'] - alt['limiting_response']),
                   finite_p_F_minus_limit_original_pp=100 * (entry['F'] - orig['limiting_response']),
                   limit_alternative_inside_finite_p_95=entry['F95']['interval'][0] <= alt['limiting_response'] <= entry['F95']['interval'][1],
                   measured_control_rate=entry['control_rate'], measured_policy_rate=entry['policy_rate'])
        for label, block, d0 in (('alternative_K', alt, d0_alt), ('original_K', orig, d0_orig)):
            r0 = 1000 * p * p * block['K0']['mean'] / d0
            r1 = 1000 * p * p * block['K1']['mean'] / d0
            row['leading_order_' + label] = dict(D0=d0, control_rate=r0, policy_rate=r1,
                                                  control_ratio_to_measured=r0 / entry['control_rate'],
                                                  policy_ratio_to_measured=r1 / entry['policy_rate'])
        row['original_K_control_rate_matches_RESULTS_md_2dp'] = round(row['leading_order_original_K']['control_rate'], 2) == RESULTS_MD_ORIGINAL_K_LEADING_CONTROL[legacy]
        out[legacy] = row
    return dict(label='POST HOC: the finite-p alternative-workflow results were seen before this protocol; not prespecified.',
                formula='leading-order rate R_h = 1000 p^2 mean(K_h) / D0 per 1,000 patient-days (D0 = mean observation patient-days)',
                inputs={str(FINITE_P_ANALYSIS): FINITE_P_ANALYSIS_SHA256, str(FINITE_P_SELECTION): FINITE_P_SELECTION_SHA256},
                by_profile=out)


def workflow_diagnostics(alt, seeds):
    counters = [alt[s]['workflow']['counters_observation'] for s in seeds]
    totals = [alt[s]['workflow']['counters_total'] for s in seeds]
    keys = sorted(counters[0])
    out = dict(per_record_observation_mean={k: statistics.fmean(c[k] for c in counters) for k in keys},
               per_record_total_mean={k: statistics.fmean(c[k] for c in totals) for k in keys},
               max_over_records={k: max(c[k] for c in totals) for k in keys})
    out['deferral_patient_days_per_observation_window_mean'] = out['per_record_observation_mean']['workflow_events.isolation_deferred_no_empty_bed'] / TICKS_PER_DAY
    out['no_swap_all_records'] = all(c['audit.swaps'] == 0 and c['audit.moves_into_occupied_bed'] == 0
                                     and c['audit.moves_from_isolation_bed'] == 0 for c in totals)
    mixes = [alt[s]['workflow']['patient_contacts_by_period_status_bed'] for s in seeds]
    labels = sorted({k for m in mixes for k in m})
    out['patient_contacts_by_period_status_bed_mean'] = {k: statistics.fmean(m.get(k, 0) for m in mixes) for k in labels}
    obs_carrier = [k for k in labels if k.startswith('observation|carrier|')]
    total_obs_carrier = sum(out['patient_contacts_by_period_status_bed_mean'][k] for k in obs_carrier)
    out['observation_carrier_contacts_share_in_shared_beds'] = None if total_obs_carrier == 0 else (
        out['patient_contacts_by_period_status_bed_mean'].get('observation|carrier|shared_bed', 0) / total_obs_carrier)
    return out


# --------------------------------------------------------------- main
def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w') as handle:
        handle.write(json.dumps(value, indent=1, sort_keys=True, allow_nan=False) + '\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--records', type=Path, default=run_batch.RECORDS)
    parser.add_argument('--output', type=Path, default=None)
    args = parser.parse_args(argv)
    started = time.time()

    all_seeds, banks = seedbook()
    original = load_original()
    extended = load_extended()
    if set(original) != set(all_seeds) or set(extended) != set(all_seeds):
        raise RuntimeError('Archived original/extended records do not cover the 256-seed seedbook')
    checks_pass, checks = self_checks(original, extended, all_seeds)

    alt, manifest, failed, binding = load_alternative(args.records)
    seeds = [s for s in all_seeds if s in alt]
    if not seeds:
        raise SystemExit('No alternative-workflow records in ' + str(args.records))
    surface_seeds = [s for s in seeds if alt[s]['surface'] == 'full22']
    plan_surface, tasks = planned_tasks()
    if tasks is not None and args.records.resolve() == run_batch.RECORDS.resolve():
        missing = [seed for seed, _, kind in tasks if seed not in alt or (kind == 'full22' and alt[seed]['surface'] != 'full22')]
        status = 'complete' if not missing and not failed else 'partial'
    else:
        missing, status = None, 'partial_or_test_bank'
    _, config = wb.archived_scientific_config()

    boot = Bootstrap(seeds)
    primary = dict(n_seeds=len(seeds), bootstrap=boot.describe(), paired={}, original_all_256={}, banks={})
    for legacy, pid in PRIMARY:
        o = vectors(original, seeds, lambda r, k=legacy: r['profiles'][k])
        a = vectors(alt, seeds, lambda r, k=legacy: r['profiles'][k])
        primary['paired'][legacy] = paired_block(seeds, o, a, boot)
        primary['paired'][legacy]['pair_id'] = pid
        allv = vectors(original, all_seeds, lambda r, k=legacy: r['profiles'][k])
        primary['original_all_256'][legacy] = ratio_block(allv['K0'], allv['K1'])
        primary['banks'][legacy] = dict(
            alternative=bank_split(alt, seeds, banks, lambda r, k=legacy: r['profiles'][k]),
            original_same_seeds=bank_split(original, seeds, banks, lambda r, k=legacy: r['profiles'][k]))
        primary['paired'][legacy]['K_by_opportunities_mean'] = dict(
            alternative=histogram_means(alt, seeds, legacy), original=histogram_means(original, seeds, legacy))
        primary['paired'][legacy]['stationary_q2_response_for_comparison_only'] = q2_response(
            config['h0'], config['h1'], *next((p['rho'], p['e']) for p in config['profiles'] if p['id'] == legacy))
    # Frozen summarizer agreement on the alternative rows (same conventions as the original study).
    frozen_alt = wb.bridge_run.summarize([alt[s] for s in seeds], config)
    primary['frozen_run_summarize_alternative'] = frozen_alt
    frozen_agree = all(frozen_alt['profiles'][legacy][k] == primary['paired'][legacy]['alternative'][k]
                       for legacy, _ in PRIMARY
                       for k in ('limiting_response', 'limiting_response_delta_monte_carlo_se', 'K1_over_K0'))
    checks['alternative_ratio_block_equals_frozen_run_summarize'] = frozen_agree
    checks_pass = checks_pass and frozen_agree

    surface = None
    surface_boot = None
    if surface_seeds:
        surface_boot = boot if surface_seeds == seeds else Bootstrap(surface_seeds)
        surface = dict(n_seeds=len(surface_seeds), banks=sorted({banks[s] for s in surface_seeds}),
                       bootstrap=surface_boot.describe(), pairs={})
        for pair in wb.SURFACE_PAIRS:
            pid = pair['id']
            o = vectors(extended, surface_seeds, lambda r, k=pid: r['pairs'][k])
            a = vectors(alt, surface_seeds, lambda r, k=pid: r['pairs'][k])
            block = paired_block(surface_seeds, o, a, surface_boot)
            allv = vectors(extended, all_seeds, lambda r, k=pid: r['pairs'][k])
            block.update(parameters={k: pair[k] for k in ('h0', 'h1', 'rho', 'e')}, legacy_id=pair.get('legacy_id'),
                         original_all_256=ratio_block(allv['K0'], allv['K1']),
                         stationary_q2_response_for_comparison_only=q2_response(pair['h0'], pair['h1'], pair['rho'], pair['e']))
            surface['pairs'][pid] = block

    # Seeds whose alternative-workflow history is event-for-event the original one
    # (no swap ever needed in the original, no deferral in the alternative). Their
    # K must then be bitwise identical: a consistency check on the adapter.
    identical = [s for s in seeds if alt[s]['ordered_bridge_events_sha256'] == original[s]['ordered_bridge_events_sha256']]
    identical_k_equal = all(alt[s]['profiles'][k] == original[s]['profiles'][k] for s in identical for k, _ in PRIMARY)
    checks['identical_history_seeds_have_identical_primary_K'] = identical_k_equal
    checks_pass = checks_pass and identical_k_equal
    history_identity = dict(
        n_seeds_identical_ordered_event_hash=len(identical), n_seeds=len(seeds), seeds=identical,
        primary_K_bitwise_equal_on_identical_seeds=identical_k_equal,
        note='Identical hashes mean the alternative rules never changed this background; differences in K come only from the other seeds.')

    ex_alt, ex_orig = exposure(alt, seeds), exposure(original, seeds)
    exposure_block = {k: dict(alternative=mean_uncertainty(ex_alt[k]), original=mean_uncertainty(ex_orig[k]),
                              paired_alt_minus_orig=mean_uncertainty([a - o for a, o in zip(ex_alt[k], ex_orig[k])]))
                      for k in ex_alt}
    d0_alt = statistics.fmean(ex_alt['observation_patient_days'])
    d0_orig_all = statistics.fmean(original[s]['observation_patient_days'] for s in all_seeds)

    summary = dict(
        schema_version=1, study='workflow_coefficient_2026-10-21', status=status,
        records_dir=str(args.records), n_alternative_records=len(seeds), failed_records=failed,
        planned_surface=plan_surface, missing_planned_seeds=missing,
        alternative_binding=binding,
        seeds_by_bank={b: sum(1 for s in seeds if banks[s] == b) for b in 'AB'},
        estimand='F* = 1 - mean(K1)/mean(K0) over backgrounds; paired workflow differences on the same seeds.',
        self_checks_pass=checks_pass, self_checks=checks,
        primary=primary, surface=surface,
        predictions=predictions(primary, surface_boot, surface_seeds, alt, extended),
        post_hoc_finite_p=post_hoc({k: primary['paired'][k]['alternative'] for k, _ in PRIMARY},
                                   primary['original_all_256'], d0_alt, d0_orig_all),
        history_identity=history_identity,
        exposure=exposure_block,
        workflow_diagnostics=workflow_diagnostics(alt, seeds),
        limitation='The coefficient is a fixed-horizon p->0 limit; it gives no finite-p error bound.',
        inputs=dict(original_records=str(ORIGINAL_RESULTS), extended_records=str(wb.EXTENDED_RESULTS),
                    original_aggregate_sha256=sha256(ORIGINAL_AGGREGATE), supplement_sha256=sha256(SUPPLEMENT),
                    frozen_sha256=wb.verify_frozen(), adapter_sha256=wb.adapter_sha256(),
                    summarize_sha256=sha256(Path(__file__)), run_batch_sha256=sha256(CODE / 'run_batch.py')),
        runtime=dict(python=sys.version, numpy=np.__version__),
        computed_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        compute_seconds=None)
    summary['compute_seconds'] = time.time() - started
    if args.output is not None:
        output = args.output
    else:
        output = STUDY / ('summary.json' if status == 'complete' else 'progress_summary.json')
    atomic_json(output, summary)
    print(json.dumps(dict(output=str(output), status=status, n=len(seeds), surface_n=len(surface_seeds),
                          self_checks_pass=checks_pass,
                          F_alt={k: round(primary['paired'][k]['alternative']['limiting_response'], 6) for k, _ in PRIMARY},
                          F_orig_same_seeds={k: round(primary['paired'][k]['original']['limiting_response'], 6) for k, _ in PRIMARY}),
                     indent=1))
    return 0 if checks_pass else 1


if __name__ == '__main__':
    sys.exit(main())
