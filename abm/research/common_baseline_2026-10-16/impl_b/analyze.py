#!/usr/bin/env python3
"""Common-baseline prevention (post hoc, secondary): independent implementation B.

Implements PROTOCOL.md (abm/research/common_baseline_2026-10-16/PROTOCOL.md):
reproduction gate, primary log-log interpolation estimator at B*=3.0, whole-seed
within-bank bootstrap (50,000 draws, seed 20261016), pairwise contrasts with a
Bonferroni-3 level, and descriptive sensitivity analyses S1-S4. D1 is not done here.

All inputs are read-only. The only writes are results.json (and the caller's
run log) inside this script's own folder. No archived analysis code is imported.
"""
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ABM = Path(__file__).resolve().parents[3]
OUTER = (Path(__file__).resolve().parents[4] / 'research')

# ----------------------------------------------------------------- input paths
RS_STUDY = ABM / 'research/rule_surface_followup_2026-10-11'
RS_REPORT = ABM / 'research/rule_surface_results_2026-10-12/completed_analysis/report.json'
OLD_INDEX = RS_STUDY / 'old_input_index.json'
UPPER_COLLECTION = OUTER / 'rule_surface_followup_2026-10-11/collections/independent_upper_20261012_0125_utc'
UPPER_ROOT = UPPER_COLLECTION / 'restored/repo/results/rule_surface_followup_2026-10-11/independent_upper'
UPPER_CMANIFEST = UPPER_COLLECTION / 'received/collection_manifest.json'
HOLDOUT_STUDY = ABM / 'research/full_icu_confirmation_2026-10-08'
HOLDOUT_REPORT = HOLDOUT_STUDY / 'holdout_report/report.json'
SENS_STUDY = ABM / 'research/author_notes_2026-10-10'
SENS_COLLECTION = OUTER / 'author_notes_2026-10-10/collections/sensitivity_20261011_101348_utc'
SENS_ROOT = SENS_COLLECTION / 'results/sensitivity'
SENS_SHA_MANIFEST = SENS_COLLECTION / 'SHA256_MANIFEST.json'
SENS_REPORT = OUTER / 'author_notes_2026-10-10/analysis_2026-10-11/sensitivity_replay_20261011_1015_utc/report.json'

# ------------------------------------------------------------ fixed constants
TARGET = 3.0
DRAWS = 50000
BOOT_SEED = 20261016
CHUNK = 1000
CI_PRIMARY = 0.95
CI_CONTRAST = 1 - 0.05 / 3
GATE_TOL = 1e-12
PROFILES = ('independent', 'clustered', 'imperfect')
PAIRS = (('independent', 'clustered'), ('independent', 'imperfect'), ('clustered', 'imperfect'))
LOG = []


def log(*args):
    text = ' '.join(str(a) for a in args)
    LOG.append(text)
    print(text, flush=True)


def need(cond, msg):
    if not cond:
        raise RuntimeError(msg)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load_json_bytes(path):
    raw = Path(path).read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def rel_err(a, b):
    if a == b:
        return 0.0
    return abs(a - b) / max(abs(b), 1e-300)


# ---------------------------------------------------------------- record I/O
def extract(row, path):
    """Return the per-seed observation vector and identity fields of one record."""
    need(row['status'] == 'complete', f'incomplete record {path}')
    out = {}
    for arm, h_key in (('baseline', 'h0'), ('policy', 'h1')):
        a = row['arms'][arm]
        need(a['h'] == row['scenario'][h_key], f'arm h mismatch {path}')
        ob = a['observation']
        need(ob['days'] == 56 and ob['ticks'] == 2592 * 56, f'observation window not 56 days {path}')
        c, d = ob['acquired_colonizations'], ob['patient_days']
        need(type(c) is int and c >= 0 and d > 0, f'bad counter {path}')
        # Field semantics: the record's own rate must equal 1000*acquired/patient_days.
        need(math.isclose(ob['acquired_per_1000_patient_days'], 1000 * c / d, rel_tol=1e-12, abs_tol=0),
             f'rate semantics differ {path}')
        need(sum(b['acquired_colonizations'] for b in a['observation_blocks']) == c and
             len(a['observation_blocks']) == 56, f'daily blocks do not sum {path}')
        out[arm] = (c, d)
    need(row['model_kwargs']['seed'] == row['seed'] and
         row['model_kwargs']['prob_transmission'] == row['scenario']['p'], f'kwargs seed/p mismatch {path}')
    return dict(seed=row['seed'], scenario=row['scenario'], kwargs=row['model_kwargs'],
                execution=row.get('execution', {}), C0=out['baseline'][0], D0=out['baseline'][1],
                C1=out['policy'][0], D1=out['policy'][1])


# ------------------------------------------------------------------ loaders
def load_independent_upper(hash_report):
    rs, _ = load_json_bytes(RS_REPORT)
    bind = rs['new_bank_bindings']['independent_upper']
    cman, _ = load_json_bytes(UPPER_CMANIFEST)
    cfg_path = RS_STUDY / 'configs/independent_upper.json'
    cfg, cfg_sha = load_json_bytes(cfg_path)
    freeze, _ = load_json_bytes(RS_STUDY / 'freeze.json')
    checks = {}
    checks['config_sha_vs_rule_surface_freeze'] = cfg_sha == freeze['files']['research/rule_surface_followup_2026-10-11/configs/independent_upper.json']
    checks['config_sha_vs_collection_member'] = cfg_sha == cman['member_sha256']['repo/research/rule_surface_followup_2026-10-11/configs/independent_upper.json']['sha256']
    checks['manifest_sha_vs_report'] = sha256_file(UPPER_ROOT / 'manifest.json') == bind['manifest_sha256']
    checks['binding_sha_vs_report'] = sha256_file(UPPER_ROOT / 'study_binding.json') == bind['binding_sha256']
    files = sorted((UPPER_ROOT / 'completed').iterdir())
    need(len(files) == 256 and len(bind['input_sha256']) == 256, 'independent_upper inventory')
    rec_vs_report = rec_vs_collection = 0
    rows = []
    for p in files:
        raw = p.read_bytes(); h = hashlib.sha256(raw).hexdigest()
        rec_vs_report += h == bind['input_sha256']['completed/' + p.name]
        rec_vs_collection += h == cman['member_sha256']['repo/results/rule_surface_followup_2026-10-11/independent_upper/completed/' + p.name]['sha256']
        rows.append(extract(json.loads(raw), p))
    checks['records_match_report_input_sha256'] = f'{rec_vs_report}/256'
    checks['records_match_collection_manifest'] = f'{rec_vs_collection}/256'
    ok = all(v is True for v in checks.values() if isinstance(v, bool)) and rec_vs_report == 256 and rec_vs_collection == 256
    hash_report['independent_upper'] = dict(passed=ok, **checks)
    need(ok, 'independent_upper hash check failed')
    seeds = cfg['seeds']
    by_seed = {r['seed']: r for r in rows}
    need(set(by_seed) == set(seeds) and len(by_seed) == 256, 'upper seeds differ from config')
    need(all(r['scenario'] == cfg['scenarios'][0] for r in rows), 'upper scenario differs from config')
    return seeds, [by_seed[s] for s in seeds], cfg


def load_old_holdout(hash_report):
    index, index_sha = load_json_bytes(OLD_INDEX)
    rs, _ = load_json_bytes(RS_REPORT)
    hr, _ = load_json_bytes(HOLDOUT_REPORT)
    freeze, _ = load_json_bytes(RS_STUDY / 'freeze.json')
    checks = {}
    checks['old_index_sha_vs_rule_surface_report'] = index_sha == rs['old_input_index_sha256']
    checks['old_index_sha_vs_rule_surface_freeze'] = index_sha == freeze['files']['research/rule_surface_followup_2026-10-11/old_input_index.json']
    for key in ('config', 'seedbook', 'report'):
        checks[f'{key}_sha_vs_index'] = sha256_file(ABM / index[key]['path']) == index[key]['sha256']
    for name, h in index['manifests'].items():
        checks['manifest_sha_vs_index:' + name.split('/')[-2]] = sha256_file(ABM / name) == h
    # Holdout report also binds each record per site manifest.
    site_hashes = {}
    for bank in hr['banks']:
        site = {'68527ef19490ad27a67dbf96333b3cc0646ab9d6def1fe0eb854c3919b34d565': 'holdout_cluster',
                '59c2e7978aff6a683952156e809e0eae365617d252438c3e614647a1d7e6db62': 'holdout_local'}[bank['manifest_sha256']]
        site_hashes[site] = bank['completed_records_sha256']
    need(len(index['records']) == 1536, 'old index record count')
    on_disk = {f'results/full_icu_confirmation_2026-10-08/{site}/completed/{p.name}'
               for site in ('holdout_cluster', 'holdout_local')
               for p in (ABM / f'results/full_icu_confirmation_2026-10-08/{site}/completed').iterdir()}
    checks['completed_inventory_equals_index'] = on_disk == set(index['records'])
    n_idx = n_hr = 0
    rows = []
    for name, h in sorted(index['records'].items()):
        raw = (ABM / name).read_bytes(); hh = hashlib.sha256(raw).hexdigest()
        n_idx += hh == h
        site = name.split('/')[2]
        n_hr += hh == site_hashes[site][name.split('/')[-1]]
        r = extract(json.loads(raw), name)
        meta = index['record_metadata'][name]
        need(meta['seed'] == r['seed'] and meta['scenario'] == r['scenario']['id'], 'old index metadata mismatch')
        r['site'] = meta['execution_site']; r['index_bank'] = meta['statistical_bank']
        rows.append(r)
    checks['records_match_old_input_index'] = f'{n_idx}/1536'
    checks['records_match_holdout_report_bank_hashes'] = f'{n_hr}/1536'
    ok = all(v is True for v in checks.values() if isinstance(v, bool)) and n_idx == 1536 and n_hr == 1536
    hash_report['old_holdout'] = dict(passed=ok, **checks)
    need(ok, 'old holdout hash check failed')
    seedbook = json.loads((ABM / index['seedbook']['path']).read_text())
    banks = seedbook['seed_blocks']['holdout']
    need(banks == index['ordered_statistical_banks'], 'old bank identity')
    cfg = json.loads((ABM / index['config']['path']).read_text())
    need(cfg['seeds'] == banks['A'] + banks['B'], 'old config seed order is not A+B')
    need(index['independent_upper_comparison_seeds'] == banks['A'][:128] + banks['B'][:128], 'upper subset definition')
    table = {}
    for r in rows:
        key = (r['scenario']['id'], r['seed'])
        need(key not in table, 'duplicate old record'); table[key] = r
        need(r['index_bank'] == ('A' if r['seed'] in set(banks['A']) else 'B'), 'index bank label disagrees with seedbook')
    return cfg, banks, table, index


def load_sensitivity(hash_report):
    rep, _ = load_json_bytes(SENS_REPORT)
    shaman, shaman_sha = load_json_bytes(SENS_SHA_MANIFEST)
    receipt = json.loads((SENS_COLLECTION / 'collection_receipt.json').read_text())
    cfg_path = SENS_STUDY / 'configs/sensitivity.json'
    cfg, cfg_sha = load_json_bytes(cfg_path)
    sb_path = SENS_STUDY / 'seedbook.json'
    seedbook, sb_sha = load_json_bytes(sb_path)
    freeze, _ = load_json_bytes(SENS_STUDY / 'freeze.json')
    checks = {}
    checks['config_sha_vs_report'] = cfg_sha == rep['config_sha256']
    checks['config_sha_vs_freeze'] = cfg_sha == freeze['files']['research/author_notes_2026-10-10/configs/sensitivity.json']
    checks['seedbook_sha_vs_freeze'] = sb_sha == freeze['files']['research/author_notes_2026-10-10/seedbook.json']
    checks['manifest_sha_vs_report'] = sha256_file(SENS_ROOT / 'manifest.json') == rep['manifest_sha256']
    checks['manifest_sha_vs_SHA256_MANIFEST'] = sha256_file(SENS_ROOT / 'manifest.json') == shaman['files']['manifest.json']['sha256']
    checks['study_binding_sha_vs_SHA256_MANIFEST'] = sha256_file(SENS_ROOT / 'study_binding.json') == shaman['files']['study_binding.json']['sha256']
    checks['SHA256_MANIFEST_sha_vs_collection_receipt'] = shaman_sha == receipt['manifest']['sha256']
    files = sorted((SENS_ROOT / 'completed').iterdir())
    need(len(files) == 1536 and len(rep['input_sha256']) == 1536, 'sensitivity inventory')
    n_rep = n_man = 0
    rows = []
    for p in files:
        raw = p.read_bytes(); h = hashlib.sha256(raw).hexdigest()
        n_rep += h == rep['input_sha256'][p.name]
        n_man += h == shaman['files']['completed/' + p.name]['sha256']
        rows.append(extract(json.loads(raw), p))
    checks['records_match_report_input_sha256'] = f'{n_rep}/1536'
    checks['records_match_SHA256_MANIFEST'] = f'{n_man}/1536'
    ok = all(v is True for v in checks.values() if isinstance(v, bool)) and n_rep == 1536 and n_man == 1536
    hash_report['sensitivity'] = dict(passed=ok, **checks)
    need(ok, 'sensitivity hash check failed')
    book = seedbook['sensitivity']
    need(book['A'] + book['B'] == cfg['seeds'] and len(book['A']) == len(book['B']) == 128, 'sensitivity banks')
    table = {}
    for r in rows:
        key = (r['scenario']['id'], r['seed'])
        need(key not in table, 'duplicate sensitivity record'); table[key] = r
    need(len(table) == 6 * 256, 'sensitivity coverage')
    return cfg, book, table


# ------------------------------------------------------------ condition data
def condition(name, rows, scenario):
    """rows ordered by the module's seed list; returns arrays and pooled point values."""
    need(all(r['scenario'] == scenario for r in rows), f'scenario mismatch in {name}')
    V = np.array([[r['C0'], r['D0'], r['C1'], r['D1']] for r in rows], dtype=np.float64)
    C0 = sum(r['C0'] for r in rows); C1 = sum(r['C1'] for r in rows)
    D0 = math.fsum(r['D0'] for r in rows); D1 = math.fsum(r['D1'] for r in rows)
    return dict(name=name, scenario=scenario, p=scenario['p'], seeds=[r['seed'] for r in rows], rows=rows, V=V,
                C0=C0, D0=D0, C1=C1, D1=D1, **pooled(C0, D0, C1, D1))


def pooled(C0, D0, C1, D1):
    B = 1000.0 * C0 / D0; P = 1000.0 * C1 / D1
    return dict(control_rate=B, policy_rate=P, F=1.0 - P / B)


def bank_point(cond, idx):
    rr = [cond['rows'][i] for i in idx]
    return pooled(sum(r['C0'] for r in rr), math.fsum(r['D0'] for r in rr),
                  sum(r['C1'] for r in rr), math.fsum(r['D1'] for r in rr))


# ------------------------------------------------------------ estimators
def interp(B_lo, P_lo, B_hi, P_hi, p_lo, p_hi, target):
    """Primary log-log interpolation; arrays or scalars. Returns dict incl. NaN-guard mask."""
    B_lo, P_lo, B_hi, P_hi = map(np.asarray, (B_lo, P_lo, B_hi, P_hi))
    undefined = (B_lo <= 0) | (B_hi <= 0) | (P_lo <= 0) | (P_hi <= 0) | (B_hi == B_lo)
    with np.errstate(divide='ignore', invalid='ignore'):
        lB_lo, lB_hi, lP_lo, lP_hi = np.log(B_lo), np.log(B_hi), np.log(P_lo), np.log(P_hi)
        w = (math.log(target) - lB_lo) / (lB_hi - lB_lo)
        F_lo = 1 - P_lo / B_lo; F_hi = 1 - P_hi / B_hi
        # Protocol formula verbatim: F* = 1 - (1-F_lo)^(1-w) (1-F_hi)^w
        Fstar = 1 - np.power(1 - F_lo, 1 - w) * np.power(1 - F_hi, w)
        # S1: arithmetic interpolation of F in ln B (same weight w)
        F_s1 = (1 - w) * F_lo + w * F_hi
        pstar = np.exp((1 - w) * math.log(p_lo) + w * math.log(p_hi))
    Fstar = np.where(undefined, np.nan, Fstar)
    F_s1 = np.where(undefined, np.nan, F_s1)
    w = np.where(undefined, np.nan, w)
    extrap = (B_lo >= target) | (B_hi <= target)
    return dict(w=w, Fstar=Fstar, F_s1=F_s1, pstar=np.where(undefined, np.nan, pstar),
                undefined=undefined, extrapolated=extrap)


def quad_eval_root(x, yB, yP, target, x_ref):
    """S2: exact quadratic through three (ln p, ln rate) points, per draw (vectorised).

    x: (3,) ln p for low, center, high. yB, yP: (..., 3). Solves q_B(x*) = ln target.
    Root rule: real roots inside [x_low, x_high]; if two, the one nearest x_ref
    (the primary center-high log-linear solution of the same draw); none -> undefined.
    """
    xl, xc, xh = x
    ul, uh = xl - xc, xh - xc  # centred abscissae (center at 0)
    def coeffs(y):
        yl, yc, yh = y[..., 0], y[..., 1], y[..., 2]
        s_l = (yc - yl) / (0 - ul); s_h = (yh - yc) / (uh - 0)
        a = (s_h - s_l) / (uh - ul)
        b = s_h - a * uh
        return a, b, yc  # y(u) = a u^2 + b u + c
    aB, bB, cB = coeffs(yB)
    aP, bP, cP = coeffs(yP)
    c0 = cB - math.log(target)
    disc = bB * bB - 4 * aB * c0
    with np.errstate(divide='ignore', invalid='ignore'):
        sq = np.sqrt(np.where(disc >= 0, disc, np.nan))
        # numerically stable roots
        qq = -0.5 * (bB + np.copysign(sq, bB))
        r1 = qq / aB
        r2 = c0 / qq
        lin = -c0 / bB
    small = np.abs(aB) < 1e-14
    r1 = np.where(small, lin, r1); r2 = np.where(small, np.nan, r2)
    u_ref = x_ref - xc
    lo, hi = ul - 1e-15, uh + 1e-15
    in1 = np.isfinite(r1) & (r1 >= lo) & (r1 <= hi)
    in2 = np.isfinite(r2) & (r2 >= lo) & (r2 <= hi)
    pick2 = in2 & (~in1 | (np.abs(r2 - u_ref) < np.abs(r1 - u_ref)))
    u = np.where(pick2, r2, np.where(in1, r1, np.nan))
    two_roots = in1 & in2
    lnP = aP * u * u + bP * u + cP
    F = 1 - np.exp(lnP - math.log(target))
    return dict(Fstar=F, pstar=np.exp(u + xc), undefined=~np.isfinite(u), two_in_range=two_roots,
                curvature_B=aB, curvature_P=aP)


def pct(values, level):
    v = np.asarray(values, dtype=float)
    finite = v[np.isfinite(v)]
    lo_q, hi_q = (1 - level) / 2, 1 - (1 - level) / 2
    if finite.size == 0:
        return dict(level=level, interval=None, undefined_draws=int(v.size))
    lo, hi = np.quantile(finite, [lo_q, hi_q], method='linear')
    return dict(level=level, interval=[float(lo), float(hi)], undefined_draws=int(v.size - finite.size),
                conditional_on_defined_draws=bool(v.size != finite.size))


# ------------------------------------------------------------ bootstrap
def resample_counts(rng, n_bank, size):
    """Whole-seed resampling with replacement of n_bank seeds: returns multiplicity counts."""
    idx = rng.integers(0, n_bank, size=(size, n_bank))
    flat = idx + n_bank * np.arange(size)[:, None]
    return np.bincount(flat.ravel(), minlength=size * n_bank).reshape(size, n_bank).astype(np.float64)


def module_weights(rng, banks, n, size):
    """banks: list of index lists partitioning range(n). Independent resampling within each bank."""
    W = np.zeros((size, n))
    for bank in banks:
        W[:, bank] = resample_counts(rng, len(bank), size)
    return W


def boot_module_sums(rng, banks, V, draws, order_log, label):
    """V: (n_seeds, n_cond, 4). Returns (draws, n_cond, 4) bootstrap sums."""
    n, k, _ = V.shape
    flatV = V.reshape(n, k * 4)
    out = np.empty((draws, k, 4))
    wh = hashlib.sha256()
    for start in range(0, draws, CHUNK):
        size = min(CHUNK, draws - start)
        W = module_weights(rng, banks, n, size)
        wh.update(W.astype('<i8').tobytes())
        out[start:start + size] = (W @ flatV).reshape(size, k, 4)
    order_log.append(dict(stream=label, draws=draws, chunk=CHUNK, banks=[len(b) for b in banks],
                          weight_sha256_int64LE=wh.hexdigest()))
    return out


def rates_from_sums(S):
    with np.errstate(divide='ignore', invalid='ignore'):
        B = 1000 * S[..., 0] / S[..., 1]
        P = 1000 * S[..., 2] / S[..., 3]
    return B, P


# ------------------------------------------------------------ main
def main():
    t0 = time.time()
    log('impl_b analyze.py start', time.strftime('%Y-%m-%dT%H:%M:%S%z'))
    log('python', sys.version.split()[0], 'numpy', np.__version__)
    need(np.__version__ == '1.26.4', 'numpy 1.26.4 required')
    hashes = {}

    # ---------------- load and hash-verify inputs
    up_seeds, up_rows, up_cfg = load_independent_upper(hashes)
    old_cfg, old_banks, old_table, index = load_old_holdout(hashes)
    s_cfg, s_book, s_table = load_sensitivity(hashes)
    log('hash checks:', json.dumps(hashes, sort_keys=True))

    # ---------------- assemble conditions
    old_sc = {s['id']: s for s in old_cfg['scenarios']}
    s_sc = {s['id']: s for s in s_cfg['scenarios']}
    # independent module: seed order = upper config = old A[:128] + B[:128]
    need(up_seeds == old_banks['A'][:128] + old_banks['B'][:128], 'upper seed order')
    ind_lo = condition('independent_center_matched256',
                       [old_table[('independent_target_candidate', s)] for s in up_seeds],
                       old_sc['independent_target_candidate'])
    ind_hi = condition('independent_upper', up_rows, up_cfg['scenarios'][0])
    ind_banks = [list(range(0, 128)), list(range(128, 256))]
    # sensitivity module: seed order = cfg seeds = book A + book B
    s_seeds = s_cfg['seeds']
    sens = {sid: condition(sid, [s_table[(sid, s)] for s in s_seeds], s_sc[sid]) for sid in s_sc}
    sens_banks = [list(range(0, 128)), list(range(128, 256))]
    need([s_seeds[i] for i in sens_banks[0]] == s_book['A'] and [s_seeds[i] for i in sens_banks[1]] == s_book['B'], 'sens banks')
    # old 512 holdout (S3): seed order = A + B (256 each)
    o_seeds = old_cfg['seeds']
    old = {sid: condition('old_full_' + sid, [old_table[(sid, s)] for s in o_seeds], old_sc[sid]) for sid in old_sc}
    old_bank_idx = [list(range(0, 256)), list(range(256, 512))]

    # ---------------- reproduction gate
    rs, _ = load_json_bytes(RS_REPORT)
    srep, _ = load_json_bytes(SENS_REPORT)
    hrep, _ = load_json_bytes(HOLDOUT_REPORT)
    rs_entries = {e['id']: e for e in rs['scenarios']}
    s_entries = {e['scenario']['id']: e for e in srep['scenarios']}
    h_entries = {e['scenario']['id']: e for e in hrep['scenarios']}
    gate = {}

    def gate_check(label, cond, arch, bank_idx, arch_banks, counts=None, scenario_arch=None):
        res = dict(archived={}, recomputed={}, rel_err={})
        for key in ('control_rate', 'policy_rate', 'F'):
            res['archived'][key] = arch[key]; res['recomputed'][key] = cond[key]
            res['rel_err'][key] = rel_err(cond[key], arch[key])
        for bname, idx in zip(('A', 'B'), bank_idx):
            bp = bank_point(cond, idx)
            ab = arch_banks[bname]
            need(ab['n'] == len(idx), f'bank size {label}/{bname}')
            for key in ('control_rate', 'policy_rate', 'F'):
                k = f'bank{bname}_{key}'
                res['archived'][k] = ab[key]; res['recomputed'][k] = bp[key]
                res['rel_err'][k] = rel_err(bp[key], ab[key])
        if counts is not None:
            exact_counts = (cond['C0'] == counts[0] and cond['C1'] == counts[2])
            res['counts_exact'] = bool(exact_counts)
            res['rel_err']['control_patient_days'] = rel_err(cond['D0'], counts[1])
            res['rel_err']['policy_patient_days'] = rel_err(cond['D1'], counts[3])
        if scenario_arch is not None:
            res['scenario_matches_archived'] = scenario_arch == cond['scenario']
        res['max_rel_err'] = max(res['rel_err'].values())
        res['passed'] = res['max_rel_err'] <= GATE_TOL and res.get('counts_exact', True) and res.get('scenario_matches_archived', True)
        gate[label] = res

    for label, cond, rid, bidx in (('independent_center_matched256', ind_lo, 'independent_center_matched256', ind_banks),
                                   ('independent_upper', ind_hi, 'independent_upper', ind_banks)):
        e = rs_entries[rid]
        gate_check(label, cond, e, bidx, e['banks'], (e['control_count'], e['control_patient_days'],
                                                       e['policy_count'], e['policy_patient_days']), e['scenario'])
    for sid, cond in old.items():
        e = rs_entries['old_full_' + sid]
        gate_check('old_full_' + sid, cond, e, old_bank_idx, e['banks'],
                   (e['control_count'], e['control_patient_days'], e['policy_count'], e['policy_patient_days']), e['scenario'])
        # Also against the original holdout report itself (pooled only; its own two-bank diagnostics too).
        h = h_entries[sid]
        hb = hrep['two_bank_diagnostics']['banks']
        arch = dict(control_rate=h['arms']['baseline']['rate_per_1000_patient_days'],
                    policy_rate=h['arms']['policy']['rate_per_1000_patient_days'], F=h['rate_prevented_fraction'])
        arch_banks = {b: dict(n=hb[b]['n'], control_rate=hb[b]['conditions'][sid]['R0']['estimate'],
                              policy_rate=hb[b]['conditions'][sid]['R1']['estimate'],
                              F=hb[b]['conditions'][sid]['F']['estimate']) for b in ('A', 'B')}
        need(hb['A']['seeds'] == old_banks['A'] and hb['B']['seeds'] == old_banks['B'], 'holdout report bank seeds')
        gate_check('holdout_report_' + sid, cond, arch, old_bank_idx, arch_banks,
                   (h['arms']['baseline']['acquired_total'], h['arms']['baseline']['patient_days_total'],
                    h['arms']['policy']['acquired_total'], h['arms']['policy']['patient_days_total']), h['scenario'])
    for sid, cond in sens.items():
        e = s_entries[sid]
        gate_check(sid, cond, e, sens_banks, e['banks'], None, e['scenario'])
    gate_passed = all(g['passed'] for g in gate.values())
    gate_max = max(g['max_rel_err'] for g in gate.values())
    for k, g in gate.items():
        log(f'gate {k}: passed={g["passed"]} max_rel_err={g["max_rel_err"]:.3e} '
            f'B={g["recomputed"]["control_rate"]:.15g} P={g["recomputed"]["policy_rate"]:.15g} F={g["recomputed"]["F"]:.15g}')
    hash_ok = all(v['passed'] for v in hashes.values())
    log('reproduction gate passed:', gate_passed and hash_ok, 'max rel err', gate_max)

    results = dict(
        schema='common_baseline_impl_b_v1', protocol=str(ABM / 'research/common_baseline_2026-10-16/PROTOCOL.md'),
        protocol_sha256=sha256_file(ABM / 'research/common_baseline_2026-10-16/PROTOCOL.md'),
        script_sha256=sha256_file(Path(__file__)), numpy=np.__version__, python=sys.version,
        label='post hoc, secondary analysis', hash_checks=hashes,
        reproduction_gate=dict(passed=bool(gate_passed and hash_ok), tolerance=GATE_TOL, max_rel_err=gate_max,
                               conditions=gate))
    if not (gate_passed and hash_ok):
        (HERE / 'results.json').write_text(json.dumps(results, indent=2, sort_keys=True, default=float) + '\n')
        log('STOP: reproduction gate failed; no estimates computed.')
        return 1

    # ---------------- sanity: model_kwargs and scenario identity across bracket sets
    def kwargs_diff(conds, matched):
        diffs = {}
        if matched:
            for i in range(len(conds[0]['rows'])):
                ks = [c['rows'][i]['kwargs'] for c in conds]
                need(len({c['rows'][i]['seed'] for c in conds}) == 1, 'matched seed mismatch')
                for k in set().union(*ks):
                    vals = [kw.get(k, '<absent>') for kw in ks]
                    if any(v != vals[0] for v in vals):
                        diffs.setdefault(k, set()).add(json.dumps(vals))
        else:
            for k in set().union(*(r['kwargs'] for c in conds for r in c['rows'])):
                if k == 'seed':
                    continue
                vals = [sorted({json.dumps(r['kwargs'].get(k, '<absent>')) for r in c['rows']}) for c in conds]
                if any(len(v) != 1 for v in vals):
                    diffs.setdefault(k, set()).add('varies within a condition')
                if len({v[0] for v in vals}) > 1:
                    diffs.setdefault(k, set()).add(json.dumps([json.loads(v[0]) for v in vals]))
        scen_diff = sorted({k for k in set().union(*(c['scenario'] for c in conds))
                            if len({json.dumps(c['scenario'].get(k)) for c in conds}) > 1})
        return dict(kwargs_keys_differing=sorted(diffs), only_prob_transmission=sorted(diffs) == ['prob_transmission'],
                    differing_values={k: sorted(v)[:3] for k, v in diffs.items()}, scenario_fields_differing=scen_diff)

    sanity = {
        'independent_primary(matched256 vs upper, per seed)': kwargs_diff([ind_lo, ind_hi], True),
        'clustered_primary(center vs high, per seed)': kwargs_diff([sens['clustered_center'], sens['clustered_high']], True),
        'imperfect_primary(center vs high, per seed)': kwargs_diff([sens['imperfect_center'], sens['imperfect_high']], True),
        'clustered_S2(low,center,high per seed)': kwargs_diff([sens['clustered_low'], sens['clustered_center'], sens['clustered_high']], True),
        'imperfect_S2(low,center,high per seed)': kwargs_diff([sens['imperfect_low'], sens['imperfect_center'], sens['imperfect_high']], True),
        'independent_S3(old512 vs upper, unmatched, excluding seed)': kwargs_diff([old['independent_target_candidate'], ind_hi], False),
        'clustered_S3(old512 vs sens high, unmatched, excluding seed)': kwargs_diff([old['persistent_target_candidate'], sens['clustered_high']], False),
        'imperfect_S3(old512 vs sens high, unmatched, excluding seed)': kwargs_diff([old['imperfect_target_candidate'], sens['imperfect_high']], False),
    }
    # execution environment (not model input, recorded for transparency)
    def env_summary(c):
        out = {}
        for r in c['rows']:
            key = r['execution'].get('platform', '?') + ' | python ' + r['execution'].get('python', '?').split()[0]
            out[key] = out.get(key, 0) + 1
        return out
    execution_env = {c['name']: env_summary(c) for c in [ind_lo, ind_hi, sens['clustered_center'], sens['clustered_high'],
                                                       sens['imperfect_center'], sens['imperfect_high'],
                                                       old['persistent_target_candidate'], old['imperfect_target_candidate']]}
    for k, v in sanity.items():
        log('sanity', k, json.dumps(v))
    log('execution environments', json.dumps(execution_env))

    # ---------------- brackets
    brackets = {
        'independent': (ind_lo, ind_hi, ind_banks, 'independent_module'),
        'clustered': (sens['clustered_center'], sens['clustered_high'], sens_banks, 'sensitivity_module'),
        'imperfect': (sens['imperfect_center'], sens['imperfect_high'], sens_banks, 'sensitivity_module'),
    }

    def point_estimates(lo, hi, target):
        r = interp(lo['control_rate'], lo['policy_rate'], hi['control_rate'], hi['policy_rate'], lo['p'], hi['p'], target)
        # log-form cross-check of the protocol formula
        w = float(r['w'])
        logform = 1 - math.exp((1 - w) * math.log(lo['policy_rate']) + w * math.log(hi['policy_rate']) - math.log(target))
        need(abs(logform - float(r['Fstar'])) < 1e-13, 'F* formula cross-check')
        return {k: (float(v) if np.ndim(v) == 0 else v) for k, v in r.items()}

    # ---------------- bootstrap (primary streams, then S3 stream)
    rng = np.random.default_rng(BOOT_SEED)
    order = []
    V_ind = np.stack([ind_lo['V'], ind_hi['V']], axis=1)  # (256, 2, 4)
    sens_ids = ['clustered_low', 'clustered_center', 'clustered_high', 'imperfect_low', 'imperfect_center', 'imperfect_high']
    V_sens = np.stack([sens[s]['V'] for s in sens_ids], axis=1)  # (256, 6, 4)
    old_ids = ['independent_target_candidate', 'persistent_target_candidate', 'imperfect_target_candidate']
    V_old = np.stack([old[s]['V'] for s in old_ids], axis=1)  # (512, 3, 4)
    S_ind = boot_module_sums(rng, ind_banks, V_ind, DRAWS, order, 'independent_module (A=0..127, B=128..255 of upper seed order)')
    S_sens = boot_module_sums(rng, sens_banks, V_sens, DRAWS, order, 'sensitivity_module (A=0..127, B=128..255 of config seed order), all 6 conditions synchronous')
    S_old = boot_module_sums(rng, old_bank_idx, V_old, DRAWS, order, 'S3 old 512 holdout (A=0..255, B=256..511), 3 conditions synchronous')
    log('bootstrap sums done', f'{time.time() - t0:.1f}s')
    Bi, Pi = rates_from_sums(S_ind)
    Bs, Ps = rates_from_sums(S_sens)
    Bo, Po = rates_from_sums(S_old)
    zero_counts = dict(independent_module=int(((S_ind[..., 0] == 0) | (S_ind[..., 2] == 0)).any(axis=1).sum()),
                       sensitivity_module=int(((S_sens[..., 0] == 0) | (S_sens[..., 2] == 0)).any(axis=1).sum()),
                       old_holdout=int(((S_old[..., 0] == 0) | (S_old[..., 2] == 0)).any(axis=1).sum()))
    sidx = {s: i for i, s in enumerate(sens_ids)}
    draws_bp = {
        'independent': (Bi[:, 0], Pi[:, 0], Bi[:, 1], Pi[:, 1]),
        'clustered': (Bs[:, sidx['clustered_center']], Ps[:, sidx['clustered_center']], Bs[:, sidx['clustered_high']], Ps[:, sidx['clustered_high']]),
        'imperfect': (Bs[:, sidx['imperfect_center']], Ps[:, sidx['imperfect_center']], Bs[:, sidx['imperfect_high']], Ps[:, sidx['imperfect_high']]),
    }
    oidx = {'independent': 0, 'clustered': 1, 'imperfect': 2}
    draws_s3 = {h: (Bo[:, oidx[h]], Po[:, oidx[h]], draws_bp[h][2], draws_bp[h][3]) for h in PROFILES}

    def analyse(target, draw_source, lo_override=None):
        out, dr = {}, {}
        for h in PROFILES:
            lo, hi, _, module = brackets[h]
            if lo_override is not None:
                lo = lo_override[h]
            pt = point_estimates(lo, hi, target)
            bd = interp(*draw_source[h], lo['p'], hi['p'], target)
            dr[h] = bd
            undefined = int(bd['undefined'].sum())
            ci = pct(bd['Fstar'], CI_PRIMARY)
            out[h] = dict(profile=h, module=module, target=target, lower=lo['name'], upper=hi['name'],
                          p_lo=lo['p'], p_hi=hi['p'], B_lo=lo['control_rate'], B_hi=hi['control_rate'],
                          P_lo=lo['policy_rate'], P_hi=hi['policy_rate'], F_lo=lo['F'], F_hi=hi['F'],
                          n_lo=len(lo['rows']), n_hi=len(hi['rows']),
                          w=pt['w'], p_star=pt['pstar'], F_star=pt['Fstar'],
                          point_extrapolates=bool(pt['extrapolated']),
                          ci95=ci['interval'], ci=ci,
                          frac_extrapolated=float(bd['extrapolated'].mean()),
                          n_extrapolated=int(bd['extrapolated'].sum()),
                          undefined_draws=undefined,
                          w_ci95=pct(bd['w'], .95)['interval'], p_star_ci95=pct(bd['pstar'], .95)['interval'],
                          B_lo_ci95=pct(draw_source[h][0], .95)['interval'], B_hi_ci95=pct(draw_source[h][2], .95)['interval'],
                          S1_arith_F_in_lnB=dict(F_star=pt['F_s1'], ci95=pct(bd['F_s1'], .95)['interval'],
                                                 undefined_draws=int((~np.isfinite(bd['F_s1'])).sum())))
        contrasts = []
        for a, b in PAIRS:
            d = dr[a]['Fstar'] - dr[b]['Fstar']
            ci = pct(d, CI_CONTRAST)
            est = out[a]['F_star'] - out[b]['F_star']
            contrasts.append(dict(pair=f'{a}-{b}', estimate=est, level=CI_CONTRAST, ci=ci['interval'],
                                  undefined_draws=ci['undefined_draws'],
                                  resolved=bool(ci['interval'] is not None and (ci['interval'][0] > 0 or ci['interval'][1] < 0))))
        return out, contrasts

    # ---------------- bootstrap mechanics check: per-condition 95% intervals vs archived ones
    # (same bank structure for the independent and sensitivity modules; Monte Carlo differences only).
    def cond_ci(B, P):
        return dict(control_rate=pct(B, .95)['interval'], policy_rate=pct(P, .95)['interval'],
                    F=pct(1 - P / B, .95)['interval'])
    boot_check = {}
    arch_ci = {'independent_center_matched256': (rs_entries['independent_center_matched256'], 'rs'),
               'independent_upper': (rs_entries['independent_upper'], 'rs')}
    for sid in sens_ids:
        arch_ci[sid] = (s_entries[sid], 'sens')
    for sid in old_ids:
        arch_ci['old_full_' + sid] = (rs_entries['old_full_' + sid], 'rs')
    mine = {'independent_center_matched256': cond_ci(Bi[:, 0], Pi[:, 0]), 'independent_upper': cond_ci(Bi[:, 1], Pi[:, 1])}
    for sid in sens_ids:
        mine[sid] = cond_ci(Bs[:, sidx[sid]], Ps[:, sidx[sid]])
    for j, sid in enumerate(old_ids):
        mine['old_full_' + sid] = cond_ci(Bo[:, j], Po[:, j])
    for name, (e, kind) in arch_ci.items():
        arch = {k: e[k + '95']['interval'] for k in ('control_rate', 'policy_rate', 'F')}
        diff = max(abs(a - b) for k in arch for a, b in zip(arch[k], mine[name][k]))
        f_diff = max(abs(a - b) for a, b in zip(arch['F'], mine[name]['F']))
        boot_check[name] = dict(impl_b=mine[name], archived=arch, max_abs_bound_diff=diff, F_max_abs_bound_diff=f_diff,
                                note=('archived old_full_* used four 128-seed strata (rule-surface); impl_b S3 uses two 256-seed banks'
                                      if name.startswith('old_full_') else 'same bank structure; Monte Carlo difference only'))
        log(f'bootstrap check {name}: max |bound diff| (B,P,F) = {diff:.4g}; F bounds diff = {f_diff:.4g}')

    primary, contrasts = analyse(TARGET, draws_bp)
    for h in PROFILES:
        r = primary[h]
        log(f'PRIMARY {h}: w={r["w"]:.6f} p*={r["p_star"]:.8f} F*={r["F_star"]:.6f} CI95={r["ci95"]} '
            f'frac_extrap={r["frac_extrapolated"]:.5f} undefined={r["undefined_draws"]}')
    for c in contrasts:
        log(f'CONTRAST {c["pair"]}: {c["estimate"]:.6f} CI{100 * c["level"]:.3f}={c["ci"]} resolved={c["resolved"]}')

    # ---------------- S1 (already computed alongside primary, same draws)
    S1 = {h: dict(F_star=primary[h]['S1_arith_F_in_lnB']['F_star'], ci95=primary[h]['S1_arith_F_in_lnB']['ci95'],
                  undefined_draws=primary[h]['S1_arith_F_in_lnB']['undefined_draws'],
                  minus_primary=primary[h]['S1_arith_F_in_lnB']['F_star'] - primary[h]['F_star']) for h in PROFILES}
    for h in PROFILES:
        log(f'S1 {h}: F*={S1[h]["F_star"]:.6f} CI95={S1[h]["ci95"]} (S1-primary={S1[h]["minus_primary"]:.6f})')

    # ---------------- S2 quadratic in ln p (clustered, imperfect)
    S2 = {}
    for h in ('clustered', 'imperfect'):
        ids3 = [f'{h}_low', f'{h}_center', f'{h}_high']
        x = np.array([math.log(sens[s]['p']) for s in ids3])
        yB = np.array([math.log(sens[s]['control_rate']) for s in ids3])
        yP = np.array([math.log(sens[s]['policy_rate']) for s in ids3])
        x_ref = math.log(primary[h]['p_star'])
        pt = quad_eval_root(x, yB[None, :], yP[None, :], TARGET, np.array([x_ref]))
        # point via an independent route: numpy polyfit on raw ln p, then np.roots
        cb = np.polyfit(x, yB, 2); cp = np.polyfit(x, yP, 2)
        roots = np.roots(cb - np.array([0, 0, math.log(TARGET)]))
        real = [float(r.real) for r in roots if abs(r.imag) < 1e-12 and x[0] - 1e-12 <= r.real <= x[2] + 1e-12]
        need(len(real) >= 1, f'S2 no in-range root at point estimate for {h}')
        xr = min(real, key=lambda r: abs(r - x_ref))
        F_poly = 1 - math.exp(float(np.polyval(cp, xr)) - math.log(TARGET))
        need(abs(F_poly - float(pt['Fstar'][0])) < 1e-9, f'S2 cross-check {h}: {F_poly} vs {pt["Fstar"][0]}')
        yBd = np.log(np.stack([Bs[:, sidx[s]] for s in ids3], axis=1))
        yPd = np.log(np.stack([Ps[:, sidx[s]] for s in ids3], axis=1))
        primary_draws = interp(*draws_bp[h], brackets[h][0]['p'], brackets[h][1]['p'], TARGET)
        x_ref_d = np.log(primary_draws['pstar'])
        bd = quad_eval_root(x, yBd, yPd, TARGET, x_ref_d)
        S2[h] = dict(F_star=float(pt['Fstar'][0]), p_star=float(pt['pstar'][0]),
                     F_star_polyfit_crosscheck=F_poly, in_range_roots_at_point=len(real),
                     quadratic_coef_lnB_in_lnp=[float(v) for v in cb], quadratic_coef_lnP_in_lnp=[float(v) for v in cp],
                     minus_primary=float(pt['Fstar'][0]) - primary[h]['F_star'],
                     ci95=pct(bd['Fstar'], .95)['interval'],
                     undefined_draws=int(bd['undefined'].sum()),
                     draws_with_two_in_range_roots=int(bd['two_in_range'].sum()),
                     root_rule='real roots of q_B(ln p)=ln 3 within [ln p_low, ln p_high]; if two, the one nearest the same-draw primary log-linear ln p*; none -> undefined')
        log(f'S2 {h}: F*={S2[h]["F_star"]:.6f} p*={S2[h]["p_star"]:.8f} CI95={S2[h]["ci95"]} '
            f'(S2-primary={S2[h]["minus_primary"]:.6f}) undefined={S2[h]["undefined_draws"]} two_roots={S2[h]["draws_with_two_in_range_roots"]}')

    # ---------------- S3: original 512-seed holdout center as lower bracket
    s3_lo = {'independent': old['independent_target_candidate'], 'clustered': old['persistent_target_candidate'],
             'imperfect': old['imperfect_target_candidate']}
    S3_res, S3_con = analyse(TARGET, draws_s3, lo_override=s3_lo)
    S3 = {h: {k: S3_res[h][k] for k in ('lower', 'upper', 'n_lo', 'n_hi', 'p_lo', 'p_hi', 'B_lo', 'B_hi', 'F_lo', 'F_hi',
                                         'w', 'p_star', 'F_star', 'ci95', 'frac_extrapolated', 'undefined_draws')}
          for h in PROFILES}
    for h in PROFILES:
        S3[h]['minus_primary'] = S3[h]['F_star'] - primary[h]['F_star']
        log(f'S3 {h}: B_lo={S3[h]["B_lo"]:.6f} F_lo={S3[h]["F_lo"]:.6f} w={S3[h]["w"]:.6f} F*={S3[h]["F_star"]:.6f} '
            f'CI95={S3[h]["ci95"]} frac_extrap={S3[h]["frac_extrapolated"]:.5f} undefined={S3[h]["undefined_draws"]}')
    for c in S3_con:
        log(f'S3 contrast {c["pair"]}: {c["estimate"]:.6f} CI{100 * c["level"]:.3f}={c["ci"]} excl0={c["resolved"]}')

    # S3 convention check (descriptive only; drawn AFTER all other streams so they are unaffected):
    # the archived holdout primary intervals used unstratified 512-seed resampling, and the rule-surface
    # reanalysis used four 128-seed strata. Report F* intervals under those alternatives too.
    S3_alt = {}
    for label, banks_alt in (('unstratified_512', [list(range(512))]),
                             ('four_strata_128_rule_surface', [list(range(a, a + 128)) for a in (0, 128, 256, 384)])):
        S_alt = boot_module_sums(rng, banks_alt, V_old, DRAWS, order, f'S3 convention check: old 512 holdout, {label}')
        Ba, Pa = rates_from_sums(S_alt)
        alt_src = {h: (Ba[:, oidx[h]], Pa[:, oidx[h]], draws_bp[h][2], draws_bp[h][3]) for h in PROFILES}
        ra, ca = analyse(TARGET, alt_src, lo_override=s3_lo)
        S3_alt[label] = dict(profiles={h: dict(F_star=ra[h]['F_star'], ci95=ra[h]['ci95'],
                                               frac_extrapolated=ra[h]['frac_extrapolated'],
                                               undefined_draws=ra[h]['undefined_draws']) for h in PROFILES},
                             contrasts=[dict(pair=c['pair'], ci=c['ci'], resolved=c['resolved']) for c in ca])
        for h in PROFILES:
            log(f'S3-alt {label} {h}: CI95={ra[h]["ci95"]} frac_extrap={ra[h]["frac_extrapolated"]:.5f}')

    # ---------------- S4: targets 2.7 and 3.3
    S4 = {}
    for T in (2.7, 3.3):
        r, cc = analyse(T, draws_bp)
        S4[str(T)] = dict(profiles={h: {k: r[h][k] for k in ('w', 'p_star', 'F_star', 'ci95', 'frac_extrapolated',
                                                              'undefined_draws', 'point_extrapolates')} for h in PROFILES},
                          contrasts=cc)
        for h in PROFILES:
            log(f'S4 T={T} {h}: w={r[h]["w"]:.6f} F*={r[h]["F_star"]:.6f} CI95={r[h]["ci95"]} '
                f'frac_extrap={r[h]["frac_extrapolated"]:.5f} point_extrapolates={r[h]["point_extrapolates"]}')
        for c in cc:
            log(f'S4 T={T} contrast {c["pair"]}: {c["estimate"]:.6f} CI{100 * c["level"]:.3f}={c["ci"]} excl0={c["resolved"]}')

    # ---------------- assemble results
    for h in PROFILES:
        primary[h].pop('ci', None)
    results.update(
        bank_definitions=dict(
            independent_module=dict(
                seed_list='rule_surface configs/independent_upper.json seeds = old holdout ordered banks A[:128] + B[:128] '
                          '(old_input_index.json independent_upper_comparison_seeds; seedbook seed_blocks.holdout)',
                bank_A='positions 0-127 (old bank A core: first 128 of ordered holdout bank A)',
                bank_B='positions 128-255 (old bank B core: first 128 of ordered holdout bank B)',
                lower_bracket='independent_target_candidate records of the original holdout (cluster+local sites) at those 256 seeds, '
                              'reordered to the upper seed order; = rule-surface report independent_center_matched256 '
                              '(rule_surface analyze.py: core = rows 0..127 and 256..383 of old A+B order)',
                verified_by='bank A/B control, policy, F reproduced against rule-surface report banks to <=1e-12'),
            sensitivity_module=dict(
                seed_list='abm/research/author_notes_2026-10-10/configs/sensitivity.json seeds (256) = seedbook.json["sensitivity"] A + B',
                bank_A='positions 0-127 (seedbook sensitivity A)', bank_B='positions 128-255 (seedbook sensitivity B)',
                verified_by='bank A/B values of all 6 conditions reproduced against sensitivity report to <=1e-12'),
            S3_old_holdout=dict(
                seed_list='full_icu_confirmation_2026-10-08 configs/holdout_all.json seeds (512) = seedbook seed_blocks.holdout A + B',
                convention='original holdout two-bank convention: banks A and B of 256 seeds each (seedbook seed_blocks.holdout, '
                           'two_bank_diagnostics of the holdout report); resampled independently within each bank, synchronously '
                           'across the three profile centers; independent of the upper-bracket streams. (The archived holdout '
                           'primary intervals used unstratified 512-seed resampling; the rule-surface reanalysis used four 128 '
                           'strata; neither is used here.)',
                verified_by='bank A/B values reproduced against rule-surface old_full_* banks and holdout report two_bank_diagnostics to <=1e-12')),
        bootstrap=dict(draws=DRAWS, seed=BOOT_SEED, generator='numpy.random.default_rng(20261016), PCG64',
                       method='whole-seed resampling with replacement: per bank, rng.integers(0,n_bank,(chunk,n_bank)) -> '
                              'multiplicity counts; sums = counts @ per-seed (C0,D0,C1,D1); rates as ratio of sums',
                       stream_order=order, chunk=CHUNK,
                       quantiles='numpy.quantile method=linear on all draws (no undefined draws occurred)',
                       draws_with_any_zero_count=zero_counts),
        brackets={h: dict(lower=brackets[h][0]['name'], upper=brackets[h][1]['name'], module=brackets[h][3])
                  for h in PROFILES},
        primary=primary, contrasts=contrasts,
        sensitivity=dict(S1=S1, S2=S2, S3=dict(profiles=S3, contrasts_descriptive=S3_con,
                                                convention='old holdout two banks A/B of 256 (seedbook seed_blocks.holdout)',
                                                alternative_conventions_descriptive=S3_alt), S4=S4),
        sanity_model_kwargs=sanity, execution_environments=execution_env,
        bootstrap_mechanics_check_vs_archived_condition_intervals=boot_check,
        runtime_seconds=time.time() - t0)
    (HERE / 'results.json').write_text(json.dumps(results, indent=2, sort_keys=True, allow_nan=False) + '\n')
    log('wrote', HERE / 'results.json', f'{time.time() - t0:.1f}s')
    return 0


if __name__ == '__main__':
    sys.exit(main())
