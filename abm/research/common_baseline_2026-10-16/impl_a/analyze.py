#!/usr/bin/env python
"""Common-baseline prevention (post hoc, secondary): implementation A.

Implements PROTOCOL.md (common_baseline_2026-10-16): reproduction gate, primary
log-log interpolation estimator at B*=3, whole-seed bank-stratified bootstrap
(50,000 draws, default_rng(20261016)), three pairwise contrasts, and descriptive
sensitivity analyses S1-S4, plus a model_kwargs sanity check. D1 is not done here.

Written independently; archived analysis code was read only for conventions
(bank membership, record fields, ratio-of-sums), never imported.

All inputs are read-only. Outputs: results.json next to this script.
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
DM = ABM.parent

# ---------------------------------------------------------------- input paths
RS_DIR = ABM / 'research/rule_surface_followup_2026-10-11'
RS_REPORT = ABM / 'research/rule_surface_results_2026-10-12/completed_analysis/report.json'
RS_RELEASE = ABM / 'research/rule_surface_results_2026-10-12/RELEASE.json'
UPPER_ROOT = (DM / 'research/rule_surface_followup_2026-10-11/collections/independent_upper_20261012_0125_utc'
              / 'restored/repo/results/rule_surface_followup_2026-10-11/independent_upper')
HOLD_CFG = ABM / 'research/full_icu_confirmation_2026-10-08/configs/holdout_all.json'
HOLD_SEEDBOOK = ABM / 'research/full_icu_confirmation_2026-10-08/seedbook.json'
HOLD_REPORT = ABM / 'research/full_icu_confirmation_2026-10-08/holdout_report/report.json'
AN_DIR = ABM / 'research/author_notes_2026-10-10'
SENS_ROOT = DM / 'research/author_notes_2026-10-10/collections/sensitivity_20261011_101348_utc/results/sensitivity'
SENS_SHA_MANIFEST = DM / 'research/author_notes_2026-10-10/collections/sensitivity_20261011_101348_utc/SHA256_MANIFEST.json'
SENS_REPORT_DIR = DM / 'research/author_notes_2026-10-10/analysis_2026-10-11/sensitivity_replay_20261011_1015_utc'
SENS_REPORT = SENS_REPORT_DIR / 'report.json'

DRAWS = 50000
BOOT_SEED = 20261016
TARGET = 3.0
LEVEL_F = 0.95
LEVEL_C = 1 - 0.05 / 3
REL_TOL = 1e-12

LOG = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    LOG.append(s)
    print(s, flush=True)


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def load(path):
    with open(path, 'rb') as f:
        return json.loads(f.read())


class GateFailure(Exception):
    pass


HASH_CHECKS = []      # (label, n_checked, n_mismatch)
HASH_MISMATCHES = []


SINGLE_HASH_CHECKS = []


def check_hash(label, path, expected, single=True):
    got = sha(path)
    if got != expected:
        HASH_MISMATCHES.append(dict(label=label, path=str(path), expected=expected, got=got))
    if single:
        SINGLE_HASH_CHECKS.append(dict(label=label, path=str(path), sha256=got, match=got == expected))
        log(f'  hash check {label}: {"match" if got == expected else "MISMATCH"}')
    return got == expected


def check_many(label, pairs):
    bad = 0
    for path, expected in pairs:
        if not check_hash(label, path, expected, single=False):
            bad += 1
    HASH_CHECKS.append(dict(label=label, n_files=len(pairs), mismatches=bad))
    log(f'  hash check {label}: {len(pairs)} files, {bad} mismatches')


# ------------------------------------------------------------- record reading
def record_values(row):
    """(control acq, control pd, policy acq, policy pd) from the 56-day observation."""
    out = []
    for arm, h in (('baseline', 0.9), ('policy', 0.95)):
        a = row['arms'][arm]
        obs = a['observation']
        assert a['h'] == h, ('arm h', arm, a['h'])
        assert obs['days'] == 56, obs['days']
        assert obs['ticks'] == 145152, obs['ticks']
        c, d = obs['acquired_colonizations'], obs['patient_days']
        assert type(c) is int and c >= 0 and d > 0
        # field-semantics checks
        assert math.isclose(obs['acquired_per_1000_patient_days'], 1000 * c / d, rel_tol=1e-12, abs_tol=0)
        blocks = a['observation_blocks']
        assert len(blocks) == 56 and sum(b['acquired_colonizations'] for b in blocks) == c
        assert math.isclose(math.fsum(b['patient_days'] for b in blocks), d, rel_tol=1e-9)
        out += [c, d]
    assert row['status'] == 'complete'
    assert row['distinct_observation_trajectories'] == 2
    return out


def kw_strip(kw, drop):
    return {k: v for k, v in kw.items() if k not in drop}


# ------------------------------------------------------------------- loaders
def load_independent_upper(rs_report):
    cfg = load(RS_DIR / 'configs/independent_upper.json')
    seeds = cfg['seeds']
    assert len(cfg['scenarios']) == 1
    binding = rs_report['new_bank_bindings']['independent_upper']
    log('Independent upper bracket:', UPPER_ROOT)
    check_hash('independent_upper manifest.json', UPPER_ROOT / 'manifest.json', binding['manifest_sha256'])
    check_hash('independent_upper study_binding.json', UPPER_ROOT / 'study_binding.json', binding['binding_sha256'])
    files = sorted((UPPER_ROOT / 'completed').glob('job_*.json'))
    assert len(files) == 256 == len(binding['input_sha256'])
    assert {'completed/' + f.name for f in files} == set(binding['input_sha256'])
    check_many('independent_upper records vs rule-surface report input_sha256',
               [(UPPER_ROOT / k, v) for k, v in binding['input_sha256'].items()])
    manifest = load(UPPER_ROOT / 'manifest.json')
    assert manifest['config'] == cfg
    vals = {}
    meta = {}
    for f in files:
        row = load(f)
        assert row['identity'] == manifest['identity']
        assert row['scenario'] == cfg['scenarios'][0]
        s = row['seed']
        assert s not in vals
        assert row['model_kwargs']['seed'] == s
        assert row['model_kwargs']['prob_transmission'] == cfg['scenarios'][0]['p']
        vals[s] = record_values(row)
        meta[s] = dict(kw=row['model_kwargs'], scenario=row['scenario'])
    assert set(vals) == set(seeds)
    return cfg, seeds, vals, meta


def load_holdout():
    index = load(RS_DIR / 'old_input_index.json')
    cfg = load(HOLD_CFG)
    seedbook = load(HOLD_SEEDBOOK)
    book = seedbook['seed_blocks']['holdout']
    assert book == index['ordered_statistical_banks']
    assert book['A'] + book['B'] == cfg['seeds'] and len(book['A']) == len(book['B']) == 256
    assert index['independent_upper_comparison_seeds'] == book['A'][:128] + book['B'][:128]
    # manifests and records listed in the frozen old-input index
    check_many('holdout manifests vs old_input_index.json',
               [(ABM / k, v) for k, v in index['manifests'].items()])
    check_many('holdout records (1536) vs old_input_index.json',
               [(ABM / k, v) for k, v in index['records'].items()])
    # second, independent listing: the original holdout report's per-site bank hashes
    hold_report = load(HOLD_REPORT)
    pairs = []
    for b in hold_report['banks']:
        site = None
        for name in ('holdout_cluster', 'holdout_local'):
            if sha(ABM / 'results/full_icu_confirmation_2026-10-08' / name / 'manifest.json') == b['manifest_sha256']:
                site = name
        assert site is not None
        for fn, h in b['completed_records_sha256'].items():
            pairs.append((ABM / 'results/full_icu_confirmation_2026-10-08' / site / 'completed' / fn, h))
    assert len(pairs) == 1536
    check_many('holdout records (1536) vs original holdout report banks[].completed_records_sha256', pairs)
    scen_ids = [s['id'] for s in cfg['scenarios']]
    vals = {sid: {} for sid in scen_ids}
    meta = {sid: {} for sid in scen_ids}
    for rel, md in index['record_metadata'].items():
        row = load(ABM / rel)
        sid = row['scenario']['id']
        assert sid == md['scenario'] and row['seed'] == md['seed']
        assert row['scenario'] == cfg['scenarios'][scen_ids.index(sid)]
        s = row['seed']
        si = cfg['seeds'].index(s)
        assert md['statistical_bank'] == ('A' if si < 256 else 'B')
        assert s not in vals[sid]
        assert row['model_kwargs']['seed'] == s
        assert row['model_kwargs']['prob_transmission'] == row['scenario']['p']
        vals[sid][s] = record_values(row)
        meta[sid][s] = dict(kw=row['model_kwargs'], scenario=row['scenario'])
    for sid in scen_ids:
        assert set(vals[sid]) == set(cfg['seeds'])
    return cfg, book, vals, meta, hold_report, index


def load_sensitivity(sens_report):
    cfg_path = AN_DIR / 'configs/sensitivity.json'
    cfg = load(cfg_path)
    seedbook = load(AN_DIR / 'seedbook.json')['sensitivity']
    freeze = load(AN_DIR / 'freeze.json')
    for rel in ('research/author_notes_2026-10-10/configs/sensitivity.json',
                'research/author_notes_2026-10-10/seedbook.json',
                'research/author_notes_2026-10-10/analyze_positive.py'):
        check_hash('author_notes freeze ' + rel, ABM / rel, freeze['files'][rel])
    HASH_CHECKS.append(dict(label='author_notes freeze.json (config, seedbook, analysis source)', n_files=3,
                            mismatches=sum(1 for m in HASH_MISMATCHES if m['label'].startswith('author_notes freeze'))))
    assert sha(cfg_path) == sens_report['config_sha256']
    assert sha(AN_DIR / 'analyze_positive.py') == sens_report['source_sha256']
    binding = load(SENS_ROOT / 'study_binding.json')
    assert binding['config_sha256'] == sha(cfg_path) and binding['freeze_identity'] == freeze['identity_sha256']
    check_hash('sensitivity manifest.json vs report', SENS_ROOT / 'manifest.json', sens_report['manifest_sha256'])
    shaman = load(SENS_SHA_MANIFEST)
    check_hash('sensitivity manifest.json vs SHA256_MANIFEST', SENS_ROOT / 'manifest.json', shaman['files']['manifest.json']['sha256'])
    check_hash('sensitivity study_binding.json vs SHA256_MANIFEST', SENS_ROOT / 'study_binding.json', shaman['files']['study_binding.json']['sha256'])
    assert seedbook['A'] + seedbook['B'] == cfg['seeds'] and len(seedbook['A']) == len(seedbook['B']) == 128
    files = sorted((SENS_ROOT / 'completed').glob('job_*.json'))
    assert len(files) == 1536 == len(sens_report['input_sha256'])
    assert {f.name for f in files} == set(sens_report['input_sha256'])
    check_many('sensitivity records vs sensitivity report input_sha256',
               [(SENS_ROOT / 'completed' / k, v) for k, v in sens_report['input_sha256'].items()])
    check_many('sensitivity records vs collection SHA256_MANIFEST.json',
               [(SENS_ROOT / k, v['sha256']) for k, v in shaman['files'].items() if k.startswith('completed/')])
    manifest = load(SENS_ROOT / 'manifest.json')
    assert manifest['config'] == cfg
    scen_ids = [s['id'] for s in cfg['scenarios']]
    vals = {sid: {} for sid in scen_ids}
    meta = {sid: {} for sid in scen_ids}
    n = len(cfg['seeds'])
    for f in files:
        row = load(f)
        assert row['identity'] == manifest['identity']
        sid = row['scenario']['id']
        s = row['seed']
        ci, si = divmod(row['job_id'], n)              # cross-check only
        assert scen_ids[ci] == sid and cfg['seeds'][si] == s
        assert row['scenario'] == cfg['scenarios'][ci]
        assert s not in vals[sid]
        assert row['model_kwargs']['seed'] == s
        assert row['model_kwargs']['prob_transmission'] == row['scenario']['p']
        vals[sid][s] = record_values(row)
        meta[sid][s] = dict(kw=row['model_kwargs'], scenario=row['scenario'])
    for sid in scen_ids:
        assert set(vals[sid]) == set(cfg['seeds'])
    return cfg, seedbook, vals, meta


# ---------------------------------------------------------------- estimators
def rates(sums):
    """sums[..., 4] -> B, P, F (ratio of sums, per 1000 patient-days)."""
    B = 1000.0 * sums[..., 0] / sums[..., 1]
    P = 1000.0 * sums[..., 2] / sums[..., 3]
    F = 1.0 - P / B
    return B, P, F


def interp(B_lo, B_hi, P_lo, P_hi, target=TARGET):
    """Primary estimator: log rate linear in log p between brackets."""
    with np.errstate(divide='ignore', invalid='ignore'):
        lbl, lbh = np.log(B_lo), np.log(B_hi)
        w = (math.log(target) - lbl) / (lbh - lbl)
        F_lo = 1.0 - P_lo / B_lo
        F_hi = 1.0 - P_hi / B_hi
        Fs = 1.0 - (1.0 - F_lo) ** (1.0 - w) * (1.0 - F_hi) ** w
    return w, Fs, F_lo, F_hi


def undefined_mask(c_lo, c_hi, B_lo, B_hi):
    """Zero control or policy count in either bracket, or B_hi == B_lo."""
    return ((c_lo[..., 0] <= 0) | (c_lo[..., 2] <= 0) | (c_hi[..., 0] <= 0) | (c_hi[..., 2] <= 0)
            | (B_hi == B_lo) | ~np.isfinite(B_lo) | ~np.isfinite(B_hi))


def pct_interval(x, level, undefined):
    good = x[~undefined & np.isfinite(x)]
    q = [(1 - level) / 2, 1 - (1 - level) / 2]
    lo, hi = np.quantile(good, q)
    return float(lo), float(hi)


def p_star(w, p_lo, p_hi):
    return np.exp((1 - w) * math.log(p_lo) + w * math.log(p_hi))


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    log('common_baseline_2026-10-16 implementation A')
    log('python', sys.version.split()[0], 'numpy', np.__version__)
    log('script sha256', sha(Path(__file__)))

    # ---------------- archived reports and frozen descriptors
    rs_release = load(RS_RELEASE)
    check_hash('rule-surface report.json vs RELEASE.json', RS_REPORT, rs_release['report_sha256'])
    rs_report = load(RS_REPORT)
    receipt = load(SENS_REPORT_DIR / 'replay_receipt.json')
    check_hash('sensitivity report.json vs replay_receipt.json', SENS_REPORT, receipt['replay_report_sha256'])
    sens_report = load(SENS_REPORT)
    rs_freeze = load(RS_DIR / 'freeze.json')
    for rel in ('research/rule_surface_followup_2026-10-11/old_input_index.json',
                'research/rule_surface_followup_2026-10-11/seedbook.json',
                'research/rule_surface_followup_2026-10-11/configs/independent_upper.json',
                'research/rule_surface_followup_2026-10-11/design.json',
                'research/rule_surface_followup_2026-10-11/analyze.py',
                'research/full_icu_confirmation_2026-10-08/configs/holdout_all.json',
                'research/full_icu_confirmation_2026-10-08/seedbook.json'):
        check_hash('rule-surface freeze ' + rel, ABM / rel, rs_freeze['files'][rel])
    HASH_CHECKS.append(dict(label='rule-surface freeze.json (index, seedbooks, configs, design, analysis source)', n_files=7,
                            mismatches=sum(1 for m in HASH_MISMATCHES if m['label'].startswith('rule-surface freeze'))))
    assert sha(RS_DIR / 'old_input_index.json') == rs_report['old_input_index_sha256']
    index = load(RS_DIR / 'old_input_index.json')
    check_hash('holdout report.json vs old_input_index', HOLD_REPORT, index['report']['sha256'])
    check_hash('holdout seedbook vs old_input_index', HOLD_SEEDBOOK, index['seedbook']['sha256'])
    check_hash('holdout config vs old_input_index', HOLD_CFG, index['config']['sha256'])
    rs_seedbook = load(RS_DIR / 'seedbook.json')['independent_upper']

    # ---------------- records
    log('Loading records and verifying hashes ...')
    up_cfg, up_seeds, up_vals, up_meta = load_independent_upper(rs_report)
    h_cfg, h_book, h_vals, h_meta, h_report, _ = load_holdout()
    s_cfg, s_book, s_vals, s_meta = load_sensitivity(sens_report)
    assert rs_seedbook['A'] == h_book['A'][:128] and rs_seedbook['B'] == h_book['B'][:128]
    assert up_seeds == rs_seedbook['A'] + rs_seedbook['B']
    log(f'  hash mismatches total: {len(HASH_MISMATCHES)}')

    # ---------------- arrays on fixed seed axes
    # Independent module axis: rule-surface independent_upper seed order = holdout A[:128] + B[:128]
    ind_seeds = up_seeds
    IND = np.zeros((256, 2, 4))      # [:,0]=center matched256 (lower), [:,1]=independent_upper
    for i, s in enumerate(ind_seeds):
        IND[i, 0] = h_vals['independent_target_candidate'][s]
        IND[i, 1] = up_vals[s]
    ind_banks = [np.arange(0, 128), np.arange(128, 256)]
    # Sensitivity axis: config seed order = seedbook A + B
    s_ids = [sc['id'] for sc in s_cfg['scenarios']]
    SEN = np.zeros((256, 6, 4))
    for j, sid in enumerate(s_ids):
        for i, s in enumerate(s_cfg['seeds']):
            SEN[i, j] = s_vals[sid][s]
    sen_banks = [np.arange(0, 128), np.arange(128, 256)]
    # Holdout axis: holdout config seed order = seed_blocks.holdout A + B
    h_ids = [sc['id'] for sc in h_cfg['scenarios']]
    HOL = np.zeros((512, 3, 4))
    for j, sid in enumerate(h_ids):
        for i, s in enumerate(h_cfg['seeds']):
            HOL[i, j] = h_vals[sid][s]
    hol_banks = [np.arange(0, 256), np.arange(256, 512)]

    # ---------------- reproduction gate
    log('Reproduction gate (rel err <= 1e-12) ...')
    gate_rows = []

    def gate(name, arr, archived, rows_by_bank=None, src=''):
        tot = arr.sum(axis=0)
        B, P, F = rates(tot)
        errs = {}
        for k, v in (('control_rate', B), ('policy_rate', P), ('F', F)):
            a = archived[k]
            errs[k] = abs(v - a) / abs(a)
        row = dict(condition=name, source=src, n=int(arr.shape[0]),
                   control_rate=float(B), archived_control_rate=archived['control_rate'],
                   policy_rate=float(P), archived_policy_rate=archived['policy_rate'],
                   F=float(F), archived_F=archived['F'],
                   rel_err={k: float(v) for k, v in errs.items()})
        if 'control_count' in archived:
            row['counts_exact'] = (int(tot[0]) == archived['control_count'] and int(tot[2]) == archived['policy_count'])
            row['patient_days_rel_err'] = float(max(abs(tot[1] - archived['control_patient_days']) / archived['control_patient_days'],
                                                    abs(tot[3] - archived['policy_patient_days']) / archived['policy_patient_days']))
        if rows_by_bank is not None:
            row['banks'] = {}
            for bname, rows in rows_by_bank.items():
                bt = arr[rows].sum(axis=0)
                bB, bP, bF = rates(bt)
                ab = archived['banks'][bname]
                assert ab['n'] == len(rows)
                row['banks'][bname] = dict(n=len(rows), rel_err=dict(
                    control_rate=float(abs(bB - ab['control_rate']) / abs(ab['control_rate'])),
                    policy_rate=float(abs(bP - ab['policy_rate']) / abs(ab['policy_rate'])),
                    F=float(abs(bF - ab['F']) / abs(ab['F']))))
        all_errs = list(row['rel_err'].values()) + [e for b in row.get('banks', {}).values() for e in b['rel_err'].values()]
        if 'patient_days_rel_err' in row:
            all_errs.append(row['patient_days_rel_err'])
        row['max_rel_err'] = float(max(all_errs))
        row['passed'] = bool(row['max_rel_err'] <= REL_TOL and row.get('counts_exact', True))
        gate_rows.append(row)
        log(f'  {name:42s} B={B:.15g} P={P:.15g} F={F:.15g} max_rel_err={row["max_rel_err"]:.3e} '
            f'{"PASS" if row["passed"] else "FAIL"}')

    rs_by_id = {s['id']: s for s in rs_report['scenarios']}
    gate('independent_center_matched256', IND[:, 0], rs_by_id['independent_center_matched256'],
         {'A': ind_banks[0], 'B': ind_banks[1]}, 'rule-surface report')
    gate('independent_upper', IND[:, 1], rs_by_id['independent_upper'],
         {'A': ind_banks[0], 'B': ind_banks[1]}, 'rule-surface report')
    sens_by_id = {s['scenario']['id']: s for s in sens_report['scenarios']}
    for j, sid in enumerate(s_ids):
        gate(sid, SEN[:, j], sens_by_id[sid], {'A': sen_banks[0], 'B': sen_banks[1]}, 'sensitivity report')
    for j, sid in enumerate(h_ids):
        gate('old_full_' + sid, HOL[:, j], rs_by_id['old_full_' + sid],
             {'A': hol_banks[0], 'B': hol_banks[1]}, 'rule-surface report (S3 input)')
    # cross-check against the original holdout report too
    hold_by_id = {s['scenario']['id']: s for s in h_report['scenarios']}
    for j, sid in enumerate(h_ids):
        a = hold_by_id[sid]
        gate('old_full_' + sid + ' [holdout report]', HOL[:, j],
             dict(control_rate=a['arms']['baseline']['rate_per_1000_patient_days'],
                  policy_rate=a['arms']['policy']['rate_per_1000_patient_days'],
                  F=a['rate_prevented_fraction'], control_count=a['arms']['baseline']['acquired_total'],
                  policy_count=a['arms']['policy']['acquired_total'],
                  control_patient_days=a['arms']['baseline']['patient_days_total'],
                  policy_patient_days=a['arms']['policy']['patient_days_total']),
             None, 'original holdout report (S3 input)')
    gate_passed = all(r['passed'] for r in gate_rows) and not HASH_MISMATCHES
    log('Reproduction gate', 'PASSED' if gate_passed else 'FAILED')

    results = dict(
        schema='common_baseline_2026-10-16/impl_a/results', status=None,
        label='Post hoc secondary analysis; does not replace or soften any failed matching gate.',
        protocol_sha256=sha(HERE.parent / 'PROTOCOL.md'), script_sha256=sha(Path(__file__)),
        software=dict(python=sys.version, numpy=np.__version__),
        reproduction_gate=dict(passed=gate_passed, tolerance=REL_TOL, conditions=gate_rows,
                               hash_checks=HASH_CHECKS, single_file_hash_checks=SINGLE_HASH_CHECKS,
                               hash_mismatches=HASH_MISMATCHES),
    )
    if not gate_passed:
        results['status'] = 'stopped_reproduction_gate_failed'
        (HERE / 'results.json').write_text(json.dumps(results, indent=1))
        log('STOP: reproduction gate failed.')
        return 1

    # ---------------- sanity: model_kwargs across bracket pairs
    log('Sanity: model_kwargs across bracket pairs ...')
    sanity = []

    def kw_pair(name, lo_meta, hi_meta, seeds_lo, seeds_hi, same_seeds):
        diffs = {}
        n = 0
        if same_seeds:
            pairs = [(s, s) for s in seeds_lo]
        else:
            pairs = [(seeds_lo[0], seeds_hi[0])]
        # every seed within each condition: kwargs identical except seed
        within = set()
        for m, seeds in ((lo_meta, seeds_lo), (hi_meta, seeds_hi)):
            ref = kw_strip(m[seeds[0]]['kw'], {'seed'})
            for s in seeds:
                kws = kw_strip(m[s]['kw'], {'seed'})
                if kws != ref:
                    within.update(k for k in set(ref) | set(kws) if ref.get(k) != kws.get(k))
        for a, b in pairs:
            ka, kb = lo_meta[a]['kw'], hi_meta[b]['kw']
            n += 1
            for k in set(ka) | set(kb):
                if ka.get(k) != kb.get(k):
                    diffs.setdefault(k, set()).add((repr(ka.get(k)), repr(kb.get(k))))
        sa = lo_meta[seeds_lo[0]]['scenario']
        sb = hi_meta[seeds_hi[0]]['scenario']
        sc_diff = sorted(k for k in set(sa) | set(sb) if sa.get(k) != sb.get(k))
        other = sorted(k for k in diffs if k not in ('prob_transmission', 'seed'))
        entry = dict(pair=name, seed_pairs_compared=n, matched_seeds=same_seeds,
                     differing_kwargs=sorted(diffs), other_than_p_and_seed=other,
                     within_condition_nonseed_differences=sorted(within),
                     differing_scenario_fields=sc_diff,
                     prob_transmission=[sorted({x[0] for x in diffs.get('prob_transmission', set())}),
                                        sorted({x[1] for x in diffs.get('prob_transmission', set())})],
                     ok=(not other and not within and set(sc_diff) <= {'id', 'p'}))
        sanity.append(entry)
        log(f'  {name}: differing kwargs={entry["differing_kwargs"]} scenario diffs={sc_diff} '
            f'within-condition={entry["within_condition_nonseed_differences"]} ok={entry["ok"]}')

    ind_lo_meta = h_meta['independent_target_candidate']
    kw_pair('independent: center_matched256 vs upper', ind_lo_meta, up_meta, ind_seeds, ind_seeds, True)
    kw_pair('clustered: center vs high', s_meta['clustered_center'], s_meta['clustered_high'], s_cfg['seeds'], s_cfg['seeds'], True)
    kw_pair('imperfect: center vs high', s_meta['imperfect_center'], s_meta['imperfect_high'], s_cfg['seeds'], s_cfg['seeds'], True)
    kw_pair('clustered S2: low vs center', s_meta['clustered_low'], s_meta['clustered_center'], s_cfg['seeds'], s_cfg['seeds'], True)
    kw_pair('imperfect S2: low vs center', s_meta['imperfect_low'], s_meta['imperfect_center'], s_cfg['seeds'], s_cfg['seeds'], True)
    kw_pair('S3 independent: holdout512 center vs upper', ind_lo_meta, up_meta, h_cfg['seeds'], ind_seeds, False)
    kw_pair('S3 clustered: holdout512 persistent center vs sensitivity high', h_meta['persistent_target_candidate'],
            s_meta['clustered_high'], h_cfg['seeds'], s_cfg['seeds'], False)
    kw_pair('S3 imperfect: holdout512 center vs sensitivity high', h_meta['imperfect_target_candidate'],
            s_meta['imperfect_high'], h_cfg['seeds'], s_cfg['seeds'], False)
    # identical p across modules where the protocol treats them as the same center
    p_checks = dict(
        clustered_center_equals_holdout_persistent_p=s_cfg['scenarios'][s_ids.index('clustered_center')]['p'] == h_cfg['scenarios'][1]['p'],
        imperfect_center_equals_holdout_imperfect_p=s_cfg['scenarios'][s_ids.index('imperfect_center')]['p'] == h_cfg['scenarios'][2]['p'],
        independent_upper_p_over_center=up_cfg['scenarios'][0]['p'] / h_cfg['scenarios'][0]['p'])

    # ---------------- profile definitions
    p_of = {sc['id']: sc['p'] for sc in s_cfg['scenarios']}
    p_ind_lo = h_cfg['scenarios'][0]['p']
    p_ind_hi = up_cfg['scenarios'][0]['p']
    profiles = ['independent', 'clustered', 'imperfect']
    lower_label = dict(independent='independent_center_matched256', clustered='clustered_center', imperfect='imperfect_center')
    upper_label = dict(independent='independent_upper', clustered='clustered_high', imperfect='imperfect_high')
    p_lo = dict(independent=p_ind_lo, clustered=p_of['clustered_center'], imperfect=p_of['imperfect_center'])
    p_hi = dict(independent=p_ind_hi, clustered=p_of['clustered_high'], imperfect=p_of['imperfect_high'])
    data_lo = dict(independent=IND[:, 0], clustered=SEN[:, s_ids.index('clustered_center')],
                   imperfect=SEN[:, s_ids.index('imperfect_center')])
    data_hi = dict(independent=IND[:, 1], clustered=SEN[:, s_ids.index('clustered_high')],
                   imperfect=SEN[:, s_ids.index('imperfect_high')])

    # ---------------- bootstrap weights (fixed, documented draw order)
    log('Drawing bootstrap weights ...')
    rng = np.random.default_rng(BOOT_SEED)

    def bank_weights(n_bank):
        return rng.multinomial(n_bank, np.full(n_bank, 1.0 / n_bank), size=DRAWS)

    W_ind = np.zeros((DRAWS, 256), dtype=np.int64)
    W_sen = np.zeros((DRAWS, 256), dtype=np.int64)
    W_ind[:, ind_banks[0]] = bank_weights(128)     # 1. independent module bank A
    W_ind[:, ind_banks[1]] = bank_weights(128)     # 2. independent module bank B
    W_sen[:, sen_banks[0]] = bank_weights(128)     # 3. sensitivity module bank A
    W_sen[:, sen_banks[1]] = bank_weights(128)     # 4. sensitivity module bank B
    W_hol = np.zeros((DRAWS, 512), dtype=np.int64)
    W_hol[:, hol_banks[0]] = bank_weights(256)     # 5. S3 holdout bank A (256)
    W_hol[:, hol_banks[1]] = bank_weights(256)     # 6. S3 holdout bank B (256)
    W_hol_pooled = bank_weights(512)               # 7. S3 variant: unstratified 512 (archived holdout primary convention)
    for W, n in ((W_ind, 256), (W_sen, 256), (W_hol, 512), (W_hol_pooled, 512)):
        assert (W.sum(axis=1) == n).all()
    weight_hashes = {k: hashlib.sha256(np.ascontiguousarray(v).astype('<i8').tobytes()).hexdigest()
                     for k, v in (('independent', W_ind), ('sensitivity', W_sen), ('holdout_AB256', W_hol),
                                  ('holdout_pooled512', W_hol_pooled))}
    assert (W_ind[:, :128].sum(axis=1) == 128).all() and (W_sen[:, 128:].sum(axis=1) == 128).all()

    weights_for = dict(independent=W_ind, clustered=W_sen, imperfect=W_sen)

    def boot_sums(W, arr):
        return (W.astype(np.float64) @ arr.reshape(arr.shape[0], -1)).reshape(DRAWS, *arr.shape[1:])

    # bootstrap sums per condition
    BS = {}
    BS['ind'] = boot_sums(W_ind, IND)          # (DRAWS, 2, 4)
    BS['sen'] = boot_sums(W_sen, SEN)          # (DRAWS, 6, 4)
    BS['hol'] = boot_sums(W_hol, HOL)          # (DRAWS, 3, 4)
    BS['hol_pooled'] = boot_sums(W_hol_pooled, HOL)
    del W_hol, W_hol_pooled
    boot_lo = dict(independent=BS['ind'][:, 0], clustered=BS['sen'][:, s_ids.index('clustered_center')],
                   imperfect=BS['sen'][:, s_ids.index('imperfect_center')])
    boot_hi = dict(independent=BS['ind'][:, 1], clustered=BS['sen'][:, s_ids.index('clustered_high')],
                   imperfect=BS['sen'][:, s_ids.index('imperfect_high')])

    # ---------------- primary estimator (generic for any lower/upper pairing and target)
    def estimate(lo_tot, hi_tot, blo, bhi, plo, phi, target):
        B_lo, P_lo, F_lo = rates(lo_tot)
        B_hi, P_hi, F_hi = rates(hi_tot)
        w, Fs, _, _ = interp(np.float64(B_lo), np.float64(B_hi), np.float64(P_lo), np.float64(P_hi), target)
        # cross-check: identical to direct log-rate interpolation
        lnP = (1 - w) * math.log(P_lo) + w * math.log(P_hi)
        assert abs((1 - math.exp(lnP) / target) - Fs) < 1e-12
        bB_lo, bP_lo, _ = rates(blo)
        bB_hi, bP_hi, _ = rates(bhi)
        bw, bFs, _, _ = interp(bB_lo, bB_hi, bP_lo, bP_hi, target)
        und = undefined_mask(blo, bhi, bB_lo, bB_hi)
        extrap = (bB_lo >= target) | (bB_hi <= target)
        nonincreasing = bB_hi <= bB_lo
        return dict(B_lo=float(B_lo), B_hi=float(B_hi), P_lo=float(P_lo), P_hi=float(P_hi),
                    F_lo=float(F_lo), F_hi=float(F_hi), w=float(w), F_star=float(Fs),
                    p_lo=plo, p_hi=phi, p_star=float(p_star(w, plo, phi)),
                    point_extrapolated=bool(B_lo >= target or B_hi <= target),
                    undefined_draws=int(und.sum()), frac_extrapolated=float(extrap[~und].mean()),
                    n_extrapolated=int(extrap[~und].sum()), nonincreasing_draws=int(nonincreasing.sum())), \
            dict(F_star=bFs, w=bw, und=und, p_star=p_star(bw, plo, phi))

    def run_set(target, lo_points, hi_points, lo_boot, hi_boot, label):
        out = {}
        draws = {}
        for h in profiles:
            pt, bd = estimate(lo_points[h], hi_points[h], lo_boot[h], hi_boot[h], p_lo[h], p_hi[h], target)
            lo95, hi95 = pct_interval(bd['F_star'], LEVEL_F, bd['und'])
            pt['ci95_lo'], pt['ci95_hi'] = lo95, hi95
            pt['p_star_ci95'] = list(pct_interval(bd['p_star'], LEVEL_F, bd['und']))
            pt['w_ci95'] = list(pct_interval(bd['w'], LEVEL_F, bd['und']))
            pt['F_star_boot_mean'] = float(np.mean(bd['F_star'][~bd['und']]))
            pt['profile'] = h
            pt['lower_condition'] = lower_label[h] if label != 'S3' else s3_lower_label[h]
            pt['upper_condition'] = upper_label[h]
            out[h] = pt
            draws[h] = bd
            log(f'  [{label} B*={target}] {h:12s} B_lo={pt["B_lo"]:.6f} B_hi={pt["B_hi"]:.6f} F_lo={pt["F_lo"]:.6f} '
                f'F_hi={pt["F_hi"]:.6f} w={pt["w"]:.6f} p*={pt["p_star"]:.8f} F*={pt["F_star"]:.6f} '
                f'95%=[{lo95:.6f},{hi95:.6f}] extrap={pt["frac_extrapolated"]:.5f} undef={pt["undefined_draws"]}')
        return out, draws

    def contrasts_of(points, draws, label, level=LEVEL_C):
        out = []
        for a, b in (('independent', 'clustered'), ('independent', 'imperfect'), ('clustered', 'imperfect')):
            d = draws[a]['F_star'] - draws[b]['F_star']
            und = draws[a]['und'] | draws[b]['und']
            lo, hi = pct_interval(d, level, und)
            est = points[a]['F_star'] - points[b]['F_star']
            out.append(dict(pair=f'{a}-{b}', estimate=float(est), ci_lo=lo, ci_hi=hi, level=level,
                            resolved=bool(lo > 0 or hi < 0), undefined_draws=int(und.sum())))
            log(f'  [{label}] {a}-{b}: {est:.6f} [{lo:.6f},{hi:.6f}] level={level:.6f} resolved={lo > 0 or hi < 0}')
        return out

    lo_pts = {h: data_lo[h].sum(axis=0) for h in profiles}
    hi_pts = {h: data_hi[h].sum(axis=0) for h in profiles}
    s3_lower_label = dict(independent='old_full_independent_target_candidate (512)',
                          clustered='old_full_persistent_target_candidate (512)',
                          imperfect='old_full_imperfect_target_candidate (512)')

    log('Primary (B*=3) ...')
    primary, prim_draws = run_set(TARGET, lo_pts, hi_pts, boot_lo, boot_hi, 'primary')
    contrasts = contrasts_of(primary, prim_draws, 'primary')

    # ---------------- S1: arithmetic interpolation of F in ln B (same w)
    log('S1: F linear in ln B ...')
    S1 = {}
    s1_draws = {}
    for h in profiles:
        pt = primary[h]
        F1 = (1 - pt['w']) * pt['F_lo'] + pt['w'] * pt['F_hi']
        bB_lo, bP_lo, bF_lo = rates(boot_lo[h])
        bB_hi, bP_hi, bF_hi = rates(boot_hi[h])
        bw = prim_draws[h]['w']
        bF1 = (1 - bw) * bF_lo + bw * bF_hi
        und = prim_draws[h]['und']
        lo95, hi95 = pct_interval(bF1, LEVEL_F, und)
        S1[h] = dict(F_star=float(F1), ci95_lo=lo95, ci95_hi=hi95, w=pt['w'],
                     minus_primary=float(F1 - pt['F_star']), undefined_draws=int(und.sum()))
        s1_draws[h] = dict(F_star=bF1, und=und)
        log(f'  S1 {h:12s} F*={F1:.6f} 95%=[{lo95:.6f},{hi95:.6f}] (primary {pt["F_star"]:.6f})')
    S1_contrasts = contrasts_of(S1, s1_draws, 'S1')

    # ---------------- S2: quadratic in ln p through low/center/high (clustered, imperfect)
    log('S2: quadratic in ln p (clustered, imperfect) ...')

    def quad_coefs(x, y):
        """Newton form through (x0,y0),(x1,y1),(x2,y2); returns (c0,c1,c2) for c0 + c1 x + c2 x^2."""
        x0, x1, x2 = x
        y0, y1, y2 = y
        d1 = (y1 - y0) / (x1 - x0)
        d2 = ((y2 - y1) / (x2 - x1) - d1) / (x2 - x0)
        c2 = d2
        c1 = d1 - d2 * (x0 + x1)
        c0 = y0 - d1 * x0 + d2 * x0 * x1
        return c0, c1, c2

    def s2_solve(x, yB, yP, x_ref, target=TARGET):
        c0, c1, c2 = quad_coefs(x, yB)
        a, b, c = c2, c1, c0 - math.log(target)
        with np.errstate(invalid='ignore', divide='ignore'):
            disc = b * b - 4 * a * c
            sq = np.sqrt(np.where(disc >= 0, disc, np.nan))
            r1 = (-b + sq) / (2 * a)
            r2 = (-b - sq) / (2 * a)
            lin = -c / b
        tiny = np.abs(a) < 1e-14
        r1 = np.where(tiny, lin, r1)
        r2 = np.where(tiny, lin, r2)
        pick = np.where(np.abs(r1 - x_ref) <= np.abs(r2 - x_ref), r1, r2)
        pick = np.where(np.isnan(r1), r2, np.where(np.isnan(r2), r1, pick))
        d0, d1_, d2_ = quad_coefs(x, yP)
        lnP = d0 + d1_ * pick + d2_ * pick * pick
        F = 1 - np.exp(lnP) / target
        no_root = ~np.isfinite(pick)
        outside = (pick < x[0]) | (pick > x[2])
        # local slope of ln B at the root (positive = increasing branch)
        slope = c1 + 2 * c2 * pick
        return F, pick, no_root, outside, slope, (c0, c1, c2)

    S2 = {}
    S2_draws = {}
    for h in ('clustered', 'imperfect'):
        ids3 = [f'{h}_low', f'{h}_center', f'{h}_high']
        js = [s_ids.index(i) for i in ids3]
        x = np.array([math.log(p_of[i]) for i in ids3])
        tot = SEN[:, js].sum(axis=0)
        B3, P3, F3 = rates(tot)
        x_ref = math.log(primary[h]['p_star'])
        F2, xr, nr, outside, slope, coefs = s2_solve(x, np.log(B3), np.log(P3), np.float64(x_ref))
        # bootstrap
        bt = BS['sen'][:, js]
        bB3, bP3, _ = rates(bt)
        bx_ref = np.log(prim_draws[h]['p_star'])
        with np.errstate(divide='ignore', invalid='ignore'):
            bF2, bxr, bnr, bout, bslope, _ = s2_solve(x, np.log(bB3).T, np.log(bP3).T, bx_ref)
        zero = (bt[:, :, 0] <= 0).any(axis=1) | (bt[:, :, 2] <= 0).any(axis=1)
        und = bnr | zero | ~np.isfinite(bF2)
        lo95, hi95 = pct_interval(bF2, LEVEL_F, und)
        S2_draws[h] = dict(F_star=bF2, und=und)
        S2[h] = dict(conditions=ids3, p=[p_of[i] for i in ids3], B=[float(v) for v in B3], P=[float(v) for v in P3],
                     F=[float(v) for v in F3], lnB_quadratic_coefs_in_lnp=[float(c) for c in coefs],
                     curvature_lnB=float(coefs[2]),
                     F_star=float(F2), p_star=float(math.exp(xr)), ci95_lo=lo95, ci95_hi=hi95,
                     minus_primary=float(F2 - primary[h]['F_star']),
                     root_outside_low_high=bool(outside), slope_at_root=float(slope),
                     undefined_draws=int(und.sum()), no_real_root_draws=int(bnr.sum()),
                     root_outside_low_high_draws=int((bout & ~und).sum()),
                     decreasing_branch_root_draws=int(((bslope <= 0) & ~und).sum()),
                     root_rule='root of quadratic ln B(ln p) = ln 3 nearest the primary (center-high log-linear) ln p* of the same draw')
        log(f'  S2 {h:12s} F*={F2:.6f} p*={math.exp(xr):.8f} 95%=[{lo95:.6f},{hi95:.6f}] (primary {primary[h]["F_star"]:.6f}) '
            f'undef={int(und.sum())} outside={int((bout & ~und).sum())}')
    # S2 contrast clustered - imperfect (only pair available); reuses the per-profile S2 draws
    d = S2_draws['clustered']['F_star'] - S2_draws['imperfect']['F_star']
    und = S2_draws['clustered']['und'] | S2_draws['imperfect']['und']
    lo, hi = pct_interval(d, LEVEL_C, und)
    S2_contrast = dict(pair='clustered-imperfect', estimate=float(S2['clustered']['F_star'] - S2['imperfect']['F_star']),
                       ci_lo=lo, ci_hi=hi, level=LEVEL_C, resolved=bool(lo > 0 or hi < 0), undefined_draws=int(und.sum()))
    log(f'  S2 clustered-imperfect: {S2_contrast["estimate"]:.6f} [{lo:.6f},{hi:.6f}]')

    # ---------------- S3: 512-seed holdout center as lower bracket
    log('S3: holdout 512-seed centers as lower brackets (A/B 256 stratified, independent of upper) ...')
    hol_lo_pts = {h: HOL[:, j].sum(axis=0) for h, j in zip(profiles, range(3))}
    hol_lo_boot = {h: BS['hol'][:, j] for h, j in zip(profiles, range(3))}
    p_lo_saved = dict(p_lo)
    # holdout center p values equal primary lower p for clustered/imperfect; check and use holdout p
    for h, j in zip(profiles, range(3)):
        p_lo[h] = h_cfg['scenarios'][j]['p']
        assert p_lo[h] == p_lo_saved[h]
    S3, s3_draws = run_set(TARGET, hol_lo_pts, hi_pts, hol_lo_boot, boot_hi, 'S3')
    S3_contrasts = contrasts_of(S3, s3_draws, 'S3')
    log('S3 variant: holdout resampled as one unstratified 512-seed multinomial ...')
    hol_lo_boot_pooled = {h: BS['hol_pooled'][:, j] for h, j in zip(profiles, range(3))}
    S3p, s3p_draws = run_set(TARGET, hol_lo_pts, hi_pts, hol_lo_boot_pooled, boot_hi, 'S3')
    S3p_contrasts = contrasts_of(S3p, s3p_draws, 'S3-pooled512')

    # ---------------- S4: primary estimator at B*=2.7 and 3.3
    S4 = {}
    S4_contrasts = {}
    for tgt in (2.7, 3.3):
        log(f'S4: B*={tgt} ...')
        pts, dr = run_set(tgt, lo_pts, hi_pts, boot_lo, boot_hi, 'S4')
        S4[str(tgt)] = pts
        S4_contrasts[str(tgt)] = contrasts_of(pts, dr, f'S4 B*={tgt}')

    results.update(
        status='complete',
        bank_definitions=dict(
            independent_module=dict(
                seeds='rule-surface configs/independent_upper.json seeds = holdout seed_blocks.holdout A[:128] + B[:128] '
                      '(= old_input_index.json independent_upper_comparison_seeds = rule-surface seedbook independent_upper A+B)',
                bank_A='positions 0-127 (holdout bank A core: first 128 of ordered bank A)',
                bank_B='positions 128-255 (holdout bank B core: first 128 of ordered bank B)',
                lower='holdout records with scenario.id independent_target_candidate for those 256 seeds, located via '
                      'old_input_index.json record_metadata (seed, scenario), both holdout sites',
                upper='256 independent_upper records'),
            sensitivity_module=dict(seeds='author_notes configs/sensitivity.json seeds = seedbook.json sensitivity A + B',
                                    bank_A='positions 0-127', bank_B='positions 128-255'),
            holdout_S3=dict(seeds='holdout_all.json seeds = full_icu_confirmation seedbook seed_blocks.holdout A + B',
                            bank_A='positions 0-255', bank_B='positions 256-511',
                            convention='S3 main: independent multinomial resampling within holdout banks A and B (256 each), '
                                       'synchronized across the three holdout profiles. Variant: one unstratified 512-seed '
                                       'multinomial (the archived holdout report primary bootstrap convention).'),
            draw_order='default_rng(20261016); multinomial(n_bank, uniform, size=50000) in order: independent A(128), '
                       'independent B(128), sensitivity A(128), sensitivity B(128), holdout A(256), holdout B(256), '
                       'holdout pooled(512)',
            weight_sha256_int64LE=weight_hashes),
        sanity_model_kwargs=dict(pairs=sanity, p_checks=p_checks),
        primary=dict(target=TARGET, level=LEVEL_F, profiles=primary),
        contrasts=dict(level=LEVEL_C, family_size=3, pairs=contrasts),
        S1=dict(description='F linear in ln B between the same brackets (same w)', profiles=S1, contrasts=S1_contrasts),
        S2=dict(description='quadratic in ln p through low/center/high for ln B and ln P, evaluated at B=3',
                profiles=S2, contrast=S2_contrast),
        S3=dict(description='lower bracket = original 512-seed holdout center; lower resampled independently of upper',
                main_AB256=dict(profiles=S3, contrasts=S3_contrasts),
                variant_pooled512=dict(profiles=S3p, contrasts=S3p_contrasts)),
        S4=dict(description='primary estimator at B*=2.7 and B*=3.3', targets=S4, contrasts=S4_contrasts),
        bootstrap=dict(draws=DRAWS, seed=BOOT_SEED, interval_method='numpy.quantile (linear) percentile over defined draws'),
        runtime_seconds=time.time() - t0,
    )
    (HERE / 'results.json').write_text(json.dumps(results, indent=1))
    log('wrote', HERE / 'results.json')
    log(f'done in {time.time() - t0:.1f}s')
    return 0


if __name__ == '__main__':
    try:
        rc = main()
    except AssertionError as e:
        import traceback
        traceback.print_exc()
        rc = 2
    sys.exit(rc)
