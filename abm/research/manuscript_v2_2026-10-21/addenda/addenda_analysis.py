#!/usr/bin/env python
"""Analysis-only addenda for manuscript v2 (21 October 2026). No simulation.

Reads archived, hash-checked records only and writes addenda.json next to this
script. Four parts:

  1. Baseline-mismatch elasticity: slope of prevented fraction F against
     ln(baseline rate) per profile, from seed-paired original-workflow
     conditions, with whole-seed bootstrap intervals; converted to the implied
     bias in F from each observed baseline shortfall relative to 3.0, for the
     original and alternative workflows.
  2. Selection-uncertainty propagation: archived selected-p intervals x dF/dln p.
  3. Route-label table for both workflows (descriptive labels).
  4. Data for the study-flow table and the compliance-range table (Table S5).

Bootstrap: 20,000 whole-seed draws. numpy SeedSequence(20261021) is spawned
into five independent child streams, one per seed list:
  [0] author-notes sensitivity bank (A/B 128, shared by all six conditions)
  [1] original holdout bank (four 128-seed strata: A-core, A-rem, B-core, B-rem);
      the independent upper condition uses the two cores (its seeds are
      A[:128] + B[:128]), so overlap with the holdout is preserved exactly
  [2] environmental-comparison bank (A/B 256, shared by all six conditions)
  [3] alternative-workflow policy bank (512, unstratified, as in its protocol)
  [4] alternative-workflow grid bank (128, unstratified)
Rates are ratios of pooled counts to pooled patient-days in every draw.
Percentile intervals; undefined (non-finite) draws are counted and reported.
"""
import hashlib
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ABM = Path(__file__).resolve().parents[3]
DM = ABM.parent

DRAWS = 20000
BOOT_SEED = 20261021
TARGET = 3.0
CHUNK = 2000

# ----------------------------------------------------------------- input paths
SENS_ROOT = DM / 'research/author_notes_2026-10-10/collections/sensitivity_20261011_101348_utc/results/sensitivity'
SENS_REPORT = DM / 'research/author_notes_2026-10-10/analysis_2026-10-11/sensitivity_replay_20261011_1015_utc/report.json'
AN_SEEDBOOK = ABM / 'research/author_notes_2026-10-10/seedbook.json'
UPPER_ROOT = (DM / 'research/rule_surface_followup_2026-10-11/collections/independent_upper_20261012_0125_utc'
              / 'restored/repo/results/rule_surface_followup_2026-10-11/independent_upper')
RS_REPORT = ABM / 'research/rule_surface_results_2026-10-12/completed_analysis/report.json'
RS_INDEX = ABM / 'research/rule_surface_followup_2026-10-11/old_input_index.json'
CONF_RESULTS = ABM / 'results/full_icu_confirmation_2026-10-08'
CONF_SEEDBOOK = ABM / 'research/full_icu_confirmation_2026-10-08/seedbook.json'
HOLD_REPORT = ABM / 'research/full_icu_confirmation_2026-10-08/holdout_report/report.json'
ENV_REPORT = ABM / 'research/full_icu_confirmation_2026-10-08/environment_report/report.json'
ALT_EXEC = DM / 'research/robustness_uq_2026-10-17/collection_2026-10-21/execution01'
ALT_LOCAL_SHA = DM / 'research/robustness_uq_2026-10-17/collection_2026-10-21/local_sha256.txt'
ALT_SEEDBOOK = ALT_EXEC / 'source/research/robustness_uq_2026-10-17/seedbook.json'
ALT_SELECTION = ALT_EXEC / 'source/research/robustness_uq_2026-10-17/selection.json'
ALT_ANALYSIS = ALT_EXEC / 'outputs/analysis.json'
ORIG_SEL_INDEP = ABM / 'research/author_notes_2026-10-10/saved_analysis/outputs/calibration.json'
ORIG_SEL_EXPL = ABM / 'research/editorial_polish_2026-10-11/source_revision_02/input/evidence/selection_uncertainty.json'
COMMON_BASELINE = ABM / 'research/common_baseline_2026-10-16/results_final.json'
SUPPLEMENT = ABM / 'manuscript/supplement.tex'

PENDING = '[PENDING: workflow_coefficient_2026-10-21 results not yet available]'

# Study-flow rows, transcribed from the cited archived sources (counts are records or paired jobs).
STUDY_FLOW = [
    dict(step=1, study='Antecedent exploration (before 8 Oct)',
         question='Select p for each profile at a baseline of 3.0 per 1,000 patient-days',
         p_selection='Interpolation within the lowest-p increasing adjacent bracket (archived method '
                     '"lowest_p_increasing_adjacent_bracket"); no fallback or extrapolation',
         n='27 conditions x 64 paired seeds (1,728 records)',
         outcome='Clustered p = .07602 [.07346, .07873], imperfect p = .07484 [.07243, .07785] (retrospective '
                 'conditional 95%); independent: no point selection',
         source='manuscript Methods (Calibration and fresh confirmation); supplement Table S14; '
                'editorial_polish_2026-10-11/.../selection_uncertainty.json (method field)'),
    dict(step=2, study='Independent-profile calibration grid (8 Oct)',
         question='Select p for the independent profile',
         p_selection='Declared bracket/interpolation rule on a baseline-only grid, p = .15 fallback',
         n='3 probabilities (.13497, .15, .165) x 128 seeds (384 baseline-only records)',
         outcome='p = .14328 [.14004, .14792] (retrospective conditional 95%; no failed pooled draws)',
         source='manuscript Methods; supplement Table S8'),
    dict(step=3, study='Original-workflow holdout confirmation (8 Oct)',
         question='Prevented fraction F for completion 90% to 95% at the fixed selected p',
         p_selection='Fixed from steps 1-2; policy outcomes never select p',
         n='3 profiles x 512 fresh paired seeds (1,536 paired jobs)',
         outcome='Control 2.83/2.69/2.80; F 67.6/55.9/27.0%; contrasts 11.7, 40.6, 28.9 pp, all exclude 0 '
                 '(98.333%); all six baseline-matching checks failed',
         source='manuscript Results (Direct simulation after baseline fitting)'),
    dict(step=4, study='Environmental comparison (8 Oct)',
         question='Does disabling hand-site exchange change F? (independent profile)',
         p_selection='Three fixed earlier probabilities (.0810, .1080, .1350); no fitting',
         n='6 conditions x 512 paired seeds (3,072 paired jobs)',
         outcome='Full minus disabled F: -10.2 [-22.4, 2.4], -11.7 [-20.7, -2.7], -14.3 [-21.1, -7.5] pp '
                 '(98.333%)',
         source='manuscript Results (Environmental routes at three probabilities)'),
    dict(step=5, study='Rare-transmission coefficient, original workflow (8-10 Oct)',
         question='Limit of F as p -> 0 from transmission-free histories',
         p_selection='None (p = 0 backgrounds)',
         n='256 transmission-free backgrounds; 22-row extended surface',
         outcome='F* = 75.1/53.6/28.0%; extended surface in Table S5',
         source='manuscript Appendix (Full ICU model coefficient); supplement Table S5'),
    dict(step=6, study='Inherited-clock baseline sweep (10-11 Oct)',
         question='Can the inherited per-contact washing clock reach the target?',
         p_selection='Baseline only at p = .15, .30, .60, 1',
         n='12 conditions x 128 seeds (1,536 baseline-only records)',
         outcome='No sampled match to the 2.7-3.3 band',
         source='supplement Table S1'),
    dict(step=7, study='Local p-sensitivity (10-11 Oct)',
         question='Is F locally sensitive to p? (clustered, imperfect)',
         p_selection='0.75, 1 and 1.25 x the fixed p',
         n='6 conditions x 256 fresh paired seeds (1,536 paired jobs)',
         outcome='F 60.5/55.3/53.6% (clustered) and 28.9/29.8/32.7% (imperfect); four contrasts unresolved',
         source='supplement Tables S2-S4'),
    dict(step=8, study='Opportunity-rule and matched-surface follow-up (11-12 Oct)',
         question='Opportunity rule, surface disabling and an independent upper probability',
         p_selection='Lowest exact target hit or lowest increasing adjacent bracket, interpolation within that bracket; '
                     'no fallback, extrapolation or reselection',
         n='5,760 records: 19-condition grid x 128 (2,432), 2 x 512 selected-clock pairs, 512 p = 1 pairs, '
           '3 x 512 surface-disabled pairs, 256 independent-upper pairs',
         outcome='Selected p .283/.345; both target checks failed (3.09, 2.85); six-contrast family: +16.3, +37.0, '
                 '+9.7, -4.5, +0.6, +4.7 pp (99.167%)',
         source='manuscript Methods; supplement Tables S15-S22'),
    dict(step=9, study='Common-baseline secondary analysis (16-17 Oct, post hoc)',
         question='F at a common baseline of 3.0',
         p_selection='Log-log interpolation between center and 1.25 x center; no new runs',
         n='No new runs; 256 seed-matched pairs per bracket',
         outcome='F* 66.7 [62.1, 71.2], 55.0 [50.7, 59.4], 30.0 [23.6, 36.6]%; contrasts 11.7, 36.7, 25.0 pp '
                 'exclude 0 (98.333%)',
         source='common_baseline_2026-10-16/RESULTS.md; manuscript Table (common baseline)'),
    dict(step=10, study='Alternative-workflow robustness (17-21 Oct)',
         question='Do the profile contrasts persist with vacancy cleaning plus empty-bed-only isolation?',
         p_selection='Lowest positive increasing adjacent bracket, log-log interpolation on a 0.5/1/2/4 x p grid',
         n='Grid 4 x 3 x 128 (1,536 baseline-only records); holdout 3 x 512 fresh paired seeds (1,536 paired jobs)',
         outcome='p .16511/.08063/.07614; control 2.90/2.62/2.63; F 76.5/53.7/27.3%; contrasts 22.8, 49.2, 26.4 pp '
                 'exclude 0; clustered and imperfect baselines below 2.7',
         source='robustness_uq_2026-10-17/verification_2026-10-21/RESULTS.md'),
    dict(step=11, study='Synthetic measurement UQ (17 Oct)',
         question='Does a two-opportunity survival assay narrow the predicted benefit?',
         p_selection='Not applicable (known synthetic truths)',
         n='12,000 synthetic assay experiments (4 truths x 3 sample sizes x 1,000)',
         outcome='At 8,000 per assay, continuous k = 2 widths 25.0 -> 6.3, 25.0 -> 1.1, 11.5 -> 5.7 and '
                 '11.5 -> 6.7 pp; not ICU data',
         source='robustness_uq_2026-10-17/uq/RESULTS.md; uq/continuous_k2_2026-10-17/RESULTS.md'),
    dict(step=12, study='Mathematical extension (17 Oct)',
         question='Identification of the policy response from bridge survival features',
         p_selection='Not applicable (analytic)',
         n='Analytic; 229 exact rational assertion groups',
         outcome='Identification criterion and sharp bounds; two-opportunity clustering effect changes sign at '
                 'e ~ .9235',
         source='mathematical_extension_2026-10-17/START_HERE.txt'),
    dict(step=13, study='Rare-transmission coefficient, alternative workflow (21 Oct)',
         question='p -> 0 limit under the alternative workflow, paired with the original on the same backgrounds',
         p_selection='None (p = 0 backgrounds)',
         n='256 backgrounds for the three profiles; 22-row surface on 256 or on bank A (128) per the budget rule',
         outcome=PENDING,
         source='workflow_coefficient_2026-10-21/PROTOCOL.md'),
    dict(step=14, study='These addenda (21-22 Oct)',
         question='Baseline-mismatch elasticity, selection propagation, route labels',
         p_selection='Not applicable',
         n='No new runs; {n_records} archived records re-read',
         outcome='See ADDENDA.md sections 1-3',
         source='manuscript_v2_2026-10-21/addenda'),
]

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


HASH_CHECKS = []


def check_hashes(label, pairs):
    """pairs: iterable of (path, expected sha256)."""
    bad = [str(p) for p, e in pairs if sha(p) != e]
    HASH_CHECKS.append(dict(label=label, n_files=len(pairs), mismatches=len(bad), first_mismatches=bad[:5]))
    log(f'  hash check {label}: {len(pairs)} files, {len(bad)} mismatches')
    if bad:
        raise SystemExit(f'hash mismatch in {label}')


def input_digest(paths):
    h = hashlib.sha256()
    for p in sorted(paths, key=lambda q: q.name):
        h.update(p.name.encode())
        h.update(sha(p).encode())
    return h.hexdigest()


# --------------------------------------------------------------- record reading
ROUTE_KEYS = ('patient', 'own0', 'ownpos', 'other', 'carry0', 'carrypos')


def route_bins(rc):
    out = dict.fromkeys(ROUTE_KEYS, 0)
    for k, v in rc.items():
        route, gap = k.rsplit('|', 1)
        g = int(gap)
        if route == 'patient_pickup':
            out['patient'] += v
        elif route == 'own_zone_environment_pickup':
            out['own0' if g == 0 else 'ownpos'] += v
        elif route == 'other_zone_environment_pickup':
            out['other'] += v
        elif route == 'own_zone_environment_pickup_after_readmission':
            out['carry0' if g == 0 else 'carrypos'] += v
        else:
            raise ValueError('unknown route ' + route)
    return out


def read_record(path, paired=True):
    row = load(path)
    assert row['status'] == 'complete', path
    out = dict(seed=row['seed'], sid=row['scenario']['id'], p=row['scenario']['p'], scenario=row['scenario'],
               kw=row['model_kwargs'])
    assert row['model_kwargs']['seed'] == row['seed']
    assert row['model_kwargs']['prob_transmission'] == row['scenario']['p']
    arms = (('baseline', 0.9), ('policy', 0.95)) if paired else (('baseline', 0.9),)
    vals = []
    for arm, h in arms:
        a = row['arms'][arm]
        assert a['h'] == h, (path, arm, a['h'])
        obs = a['observation']
        assert obs['days'] == 56 and obs['ticks'] == 145152
        c, d = obs['acquired_colonizations'], obs['patient_days']
        assert isinstance(c, int) and c >= 0 and d > 0
        rc = obs['response_counts']['acquisitions_by_route_and_opportunities']
        rb = route_bins(rc)
        assert sum(rb.values()) == c, (path, arm, rb, c)
        assert sum(obs['response_counts']['acquisitions_by_route'].values()) == c
        vals += [c, d]
        out['routes_' + arm] = rb
    out['vals'] = vals
    return out


def load_bank(files, seed_order, paired=True):
    """Return {sid: dict(p, X (n x 4 or n x 2), routes0, routes1 (n x 6), meta)} in seed_order."""
    by = defaultdict(dict)
    meta = {}
    for f in files:
        r = read_record(f, paired)
        assert r['seed'] not in by[r['sid']], ('duplicate', r['sid'], r['seed'])
        by[r['sid']][r['seed']] = r
        meta.setdefault(r['sid'], r['scenario'])
    out = {}
    for sid, recs in by.items():
        assert set(recs) == set(seed_order), (sid, len(recs), len(seed_order))
        rows = [recs[s] for s in seed_order]
        ps = {r['p'] for r in rows}
        assert len(ps) == 1
        X = np.array([r['vals'] for r in rows], dtype=float)
        d = dict(p=ps.pop(), X=X, scenario=meta[sid], n=len(rows))
        d['routes0'] = np.array([[r['routes_baseline'][k] for k in ROUTE_KEYS] for r in rows], dtype=float)
        if paired:
            d['routes1'] = np.array([[r['routes_policy'][k] for k in ROUTE_KEYS] for r in rows], dtype=float)
        # model settings identical across seeds apart from the seed itself
        kws = {json.dumps({k: v for k, v in r['kw'].items() if k != 'seed'}, sort_keys=True) for r in rows}
        assert len(kws) == 1, sid
        d['model_kwargs_except_seed'] = json.loads(kws.pop())
        out[sid] = d
    return out


# ---------------------------------------------------------------- bootstrap
def weights(rng, strata, draws=DRAWS):
    """Whole-seed resampling within strata. strata: list of index arrays.
    Returns int32 count matrix (draws x n_total)."""
    n = sum(len(s) for s in strata)
    W = np.zeros((draws, n), dtype=np.int32)
    for s in strata:
        s = np.asarray(s)
        idx = rng.integers(0, len(s), size=(draws, len(s)))
        cols = s[idx]
        rows = np.broadcast_to(np.arange(draws)[:, None], cols.shape)
        np.add.at(W, (rows, cols), 1)
    assert (W.sum(1) == n).all()
    return W


def boot_sums(W, X):
    out = np.empty((W.shape[0], X.shape[1]))
    for i in range(0, W.shape[0], CHUNK):
        out[i:i + CHUNK] = W[i:i + CHUNK].astype(float) @ X
    return out


def rates(S):
    """S[..., 4] = (c0, d0, c1, d1) -> B, P, F (F as fraction)."""
    with np.errstate(divide='ignore', invalid='ignore'):
        B = 1000.0 * S[..., 0] / S[..., 1]
        P = 1000.0 * S[..., 2] / S[..., 3]
        F = 1.0 - P / B
    return B, P, F


def ci(x, level=0.95):
    x = np.asarray(x, dtype=float)
    good = x[np.isfinite(x)]
    q = [(1 - level) / 2, 1 - (1 - level) / 2]
    lo, hi = np.quantile(good, q)
    return dict(lo=float(lo), hi=float(hi), level=level, undefined_draws=int(x.size - good.size))


class Cond:
    """A condition with point estimates and aligned bootstrap draws."""

    def __init__(self, name, p, X, W):
        self.name, self.p = name, p
        S0 = X.sum(0)
        self.B, self.P, self.F = (float(v) for v in rates(S0))
        self.c0, self.d0, self.c1, self.d1 = (float(v) for v in S0)
        Sb = boot_sums(W, X)
        self.Bb, self.Pb, self.Fb = rates(Sb)
        self.n = X.shape[0]

    def summary(self):
        return dict(condition=self.name, p=self.p, n_seeds=self.n, control_acquisitions=int(self.c0),
                    control_patient_days=self.d0, policy_acquisitions=int(self.c1), policy_patient_days=self.d1,
                    control_rate=self.B, control_rate_ci95=ci(self.Bb), policy_rate=self.P,
                    F_pct=100 * self.F, F_pct_ci95=ci(100 * self.Fb))


def segment(lo, hi):
    """Two-point slopes between conditions lo and hi (same seed list)."""
    dlp = math.log(hi.p / lo.p)
    with np.errstate(divide='ignore', invalid='ignore'):
        dlB = math.log(hi.B / lo.B)
        dlBb = np.log(hi.Bb) - np.log(lo.Bb)
        beta_B = 100 * (hi.F - lo.F) / dlB
        beta_Bb = 100 * (hi.Fb - lo.Fb) / dlBb
        beta_p = 100 * (hi.F - lo.F) / dlp
        beta_pb = 100 * (hi.Fb - lo.Fb) / dlp
        gamma = dlB / dlp
        gammab = dlBb / dlp
    return dict(kind='two-point', conditions=[lo.name, hi.name], p=[lo.p, hi.p],
                control_rates=[lo.B, hi.B], F_pct=[100 * lo.F, 100 * hi.F],
                beta_B=beta_B, beta_B_draws=beta_Bb, beta_p=beta_p, beta_p_draws=beta_pb,
                gamma=gamma, gamma_draws=gammab)


def ols_slope(x, y):
    """Equal-weight OLS slope along the last axis."""
    xm = x.mean(-1, keepdims=True)
    ym = y.mean(-1, keepdims=True)
    return ((x - xm) * (y - ym)).sum(-1) / ((x - xm) ** 2).sum(-1)


def ols3(conds):
    lp = np.array([math.log(c.p) for c in conds])
    lB = np.array([math.log(c.B) for c in conds])
    F = np.array([100 * c.F for c in conds])
    with np.errstate(divide='ignore', invalid='ignore'):
        lBb = np.stack([np.log(c.Bb) for c in conds], -1)
    Fb = np.stack([100 * c.Fb for c in conds], -1)
    lpb = np.broadcast_to(lp, lBb.shape)
    return dict(kind='three-point OLS', conditions=[c.name for c in conds], p=[c.p for c in conds],
                control_rates=[c.B for c in conds], F_pct=list(F),
                beta_B=float(ols_slope(lB, F)), beta_B_draws=ols_slope(lBb, Fb),
                beta_p=float(ols_slope(lp, F)), beta_p_draws=ols_slope(lpb, Fb),
                gamma=float(ols_slope(lp, lB)), gamma_draws=ols_slope(lpb, lBb))


def slope_summary(s):
    out = {k: v for k, v in s.items() if not isinstance(v, np.ndarray)}
    for k in ('beta_B', 'beta_p', 'gamma'):
        out[k + '_ci95'] = ci(s[k + '_draws'])
    out['beta_B_per_10pct'] = s['beta_B'] * math.log(1.1)
    out['beta_B_per_10pct_ci95'] = ci(s['beta_B_draws'] * math.log(1.1))
    return out


# ------------------------------------------------------------------- main
def main():
    t0 = time.time()
    log('manuscript_v2 addenda analysis', time.strftime('%Y-%m-%d %H:%M:%S'))
    log('python', sys.version.split()[0], 'numpy', np.__version__)
    script_sha = sha(Path(__file__))
    log('script sha256', script_sha)

    # ---------------------------------------------------------- load + verify
    log('Loading and verifying inputs')
    sens_report = load(SENS_REPORT)
    sens_files = sorted((SENS_ROOT / 'completed').glob('job_*.json'))
    assert len(sens_files) == 1536 and {f.name for f in sens_files} == set(sens_report['input_sha256'])
    check_hashes('sensitivity records vs archived sensitivity report input_sha256',
                 [(SENS_ROOT / 'completed' / k, v) for k, v in sens_report['input_sha256'].items()])
    sb = load(AN_SEEDBOOK)['sensitivity']
    sens_seeds = sb['A'] + sb['B']
    sens = load_bank(sens_files, sens_seeds)

    rs_report = load(RS_REPORT)
    upb = rs_report['new_bank_bindings']['independent_upper']
    up_files = sorted((UPPER_ROOT / 'completed').glob('job_*.json'))
    assert len(up_files) == 256 and {'completed/' + f.name for f in up_files} == set(upb['input_sha256'])
    check_hashes('independent upper records vs rule-surface report binding',
                 [(UPPER_ROOT / k, v) for k, v in upb['input_sha256'].items()])

    hold_report = load(HOLD_REPORT)
    env_report = load(ENV_REPORT)

    def conf_pairs(report, prefix):
        pairs = []
        for b in report['banks']:
            site = [s for s in (prefix + '_cluster', prefix + '_local')
                    if sha(CONF_RESULTS / s / 'manifest.json') == b['manifest_sha256']]
            assert len(site) == 1
            pairs += [(CONF_RESULTS / site[0] / 'completed' / fn, h) for fn, h in b['completed_records_sha256'].items()]
        return pairs

    hold_pairs = conf_pairs(hold_report, 'holdout')
    assert len(hold_pairs) == 1536
    check_hashes('holdout records vs archived holdout report', hold_pairs)
    env_pairs = conf_pairs(env_report, 'environment')
    assert len(env_pairs) == 3072
    check_hashes('environment records vs archived environment report', env_pairs)

    csb = load(CONF_SEEDBOOK)['seed_blocks']
    hold_seeds = csb['holdout']['A'] + csb['holdout']['B']
    env_seeds = csb['environment']['A'] + csb['environment']['B']
    rs_index = load(RS_INDEX)
    core_seeds = rs_index['independent_upper_comparison_seeds']
    assert core_seeds == csb['holdout']['A'][:128] + csb['holdout']['B'][:128]
    hold = load_bank([p for p, _ in hold_pairs], hold_seeds)
    env = load_bank([p for p, _ in env_pairs], env_seeds)
    upper = load_bank(up_files, core_seeds)

    alt_sha = {}
    for line in ALT_LOCAL_SHA.read_text().splitlines():
        h, rel = line.split(None, 1)
        alt_sha[rel.strip()] = h
    alt_pol_files = sorted((ALT_EXEC / 'banks/policy/completed').glob('job_*.json'))
    alt_grid_files = sorted((ALT_EXEC / 'banks/grid/completed').glob('job_*.json'))
    assert len(alt_pol_files) == 1536 and len(alt_grid_files) == 1536
    check_hashes('alternative policy records vs collection local_sha256.txt',
                 [(f, alt_sha['./banks/policy/completed/' + f.name]) for f in alt_pol_files])
    check_hashes('alternative grid records vs collection local_sha256.txt',
                 [(f, alt_sha['./banks/grid/completed/' + f.name]) for f in alt_grid_files])
    asb = load(ALT_SEEDBOOK)
    alt = load_bank(alt_pol_files, asb['policy'])
    altgrid = load_bank(alt_grid_files, asb['grid'], paired=False)

    inputs = dict(
        sensitivity=dict(root=str(SENS_ROOT), n_records=len(sens_files), digest=input_digest(sens_files)),
        independent_upper=dict(root=str(UPPER_ROOT), n_records=len(up_files), digest=input_digest(up_files)),
        holdout=dict(root=str(CONF_RESULTS) + '/holdout_{cluster,local}', n_records=len(hold_pairs)),
        environment=dict(root=str(CONF_RESULTS) + '/environment_{cluster,local}', n_records=len(env_pairs)),
        alternative_policy=dict(root=str(ALT_EXEC / 'banks/policy'), n_records=len(alt_pol_files),
                                digest=input_digest(alt_pol_files)),
        alternative_grid=dict(root=str(ALT_EXEC / 'banks/grid'), n_records=len(alt_grid_files),
                              digest=input_digest(alt_grid_files)),
        other_files_sha256={str(p): sha(p) for p in (SENS_REPORT, AN_SEEDBOOK, RS_REPORT, RS_INDEX, CONF_SEEDBOOK,
                                                     HOLD_REPORT, ENV_REPORT, ALT_LOCAL_SHA, ALT_SEEDBOOK,
                                                     ALT_SELECTION, ALT_ANALYSIS, ORIG_SEL_INDEP, ORIG_SEL_EXPL,
                                                     COMMON_BASELINE, SUPPLEMENT)},
    )

    # ------------------------------------------------------------ bootstrap
    log('Building bootstrap weights')
    ss = np.random.SeedSequence(BOOT_SEED)
    g_sens, g_hold, g_env, g_alt, g_grid = (np.random.default_rng(c) for c in ss.spawn(5))
    W_sens = weights(g_sens, [np.arange(0, 128), np.arange(128, 256)])
    W_hold = weights(g_hold, [np.arange(0, 128), np.arange(128, 256), np.arange(256, 384), np.arange(384, 512)])
    core_idx = np.concatenate([np.arange(0, 128), np.arange(256, 384)])
    assert [hold_seeds[i] for i in core_idx] == core_seeds
    W_core = W_hold[:, core_idx]
    W_env = weights(g_env, [np.arange(0, 256), np.arange(256, 512)])
    W_alt = weights(g_alt, [np.arange(512)])
    W_grid = weights(g_grid, [np.arange(128)])

    # ------------------------------------------------------------ conditions
    log('Computing condition estimates')
    C = {}
    for sid in ('clustered_low', 'clustered_center', 'clustered_high', 'imperfect_low', 'imperfect_center',
                'imperfect_high'):
        C['sens:' + sid] = Cond('sensitivity ' + sid, sens[sid]['p'], sens[sid]['X'], W_sens)
    hi_ind = hold['independent_target_candidate']
    core_rows = np.array([hold_seeds.index(s) for s in core_seeds])
    C['ind:center_matched256'] = Cond('independent center (matched 256)', hi_ind['p'], hi_ind['X'][core_rows], W_core)
    C['ind:upper'] = Cond('independent upper (matched 256)', upper['independent_upper']['p'],
                          upper['independent_upper']['X'], W_core)
    for sid in ('independent_low075_env1', 'independent_low100_env1', 'independent_low125_env1'):
        C['env:' + sid] = Cond('environment full model ' + sid, env[sid]['p'], env[sid]['X'], W_env)
    hold_map = dict(independent='independent_target_candidate', clustered='persistent_target_candidate',
                    imperfect='imperfect_target_candidate')
    for prof, sid in hold_map.items():
        C['hold:' + prof] = Cond('original holdout ' + prof, hold[sid]['p'], hold[sid]['X'], W_hold)
    for prof in ('independent', 'clustered', 'imperfect'):
        C['alt:' + prof] = Cond('alternative holdout ' + prof, alt[prof]['p'], alt[prof]['X'], W_alt)

    # reproduction checks vs archived point values
    repro = []

    def rcheck(label, got, want, tol):
        ok = abs(got - want) <= tol
        repro.append(dict(check=label, recomputed=got, archived=want, tolerance=tol, pass_=ok))
        log(f'  repro {label}: {got:.10g} vs {want:.10g} {"ok" if ok else "FAIL"}')
        return ok

    for sc in sens_report['scenarios']:
        c = C['sens:' + sc['scenario']['id']]
        rcheck('sensitivity ' + sc['scenario']['id'] + ' F', c.F, sc['F'], 1e-12)
        rcheck('sensitivity ' + sc['scenario']['id'] + ' control rate', c.B, sc['control_rate'], 1e-12)
    cb = load(COMMON_BASELINE)['primary']['profiles']
    rcheck('independent matched center B (common baseline B_lo)', C['ind:center_matched256'].B, cb['independent']['B_lo'], 1e-12)
    rcheck('independent upper B (common baseline B_hi)', C['ind:upper'].B, cb['independent']['B_hi'], 1e-12)
    rcheck('independent matched center F', C['ind:center_matched256'].F, cb['independent']['F_lo'], 1e-12)
    rcheck('independent upper F', C['ind:upper'].F, cb['independent']['F_hi'], 1e-12)
    for sc in load(ALT_ANALYSIS)['scenarios']:
        c = C['alt:' + sc['scenario']['id']]
        rcheck('alternative ' + sc['scenario']['id'] + ' F', c.F, sc['F'], 1e-12)
        rcheck('alternative ' + sc['scenario']['id'] + ' control rate', c.B, sc['control_rate'], 1e-12)
    for prof, (b, f) in dict(independent=(2.83, 67.6), clustered=(2.69, 55.9), imperfect=(2.80, 27.0)).items():
        rcheck(f'original holdout {prof} control rate (manuscript, 2 dp)', round(C['hold:' + prof].B, 2), b, 1e-9)
        rcheck(f'original holdout {prof} F (manuscript, 1 dp)', round(100 * C['hold:' + prof].F, 1), f, 1e-9)

    # ------------------------------------------------------------ elasticity
    log('Elasticity')
    S = {}
    S['clustered'] = dict(
        primary=segment(C['sens:clustered_center'], C['sens:clustered_high']),
        low_center=segment(C['sens:clustered_low'], C['sens:clustered_center']),
        ols3=ols3([C['sens:clustered_low'], C['sens:clustered_center'], C['sens:clustered_high']]))
    S['imperfect'] = dict(
        primary=segment(C['sens:imperfect_center'], C['sens:imperfect_high']),
        low_center=segment(C['sens:imperfect_low'], C['sens:imperfect_center']),
        ols3=ols3([C['sens:imperfect_low'], C['sens:imperfect_center'], C['sens:imperfect_high']]))
    e1 = [C['env:independent_low075_env1'], C['env:independent_low100_env1'], C['env:independent_low125_env1']]
    S['independent'] = dict(
        primary=segment(C['ind:center_matched256'], C['ind:upper']),
        env_low_mid=segment(e1[0], e1[1]),
        env_mid_high=segment(e1[1], e1[2]),
        env_ols3=ols3(e1))
    slope_labels = dict(
        primary='center to 1.25 x center (brackets 3.0)',
        low_center='0.75 x center to center',
        ols3='three-point OLS (0.75, 1, 1.25 x center)',
        env_low_mid='environment full model, p .0810 to .1080 (separate 512-seed bank)',
        env_mid_high='environment full model, p .1080 to .1350 (separate 512-seed bank)',
        env_ols3='environment full model, three-point OLS (B 0.40 to 2.29; below the target region)',
        piecewise='piecewise-local: low->center slope below the center, center->high slope above it '
                  '(independent: center->high on both sides, no paired lower point)')

    # consistency: primary slope applied at the sensitivity/matched center reproduces common-baseline F*
    consistency = {}
    for prof in ('independent', 'clustered', 'imperfect'):
        s = S[prof]['primary']
        Fc, Bc = s['F_pct'][0], s['control_rates'][0]
        f_arith = Fc + s['beta_B'] * math.log(TARGET / Bc)
        consistency[prof] = dict(F_center_pct=Fc, B_center=Bc, F_at_3_arithmetic_pct=f_arith,
                                 common_baseline_F_star_pct=100 * cb[prof]['F_star'],
                                 difference_pp=f_arith - 100 * cb[prof]['F_star'])

    # ------------------------------------------------------------ bias
    log('Bias from baseline shortfall')
    observed = dict(
        original=dict(independent=C['hold:independent'], clustered=C['hold:clustered'], imperfect=C['hold:imperfect']),
        alternative=dict(independent=C['alt:independent'], clustered=C['alt:clustered'], imperfect=C['alt:imperfect']))
    stated = dict(original=dict(independent=2.83, clustered=2.69, imperfect=2.80),
                  alternative=dict(independent=2.90, clustered=2.62, imperfect=2.63))
    center_cond = dict(independent=C['ind:center_matched256'], clustered=C['sens:clustered_center'],
                       imperfect=C['sens:imperfect_center'])

    def piecewise(prof, x, draws):
        """F_hat(x) - F_hat(ln 3), F_hat piecewise linear in ln B through the seed-paired points:
        low->center slope left of the center, center->high slope right of it (end segments extrapolate).
        Independent has no paired point below its matched center, so the center->high slope is used on
        both sides (identical to the primary estimator). x: ln(observed baseline), float or draws."""
        c = center_cond[prof]
        b_ch = S[prof]['primary']['beta_B_draws'] if draws else S[prof]['primary']['beta_B']
        if 'low_center' in S[prof]:
            b_lc = S[prof]['low_center']['beta_B_draws'] if draws else S[prof]['low_center']['beta_B']
        else:
            b_lc = b_ch
        lnc = np.log(c.Bb) if draws else math.log(c.B)

        def g(v):
            return np.where(v < lnc, b_lc * (v - lnc), b_ch * (v - lnc))

        with np.errstate(invalid='ignore'):
            out = g(x) - g(math.log(TARGET))
        return out if draws else float(out)

    bias = {}
    joint_draws = {}
    for wf, profs in observed.items():
        bias[wf] = {}
        for prof, cond in profs.items():
            assert round(cond.B, 2) == stated[wf][prof], (wf, prof, cond.B)
            s_pt = math.log(cond.B / TARGET)
            s_b = np.log(cond.Bb / TARGET)
            rows = {}
            # piecewise-local estimator (recommended local reading; equals primary when the observed
            # baseline and 3.0 lie on the same side of the center)
            pw_pt = piecewise(prof, math.log(cond.B), False)
            pw_beta_only = piecewise(prof, math.log(cond.B), True)
            pw_joint = piecewise(prof, np.log(cond.Bb), True)
            rows['piecewise'] = dict(
                slope_label=slope_labels['piecewise'], bias_pp=pw_pt,
                bias_pp_ci95_beta_only=ci(pw_beta_only), bias_pp_ci95_joint=ci(pw_joint),
                observed_below_center=bool(cond.B < center_cond[prof].B),
                center_rate=center_cond[prof].B,
                F_observed_pct=100 * cond.F, F_adjusted_to_3_pct=100 * cond.F - pw_pt,
                F_adjusted_to_3_pct_ci95_joint=ci(100 * cond.Fb - pw_joint))
            joint_draws[(wf, prof, 'piecewise')] = (pw_pt, pw_joint)
            for key, sl in S[prof].items():
                b_pt = sl['beta_B'] * s_pt
                rows[key] = dict(slope_label=slope_labels[key], beta_B=sl['beta_B'],
                                 bias_pp=b_pt, bias_pp_ci95_beta_only=ci(sl['beta_B_draws'] * s_pt))
                if key == 'primary':
                    # joint: shortfall resampled too. Draws are aligned where seeds are shared
                    # (original independent: matched center seeds are a subset of the holdout);
                    # otherwise the streams are independent.
                    jb = sl['beta_B_draws'] * s_b
                    rows[key]['bias_pp_ci95_joint'] = ci(jb)
                    Fadj_b = 100 * cond.Fb - jb
                    rows[key]['F_observed_pct'] = 100 * cond.F
                    rows[key]['F_observed_pct_ci95'] = ci(100 * cond.Fb)
                    rows[key]['F_adjusted_to_3_pct'] = 100 * cond.F - b_pt
                    rows[key]['F_adjusted_to_3_pct_ci95_joint'] = ci(Fadj_b)
                    joint_draws[(wf, prof, 'primary')] = (b_pt, jb)
            bias[wf][prof] = dict(control_rate=cond.B, control_rate_ci95=ci(cond.Bb), stated_rate=stated[wf][prof],
                                  shortfall_ln=s_pt, shortfall_pct=100 * (cond.B / TARGET - 1),
                                  shortfall_ln_ci95=ci(s_b), by_slope=rows,
                                  transfer_assumption=(wf == 'alternative'))
            envs = [r['bias_pp'] for r in rows.values()]
            bias[wf][prof]['point_range_across_slope_estimators_pp'] = [min(envs), max(envs)]

    # contrast-level implied bias. Draws are index-aligned: shared seed streams keep their covariance
    # (clustered/imperfect slopes share the sensitivity stream; observed shortfalls and F share the holdout or
    # alternative stream within a workflow); separate streams are independent.
    pairs3 = (('independent', 'clustered'), ('independent', 'imperfect'), ('clustered', 'imperfect'))
    contrast_bias = {}
    for wf, profs in observed.items():
        contrast_bias[wf] = {}
        for a, b in pairs3:
            name = f'{a}-{b}'
            obs_pt = 100 * (profs[a].F - profs[b].F)
            obs_b = 100 * (profs[a].Fb - profs[b].Fb)
            row = dict(observed_contrast_pp=obs_pt, observed_contrast_ci95=ci(obs_b),
                       observed_contrast_ci98333=ci(obs_b, 1 - 0.05 / 3))
            for key in ('primary', 'piecewise'):
                pa, da = joint_draws[(wf, a, key)]
                pb, db = joint_draws[(wf, b, key)]
                row[key] = dict(contrast_bias_pp=pa - pb, contrast_bias_ci95_joint=ci(da - db),
                                adjusted_contrast_pp=obs_pt - (pa - pb),
                                adjusted_contrast_ci95_joint=ci(obs_b - (da - db)),
                                adjusted_contrast_ci98333_joint=ci(obs_b - (da - db), 1 - 0.05 / 3))
            contrast_bias[wf][name] = row

    # ------------------------------------------------------------ selection propagation
    log('Selection propagation')
    si = load(ORIG_SEL_INDEP)['cohorts']['pooled']
    sx = [c for c in load(ORIG_SEL_EXPL)['cohorts'] if c['cohort'] == 'pooled64'][0]['profiles']
    asel = load(ALT_SELECTION)['profiles']

    def interval_from(prof):
        b = sx[prof]['bootstrap']
        return b['conditional_valid_draw_interval95']

    sel_inputs = dict(
        original=dict(
            independent=dict(p=si['selected_p'], lo=si['conditional_selected_p_uncertainty']['interval95'][0],
                             hi=si['conditional_selected_p_uncertainty']['interval95'][1],
                             source='Supplement Table S8, pooled 128-seed row (conditional 95%; equals fallback-inclusive)',
                             file=str(ORIG_SEL_INDEP)),
            clustered=dict(p=C['hold:clustered'].p, lo=interval_from('persistent')[0], hi=interval_from('persistent')[1],
                           source='Supplement Table S14, pooled64 row (conditional 95%; 0 failed draws)',
                           file=str(ORIG_SEL_EXPL)),
            imperfect=dict(p=C['hold:imperfect'].p, lo=interval_from('imperfect')[0], hi=interval_from('imperfect')[1],
                           source='Supplement Table S14, pooled64 row (conditional 95%; 0 failed draws)',
                           file=str(ORIG_SEL_EXPL))),
        alternative={prof: dict(p=asel[prof]['selected_p'],
                                lo=asel[prof]['selection_uncertainty']['conditional_on_success95'][0],
                                hi=asel[prof]['selection_uncertainty']['conditional_on_success95'][1],
                                source='verification RESULTS.md section 3 (selection.json conditional 95%; 0 failed draws)',
                                lower_bracket_draws=asel[prof]['selection_uncertainty']['selected_bracket_frequencies'],
                                file=str(ALT_SELECTION))
                     for prof in ('independent', 'clustered', 'imperfect')})
    # supplement display check
    disp = dict(independent=(0.14328, 0.14004, 0.14792), clustered=(0.07602, 0.07346, 0.07873),
                imperfect=(0.07484, 0.07243, 0.07785))
    for prof, (p, lo, hi) in disp.items():
        d = sel_inputs['original'][prof]
        for nm, got, want in (('p', d['p'], p), ('lo', d['lo'], lo), ('hi', d['hi'], hi)):
            rcheck(f'original selection {prof} {nm} vs supplement display', round(got, 5), want, 1e-12)
    for prof, (p, lo, hi) in dict(independent=(0.165111, 0.15999, 0.17033), clustered=(0.080633, 0.07801, 0.08332),
                                  imperfect=(0.076144, 0.07338, 0.07862)).items():
        d = sel_inputs['alternative'][prof]
        rcheck(f'alternative selection {prof} p vs RESULTS.md', round(d['p'], 6), p, 1e-12)
        rcheck(f'alternative selection {prof} lo vs RESULTS.md', round(d['lo'], 5), lo, 1e-12)
        rcheck(f'alternative selection {prof} hi vs RESULTS.md', round(d['hi'], 5), hi, 1e-12)
    assert abs(sel_inputs['original']['independent']['p'] - C['hold:independent'].p) < 1e-15

    # alternative-workflow local log-log slope of baseline on p (grid m1 -> m2)
    gamma_alt = {}
    for prof in ('independent', 'clustered', 'imperfect'):
        g1, g2 = altgrid[prof + '_m1'], altgrid[prof + '_m2']
        S1, S2 = g1['X'].sum(0), g2['X'].sum(0)
        B1, B2 = 1000 * S1[0] / S1[1], 1000 * S2[0] / S2[1]
        b1 = boot_sums(W_grid, g1['X'])
        b2 = boot_sums(W_grid, g2['X'])
        with np.errstate(divide='ignore', invalid='ignore'):
            gb = (np.log(1000 * b2[:, 0] / b2[:, 1]) - np.log(1000 * b1[:, 0] / b1[:, 1])) / math.log(g2['p'] / g1['p'])
        gamma_alt[prof] = dict(B_m1=float(B1), B_m2=float(B2), p_m1=g1['p'], p_m2=g2['p'],
                               gamma=float(math.log(B2 / B1) / math.log(g2['p'] / g1['p'])), gamma_ci95=ci(gb),
                               gamma_draws=gb)

    hw = dict(original=dict(independent=(64.1, 70.8), clustered=(52.1, 59.5), imperfect=(21.6, 32.1)),
              alternative=dict(independent=(74.1, 78.7), clustered=(49.7, 57.5), imperfect=(22.2, 31.9)))
    selection = {}
    for wf in ('original', 'alternative'):
        selection[wf] = {}
        for prof in ('independent', 'clustered', 'imperfect'):
            d = sel_inputs[wf][prof]
            dlo, dhi = math.log(d['lo'] / d['p']), math.log(d['hi'] / d['p'])
            sl = S[prof]['primary']
            bp, bpb = sl['beta_p'], sl['beta_p_draws']
            pt = sorted([bp * dlo, bp * dhi])
            lo_b = np.minimum(bpb * dlo, bpb * dhi)
            hi_b = np.maximum(bpb * dlo, bpb * dhi)
            bci = ci(bpb)
            env_vals = [b * dd for b in (bci['lo'], bci['hi'], bp) for dd in (dlo, dhi, 0.0)]
            mc_hw = (hw[wf][prof][1] - hw[wf][prof][0]) / 2
            sel_hw = max(abs(pt[0]), abs(pt[1]))
            env_hw = max(abs(min(env_vals)), abs(max(env_vals)))
            row = dict(p_selected=d['p'], p_ci95=[d['lo'], d['hi']], source=d['source'],
                       dlnp=[dlo, dhi], dlnp_pct=[100 * (d['lo'] / d['p'] - 1), 100 * (d['hi'] / d['p'] - 1)],
                       dF_dlnp=bp, dF_dlnp_ci95=bci, dF_dlnp_source='original workflow, ' + slope_labels['primary'],
                       implied_F_shift_pp_point=pt,
                       implied_F_shift_pp_with_slope_uncertainty=[float(np.quantile(lo_b[np.isfinite(lo_b)], 0.025)),
                                                                  float(np.quantile(hi_b[np.isfinite(hi_b)], 0.975))],
                       implied_F_shift_pp_envelope=[min(env_vals), max(env_vals)],
                       conditional_MC_F_ci95_pct=list(hw[wf][prof]), conditional_MC_half_width_pp=mc_hw,
                       selection_half_width_pp_point=sel_hw, selection_half_width_pp_envelope=env_hw,
                       quadrature_half_width_pp_point=math.hypot(mc_hw, sel_hw),
                       quadrature_half_width_pp_envelope=math.hypot(mc_hw, env_hw),
                       transfer_assumption=(wf == 'alternative'))
            if wf == 'alternative':
                # sensitivity: transferred dF/dlnB x the alternative workflow's own grid log-log slope
                ga = gamma_alt[prof]
                bp2 = sl['beta_B'] * ga['gamma']
                bp2b = sl['beta_B_draws'] * ga['gamma_draws']
                row['sensitivity_beta_B_times_alt_gamma'] = dict(
                    alt_grid_gamma=ga['gamma'], alt_grid_gamma_ci95=ga['gamma_ci95'],
                    alt_grid_bracket=[ga['p_m1'], ga['p_m2']], alt_grid_rates=[ga['B_m1'], ga['B_m2']],
                    dF_dlnp=bp2, dF_dlnp_ci95=ci(bp2b),
                    implied_F_shift_pp_point=sorted([bp2 * dlo, bp2 * dhi]))
            selection[wf][prof] = row

    # ------------------------------------------------------------ routes
    log('Route labels')
    route_out = {}
    for wf, bank, W, idmap in (('original', hold, W_hold, hold_map),
                               ('alternative', alt, W_alt, {k: k for k in ('independent', 'clustered', 'imperfect')})):
        route_out[wf] = {}
        for prof, sid in idmap.items():
            b = bank[sid]
            R0, R1 = b['routes0'], b['routes1']
            D0, D1 = b['X'][:, 1].sum(), b['X'][:, 3].sum()
            tot0, tot1 = R0.sum(), R1.sum()
            assert tot0 == b['X'][:, 0].sum() and tot1 == b['X'][:, 2].sum()
            R0b, R1b = boot_sums(W, R0), boot_sums(W, R1)
            D0b = boot_sums(W, b['X'][:, [1]])[:, 0]
            D1b = boot_sums(W, b['X'][:, [3]])[:, 0]
            labels = {}
            groups = dict(patient=['patient'], own_current_0=['own0'], own_current_pos=['ownpos'],
                          own_current_all=['own0', 'ownpos'], other_zone=['other'],
                          carry_over=['carry0', 'carrypos'], carry_over_0=['carry0'], carry_over_pos=['carrypos'],
                          all_except_carry=['patient', 'own0', 'ownpos', 'other'],
                          paper_metric_site_0=['own0', 'carry0'],
                          all=list(ROUTE_KEYS))
            for g, keys in groups.items():
                ix = [ROUTE_KEYS.index(k) for k in keys]
                c0, c1 = float(R0[:, ix].sum()), float(R1[:, ix].sum())
                c0b, c1b = R0b[:, ix].sum(1), R1b[:, ix].sum(1)
                with np.errstate(divide='ignore', invalid='ignore'):
                    red_count = 100 * (1 - c1 / c0) if c0 > 0 else None
                    red_count_b = 100 * (1 - c1b / c0b)
                    red_rate = 100 * (1 - (c1 / D1) / (c0 / D0)) if c0 > 0 else None
                    contrib = 100 * (c0 / D0 - c1 / D1) / (tot0 / D0)
                    contrib_b = 100 * (c0b / D0b - c1b / D1b) / (R0b.sum(1) / D0b)
                labels[g] = dict(control=int(c0), policy=int(c1),
                                 control_share_pct=100 * c0 / tot0, policy_share_pct=100 * c1 / tot1,
                                 reduction_count_pct=red_count,
                                 reduction_count_pct_ci95=(ci(red_count_b) if c0 > 0 else None),
                                 reduction_rate_pct=red_rate,
                                 contribution_to_F_pp=contrib, contribution_to_F_pp_ci95=ci(contrib_b),
                                 small_count_flag=bool(c0 < 30))
            route_out[wf][prof] = dict(records=int(b['n']), control_patient_days=float(D0),
                                       policy_patient_days=float(D1),
                                       patient_days_relative_difference=float(D1 / D0 - 1), labels=labels)
            # sanity: contributions sum to F
            Fsum = sum(labels[g]['contribution_to_F_pp'] for g in ('patient', 'own_current_all', 'other_zone', 'carry_over'))
            Fd = 100 * (1 - (tot1 / D1) / (tot0 / D0))
            assert abs(Fsum - Fd) < 1e-9, (wf, prof, Fsum, Fd)
    # reproduce the verification report's route counts
    vr = dict(original=dict(independent=(658, 143, 344, 80, 923, 228, 450, 319),
                            clustered=(1361, 590, 324, 130, 495, 222, 75, 53),
                            imperfect=(1467, 1047, 331, 231, 492, 382, 62, 58)),
              alternative=dict(independent=(824, 195, 597, 137, 952, 213, 65, 28),
                               clustered=(1342, 640, 371, 142, 466, 230, 22, 7),
                               imperfect=(1404, 1029, 347, 238, 438, 330, 16, 7)))
    for wf, d in vr.items():
        for prof, v in d.items():
            L = route_out[wf][prof]['labels']
            got = (L['patient']['control'], L['patient']['policy'], L['other_zone']['control'], L['other_zone']['policy'],
                   L['own_current_all']['control'], L['own_current_all']['policy'], L['carry_over']['control'],
                   L['carry_over']['policy'])
            rcheck(f'route counts {wf} {prof} vs verification RESULTS.md section 7', float(got == v), 1.0, 0)

    # ------------------------------------------------------------ compliance table (S5)
    log('Compliance-range table from supplement Table S5')
    tex = SUPPLEMENT.read_text()
    seg_tex = tex.split(r'\label{tab:author-coefficients}')[1].split(r'\bottomrule')[0]
    s5 = []
    for line in seg_tex.splitlines():
        m = re.match(r'^\s*([0-9.]+) & ([0-9.]+) & ([0-9.]+) & ([0-9.]+) & ([0-9.]+) & ([0-9.]+) & '
                     r'([0-9.]+) \[([0-9.]+), ([0-9.]+)\] & ([0-9.]+)\\\\', line)
        if m:
            v = [float(x) for x in m.groups()]
            s5.append(dict(h0=v[0], h1=v[1], rho=v[2], e=v[3], K0=v[4], K1=v[5], F_star_pct=v[6],
                           F_star_ci95=[v[7], v[8]], q2_response_pct=v[9]))
    assert len(s5) == 22, len(s5)
    prof_of = {(0.0, 1.0): 'independent', (0.5, 1.0): 'clustered', (0.0, 0.8): 'imperfect'}
    compliance = []
    for hp in ((0.80, 0.90), (0.90, 0.95), (0.95, 0.99)):
        row = dict(h0=hp[0], h1=hp[1])
        for r in s5:
            if (r['h0'], r['h1']) == hp and (r['rho'], r['e']) in prof_of:
                row[prof_of[(r['rho'], r['e'])]] = r
        assert all(k in row for k in prof_of.values()), hp
        row['independent_minus_clustered_pp'] = row['independent']['F_star_pct'] - row['clustered']['F_star_pct']
        row['independent_minus_imperfect_pp'] = row['independent']['F_star_pct'] - row['imperfect']['F_star_pct']
        row['clustered_minus_imperfect_pp'] = row['clustered']['F_star_pct'] - row['imperfect']['F_star_pct']
        row['alternative_workflow_F_star'] = 'PENDING: workflow_coefficient_2026-10-21'
        compliance.append(row)

    # ------------------------------------------------------------ assemble
    def strip(d):
        if isinstance(d, dict):
            return {k: strip(v) for k, v in d.items() if not isinstance(v, np.ndarray)}
        if isinstance(d, list):
            return [strip(v) for v in d]
        if isinstance(d, (np.floating,)):
            return float(d)
        return d

    out = dict(
        schema='manuscript_v2_addenda/1',
        date='2026-10-22', folder_date='2026-10-21',
        status='analysis-only addenda from saved records; no simulation',
        script=str(Path(__file__)), script_sha256=script_sha,
        runtime=dict(python=sys.version.split()[0], numpy=np.__version__),
        bootstrap=dict(draws=DRAWS, seed=BOOT_SEED, streams='SeedSequence(20261021).spawn(5): sensitivity, '
                       'holdout (4 strata; upper uses cores), environment, alternative policy, alternative grid',
                       interval='percentile 95%'),
        inputs=inputs, hash_checks=HASH_CHECKS, reproduction_checks=repro,
        conditions={k: v.summary() for k, v in C.items()},
        elasticity=dict(
            definition='beta_B = dF/dln(B) in percentage points of F per unit ln(baseline rate); '
                       'beta_p = dF/dln(p); gamma = dln(B)/dln(p). Two-point slopes are finite differences '
                       'between seed-paired conditions; OLS slopes are equal-weight.',
            slope_labels=slope_labels,
            slopes={prof: {k: slope_summary(v) for k, v in d.items()} for prof, d in S.items()},
            consistency_with_common_baseline=consistency),
        bias=dict(
            definition='bias = F(B_obs) - F(3.0) ~= beta_B * ln(B_obs/3), in percentage points. '
                       'F at 3.0 ~= F_obs - bias. Negative shortfall (B_obs < 3) with negative beta_B gives a '
                       'positive bias (observed F too high), and vice versa.',
            interval_notes='beta_only: shortfall fixed at its point value, slope resampled. joint: shortfall and '
                           'slope both resampled (aligned draws where seeds are shared).',
            results=strip(bias),
            contrast_definition='contrast bias = bias(a) - bias(b); adjusted contrast = observed contrast - '
                                'contrast bias. Post hoc and descriptive; it relies on the same local '
                                'log-linearity and (alternative workflow) slope-transfer assumptions.',
            contrasts=strip(contrast_bias)),
        selection=dict(
            definition='implied F shift = dF/dln p x ln(p_bound/p_selected), in percentage points. '
                       'Quadrature half-widths combine the archived conditional MC half-width and the selection '
                       'half-width as if independent and symmetric (descriptive only).',
            inputs=strip(sel_inputs), results=strip(selection),
            alt_grid_gamma=strip(gamma_alt)),
        routes=dict(
            definition='Passive labels of the pickup that began the infecting hand contamination episode. '
                       'own_current = own-zone site pickup, contamination from the current stay; carry_over = '
                       'own_zone_environment_pickup_after_readmission (site lineage predates the current '
                       'admission); suffix _0/_pos = zero / at least one intervening hygiene opportunity. '
                       'reduction = 1 - policy/control; contribution_to_F = (control rate - policy rate of label) '
                       '/ total control rate, which sums to F over the four labels.',
            results=strip(route_out)),
        compliance_range=dict(source='supplement.tex Table S5 (tab:author-coefficients), parsed', rows_all=s5,
                              by_h_pair=compliance),
        study_flow=[dict(r, n=r['n'].format(n_records=f"{sum(v['n_records'] for k, v in inputs.items() if isinstance(v, dict) and 'n_records' in v):,}"))
                    for r in STUDY_FLOW],
        placeholders=dict(workflow_coefficient=PENDING),
        elapsed_seconds=None,
    )
    # ------------------------------------------------------------ headline (derived; rendered in ADDENDA.md)
    PR3 = ('independent', 'clustered', 'imperfect')
    loc = ('piecewise', 'primary')
    archived_contrasts = dict(original=[11.7, 40.6, 28.9], alternative=[22.8, 49.2, 26.4])
    hws = [(hw[wf][p][1] - hw[wf][p][0]) / 2 for wf in hw for p in PR3]
    jl = [x for wf in bias for p in PR3 for k in loc
          for x in (bias[wf][p]['by_slope'][k]['bias_pp_ci95_joint']['lo'],
                    bias[wf][p]['by_slope'][k]['bias_pp_ci95_joint']['hi'])]
    cbj = [x for wf in contrast_bias for c in contrast_bias[wf] for k in loc
           for x in (contrast_bias[wf][c][k]['contrast_bias_ci95_joint']['lo'],
                     contrast_bias[wf][c][k]['contrast_bias_ci95_joint']['hi'])]
    out['headline'] = dict(
        max_abs_bias_pp_local_estimators=max(abs(bias[wf][p]['by_slope'][k]['bias_pp'])
                                             for wf in bias for p in PR3 for k in loc),
        max_abs_bias_pp_any_estimator=max(abs(v['bias_pp']) for wf in bias for p in PR3
                                          for v in bias[wf][p]['by_slope'].values()),
        joint_ci95_range_pp_local_estimators=[min(jl), max(jl)],
        archived_MC_F_half_width_pp_range=[min(hws), max(hws)],
        archived_MC_F_ci95_pct=hw,
        archived_primary_contrasts_pp=archived_contrasts,
        all_slope_ci95_include_zero=all(ci(v['beta_B_draws'])['lo'] < 0 < ci(v['beta_B_draws'])['hi']
                                        for d in S.values() for v in d.values()),
        max_abs_contrast_bias_pp_local=max(abs(contrast_bias[wf][c][k]['contrast_bias_pp'])
                                           for wf in contrast_bias for c in contrast_bias[wf] for k in loc),
        contrast_bias_joint_ci95_range_pp_local=[min(cbj), max(cbj)],
        all_adjusted_contrasts_exclude_zero_98333=all(
            contrast_bias[wf][c][k]['adjusted_contrast_ci98333_joint']['lo'] > 0
            for wf in contrast_bias for c in contrast_bias[wf] for k in loc),
        selected_p_interval_pct_range=[min(selection[wf][p]['dlnp_pct'][0] for wf in selection for p in PR3),
                                       max(selection[wf][p]['dlnp_pct'][1] for wf in selection for p in PR3)],
        max_abs_selection_F_shift_pp_point_slope=max(abs(x) for wf in selection for p in PR3
                                                     for x in selection[wf][p]['implied_F_shift_pp_point']),
        max_abs_selection_F_shift_pp_with_slope_uncertainty=max(
            abs(x) for wf in selection for p in PR3
            for x in selection[wf][p]['implied_F_shift_pp_with_slope_uncertainty']),
        holdout_F_ci95_max_abs_diff_from_published_pp=max(
            abs(C['hold:' + p].summary()['F_pct_ci95'][e] - hw['original'][p][i])
            for p in PR3 for i, e in enumerate(('lo', 'hi'))),
    )
    out['elapsed_seconds'] = round(time.time() - t0, 1)
    out['log'] = LOG
    fails = [r for r in repro if not r['pass_']]
    out['reproduction_failures'] = len(fails)
    (HERE / 'addenda.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    log('wrote', HERE / 'addenda.json', 'reproduction failures:', len(fails), 'elapsed', out['elapsed_seconds'])


if __name__ == '__main__':
    main()
