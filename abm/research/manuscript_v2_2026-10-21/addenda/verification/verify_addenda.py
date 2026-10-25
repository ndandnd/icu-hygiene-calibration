"""Independent verification of ADDENDA.md sections 1-3 (written without reading
addenda_analysis.py or render_addenda_md.py).

Reads raw records only (plus seedbooks/configs for seed order, and archived
reports only for MC half-widths and p intervals). Writes verify_results.json
next to this file. No simulation; read-only on all inputs.
"""
from pathlib import Path
import glob
import json
import math
import os
import sys

import numpy as np

ROOT = str(Path(__file__).resolve().parents[5])
ABM = ROOT + '/abm'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'verify_results.json')
D = 20000
RNG_SEED = 777001  # deliberately different from the addenda's 20261021
TARGET = 3.0

LABELS = ['patient', 'own0', 'ownpos', 'other', 'carry0', 'carrypos']


def label_of(key):
    route, k = key.split('|')
    k = int(k)
    if route == 'patient_pickup':
        return 'patient'
    if route == 'own_zone_environment_pickup':
        return 'own0' if k == 0 else 'ownpos'
    if route == 'other_zone_environment_pickup':
        return 'other'
    if route == 'own_zone_environment_pickup_after_readmission':
        return 'carry0' if k == 0 else 'carrypos'
    raise ValueError('unknown route key ' + key)


def read_dir(dirs, keep=None):
    rows = {}
    for d in dirs:
        files = sorted(glob.glob(d + '/*.json'))
        for f in files:
            r = json.load(open(f))
            sid = r['scenario']['id']
            if keep is not None and sid not in keep:
                continue
            shared = 'shared_with' in r['arms']['policy']
            rec = dict(sid=sid, p=r['scenario']['p'], seed=r['seed'],
                       h=(r['arms']['baseline']['h'], None if shared else r['arms']['policy']['h']))
            for arm in ('baseline', 'policy'):
                src = 'baseline' if (arm == 'policy' and shared) else arm
                ob = r['arms'][src]['observation']
                rao = ob['response_counts']['acquisitions_by_route_and_opportunities']
                assert sum(rao.values()) == ob['acquired_colonizations'], f
                lab = dict.fromkeys(LABELS, 0)
                for k, v in rao.items():
                    lab[label_of(k)] += v
                rec[arm] = (ob['acquired_colonizations'], ob['patient_days'], lab)
            key = (sid, rec['seed'])
            assert key not in rows, ('duplicate', key)
            rows[key] = rec
    return rows


def cond_arrays(rows, sid, seeds):
    recs = [rows[(sid, s)] for s in seeds]
    ps = {r['p'] for r in recs}
    assert len(ps) == 1
    X = np.array([[r['baseline'][0], r['baseline'][1], r['policy'][0], r['policy'][1]] for r in recs], float)
    L0 = np.array([[r['baseline'][2][l] for l in LABELS] for r in recs], float)
    L1 = np.array([[r['policy'][2][l] for l in LABELS] for r in recs], float)
    return dict(X=X, L0=L0, L1=L1, p=ps.pop(), n=len(seeds), h=recs[0]['h'])


def weights(rng, strata_sizes):
    """D x n multinomial resampling counts, independent within each stratum."""
    blocks = [rng.multinomial(m, np.full(m, 1.0 / m), size=D).astype(float) for m in strata_sizes]
    return np.concatenate(blocks, axis=1)


def BPF(X, W=None):
    t = X.sum(0) if W is None else W @ X
    t = np.asarray(t)
    B = 1000 * t[..., 0] / t[..., 1]
    P = 1000 * t[..., 2] / t[..., 3]
    return B, P, 1 - P / B


def q(a, lo=2.5, hi=97.5):
    return [float(np.percentile(a, lo)), float(np.percentile(a, hi))]


def main():
    res = {}
    rng = np.random.default_rng(RNG_SEED)
    # ---------------------------------------------------------------- loading
    sens_book = json.load(open(ABM + '/research/author_notes_2026-10-10/seedbook.json'))['sensitivity']
    sens_seeds = sens_book['A'] + sens_book['B']
    sens_cfg = json.load(open(ABM + '/research/author_notes_2026-10-10/configs/sensitivity.json'))
    assert sens_cfg['seeds'] == sens_seeds
    sens_rows = read_dir([ROOT + '/research/author_notes_2026-10-10/collections/sensitivity_20261011_101348_utc/results/sensitivity/completed'])
    assert len(sens_rows) == 1536

    fic_book = json.load(open(ABM + '/research/full_icu_confirmation_2026-10-08/seedbook.json'))['seed_blocks']
    hA, hB = fic_book['holdout']['A'], fic_book['holdout']['B']
    hold_seeds = hA[:128] + hA[128:] + hB[:128] + hB[128:]  # 4 strata
    hold_rows = read_dir([ABM + '/results/full_icu_confirmation_2026-10-08/holdout_cluster/completed',
                          ABM + '/results/full_icu_confirmation_2026-10-08/holdout_local/completed'])
    assert len(hold_rows) == 1536
    upper_cfg = json.load(open(ABM + '/research/rule_surface_followup_2026-10-11/configs/independent_upper.json'))
    upper_seeds = hA[:128] + hB[:128]
    assert upper_cfg['seeds'] == upper_seeds
    upper_rows = read_dir([ROOT + '/research/rule_surface_followup_2026-10-11/collections/independent_upper_20261012_0125_utc/restored/repo/results/rule_surface_followup_2026-10-11/independent_upper/completed'])
    assert len(upper_rows) == 256

    eA, eB = fic_book['environment']['A'], fic_book['environment']['B']
    env_seeds = eA + eB
    env_keep = {'independent_low075_env1', 'independent_low100_env1', 'independent_low125_env1'}
    env_rows = read_dir([ABM + '/results/full_icu_confirmation_2026-10-08/environment_cluster/completed',
                         ABM + '/results/full_icu_confirmation_2026-10-08/environment_local/completed'], keep=env_keep)
    assert len(env_rows) == 1536

    alt_src = ROOT + '/research/robustness_uq_2026-10-17/collection_2026-10-21/execution01/source/research/robustness_uq_2026-10-17'
    alt_book = json.load(open(alt_src + '/seedbook.json'))
    pol_rows = read_dir([ROOT + '/research/robustness_uq_2026-10-17/collection_2026-10-21/execution01/banks/policy/completed'])
    grid_rows = read_dir([ROOT + '/research/robustness_uq_2026-10-17/collection_2026-10-21/execution01/banks/grid/completed'])
    assert len(pol_rows) == 1536 and len(grid_rows) == 1536
    selection = json.load(open(alt_src + '/selection.json'))

    C = {}
    for sid in ['clustered_low', 'clustered_center', 'clustered_high', 'imperfect_low', 'imperfect_center', 'imperfect_high']:
        C['sens_' + sid] = cond_arrays(sens_rows, sid, sens_seeds)
    hold_ids = {'independent': 'independent_target_candidate', 'clustered': 'persistent_target_candidate',
                'imperfect': 'imperfect_target_candidate'}
    for prof, sid in hold_ids.items():
        C['hold_' + prof] = cond_arrays(hold_rows, sid, hold_seeds)
    C['ind_center256'] = cond_arrays(hold_rows, 'independent_target_candidate', upper_seeds)
    C['ind_upper'] = cond_arrays(upper_rows, 'independent_upper', upper_seeds)
    for sid in sorted(env_keep):
        C['env_' + sid] = cond_arrays(env_rows, sid, env_seeds)
    for prof in ('independent', 'clustered', 'imperfect'):
        C['alt_' + prof] = cond_arrays(pol_rows, prof, alt_book['policy'])
        for m in ('m1', 'm2'):
            C[f'grid_{prof}_{m}'] = cond_arrays(grid_rows, f'{prof}_{m}', alt_book['grid'])
    for k, c in C.items():
        if not k.startswith('grid'):
            assert c['h'] == (0.9, 0.95), (k, c['h'])
    # grid counts vs selection.json
    for prof in ('independent', 'clustered', 'imperfect'):
        gp = {g['scenario_id']: g for g in selection['profiles'][prof]['grid_points']}
        for m in ('m1', 'm2'):
            X = C[f'grid_{prof}_{m}']['X']
            assert int(X[:, 0].sum()) == gp[f'{prof}_{m}']['control_count']
            assert abs(X[:, 1].sum() - gp[f'{prof}_{m}']['control_patient_days']) < 1e-6

    # ---------------------------------------------------------------- weights
    W_sens = weights(rng, [128, 128])
    W_hold = weights(rng, [128, 128, 128, 128])
    W_c256 = np.concatenate([W_hold[:, 0:128], W_hold[:, 256:384]], axis=1)  # A-first + B-first
    W_env = weights(rng, [256, 256])
    W_alt = weights(rng, [512])
    W_grid = weights(rng, [128])
    Wof = lambda k: (W_sens if k.startswith('sens') else W_hold if k.startswith('hold') else
                     W_c256 if k in ('ind_center256', 'ind_upper') else W_env if k.startswith('env') else
                     W_alt if k.startswith('alt') else W_grid)

    pt, bt = {}, {}
    for k, c in C.items():
        B, P, F = BPF(c['X'])
        Bb, Pb, Fb = BPF(c['X'], Wof(k))
        pt[k] = dict(B=float(B), P=float(P), F=float(F), p=c['p'], n=c['n'])
        bt[k] = dict(B=Bb, P=Pb, F=Fb)
        pt[k]['B95'] = q(Bb)
        pt[k]['F95'] = q(100 * Fb)
    res['conditions'] = pt

    # ---------------------------------------------------------------- slopes
    def seg(lo, hi, src):
        out = {}
        for s, d in (('pt', pt), ('bt', bt)):
            dlnB = np.log(d[hi]['B']) - np.log(d[lo]['B'])
            dF = 100 * (d[hi]['F'] - d[lo]['F'])
            dlnp = math.log(pt[hi]['p'] / pt[lo]['p'])
            out[s] = dict(dFdlnB=dF / dlnB, dFdlnp=dF / dlnp, dlnBdlnp=dlnB / dlnp)
        return out

    def ols(keys):
        out = {}
        lnp = np.array([math.log(pt[k]['p']) for k in keys])
        for s, d in (('pt', pt), ('bt', bt)):
            lnB = np.stack([np.log(np.asarray(d[k]['B'])) for k in keys], -1)
            F = np.stack([100 * np.asarray(d[k]['F']) for k in keys], -1)

            def slope(x, y):
                xc = x - x.mean(-1, keepdims=True)
                return (xc * (y - y.mean(-1, keepdims=True))).sum(-1) / (xc * xc).sum(-1)
            lnp_b = np.broadcast_to(lnp, lnB.shape)
            out[s] = dict(dFdlnB=slope(lnB, F), dFdlnp=slope(lnp_b, F), dlnBdlnp=slope(lnp_b, lnB))
        return out

    S = {
        ('independent', 'primary'): seg('ind_center256', 'ind_upper', ''),
        ('independent', 'env_seg1'): seg('env_independent_low075_env1', 'env_independent_low100_env1', ''),
        ('independent', 'env_seg2'): seg('env_independent_low100_env1', 'env_independent_low125_env1', ''),
        ('independent', 'env_ols'): ols(['env_independent_low075_env1', 'env_independent_low100_env1', 'env_independent_low125_env1']),
        ('clustered', 'primary'): seg('sens_clustered_center', 'sens_clustered_high', ''),
        ('clustered', 'low_center'): seg('sens_clustered_low', 'sens_clustered_center', ''),
        ('clustered', 'ols'): ols(['sens_clustered_low', 'sens_clustered_center', 'sens_clustered_high']),
        ('imperfect', 'primary'): seg('sens_imperfect_center', 'sens_imperfect_high', ''),
        ('imperfect', 'low_center'): seg('sens_imperfect_low', 'sens_imperfect_center', ''),
        ('imperfect', 'ols'): ols(['sens_imperfect_low', 'sens_imperfect_center', 'sens_imperfect_high']),
    }
    res['slopes'] = {}
    for (prof, est), v in S.items():
        row = {}
        for m in ('dFdlnB', 'dFdlnp', 'dlnBdlnp'):
            row[m] = float(v['pt'][m])
            row[m + '95'] = q(v['bt'][m])
        row['per10'] = row['dFdlnB'] * math.log(1.1)
        row['per10_95'] = [x * math.log(1.1) for x in row['dFdlnB95']]
        res['slopes'][f'{prof}/{est}'] = row

    # consistency with common-baseline F*
    centers = {'independent': 'ind_center256', 'clustered': 'sens_clustered_center', 'imperfect': 'sens_imperfect_center'}
    cb = {'independent': 66.69796229001121, 'clustered': 54.99471289475639, 'imperfect': 30.03674227677826}
    res['consistency'] = {}
    for prof, ck in centers.items():
        val = 100 * pt[ck]['F'] + S[(prof, 'primary')]['pt']['dFdlnB'] * math.log(TARGET / pt[ck]['B'])
        res['consistency'][prof] = dict(linear=val, common_baseline=cb[prof], diff=val - cb[prof])

    # ---------------------------------------------------------------- bias
    obs = {('original', p): 'hold_' + p for p in ('independent', 'clustered', 'imperfect')}
    obs.update({('alternative', p): 'alt_' + p for p in ('independent', 'clustered', 'imperfect')})
    res['bias'] = {}
    bias_draws = {}
    for (wf, prof), ok in obs.items():
        ck = centers[prof]
        bhi_pt = S[(prof, 'primary')]['pt']['dFdlnB']
        bhi_bt = S[(prof, 'primary')]['bt']['dFdlnB']
        Bo_pt, Bo_bt = pt[ok]['B'], bt[ok]['B']
        Bc_pt, Bc_bt = pt[ck]['B'], bt[ck]['B']
        prim_pt = bhi_pt * math.log(Bo_pt / TARGET)
        prim_bt = bhi_bt * np.log(Bo_bt / TARGET)
        if prof == 'independent':
            pw_pt, pw_bt = prim_pt, prim_bt
        else:
            blo_pt = S[(prof, 'low_center')]['pt']['dFdlnB']
            blo_bt = S[(prof, 'low_center')]['bt']['dFdlnB']
            pw_pt = (blo_pt * math.log(Bo_pt / Bc_pt) + bhi_pt * math.log(Bc_pt / TARGET)) if Bo_pt < Bc_pt else prim_pt
            below = Bo_bt < Bc_bt
            pw_bt = np.where(below, blo_bt * np.log(Bo_bt / Bc_bt) + bhi_bt * np.log(Bc_bt / TARGET), prim_bt)
        beta_only = sorted(x * math.log(Bo_pt / TARGET) for x in q(bhi_bt))
        allest = [float(S[(prof, e)]['pt']['dFdlnB']) * math.log(Bo_pt / TARGET)
                  for (pp, e) in S if pp == prof]
        Fo_pt, Fo_bt = 100 * pt[ok]['F'], 100 * bt[ok]['F']
        res['bias'][f'{wf}/{prof}'] = dict(
            B_obs=Bo_pt, pct_vs_3=100 * (Bo_pt / TARGET - 1), B_center=Bc_pt,
            piecewise=pw_pt, piecewise95=q(pw_bt), primary=prim_pt, primary95=q(prim_bt),
            beta_only95=beta_only, range=[min(allest), max(allest)],
            F_obs=Fo_pt, F_at3_pw=Fo_pt - pw_pt, F_at3_pw95=q(Fo_bt - pw_bt),
            frac_draws_below_center=float(np.mean(Bo_bt < Bc_bt)))
        bias_draws[(wf, prof)] = dict(pw=pw_bt, prim=prim_bt, F=Fo_bt, pw_pt=pw_pt, prim_pt=prim_pt, F_pt=Fo_pt)

    lv = 100 * (1 - 0.05 / 3)
    lo3, hi3 = (100 - lv) / 2, 100 - (100 - lv) / 2
    res['contrast_bias'] = {}
    for wf in ('original', 'alternative'):
        for a, b in (('independent', 'clustered'), ('independent', 'imperfect'), ('clustered', 'imperfect')):
            x, y = bias_draws[(wf, a)], bias_draws[(wf, b)]
            obs_pt = x['F_pt'] - y['F_pt']
            obs_bt = x['F'] - y['F']
            cb_pw_pt, cb_pw_bt = x['pw_pt'] - y['pw_pt'], x['pw'] - y['pw']
            cb_pr_pt, cb_pr_bt = x['prim_pt'] - y['prim_pt'], x['prim'] - y['prim']
            res['contrast_bias'][f'{wf}/{a}-{b}'] = dict(
                observed=obs_pt, observed98=q(obs_bt, lo3, hi3),
                bias_pw=cb_pw_pt, bias_pw95=q(cb_pw_bt), bias_prim=cb_pr_pt, bias_prim95=q(cb_pr_bt),
                adjusted_pw=obs_pt - cb_pw_pt, adjusted_pw98=q(obs_bt - cb_pw_bt, lo3, hi3))

    # ---------------------------------------------------------------- routes
    res['routes'] = {}
    for wf, pre in (('original', 'hold_'), ('alternative', 'alt_')):
        for prof in ('independent', 'clustered', 'imperfect'):
            c = C[pre + prof]
            W = Wof(pre + prof)
            L0, L1 = c['L0'], c['L1']
            pd0, pd1 = c['X'][:, 1].sum(), c['X'][:, 3].sum()
            tot0, tot1 = c['X'][:, 0].sum(), c['X'][:, 2].sum()
            assert tot0 == L0.sum() and tot1 == L1.sum()
            idx = {l: i for i, l in enumerate(LABELS)}
            groups = {
                'patient': ['patient'], 'own0': ['own0'], 'ownpos': ['ownpos'], 'other': ['other'],
                'carry_all': ['carry0', 'carrypos'], 'carry0': ['carry0'], 'carrypos': ['carrypos'],
                'all': LABELS, 'all_except_carry': ['patient', 'own0', 'ownpos', 'other'],
                'manuscript_metric': ['own0', 'carry0'], 'own_any': ['own0', 'ownpos']}
            out = dict(pd_rel_diff=float(pd1 / pd0 - 1))
            T0 = W @ L0
            T1 = W @ L1
            for g, ls in groups.items():
                ii = [idx[l] for l in ls]
                n0, n1 = float(L0[:, ii].sum()), float(L1[:, ii].sum())
                r = dict(ctrl=int(n0), pol=int(n1), share_ctrl=100 * n0 / tot0, share_pol=100 * n1 / tot1,
                         contribution=100 * (n0 / pd0 - n1 / pd1) / (tot0 / pd0))
                if n0 > 0:
                    r['reduction'] = 100 * (1 - n1 / n0)
                    b0, b1 = T0[:, ii].sum(1), T1[:, ii].sum(1)
                    with np.errstate(divide='ignore', invalid='ignore'):
                        red = 100 * (1 - b1 / b0)
                    good = np.isfinite(red)
                    r['reduction95'] = q(red[good])
                    r['undefined_draws'] = int((~good).sum())
                    r['rate_reduction'] = 100 * (1 - (n1 / pd1) / (n0 / pd0))
                r['dagger'] = n0 < 30
                out[g] = r
            out['contrib_sum'] = sum(out[g]['contribution'] for g in ('patient', 'own0', 'ownpos', 'other', 'carry_all'))
            out['F'] = 100 * pt[pre + prof]['F']
            res['routes'][f'{wf}/{prof}'] = out

    # ---------------------------------------------------------------- selection arithmetic
    # archived MC half-widths
    hr = json.load(open(ABM + '/research/full_icu_confirmation_2026-10-08/holdout_report/report.json'))
    mc = {}
    for s in hr['scenarios']:
        sid = s['scenario']['id'] if isinstance(s['scenario'], dict) else s['scenario']
        prof = {v: k for k, v in hold_ids.items()}[sid]
        iv = s['rate_prevented_fraction_uncertainty']['interval95']
        mc[('original', prof)] = dict(F=100 * s['rate_prevented_fraction'], iv=[100 * x for x in iv],
                                      half=50 * (iv[1] - iv[0]))
    an = json.load(open(ROOT + '/research/robustness_uq_2026-10-17/collection_2026-10-21/execution01/outputs/analysis.json'))
    for s in an['scenarios']:
        prof = s['scenario']['id'] if isinstance(s['scenario'], dict) else s['scenario']
        iv = s['F95']['interval']
        mc[('alternative', prof)] = dict(F=100 * s['F'], iv=[100 * x for x in iv], half=50 * (iv[1] - iv[0]))
    su = json.load(open(ABM + '/research/editorial_polish_2026-10-11/source_revision_02/input/evidence/selection_uncertainty.json'))
    pooled64 = [c for c in su['cohorts'] if c['cohort'] == 'pooled64'][0]['profiles']
    pint = {
        ('original', 'independent'): (0.14327558866543832, [0.14004, 0.14792], 'S8 pooled (rounded, supplement)'),
        ('original', 'clustered'): (pooled64['persistent']['point']['p'], pooled64['persistent']['bootstrap']['conditional_valid_draw_interval95'], 'S14 pooled64 (selection_uncertainty.json)'),
        ('original', 'imperfect'): (pooled64['imperfect']['point']['p'], pooled64['imperfect']['bootstrap']['conditional_valid_draw_interval95'], 'S14 pooled64 (selection_uncertainty.json)'),
    }
    for prof in ('independent', 'clustered', 'imperfect'):
        sp = selection['profiles'][prof]
        pint[('alternative', prof)] = (sp['selected_p'], sp['selection_uncertainty']['conditional_on_success95'], 'selection.json')
    res['selection'] = {}
    for (wf, prof), (p0, (pl, ph), src) in pint.items():
        sl = S[(prof, 'primary')]
        b_pt = float(sl['pt']['dFdlnp'])
        b_bt = sl['bt']['dFdlnp']
        b_lo, b_hi = q(b_bt)
        a_lo, a_hi = math.log(pl / p0), math.log(ph / p0)
        lo_draws, hi_draws = b_bt * a_lo, b_bt * a_hi
        env_min = np.minimum(lo_draws, hi_draws)
        env_max = np.maximum(lo_draws, hi_draws)
        Smax = max(abs(b * a) for b in (b_pt, b_lo, b_hi) for a in (a_lo, a_hi))
        half = mc[(wf, prof)]['half']
        r = dict(p=p0, p_int=[pl, ph], source=src, pct=[100 * (pl / p0 - 1), 100 * (ph / p0 - 1)],
                 slope=b_pt, slope95=[b_lo, b_hi],
                 shift_point=[b_pt * a_lo, b_pt * a_hi],
                 shift_slopeunc_endpoints=[min(b_lo * a_lo, b_hi * a_lo, b_lo * a_hi, b_hi * a_hi),
                                           max(b_lo * a_lo, b_hi * a_lo, b_lo * a_hi, b_hi * a_hi)],
                 shift_slopeunc_lowend_highend=[float(np.percentile(lo_draws, 2.5)), float(np.percentile(hi_draws, 97.5))],
                 shift_slopeunc_envelope=[float(np.percentile(env_min, 2.5)), float(np.percentile(env_max, 97.5))],
                 S=Smax, mc_half=half, mc_iv=mc[(wf, prof)]['iv'],
                 mc_half_rounded_iv=(round(mc[(wf, prof)]['iv'][1], 1) - round(mc[(wf, prof)]['iv'][0], 1)) / 2,
                 combined=math.sqrt(half ** 2 + Smax ** 2))
        r['combined_rounded_mc'] = math.sqrt(r['mc_half_rounded_iv'] ** 2 + Smax ** 2)
        if wf == 'alternative':
            g1, g2 = f'grid_{prof}_m1', f'grid_{prof}_m2'
            gs_pt = math.log(pt[g2]['B'] / pt[g1]['B']) / math.log(2)
            gs_bt = np.log(bt[g2]['B'] / bt[g1]['B']) / math.log(2)
            bB_pt = float(sl['pt']['dFdlnB'])
            bB_bt = sl['bt']['dFdlnB']
            prod_pt = bB_pt * gs_pt
            prod_bt = bB_bt * gs_bt
            r['alt_grid'] = dict(B_m1=pt[g1]['B'], B_m2=pt[g2]['B'], dlnBdlnp=gs_pt, dlnBdlnp95=q(gs_bt),
                                 orig_dlnBdlnp=float(sl['pt']['dlnBdlnp']), dFdlnp=prod_pt, dFdlnp95=q(prod_bt),
                                 shift_point=[prod_pt * a_lo, prod_pt * a_hi])
        res['selection'][f'{wf}/{prof}'] = r

    res['meta'] = dict(draws=D, rng_seed=RNG_SEED, python=sys.version, numpy=np.__version__,
                       n_records=len(sens_rows) + len(hold_rows) + len(upper_rows) + len(env_rows) + len(pol_rows) + len(grid_rows),
                       note='Environment records filtered to env1 at read time; env_rows count is env1 only.')
    json.dump(res, open(OUT, 'w'), indent=1, default=float)
    print('wrote', OUT)


if __name__ == '__main__':
    main()
