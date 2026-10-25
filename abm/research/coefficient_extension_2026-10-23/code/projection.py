"""Projection of the saved measurement-study assay limits through an ICU surface.

PROTOCOL.md "Secondary, specified now: exact ICU projection of the measurement
study". The feasible-set construction is the one of supplement S11 and of
``round1_analysis.py`` section 9, ported line for line:

  for each saved experiment at 8,000 observations per assay, with confidence
  limits (q1_lo, q1_hi, q2_lo, q2_hi) and baseline h = .9:
    e in [max((1-q1_hi)/h, 1e-9), min((1-q1_lo)/h, 1)], 81 values (empty if reversed);
    at each e, rho in [max((q2_lo-b)/v, 0), min((q2_hi-b)/v, 1)] with
      b = (1-h e)^2, v = e^2 h (1-h), 61 values (slice skipped if reversed);
    width = max - min of the surface over every evaluated (rho, e); a set
      with no nonempty slice is empty.

Three surfaces go through the same loop:
  * isolated two-opportunity bridge (the code check: must reproduce the
    round-1 values ``median_width_pp_isolated_k2_same_grid_check``);
  * the round-1 approximate surface G = F_geo + bilinear residual on the
    4 x 4 archived grid, held constant beyond it (reproduces 6.2/0.9/6.7/11.7;
    used only as the labelled comparison);
  * the exact surface: bilinear interpolation of the computed ICU F* on the
    11 x 12 fine grid. Its domain is rho in [0, .95], e in [.5, 1]. A set with
    any evaluated feasible point outside the domain is counted as out of domain
    and excluded from the widths; nothing is extrapolated or clipped.
No new measurement, no coverage claim.
"""
import sys

sys.dont_write_bytecode = True

import hashlib  # noqa: E402
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import closed_forms  # noqa: E402

CODE = Path(__file__).resolve().parent
ABM = CODE.parents[2]
UQ = ABM / 'research/robustness_uq_2026-10-17/uq'
RAW_UQ = UQ / 'run_001/raw_benchmark.npz'
RAW_UQ_SHA256 = 'b5468a0f2a6a16929cc86655cbecf8b68656570d15bdded0425b624a8191d219'
CONTINUOUS = UQ / 'continuous_k2_2026-10-17/run_001/results.json'
CONTINUOUS_SHA256 = 'c10e34b2cd49b5fd5336d589669925bab0fa5c7c54754a225d0ef985bb3eefe2'
ROUND1 = ABM / 'research/manuscript_v2_2026-10-21/reviews_round1/analysis/round1_numbers.json'
ROUND1_SHA256 = 'c4a23e0a673a23711f3e453cb3194c6e4e66aa35d685df1bcbfd0569170a7e0b'
WC_SUMMARY = ABM / 'research/workflow_coefficient_2026-10-21/summary.json'
WC_SUMMARY_SHA256 = '6cd2ff7e768f6ba375e3f5c835af0d3d81dad4de5a86ba32170871476579bad2'

TRUTHS = ('independent', 'clustered', 'imperfect', 'combined')
N_PER_ASSAY = 8000
H_BASE = 0.9
N_E, N_RHO = 81, 61
APPROXIMATE_MEDIANS_PP_PROTOCOL = {'independent': 6.2, 'clustered': .9, 'imperfect': 6.7, 'combined': 11.7}
EXACT_DOMAIN = dict(rho_max=.95, e_min=.5)
# round1 approximate grid
RHO_G = np.array([0.0, 0.25, 0.5, 0.75])
E_G = np.array([0.6, 0.8, 0.9, 1.0])


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checked_json(path, expected):
    if sha256(path) != expected:
        raise RuntimeError('Input changed: ' + str(path))
    return json.loads(Path(path).read_text())


def load_limits():
    if sha256(RAW_UQ) != RAW_UQ_SHA256:
        raise RuntimeError('Saved assay limits changed: ' + str(RAW_UQ))
    z = np.load(RAW_UQ)
    return {truth: z[f'{truth}_n{N_PER_ASSAY}_confidence_limits'] for truth in TRUTHS}


# ------------------------------------------------------------------ surfaces
def isolated_k2(rg, e):
    """round1 section 9 'iso' check, expression for expression."""
    q2a = (1 - .95 * e) ** 2 + rg * e * e * .95 * .05
    q2b = (1 - .90 * e) ** 2 + rg * e * e * .90 * .10
    return 1 - q2a / q2b


def approximate_surface(original_surface_F):
    """round1 G('original', rho, e) from the 16 archived original-workflow F* values.

    ``original_surface_F[(rho, e)]`` = F* on rho in RHO_G x e in E_G.
    """
    res = np.array([[original_surface_F[(float(r), float(e))] - closed_forms.Fgeo(.9, .95, r, e) for e in E_G]
                    for r in RHO_G])

    def interp_res(rho, e):
        """Bilinear interpolation of the residual on the 4x4 grid; constant beyond the grid."""
        rho = np.clip(rho, RHO_G[0], RHO_G[-1])
        e = np.clip(e, E_G[0], E_G[-1])
        i = np.clip(np.searchsorted(RHO_G, rho, side='right') - 1, 0, len(RHO_G) - 2)
        j = np.clip(np.searchsorted(E_G, e, side='right') - 1, 0, len(E_G) - 2)
        tr = (rho - RHO_G[i]) / (RHO_G[i + 1] - RHO_G[i])
        te = (e - E_G[j]) / (E_G[j + 1] - E_G[j])
        return ((1 - tr) * (1 - te) * res[i, j] + tr * (1 - te) * res[i + 1, j]
                + (1 - tr) * te * res[i, j + 1] + tr * te * res[i + 1, j + 1])

    def G(rho, e):
        return closed_forms.Fgeo_vec(.9, .95, rho, e) + interp_res(rho, e)
    return G, res


def wc_original_surface():
    """The 16 original-workflow F* values used by round1 (workflow-coefficient summary.json)."""
    summary = checked_json(WC_SUMMARY, WC_SUMMARY_SHA256)
    icu = summary['surface']['pairs']
    out = {}
    for r in RHO_G:
        for e in E_G:
            key = f"h0.9_to_0.95_rho{('0' if r == 0 else ('%g' % r))}_e{('1' if e == 1 else ('%g' % e))}"
            out[(float(r), float(e))] = icu[key]['original']['limiting_response']
    return out


class GridSurface:
    """Bilinear interpolation of F* on a rectilinear (rho, e) grid; refuses to extrapolate."""

    def __init__(self, rho_nodes, e_nodes, values):
        self.R = np.asarray(rho_nodes, float)
        self.E = np.asarray(e_nodes, float)
        self.Z = np.asarray(values, float)
        if self.Z.shape != (len(self.R), len(self.E)) or not np.all(np.diff(self.R) > 0) or not np.all(np.diff(self.E) > 0):
            raise ValueError('Grid must be strictly increasing with values of shape (len(rho), len(e))')
        if not np.all(np.isfinite(self.Z)):
            raise ValueError('Grid values must be finite')

    def inside(self, rho, e):
        rho = np.asarray(rho, float)
        return bool(np.all(rho >= self.R[0]) and np.all(rho <= self.R[-1]) and self.E[0] <= e <= self.E[-1])

    def __call__(self, rho, e):
        rho = np.asarray(rho, float)
        if not self.inside(rho, e):
            raise ValueError('Point outside the computed grid: no extrapolation')
        R, E, Z = self.R, self.E, self.Z
        i = np.clip(np.searchsorted(R, rho, side='right') - 1, 0, len(R) - 2)
        j = int(np.clip(np.searchsorted(E, e, side='right') - 1, 0, len(E) - 2))
        tr = (rho - R[i]) / (R[i + 1] - R[i])
        te = (e - E[j]) / (E[j + 1] - E[j])
        return ((1 - tr) * (1 - te) * Z[i, j] + tr * (1 - te) * Z[i + 1, j]
                + (1 - tr) * te * Z[i, j + 1] + tr * te * Z[i + 1, j + 1])


# ------------------------------------------------------------- feasible sets
def feasible_slices(limit, h=H_BASE):
    """(e, rho-grid) slices of one experiment's joint feasible set; None if the e interval is empty."""
    l1, u1, l2, u2 = limit
    e_lo = max((1 - u1) / h, 1e-9)
    e_hi = min((1 - l1) / h, 1.0)
    if e_lo > e_hi:
        return None
    slices = []
    for e in np.linspace(e_lo, e_hi, N_E):
        base = (1 - h * e) ** 2
        var = e * e * h * (1 - h)
        r_lo = max((l2 - base) / var, 0.0)
        r_hi = min((u2 - base) / var, 1.0)
        if r_lo <= r_hi:
            slices.append((e, np.linspace(r_lo, r_hi, N_RHO)))
    return slices


def reaches_rho_above(limit, threshold, h=H_BASE):
    """round1 'joint_sets_reaching_rho_gt_0p75' rule (unclipped q2 upper bound over the e grid)."""
    l1, u1, l2, u2 = limit
    e_lo = max((1 - u1) / h, 1e-9)
    e_hi = min((1 - l1) / h, 1.0)
    return any(((u2 - (1 - h * e) ** 2) / (e * e * h * (1 - h))) > threshold for e in np.linspace(e_lo, e_hi, N_E))


def project(limits, surface, domain=None):
    """Widths (fraction units) of ``surface`` over every experiment's feasible set.

    ``domain`` (GridSurface or None): sets with any evaluated point outside it are
    counted out of domain and excluded; never evaluated outside the domain.
    """
    widths, empty, out = [], 0, dict(total=0, rho_above_max=0, e_below_min=0)
    for limit in limits:
        slices = feasible_slices(limit)
        if not slices:
            empty += 1
            continue
        if domain is not None:
            rho_out = any(rg[-1] > domain.R[-1] or rg[0] < domain.R[0] for _, rg in slices)
            e_out = any(e < domain.E[0] or e > domain.E[-1] for e, _ in slices)
            if rho_out or e_out:
                out['total'] += 1
                out['rho_above_max'] += rho_out
                out['e_below_min'] += e_out
                continue
        values = np.concatenate([surface(rg, e) for e, rg in slices])
        widths.append(values.max() - values.min())
    return dict(widths=widths, empty=empty, out_of_domain=out, defined=len(widths), experiments=len(limits))


def describe(result):
    w = np.asarray(result['widths'], float)
    return dict(median_width_pp=100 * float(np.median(w)) if w.size else None,
                quantiles_10_50_90_pp=[100 * float(q) for q in np.percentile(w, [10, 50, 90])] if w.size else None,
                defined=result['defined'], empty=result['empty'], experiments=result['experiments'],
                out_of_domain=result['out_of_domain'])


def code_check(limits=None):
    """Reproduce round-1 isolated-bridge and approximate widths; compare with certified continuous widths."""
    limits = load_limits() if limits is None else limits
    round1 = checked_json(ROUND1, ROUND1_SHA256)['icu_projection_approximate']['results']
    continuous = {c['truth']: c for c in checked_json(CONTINUOUS, CONTINUOUS_SHA256)['cases'] if c['n_per_assay'] == N_PER_ASSAY}
    G, _ = approximate_surface(wc_original_surface())
    out = {}
    for truth in TRUTHS:
        iso = project(limits[truth], isolated_k2)
        approx = project(limits[truth], G)
        reference = round1[truth][str(N_PER_ASSAY)]
        iso_median = 100 * float(np.median(iso['widths']))
        approx_median = 100 * float(np.median(approx['widths']))
        certified = 100 * continuous[truth]['information_sets']['q1_and_q2']['continuous_width']['quantiles_10_50_90'][1]
        beyond = sum(reaches_rho_above(lim, .75) for lim in limits[truth] if feasible_slices(lim))
        out[truth] = dict(
            isolated_median_pp=iso_median,
            isolated_round1_reference_pp=reference['median_width_pp_isolated_k2_same_grid_check'],
            isolated_equals_round1=iso_median == reference['median_width_pp_isolated_k2_same_grid_check'],
            isolated_certified_continuous_median_pp=certified,
            isolated_minus_certified_pp=iso_median - certified,
            approximate_median_pp=approx_median,
            approximate_round1_reference_pp=reference['median_width_pp_q1_and_q2'],
            approximate_equals_round1=approx_median == reference['median_width_pp_q1_and_q2'],
            approximate_rounded_matches_protocol=round(approx_median, 1) == APPROXIMATE_MEDIANS_PP_PROTOCOL[truth],
            defined=approx['defined'], empty=approx['empty'],
            counts_equal_round1=(approx['defined'], approx['empty']) == (reference['defined_joint'], reference['empty_joint'])
            and iso['defined'] == reference['defined_joint'],
            joint_sets_reaching_rho_gt_0p75=beyond,
            beyond_equals_round1=beyond == reference['joint_sets_reaching_rho_gt_0p75'])
    passed = all(v['isolated_equals_round1'] and v['approximate_equals_round1'] and v['counts_equal_round1']
                 and v['approximate_rounded_matches_protocol'] and v['beyond_equals_round1']
                 and abs(v['isolated_minus_certified_pp']) < 1e-6 for v in out.values())
    return dict(passed=passed, by_truth=out,
                note=('The same feasible-set loop reproduces the round-1 isolated-bridge grid widths exactly (the '
                      'code check) and the round-1 approximate ICU widths exactly; the isolated grid widths agree '
                      'with the certified continuous widths to < 1e-6 pp.'),
                inputs={str(RAW_UQ): RAW_UQ_SHA256, str(CONTINUOUS): CONTINUOUS_SHA256,
                        str(ROUND1): ROUND1_SHA256, str(WC_SUMMARY): WC_SUMMARY_SHA256})


def exact_projection(surface, limits=None):
    """Exact ICU projection through ``surface`` (a GridSurface over the fine grid)."""
    limits = load_limits() if limits is None else limits
    G, _ = approximate_surface(wc_original_surface())
    out = {}
    for truth in TRUTHS:
        exact = project(limits[truth], surface, domain=surface)
        block = describe(exact)
        # Like-for-like extra: the approximate surface on the same in-domain sets.
        same = []
        for limit in limits[truth]:
            slices = feasible_slices(limit)
            if not slices:
                continue
            if any(rg[-1] > surface.R[-1] or rg[0] < surface.R[0] for _, rg in slices) or \
                    any(e < surface.E[0] or e > surface.E[-1] for e, _ in slices):
                continue
            values = np.concatenate([G(rg, e) for e, rg in slices])
            same.append(values.max() - values.min())
        block.update(
            approximate_median_pp_comparison=APPROXIMATE_MEDIANS_PP_PROTOCOL[truth],
            exact_minus_approximate_pp=None if block['median_width_pp'] is None else block['median_width_pp'] - APPROXIMATE_MEDIANS_PP_PROTOCOL[truth],
            approximate_surface_median_pp_on_same_in_domain_sets=100 * float(np.median(same)) if same else None)
        out[truth] = block
    return dict(
        label='Projection of saved measurements (comparison with the approximate values is labelled as a comparison); no new measurement or coverage claim.',
        method=('Feasible sets as in supplement S11 / round1 section 9 (81 e x 61 rho, n = 8,000 per assay); surface = '
                'bilinear interpolation of the computed ICU F* on the fine grid; domain rho <= .95, e >= .5; sets with any '
                'evaluated feasible point outside the domain are counted (out_of_domain) and excluded, never extrapolated.'),
        domain=dict(rho=[float(surface.R[0]), float(surface.R[-1])], e=[float(surface.E[0]), float(surface.E[-1])]),
        by_truth=out)


__all__ = ['APPROXIMATE_MEDIANS_PP_PROTOCOL', 'GridSurface', 'TRUTHS', 'approximate_surface', 'code_check',
           'describe', 'exact_projection', 'feasible_slices', 'isolated_k2', 'load_limits', 'project']
