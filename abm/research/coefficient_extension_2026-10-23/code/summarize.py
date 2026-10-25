#!/usr/bin/env python3
"""Protocol outputs of coefficient_extension_2026-10-23 (summary JSON).

Modelled on research/workflow_coefficient_2026-10-21/code/summarize.py (not edited).

Inputs (read-only):
  <records>/<seed>.json + manifest.json, each validated by run_batch.validate_record
  (payload digest, manifest identity, seed, horizon, required pairs) and its gate
  recomputed here, independently, against the payload-verified archived
  extended record of the same seed. A seed whose stored or recomputed gate
  fails is excluded from every estimate and reported.

Statistics (PROTOCOL.md "Inference"):
  * F* = 1 - mean(K1)/mean(K0) over seeds (ratio of means, no seed dropped).
  * Delta-method normal intervals as in the earlier studies (frozen
    run.mean_uncertainty; SE of F* = sd((K1_i - r K0_i)/mean K0)/sqrt(n)); a
    difference D = F_a - F_b on the same seeds uses
    psi_i = -(y_a - r_a x_a)/mean x_a + (y_b - r_b x_b)/mean x_b.
  * Whole-seed bootstrap: 20,000 draws, one index matrix
    numpy.random.default_rng(20261023).integers(0, n, (20000, n)) per seed set
    (seeds in seedbook order), shared by every statistic on that set; percentile
    95% intervals; undefined draws counted and excluded.
  * Seed sets: pairs carried by every record (archived, crossing, completion,
    transition) use all passing seeds; fine-grid pairs use the passing seeds
    that carry the fine grid (all of them under the full plan, bank A under
    the fallback plan).

Predictions (PROTOCOL.md "Predictions written before any outcome"), reported
whatever the outcome, with the protocol's criteria verbatim:
  1. crossing e* in [.91, .935];
  2. completion-range signs;
  3. transition-defined |D| < 5 pp at e = 1;
  4. max over the 132 fine-grid rows of |F_geo - F*| <= 3.5 pp.
Secondary: the exact ICU projection of the measurement study (projection.py),
run only after its code check reproduces the round-1 values.

summary.json is written only when every planned record (plan from
../BENCHMARK.json) exists, passes the gate (stored and recomputed) and every
self-check passes; otherwise the same content goes to progress_summary.json.

Usage (from the study folder):
  python -B code/summarize.py
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
from adapter import wb  # noqa: E402
import closed_forms  # noqa: E402
import gate  # noqa: E402
import pairs  # noqa: E402
import projection  # noqa: E402
import run_batch  # noqa: E402

STUDY = CODE.parent
BOOTSTRAP_DRAWS = 20000
BOOTSTRAP_SEED = 20261023
H_CROSS = (.9, .95)
# PROTOCOL.md predictions (verbatim criteria).
P1_INTERVAL = (.91, .935)
P1_ISOLATED_REFERENCE = .9235
P2_E08_EXPECTED_SIGN = {(.6, .7): -1, (.7, .8): -1, (.8, .9): +1}
P2_ISOLATED_PP_PROTOCOL = {(.6, .7): -3.5, (.7, .8): -2.4, (.8, .9): 1.2}
P2_RANGES = ((.6, .7), (.7, .8), (.8, .9), (.9, .95), (.95, .99))
P3_TOLERANCE = .05
P3_PER_BOUNDARY_PP_PROTOCOL = 21.5
P4_TOLERANCE = .035
P4_ARCHIVED_REFERENCE_PP = 2.78
ROUND1_ARCHIVED_GEO_MAX_PP = 2.7802003649699536


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


canonical = gate.canonical


# ----------------------------------------------------------------- loading
def seedbook_order():
    book = json.loads((wb.BRIDGE / 'seedbook.json').read_text())
    return list(book['all_seeds']), wb.seed_banks()


def load_records(records):
    """Validate every record, recompute its gate, keep only what the analysis needs."""
    manifest_path = records / 'manifest.json'
    if not manifest_path.exists():
        raise SystemExit(f'{manifest_path} missing: no run in {records}')
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('study') != run_batch.STUDY_NAME:
        raise SystemExit('Records belong to another study')
    data, problems = {}, []
    legacy_ids = set(gate.PRIMARY_PAIR_IDS)
    for path in sorted(records.glob('[0-9]*.json')):
        seed = int(path.stem)
        try:
            kind = json.loads(path.read_text())['surface']
            row = run_batch.validate_record(path, manifest, seed, kind)
        except Exception as error:  # reported; the seed counts as missing
            problems.append(dict(seed=seed, problem='invalid record (not used)', error=repr(error)))
            continue
        extended = gate.load_extended(seed)
        recomputed = gate.run_gate(row, extended['row'])
        stored = {k: v for k, v in row['gate'].items() if k != 'archived_extended_record'}
        recheck = dict(passed=recomputed['passed'], equals_stored=canonical(recomputed) == canonical(stored),
                       archived_file_sha256_unchanged=row['gate']['archived_extended_record']['file_sha256'] == extended['file_sha256'])
        if not (recheck['equals_stored'] and recheck['archived_file_sha256_unchanged']):
            problems.append(dict(seed=seed, problem='stored gate disagrees with recomputation', recheck=recheck))
        K = {pid: (entry['baseline']['K'], entry['policy']['K']) for pid, entry in row['pairs'].items()}
        archived_K = {pid: (extended['row']['pairs'][pid]['baseline']['K'], extended['row']['pairs'][pid]['policy']['K'])
                      for pid in pairs.ARCHIVED_IDS}
        data[seed] = dict(
            seed=seed, bank=row['bank'], kind=kind, K=K, archived_K=archived_K,
            gate_ok=bool(row['gate']['passed'] and recomputed['passed'] and recheck['equals_stored']
                         and recheck['archived_file_sha256_unchanged']),
            gate=gate.gate_summary(row['gate']), gate_recheck=recheck,
            profiles={k: row['profiles'][k] for k in legacy_ids},
            observation_patient_days=row['observation_patient_days'], event_counts=row['event_counts'],
            cpu_seconds=row['cpu_seconds'], process_cpu_seconds_total=row.get('process_cpu_seconds_total'),
            wall_seconds=row['wall_seconds'], record_sha256=sha256(path))
        del row, extended
    failed = sorted(p.name for p in records.glob('failed_*.json') if re.fullmatch(r'failed_\d+\.json', p.name))
    current = run_batch.build_manifest()
    binding = dict(manifest_identity_sha256=manifest['identity_sha256'], manifest_sha256=sha256(manifest_path),
                   current_code_identity_sha256=current['identity_sha256'],
                   current_code_matches_manifest=current['identity_sha256'] == manifest['identity_sha256'])
    return data, manifest, failed, problems, binding


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
    return dict(level=level, interval=[float(low), float(high)], defined_draws=int(defined.size), undefined_draws=undefined)


def boot_response(boot, x, y):
    mx, my = boot.means(x), boot.means(y)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(mx == 0, np.nan, 1 - my / mx)


class Surface:
    """F* of every pair on its seed set, with bootstrap draws cached."""

    def __init__(self, data, all_seeds, fine_seeds):
        self.data = data
        self.sets = {'all': list(all_seeds), 'fine': list(fine_seeds)}
        self.boots = {'all': Bootstrap(all_seeds) if all_seeds else None}
        self.boots['fine'] = self.boots['all'] if fine_seeds == all_seeds else (Bootstrap(fine_seeds) if fine_seeds else None)
        self.common = set(p['id'] for p in pairs.KINDS['no_fine_grid'])
        self._draws = {}
        self._blocks = {}

    def set_name(self, pid):
        if pid in self.common or self.sets['fine'] == self.sets['all']:
            return 'all'
        return 'fine'

    def seeds(self, pid):
        return self.sets[self.set_name(pid)]

    def vectors(self, pid, seeds=None):
        seeds = self.seeds(pid) if seeds is None else seeds
        return (np.array([self.data[s]['K'][pid][0] for s in seeds], float),
                np.array([self.data[s]['K'][pid][1] for s in seeds], float))

    def draws(self, pid):
        if pid not in self._draws:
            x, y = self.vectors(pid)
            self._draws[pid] = boot_response(self.boots[self.set_name(pid)], x, y)
        return self._draws[pid]

    def block(self, pid):
        if pid not in self._blocks:
            x, y = self.vectors(pid)
            b = ratio_block(x.tolist(), y.tolist())
            b['F_bootstrap_95'] = percentile_interval(self.draws(pid))
            b['seed_set'] = self.set_name(pid)
            self._blocks[pid] = b
        return self._blocks[pid]

    def F(self, pid):
        return self.block(pid)['limiting_response']

    def difference(self, pid_a, pid_b):
        """D = F*(a) - F*(b) on the common seed set (both pairs must share it)."""
        if self.set_name(pid_a) != self.set_name(pid_b):
            raise ValueError('Difference across different seed sets')
        xa, ya = self.vectors(pid_a)
        xb, yb = self.vectors(pid_b)
        fa, fb = self.F(pid_a), self.F(pid_b)
        n = len(xa)
        ra, rb = self.block(pid_a)['K1_over_K0'], self.block(pid_b)['K1_over_K0']
        psi = -(ya - ra * xa) / statistics.fmean(xa.tolist()) + (yb - rb * xb) / statistics.fmean(xb.tolist())
        se = statistics.stdev(psi.tolist()) / math.sqrt(n) if n > 1 else None
        d = fa - fb
        draws = self.draws(pid_a) - self.draws(pid_b)
        finite = draws[np.isfinite(draws)]
        return dict(pair_a=pid_a, pair_b=pid_b, F_a=fa, F_b=fb, D=d, n=n, seed_set=self.set_name(pid_a),
                    delta_monte_carlo_se=se, normal_95_interval=None if se is None else [d - 1.96 * se, d + 1.96 * se],
                    bootstrap_95=percentile_interval(draws),
                    fraction_draws_positive=float(np.mean(finite > 0)) if finite.size else None,
                    fraction_draws_negative=float(np.mean(finite < 0)) if finite.size else None,
                    _draws=draws)


def strip_draws(value):
    if isinstance(value, dict):
        return {k: strip_draws(v) for k, v in value.items() if not k.startswith('_')}
    if isinstance(value, list):
        return [strip_draws(v) for v in value]
    return value


def pair_lookup(h0, h1, rho, e, prefer=('archived', 'crossing', 'completion', 'fine_grid')):
    """Id of the computed pair with these parameters (duplicates are bitwise equal; prefer the archived one)."""
    matches = [p['id'] for p in pairs.ALL if (p['h0'], p['h1'], p['rho'], p['e']) == (h0, h1, rho, e)]
    for group in prefer:
        for pid in matches:
            if pairs.PAIR_INFO[pid]['group'] == group:
                return pid
    return None


# --------------------------------------------------------------- predictions
def interpolated_root(e_values, d_values):
    """All adjacent sign changes of D (sign: D > 0 vs D <= 0) with linear-interpolation roots; first = lowest e."""
    roots = []
    for k in range(len(e_values) - 1):
        lo, hi = d_values[k], d_values[k + 1]
        if (lo > 0) != (hi > 0):
            roots.append(dict(e_low=e_values[k], e_high=e_values[k + 1], root=e_values[k] + (e_values[k + 1] - e_values[k]) * lo / (lo - hi)))
    return roots


def bootstrap_roots(e_values, draws):
    """Per-draw first sign change (same rule as the point estimate); counts of 0 and >1 sign changes."""
    D = np.vstack(draws)  # (n_e, draws)
    e = np.asarray(e_values, float)
    finite = np.all(np.isfinite(D), axis=0)
    pos = D > 0
    change = pos[:-1] != pos[1:]
    count = change.sum(axis=0)
    first = np.argmax(change, axis=0)
    cols = np.arange(D.shape[1])
    lo, hi = D[first, cols], D[first + 1, cols]
    with np.errstate(divide='ignore', invalid='ignore'):
        root = e[first] + (e[first + 1] - e[first]) * lo / (lo - hi)
    root = np.where((count >= 1) & finite, root, np.nan)
    return root, dict(draws=int(D.shape[1]), draws_nonfinite=int((~finite).sum()),
                      draws_no_sign_change=int(((count == 0) & finite).sum()),
                      draws_one_sign_change=int(((count == 1) & finite).sum()),
                      draws_several_sign_changes=int(((count > 1) & finite).sum()))


def prediction_1(surface):
    e_all = sorted({p['e'] for p in pairs.ALL if (p['h0'], p['h1']) == H_CROSS and p['rho'] in (0., .5)})
    grid, sources = [], {}
    for e in e_all:
        p0, p5 = pair_lookup(*H_CROSS, 0., e), pair_lookup(*H_CROSS, .5, e)
        if p0 and p5 and surface.set_name(p0) == 'all' and surface.set_name(p5) == 'all':
            grid.append(e)
            sources[e] = (p5, p0)
    diffs = {e: surface.difference(*sources[e]) for e in grid}
    d = [diffs[e]['D'] for e in grid]
    roots = interpolated_root(grid, d)
    estar = roots[0]['root'] if roots else None
    boot_root, counts = bootstrap_roots(grid, [diffs[e]['_draws'] for e in grid])
    iso = [closed_forms.Fk(2, *H_CROSS, .5, e) - closed_forms.Fk(2, *H_CROSS, 0., e) for e in grid]
    iso_roots = interpolated_root(grid, iso)
    return dict(
        statement=('Crossing. The ICU crossing e* is the root of D(e) = F*(rho=.5) - F*(rho=0) at h .90->.95, found by '
                   'linear interpolation between adjacent grid e values. It lies in [.91, .935]. (Isolated k=2 root: .9235.)'),
        criterion='holds iff the interpolated root e* (lowest-e sign change of D on the grid) lies in [.91, .935]',
        e_grid=grid, n_seeds=len(surface.sets['all']),
        grid_note=('Every e at which both rho=0 and rho=.5 rows at h .90->.95 are computed on all inference seeds '
                   '(archived, crossing and, under the full plan, fine-grid rows; duplicates are bitwise equal and the '
                   'archived/crossing row is used).'),
        D_by_e={f'{e:g}': diffs[e] for e in grid},
        sign_changes=roots, n_sign_changes=len(roots), e_star=estar,
        e_star_bootstrap_95=percentile_interval(boot_root), bootstrap_sign_change_counts=counts,
        isolated_k2_exact_root=closed_forms.estar(*H_CROSS),
        isolated_k2_linear_interpolation_on_same_grid=iso_roots[0]['root'] if iso_roots else None,
        holds=estar is not None and P1_INTERVAL[0] <= estar <= P1_INTERVAL[1],
        unique_sign_change=len(roots) == 1)


def prediction_2(surface):
    rows = {}
    for h0, h1 in P2_RANGES:
        for e in (.8, 1.):
            p5, p0 = pair_lookup(h0, h1, .5, e), pair_lookup(h0, h1, 0., e)
            if p5 is None or p0 is None:
                continue
            block = surface.difference(p5, p0)
            iso = closed_forms.Fk(2, h0, h1, .5, e) - closed_forms.Fk(2, h0, h1, 0., e)
            block.update(h0=h0, h1=h1, e=e, isolated_k2_D=iso)
            rows[f'h{h0:g}_to_{h1:g}_e{e:g}'] = block
    components = {}
    for (h0, h1), sign in P2_E08_EXPECTED_SIGN.items():
        key = f'h{h0:g}_to_{h1:g}_e0.8'
        d = rows[key]['D']
        components['e0.8_' + key] = dict(expected='D < 0' if sign < 0 else 'D > 0', D=d, holds=(d < 0) if sign < 0 else (d > 0),
                                         isolated_k2_D_pp=100 * rows[key]['isolated_k2_D'],
                                         isolated_protocol_pp=P2_ISOLATED_PP_PROTOCOL[(h0, h1)],
                                         isolated_matches_protocol_1dp=round(100 * rows[key]['isolated_k2_D'], 1) == P2_ISOLATED_PP_PROTOCOL[(h0, h1)])
    for h0, h1 in P2_RANGES:
        key = f'h{h0:g}_to_{h1:g}_e1'
        components['e1_' + key] = dict(expected='D < 0', D=rows[key]['D'], holds=rows[key]['D'] < 0,
                                       isolated_k2_D_pp=100 * rows[key]['isolated_k2_D'])
    return dict(
        statement=('Completion ranges. At e = .8, D < 0 for .60->.70 and .70->.80 (isolated -3.5 and -2.4 pp), and D > 0 '
                   'for .80->.90 (isolated +1.2 pp). At e = 1, D < 0 in every range.'),
        criterion='holds iff every component holds (point estimates); bootstrap intervals and sign fractions are reported',
        every_range_at_e1=[f'{h0:g}->{h1:g}' for h0, h1 in P2_RANGES],
        D_rows=rows, components=components, holds=all(c['holds'] for c in components.values()))


def prediction_3(surface):
    out = {}
    for e in (1., .8):
        out[f'e{e:g}'] = surface.difference(f'x_trans_c0.9_to_0.95_rho0.5_e{e:g}', f'x_trans_c0.9_to_0.95_rho0_e{e:g}')
        roots = pairs.TRANSITION_ROOTS
        h = {r: (roots[f'rho{r:g}']['h0']['h_double'], roots[f'rho{r:g}']['h1']['h_double']) for r in (0., .5)}
        out[f'e{e:g}']['isolated_k2_D'] = closed_forms.Fk(2, *h[.5], .5, e) - closed_forms.Fk(2, *h[0.], 0., e)
    per_boundary = surface.difference(pair_lookup(*H_CROSS, .5, 1.), pair_lookup(*H_CROSS, 0., 1.))
    d = out['e1']['D']
    boot = out['e1']['_draws']
    finite = boot[np.isfinite(boot)]
    return dict(
        statement=('Transition-defined policy. At e = 1, |F*(rho=.5) - F*(rho=0)| < 5 pp, against 21.5 pp under the '
                   'per-boundary policy. (Isolated k=2: exactly 0.)'),
        criterion='holds iff |D_transition(e=1)| < .05 (point estimate)',
        transition_roots=pairs.TRANSITION_ROOTS, D=out,
        abs_D_e1=abs(d), abs_D_e1_bootstrap_95=percentile_interval(np.abs(boot)),
        fraction_draws_abs_D_below_5pp=float(np.mean(np.abs(finite) < P3_TOLERANCE)) if finite.size else None,
        per_boundary_D_e1=per_boundary, per_boundary_abs_D_e1_pp=100 * abs(per_boundary['D']),
        per_boundary_matches_protocol_21p5=round(100 * abs(per_boundary['D']), 1) == P3_PER_BOUNDARY_PP_PROTOCOL,
        holds=abs(d) < P3_TOLERANCE)


def prediction_4(surface):
    rows = {}
    for p in pairs.FINE_GRID:
        fstar = surface.F(p['id'])
        fgeo = closed_forms.Fgeo(p['h0'], p['h1'], p['rho'], p['e'])
        rows[p['id']] = dict(rho=p['rho'], e=p['e'], F_star=fstar, F_geo=fgeo, F_geo_minus_F_star=fgeo - fstar,
                             F_star_bootstrap_95=surface.block(p['id'])['F_bootstrap_95']['interval'])
    errors = {pid: abs(r['F_geo_minus_F_star']) for pid, r in rows.items()}
    worst = max(errors, key=errors.get)
    archived = {}
    for p in pairs.ARCHIVED:
        fstar = surface.F(p['id'])
        fgeo = closed_forms.Fgeo(p['h0'], p['h1'], p['rho'], p['e'])
        archived[p['id']] = dict(F_star=fstar, F_geo=fgeo, F_geo_minus_F_star=fgeo - fstar)
    archived_max = max(abs(r['F_geo_minus_F_star']) for r in archived.values())
    return dict(
        statement=('Renewal comparator. Over the 132-row fine grid, the largest |F_geo - F*| is at most 3.5 pp. F_geo is '
                   'the closed form defined in the manuscript (supplement S4.2); its largest deviation on the 22 archived '
                   'rows was 2.78 pp.'),
        criterion='holds iff max over the 132 fine-grid rows of |F_geo - F*| <= .035',
        n_rows=len(rows), n_seeds=len(surface.sets[surface.set_name(pairs.FINE_GRID[0]['id'])]),
        max_abs_error=errors[worst], max_abs_error_pp=100 * errors[worst], worst_row=worst,
        rows_within_0p5pp=sum(v < .005 for v in errors.values()), rows_within_1pp=sum(v < .01 for v in errors.values()),
        F_geo_exceeds_F_star_rows=sum(r['F_geo_minus_F_star'] > 0 for r in rows.values()),
        rows=rows,
        archived_22=dict(max_abs_error_pp=100 * archived_max, protocol_reference_pp=P4_ARCHIVED_REFERENCE_PP,
                         n_seeds=len(surface.sets['all']), rows=archived),
        holds=errors[worst] <= P4_TOLERANCE)


def secondary_projection(surface):
    check = projection.code_check()
    if not check['passed']:
        return dict(status='not_run', reason='projection code check failed', code_check=check)
    Z = np.array([[surface.F('x_fine_' + pairs.pair_id(*H_CROSS, r, e)) for e in pairs.FINE_E] for r in pairs.FINE_RHO])
    grid = projection.GridSurface(pairs.FINE_RHO, pairs.FINE_E, Z)
    result = projection.exact_projection(grid)
    result.update(status='computed', code_check=check,
                  fine_grid_seed_set=surface.set_name(pairs.FINE_GRID[0]['id']),
                  fine_grid_n_seeds=len(surface.seeds(pairs.FINE_GRID[0]['id'])))
    return result


# ------------------------------------------------------------- self-checks
def self_checks(surface, data, all_seeds_256, manifest):
    checks = {}
    seeds = surface.sets['all']
    # Archived pairs reproduce the archived extended F* on the same seeds (bitwise K).
    equal = all(surface.F(pid) == ratio_block([data[s]['archived_K'][pid][0] for s in seeds],
                                              [data[s]['archived_K'][pid][1] for s in seeds])['limiting_response']
                for pid in pairs.ARCHIVED_IDS) if seeds else None
    checks['archived_pairs_F_equal_extended_archive_same_seeds'] = equal
    # Full 256: equal to the workflow-coefficient study's original_all_256 and round-1 F_geo max.
    if sorted(seeds) == sorted(all_seeds_256):
        wc = projection.checked_json(projection.WC_SUMMARY, projection.WC_SUMMARY_SHA256)['surface']['pairs']
        checks['archived_pairs_F_equal_workflow_coefficient_original_all_256'] = all(
            surface.F(pid) == wc[pid]['original_all_256']['limiting_response'] for pid in pairs.ARCHIVED_IDS)
        geo = max(abs(closed_forms.Fgeo(p['h0'], p['h1'], p['rho'], p['e']) - surface.F(p['id'])) for p in pairs.ARCHIVED)
        checks['archived_F_geo_max_error_equals_round1'] = abs(100 * geo - ROUND1_ARCHIVED_GEO_MAX_PP) < 1e-9
    else:
        checks['archived_pairs_F_equal_workflow_coefficient_original_all_256'] = 'not_applicable_partial_seed_set'
        checks['archived_F_geo_max_error_equals_round1'] = 'not_applicable_partial_seed_set'
    # Frozen run.summarize agreement on the primary profiles.
    if len(seeds) >= 2:
        _, config = wb.archived_scientific_config()
        rows = [dict(profiles=data[s]['profiles'], observation_patient_days=data[s]['observation_patient_days'],
                     event_counts=data[s]['event_counts']) for s in seeds]
        frozen = wb.bridge_run.summarize(rows, config)
        checks['primary_ratio_block_equals_frozen_run_summarize'] = all(
            frozen['profiles'][legacy][k] == surface.block(pid)[k]
            for legacy, pid in gate.PRIMARY_PAIR_IDS.items()
            for k in ('K0', 'K1', 'K1_over_K0', 'limiting_response', 'limiting_response_delta_monte_carlo_se',
                      'limiting_response_normal_95_interval'))
    else:
        checks['primary_ratio_block_equals_frozen_run_summarize'] = 'not_applicable_fewer_than_2_seeds'
    # Transition roots in the manifest equal a fresh computation.
    checks['transition_roots_equal_manifest'] = canonical(manifest['science']['transition_roots']) == canonical(pairs.TRANSITION_ROOTS)
    # F_geo closed form against truncated direct sums (as round1_analysis.py).
    direct = []
    for h0, h1, rho, e in ((.9, .95, .5, .8), (.9, .95, 0, 1), (.8, .9, .5, 1), (.9, .95, .95, .5)):
        value = 1 - (sum(closed_forms.qk(2 * j, h1, rho, e) for j in range(1, 4000))
                     / sum(closed_forms.qk(2 * j, h0, rho, e) for j in range(1, 4000)))
        direct.append(abs(value - closed_forms.Fgeo(h0, h1, rho, e)) < 1e-10)
    checks['F_geo_equals_truncated_direct_sums'] = all(direct)
    passed = all(v is True or (isinstance(v, str) and v.startswith('not_applicable')) for v in checks.values())
    return passed, checks


# --------------------------------------------------------------- main
def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name('.' + path.name + '.tmp')
    with temporary.open('w') as handle:
        handle.write(json.dumps(value, indent=1, sort_keys=True, allow_nan=False) + '\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def planned_tasks(records):
    if records.resolve() != run_batch.RECORDS.resolve() or not run_batch.BENCHMARK.exists():
        return None, None
    surface = run_batch.resolve_surface('auto')
    return surface, run_batch.plan(surface)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--records', type=Path, default=run_batch.RECORDS)
    parser.add_argument('--output', type=Path, default=None)
    args = parser.parse_args(argv)
    started = time.time()
    records = args.records.resolve()
    all_seeds, banks = seedbook_order()
    data, manifest, failed, problems, binding = load_records(records)
    if not data:
        raise SystemExit('No records in ' + str(records))
    excluded = [dict(seed=s, bank=data[s]['bank'], gate=data[s]['gate'], recheck=data[s]['gate_recheck'])
                for s in all_seeds if s in data and not data[s]['gate_ok']]
    passing = [s for s in all_seeds if s in data and data[s]['gate_ok']]
    fine = [s for s in passing if data[s]['kind'] == 'all']
    if not passing:
        raise SystemExit('No record passes the gate')
    surface = Surface(data, passing, fine)

    plan_surface, tasks = planned_tasks(records)
    if tasks is not None:
        required = {seed: kind for seed, _, kind in tasks}
        missing = [seed for seed, kind in required.items()
                   if seed not in data or not set(p['id'] for p in pairs.KINDS[kind]) <= set(data[seed]['K'])]
        complete_records = not missing and not excluded and not problems
    else:
        missing, complete_records = None, False

    pair_blocks = {}
    for p in pairs.ALL:
        if p['id'] in surface.common or fine:
            block = dict(surface.block(p['id']))
            block.update(parameters={k: p[k] for k in ('h0', 'h1', 'rho', 'e')}, **pairs.PAIR_INFO[p['id']],
                         isolated_k2_response_for_comparison_only=closed_forms.Fk(2, p['h0'], p['h1'], p['rho'], p['e']),
                         F_geo=closed_forms.Fgeo(p['h0'], p['h1'], p['rho'], p['e']))
            pair_blocks[p['id']] = block
    predictions = dict(prediction_1_crossing=prediction_1(surface), prediction_2_completion_ranges=prediction_2(surface),
                       prediction_3_transition_policy=prediction_3(surface))
    if fine:
        predictions['prediction_4_renewal_comparator'] = prediction_4(surface)
        secondary = secondary_projection(surface)
    else:
        predictions['prediction_4_renewal_comparator'] = dict(status='not_evaluable', reason='no passing record carries the fine grid')
        secondary = dict(status='not_run', reason='no passing record carries the fine grid')
    checks_pass, checks = self_checks(surface, data, all_seeds, manifest)
    checks['records_valid_and_stored_gates_agree_with_recomputation'] = not problems
    checks_pass = checks_pass and not problems
    if secondary.get('status') == 'not_run' and fine:
        checks_pass = False
    complete = bool(complete_records and checks_pass and binding['current_code_matches_manifest'])
    status = 'complete' if complete else ('partial' if tasks is not None else 'partial_or_test_bank')

    summary = dict(
        schema_version=1, study=run_batch.STUDY_NAME, status=status, records_dir=str(records),
        planned_surface=plan_surface, planned_seeds=None if tasks is None else len(tasks),
        missing_planned_seeds=missing, failed_records=failed,
        n_records=len(data), n_passing=len(passing), n_fine_grid_passing=len(fine),
        seeds_by_bank={b: sum(1 for s in passing if banks[s] == b) for b in 'AB'},
        excluded_gate_failures=excluded, record_problems=problems,
        gate=dict(definition=gate.run_gate.__doc__, all_records_pass=not excluded and not problems,
                  per_seed={str(s): data[s]['gate'] for s in all_seeds if s in data}),
        binding=binding,
        estimand='F* = 1 - mean(K1)/mean(K0) over backgrounds (original workflow, p -> 0, fixed horizon).',
        bootstrap=dict(all=surface.boots['all'].describe(),
                       fine=None if surface.boots['fine'] is None else surface.boots['fine'].describe()),
        seed_sets=dict(all=dict(n=len(passing), note='archived, crossing, completion and transition pairs'),
                       fine=dict(n=len(fine), same_as_all=fine == passing, note='fine-grid pairs')),
        self_checks_pass=checks_pass, self_checks=checks,
        pairs=pair_blocks,
        predictions=predictions,
        secondary_exact_icu_projection=secondary,
        exposure=dict(observation_patient_days=mean_uncertainty([data[s]['observation_patient_days'] for s in passing])),
        compute=dict(cpu_hours_records=sum(data[s]['cpu_seconds'] for s in data) / 3600,
                     process_cpu_hours_records=sum(data[s]['process_cpu_seconds_total'] or 0 for s in data) / 3600,
                     mean_wall_seconds=statistics.fmean(data[s]['wall_seconds'] for s in data)),
        limitations=['These are fixed-horizon p -> 0 limits. They give no finite-p error bound, and the finite-p reversal '
                     'remains untested (O-3).',
                     "Transition-defined compliance assumes the stationary Markov law and the model's exit/entry opportunity rule."],
        inputs=dict(records_sha256={str(s): data[s]['record_sha256'] for s in all_seeds if s in data},
                    frozen_sha256=wb.verify_frozen(), adapter_sha256=wb.adapter_sha256(),
                    extended_manifest_sha256=wb.sha256(wb.EXTENDED_MANIFEST),
                    code_sha256={name: sha256(CODE / name) for name in
                                 ('summarize.py', 'projection.py', *run_batch.CODE_FILES)},
                    protocol_sha256=sha256(STUDY / 'PROTOCOL.md')),
        runtime=dict(python=sys.version, numpy=np.__version__),
        computed_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), compute_seconds=None)
    summary = strip_draws(summary)
    summary['compute_seconds'] = time.time() - started
    final = (STUDY / 'summary.json').resolve()
    output = args.output.resolve() if args.output is not None else (
        STUDY / ('summary.json' if complete else 'progress_summary.json')).resolve()
    if output == final and not complete:
        raise SystemExit('Refusing to write summary.json: not every planned record exists and passes the gate '
                         '(or a self-check failed). Use progress_summary.json.')
    atomic_json(output, summary)
    p = summary['predictions']
    print(json.dumps(dict(output=str(output), status=status, n_passing=len(passing), n_fine=len(fine),
                          excluded=[e['seed'] for e in excluded], self_checks_pass=checks_pass,
                          prediction_holds={k: v.get('holds') for k, v in p.items()},
                          e_star=p['prediction_1_crossing']['e_star']), indent=1))
    return 0 if checks_pass and not excluded and not problems else 1


if __name__ == '__main__':
    sys.exit(main())
