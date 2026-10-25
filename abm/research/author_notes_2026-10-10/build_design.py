#!/usr/bin/env python3
"""Construct fixed positive-p follow-ups with explicit audited seed lists."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n'
    if path.exists() and path.read_text() != text:
        raise ValueError('Refusing to change declared design: ' + str(path))
    path.write_text(text)


def collect_seeds(value, found):
    if isinstance(value, dict):
        for key, child in value.items():
            if key in ('seed', 'bootstrap_seed', 'base_seed') and type(child) is int:
                found.add(child)
            elif key == 'seeds' and isinstance(child, list):
                found.update(x for x in child if type(x) is int)
            collect_seeds(child, found)
    elif isinstance(value, list):
        for child in value:
            collect_seeds(child, found)


def main():
    found = set()
    scanned = 0
    for path in (ABM / 'research').rglob('*.json'):
        if HERE in path.parents or path.stat().st_size > 8_000_000:
            continue
        try:
            value = json.loads(path.read_text())
        except (ValueError, UnicodeError):
            continue
        collect_seeds(value, found)
        scanned += 1
    old_count = len(found)

    def seeds(module, n):
        values = []
        counter = 0
        while len(values) < n:
            value = int.from_bytes(hashlib.sha256(
                f'author_notes_2026-10-10/{module}/{counter}'.encode()).digest()[:4], 'big')
            counter += 1
            if value not in found:
                values.append(value)
                found.add(value)
        return values

    clock_seeds = seeds('clock', 128)
    sensitivity_seeds = seeds('sensitivity', 256)
    profiles = [('independent', 0., 1.), ('clustered', .5, 1.), ('imperfect', 0., .8)]
    clock = [dict(id=f'{name}_p{int(p*100):03d}', p=p, h0=.9, h1=.9,
                  rho=rho, e=e, opportunities='separate', environment_multiplier=1.)
             for name, rho, e in profiles for p in (.15, .30, .60, 1.)]
    sensitivity = [dict(id=f'{name}_{label}', p=p*factor, h0=.9, h1=.95,
                        rho=rho, e=e, opportunities='zone_separate', environment_multiplier=1.)
                   for name, rho, e, p in [('clustered', .5, 1., .07602158239372622),
                                          ('imperfect', 0., .8, .07483995991409355)]
                   for label, factor in [('low', .75), ('center', 1.), ('high', 1.25)]]
    common = dict(warmup_days=90., observation_days=56., block_days=1.,
                  model_overrides=dict(prob_patient_sick=.03, initial_environment_contaminated=False))
    for name, scenarios, selected in [('clock', clock, clock_seeds),
                                       ('sensitivity', sensitivity, sensitivity_seeds)]:
        write(HERE / 'configs' / f'{name}.json', dict(common, scenarios=scenarios, seeds=selected,
             description='New fixed author-note follow-up; no refitting or budget extension.'))
    write(HERE / 'seedbook.json', dict(
        schema_version=1, namespace=HERE.name,
        clock=dict(A=clock_seeds[:64], B=clock_seeds[64:]),
        sensitivity=dict(A=sensitivity_seeds[:128], B=sensitivity_seeds[128:]),
        audit=dict(accessible_public_json_files=scanned, unique_explicit_old_seeds=old_count,
                   excluded_json_larger_than_bytes=8_000_000,
                   external_absent_archives_certified=False,
                   within_module_shared_seeds=True, cross_module_disjoint=True)))
    write(HERE / 'output_scopes.json', dict(banks={
        name: f'results/{HERE.name}/{name}' for name in ('clock', 'sensitivity')}))
    write(HERE / 'design.json', dict(schema_version=1, study=HERE.name,
        clock_jobs=1536, sensitivity_pairs=1536, coefficient_background_replays=256,
        coefficient_pairs=22, bootstrap_draws=50000, bootstrap_seed=20261010,
        interval_levels=dict(descriptive=.95, clock_rate_family=1-.05/12,
                             sensitivity_contrast_family=1-.05/4),
        sampled_clock_target_band=[2.7,3.3], no_selection_or_grid_extension=True,
        interpretation='New follow-up after old outcomes, not retrospective registration.'))
    print('Declared 1536 baseline jobs, 1536 paired jobs, and 256 p=0 background replays.')


if __name__ == '__main__':
    main()
