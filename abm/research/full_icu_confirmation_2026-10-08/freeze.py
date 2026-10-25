#!/usr/bin/env python3
"""Create a write-once source/configuration freeze before a declared stage."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]
OLD = ABM / 'research/full_icu_response_2026-10-08'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--phase', choices=('initial', 'holdout', 'bridge'), required=True)
    args = p.parse_args()
    paths = [ABM / 'src' / n for n in ('CPEmodel.py', 'agents.py', 'routes.py')]
    paths += [ABM / 'research/full_icu_calibration_2026-10-07/run.py']
    paths += [OLD / n for n in ('run.py', 'instrumented_model.py', 'array_worker.py', 'PROTOCOL.md')]
    paths += [HERE / n for n in ('PROTOCOL.md', 'DESIGN_REVIEW.md', 'candidate_design.json',
               'build_design.py', 'seedbook.json', 'seed_audit.json', 'output_scopes.json',
               'array_worker.py', 'array_resume.sbatch', 'test_array_binding.py',
               'test_build_design.py', 'runtime_probe.py', 'runtime_validation.json', 'freeze.py')]
    paths += sorted((HERE / 'configs').glob('calibration_*.json'))
    paths += sorted((HERE / 'configs').glob('environment_*.json'))
    if args.phase == 'holdout':
        paths += [HERE / 'freeze.json', HERE / 'calibration_selection.json',
                  HERE / 'calibration_report/report.json']
        paths += sorted((HERE / 'configs').glob('holdout_*.json'))
        if len(list((HERE / 'configs').glob('holdout_*.json'))) != 3:
            raise ValueError('All three locked holdout configs required')
    elif args.phase == 'bridge':
        paths += [HERE / 'freeze.json']
        paths += [HERE / 'bridge' / n for n in ('PROTOCOL.md', 'coefficient.py',
                  'background.py', 'run.py', 'array_worker.py', 'test_bridge.py',
                  'INDEPENDENT_REVIEW.md', 'seedbook.json')]
        paths += sorted((HERE / 'bridge/configs').glob('*.json'))
    for path in paths:
        if not path.is_file() or path.is_symlink():
            raise ValueError('Missing regular freeze input: ' + str(path.relative_to(ABM)))
    value = dict(schema_version=1, study=HERE.name, phase=args.phase,
                 created_utc=datetime.now(timezone.utc).isoformat(),
                 files={str(path.relative_to(ABM)): digest(path) for path in sorted(set(paths))},
                 scientific_backend='Frozen response-study backend, unchanged',
                 protocol_scope='Fixed prospective follow-up to known results; not a retrospective preregistration')
    value['identity_sha256'] = hashlib.sha256(canonical(value).encode()).hexdigest()
    target = HERE / ('freeze.json' if args.phase == 'initial' else 'freeze_' + args.phase + '.json')
    with target.open('x') as handle:
        handle.write(json.dumps(value, indent=2, sort_keys=True) + '\n')
    print(json.dumps(dict(phase=args.phase, identity=value['identity_sha256'], files=len(value['files']))))


if __name__ == '__main__':
    main()
