#!/usr/bin/env python3
"""Bind this new source/configuration closure without modifying old releases."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(path=HERE / 'freeze.json', root=ABM):
    frozen = json.loads(path.read_text())
    core = {k:v for k,v in frozen.items() if k != 'identity_sha256'}
    if digest_bytes(canonical(core).encode()) != frozen['identity_sha256']:
        raise ValueError('Invalid freeze identity')
    if frozen['study'] != HERE.name:
        raise ValueError('Wrong study identity')
    for name, expected in frozen['files'].items():
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Unsafe frozen path')
        target = root / relative
        if target.is_symlink() or digest(target) != expected:
            raise ValueError('Frozen file changed: ' + name)
    return frozen


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def main():
    files = [p for p in HERE.rglob('*') if p.is_file() and p.suffix in ('.py','.md','.json','.sbatch')
             and p.name != 'freeze.json' and not any(x in p.parts for x in ('outputs','qa','__pycache__'))]
    dependencies = ['src/CPEmodel.py','src/agents.py','src/routes.py',
        'research/full_icu_calibration_2026-10-07/run.py',
        'research/full_icu_response_2026-10-08/run.py',
        'research/full_icu_response_2026-10-08/array_worker.py',
        'research/full_icu_response_2026-10-08/instrumented_model.py',
        'research/full_icu_response_2026-10-08/PROTOCOL.md',
        'research/full_icu_confirmation_2026-10-08/bridge/background.py',
        'research/full_icu_confirmation_2026-10-08/bridge/coefficient.py',
        'research/full_icu_confirmation_2026-10-08/bridge/run.py',
        'research/full_icu_confirmation_2026-10-08/bridge/test_bridge.py']
    files += [ABM / p for p in dependencies]
    bridge_closure = json.loads((HERE / 'bridge/dependency_paths.json').read_text())['files']
    for name, expected in bridge_closure.items():
        path=ABM/name
        if digest(path)!=expected:raise ValueError('Bridge dependency drift: '+name)
        files.append(path)
    value = dict(schema_version=1, study=HERE.name,
                 files={str(p.relative_to(ABM)):digest(p) for p in sorted(set(files))},
                 preservation='All older study freezes remain unchanged.')
    value['identity_sha256'] = digest_bytes(canonical(value).encode())
    path = HERE / 'freeze.json'
    if path.exists() and json.loads(path.read_text()) != value:
        raise ValueError('Freeze already exists with different bytes')
    path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
    verify()
    print(value['identity_sha256'])


if __name__ == '__main__':
    main()
