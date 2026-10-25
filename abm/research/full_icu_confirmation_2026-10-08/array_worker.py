#!/usr/bin/env python3
"""Bind new confirmation banks to a frozen plan; reuse the unchanged backend."""
import argparse
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]
OLD = ABM / 'research/full_icu_response_2026-10-08'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def verify_freeze(path, root=ABM):
    freeze = json.loads(path.read_text())
    expected = freeze['identity_sha256']
    payload = {k: v for k, v in freeze.items() if k != 'identity_sha256'}
    if hashlib.sha256(canonical(payload).encode()).hexdigest() != expected:
        raise ValueError('Invalid freeze identity')
    if freeze.get('study') != HERE.name:
        raise ValueError('Wrong study')
    for name, sha in freeze['files'].items():
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Unsafe source path')
        source = root / relative
        if source.is_symlink() or digest(source) != sha:
            raise ValueError('Frozen file changed: ' + name)
    return freeze


def binding(args):
    freeze = verify_freeze(args.freeze)
    config_name = str(args.config.resolve().relative_to(ABM))
    if freeze['files'].get(config_name) != digest(args.config):
        raise ValueError('Config is not frozen')
    scope = json.loads((HERE / 'output_scopes.json').read_text())
    name = args.config.stem
    output_name = str(args.output.resolve().relative_to(ABM))
    if scope['banks'].get(name) != output_name:
        raise ValueError('Output is outside its declared bank')
    return dict(schema_version=1, study=HERE.name,
                freeze_sha256=digest(args.freeze), freeze_identity=freeze['identity_sha256'],
                config_path=config_name, config_sha256=digest(args.config),
                seedbook_sha256=digest(HERE / 'seedbook.json'),
                output_scope=output_name, task_count=args.task_count,
                backend='research/full_icu_response_2026-10-08/array_worker.py',
                note='Unchanged backend; this receipt identifies the new study.')


def work(args):
    receipt = binding(args)
    args.output.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(OLD))
    spec = importlib.util.spec_from_file_location('_confirmation_frozen_worker', OLD / 'array_worker.py')
    worker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(worker)
    with (args.output / '.confirmation.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        path = args.output / 'confirmation_binding.json'
        if path.exists():
            if json.loads(path.read_text()) != receipt:
                raise ValueError('Changed confirmation binding: resume refused')
        else:
            worker.atomic(path, receipt)
    result = worker.work(args)
    if binding(args) != receipt:
        raise ValueError('Confirmation inputs changed during execution')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--freeze', type=Path, default=Path(os.environ.get('ICU_FREEZE', HERE / 'freeze.json')))
    p.add_argument('--task-index', type=int, required=True)
    p.add_argument('--task-count', type=int, required=True)
    p.add_argument('--wall-seconds', type=float, default=3300)
    args = p.parse_args()
    if not 0 <= args.task_index < args.task_count or not 0 < args.wall_seconds < 86400:
        p.error('Invalid partition or wall budget')
    try:
        return work(args)
    except TimeoutError:
        return 75


if __name__ == '__main__':
    sys.exit(main())
