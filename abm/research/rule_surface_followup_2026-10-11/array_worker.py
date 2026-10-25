#!/usr/bin/env python3
"""Run new isolated positive-p banks through the unchanged validated backend."""
import argparse
import fcntl
import importlib.util
import json
from pathlib import Path
import sys
from freeze import canonical, digest, verify, verify_config

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]
OLD = ABM / 'research/full_icu_response_2026-10-08'


def binding(args):
    frozen = verify_config(args.config)
    config = str(args.config.resolve().relative_to(ABM))
    if frozen['files'].get(config) != digest(args.config):
        raise ValueError('Configuration is not frozen')
    scope = json.loads((HERE / 'output_scopes.json').read_text())
    output = str(args.output.resolve().relative_to(ABM))
    if scope['banks'].get(args.config.stem) != output:
        raise ValueError('Output is outside its declared bank')
    return dict(study=HERE.name, freeze_identity=frozen['identity_sha256'],
                config_path=config, config_sha256=digest(args.config),
                output_scope=output, task_count=args.task_count)


def work(args):
    receipt = binding(args)
    args.output.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(OLD))
    spec = importlib.util.spec_from_file_location('_author_note_response_worker', OLD / 'array_worker.py')
    worker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(worker)
    with (args.output / '.rule_surface.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        path = args.output / 'study_binding.json'
        if path.exists():
            if json.loads(path.read_text()) != receipt:
                raise ValueError('Study binding changed; resume refused')
        else:
            worker.atomic(path, receipt)
    code = worker.work(args)
    if binding(args) != receipt:
        raise ValueError('Inputs changed during execution')
    return code


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
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
