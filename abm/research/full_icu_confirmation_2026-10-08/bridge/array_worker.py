#!/usr/bin/env python3
"""Frozen outer binding for one portable p=0 partition; incomplete seeds restart."""
import argparse
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
import sys

HERE = Path(__file__).resolve().parent
NEW = HERE.parent
ABM = HERE.parents[2]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def validate_partition(manifest, audit, site, index):
    """Pure exact sample/config/scope validation; no model import or simulation."""
    if manifest.get('study') != NEW.name or manifest.get('module') != 'bridge':
        raise ValueError('Wrong partition study/module')
    banks = audit['bridge_seed_blocks']
    if set(banks) != {'A', 'B'} or any(len(banks[b]) != 128 for b in ('A', 'B')):
        raise ValueError('Expected two 128-seed bridge banks')
    all_seeds = banks['A']+banks['B']
    if any(type(s) is not int or not 0 < s < 2**32 for s in all_seeds) or len(set(all_seeds)) != 256:
        raise ValueError('Invalid or overlapping bridge seeds')
    if manifest['seed_blocks'] != banks or manifest['all_seeds'] != all_seeds:
        raise ValueError('Bridge seedbook disagrees with audit')
    allocations = {'local': banks['A'][:16]+banks['B'][:16],
                   'cluster': banks['A'][16:]+banks['B'][16:]}
    if manifest['site_allocations'] != allocations:
        raise ValueError('Changed 1/8 local, 7/8 cluster allocation')
    expected_science = {'description': 'Fixed 256 fresh p=0 backgrounds for exact finite-horizon ICU coefficient; no positive-p runs.',
        'warmup_ticks': 90*2592, 'observation_ticks': 56*2592, 'h0': .9, 'h1': .95,
        'opportunities': 'zone_separate', 'histogram_cap': 2,
        'profiles': [{'id': 'iid_perfect', 'rho': 0., 'e': 1.},
                     {'id': 'persistent_perfect', 'rho': .5, 'e': 1.},
                     {'id': 'iid_imperfect', 'rho': 0., 'e': .8}],
        'model_overrides': {'prob_patient_sick': .03}}
    if manifest['scientific_config'] != expected_science:
        raise ValueError('Changed bridge scientific settings')
    counts = {'local': 6, 'cluster': 32}
    for name, count in counts.items():
        parts = manifest['partitions'][name]
        if len(parts) != count: raise ValueError('Incorrect partition count')
        expected = [{'index': i,
            'config_path': f'research/{NEW.name}/bridge/configs/{name}_part_{i:02d}.json',
            'seeds': allocations[name][i::count],
            'output_scope': f'results/{NEW.name}/bridge_{name}/part_{i:02d}'} for i in range(count)]
        if parts != expected: raise ValueError('Partition membership/order/scope changed')
    if site not in counts or type(index) is not int or not 0 <= index < counts[site]:
        raise ValueError('Invalid partition')
    return manifest['partitions'][site][index], expected_science


def binding(args):
    primary = load_module('_bridge_outer_freeze_verifier', NEW/'array_worker.py')
    freeze = primary.verify_freeze(args.freeze, ABM)
    if freeze.get('phase') != 'bridge': raise ValueError('Bridge phase freeze required')
    prefix = f'research/{NEW.name}'
    required = [f'{prefix}/array_worker.py', f'{prefix}/seed_audit.json', f'{prefix}/seedbook.json',
                f'{prefix}/bridge/PROTOCOL.md', f'{prefix}/bridge/seedbook.json',
                *[f'{prefix}/bridge/{name}' for name in ('array_worker.py', 'run.py', 'coefficient.py',
                                                        'background.py', 'test_bridge.py')],
                *[f'src/{name}' for name in ('CPEmodel.py', 'agents.py', 'routes.py')]]
    if any(name not in freeze['files'] for name in required):
        raise ValueError('Bridge freeze omits a required binding/scientific file')
    manifest_path = HERE/'seedbook.json'
    manifest = json.loads(manifest_path.read_text())
    audit_path = NEW/'seed_audit.json'
    if manifest['seed_audit_path'] != f'{prefix}/seed_audit.json' or manifest['seed_audit_sha256'] != digest(audit_path):
        raise ValueError('Partition seed audit binding changed')
    part, science = validate_partition(manifest, json.loads(audit_path.read_text()), args.site, args.task_index)
    primary_seedbook = json.loads((NEW/'seedbook.json').read_text())
    primary_seeds = {seed for banks in primary_seedbook['seed_blocks'].values()
                     for values in banks.values() for seed in values}
    if primary_seeds.intersection(manifest['all_seeds']):
        raise ValueError('Bridge/primary seed overlap')
    config_path = ABM/part['config_path']
    if freeze['files'].get(part['config_path']) != digest(config_path):
        raise ValueError('Partition config is not frozen')
    config = json.loads(config_path.read_text())
    if config != dict(science, seeds=part['seeds']):
        raise ValueError('Partition config disagrees with exact science/seed membership')
    dispatch_name = f'{prefix}/bridge/configs/bridge_{args.site}.json'
    dispatch_path = ABM/dispatch_name
    if freeze['files'].get(dispatch_name) != digest(dispatch_path):
        raise ValueError('Site dispatcher config is not frozen')
    if json.loads(dispatch_path.read_text()) != dict(science, seeds=manifest['site_allocations'][args.site]):
        raise ValueError('Site dispatcher config membership changed')
    supplied_config = getattr(args, 'config', None)
    if supplied_config is not None and supplied_config.resolve() != dispatch_path.resolve():
        raise ValueError('Dispatcher must use the exact declared site config')
    supplied_count = getattr(args, 'task_count', None)
    if supplied_count is not None and supplied_count != len(manifest['partitions'][args.site]):
        raise ValueError('Dispatcher task-count changed')
    supplied_output = getattr(args, 'output', None)
    dispatch_output = ABM/f'results/{NEW.name}/bridge_{args.site}'
    if supplied_output is not None and supplied_output.resolve() != dispatch_output.resolve():
        raise ValueError('Dispatcher output is outside its declared site scope')
    output = ABM/part['output_scope']
    if str(output.resolve().relative_to(ABM.resolve())) != part['output_scope']:
        raise ValueError('Unsafe output scope')
    receipt = {'schema_version': 1, 'study': NEW.name, 'module': 'bridge',
        'freeze_identity': freeze['identity_sha256'], 'freeze_sha256': digest(args.freeze),
        'site': args.site, 'partition_index': args.task_index,
        'dispatch_config_path': dispatch_name, 'dispatch_config_sha256': digest(dispatch_path),
        'dispatch_output_scope': f'results/{NEW.name}/bridge_{args.site}',
        'partition_count': len(manifest['partitions'][args.site]),
        'config_path': part['config_path'], 'config_sha256': digest(config_path),
        'protocol_path': f'{prefix}/bridge/PROTOCOL.md', 'protocol_sha256': digest(HERE/'PROTOCOL.md'),
        'seed_audit_path': f'{prefix}/seed_audit.json', 'seed_audit_sha256': digest(audit_path),
        'bridge_seedbook_path': f'{prefix}/bridge/seedbook.json',
        'bridge_seedbook_sha256': digest(manifest_path),
        'seeds': part['seeds'], 'output_scope': part['output_scope'],
        'backend': f'{prefix}/bridge/run.py',
        'scientific_sources': {name: freeze['files'][name] for name in required
                               if name.startswith('src/') or name.endswith(('coefficient.py', 'background.py', '/run.py'))},
        'resume_note': 'Existing backend immutable seed outputs/run lock; incomplete seed restarts. No eventwise positive-p coupling.'}
    return receipt, config_path, output


def write_receipt(output, receipt):
    output.mkdir(parents=True, exist_ok=True)
    with (output/'.bridge_binding.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        path = output/'bridge_binding.json'
        if path.exists():
            if json.loads(path.read_text()) != receipt:
                raise ValueError('Changed outer bridge binding; resume refused')
        else:
            temporary = path.with_suffix('.json.tmp')
            with temporary.open('w') as stream:
                stream.write(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False)+'\n')
                stream.flush(); os.fsync(stream.fileno())
            temporary.replace(path)


def work(args):
    receipt, config_path, output = binding(args)
    write_receipt(output, receipt)
    # Recheck frozen bytes immediately before invoking the reviewed backend.
    if binding(args)[0] != receipt: raise ValueError('Bridge binding changed before execution')
    sys.path.insert(0, str(HERE))
    backend = load_module('_bridge_reviewed_partition_backend', HERE/'run.py')
    result = backend.main(['--config', str(config_path), '--output', str(output),
                           '--abm-root', str(ABM), '--wall-seconds', str(args.wall_seconds)])
    if binding(args)[0] != receipt: raise ValueError('Bridge inputs changed during execution')
    return result


def interrupted(number, frame):
    raise TimeoutError('Scheduler interruption; incomplete bridge seed will restart')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--site', choices=('local', 'cluster'))
    parser.add_argument('--config', type=Path, help='Exact bridge_local.json or bridge_cluster.json dispatcher config')
    parser.add_argument('--output', type=Path, help='Declared results/NEW/bridge_SITE dispatcher scope')
    parser.add_argument('--task-count', type=int, help='Must match the fixed six local or 32 cluster partitions')
    parser.add_argument('--task-index', type=int, required=True)
    parser.add_argument('--freeze', type=Path, default=Path(os.environ.get('ICU_BRIDGE_FREEZE', NEW/'freeze_bridge.json')))
    parser.add_argument('--wall-seconds', type=float, default=3300.)
    args = parser.parse_args(argv)
    if args.site is None:
        names = {'bridge_local.json': 'local', 'bridge_cluster.json': 'cluster'}
        if args.config is None or args.config.name not in names:
            parser.error('Supply --site or an exact bridge_local/cluster.json config')
        args.site = names[args.config.name]
    if not math.isfinite(args.wall_seconds) or not 0 < args.wall_seconds <= 3300:
        parser.error('Wall budget must be finite in (0,3300] seconds')
    try:
        signal.signal(signal.SIGUSR1, interrupted)
        signal.signal(signal.SIGTERM, interrupted)
        return work(args)
    except TimeoutError:
        return 75


if __name__ == '__main__': sys.exit(main())
