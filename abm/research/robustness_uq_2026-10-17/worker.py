"""Restartable partitions using the isolated, explicitly frozen workflow variant."""
import argparse
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import sys
import time
import uuid
from freeze import ABM,HERE,canonical,digest,verify,verify_config

spec=importlib.util.spec_from_file_location('_robustness_backend',ABM/'research/full_icu_response_2026-10-08/run.py')
backend=importlib.util.module_from_spec(spec);sys.modules[spec.name]=backend;spec.loader.exec_module(backend)

def atomic(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
    with tmp.open('x') as f:f.write(canonical(value)+'\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)

def runtime(production):
    versions=backend.versions()
    for name,want in {'Mesa':'1.2.1','numpy':'1.26.4','pandas':'3.0.6','networkx':'3.7'}.items():
        if versions[name]!=want:raise ValueError('Runtime mismatch: '+name)
    if production and (sys.version_info[:3]!=(3,12,3) or not sys.platform.startswith('linux')):
        raise ValueError('Production requires pinned Python3.12.3/Linux')
    if any(os.environ.get(k)!='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')):
        raise ValueError('All numerical thread limits must equal1')
    return versions

def python_executable():
    declared=Path(os.path.abspath(sys.executable));resolved=declared.resolve(strict=True)
    return dict(declared=str(declared),resolved=str(resolved),sha256=digest(resolved))

def production_anchor(path,parent_identity,software):
    """Bind every production partition to the successful compute-node replay."""
    payload=path.read_bytes();probe=json.loads(payload)
    if (probe.get('status')!='passed' or probe.get('production') is not True or
        probe.get('reference_verified') is not True or probe.get('parent_freeze_identity')!=parent_identity):
        raise ValueError('Successful parent-bound production probe required')
    versions=probe.get('versions',{})
    if not versions.get('platform','').startswith('Linux'):
        raise ValueError('Production probe must use Linux')
    if {k:v for k,v in versions.items() if k!='platform'}!={k:v for k,v in software.items() if k!='platform'}:
        raise ValueError('Production runtime differs from probe, including full Python build')
    return hashlib.sha256(payload).hexdigest()

def work(args):
    freeze=verify_config(args.config);parent=verify()
    initial_python=python_executable();probe_path=None;probe_sha=None
    def unchanged():
        if verify()!=parent or verify_config(args.config)!=freeze:
            raise ValueError('Frozen study changed during execution')
        if python_executable()!=initial_python:raise ValueError('Python executable changed within worker stage')
        if probe_path is not None and digest(probe_path)!=probe_sha:raise ValueError('Production probe changed during execution')
    raw=json.loads(args.config.read_text());cfg=backend.normalize_config(raw)
    jobs=backend.build_jobs(cfg)
    is_probe=args.config.resolve()==(HERE/'configs/probe.json').resolve()
    software=runtime(not is_probe)
    if not is_probe:
        if os.environ.get('ROBUST_PYTHON')!=initial_python['declared']:
            raise ValueError('Worker interpreter differs from declared entry path')
        output=Path(os.path.abspath(args.output))
        if output.parent.name!='banks' or output.name not in ('grid','policy') or args.config.stem!=output.name:
            raise ValueError('Production bank must be ROOT/banks/grid or ROOT/banks/policy')
        probe_path=output.parent.parent/'probe/replay.json'
        probe_sha=production_anchor(probe_path,parent['identity_sha256'],software)
    unchanged()
    source_files=dict(parent['files']);source_files.update(freeze['files'])
    manifest=dict(schema_version=1,study=HERE.name,config=raw,resolved_config=cfg,jobs=jobs,
        sources=source_files,software={k:v for k,v in software.items() if k!='platform'},
        task_count=args.task_count,freeze_identity=freeze['identity_sha256'],
        parent_freeze_identity=parent['identity_sha256'],model_factory='WorkflowRobustnessModel',
        production_probe_sha256=probe_sha,
        checkpoint='Atomic complete records; incomplete record restarts from seed; no biological checkpoint.')
    manifest['identity']=hashlib.sha256(canonical(manifest).encode()).hexdigest()
    args.output.mkdir(parents=True,exist_ok=True)
    with (args.output/'.manifest.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        path=args.output/'manifest.json'
        if path.exists():
            if json.loads(path.read_text())!=manifest:raise ValueError('Source/config/runtime/partition changed on resume')
        else:atomic(path,manifest)
    for d in ('completed','progress','attempts','locks'):(args.output/d).mkdir(exist_ok=True)
    with (args.output/'locks'/f'{args.task_index}.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        from workflow_variant import WorkflowRobustnessModel
        unchanged()
        start=time.monotonic();deadline=start+args.wall_seconds
        stop={'requested':False};done=skipped=0
        def request(*unused):stop['requested']=True
        def terminate(*unused):
            stop['requested']=True
            raise TimeoutError('Interrupted; incomplete record will restart')
        signal.signal(signal.SIGUSR1,request);signal.signal(signal.SIGTERM,terminate)
        def progress(status,job=None):
            atomic(args.output/'progress'/f'{args.task_index}.json',dict(identity=manifest['identity'],
                task=args.task_index,status=status,done=done,skipped=skipped,current_job=job,
                python_executable=initial_python,production_probe_sha256=probe_sha,
                elapsed_seconds=time.monotonic()-start))
        progress('running')
        for job in jobs:
            if job['job_id']%args.task_count!=args.task_index:continue
            path=args.output/'completed'/f"job_{job['job_id']:08d}.json"
            if path.exists():
                old=json.loads(path.read_text())
                if old.get('identity')!=manifest['identity'] or old.get('status')!='complete' or any(old.get(k)!=v for k,v in job.items()):
                    raise ValueError('Completed record does not match expected identity')
                skipped+=1;continue
            if stop['requested'] or time.monotonic()>=deadline:
                progress('restart_required')
                return 75
            unchanged()
            progress('running',job['job_id'])
            row=backend.simulate_pair(job,cfg,deadline=deadline,model_factory=WorkflowRobustnessModel)
            current_runtime=runtime(not is_probe)
            if current_runtime!=software:raise ValueError('Runtime changed within worker stage')
            row.update(identity=manifest['identity'],execution=current_runtime,python_executable=initial_python)
            unchanged()
            if row.get('status')!='complete':
                atomic(args.output/'attempts'/f"job_{job['job_id']:08d}_{uuid.uuid4().hex}.json",row)
                progress('restart_required' if row.get('status')=='timeout' else 'error',job['job_id'])
                return 75 if row.get('status')=='timeout' else 2
            unchanged()
            atomic(path,row);done+=1;progress('running')
            print(canonical(dict(task=args.task_index,done=done,skipped=skipped)),flush=True)
        unchanged();progress('complete');return 0

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--task-index',type=int,required=True);p.add_argument('--task-count',type=int,required=True)
    p.add_argument('--wall-seconds',type=float,default=3300)
    a=p.parse_args()
    if not 0<=a.task_index<a.task_count or not 0<a.wall_seconds<86400:p.error('Invalid partition or wall budget')
    try:return work(a)
    except TimeoutError:return 75

if __name__=='__main__':sys.exit(main())
