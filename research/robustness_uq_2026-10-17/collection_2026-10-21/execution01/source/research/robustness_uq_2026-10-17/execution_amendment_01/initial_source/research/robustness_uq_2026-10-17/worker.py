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

def work(args):
    freeze=verify_config(args.config);parent=verify()
    def unchanged():
        if verify()!=parent or verify_config(args.config)!=freeze:
            raise ValueError('Frozen study changed during execution')
    raw=json.loads(args.config.read_text());cfg=backend.normalize_config(raw)
    jobs=backend.build_jobs(cfg)
    is_probe=args.config.resolve()==(HERE/'configs/probe.json').resolve()
    software=runtime(not is_probe)
    source_files=dict(parent['files']);source_files.update(freeze['files'])
    manifest=dict(schema_version=1,study=HERE.name,config=raw,resolved_config=cfg,jobs=jobs,
        sources=source_files,software={k:v for k,v in software.items() if k!='platform'},
        task_count=args.task_count,freeze_identity=freeze['identity_sha256'],
        parent_freeze_identity=parent['identity_sha256'],model_factory='WorkflowRobustnessModel',
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
            row.update(identity=manifest['identity'],execution=runtime(not is_probe))
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
