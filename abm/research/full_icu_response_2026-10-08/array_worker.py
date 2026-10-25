#!/usr/bin/env python3
"""Hash-bound, restartable serial partitions of paired ICU simulations."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import time
import uuid
import run

def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)

def atomic(path,value):
    temp=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
    with temp.open('x') as stream:
        stream.write(canonical(value)+'\n');stream.flush();os.fsync(stream.fileno())
    temp.replace(path)

def work(args):
    raw=json.loads(args.config.read_text()); config=run.normalize_config(raw)
    jobs=run.build_jobs(config)
    sources=run.source_hashes()
    sources['array_worker.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    manifest=dict(schema_version=1,config=raw,resolved_config=config,jobs=jobs,
        sources=sources,software={k:v for k,v in run.versions().items() if k!='platform'},
        task_count=args.task_count,checkpoint='Atomic completed pairs; restart incomplete pairs from seed.')
    identity=hashlib.sha256(canonical(manifest).encode()).hexdigest()
    manifest['identity']=identity
    args.output.mkdir(parents=True,exist_ok=True)
    with (args.output/'.manifest.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        path=args.output/'manifest.json'
        if path.exists():
            if json.loads(path.read_text())!=manifest:raise ValueError('Changed config/source/software/partition: resume refused')
        else:atomic(path,manifest)
    for d in ('completed','attempts','locks','progress'):(args.output/d).mkdir(exist_ok=True)
    lock=(args.output/'locks'/f'{args.task_index}.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    stop={'requested':False};start=time.monotonic();deadline=start+args.wall_seconds
    def request(number,frame):stop['requested']=True
    def terminate(number,frame):
        stop['requested']=True
        raise TimeoutError('Scheduler interruption; incomplete pair will restart')
    signal.signal(signal.SIGUSR1,request);signal.signal(signal.SIGTERM,terminate)
    done=skipped=0
    def progress(status,job=None):
        atomic(args.output/'progress'/f'{args.task_index}.json',dict(identity=identity,
            task=args.task_index,status=status,done=done,skipped=skipped,current_job=job,
            elapsed_seconds=time.monotonic()-start))
    progress('running')
    for job in jobs:
        if job['job_id']%args.task_count!=args.task_index:continue
        if stop['requested'] or time.monotonic()>=deadline:
            progress('restart_required');return 75
        path=args.output/'completed'/f"job_{job['job_id']:08d}.json"
        if path.exists():
            old=json.loads(path.read_text())
            if old.get('identity')!=identity or old.get('status')!='complete' or any(old.get(k)!=job[k] for k in job):
                raise ValueError('Existing completed pair does not match this job')
            skipped+=1;continue
        progress('running',job['job_id'])
        try:
            row=run.simulate_pair(job,config,deadline=deadline)
        except TimeoutError:
            progress('restart_required',job['job_id']);return 75
        row.update(identity=identity,execution=run.versions())
        if row.get('status')!='complete':
            atomic(args.output/'attempts'/f"job_{job['job_id']:08d}_{uuid.uuid4().hex}.json",row)
            if row.get('status')=='timeout':progress('restart_required',job['job_id']);return 75
            progress('error',job['job_id']);return 2
        # Also verify identity after execution, before publishing a completed row.
        current=run.source_hashes();current['array_worker.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        if current!=sources:raise ValueError('Scientific source changed during execution')
        atomic(path,row);done+=1;progress('running')
        print(canonical(dict(task=args.task_index,complete=done,skipped=skipped)),flush=True)
    progress('complete');return 0

def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--task-index',type=int,required=True)
    p.add_argument('--task-count',type=int,required=True);p.add_argument('--wall-seconds',type=float,default=3300)
    a=p.parse_args()
    if not 0<=a.task_index<a.task_count or not 0<a.wall_seconds<86400:p.error('Invalid partition or wall budget')
    try:
        return work(a)
    except TimeoutError:
        # A scheduler signal can arrive between pairs or during an atomic write.
        # Complete records remain immutable; the next process validates and skips them.
        return 75

if __name__=='__main__':sys.exit(main())
