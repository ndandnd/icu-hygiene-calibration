"""Bounded one-time probe→grid→selection→policy→analysis dispatch.

Deployment is ROOT/source (the ABM directory); no private host/path is hardcoded.
Initial --stage start requires --root, --python, --sbatch, --scontrol, --squeue,
--sacct, --account, --partition and --reference; --nodelist is optional.
Only start submits initial jobs. Auxiliary stages use immutable launch_config.
This module never changes another project's jobs or recovers unknown submissions.
"""
import argparse
import datetime
import fcntl
import getpass
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from freeze import ABM,HERE,check,digest,identity,verify,write_new
from submit_once import submit


def path(value):
    p=Path(os.path.abspath(value))
    if any(x.is_symlink() for x in (p,*p.parents)):raise ValueError('Symlink path refused')
    if any(c in str(p) for c in (',','\n','\r','\x00')):raise ValueError('Unsupported path character')
    return p


def read(p):return json.loads(Path(p).read_text())


def executable(value):
    """Retain a venv entry path, while separately binding its real binary."""
    declared=Path(os.path.abspath(value))
    if any(c in str(declared) for c in (',','\n','\r','\x00')):raise ValueError('Unsupported executable path')
    resolved=declared.resolve(strict=True)
    if not resolved.is_file() or not os.access(declared,os.X_OK):raise ValueError('Executable is unavailable')
    return dict(declared=str(declared),resolved=str(resolved),sha256=digest(resolved))


def guard(root,cfg,stage_python=None):
    if path(root/'source')!=ABM or path(root/'source/research'/HERE.name)!=HERE:
        raise ValueError('Use this deployment under ROOT/source')
    parent=verify()
    if parent['identity_sha256']!=cfg['parent_freeze_identity']:raise ValueError('Parent freeze differs from launch declaration')
    if digest(path(cfg['reference']))!=cfg['reference_sha256']:raise ValueError('Replay reference changed')
    for name,binding in cfg['executables'].items():
        actual=executable(cfg[name])
        if name=='python':
            if actual['declared']!=binding['declared']:raise ValueError('Python entry path changed')
            if stage_python is not None and actual!=stage_python:
                raise ValueError('Python executable changed within stage')
        elif actual!=binding:raise ValueError('Executable identity changed: '+name)
    return parent


def priority_guard(root,cfg):
    result=subprocess.run([cfg['squeue'],'--noheader','--user',getpass.getuser(),'--format','%i|%j|%T'],
                          check=True,text=True,capture_output=True,timeout=30)
    rows=[line.split('|') for line in result.stdout.splitlines() if line.strip()]
    if any(len(row)!=3 for row in rows):raise ValueError('Unparseable priority queue check')
    blocked=[row for row in rows if 'opte' in row[1].lower()]
    receipt=dict(checked_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 scope='Current user OPTE job-name guard; no scheduler mutation',active_opte=blocked,
                 scanned_job_count=len(rows),squeue_binary=cfg['executables']['squeue'])
    # Each attempt is retained; successful dispatch remains one-time.
    from submit_once import atomic
    import uuid
    atomic(root/'priority_checks'/('check-'+uuid.uuid4().hex+'.json'),receipt)
    if blocked:raise ValueError('Active OPTE work defers new ABM submissions')


def await_stage_receipt(root,name,job_id,timeout=30):
    """A fast-starting job may precede sbatch's durable returned-ID receipt."""
    if not job_id.isdigit():raise ValueError('Numeric SLURM_JOB_ID required')
    receipt=root/'submissions'/(name+'.jobid');intent=root/'submissions'/(name+'.intent.json')
    deadline=time.monotonic()+timeout
    while True:
        state=read(intent) if intent.exists() else None
        if receipt.exists():
            # Submission publishes the submitted intent before its job-ID file.
            # The receipt may have appeared after the first intent read above.
            state=read(intent) if intent.exists() else None
            if receipt.read_text().strip()!=job_id or state is None or state.get('state')!='submitted' or state.get('job_id')!=job_id:
                raise ValueError('Auxiliary stage receipt/intent differs from its SLURM_JOB_ID')
            return
        if state is not None and state.get('state')=='submitted' and state.get('job_id')!=job_id:
            raise ValueError('Auxiliary stage intent belongs to another job')
        remaining=deadline-time.monotonic()
        if remaining<=0:raise ValueError('Auxiliary stage receipt not durable within30seconds; no resubmission')
        time.sleep(min(.1,remaining))


def run(cfg,args):
    subprocess.run([cfg['python'],*map(str,args)],check=True)


def exported(items):
    for value in items.values():
        if any(c in str(value) for c in (',','\n','\r','\x00')):raise ValueError('Unsupported exported value')
    return '--export='+','.join(k+'='+str(v) for k,v in sorted(items.items()))


def dispatch(root,cfg,name,dependency=None):
    guard(root,cfg)
    common=['--account='+cfg['account'],'--partition='+cfg['partition'],
        '--nodes=1','--ntasks=1','--cpus-per-task=1','--mem=1G','--nice=10000',
        '--job-name=icu-robust-'+name,'--output='+str(root/'logs'/(name+'-%A_%a.log'))]
    if cfg.get('nodelist'):common+=['--nodelist='+cfg['nodelist']]
    if dependency is not None:
        if not str(dependency).isdigit():raise ValueError('Numeric dependency required')
        common+=['--dependency=afterok:'+str(dependency)]
    env=dict(ROBUST_PYTHON=cfg['python'],ROBUST_STUDY=str(HERE),ROBUST_ROOT=str(root))
    if name in ('grid','policy'):
        env.update(ROBUST_CONFIG=str(HERE/'configs'/(name+'.json')),ROBUST_OUTPUT=str(root/'banks'/name),
                   ROBUST_TASKS='64',ROBUST_SCONTROL=cfg['scontrol'])
        arguments=[*common,'--array=0-63%16','--time=03:00:00',exported(env),str(HERE/'array.sbatch')]
    else:
        stage={'probe':'probe','selection':'select-and-dispatch','analysis':'analyze'}[name]
        env['ROBUST_STAGE']=stage
        arguments=[*common,'--no-requeue','--time='+('00:05:00' if name=='probe' else '00:30:00'),
                   exported(env),str(HERE/'auxiliary.sbatch')]
    return submit(root,name,arguments,cfg['sbatch'])


def start(args,root):
    for key in ('python','sbatch','scontrol','squeue','sacct','reference','account','partition'):
        if not getattr(args,key):raise ValueError('Initial launch requires --'+key)
    if sys.version_info[:3]!=(3,12,3) or not sys.platform.startswith('linux'):
        raise ValueError('Start under the pinned Linux Python3.12.3 runtime')
    parent=verify()
    bindings={key:executable(getattr(args,key)) for key in ('python','sbatch','scontrol','squeue','sacct')}
    cfg=dict(schema_version=1,study=HERE.name,root=str(root),parent_freeze_identity=parent['identity_sha256'],
        **{key:value['declared'] for key,value in bindings.items()},executables=bindings,reference=str(path(args.reference)),
        account=args.account,partition=args.partition,nodelist=args.nodelist)
    if os.path.abspath(sys.executable)!=cfg['python']:raise ValueError('Launch interpreter differs from declared runtime')
    for key in ('account','partition','nodelist'):
        if cfg[key] is not None and (not cfg[key] or any(c.isspace() for c in cfg[key])):raise ValueError('Invalid scheduler field')
    cfg['reference_sha256']=digest(path(cfg['reference']))
    reference=read(cfg['reference'])
    if reference.get('status')!='passed' or len(reference.get('scientific_record_sha256',[]))!=9:
        raise ValueError('Nine-fixture successful local reference required')
    guard(root,cfg)
    root.mkdir(parents=True,exist_ok=True)
    for name in ('logs','submissions','banks','probe','outputs'):(root/name).mkdir(exist_ok=True)
    write_new(root/'launch_config.json',cfg)
    with (root/'launch.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        priority_guard(root,cfg)
        probe=dispatch(root,cfg,'probe');grid=dispatch(root,cfg,'grid',probe)
        selection=dispatch(root,cfg,'selection',grid)
        write_new(root/'initial_dispatch.json',dict(parent_freeze_identity=parent['identity_sha256'],
                  jobs=dict(probe=probe,grid=grid,selection=selection)))
        guard(root,cfg,stage_python=bindings['python'])
    print(json.dumps(dict(probe=probe,grid=grid,selection=selection)))


def auxiliary(args,root):
    cfg=read(root/'launch_config.json')
    if cfg['root']!=str(root) or cfg['study']!=HERE.name:raise ValueError('Launch configuration scope differs')
    if os.path.abspath(sys.executable)!=cfg['python']:raise ValueError('Auxiliary interpreter differs from declaration')
    parent=guard(root,cfg)
    from worker import runtime
    actual_runtime=runtime(True)
    stage_name={'probe':'probe','select-and-dispatch':'selection','analyze':'analysis'}[args.stage]
    job_id=os.environ.get('SLURM_JOB_ID','')
    await_stage_receipt(root,stage_name,job_id)
    stage_binding=dict(stage=args.stage,job_id=job_id,python=sys.version,runtime=actual_runtime,
        python_executable=executable(cfg['python']),launch_config_sha256=digest(root/'launch_config.json'),
        parent_freeze_identity=parent['identity_sha256'])
    write_new(root/'stages'/(stage_name+'.start.json'),stage_binding)
    def finish(**details):
        guard(root,cfg,stage_python=stage_binding['python_executable'])
        write_new(root/'stages'/(stage_name+'.finish.json'),dict(stage_binding,status='complete',**details))
    if args.stage=='probe':
        output=root/'probe/replay.json'
        run(cfg,[HERE/'probe.py','--output',output,'--reference',cfg['reference'],'--production'])
        receipt=read(output)
        if receipt.get('status')!='passed' or receipt.get('reference_verified') is not True:
            raise ValueError('Production replay did not pass')
        finish(probe_receipt_sha256=digest(output))
        return
    if args.stage=='select-and-dispatch':
        with (root/'selection_dispatch.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            selection=HERE/'selection.json';policy=HERE/'configs/policy.json'
            run(cfg,[HERE/'analyze.py','--stage','select','--bank',root/'banks/grid',
                '--config',HERE/'configs/grid.json','--output',selection,'--policy-seeds',HERE/'seedbook.json',
                '--policy-config',policy,'--probe',root/'probe/replay.json'])
            chosen=read(selection)
            if chosen['parent_freeze_identity']!=parent['identity_sha256'] or chosen['policy_outcomes_used']:
                raise ValueError('Selection parent/scope mismatch')
            files={str(selection.relative_to(ABM)):digest(selection)}
            if chosen['policy_config'] is not None:
                if read(policy)!=chosen['policy_config']:raise ValueError('Derived policy differs from selected configuration')
                if len(chosen['policy_config']['scenarios'])>3 or len(chosen['policy_config']['seeds'])!=512:
                    raise ValueError('Selected workload exceeds fixed ceiling')
                files[str(policy.relative_to(ABM))]=digest(policy)
            elif policy.exists():raise ValueError('Policy config exists after complete selection failure')
            derived=dict(schema_version=1,study=HERE.name,parent_freeze_identity=parent['identity_sha256'],
                         files=files,grid_input_receipt=chosen['input_receipt'],selection_runtime=chosen['analysis_runtime'])
            derived['identity_sha256']=identity(derived)
            write_new(HERE/'policy_freeze.json',derived);check(derived);guard(root,cfg)
            if chosen['policy_config'] is None:
                write_new(root/'outputs/final_status.json',dict(status='complete_selection_failure_no_policy',
                    parent_freeze_identity=parent['identity_sha256'],selection_sha256=digest(selection),
                    reason='All prespecified selections failed; no fallback, no policy or analysis jobs.',
                    primary_contrasts=[dict(name=x,estimate=None,interval=None,reason='Selection failed; family remains three')
                        for x in ('independent-clustered','independent-imperfect','clustered-imperfect')]))
                finish(selection_sha256=digest(selection),policy_pairs=0)
                print('All selections failed; policy not submitted.');return
            policy_job=dispatch(root,cfg,'policy');analysis=dispatch(root,cfg,'analysis',policy_job)
            write_new(root/'policy_dispatch.json',dict(policy=policy_job,analysis=analysis,
                selected_profiles=[s['id'] for s in chosen['policy_config']['scenarios']],
                policy_pairs=512*len(chosen['policy_config']['scenarios']),derived_freeze_identity=derived['identity_sha256']))
            finish(selection_sha256=digest(selection),derived_freeze_identity=derived['identity_sha256'],policy=policy_job,analysis=analysis)
            print(json.dumps(dict(policy=policy_job,analysis=analysis)))
        return
    if args.stage=='analyze':
        output=root/'outputs/analysis.json'
        run(cfg,[HERE/'analyze.py','--stage','analyze','--bank',root/'banks/policy',
            '--config',HERE/'configs/policy.json','--output',output,'--selection',HERE/'selection.json',
            '--derived-freeze',HERE/'policy_freeze.json','--probe',root/'probe/replay.json'])
        guard(root,cfg)
        write_new(root/'outputs/final_status.json',dict(status='execution_and_analysis_complete_awaiting_independent_verification',
            analysis_sha256=digest(output),parent_freeze_identity=parent['identity_sha256'],runtime=actual_runtime))
        finish(analysis_sha256=digest(output))
        return
    raise ValueError('Unknown auxiliary stage')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=('start','probe','select-and-dispatch','analyze'),required=True)
    p.add_argument('--root',required=True)
    for key in ('python','sbatch','scontrol','squeue','sacct','reference','account','partition','nodelist'):p.add_argument('--'+key)
    args=p.parse_args();root=path(args.root)
    if args.stage=='start':start(args,root)
    else:auxiliary(args,root)


if __name__=='__main__':main()
