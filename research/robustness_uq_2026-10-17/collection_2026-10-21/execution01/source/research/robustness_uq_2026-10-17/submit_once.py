"""Durable, fail-closed submission for this study's five named jobs only."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import uuid

ALLOWED={'probe','grid','selection','policy','analysis'}


def atomic(path,value):
    payload=value if isinstance(value,str) else json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
    with temp.open('x') as f:f.write(payload);f.flush();os.fsync(f.fileno())
    os.replace(temp,path)
    fd=os.open(path.parent,os.O_RDONLY)
    try:os.fsync(fd)
    finally:os.close(fd)


def submit(root,name,arguments,sbatch):
    if name not in ALLOWED:raise ValueError('Unapproved job name')
    if not Path(sbatch).is_absolute():raise ValueError('Absolute sbatch path required')
    folder=Path(root)/'submissions';folder.mkdir(parents=True,exist_ok=True)
    receipt=folder/(name+'.jobid');intent=folder/(name+'.intent.json')
    script=Path(arguments[-1])
    if not script.is_file():raise ValueError('Submission script absent')
    declaration=dict(name=name,sbatch=str(sbatch),arguments=arguments,
        script_sha256=hashlib.sha256(script.read_bytes()).hexdigest())
    with (folder/(name+'.submit.lock')).open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if receipt.exists():
            job=receipt.read_text().strip()
            old=json.loads(intent.read_text()) if intent.exists() else {}
            if not job.isdigit() or old.get('state')!='submitted' or old.get('job_id')!=job or old.get('declaration')!=declaration:
                raise ValueError('Submission receipt/declaration mismatch')
            return job
        if intent.exists():
            old=json.loads(intent.read_text())
            if old.get('declaration')==declaration and old.get('state')=='submitted' and str(old.get('job_id','')).isdigit():
                atomic(receipt,old['job_id']+'\n');return old['job_id']
            raise ValueError('Unknown prior submission outcome; diagnose receipt/accounting before recovery')
        state=dict(state='pending_unknown_until_receipted',declaration=declaration,
                   comment_token='icu-robust-'+uuid.uuid4().hex)
        atomic(intent,state)
        command=[str(sbatch),'--parsable','--comment='+state['comment_token'],*arguments]
        try:
            result=subprocess.run(command,text=True,capture_output=True,timeout=60)
        except BaseException as error:
            atomic(folder/(name+'.submission-error.json'),dict(type=type(error).__name__,message=str(error),
                outcome='unknown; never auto-resubmit'))
            raise
        atomic(folder/(name+'.submission-output.json'),dict(returncode=result.returncode,stdout=result.stdout,
                                                          stderr=result.stderr,command=command))
        match=re.fullmatch(r'([1-9][0-9]*)(?:;[A-Za-z0-9_.-]+)?',result.stdout.strip())
        if result.returncode or match is None:
            raise ValueError('No confirmed numeric submission result; intent retained, no retry')
        job=match.group(1);state.update(state='submitted',job_id=job)
        atomic(intent,state);atomic(receipt,job+'\n')
        return job
