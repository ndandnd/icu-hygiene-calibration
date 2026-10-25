"""Read only this deployment's receipted scheduler IDs and completion metadata."""
import argparse
import datetime
import json
from pathlib import Path
import subprocess
from submit_once import ALLOWED


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    root=a.root.resolve();cfg=json.loads((root/'launch_config.json').read_text())
    if cfg['root']!=str(root):raise ValueError('Launch root mismatch')
    jobs={};unknown=[]
    for file in (root/'submissions').glob('*.intent.json'):
        name=file.name.removesuffix('.intent.json');intent=json.loads(file.read_text())
        if name not in ALLOWED:raise ValueError('Unexpected submission scope')
        receipt=root/'submissions'/(name+'.jobid')
        if intent.get('state')!='submitted' or not receipt.exists():
            unknown.append(dict(name=name,state=intent.get('state'),comment_token=intent.get('comment_token')));continue
        job=receipt.read_text().strip()
        if not job.isdigit() or job!=intent.get('job_id'):raise ValueError('Receipt/intent mismatch')
        jobs[name]=job
    queue=[];accounting=[]
    if jobs:
        ids=','.join(jobs.values())
        queue=subprocess.run([cfg['squeue'],'--noheader','--array','--jobs',ids,'--format','%i|%j|%T|%R'],
            check=True,text=True,capture_output=True,timeout=30).stdout.splitlines()
        accounting=subprocess.run([cfg['sacct'],'--noheader','--parsable2','-X','--jobs',ids,'--format','JobID,JobName,State,ExitCode'],
            check=True,text=True,capture_output=True,timeout=30).stdout.splitlines()
    counts={name:dict(completed_records=len(list((root/'banks'/name/'completed').glob('job_*.json'))),
                     progress_files=len(list((root/'banks'/name/'progress').glob('*.json')))) for name in ('grid','policy')}
    print(json.dumps(dict(checked_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),jobs=jobs,
        unknown_submission_outcomes=unknown,queue=queue,accounting=accounting,counts=counts),indent=2))


if __name__=='__main__':main()
