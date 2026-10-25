#!/usr/bin/env python3
"""Stop retryable infrastructure stalls; never replace an unfinished seed."""
import argparse
import json
from pathlib import Path


def update(output,task,module):
    if module=='bridge':
        done=len(list((output/f'part_{task:02d}'/'completed').glob('*.json')))
    else:
        progress=json.loads((output/'progress'/f'{task}.json').read_text())
        done=progress['done']+progress['skipped']
    root=output/'retry_state';root.mkdir(exist_ok=True)
    path=root/f'{task}.json'
    previous=json.loads(path.read_text()) if path.exists() else dict(completed=0,consecutive_stalls=0)
    stalls=previous['consecutive_stalls']+1 if done<=previous['completed'] else 0
    value=dict(completed=done,consecutive_stalls=stalls,status='stop_stalled' if stalls>=3 else 'retry_allowed')
    temp=path.with_suffix('.json.tmp');temp.write_text(json.dumps(value,sort_keys=True)+'\n');temp.replace(path)
    return value


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--task',type=int,required=True);p.add_argument('--module',choices=['clock_grid','surfaces_off','independent_upper','clock_boundary','clock_policy'],required=True)
    args=p.parse_args();value=update(args.output,args.task,args.module);print(json.dumps(value,sort_keys=True))
    return 2 if value['status']=='stop_stalled' else 0


if __name__=='__main__':raise SystemExit(main())
