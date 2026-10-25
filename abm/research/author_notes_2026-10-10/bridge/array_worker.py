#!/usr/bin/env python3
"""One fixed strided subset of the original seed bank; no submission commands."""
import argparse
import json
from pathlib import Path
import sys
try:
    from .run import atomic_json, main as run_main, normalize
except ImportError:
    from run import atomic_json, main as run_main, normalize


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True,help='Parent output scope; appends part_NN')
    parser.add_argument('--task-index',type=int,required=True)
    parser.add_argument('--task-count',type=int,required=True)
    parser.add_argument('--wall-seconds',type=float,default=3300.)
    parser.add_argument('--abm-root',type=Path)
    parser.add_argument('--freeze-only',action='store_true')
    args=parser.parse_args(argv)
    raw=json.loads(args.config.read_text())
    config=normalize(raw)
    if not 1<=args.task_count<=len(config['seeds']) or not 0<=args.task_index<args.task_count:
        parser.error('Task count/index outside original declared bank')
    output=args.output/f'part_{args.task_index:02d}'
    output.mkdir(parents=True,exist_ok=True)
    task=dict(seeds=config['seeds'][args.task_index::args.task_count],
              description=f'Original-bank strided partition {args.task_index} of {args.task_count}')
    path=output/'config.json'
    atomic_json(path,task,immutable=True)
    command=['--config',str(path),'--output',str(output),'--wall-seconds',str(args.wall_seconds)]
    if args.abm_root is not None:command+=['--abm-root',str(args.abm_root)]
    if args.freeze_only:command+=['--freeze-only']
    return run_main(command)


if __name__=='__main__':
    try:
        sys.exit(main())
    except (ValueError,RuntimeError,OSError,KeyError,TypeError) as error:
        print(f'Fatal bridge validation: {error}',file=sys.stderr)
        sys.exit(2)
