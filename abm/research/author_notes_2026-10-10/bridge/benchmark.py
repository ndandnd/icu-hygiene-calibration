#!/usr/bin/env python3
"""Bounded adapter benchmark only; never runs the scientific full horizon."""
import argparse
import json
from pathlib import Path
import resource
import sys
try:
    from .background import make_background_model
    from .run import atomic_json, simulate_seed, PAIRS
except ImportError:
    from background import make_background_model
    from run import atomic_json, simulate_seed, PAIRS


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ticks',type=int,default=1024)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if not 129<=args.ticks<=4096:parser.error('Bounded benchmark accepts 129..4096 ticks')
    row=simulate_seed(10801,make_background_model(),warmup_ticks=128,observation_ticks=args.ticks-128)
    seconds=row['elapsed_seconds']
    report=dict(description='Tiny implementation benchmark; no study seed/horizon executed',ticks=args.ticks,
                pair_count=len(PAIRS),arm_count=2*len(PAIRS),elapsed_seconds=seconds,
                seconds_per_tick=seconds/args.ticks,peak_rss_native_units=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                projected_full_seed_seconds=seconds/args.ticks*(90+56)*2592,
                projected_256_seed_cpu_hours=seconds/args.ticks*(90+56)*2592*256/3600,
                extrapolation_note='Tiny cold-run extrapolation, not a full-horizon timing guarantee.',
                event_counts=row['event_counts'],ordered_bridge_events_sha256=row['ordered_bridge_events_sha256'])
    atomic_json(args.output,report,immutable=True)
    print(json.dumps(report,indent=2))


if __name__=='__main__':sys.exit(main())
