#!/usr/bin/env python3
"""Short extreme-parameter implementation/replay check; not scientific data."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ABM=HERE.parents[1]
OLD=ABM/'research/full_icu_response_2026-10-08'


def main():
    p=argparse.ArgumentParser();p.add_argument('--write-reference',action='store_true');args=p.parse_args()
    sys.path.insert(0,str(OLD))
    spec=importlib.util.spec_from_file_location('_author_notes_probe_backend',OLD/'run.py')
    run=importlib.util.module_from_spec(spec);spec.loader.exec_module(run)
    scenarios=[]
    for module in ('clock_grid','surfaces_off','independent_upper','clock_boundary'):
        values=json.loads((HERE/'configs'/(module+'.json')).read_text())['scenarios']
        if module=='clock_grid': scenarios.extend(values[i] for i in (0,8,9,18))
        else: scenarios.extend(values)
    config=run.normalize_config(dict(scenarios=scenarios,seeds=[231023],warmup_days=1.,observation_days=1.,block_days=1.,
              model_overrides=dict(prob_patient_sick=.03,initial_environment_contaminated=False)))
    hashes={}
    for job in run.build_jobs(config):
        row=run.simulate_pair(job,config,expected_hashes=run.source_hashes())
        if row['status']!='complete':raise RuntimeError(row)
        row.pop('elapsed_seconds',None)
        hashes[job['scenario']['id']]=hashlib.sha256(run.canonical(row).encode()).hexdigest()
    path=HERE/'probe_reference.json'
    if args.write_reference:
        if path.exists():raise ValueError('Reference already exists')
        path.write_text(json.dumps(hashes,indent=2,sort_keys=True)+'\n')
    elif json.loads(path.read_text())!=hashes:
        raise ValueError('Short numerical replay differs from local reference')
    print(json.dumps(dict(status='pass',scope='Nine one-day warmup/one-day observation checks only',hashes=hashes,versions=run.versions()),sort_keys=True))


if __name__=='__main__':main()
