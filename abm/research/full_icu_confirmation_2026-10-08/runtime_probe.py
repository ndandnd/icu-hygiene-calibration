#!/usr/bin/env python3
"""Short deterministic portability check, not a scientific experiment."""
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]
OLD = ABM / 'research/full_icu_response_2026-10-08'
spec = importlib.util.spec_from_file_location('_confirmation_probe', OLD / 'run.py')
run = importlib.util.module_from_spec(spec)
spec.loader.exec_module(run)

rows = []
for i, (rho, efficacy) in enumerate(((0., 1.), (.5, 1.), (0., .8))):
    raw = dict(warmup_days=1., observation_days=1., block_days=1,
               model_overrides=dict(prob_patient_sick=.03, initial_environment_contaminated=False),
               seeds=[10801], scenarios=[dict(id='runtime_validation_'+str(i),
               p=.2, h0=.9, h1=.95, rho=rho, e=efficacy,
               opportunities='zone_separate', environment_multiplier=1.)])
    config = run.normalize_config(raw)
    result = run.simulate_pair(run.build_jobs(config)[0], config)
    if result['status'] != 'complete':
        raise RuntimeError(result)
    result.pop('elapsed_seconds')
    rows.append(result)
rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print(json.dumps(dict(versions=run.versions(),
                     short_replay_sha256=hashlib.sha256(run.canonical(rows).encode()).hexdigest(),
                     max_rss_bytes=rss if sys.platform == 'darwin' else rss*1024,
                     source_sha256=run.source_hashes(),
                     scope='Three one-day warmup/one-day observation implementation checks only.'),
                 sort_keys=True))
