"""Nine fixed short fixtures: deterministic replay, runtime and counter checks."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import time
from copy import deepcopy

HERE=Path(__file__).resolve().parent;ABM=HERE.parents[1]
spec=importlib.util.spec_from_file_location('_robustness_probe_runner',ABM/'research/full_icu_response_2026-10-08/run.py')
backend=importlib.util.module_from_spec(spec);sys.modules[spec.name]=backend;spec.loader.exec_module(backend)
from workflow_variant import WorkflowRobustnessModel
from freeze import verify,digest
from worker import runtime

def scientific(row):
    row=deepcopy(row);row.pop('elapsed_seconds',None)
    return hashlib.sha256(backend.canonical(row).encode()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--reference',type=Path)
    p.add_argument('--production',action='store_true');a=p.parse_args()
    if a.output.exists():raise ValueError('Probe output already exists')
    if a.production and a.reference is None:raise ValueError('Production probe requires a frozen local reference')
    frozen=verify();software=runtime(a.production)
    reference_sha=digest(a.reference) if a.reference else None
    reference=json.loads(a.reference.read_text()) if a.reference else None
    if reference is not None and reference.get('parent_freeze_identity')!=frozen['identity_sha256']:
        raise ValueError('Reference belongs to another freeze')
    raw=json.loads((HERE/'configs/probe.json').read_text());cfg=backend.normalize_config(raw)
    assert cfg['warmup_ticks']==cfg['observation_ticks']==2592
    started=time.monotonic();records=[]; hashes=[]
    for job in backend.build_jobs(cfg):
        row=backend.simulate_pair(job,cfg,model_factory=WorkflowRobustnessModel)
        if row['status']!='complete':raise ValueError(row)
        records.append(row);hashes.append(scientific(row))
    if reference is not None and hashes!=reference['scientific_record_sha256']:
        raise ValueError('Cross-runtime replay disagrees')
    if verify()!=frozen or runtime(a.production)!=software:raise ValueError('Probe sources/runtime changed')
    if a.reference and digest(a.reference)!=reference_sha:raise ValueError('Probe reference changed')
    obj=dict(status='passed',scope='Nine repeated short implementation fixtures; not study evidence',
        versions=software,scientific_record_sha256=hashes,
        parent_freeze_identity=frozen['identity_sha256'],reference_sha256=reference_sha,
        config_sha256=digest(HERE/'configs/probe.json'),probe_sha256=digest(Path(__file__)),production=a.production,
        elapsed_seconds=time.monotonic()-started,record_seconds=[r['elapsed_seconds'] for r in records],
        max_rss_native_units=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        source_sha256={**backend.source_hashes(),
            str((HERE/'workflow_variant.py').relative_to(ABM)):hashlib.sha256((HERE/'workflow_variant.py').read_bytes()).hexdigest()},
        reference_verified=a.reference is not None,records=records)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(obj,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps({k:obj[k] for k in ('status','elapsed_seconds','record_seconds','reference_verified')}))

if __name__=='__main__':main()
