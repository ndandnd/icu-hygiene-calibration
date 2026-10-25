#!/usr/bin/env python3
"""Immutable first-phase closure and separately bound derived clock selection."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ABM=HERE.parents[1]
def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def identity(value):
    return hashlib.sha256(canonical({k:v for k,v in value.items() if k!='identity_sha256'}).encode()).hexdigest()
def check(value,root=ABM):
    if value.get('study')!=HERE.name or value.get('identity_sha256')!=identity(value):
        raise ValueError('Invalid freeze identity')
    for name,sha in value['files'].items():
        p=Path(name)
        if p.is_absolute() or '..' in p.parts or (root/p).is_symlink() or digest(root/p)!=sha:
            raise ValueError('Frozen file changed: '+name)
    return value
def verify(): return check(json.loads((HERE/'freeze.json').read_text()))
def verify_config(path):
    base=verify();name=str(path.resolve().relative_to(ABM))
    if base['files'].get(name)==digest(path): return base
    derived=check(json.loads((HERE/'clock_policy_freeze.json').read_text()))
    if derived.get('parent_freeze_identity')!=base['identity_sha256'] or derived['files'].get(name)!=digest(path):
        raise ValueError('Derived config is not bound to this study')
    return derived
def main():
    own=['PROTOCOL.md','build_design.py','freeze.py','array_worker.py','array_resume.sbatch',
         'retry_state.py','probe.py','probe_reference.json','select_clock.py','analyze.py','test_analysis.py',
         'design.json','seedbook.json','output_scopes.json','old_input_index.json']
    own += ['configs/'+p.name for p in sorted((HERE/'configs').glob('*.json'))]
    dependencies=['src/CPEmodel.py','src/agents.py','src/routes.py',
        'research/full_icu_calibration_2026-10-07/run.py',
        'research/full_icu_response_2026-10-08/run.py',
        'research/full_icu_response_2026-10-08/array_worker.py',
        'research/full_icu_response_2026-10-08/instrumented_model.py',
        'research/full_icu_response_2026-10-08/PROTOCOL.md',
        'research/full_icu_confirmation_2026-10-08/configs/holdout_all.json',
        'research/full_icu_confirmation_2026-10-08/seedbook.json',
        'research/reviewer_comments_2026-10-11/saved_data/build_old_index.py']
    paths=[HERE/p for p in own]+[ABM/p for p in dependencies]
    value=dict(schema_version=1,study=HERE.name,files={str(p.relative_to(ABM)):digest(p) for p in sorted(paths)})
    value['identity_sha256']=identity(value)
    text=json.dumps(value,indent=2,sort_keys=True)+'\n';target=HERE/'freeze.json'
    if target.exists() and target.read_text()!=text: raise ValueError('Existing freeze differs')
    target.write_text(text);verify();print(value['identity_sha256'])
if __name__=='__main__': main()
