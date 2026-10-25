"""Immutable source/design closure and separately bound selected-policy files."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ABM=HERE.parents[1]

def canonical(x):return json.dumps(x,sort_keys=True,separators=(",",":"),allow_nan=False)
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def identity(x):return hashlib.sha256(canonical({k:v for k,v in x.items() if k!='identity_sha256'}).encode()).hexdigest()
def write_new(path,value):
    data=json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+"\n"
    if path.exists():
        if path.read_text()!=data:raise ValueError("Existing immutable file differs: "+str(path))
    else:
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('x') as f:f.write(data)

def check(value):
    if value.get('study')!=HERE.name or value.get('identity_sha256')!=identity(value):
        raise ValueError('Invalid study freeze identity')
    for name,sha in value['files'].items():
        p=Path(name)
        if (p.is_absolute() or '..' in p.parts or
            any(q.is_symlink() for q in (ABM/p,*(ABM/p).parents)) or digest(ABM/p)!=sha):
            raise ValueError('Frozen file changed: '+name)
    return value

def verify():return check(json.loads((HERE/'freeze.json').read_text()))

def verify_config(path):
    base=verify();name=str(path.resolve().relative_to(ABM))
    if base['files'].get(name)==digest(path):return base
    selected=check(json.loads((HERE/'policy_freeze.json').read_text()))
    if selected.get('parent_freeze_identity')!=base['identity_sha256'] or selected['files'].get(name)!=digest(path):
        raise ValueError('Selected policy is not bound to parent study')
    return selected

def main():
    own=['PROTOCOL.md','EXECUTION_PLAN.md','build_design.py','design.json','seedbook.json',
         'configs/grid.json','configs/probe.json','workflow_variant.py','test_workflow_variant.py',
         'analyze.py','test_analysis.py','worker.py','freeze.py','probe.py','array.sbatch',
         'launch.py','submit_once.py','auxiliary.sbatch','status.py','test_launch.py']
    dependencies=['src/CPEmodel.py','src/agents.py','src/routes.py',
        'research/full_icu_calibration_2026-10-07/run.py',
        'research/full_icu_response_2026-10-08/run.py',
        'research/full_icu_response_2026-10-08/instrumented_model.py',
        'research/full_icu_response_2026-10-08/PROTOCOL.md',
        'research/rule_surface_followup_2026-10-11/analyze.py',
        'research/rule_surface_followup_2026-10-11/freeze.py',
        'research/full_icu_confirmation_2026-10-08/configs/holdout_all.json']
    paths=[HERE/x for x in own]+[ABM/x for x in dependencies]
    obj=dict(schema_version=1,study=HERE.name,files={str(x.relative_to(ABM)):digest(x) for x in sorted(paths)})
    obj['identity_sha256']=identity(obj)
    write_new(HERE/'freeze.json',obj);verify();print(obj['identity_sha256'])

if __name__=='__main__':main()
