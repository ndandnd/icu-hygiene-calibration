"""Freeze this one execution repair after independently reviewed source changes."""
import ast
import difflib
import hashlib
import importlib.util
import json
import os
from pathlib import Path

HERE=Path(__file__).resolve().parent
STUDY=HERE.parent
ABM=STUDY.parents[1]
spec=importlib.util.spec_from_file_location('_repair_freeze',STUDY/'freeze.py')
freeze=importlib.util.module_from_spec(spec);spec.loader.exec_module(freeze)
original_root=HERE/'initial_source'
original_path=original_root/'research'/STUDY.name/'freeze.json'
old=json.loads(original_path.read_text())
assert old['identity_sha256']==freeze.identity(old)
allowed={'launch.py','worker.py','analyze.py','test_launch.py','test_analysis.py'}
allowed={str((STUDY/p).relative_to(ABM)) for p in allowed}
changed={};patch=[]
for name,sha in old['files'].items():
    original=original_root/name;current=ABM/name
    assert freeze.digest(original)==sha
    current_sha=freeze.digest(current)
    if current_sha!=sha:
        assert name in allowed, 'Undeclared source change: '+name
        changed[name]={'original':sha,'amended':current_sha}
        patch.extend(difflib.unified_diff(original.read_text().splitlines(True),
            current.read_text().splitlines(True),fromfile='original/'+name,tofile='amended/'+name))
assert changed
functions=('reconcile','profile','jobs_for','validate_config','validate_records','validate_software',
    'estimates','interval','finite','bootstrap_totals','choose','select_values','select',
    'summarize_values','analyze_bank')
def bodies(path):
    return {node.name:ast.dump(node,include_attributes=False) for node in ast.parse(path.read_text()).body
            if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef)) and node.name in functions}
before=bodies(original_root/'research'/STUDY.name/'analyze.py')
after=bodies(STUDY/'analyze.py')
assert set(before)==set(functions) and before==after, 'Scientific analysis functions changed'
science={name:hashlib.sha256(body.encode()).hexdigest() for name,body in before.items()}
patch_path=HERE/'execution_diff.patch'
patch_text=''.join(patch)
if patch_path.exists():assert patch_path.read_text()==patch_text
else:patch_path.write_text(patch_text)
details=dict(original_freeze_identity=old['identity_sha256'],changed_files=changed,
    unchanged_scientific_analysis_ast_sha256=science,scientific_records_before_repair=0,
    replacement_probe_limit=1,extra_probe_allocation_seconds=300,
    full_scientific_record_pair_ceiling=3072,maximum_allocated_cpu_seconds=769*3600+10*60)
freeze.write_new(HERE/'SOURCE_RECONCILIATION.json',details)
files={name:freeze.digest(ABM/name) for name in old['files']}
for p in sorted(HERE.rglob('*')):
    if p.is_file() and '__pycache__' not in p.parts and p.name!='amended_freeze.json':
        files[str(p.relative_to(ABM))]=freeze.digest(p)
obj=dict(schema_version=1,study=STUDY.name,files=files,execution_amendment=details)
obj['identity_sha256']=freeze.identity(obj)
freeze.write_new(HERE/'amended_freeze.json',obj)
current=json.loads((STUDY/'freeze.json').read_text())
assert current in (old,obj), 'Unexpected existing freeze'
if current!=obj:
    temporary=STUDY/'freeze.execution01.tmp'
    with temporary.open('x') as f:
        f.write(json.dumps(obj,indent=2,sort_keys=True)+'\n');f.flush();os.fsync(f.fileno())
    os.replace(temporary,STUDY/'freeze.json')
freeze.verify()
print(json.dumps({'identity_sha256':obj['identity_sha256'],'changed_files':list(changed),
    'scientific_functions_unchanged':len(science),'files':len(files)}))
