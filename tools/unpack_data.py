#!/usr/bin/env python3
"""Verify all compressed bytes and member hashes before safe extraction."""
import argparse, hashlib, io, json, tarfile
from pathlib import Path, PurePosixPath
ROOT=Path(__file__).resolve().parents[1]
def digest(b): return hashlib.sha256(b).hexdigest()
p=argparse.ArgumentParser();p.add_argument('--study',action='append');p.add_argument('--verify-only',action='store_true');a=p.parse_args()
index=json.loads((ROOT/'DATA_INDEX.json').read_text())
selected=[s for s in index['studies'] if not a.study or s['study'] in a.study]
if a.study and set(a.study)-{s['study']for s in selected}:raise SystemExit('Unknown study')
for study in selected:
 expected={m['path']:m for m in study['members']};seen=set()
 for chunk in study['archives']:
  raw=(ROOT/chunk['path']).read_bytes();assert digest(raw)==chunk['sha256'],chunk['path']
  with tarfile.open(fileobj=io.BytesIO(raw),mode='r:gz')as t:
   for m in t:
    n=PurePosixPath(m.name)
    assert m.isfile() and not n.is_absolute() and '..'not in n.parts and m.name in expected,m.name
    assert m.name not in seen,m.name;seen.add(m.name)
    b=t.extractfile(m).read();assert digest(b)==expected[m.name]['sha256'] and len(b)==expected[m.name]['bytes'],m.name
    if not a.verify_only:
     out=ROOT/m.name;assert out.resolve().is_relative_to(ROOT)
     if out.exists():assert out.is_file()and digest(out.read_bytes())==digest(b),'Existing file differs: '+m.name
     else:out.parent.mkdir(parents=True,exist_ok=True);out.write_bytes(b)
 assert seen==set(expected),study['study']
 print(study['study'],len(seen),'members verified'+(' and extracted'if not a.verify_only else ''))
