#!/usr/bin/env python3
"""Bind all delivered files; run after adding the final paper artifacts."""
import hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
rows=[]
for p in sorted(root.rglob('*')):
 if '.git'in p.relative_to(root).parts:continue
 if p.is_symlink():raise SystemExit('Symlinks are not accepted: '+str(p))
 if p.is_file()and p.name!='SHA256SUMS':
  if p.stat().st_size>=90_000_000:raise SystemExit('File exceeds release cap: '+str(p))
  rows.append(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(root).as_posix())
(root/'SHA256SUMS').write_text('\n'.join(rows)+'\n')
print(len(rows),'regular files bound in SHA256SUMS')
