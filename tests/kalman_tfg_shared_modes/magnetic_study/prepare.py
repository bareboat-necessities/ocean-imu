"""Construct a research-only include tree. Patch checked copies, never originals."""
import shutil
import subprocess
import sys
from pathlib import Path

here=Path(__file__).resolve().parent
root=here.parents[2]
dest=Path(sys.argv[1]).resolve()
patch=here/'instrumentation.patch'
paths=[line.split(' b/',1)[1] for line in patch.read_text().splitlines() if line.startswith('diff --git')]
if not paths or dest in (root,here):
    raise ValueError('invalid overlay destination or patch')
dest.mkdir(parents=True,exist_ok=True)
for relative in paths:
    target=dest/relative
    if target.is_symlink() or not target.resolve().is_relative_to(dest):
        raise ValueError(f'unsafe overlay target {relative}')
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(root/relative,target)
# Quoted same-directory includes and ../ in the test adapter need their siblings.
for parent in sorted({str(Path(p).parent) for p in paths}):
    base=root/parent
    for item in base.iterdir():
        target=dest/parent/item.name
        if not target.exists() and not target.is_symlink():
            target.symlink_to(item,target_is_directory=item.is_dir())
for item in (root/'tests').iterdir():
    target=dest/'tests'/item.name
    if item.is_dir() and not target.exists():
        target.symlink_to(item,target_is_directory=True)
subprocess.run(['patch','--batch','--fuzz=0','--no-backup-if-mismatch','-s','-p1','-i',str(patch)],cwd=dest,check=True)
