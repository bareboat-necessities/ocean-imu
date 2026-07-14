#!/bin/bash -e

# The corrected telemetry names the axial field explicitly. Execute the existing
# plotting script with an in-memory compatibility substitution so both new
# `dir_axis_deg` CSVs and older `dir_deg` CSVs remain plottable.
python3 - <<'PY'
from pathlib import Path

path = Path("kalman_ou_iii-plots.py")
source = path.read_text(encoding="utf-8")
old = '        ("dir_deg",        r"Dir (deg, axial)"),'
new = '        (("dir_axis_deg" if "dir_axis_deg" in df.columns else "dir_deg"), r"Dir (deg, axial)"),'
if source.count(old) != 1:
    raise RuntimeError("OU-III direction plot column anchor not found exactly once")
namespace = {"__name__": "__main__", "__file__": str(path)}
exec(compile(source.replace(old, new, 1), str(path), "exec"), namespace)
PY

python3 ../spectrum/spectrum-plots.py  # .tex needs some spectrum .pgf later
