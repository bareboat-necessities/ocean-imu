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

# Run the existing VVR-only time-varying-gain NLO on the same source records
# used by this OU-III job.  Keeping both runners in one CI workspace guarantees
# paired timestamps, wave realization, and sensor-error realization.
NLO_TEST_DIR="../../tests/nlo"
find . -maxdepth 1 -type f -name '*.csv' ! -name '*_fusion_ou3.csv' ! -name '*_tvg_nlo_nomag_nognss.csv' \
  -exec cp -f {} "${NLO_TEST_DIR}/" \;
(
  cd "${NLO_TEST_DIR}"
  make build
  chmod +x ./*.sh
  ./run_tests.sh
  cp -f ./*_tvg_nlo_nomag_nognss.csv ../../plots/kalman_ou_iii/
  make clean
)

python3 baseline-comparison.py

# The table is generated from CSVs rather than transcribed by hand.  Copy the
# generated include next to the article and insert the stable section include
# immediately before the hardware-validation section in the CI workspace.
DOC_DIR="../../doc/kalman_ou_iii"
cp -f w3d-baseline-results-generated.tex-part "${DOC_DIR}/"
cp -f w3d_ou3_vs_tvg_nlo_jonswap_medium.pgf "${DOC_DIR}/"
cp -f w3d_ou3_vs_tvg_nlo_jonswap_medium.svg "${DOC_DIR}/"

python3 - <<'PY'
from pathlib import Path

path = Path("../../doc/kalman_ou_iii/w3d-sim-charts.tex-part")
source = path.read_text(encoding="utf-8")
include = r"\input{w3d-baseline-comparison.tex-part}"
anchor = r"\section{Real-Hardware Validation Platform}"
if include not in source:
    if source.count(anchor) != 1:
        raise RuntimeError("hardware-section insertion anchor not found exactly once")
    source = source.replace(anchor, include + "\n\n" + anchor, 1)
    path.write_text(source, encoding="utf-8")
PY
