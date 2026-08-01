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

python3 ../spectrum/spectrum-plots.py

# Re-run all comparison observers in this workspace from the same source records.
# This guarantees identical timestamps, wave realization, and sensor-error realization.
run_comparison() {
  local test_dir="$1"
  local output_glob="$2"
  find . -maxdepth 1 -type f -name '*.csv' \
    ! -name '*_fusion_ou3.csv' \
    ! -name '*_fusion_ou2.csv' \
    ! -name '*_nonkalman_fusion.csv' \
    ! -name '*_tvg_nlo_nomag_nognss.csv' \
    -exec cp -f {} "${test_dir}/" \;
  (
    cd "${test_dir}"
    make build
    chmod +x ./*.sh
    ./run_tests.sh
    cp -f ${output_glob} ../../plots/kalman_ou_iii/
    make clean
  )
}

run_comparison "../../tests/nlo" "*_tvg_nlo_nomag_nognss.csv"
run_comparison "../../tests/kalman_ou_ii" "*_fusion_ou2.csv"
run_comparison "../../tests/pii_observer" "*_nonkalman_fusion.csv"

python3 baseline-comparison.py

DOC_DIR="../../doc/kalman_ou_iii"
cp -f w3d-baseline-results-generated.tex-part "${DOC_DIR}/"
cp -f w3d_multi_observer_jonswap_medium.pgf "${DOC_DIR}/"
cp -f w3d_multi_observer_jonswap_medium.svg "${DOC_DIR}/"
cp -f ../../reports/results/ou_validation/ou_validation_publication.tex \
  "${DOC_DIR}/w3d-ou-validation-results-generated.tex-part"
cp -f ../../reports/results/ou_validation/ou_validation_vertical.svg \
  "${DOC_DIR}/ou_validation_vertical.svg"
cp -f ../../reports/results/ou_robustness/ou_robustness_publication.tex \
  "${DOC_DIR}/w3d-ou-robustness-results-generated.tex-part"
cp -f ../../reports/results/ou_robustness/ou_robustness_sensitivity.svg \
  "${DOC_DIR}/ou_robustness_sensitivity.svg"
cp -f ../../reports/results/ou_robustness/ou_robustness_stress.svg \
  "${DOC_DIR}/ou_robustness_stress.svg"

python3 - <<'PY'
from pathlib import Path

sim_path = Path("../../doc/kalman_ou_iii/w3d-sim-charts.tex-part")
source = sim_path.read_text(encoding="utf-8")
# Use one manuscript name for the nonlinear Pierson--Moskowitz cases.
source = source.replace("PM+Stokes", "PM--Stokes")
include = r"\input{w3d-baseline-comparison.tex-part}"
anchor = r"\section{Real-Hardware Validation Platform}"
if include not in source:
    if source.count(anchor) != 1:
        raise RuntimeError("hardware-section insertion anchor not found exactly once")
    source = source.replace(anchor, include + "\n\n" + anchor, 1)
sim_path.write_text(source, encoding="utf-8")

main_path = Path("../../doc/kalman_ou_iii/kalman_ou-w3d.tex")
main = main_path.read_text(encoding="utf-8")
old_bib = r"\bibliography{w3d,w3d-iss}"
new_bib = r"\bibliography{w3d,w3d-iss,w3d-baselines}"
if new_bib not in main:
    if main.count(old_bib) != 1:
        raise RuntimeError("article bibliography anchor not found exactly once")
    main_path.write_text(main.replace(old_bib, new_bib, 1), encoding="utf-8")
PY
