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
cp -f ../../reports/results/ou_validation/ou_validation_macros.tex \
  "${DOC_DIR}/w3d-ou-validation-macros-generated.tex-part"
cp -f ../../reports/results/ou_validation/ou_validation_vertical.svg \
  "${DOC_DIR}/ou_validation_vertical.svg"
cp -f ../../reports/results/ou_robustness/ou_robustness_publication.tex \
  "${DOC_DIR}/w3d-ou-robustness-results-generated.tex-part"
cp -f ../../reports/results/ou_robustness/ou_robustness_macros.tex \
  "${DOC_DIR}/w3d-ou-robustness-macros-generated.tex-part"
cp -f ../../reports/results/ou_robustness/ou_robustness_sensitivity.svg \
  "${DOC_DIR}/ou_robustness_sensitivity.svg"
cp -f ../../reports/results/ou_robustness/ou_robustness_stress.svg \
  "${DOC_DIR}/ou_robustness_stress.svg"

# Keep the complete machine-generated validation tables as evidence, but make a
# compact publication view. Detailed covariance-control and per-scenario
# direction tables remain in the full generated study. The adaptive/fixed
# vertical table is also omitted from the publication because the same values
# are already presented in Fig. ou_mc_vertical.
python3 - <<'PY'
from pathlib import Path


def strip_table(text: str, label: str) -> str:
    pos = text.find(label)
    if pos < 0:
        raise RuntimeError(f"supplemental table label not found: {label}")
    start = text.rfind(r"\begin{table*}", 0, pos)
    end_token = r"\end{table*}"
    end = text.find(end_token, pos)
    if start < 0 or end < 0:
        raise RuntimeError(f"could not delimit supplemental table: {label}")
    end += len(end_token)
    return text[:start].rstrip() + "\n\n" + text[end:].lstrip()

src = Path("../../doc/kalman_ou_iii/w3d-ou-validation-results-generated.tex-part")
dst = Path("../../doc/kalman_ou_iii/w3d-ou-validation-results-publication.tex-part")
text = src.read_text(encoding="utf-8")
for label in (
    r"\label{tab:ou_mc_adaptation}",
    r"\label{tab:ou_mc_covsync}",
    r"\label{tab:ou_mc_direction}",
):
    text = strip_table(text, label)
dst.write_text(text, encoding="utf-8")
PY

# The robustness publication view omits the two full-width tables because both
# datasets are already shown by the sensitivity and degradation figures. The
# detailed-results document continues to consume the unfiltered generated file.
# The scalar macros quoted by the prose are not in either file: they are
# mirrored separately above and input from the preamble, because the Riccati
# analysis quotes them before these tables are read.
python3 - <<'PY'
from pathlib import Path


def strip_table(text: str, label: str) -> str:
    pos = text.find(label)
    if pos < 0:
        raise RuntimeError(f"supplemental table label not found: {label}")
    start = text.rfind(r"\begin{table*}", 0, pos)
    end_token = r"\end{table*}"
    end = text.find(end_token, pos)
    if start < 0 or end < 0:
        raise RuntimeError(f"could not delimit supplemental table: {label}")
    end += len(end_token)
    return text[:start].rstrip() + "\n\n" + text[end:].lstrip()

src = Path("../../doc/kalman_ou_iii/w3d-ou-robustness-results-generated.tex-part")
dst = Path("../../doc/kalman_ou_iii/w3d-ou-robustness-results-publication.tex-part")
text = src.read_text(encoding="utf-8")
for label in (
    r"\label{tab:ou_robustness_sensitivity}",
    r"\label{tab:ou_robustness_stress}",
):
    text = strip_table(text, label)
dst.write_text(text, encoding="utf-8")
PY

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
