#!/bin/bash -e

# Execute the existing plotting script with publication-only substitutions.
# OU-III adaptation is driven by the wave-band period, so show that schedule
# rather than the acceleration-band carrier used by direction diagnostics.
python3 - <<'PY'
from pathlib import Path

path = Path("kalman_ou_iii-plots.py")
source = path.read_text(encoding="utf-8")

old_dir = '        ("dir_deg",        r"Dir (deg, axial)"),'
new_dir = '        (("dir_axis_deg" if "dir_axis_deg" in df.columns else "dir_deg"), r"Dir (deg, axial)"),'
if source.count(old_dir) != 1:
    raise RuntimeError("OU-III direction plot column anchor not found exactly once")
source = source.replace(old_dir, new_dir, 1)

old_panel = '        ("freq_tracker_hz", "Frequency (Hz)"),'
new_panel = '        ("wave_tuning_freq_hz", "Applied wave-band frequency (Hz)"),'
if source.count(old_panel) != 1:
    raise RuntimeError("OU-III tuner frequency panel anchor not found exactly once")
source = source.replace(old_panel, new_panel, 1)

anchor = '    # === Frequency / Tuner ===\n'
inject = '''    # Publication diagnostic for the deployed WaveBand tuning path.\n    tau_for_plot = pd.to_numeric(df["tau_applied"], errors="coerce").to_numpy()\n    df["wave_tuning_freq_hz"] = np.where(\n        np.isfinite(tau_for_plot) & (tau_for_plot > 0.0),\n        1.0 / (2.0 * tau_for_plot),\n        np.nan,\n    )\n\n'''
if source.count(anchor) != 1:
    raise RuntimeError("OU-III tuner section anchor not found exactly once")
source = source.replace(anchor, inject + anchor, 1)

namespace = {"__name__": "__main__", "__file__": str(path)}
exec(compile(source, str(path), "exec"), namespace)
PY

python3 ../spectrum/spectrum-plots.py

# Re-run comparison observers from the same source records.
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
cp -f ../../reports/results/ou_validation/ou_validation_macros.tex \
  "${DOC_DIR}/w3d-ou-validation-macros-generated.tex-part"
cp -f ../../reports/results/ou_validation/ou_validation_vertical.svg \
  "${DOC_DIR}/ou_validation_vertical.svg"

# Keep full archives for provenance, curate the validation article view, and
# publish sensitivity evidence only when its coupled rows were regenerated with
# the deployed SpectralMSE law.  The permanent OFAT table is not touched here.
python3 ../../tools/ou_publication_sync.py \
  --validation-dir ../../reports/results/ou_validation
python3 ../../tools/ou_publication_robustness_sync.py \
  --robustness-dir ../../reports/results/ou_robustness \
  --doc-dir "${DOC_DIR}"
rm -f "${DOC_DIR}/ou_validation_transition.svg"

# Regenerate the only transition evidence presented by the article: the
# deployed SpectralMSE low--high--low protocol.  Build the simulator here
# because the main workflow cleans it before entering the plotting stage.
find . -maxdepth 1 -type f -name 'wave_data_*.csv' \
  -exec cp -f {} ../../tests/kalman_ou_iii/ \;
(
  cd ../../tests/kalman_ou_iii
  make build
)
python3 ../../tools/ou_roundtrip_transition.py \
  --output-dir ../../reports/results/ou_rs_law
cp -f ../../reports/results/ou_rs_law/ou_rs_roundtrip_transition.svg \
  "${DOC_DIR}/ou_rs_roundtrip_transition.svg"
cp -f ../../reports/results/ou_rs_law/ou_rs_roundtrip_scores.tex \
  "${DOC_DIR}/w3d-roundtrip-transition-scores-generated.tex-part"
(
  cd ../../tests/kalman_ou_iii
  make clean
)

# Defensive publication gates.
test ! -e "${DOC_DIR}/ou_validation_transition.svg"
if test -e "${DOC_DIR}/w3d-ou-robustness-sensitivity-current-generated.tex-part"; then
  grep -q 'c^{6/7}' "${DOC_DIR}/w3d-ou-robustness-sensitivity-current-generated.tex-part"
  grep -q 'c^{41/14}' "${DOC_DIR}/w3d-ou-robustness-sensitivity-current-generated.tex-part"
  test -e "${DOC_DIR}/ou_robustness_sensitivity.svg"
fi
grep -q 'Rise crossfade' "${DOC_DIR}/w3d-roundtrip-transition-scores-generated.tex-part"
grep -q 'Fall crossfade' "${DOC_DIR}/w3d-roundtrip-transition-scores-generated.tex-part"
