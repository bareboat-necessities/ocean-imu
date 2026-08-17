#!/bin/bash -e

# Execute the existing plotting script with two publication-only compatibility
# substitutions.  The simulator still records the acceleration-band carrier for
# direction diagnostics, but OU-III adaptation is driven by the wave-band period.
# Plot the smoothed wave-band schedule actually applied to the filter instead:
# with deployed c_tau=1, tau_applied = Tz_applied/2 and therefore
# f_wave_applied = 1/(2*tau_applied).  This avoids changing simulator telemetry
# (and therefore the evidence contract) solely for a publication figure.
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
inject = '''    # Publication diagnostic for the deployed WaveBand tuning path.\n    # The acceleration carrier in freq_tracker_hz belongs to wave direction and\n    # is intentionally not shown in the adaptation figure.\n    tau_for_plot = pd.to_numeric(df["tau_applied"], errors="coerce").to_numpy()\n    df["wave_tuning_freq_hz"] = np.where(\n        np.isfinite(tau_for_plot) & (tau_for_plot > 0.0),\n        1.0 / (2.0 * tau_for_plot),\n        np.nan,\n    )\n\n'''
if source.count(anchor) != 1:
    raise RuntimeError("OU-III tuner section anchor not found exactly once")
source = source.replace(anchor, inject + anchor, 1)

namespace = {"__name__": "__main__", "__file__": str(path)}
exec(compile(source, str(path), "exec"), namespace)
PY

# Directional frequency--direction spectra used by the evaluation-methodology
# section are generated from the same committed source records.
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
# direction tables remain in the full generated study.
python3 - <<'PY'
from pathlib import Path
import re


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


def replace_table_caption(text: str, label: str, caption: str) -> str:
    pos = text.find(label)
    if pos < 0:
        raise RuntimeError(f"table label not found: {label}")
    start = text.rfind(r"\begin{table*}", 0, pos)
    cap = text.find(r"\caption{", start, pos)
    if start < 0 or cap < 0:
        raise RuntimeError(f"caption not found for table: {label}")
    body_start = cap + len(r"\caption{")
    depth = 1
    i = body_start
    while i < pos and depth:
        if text[i] == "{" and (i == 0 or text[i - 1] != "\\"):
            depth += 1
        elif text[i] == "}" and (i == 0 or text[i - 1] != "\\"):
            depth -= 1
        i += 1
    if depth != 0:
        raise RuntimeError(f"unbalanced caption for table: {label}")
    return text[:cap] + r"\caption{" + caption + "}" + text[i:]


def compact_axes_table(text: str) -> str:
    label = r"\label{tab:ou_mc_axes}"
    pos = text.find(label)
    if pos < 0:
        raise RuntimeError("3-D axis table label not found")
    start = text.rfind(r"\begin{table*}", 0, pos)
    end_token = r"\end{table*}"
    end = text.find(end_token, pos)
    if start < 0 or end < 0:
        raise RuntimeError("could not delimit 3-D axis table")
    end += len(end_token)
    block = text[start:end]

    spacing = r"\setlength{\tabcolsep}{3.0pt}"
    if block.count(spacing) != 1:
        raise RuntimeError("3-D axis table spacing anchor not found exactly once")
    block = block.replace(spacing, r"\setlength{\tabcolsep}{1.6pt}", 1)

    block, scenario_count = re.subn(r"\$H_s=([0-9.]+)\$ m", r"\1 m", block)
    if scenario_count != 4:
        raise RuntimeError(
            f"expected four stationary H_s labels in 3-D axis table, found {scenario_count}"
        )

    block = re.sub(
        r"([0-9]+\.[0-9]+) \$\\pm\$ ([0-9]+\.[0-9]+)",
        r"$\1{\\pm}\2$",
        block,
    )
    return text[:start] + block + text[end:]


src = Path("../../doc/kalman_ou_iii/w3d-ou-validation-results-generated.tex-part")
dst = Path("../../doc/kalman_ou_iii/w3d-ou-validation-results-publication.tex-part")
text = src.read_text(encoding="utf-8")
for label in (
    r"\label{tab:ou_mc_adaptation}",
    r"\label{tab:ou_mc_covsync}",
    r"\label{tab:ou_mc_direction}",
):
    text = strip_table(text, label)
text = compact_axes_table(text)
for label, caption in (
    (r"\label{tab:ou_mc_family}",
     r"Ten-seed paired OU-family comparison over the final \SI{900}{s}."),
    (r"\label{tab:ou_mc_axes}",
     r"Ten-seed adaptive displacement RMS by axis over the final \SI{900}{s}."),
    (r"\label{tab:ou_mc_channels}",
     r"OU--III adaptation-channel ablation over the final \SI{900}{s}."),
    (r"\label{tab:ou_transition_segments}",
     r"Controlled transition scored by interval (Adaptive mode)."),
    (r"\label{tab:ou_mc_pmstokes}",
     r"Ten-seed paired OU-family comparison on PM--Stokes seas."),
):
    text = replace_table_caption(text, label, caption)
dst.write_text(text, encoding="utf-8")
PY

# The robustness publication view omits the two full-width tables because both
# datasets are already shown by the sensitivity and degradation figures.
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