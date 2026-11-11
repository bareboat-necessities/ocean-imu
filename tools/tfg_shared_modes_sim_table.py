#!/usr/bin/env python3
"""Regenerate the TFG shared-mode deterministic replay table using the common formatter."""
import json
from pathlib import Path
from carma_sim_table import main


def write_gate_table() -> None:
    root = Path(__file__).resolve().parents[1]
    baseline = json.loads((root / "reports/results/tfg_shared_modes/deterministic_regression.json").read_text())
    names = {
        "err_limit_percent_z_jonswap": "JONSWAP heave [\\% $H_s$]",
        "err_limit_percent_z_pmstokes": "PM--Stokes heave [\\% $H_s$]",
        "err_limit_yaw_deg": "Heading [deg]", "err_limit_roll_deg": "Roll [deg]",
        "err_limit_pitch_deg": "Pitch [deg]",
        "err_limit_percent_3d_jonswap": "JONSWAP 3-D position [\\%]",
        "err_limit_percent_3d_pmstokes": "PM--Stokes 3-D position [\\%]",
        "acc_z_bias_percent": "Accelerometer Z bias [\\%]",
        "bias_3d_percent": "Accelerometer 3-D bias [\\%]",
        "gyro_bias_3d_percent": "Gyroscope 3-D bias [\\%]"}
    rows = [r"% Generated from the locked deterministic regression baseline.",
            r"\begin{table}[!t]\centering\footnotesize",
            r"\caption{TFG shared-mode deterministic quality gates. Bias and 3-D position percentages use the reference maxima defined in the text.}",
            r"\label{tab:tfg_shared_gates}",
            r"\begin{tabular}{lrr}\toprule Quantity&Measured worst&Limit\\\midrule"]
    for key, name in names.items():
        rows.append(f"{name}&{baseline['measured_worst'][key]['value']:.4f}&{baseline['limits'][key]:.4f}" + r"\\")
    rows.extend([r"\bottomrule\end{tabular}\end{table}", ""])
    (root / "doc/kalman_tfg_shared_modes/tfg-shared-gates-generated.tex-part").write_text("\n".join(rows))

if __name__ == "__main__":
    status = main(
        table={"tool": "tools/tfg_shared_modes_sim_table.py",
               "variant": "TFG shared modes (estimated coupling with confidence)",
               "caption": "TFG shared-mode CARMA simulation results (last \\SI{900}{s} RMS window).",
               "label": "tab:tfg_shared_sim_results"},
        binary_default="tests/kalman_tfg_shared_modes/kalman_tfg_shared_modes-sim",
        data_default="plots/kalman_tfg_shared_modes",
        output_default="doc/kalman_tfg_shared_modes/tfg-shared-sim-results-generated.tex-part")
    if status == 0:
        write_gate_table()
    raise SystemExit(status)
