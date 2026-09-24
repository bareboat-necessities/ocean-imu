#!/usr/bin/env python3
from pathlib import Path
import csv

SCRIPT_DIR = Path(__file__).resolve().parent
TEST_DIR = SCRIPT_DIR.parent / "../tests/imu_calibrate"
OUT_PATH = SCRIPT_DIR / "calibrate_imu_quality.pgf"
SAMPLES_CSV = TEST_DIR / "calibrate_imu_test_output.csv"
SUMMARY_CSV = TEST_DIR / "calibrate_imu_test_summary.csv"
CAMPAIGN_CSV = TEST_DIR / "calibrate_accel_campaign_summary.csv"
CAMPAIGN_OUT = SCRIPT_DIR / "calibrate_accel_campaign_table.pgf"

SCENARIO_LABELS = {
    "typical": "typical",
    "adverse": "adverse",
    "cold_start_thermal": "cold start",
    "thermal_scale_mismatch": "scale drift $10^{-4}$/K",
    "thermal_scale_stress": "scale drift $3{\\times}10^{-4}$/K",
    "accel_only": "accel only",
    "local_g": "local $g = 9.780$",
    "temp_nan": "no temperature",
    "temp_confounded": "confounded T",
}
METHOD_LABELS = {"OLD": "deployed", "RICH+OLD": "new capture, old fit", "RICH+NEW": "new"}


def load_rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def coords(rows, x_key, y_key):
    return "\n".join(f"({float(r[x_key]):.3f},{float(r[y_key]):.6f})" for r in rows)


def gate_value(summary_rows, sensor, metric):
    for row in summary_rows:
        if row["sensor"] == sensor and row["metric"] == metric:
            return float(row["gate"]), float(row["value"])
    raise RuntimeError(f"missing gate for {sensor}/{metric}")


def main() -> None:
    sample_rows = load_rows(SAMPLES_CSV)
    summary_rows = load_rows(SUMMARY_CSV)

    accel_rows = [r for r in sample_rows if r["sensor"] == "accel"]
    mag_rows = [r for r in sample_rows if r["sensor"] == "mag"]
    gyro_rows = [r for r in sample_rows if r["sensor"] == "gyro"]

    accel_gate, accel_val = gate_value(summary_rows, "accel", "norm_rms_mps2")
    mag_gate, mag_val = gate_value(summary_rows, "mag", "norm_rms_uT")
    gyro_gate, _ = gate_value(summary_rows, "gyro", "bias_fit_rms_rads")

    tex = f"""\\begin{{tikzpicture}}
\\begin{{groupplot}}[
  group style={{group size=1 by 3, vertical sep=1.4cm}},
  width=14.5cm,
  height=3.9cm,
  grid=both,
  grid style={{line width=.1pt, draw=gray!25}},
  major grid style={{line width=.2pt, draw=gray!40}},
  xlabel={{}},
  ylabel near ticks,
]
\\nextgroupplot[
  title={{Accelerometer calibration error norm (RMS={accel_val:.4f}, gate={accel_gate:.4f})}},
  ylabel={{$|\\|a_{{cal}}\\|-g|$ [m/s$^2$]}},
]
\\addplot[line width=0.8pt, color=blue!70!black] coordinates {{
{coords(accel_rows, 'sample', 'error_norm')}
}};
\\addplot[dashed, line width=0.8pt, color=red!70!black] coordinates {{(0,{accel_gate:.6f}) ({len(accel_rows)-1},{accel_gate:.6f})}};
\\addplot[dashed, line width=0.8pt, color=red!70!black] coordinates {{(0,{-accel_gate:.6f}) ({len(accel_rows)-1},{-accel_gate:.6f})}};

\\nextgroupplot[
  title={{Magnetometer calibrated norm error (RMS={mag_val:.4f}, gate={mag_gate:.4f})}},
  ylabel={{$\\|m_{{cal}}\\|-B$ [$\\mu$T]}},
]
\\addplot[line width=0.8pt, color=teal!70!black] coordinates {{
{coords(mag_rows, 'sample', 'error_norm')}
}};
\\addplot[dashed, line width=0.8pt, color=red!70!black] coordinates {{(0,{mag_gate:.6f}) ({len(mag_rows)-1},{mag_gate:.6f})}};
\\addplot[dashed, line width=0.8pt, color=red!70!black] coordinates {{(0,{-mag_gate:.6f}) ({len(mag_rows)-1},{-mag_gate:.6f})}};

\\nextgroupplot[
  title={{Gyroscope residual norm after bias compensation}},
  ylabel={{$\\|\\omega_{{cal}}\\|$ [rad/s]}},
]
\\addplot[line width=0.8pt, color=purple!70!black] coordinates {{
{coords(gyro_rows, 'sample', 'error_norm')}
}};
\\addplot[dashed, line width=0.8pt, color=red!70!black] coordinates {{(0,{gyro_gate:.6f}) ({len(gyro_rows)-1},{gyro_gate:.6f})}};
\\end{{groupplot}}
\\end{{tikzpicture}}
"""
    OUT_PATH.write_text(tex)
    print(f"saved {OUT_PATH}")
    write_campaign_table()


def fmt(value: str, scale: float = 1.0, digits: int = 1) -> str:
    try:
        v = float(value)
    except ValueError:
        return "--"
    if v != v:  # NaN
        return "--"
    return f"{v * scale:.{digits}f}"


def write_campaign_table() -> None:
    """Accelerometer campaign summary as a LaTeX table (mm/s^2, s, %)."""
    order = {m: i for i, m in enumerate(METHOD_LABELS)}
    scenarios = []
    for r in load_rows(CAMPAIGN_CSV):
        if r["scenario"] not in scenarios:
            scenarios.append(r["scenario"])
    rows = sorted(load_rows(CAMPAIGN_CSV),
                  key=lambda r: (scenarios.index(r["scenario"]), order.get(r["method"], 99)))
    lines = [
        "\\begin{tabular}{@{}llrrrrrrrr@{}}",
        "\\toprule",
        "Scenario & Method & ok\\,\\% & \\multicolumn{2}{c}{bias err.\\ [mm/s$^2$]} & cross & "
        "vector & ident.\\,\\% & learned\\,\\% & time \\\\",
        " & & & median & p90 & [$10^{-4}$] & [mm/s$^2$] & & & [s] \\\\",
        "\\midrule",
    ]
    last = None
    for r in rows:
        scn = SCENARIO_LABELS.get(r["scenario"], r["scenario"])
        if last is not None and scn != last:
            lines.append("\\addlinespace")
        lines.append(
            f"{scn if scn != last else ''} & {METHOD_LABELS.get(r['method'], r['method'])} & "
            f"{fmt(r['success_rate'], 100, 0)} & {fmt(r['bias_err_median'], 1000)} & {fmt(r['bias_err_p90'], 1000)} & "
            f"{fmt(r['cross_err_median'], 1e4)} & {fmt(r['vec_rms_median'], 1000)} & "
            f"{fmt(r['identity_accepted_rate'], 100, 0)} & {fmt(r['learned_rate'], 100, 0)} & "
            f"{fmt(r['time_s_median'], 1, 0)} \\\\"
        )
        last = scn
    lines += ["\\bottomrule", "\\end{tabular}", ""]
    CAMPAIGN_OUT.write_text("\n".join(lines))
    print(f"saved {CAMPAIGN_OUT}")


if __name__ == "__main__":
    main()