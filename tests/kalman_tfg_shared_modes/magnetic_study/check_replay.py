"""Reproduce two held-out navigation scores from a released physical record."""
import csv
import hashlib
import math
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DRIVER = Path(__file__).resolve().parent / "build/magnetic-study"


def check_run(data, case, *, sustained=False):
    evidence = ROOT / "reports/results/marine_mag"
    # The 4 m PM-Stokes physical record is identical for both interventions.
    with (evidence / "validation-inputs.csv").open() as stream:
        source = next(r for r in csv.DictReader(stream) if r["case"] == case and r["scenario"] == "burst")
    path = data / source["file"]
    with path.open("rb") as stream:
        assert hashlib.file_digest(stream, "sha256").hexdigest() == source["sha256"]
    if sustained:
        with (evidence / "yaw-weight-validation.csv").open() as stream:
            rows = [r for r in csv.DictReader(stream) if float(r["height"]) == 4
                    and r["case"] == case and r["seed"] == "421"]
    else:
        with (evidence / "validation.csv").open() as stream:
            rows = [r for r in csv.DictReader(stream) if r["scenario"] == "burst"
                and r["case"] == case and r["variant"] in ("0", "8") and r["seed"] == "101"]
    assert len(rows) == 2
    for row in rows:
        env = {k: v for k, v in os.environ.items() if not k.startswith(("W3D_", "OCEAN_IMU_SIM_CASE"))}
        env.update(OCEAN_IMU_SIM_CASE=case, W3D_SEED=row["seed"])
        driver = DRIVER.with_name("magnetic-heading-study") if sustained else DRIVER
        variant = row.get("variant", "0")
        label = f"sigma x{row['sigma_factor']}" if sustained else f"V{variant}"
        if sustained:
            env["MARINE_MAG_SIGMA_FACTOR"] = row["sigma_factor"]
        output = subprocess.check_output([driver, path, variant, row["scenario"], row["dip"]], env=env, text=True)
        line = next(line for line in output.splitlines() if line.startswith("VALIDATION_METRICS ") and "segment=window" in line)
        fields = dict(f.split("=", 1) for f in line.split()[1:])
        assert int(fields["samples"]) == 180000
        for name in ("yaw_rms_deg", "roll_rms_deg", "pitch_rms_deg", "disp_3d_rms_m"):
            assert math.isclose(float(fields[name]), float(row[name]), rel_tol=1e-5, abs_tol=1e-6), (
                case, label, name, fields[name], row[name])
        print(f'{case} {label} reproduced yaw RMS {fields["yaw_rms_deg"]} degrees')


if __name__ == "__main__":
    data = Path(sys.argv[1])
    direction = sys.argv[2] if len(sys.argv) > 2 else "following"
    check_run(data, direction)
    check_run(data, direction, sustained=True)
