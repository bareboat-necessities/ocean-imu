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


def main(data, case):
    evidence = ROOT / "reports/results/marine_mag"
    with (evidence / "validation-inputs.csv").open() as stream:
        source = next(r for r in csv.DictReader(stream) if r["case"] == case and r["scenario"] == "burst")
    path = data / source["file"]
    with path.open("rb") as stream:
        assert hashlib.file_digest(stream, "sha256").hexdigest() == source["sha256"]
    with (evidence / "validation.csv").open() as stream:
        rows = [r for r in csv.DictReader(stream) if r["scenario"] == "burst"
                and r["case"] == case and r["variant"] in ("0", "8") and r["seed"] == "101"]
    assert len(rows) == 2
    for row in rows:
        env = {k: v for k, v in os.environ.items() if not k.startswith(("W3D_", "OCEAN_IMU_SIM_CASE"))}
        env.update(OCEAN_IMU_SIM_CASE=case, W3D_SEED=row["seed"])
        output = subprocess.check_output([DRIVER, path, row["variant"], "burst", row["dip"]], env=env, text=True)
        line = next(line for line in output.splitlines() if line.startswith("VALIDATION_METRICS ") and "segment=window" in line)
        fields = dict(f.split("=", 1) for f in line.split()[1:])
        assert int(fields["samples"]) == 180000
        for name in ("yaw_rms_deg", "roll_rms_deg", "pitch_rms_deg", "disp_3d_rms_m"):
            assert math.isclose(float(fields[name]), float(row[name]), rel_tol=1e-5, abs_tol=1e-6), (
                case, row["variant"], name, fields[name], row[name])
        print(f'{case} V{row["variant"]} reproduced yaw RMS {fields["yaw_rms_deg"]} degrees')


if __name__ == "__main__":
    main(Path(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else "following")
