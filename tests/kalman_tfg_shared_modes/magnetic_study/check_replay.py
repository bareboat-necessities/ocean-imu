"""Reproduce held-out pairs from a pinned physical record in each sea direction."""
import csv
import hashlib
import math
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DRIVER = Path(__file__).resolve().parent / "build/magnetic-study"


def study_input(data, case, height):
    evidence = ROOT / "reports/results/marine_mag"
    # Both panels use the release record at this wave height and direction.
    with (evidence / "validation-inputs.csv").open() as stream:
        source = next(r for r in csv.DictReader(stream) if r["case"] == case and float(r["height"]) == height)
    path = data / source["file"]
    with path.open("rb") as stream:
        assert hashlib.file_digest(stream, "sha256").hexdigest() == source["sha256"]
    return path


def check_run(path, case, study, parameter, environment, seed, scenario, levels, height=4, dip=None):
    evidence = ROOT / "reports/results/marine_mag"
    with (evidence / f"{study}.csv").open() as stream:
        rows = [r for r in csv.DictReader(stream) if float(r["height"]) == height and (dip is None or float(r["dip"]) == dip)
                and r["case"] == case and r["seed"] == seed
                and r["scenario"] == scenario and float(r[parameter]) in levels]
    assert len(rows) == len(levels) and {float(row[parameter]) for row in rows} == levels
    for row in rows:
        env = {k: v for k, v in os.environ.items() if not k.startswith(("W3D_", "OCEAN_IMU_SIM_CASE", "MARINE_MAG_"))}
        env.update(OCEAN_IMU_SIM_CASE=case, W3D_SEED=row["seed"])
        driver = DRIVER.with_name("magnetic-heading-study") if environment else DRIVER
        variant = row.get("variant", "0")
        label = f"{scenario} {parameter}={row[parameter]}"
        if environment:
            env[environment] = row[parameter]
        output = subprocess.check_output([driver, path, variant, row["scenario"], row["dip"]], env=env, text=True)
        line = next(line for line in output.splitlines() if line.startswith("VALIDATION_METRICS ") and "segment=window" in line)
        fields = dict(f.split("=", 1) for f in line.split()[1:])
        assert int(fields["samples"]) == 180000
        for name in ("yaw_rms_deg", "roll_rms_deg", "pitch_rms_deg", "disp_3d_rms_m"):
            assert math.isclose(float(fields[name]), float(row[name]), rel_tol=1e-5, abs_tol=1e-6), (
                case, label, name, fields[name], row[name])
        if parameter == "variant" and variant == "11":
            events = [int(s.split()[1]) for s in output.splitlines() if s.startswith("INSTALLATION_RELEASE ")]
            assert events == ([] if int(row["release_s"]) == -1 else [int(row["release_s"])]), (case, label, events)
        print(f'{case} {label} reproduced yaw RMS {fields["yaw_rms_deg"]} degrees')


if __name__ == "__main__":
    data = Path(sys.argv[1])
    direction = sys.argv[2] if len(sys.argv) > 2 else "following"
    path = study_input(data, direction, 4)
    check_run(path, direction, "validation", "variant", None, "101", "burst", {0, 8})
    check_run(path, direction, "yaw-weight-validation", "sigma_factor", "MARINE_MAG_SIGMA_FACTOR",
              "421", "combined", {.5, 1})
    check_run(path, direction, "gyro-bias-validation", "gyro_bias_factor", "MARINE_MAG_GYRO_BIAS_FACTOR",
              "547", "combined", {1, 4})
    for scenario in ("ideal", "clean", "combined"):
        check_run(path, direction, "installation-screen", "variant", None, "6011", scenario, {-1, 0, 11})
    low_wave = study_input(data, direction, .27)
    for dip in (45, 76):
        check_run(low_wave, direction, "installation-screen-weak", "variant", None, "7001",
                  "combined", {0, 11}, height=.27, dip=dip)
