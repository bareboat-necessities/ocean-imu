"""Compare production filter configurations on the pinned release records."""
import argparse
import csv
import math
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from sim_dataset import REFERENCE_CSV_SHA256_BY_CASE, ROOT, digest

BINARIES = {
    "tfg_shared": "kalman_tfg_shared_modes/kalman_tfg_shared_modes-sim",
    "q_shared": "kalman_carma_shared_modes/kalman_carma_shared_modes-sim",
    "tangent_shared": "kalman_tangent_shared_modes/tangent-replay",
    "tfg_carma": "kalman_tfg_carma/tfg_carma3-sim",
    "carma": "kalman_carma/carma3-sim",
    "ou2": "kalman_ou_ii/kalman_ou_ii-sim",
    "ou3": "kalman_ou_iii/kalman_ou_iii-sim",
    "tfg": "kalman_tfg/kalman_tfg-sim",
}
METRICS = ("yaw_rms_deg", "roll_rms_deg", "pitch_rms_deg", "disp_3d_rms_m",
           "accel_bias_3d_rms_mps2", "gyro_bias_3d_rms_radps")


def replay(task, *, seed=None, window_seconds=900):
    family, case, path, mode = task
    env = {k: v for k, v in os.environ.items() if not k.startswith(("W3D_", "TFG_", "SF_", "OCEAN_IMU_SIM_CASE"))}
    env.update(OCEAN_IMU_SIM_CASE=case, W3D_WRITE_TIMESERIES="0", W3D_COLLECT_ALL_GATES="1",
               W3D_REFERENCE_MOTION="underway", W3D_VALIDATION_WINDOW_SEC=str(window_seconds), W3D_MAG_CALIBRATION=mode)
    if seed is not None:
        env["W3D_SEED"] = str(seed)
    binary = ROOT / "tests" / BINARIES[family]
    with tempfile.TemporaryDirectory() as tmp:
        proc = subprocess.run([binary, "--input", path], env=env, cwd=tmp, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    rows = [dict(f.split("=", 1) for f in s.split()[1:]) for s in proc.stdout.splitlines()
            if s.startswith("VALIDATION_METRICS ") and "segment=window" in s]
    if len(rows) != 1 or int(rows[0]["samples"]) != 200 * window_seconds or proc.returncode not in (0, 1):
        raise RuntimeError(f"{task}: replay failed ({proc.returncode})\n{proc.stdout[-2000:]}")
    gates = [dict(item.split("=", 1) for item in line.split()[1:]) for line in proc.stdout.splitlines()
             if line.startswith("QUALITY_GATE: PASS=")]
    if len(gates) != 1 or gates[0]["PASS"] not in ("0", "1") or (proc.returncode and gates[0]["PASS"] == "1"):
        raise RuntimeError(f"{task}: missing or inconsistent gate result\n{proc.stdout[-2000:]}")
    values = {k: float(rows[0][k]) for k in METRICS}
    if not all(math.isfinite(v) for v in values.values()):
        raise ValueError(f"Nonfinite result: {task}")
    return dict(family=family, case=case, record=path.name, mode=mode, **values, gate_pass=int(gates[0]["PASS"]), gate_reason=gates[0]["REASON"], exit_code=proc.returncode)


def combine(folder, output):
    rows = []
    for family in BINARIES:
        with (folder / f"{family}.csv").open() as stream:
            arm = list(csv.DictReader(stream))
        expected = {(case, name, mode) for case, names in REFERENCE_CSV_SHA256_BY_CASE.items()
                    for name in names for mode in ("off", "native", "screened")}
        if (len(arm) != len(expected) or {(r["case"], r["record"], r["mode"]) for r in arm} != expected
                or any(r["family"] != family or r["gate_pass"] not in ("0", "1")
                       or not all(math.isfinite(float(r[k])) for k in METRICS) for r in arm)):
            raise ValueError(f"Invalid comparison panel for {family}")
        rows.extend(arm)
    write(rows, output)


def write(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("family", choices=BINARIES, nargs="?")
    parser.add_argument("--combine", type=Path, help="assemble all eight per-family CSVs")
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if (args.family is None) == (args.combine is None):
        parser.error("Supply a filter family or --combine DIRECTORY")
    if args.combine is not None:
        combine(args.combine, args.out)
        return
    tasks = []
    for case in ("following", "head"):
        records = REFERENCE_CSV_SHA256_BY_CASE[case]
        if len(records) != 8:
            raise ValueError("Expected eight scored wave records")
        for record, sha in sorted(records.items()):
            path = ROOT / ".sim-data/cases" / case / record
            if digest(path) != sha:
                raise ValueError(f"Unexpected physical input: {path}")
            tasks.extend((args.family, case, path, mode) for mode in ("off", "native", "screened"))
    with ThreadPoolExecutor(args.jobs) as pool:
        rows = list(pool.map(replay, tasks))
    write(rows, args.out)


if __name__ == "__main__":
    main()
