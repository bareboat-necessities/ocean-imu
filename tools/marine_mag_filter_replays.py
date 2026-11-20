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


def replay(task):
    family, case, path, mode = task
    env = {k: v for k, v in os.environ.items() if not k.startswith(("W3D_", "TFG_", "SF_", "OCEAN_IMU_SIM_CASE"))}
    env.update(OCEAN_IMU_SIM_CASE=case, W3D_WRITE_TIMESERIES="0", W3D_COLLECT_ALL_GATES="1",
               W3D_REFERENCE_MOTION="underway", W3D_MAG_CALIBRATION=mode)
    binary = ROOT / "tests" / BINARIES[family]
    with tempfile.TemporaryDirectory() as tmp:
        proc = subprocess.run([binary, "--input", path], env=env, cwd=tmp, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    rows = [dict(f.split("=", 1) for f in s.split()[1:]) for s in proc.stdout.splitlines()
            if s.startswith("VALIDATION_METRICS ") and "segment=window" in s]
    if len(rows) != 1 or int(rows[0]["samples"]) != 180000 or proc.returncode not in (0, 1):
        raise RuntimeError(f"{task}: replay failed ({proc.returncode})\n{proc.stdout[-2000:]}")
    if proc.returncode and not any(s in proc.stdout.lower() for s in ("fail", "exceed", "limit")):
        raise RuntimeError(proc.stdout[-2000:])
    values = {k: float(rows[0][k]) for k in METRICS}
    if not all(math.isfinite(v) for v in values.values()):
        raise ValueError(f"Nonfinite result: {task}")
    return dict(family=family, case=case, record=path.name, mode=mode, **values, exit_code=proc.returncode)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("family", choices=BINARIES)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
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
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
