"""Post-hoc isolation of V8's quarantine using the same held-out physical packets."""
import argparse
import csv
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

from marine_mag_study import OUT, ROOT, parse, write
from sim_dataset import REFERENCE_CSV_SHA256_BY_CASE, digest

DRIVER = ROOT / "tests/kalman_tfg_shared_modes/magnetic_study/build/magnetic-quarantine-study"


def run(pair):
    row, path = pair
    env = {k: v for k, v in os.environ.items() if not k.startswith(("W3D_", "OCEAN_IMU_SIM_CASE"))}
    env.update(OCEAN_IMU_SIM_CASE=row["case"], W3D_SEED=row["seed"])
    output = subprocess.check_output([DRIVER, path, "8", row["scenario"], row["dip"]], env=env, text=True)
    scores, _ = parse(output)
    return {**{k: row[k] for k in ("case", "spectrum", "height", "seed", "scenario", "dip")}, **scores}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    with (OUT / "validation.csv").open() as stream:
        rows = [row for row in csv.DictReader(stream) if row["variant"] == "8"]
    assert len(rows) == 32
    inputs = {}
    pairs = []
    for row in rows:
        case, height, spectrum = row["case"], row["height"], row["spectrum"]
        key = case, spectrum, height
        if key not in inputs:
            source, = (ROOT / ".sim-data/cases" / case).glob(f"wave_data_{spectrum}_H{float(height):.3f}_*.csv")
            assert digest(source) == REFERENCE_CSV_SHA256_BY_CASE[case][source.name]
            inputs[key] = source
        pairs.append((row, inputs[key]))
    with ThreadPoolExecutor(args.jobs) as pool:
        write(OUT / "quarantine-only.csv", list(pool.map(run, pairs)))


if __name__ == "__main__":
    main()
