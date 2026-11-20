"""Score external-map and attitude interventions at the actual TFG update boundary."""
import csv
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "tests/kalman_tfg_shared_modes/magnetic_study/build/magnetic-reference-study"


def main():
    rows = []
    for case in ("following", "head"):
        source, = (ROOT / ".sim-data/cases" / case).glob("wave_data_jonswap_H1.500_*.csv")
        env = {k: v for k, v in os.environ.items() if not k.startswith(("W3D_", "OCEAN_IMU_SIM_CASE"))}
        env["OCEAN_IMU_SIM_CASE"] = case
        for intervention in ("native", "map", "wrong_dip", "handoff_tilt", "handoff_yaw"):
            output = subprocess.check_output([DRIVER, source, intervention], env=env, text=True)
            if intervention.startswith("handoff"):
                assert output.count("HANDOFF_INJECTED") == 1
            fields = {}
            for line in output.splitlines():
                if line.startswith(("VALIDATION_METRICS ", "VALIDATION_METRICS_MEANS ")) and "segment=window" in line:
                    fields.update(dict(x.split("=", 1) for x in line.split()[1:]))
            assert int(fields["samples"]) == 180000
            rows.append(dict(case=case, intervention=intervention,
                             **{k: fields[k] for k in ("yaw_rms_deg", "yaw_mean_deg",
                                "roll_rms_deg", "pitch_rms_deg", "disp_3d_rms_m")}))
    with (ROOT / "reports/results/marine_mag/reference.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
