#!/usr/bin/env python3
"""Why the vertical %H_s error rises on the small seas, and what it is not.

All three deployed families report a higher vertical RMS error *as a fraction of
H_s* on the smallest records than on the largest.  The obvious suspect is the
r_S schedule, since that is the only sea-state-dependent knob between the tuner
and the linear block.  This driver runs the three measurements that decide the
question.

    stage `cj`     Sweep C_J -- the whole r_S schedule, multiplicatively -- one
                   record at a time, with the r_S floor released so the low end
                   of the sweep is not clipped by MIN_R_S.  If the schedule were
                   mistuned for small seas, the per-record argmin would drift
                   with H_s and the deployed value would sit off-centre on the
                   small records.

    stage `amp`    Rescale one record's amplitude over a 32x range with its
                   period and spectrum held fixed, and score each copy.  The
                   deployed record set moves H_s and T_z together, so it cannot
                   separate the two; this stage varies amplitude alone.

    stage `noise`  Score the deployed records with and without the simulated IMU
                   noise, which splits the error into a sensor-driven part and a
                   sensor-independent one.

    stage `screen` Sweep every other operating-point knob one at a time on the
                   four JONSWAP records, to say which of them has any vertical
                   headroom left at all.  Most have none.

    stage `xy`     Paired multi-seed sweep of the horizontal-anisotropy knob --
                   the one thing the screen did turn up.  It buys 3D RMS, not
                   vertical; see docs/ou-iii-horizontal-anisotropy-retune.md.

The amplitude stage rescales `roll_deg`/`pitch_deg` along with the linear
columns, exactly as `ou_validation.scale_wave_motion` does, so scale factors
much above 2 produce attitude excursions the sensor model was never meant to
carry.  Scores above `--amp-trust-max` are printed but flagged.

Typical use:

    python3 tools/ou_low_sea_error_study.py cj
    python3 tools/ou_low_sea_error_study.py amp --scales 0.0625,0.125,0.25,0.5,1,2
    python3 tools/ou_low_sea_error_study.py noise
    python3 tools/ou_low_sea_error_study.py screen
    python3 tools/ou_low_sea_error_study.py xy --family OU_II
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "tools"))

SIM_DIR = REPO_ROOT / "tests" / "kalman_ou_iii"
SIM_BIN = SIM_DIR / "kalman_ou_iii-sim"

# The `xy` stage is the only one that runs against either wrapper, because the
# horizontal-anisotropy knob is the only one both expose under its own name.
FAMILIES = {
    "OU_III": {
        "dir": REPO_ROOT / "tests" / "kalman_ou_iii",
        "bin": "kalman_ou_iii-sim",
        "knob": "OU_III_R_S_XY_FACTOR",
        "grid": (1.0, 0.8, 0.72, 0.6, 0.5),
    },
    "OU_II": {
        "dir": REPO_ROOT / "tests" / "kalman_ou_ii",
        "bin": "kalman_ou_ii-sim",
        "knob": "OU_II_R_P0_XY_FACTOR",
        "grid": (1.0, 0.8, 0.65, 0.55),
    },
}

# The eight scored records, smallest sea first within each spectrum.
ALL_RECORDS = tuple(
    f"wave_data_{fam}_H{h}.csv"
    for fam in ("jonswap", "pmstokes")
    for h in ("0.270_L14.047_A30.00_P60.00",
              "1.500_L50.710_A-30.00_P120.00",
              "4.000_L112.766_A30.00_P30.00",
              "8.500_L202.839_A-30.00_P72.00")
)

# One-factor-at-a-time screen of section 1: knob -> (grid, deployed value).
SCREEN = {
    "OU_III_RS_MSE_COEFF":      ((0.027, 0.0538, 0.108), 0.0538),
    "OU_ACC_NOISE_FLOOR_SIGMA": ((0.0148, 0.04, 0.12, 0.18), 0.12),
    "OU_III_TAU_COEFF":         ((0.9, 1.0, 1.15), 1.0),
    "OU_III_ACC_BIAS_INIT_STD": ((0.004, 0.025, 0.08), 0.004),
    "OU_III_SIGMA_COEFF":       ((0.4, 0.6, 0.9, 1.4), 0.9),
    "OU_III_R_S_XY_FACTOR":     ((0.6, 0.72, 1.0), 1.0),
}

# Paired metrics reported by the xy stage, in the order the note tabulates them.
PAIRED_METRICS = ("disp_x_rms_m", "disp_y_rms_m", "disp_z_rms_m", "disp_3d_rms_m",
                  "roll_rms_deg", "pitch_rms_deg", "yaw_rms_deg",
                  "accel_bias_3d_rms_mps2", "dir_travel_correct_pct")

SEEDS = ("1", "7", "99")

# The four JONSWAP records, smallest sea first.  H_s and T_z rise together
# across them, which is exactly the confound the `amp` stage removes.
RECORDS = (
    "wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv",
    "wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv",
    "wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv",
    "wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv",
)

# Deployed C_J of the SpectralMSE law, and a geometric grid around it.
C_J_DEPLOYED = 0.0538
C_J_GRID = (0.0135, 0.027, 0.0538, 0.108, 0.215)

# Released far below the deployed MIN_R_S = 0.15 so that the bottom of the C_J
# sweep measures the schedule rather than the safety clamp.
R_S_MIN_RELEASED = 0.001

SCORE = "disp_z_pct_refrms"


def run_sim(record: Path, env_extra: dict[str, str], args: list[str],
            binary: Path | None = None) -> dict[str, str]:
    env = dict(os.environ)
    env["W3D_VALIDATION_WINDOW_SEC"] = "900"
    env["W3D_WRITE_TIMESERIES"] = "0"
    env.update(env_extra)
    # Deliberately not check=True.  The simulator's deterministic FAIL_LIMITS are
    # cut against the deployed records at the deployed operating point, so both
    # an off-grid C_J and a rescaled record trip them by construction.  The exit
    # status carries no information here; the metrics line does.
    proc = subprocess.run(
        [str(binary or SIM_BIN), *args, "--input", record.name],
        cwd=record.parent,
        env=env,
        capture_output=True,
        text=True,
    )
    for line in proc.stdout.splitlines():
        if line.startswith("VALIDATION_METRICS"):
            return dict(re.findall(r"(\w+)=(\S+)", line))
    raise RuntimeError(
        f"no VALIDATION_METRICS for {record.name} (exit {proc.returncode}); "
        f"is the record name parseable?\n{proc.stderr[-500:]}")


def parallel(jobs, workers: int):
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(lambda job: (job[0], run_sim(*job[1:])), jobs))


def short(name: str) -> str:
    return name.replace("wave_data_", "").replace(".csv", "")[:34]


def stage_cj(args) -> None:
    jobs = [
        ((rec, cj), SIM_DIR / rec,
         {"OU_III_RS_MSE_COEFF": repr(cj), "OU_III_R_S_MIN": repr(R_S_MIN_RELEASED)},
         [])
        for rec in RECORDS
        for cj in C_J_GRID
    ]
    got = dict(parallel(jobs, args.jobs))

    print(f"Vertical RMS as % of reference RMS, r_S floor released to "
          f"{R_S_MIN_RELEASED} m*s.\n")
    print(f"{'record':36s}" + "".join(f"{cj:>9g}" for cj in C_J_GRID)
          + f"{'argmin':>9s}{'r_S':>9s}")
    for rec in RECORDS:
        vals = [float(got[(rec, cj)][SCORE]) for cj in C_J_GRID]
        argmin = C_J_GRID[min(range(len(vals)), key=vals.__getitem__)]
        rs = float(got[(rec, C_J_DEPLOYED)]["tuning_applied"])
        print(f"{short(rec):36s}" + "".join(f"{v:9.3f}" for v in vals)
              + f"{argmin:9g}{rs:9.3f}")
    print(f"\nDeployed C_J is {C_J_DEPLOYED}.  A schedule mistuned for the small "
          "seas would put\ntheir argmin somewhere else.")


def stage_amp(args) -> None:
    import ou_validation as ov

    source = SIM_DIR / args.amp_source
    columns, data = ov.read_wave_csv(source)
    stem = re.match(r"wave_data_(\w+)_H[\d.]+_(L.*)\.csv", source.name)
    if stem is None:
        raise SystemExit(f"cannot parse source record name: {source.name}")
    family, tail = stem.group(1), stem.group(2)
    base_h = float(re.search(r"_H([\d.]+)_", source.name).group(1))

    work = Path(tempfile.mkdtemp(prefix="ou_low_sea_amp_"))
    try:
        shutil.copy2(SIM_BIN, work / SIM_BIN.name)
        scales = [float(s) for s in args.scales.split(",")]
        made = []
        for scale in scales:
            name = f"wave_data_{family}_H{base_h * scale:.3f}_{tail}.csv"
            ov.write_wave_csv(work / name, columns,
                              ov.scale_wave_motion(columns, data, scale))
            made.append((scale, work / name))

        got = dict(parallel([((scale,), path, {}, []) for scale, path in made],
                            args.jobs))

        print(f"Amplitude swept over {min(scales):g}x..{max(scales):g}x on "
              f"{source.name},\nwith its period and spectrum held fixed.\n")
        print(f"{'scale':>8s}{'H_s':>9s}{'z RMS':>10s}{'z % ref':>10s}"
              f"{'roll':>8s}{'slope':>9s}")
        prev = None
        for scale, _ in made:
            m = got[(scale,)]
            pct, roll = float(m[SCORE]), float(m["roll_rms_deg"])
            hs = 4.0 * float(m["disp_z_ref_rms_m"])
            slope = ""
            if prev is not None:
                d_scale, d_pct = prev
                slope = f"{math.log(pct / d_pct) / math.log(scale / d_scale):9.3f}"
            flag = "  (attitude beyond model)" if roll > args.amp_trust_roll else ""
            print(f"{scale:8g}{hs:9.3f}{float(m['disp_z_rms_m']):10.4f}"
                  f"{pct:10.3f}{roll:8.3f}{slope:>9s}{flag}")
            prev = (scale, pct)
        print("\n`slope` is d ln(relative error) / d ln(amplitude) against the row "
              "above.\nThe MSE-optimal, drift-limited prediction is -4/7 = -0.571.")
    finally:
        shutil.rmtree(work, ignore_errors=True)


def stage_noise(args) -> None:
    jobs = [((rec, tag), SIM_DIR / rec, {}, extra)
            for rec in RECORDS
            for tag, extra in (("noise", []), ("clean", ["--no-noise"]))]
    got = dict(parallel(jobs, args.jobs))

    print("Deployed records with and without the simulated IMU noise.\n")
    print(f"{'record':36s}{'with noise':>12s}{'no noise':>10s}{'noise part':>12s}")
    for rec in RECORDS:
        noisy = float(got[(rec, "noise")][SCORE])
        clean = float(got[(rec, "clean")][SCORE])
        part = math.sqrt(max(0.0, noisy * noisy - clean * clean))
        print(f"{short(rec):36s}{noisy:12.3f}{clean:10.3f}{part:12.3f}")
    print("\n`noise part` subtracts the two in quadrature: the share of the error "
          "the\nsensor model contributes, which no filter tuning can remove.")


def geo_mean(values) -> float:
    return math.exp(sum(map(math.log, values)) / len(values))


def stage_screen(args) -> None:
    """One knob at a time on the four JONSWAP records, scored vertically."""
    records = [r for r in ALL_RECORDS if "jonswap" in r]
    jobs = []
    for knob, (grid, _deployed) in SCREEN.items():
        for value in grid:
            for rec in records:
                jobs.append(((knob, value, rec), SIM_DIR / rec, {knob: repr(value)}, []))
    got = dict(parallel(jobs, args.jobs))

    print("Vertical RMS as % of reference RMS, geometric mean over the four "
          "JONSWAP records.\n")
    print(f"{'knob':30s}{'deployed':>10s}{'best':>10s}{'at':>10s}{'gain %':>9s}")
    for knob, (grid, deployed) in SCREEN.items():
        scores = {v: geo_mean([float(got[(knob, v, r)][SCORE]) for r in records])
                  for v in grid}
        base = scores[deployed]
        best = min(scores, key=scores.__getitem__)
        gain = 100.0 * (base - scores[best]) / base
        print(f"{knob:30s}{base:10.3f}{scores[best]:10.3f}{best:10g}{gain:9.2f}")
    print("\nA knob whose best is its deployed value has no headroom on this "
          "ensemble.")


def stage_xy(args) -> None:
    """Paired multi-seed sweep of the horizontal-anisotropy knob."""
    fam = FAMILIES[args.family]
    sim_dir, binary, knob = fam["dir"], fam["dir"] / fam["bin"], fam["knob"]
    if not binary.exists():
        raise SystemExit(f"{binary} is missing; run `make -C {sim_dir} build`")
    grid = ([float(g) for g in args.grid.split(",")] if args.grid else list(fam["grid"]))
    base = grid[0]

    jobs = [((seed, value, rec), sim_dir / rec,
             {knob: repr(value), "W3D_SEED": seed}, [], binary)
            for seed in SEEDS for value in grid for rec in ALL_RECORDS]
    got = dict(parallel(jobs, args.jobs))

    cells = [(s, r) for s in SEEDS for r in ALL_RECORDS]
    print(f"{args.family} {knob}: paired % change against {base}, "
          f"n = {len(cells)} record x seed cells.")
    print("Negative is better; * marks a delta with the same sign in every cell.\n")
    print(f"{'metric':26s}" + "".join(f"{v:>15g}" for v in grid[1:]))
    for metric in PAIRED_METRICS:
        row = f"{metric:26s}"
        for value in grid[1:]:
            deltas = []
            for seed, rec in cells:
                b = float(got[(seed, base, rec)][metric])
                v = float(got[(seed, value, rec)][metric])
                if b > 0.0:
                    deltas.append(100.0 * (v - b) / b)
            unanimous = all(d < 0 for d in deltas) or all(d > 0 for d in deltas)
            row += f"{sum(deltas)/len(deltas):+13.2f}{'*' if unanimous else ' '} "
        print(row)
    print("\nTake the smallest value at which BOTH horizontal axes still gain: "
          "below it\nthe knob trades x against y, which depends on the records' "
          "fixed heading.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stage", choices=("cj", "amp", "noise", "screen", "xy"))
    parser.add_argument("--jobs", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--scales", default="0.0625,0.125,0.25,0.5,1,2")
    parser.add_argument("--amp-source", default=RECORDS[1])
    parser.add_argument("--family", choices=tuple(FAMILIES), default="OU_III",
                        help="xy stage only: which wrapper to sweep")
    parser.add_argument("--grid", default=None,
                        help="xy stage only: comma-separated values, baseline first")
    parser.add_argument("--amp-trust-roll", type=float, default=0.35,
                        help="roll RMS in deg above which a rescaled record is "
                             "flagged as outside the sensor model")
    args = parser.parse_args()

    if not SIM_BIN.exists():
        raise SystemExit(f"{SIM_BIN} is missing; run `make -C tests/kalman_ou_iii build`")
    if not (SIM_DIR / RECORDS[0]).exists():
        raise SystemExit("simulation records are missing; run `make fetch-sim-data`")

    {"cj": stage_cj, "amp": stage_amp, "noise": stage_noise,
     "screen": stage_screen, "xy": stage_xy}[args.stage](args)


if __name__ == "__main__":
    main()
