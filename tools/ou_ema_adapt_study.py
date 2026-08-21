#!/usr/bin/env python3
"""Tune the dynamically scaled EMA smoothers on the drift-correction channels.

Three channels smooth their target with an exponential moving average whose
time constant is proportional to the current OU operating point:

    OU-III  r_S    tau_ema = OU_ADAPT_RS_MULT    * tau_target
    OU-II   r_p0   tau_ema = OU_ADAPT_R_P0_MULT  * tau_target
    OU-II   r_v0   tau_ema = OU_ADAPT_R_V0_MULT  * tau_target

All three multipliers ship at 5.0, which was never measured against an
alternative.  This tool sweeps them and reports, for every setting, the error
of the resulting estimate relative to the shipped default.

Instrument.  A smoothing horizon is only observable while its target moves, so
the stationary records answer almost nothing: after warmup the target sits
still and every horizon converges to the same value.  The `transition` stage
therefore synthesizes a C2 round trip: low sea -> high sea over 400--520 s,
holds the high sea, then high sea -> the same low realization over 800--920 s.
Both crossfades retain the exact displacement/velocity/acceleration derivative
cross terms used by `tools/ou_validation.py`, so the generated record remains
kinematically closed in both directions.  Legacy one-way `up` and `down`
records remain available explicitly through `--transition-dir`.

Scoring.  Per record the simulator reports RMS errors over a window.  Each
candidate is scored by the ratio of its error to the default's error on the
same record, seed and interval, so records of wildly different sea states can
be pooled.  Ratios are pooled with a geometric mean -- a scale-free average of
relative change -- and reported alongside the worst record, because a horizon
that wins on average by losing badly somewhere is not an improvement.

Usage:
  tools/ou_ema_adapt_study.py stage1     --family OU_III
  tools/ou_ema_adapt_study.py transition --family OU_III
  tools/ou_ema_adapt_study.py grid2d     --family OU_II
  tools/ou_ema_adapt_study.py confirm    --family OU_II --point 'a:OU_ADAPT_R_P0_MULT=8'
"""

import argparse
import csv
import itertools
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

FAMILIES = {
    "OU_II": {
        "subdir": "kalman_ou_ii",
        "binary": "kalman_ou_ii-sim",
        "knobs": ["OU_ADAPT_R_P0_MULT", "OU_ADAPT_R_V0_MULT"],
        "short": ["p0", "v0"],
        "shipped_mult": 3.0,
    },
    "OU_III": {
        "subdir": "kalman_ou_iii",
        "binary": "kalman_ou_iii-sim",
        "knobs": ["OU_ADAPT_RS_MULT"],
        "short": ["rs"],
        "shipped_mult": 1.5,
    },
    # TFG carries the same r_S channel and the same tau-proportional smoothing
    # horizon, so `stage1` and `transition` mean the same thing here.  It is
    # also the third family whose adaptation coefficients are fitted with this
    # tool's `confirm` stage, which takes arbitrary environment knobs and does
    # not care that the tool is named for the smoothers.
    "TFG": {
        "subdir": "kalman_tfg",
        "binary": "kalman_tfg-sim",
        "knobs": ["TFG_ADAPT_RS_MULT"],
        "short": ["rs"],
        "shipped_mult": 3.0,
    },
}

# `shipped_mult` above is what each family's multipliers are compiled with.  A
# candidate that sets no knob is the baseline and reproduces that point exactly,
# so the only use of the number here is to keep the sweep grid from carrying a
# second, redundant copy of it.  It is per family because they no longer agree:
# OU-III's r_S horizon was refitted to 1.5 against a faster sea-state
# transition (docs/ou-ema-adaptation-tuning.md section 7) while OU-II's two
# channels and TFG's stayed at 3.0.  A single shared constant would prune the
# wrong point from the grid -- dropping the real comparison and keeping a
# candidate that merely re-runs the baseline under an environment override.

# Scoring window, matching the committed 900 s convention.
WINDOW_SEC = 900.0
RECORD_SEC = 1200.0
DT_SEC = 1.0 / 200.0

# Preserve the 120 s crossfade duration of the existing instrument, but move
# the build transition forward and add a symmetric decay transition later in
# the same record.  Full-record defaults are therefore:
#
#   400--520 s  low -> high
#   520--800 s  high sea
#   800--920 s  high -> the same low realization
#   920--1200 s low sea
#
# The score still opens at t=300 s, giving a clean low-sea baseline before the
# first transition and settled intervals on both sides of the high-sea plateau.
TRANSITION_START_SEC = 400.0
TRANSITION_END_SEC = 520.0
TRANSITION_RETURN_START_SEC = 800.0
TRANSITION_RETURN_END_SEC = 920.0
TRANSITION_DURATION_SEC = TRANSITION_END_SEC - TRANSITION_START_SEC
TRANSITION_RECOVER_END_SEC = TRANSITION_END_SEC + TRANSITION_DURATION_SEC
TRANSITION_RETURN_RECOVER_END_SEC = (
    TRANSITION_RETURN_END_SEC + TRANSITION_DURATION_SEC
)

METRIC_RENAME = {
    "disp_3d_rms_m": "rms_3d",
    "disp_z_rms_m": "z_rms",
    "disp_x_rms_m": "x_rms",
    "disp_y_rms_m": "y_rms",
    "roll_rms_deg": "roll_rms",
    "pitch_rms_deg": "pitch_rms",
    "yaw_rms_deg": "yaw_rms",
    "accel_bias_3d_rms_mps2": "acc_bias_rms_3d",
    "tau_applied_s": "tau_applied",
    "tuning_applied": "channel_applied",
}

METRICS = ("rms_3d", "z_rms", "xy_rms", "roll_pitch_rms", "yaw_rms", "acc_bias_rms_3d")

RX_GATE = re.compile(r"QUALITY_GATE:\s+PASS=(\d)\s+REASON=(.*)")


def log(msg):
    print(msg, flush=True)


def records(family):
    tdir = ROOT / "tests" / FAMILIES[family]["subdir"]
    return sorted(str(p) for p in tdir.glob("wave_data_*.csv"))


def parse(text):
    """One row per VALIDATION_METRICS line, plus the record's quality gate."""
    rows = []
    gate_pass, gate_reason = False, ""

    for line in text.splitlines():
        m = RX_GATE.search(line)
        if m:
            gate_pass = m.group(1) == "1"
            gate_reason = m.group(2).strip()
            continue
        # The trailing window is tagged VALIDATION_METRICS; each extra scoring
        # interval requested through W3D_VALIDATION_SEGMENTS is tagged
        # VALIDATION_SEGMENT.  Both carry the same key=value payload.
        if not (line.startswith("VALIDATION_METRICS ")
                or line.startswith("VALIDATION_SEGMENT ")):
            continue

        row = {}
        for token in line.split()[1:]:
            key, _, value = token.partition("=")
            key = METRIC_RENAME.get(key, key)
            try:
                row[key] = float(value)
            except ValueError:
                row[key] = value
        row["wave_dataset"] = os.path.basename(str(row.get("input", "")))
        row["xy_rms"] = math.hypot(row.get("x_rms", 0.0), row.get("y_rms", 0.0))
        row["roll_pitch_rms"] = math.hypot(row.get("roll_rms", 0.0),
                                           row.get("pitch_rms", 0.0))
        rows.append(row)

    for row in rows:
        row["gate_pass"] = gate_pass
        row["gate_reason"] = gate_reason
    return rows


def run_one(job):
    """Run the simulator on a single record with a single knob setting."""
    family, label, knobs, seed, record, segments, timeout_sec = job
    cfg = FAMILIES[family]
    tdir = ROOT / "tests" / cfg["subdir"]

    env = os.environ.copy()
    # "default" reproduces the shipped deterministic realization, which every
    # committed number in the repository is scored on.  Any other value expands
    # into independent sensor streams.
    if seed == "default":
        env.pop("W3D_SEED", None)
    else:
        env["W3D_SEED"] = str(seed)
    # The study reads only the RMS summaries.  Skipping the per-sample CSV both
    # speeds the run up and lets candidates run in parallel in one directory
    # without fighting over output filenames.
    env["W3D_WRITE_TIMESERIES"] = "0"
    env["W3D_COLLECT_ALL_GATES"] = "1"
    env["W3D_VALIDATION_WINDOW_SEC"] = f"{WINDOW_SEC:g}"
    if segments:
        env["W3D_VALIDATION_SEGMENTS"] = segments
    for k, v in knobs.items():
        env[k] = f"{v:.9g}"

    # Keyed by the full path, not the basename: synthesized transition records
    # all carry the same conforming filename in separate directories, and a
    # basename key would collapse them onto one baseline.
    stub = {"family": family, "label": label, "seed": seed,
            "wave_dataset": record, **knobs}

    try:
        pr = subprocess.run(
            [str(tdir / cfg["binary"]), "--input", record],
            cwd=tdir, text=True, capture_output=True, env=env, timeout=timeout_sec,
        )
    except subprocess.TimeoutExpired:
        return [{**stub, "segment": "window", "gate_pass": False,
                 "gate_reason": "timeout", "returncode": -2}]

    rows = parse(pr.stdout + "\n" + pr.stderr)
    if not rows:
        return [{**stub, "segment": "window", "gate_pass": False,
                 "gate_reason": f"no_output(rc={pr.returncode})",
                 "returncode": pr.returncode}]

    for r in rows:
        r.update({**stub, "returncode": pr.returncode})
    return rows


def run_candidates(family, candidates, seeds, recs, jobs, timeout_sec, segments=""):
    work = [
        (family, label, knobs, seed, rec, segments, timeout_sec)
        for label, knobs in candidates
        for seed in seeds
        for rec in recs
    ]
    log(f"RUN family={family} candidates={len(candidates)} seeds={len(seeds)} "
        f"records={len(recs)} runs={len(work)}")

    out, done = [], 0
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        for rows in ex.map(run_one, work, chunksize=1):
            out.extend(rows)
            done += 1
            if done % 25 == 0 or done == len(work):
                log(f"  progress {done}/{len(work)}")
    return out


def geomean(values):
    vals = [v for v in values if math.isfinite(v) and v > 0.0]
    if not vals:
        return float("nan")
    return math.exp(sum(math.log(v) for v in vals) / len(vals))


def summarize(rows, baseline_label="baseline", segment="window"):
    """Ratio every candidate against the baseline on the same record+seed."""
    rows = [r for r in rows if r.get("segment", "window") == segment]

    base = {}
    for r in rows:
        if r.get("label") == baseline_label:
            base[(r["family"], r["wave_dataset"], r["seed"])] = r

    by_label = defaultdict(list)
    for r in rows:
        by_label[(r["family"], r["label"])].append(r)

    summary = []
    for (family, label), items in sorted(by_label.items()):
        ratios = defaultdict(list)
        failures = regressions = 0
        for r in items:
            b = base.get((family, r["wave_dataset"], r["seed"]))
            if not (r.get("gate_pass", False) and r.get("returncode", 0) == 0):
                failures += 1
                # Several records fail the acceleration-bias gate for the
                # default too.  Only a gate the default passes and the
                # candidate breaks is evidence against the candidate.
                if b is not None and b.get("gate_pass", False):
                    regressions += 1
            if b is None:
                continue
            for m in METRICS:
                num, den = r.get(m), b.get(m)
                if num is None or den is None or den <= 0.0:
                    continue
                ratios[m].append(num / den)

        row = {"family": family, "label": label, "segment": segment,
               "n": len(items), "gate_failures": failures,
               "gate_regressions": regressions}
        for m in METRICS:
            row[f"{m}_gm"] = geomean(ratios[m])
            row[f"{m}_max"] = max(ratios[m]) if ratios[m] else float("nan")
        # One number to order candidates by.  Displacement dominates because
        # these channels exist to suppress double-integration drift; the
        # worst-record terms keep an average win from hiding a local loss.
        row["score"] = (
            3.0 * row["rms_3d_gm"] + 1.0 * row["rms_3d_max"]
            + 2.0 * row["z_rms_gm"] + 0.7 * row["z_rms_max"]
            + 0.8 * row["roll_pitch_rms_gm"] + 0.3 * row["yaw_rms_gm"]
            + 0.2 * row["acc_bias_rms_3d_gm"]
            + 100.0 * regressions
        )
        summary.append(row)

    summary.sort(key=lambda r: (r["family"], r["score"]))
    return summary


def write_rows(path, rows):
    if not rows:
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    log(f"WROTE {path}")


def print_summary(summary, title=""):
    hdr = (f"{'family':7s} {'candidate':22s} {'seg':7s} {'3D gm':>7s} {'3D max':>7s} "
           f"{'Z gm':>7s} {'Z max':>7s} {'RP gm':>7s} {'yaw gm':>7s} "
           f"{'bias gm':>7s} {'regr':>5s} {'score':>7s}")
    log("")
    if title:
        log(title)
    log(hdr)
    log("-" * len(hdr))
    for r in summary:
        log(f"{r['family']:7s} {r['label']:22s} {r['segment']:7s} "
            f"{r['rms_3d_gm']:7.4f} {r['rms_3d_max']:7.4f} "
            f"{r['z_rms_gm']:7.4f} {r['z_rms_max']:7.4f} "
            f"{r['roll_pitch_rms_gm']:7.4f} {r['yaw_rms_gm']:7.4f} "
            f"{r['acc_bias_rms_3d_gm']:7.4f} {r['gate_regressions']:5d} "
            f"{r['score']:7.4f}")
    log("")


#  Transition records


# Endpoint seas of the synthesized transitions.  "large" is the pair the
# committed validation uses.  "small" runs the same relative sea-state change
# one octave down in wave period, which is what separates the two candidate
# scaling laws: if the best horizon is the same *multiple of tau* on both
# pairs, the deployed tau-proportional law is right; if it is the same number
# of *seconds*, the scaling is not earning its place.
TRANSITION_PAIRS = {
    #          low sea (H_s, source tokens)      high sea, rescaled to H_s
    "large": (("1.500", "50.710", 1.500), ("8.500", "202.839", 4.000)),
    "small": (("0.270", "14.047", 0.270), ("1.500", "50.710", 0.700)),
}


def roundtrip_profile(ov, times):
    """High-sea weight and derivatives for the low->high->same-low record."""
    up, up_rate, up_accel = ov.smoothstep_profile(
        times, TRANSITION_START_SEC, TRANSITION_END_SEC
    )
    down, down_rate, down_accel = ov.smoothstep_profile(
        times, TRANSITION_RETURN_START_SEC, TRANSITION_RETURN_END_SEC
    )
    return up - down, up_rate - down_rate, up_accel - down_accel


def make_roundtrip_transition(ov, columns, low_data, high_data, seed):
    """Build a C2 low->high->same-low record with exact kinematic cross terms."""
    count = min(low_data.shape[0], high_data.shape[0])
    low = ov.phase_randomize_wave(columns, low_data[:count], seed * 2 + 1)
    high = ov.phase_randomize_wave(columns, high_data[:count], seed * 2 + 2)
    times = low[:, columns.index("time")]
    weight, weight_rate, weight_accel = roundtrip_profile(ov, times)

    displacement = [columns.index(name) for name in ("disp_x", "disp_y", "disp_z")]
    velocity = [columns.index(name) for name in ("vel_x", "vel_y", "vel_z")]
    acceleration = [columns.index(name) for name in ("acc_x", "acc_y", "acc_z")]
    attitude = [columns.index(name) for name in ("roll_deg", "pitch_deg", "yaw_deg")]

    result = low.copy()
    w = weight[:, None]
    dw = weight_rate[:, None]
    ddw = weight_accel[:, None]
    displacement_delta = high[:, displacement] - low[:, displacement]
    velocity_delta = high[:, velocity] - low[:, velocity]

    result[:, displacement] = (1.0 - w) * low[:, displacement] + w * high[:, displacement]
    result[:, velocity] = (
        (1.0 - w) * low[:, velocity]
        + w * high[:, velocity]
        + dw * displacement_delta
    )
    result[:, acceleration] = (
        (1.0 - w) * low[:, acceleration]
        + w * high[:, acceleration]
        + 2.0 * dw * velocity_delta
        + ddw * displacement_delta
    )
    result[:, attitude] = (1.0 - w) * low[:, attitude] + w * high[:, attitude]
    return ov.rebuild_body_imu(columns, result)


def transition_segment_spec(directions):
    """Return simulator segment specification and display order."""
    if "roundtrip" in directions:
        if len(directions) != 1:
            raise SystemExit("roundtrip cannot be combined with one-way transition directions")
        window_start = RECORD_SEC - WINDOW_SEC
        segments = (
            f"low_start:{window_start:g}:{TRANSITION_START_SEC:g},"
            f"rise:{TRANSITION_START_SEC:g}:{TRANSITION_END_SEC:g},"
            f"high_recover:{TRANSITION_END_SEC:g}:{TRANSITION_RECOVER_END_SEC:g},"
            f"high:{TRANSITION_RECOVER_END_SEC:g}:{TRANSITION_RETURN_START_SEC:g},"
            f"fall:{TRANSITION_RETURN_START_SEC:g}:{TRANSITION_RETURN_END_SEC:g},"
            f"low_recover:{TRANSITION_RETURN_END_SEC:g}:"
            f"{TRANSITION_RETURN_RECOVER_END_SEC:g},"
            f"low_return:{TRANSITION_RETURN_RECOVER_END_SEC:g}:{RECORD_SEC:g}"
        )
        report = [
            "window", "low_start", "rise", "high_recover", "high",
            "fall", "low_recover", "low_return",
        ]
        return segments, report

    window_start = RECORD_SEC - WINDOW_SEC
    segments = (
        f"start:{window_start:g}:{TRANSITION_START_SEC:g},"
        f"blend:{TRANSITION_START_SEC:g}:{TRANSITION_END_SEC:g},"
        f"recover:{TRANSITION_END_SEC:g}:{TRANSITION_RECOVER_END_SEC:g},"
        f"end:{TRANSITION_RECOVER_END_SEC:g}:{RECORD_SEC:g}"
    )
    return segments, ["window", "blend", "recover", "end", "start"]


def build_transition_records(family, seeds, out_dir, directions, pairs):
    """Synthesize round-trip (default) or legacy one-way transition records."""
    sys.path.insert(0, str(ROOT / "tools"))
    import ou_validation as ov

    data_dir = ROOT / "tests" / FAMILIES[family]["subdir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []

    for pair in pairs:
        (lo_h, lo_l, _), (hi_h, hi_l, hi_target) = TRANSITION_PAIRS[pair]
        lo_path = ov.find_default_input(data_dir, lo_h, lo_l)
        hi_path = ov.find_default_input(data_dir, hi_h, hi_l)
        columns, lo_data = ov.read_wave_csv(lo_path, RECORD_SEC)
        _, hi_data = ov.read_wave_csv(hi_path, RECORD_SEC)
        # The high sea keeps the source record's period and is rescaled in
        # height, exactly as the committed validation builds its transition.
        hi_scaled = ov.scale_wave_motion(columns, hi_data, hi_target / float(hi_h))

        for direction in directions:
            if direction == "roundtrip":
                start = end = None
            elif direction == "up":
                start, end = lo_data, hi_scaled
            elif direction == "down":
                start, end = hi_scaled, lo_data
            else:
                raise SystemExit(f"unknown transition direction {direction}")
            for seed in seeds:
                if direction == "roundtrip":
                    generated = make_roundtrip_transition(
                        ov, columns, lo_data, hi_scaled, seed
                    )
                else:
                    generated = ov.make_nonstationary_wave(
                        columns, start, end, seed, 1.0,
                        TRANSITION_START_SEC, TRANSITION_END_SEC,
                    )
                # The simulator infers the spectrum and the nominal H_s and
                # azimuth from the filename, so a synthesized record has to
                # carry a conforming name.  A per-record directory keeps the
                # names unique without inventing sea-state tokens that do not
                # describe the file.
                sub = out_dir / f"{pair}_{direction}_{seed}"
                sub.mkdir(parents=True, exist_ok=True)
                name = (f"wave_data_jonswap_H{hi_target:.3f}_L{hi_l}"
                        f"_A-30.00_P120.00.csv")
                path = sub / name
                ov.write_wave_csv(path, columns, generated)
                paths.append(str(path))
                log(f"  generated {pair}/{direction} seed={seed} -> {path}")
    return paths


#  Candidate construction


def parse_grid(text):
    return [float(v) for v in text.split(",") if v.strip()]


def stage1_candidates(family, grid):
    """One knob at a time, every other knob at its shipped value."""
    cfg = FAMILIES[family]
    shipped = cfg["shipped_mult"]
    cands = [("baseline", {})]
    for knob, short in zip(cfg["knobs"], cfg["short"]):
        for v in grid:
            if abs(v - shipped) < 1e-9:
                continue
            cands.append((f"{short}={v:g}", {knob: v}))
    return cands


def grid2d_candidates(family, grid_a, grid_b):
    cfg = FAMILIES[family]
    if len(cfg["knobs"]) != 2:
        raise SystemExit(f"grid2d needs a two-knob family; {family} has "
                         f"{len(cfg['knobs'])}")
    ka, kb = cfg["knobs"]
    sa, sb = cfg["short"]
    shipped = cfg["shipped_mult"]
    cands = [("baseline", {})]
    for a, b in itertools.product(grid_a, grid_b):
        if abs(a - shipped) < 1e-9 and abs(b - shipped) < 1e-9:
            continue
        cands.append((f"{sa}={a:g},{sb}={b:g}", {ka: a, kb: b}))
    return cands


def confirm_candidates(points):
    """points: list of "name:KNOB=value[:KNOB=value]" strings."""
    cands = [("baseline", {})]
    for spec in points:
        name, _, rest = spec.partition(":")
        knobs = {}
        for item in rest.split(":"):
            if not item.strip():
                continue
            k, _, v = item.partition("=")
            knobs[k.strip()] = float(v)
        cands.append((name, knobs))
    return cands


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["stage1", "grid2d", "confirm", "transition"])
    ap.add_argument("--family", default="OU_III", choices=sorted(FAMILIES))
    ap.add_argument("--grid", default="0.5,1,1.5,2,3,5,8,12,20,32,50",
                    help="stage1/transition multiplier grid")
    ap.add_argument("--grid-a", default="2,3,5,8,12", help="grid2d: first knob")
    ap.add_argument("--grid-b", default="2,3,5,8,12", help="grid2d: second knob")
    ap.add_argument("--point", action="append", default=[],
                    help="confirm stage: name:KNOB=value[:KNOB=value]")
    ap.add_argument("--seeds", default="default",
                    help="comma-separated; 'default' keeps the shipped realization")
    ap.add_argument("--records", default="",
                    help="comma-separated substrings selecting stationary records")
    ap.add_argument("--transition-seeds", default="11,12,13",
                    help="wave phase seeds for the synthesized transition records")
    ap.add_argument("--transition-dir", default="roundtrip",
                    help="'roundtrip' (default), or legacy 'up', 'down', or 'up,down'")
    ap.add_argument("--transition-pair", default="large",
                    help=f"endpoint seas: {','.join(TRANSITION_PAIRS)}")
    ap.add_argument("--transition-stage", default="stage1",
                    choices=["stage1", "grid2d", "confirm"],
                    help="which candidate set to run on the transition records")
    ap.add_argument("--keep-generated", default="",
                    help="directory to keep the synthesized records in")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2)))
    ap.add_argument("--timeout", type=float, default=1800.0)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    seeds = [s.strip() if s.strip() == "default" else int(s)
             for s in args.seeds.split(",") if s.strip()]

    def candidates_for(stage):
        if stage == "stage1":
            return stage1_candidates(args.family, parse_grid(args.grid))
        if stage == "grid2d":
            return grid2d_candidates(args.family, parse_grid(args.grid_a),
                                     parse_grid(args.grid_b))
        if not args.point:
            raise SystemExit("confirm stage needs at least one --point")
        return confirm_candidates(args.point)

    tmp_dir = None
    segments = ""
    if args.stage == "transition":
        cands = candidates_for(args.transition_stage)
        gen_root = (Path(args.keep_generated) if args.keep_generated
                    else Path(tempfile.mkdtemp(prefix="ou_ema_transition_")))
        if not args.keep_generated:
            tmp_dir = gen_root
        directions = [d.strip() for d in args.transition_dir.split(",") if d.strip()]
        segments, report_segments = transition_segment_spec(directions)
        log(f"GENERATE transition records in {gen_root}")
        recs = build_transition_records(
            args.family,
            [int(s) for s in args.transition_seeds.split(",") if s.strip()],
            gen_root,
            directions,
            [p.strip() for p in args.transition_pair.split(",") if p.strip()],
        )
    else:
        cands = candidates_for(args.stage)
        recs = records(args.family)
        if args.records:
            want = [s.strip() for s in args.records.split(",") if s.strip()]
            recs = [r for r in recs if any(w in os.path.basename(r) for w in want)]
        report_segments = ["window"]

    if not recs:
        raise SystemExit(f"no records for {args.family}; run `make fetch-sim-data`")

    try:
        rows = run_candidates(args.family, cands, seeds, recs, args.jobs,
                              args.timeout, segments)
    finally:
        if tmp_dir:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    all_summary = []
    for seg in report_segments:
        summary = summarize(rows, "baseline", seg)
        if not summary:
            continue
        print_summary(summary, title=f"segment: {seg}")
        all_summary.extend(summary)

    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        write_rows(out, rows)
        write_rows(out.with_name(out.stem + "_summary.csv"), all_summary)

    return 0


if __name__ == "__main__":
    sys.exit(main())
