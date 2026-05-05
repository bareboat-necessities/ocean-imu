#!/usr/bin/env python3
import argparse
import csv
import math
import os
import random
import re
import subprocess
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MAX_BOOT_TIMEOUT_SEC = 15.0
MIN_MARGIN_SEC = 1.0
MAX_MARGIN_SEC = 4.0

REJECTED_KEYS = {"OU_TAU_COEFF", "OU_SIGMA_COEFF"}
REJECTED_PREFIXES = ("SF_MAG_",)

PARAM_RE = re.compile(r"^(SF_|OU_)")

GRAVITY_RANGES = {
    "OU_II": {
        "SF_BOOT_TILT_ACC_TAU": (1.0, 4.0),
        "SF_BOOT_GRAV_SLOW_TAU": (2.5, 8.0),
        "SF_BOOT_GRAV_ALIGN_MAX_SIN": (0.08, 0.25),
        "SF_BOOT_GRAV_HOLD_SEC": (0.4, 2.2),
        "SF_BOOT_GRAV_MIN_SEC": (3.0, 10.0),
        "SF_BOOT_GRAV_NORM_FRAC": (0.18, 0.55),
        "SF_ONLINE_TUNE_WARMUP_SEC": (3.0, 15.0),
    },
    "OU_III": {
        "SF_BOOT_TILT_ACC_TAU": (1.0, 4.0),
        "SF_BOOT_GRAV_SLOW_TAU": (2.5, 8.0),
        "SF_BOOT_GRAV_ALIGN_MAX_SIN": (0.08, 0.25),
        "SF_BOOT_GRAV_HOLD_SEC": (0.4, 2.2),
        "SF_BOOT_GRAV_MIN_SEC": (3.0, 10.0),
        "SF_BOOT_GRAV_NORM_FRAC": (0.18, 0.55),
        "SF_ONLINE_TUNE_WARMUP_SEC": (3.0, 15.0),
    },
}

OU_COMMON = {
    "OU_ACC_NOISE_FLOOR_SIGMA": (0.06, 0.25),
    "OU_ADAPT_TAU_SEC": (0.8, 6.0),
    "OU_ADAPT_EVERY_SECS": (0.04, 0.35),
}

OU_II_ONLY = {
    "OU_P_FACTOR": (0.8, 2.5),
    "OU_R_P0_XY_FACTOR": (0.10, 0.70),
    "OU_R_P0_COEFF": (0.6, 3.5),
    "OU_R_V0_COEFF": (0.6, 3.0),
}

RX = {
    "ds": re.compile(r"Processing\s+(.+?\.csv)"),
    "ang": re.compile(
        r"Angles RMS \(deg\): Roll=([0-9.eE+-]+) Pitch=([0-9.eE+-]+) Yaw=([0-9.eE+-]+)"
    ),
    "xyz": re.compile(
        r"XYZ RMS \(m\): X=([0-9.eE+-]+) Y=([0-9.eE+-]+) Z=([0-9.eE+-]+)"
    ),
    "r3d": re.compile(r"3D RMS \(m\):\s*([0-9.eE+-]+)"),
    "accb": re.compile(
        r"Bias error RMS \(acc, m/s\^2\): X=([0-9.eE+-]+) Y=([0-9.eE+-]+) Z=([0-9.eE+-]+) \|3D\|=([0-9.eE+-]+)"
    ),
    "gate": re.compile(r"QUALITY_GATE:\s+PASS=(\d)\s+REASON=(.*)"),
}


def clean_params(p):
    out = {}
    for k, v in p.items():
        if not PARAM_RE.match(k):
            continue
        if k in REJECTED_KEYS:
            continue
        if k.startswith(REJECTED_PREFIXES):
            continue
        if isinstance(v, float) and not math.isfinite(v):
            continue
        out[k] = v
    return out


def format_duration(sec):
    sec = max(0, int(sec))
    h, rem = divmod(sec, 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def lhs(rng, n, span):
    lo, hi = span
    vals = [lo + (i + rng.random()) / n * (hi - lo) for i in range(n)]
    rng.shuffle(vals)
    return vals


def minimal_candidates():
    cands = [
        {"SF_BOOT_GRAV_HOLD_SEC": 1.0},
        {"SF_BOOT_GRAV_MIN_SEC": 6.0},
        {"SF_BOOT_GRAV_NORM_FRAC": 0.35},
        {"SF_BOOT_GRAV_SLOW_TAU": 5.0},
        {"SF_BOOT_TILT_ACC_TAU": 2.0},
        {"SF_ONLINE_TUNE_WARMUP_SEC": 8.0},
        {"SF_BOOT_GRAV_HOLD_SEC": 1.0, "SF_BOOT_GRAV_MIN_SEC": 6.0},
        {"SF_BOOT_GRAV_NORM_FRAC": 0.35, "SF_BOOT_GRAV_SLOW_TAU": 5.0},
        {"SF_BOOT_TILT_ACC_TAU": 2.0, "SF_BOOT_GRAV_HOLD_SEC": 1.0},
        {"SF_ONLINE_TUNE_WARMUP_SEC": 8.0, "SF_BOOT_GRAV_MIN_SEC": 6.0},
    ]
    return [(f"minimal_{i:03d}", p) for i, p in enumerate(cands, 1)]


def explicit_practical_candidates():
    base = {
        "SF_BOOT_TILT_ACC_TAU": 2.5,
        "SF_BOOT_GRAV_SLOW_TAU": 5.5,
        "SF_BOOT_GRAV_ALIGN_MAX_SIN": 0.12,
        "SF_BOOT_GRAV_HOLD_SEC": 1.2,
        "SF_BOOT_GRAV_MIN_SEC": 6.0,
        "SF_BOOT_GRAV_NORM_FRAC": 0.35,
        "SF_ONLINE_TUNE_WARMUP_SEC": 8.0,
    }
    out = []
    for i, timeout in enumerate((8.0, 10.0, 12.0, 15.0), 1):
        p = dict(base)
        p["SF_BOOT_GRAV_TIMEOUT_SEC"] = timeout
        out.append((f"practical_{i:03d}", p))
    return out


def apply_timeout_constraint(p, rng):
    min_sec = p["SF_BOOT_GRAV_MIN_SEC"]
    hold_sec = p["SF_BOOT_GRAV_HOLD_SEC"]
    margin_max = min(MAX_MARGIN_SEC, MAX_BOOT_TIMEOUT_SEC - min_sec - hold_sec)
    if margin_max < MIN_MARGIN_SEC:
        return False
    margin = MIN_MARGIN_SEC + rng.random() * (margin_max - MIN_MARGIN_SEC)
    p["SF_BOOT_GRAV_TIMEOUT_SEC"] = min_sec + hold_sec + margin
    return p["SF_BOOT_GRAV_TIMEOUT_SEC"] <= MAX_BOOT_TIMEOUT_SEC


def sample_params(space, n, rng, gravity=False):
    cols = {k: lhs(rng, n, span) for k, span in space.items()}
    out = []
    i = 0
    while i < n:
        p = {k: cols[k][i] for k in cols}
        if gravity and not apply_timeout_constraint(p, rng):
            continue
        out.append(p)
        i += 1
    return out


def validate_candidate(p):
    bad = []
    timeout = p.get("SF_BOOT_GRAV_TIMEOUT_SEC")
    if timeout is not None and timeout > MAX_BOOT_TIMEOUT_SEC:
        bad.append(f"SF_BOOT_GRAV_TIMEOUT_SEC>{MAX_BOOT_TIMEOUT_SEC}")
    if timeout is not None:
        min_sec = p.get("SF_BOOT_GRAV_MIN_SEC", 0.0)
        hold_sec = p.get("SF_BOOT_GRAV_HOLD_SEC", 0.0)
        if timeout <= min_sec + hold_sec + MIN_MARGIN_SEC:
            bad.append("SF_BOOT_GRAV_TIMEOUT_SEC<=min+hold+margin")
    for key in p:
        if key in REJECTED_KEYS:
            bad.append(f"forbidden_key:{key}")
        if key.startswith(REJECTED_PREFIXES):
            bad.append(f"forbidden_key:{key}")
    return (len(bad) == 0), ";".join(bad)


def parse(text):
    ds = RX["ds"].findall(text)
    ang = [m.groups() for m in RX["ang"].finditer(text)]
    xyz = [m.groups() for m in RX["xyz"].finditer(text)]
    r3d = [m.group(1) for m in RX["r3d"].finditer(text)]
    accb = [m.groups() for m in RX["accb"].finditer(text)]
    gates = [m.groups() for m in RX["gate"].finditer(text)]

    n = max(len(ds), len(ang), len(xyz), len(accb), 1)
    out = []

    def g(arr, i, j):
        return float(arr[i][j]) if i < len(arr) else float("nan")

    for i in range(n):
        gp = (gates[i][0] == "1") if i < len(gates) else False
        reason = gates[i][1].strip() if i < len(gates) else "missing_gate_status"
        out.append(
            {
                "wave_dataset": ds[i] if i < len(ds) else f"dataset_{i}",
                "roll_rms": g(ang, i, 0),
                "pitch_rms": g(ang, i, 1),
                "yaw_rms": g(ang, i, 2),
                "x_rms": g(xyz, i, 0),
                "y_rms": g(xyz, i, 1),
                "z_rms": g(xyz, i, 2),
                "rms_3d": float(r3d[i]) if i < len(r3d) else float("nan"),
                "acc_bias_rms_3d": g(accb, i, 3),
                "quality_gate_pass": gp,
                "fail_reason": reason,
            }
        )
    return out


def run_candidate(fam, cid, p, tier, seed, collect, stage):
    tdir = ROOT / "tests" / ("kalman_ou_ii" if fam == "OU_II" else "kalman_ou_iii")
    bin_name = "./kalman_ou_ii-sim" if fam == "OU_II" else "./kalman_ou_iii-sim"

    env = os.environ.copy()
    env["W3D_SEED"] = str(seed)
    env["W3D_TIER"] = tier
    if collect:
        env["W3D_COLLECT_ALL_GATES"] = "1"

    input_params = clean_params(p)
    env.update(
        {
            k: (f"{v:.9g}" if isinstance(v, float) else str(v))
            for k, v in input_params.items()
        }
    )

    # IMPORTANT: no --nomag
    pr = subprocess.run([bin_name], cwd=tdir, text=True, capture_output=True, env=env)

    rows = parse(pr.stdout + "\n" + pr.stderr)
    for r in rows:
        r.update(
            {
                "family": fam,
                "candidate": cid,
                "seed": seed,
                "tier": tier,
                "stage": stage,
                "returncode": pr.returncode,
                "roll_pitch_rms_norm": math.hypot(r["roll_rms"], r["pitch_rms"]),
                "xy_rms": math.hypot(r["x_rms"], r["y_rms"]),
                "_params": dict(input_params),
                "input_params": dict(input_params),
                "rejected_before_run": 0,
                "reject_reason": "",
            }
        )
        timeout = input_params.get("SF_BOOT_GRAV_TIMEOUT_SEC")
        r["timeout_valid"] = int(
            timeout is not None
            and isinstance(timeout, (float, int))
            and math.isfinite(float(timeout))
            and float(timeout) <= MAX_BOOT_TIMEOUT_SEC
        )
        for k, v in input_params.items():
            r[k] = v
    return rows


def percentile(vals, p):
    if not vals:
        return float("inf")
    s = sorted(vals)
    i = (len(s) - 1) * p
    lo = int(i)
    hi = min(len(s) - 1, lo + 1)
    f = i - lo
    return s[lo] * (1 - f) + s[hi] * f


def score_rows(rows):
    by = defaultdict(list)
    for r in rows:
        if not r.get("rejected_before_run"):
            by[(r["family"], r["candidate"], r["seed"], r["tier"])].append(r)

    base = {
        (r["family"], r["wave_dataset"], r["seed"], r["tier"]): r
        for r in rows
        if r.get("candidate") == "baseline" and not r.get("rejected_before_run")
    }

    def rf(x):
        return x if math.isfinite(x) and abs(x) > 1e-9 else 1e-9

    for (_, _, _, _), it in by.items():
        rp, z, r3, yaw, acc = [], [], [], [], []
        any_fail = False

        for r in it:
            b = base.get((r["family"], r["wave_dataset"], r["seed"], r["tier"]))

            if r["candidate"] != "baseline" and b is None:
                r["score"] = float("inf")
                fr = (r.get("fail_reason") or "").strip(";")
                r["fail_reason"] = (fr + ";missing_matching_baseline").strip(";")
                r["quality_gate_pass"] = False
                any_fail = True
                continue

            if b is None:
                # baseline self-compare
                b = r

            r["roll_pitch_ratio"] = rf(r["roll_pitch_rms_norm"]) / rf(
                b.get("roll_pitch_rms_norm", float("nan"))
            )
            r["z_ratio"] = rf(r["z_rms"]) / rf(b.get("z_rms", float("nan")))
            r["rms3d_ratio"] = rf(r["rms_3d"]) / rf(b.get("rms_3d", float("nan")))
            r["yaw_ratio"] = rf(r["yaw_rms"]) / rf(b.get("yaw_rms", float("nan")))
            r["acc_bias_ratio"] = rf(r["acc_bias_rms_3d"]) / rf(
                b.get("acc_bias_rms_3d", float("nan"))
            )

            rp.append(r["roll_pitch_ratio"])
            z.append(r["z_ratio"])
            r3.append(r["rms3d_ratio"])
            yaw.append(r["yaw_ratio"])
            acc.append(r["acc_bias_ratio"])
            any_fail = any_fail or (not r.get("quality_gate_pass", False))

        if not rp or not z or not r3 or not yaw or not acc:
            s = float("inf")
            for r in it:
                fr = (r.get("fail_reason") or "").strip(";")
                r["fail_reason"] = (fr + ";no_scoreable_rows").strip(";")
                r["score"] = s
            continue

        s = (
            4.0 * percentile(rp, 0.75)
            + 2.0 * max(rp)
            + 2.0 * max(0.0, max(z) - 1.02)
            + 2.0 * max(0.0, max(r3) - 1.02)
            + 0.5 * percentile(yaw, 0.75)
            + 0.5 * percentile(acc, 0.75)
            + (100.0 if any_fail else 0.0)
        )
        for r in it:
            r["score"] = s


def aggregate_candidates(rows, family, tier):
    out = []
    by = defaultdict(list)
    for r in rows:
        if r.get("family") == family and r.get("tier") == tier:
            by[r["candidate"]].append(r)

    for cand, it in by.items():
        valid = all(
            rr.get("returncode") == 0
            and rr.get("quality_gate_pass")
            and not rr.get("rejected_before_run")
            and math.isfinite(float(rr.get("score", float("inf"))))
            for rr in it
        )
        params = clean_params(it[0].get("_params", it[0].get("input_params", {})))
        obj = {
            "family": family,
            "candidate": cand,
            "tier": tier,
            "params": params,
            "valid": valid,
            "n_rows": len(it),
            "stage": it[0].get("stage", "unknown"),
        }
        if valid:
            obj.update(
                {
                    "mean_score": sum(r["score"] for r in it) / len(it),
                    "max_score": max(r["score"] for r in it),
                    "mean_rms3d": sum(r["rms_3d"] for r in it) / len(it),
                    "mean_z_rms": sum(r["z_rms"] for r in it) / len(it),
                    "mean_roll_pitch_rms": sum(r["roll_pitch_rms_norm"] for r in it)
                    / len(it),
                    "mean_yaw_rms": sum(r["yaw_rms"] for r in it) / len(it),
                    "max_rms3d_ratio": max(r.get("rms3d_ratio", float("inf")) for r in it),
                    "max_z_ratio": max(r.get("z_ratio", float("inf")) for r in it),
                    "max_roll_pitch_ratio": max(
                        r.get("roll_pitch_ratio", float("inf")) for r in it
                    ),
                    "max_yaw_ratio": max(r.get("yaw_ratio", float("inf")) for r in it),
                }
            )
        else:
            obj.update(
                {
                    "mean_score": float("inf"),
                    "max_score": float("inf"),
                    "mean_rms3d": float("inf"),
                    "mean_z_rms": float("inf"),
                    "mean_roll_pitch_rms": float("inf"),
                    "mean_yaw_rms": float("inf"),
                    "max_rms3d_ratio": float("inf"),
                    "max_z_ratio": float("inf"),
                    "max_roll_pitch_ratio": float("inf"),
                    "max_yaw_ratio": float("inf"),
                }
            )
        out.append(obj)
    return out


def rank_key(a):
    return (
        a["mean_score"],
        a["max_score"],
        a["max_rms3d_ratio"],
        a["max_z_ratio"],
        a["mean_rms3d"],
        a["mean_z_rms"],
        a["mean_roll_pitch_rms"],
        a["mean_yaw_rms"],
    )


def print_stage_summary(stage_name, agg):
    valid = [x for x in agg if x["valid"]]
    print(f"STAGE_{stage_name}_VALID_COUNT {len(valid)}")
    top = sorted(valid, key=rank_key)[:5]
    print(f"STAGE_{stage_name}_TOP {','.join(t['candidate'] for t in top) if top else 'none'}")
    for t in top:
        print(
            f"STAGE_{stage_name}_TOP_CAND candidate={t['candidate']} score={t['mean_score']:.6g} "
            f"rms3d={t['mean_rms3d']:.6g} z={t['mean_z_rms']:.6g} "
            f"rp={t['mean_roll_pitch_rms']:.6g} yaw={t['mean_yaw_rms']:.6g} params={t['params']}"
        )


def print_param_sensitivity(rows, fam):
    valid = [
        r
        for r in rows
        if r.get("family") == fam
        and r.get("returncode") == 0
        and r.get("quality_gate_pass")
        and not r.get("rejected_before_run")
        and r.get("candidate") != "baseline"
        and math.isfinite(float(r.get("score", float("nan"))))
    ]

    param_data = defaultdict(lambda: {"x": [], "y": []})
    for r in valid:
        for k, v in r.get("input_params", {}).items():
            if (
                PARAM_RE.match(k)
                and k not in REJECTED_KEYS
                and not k.startswith(REJECTED_PREFIXES)
                and isinstance(v, (float, int))
                and math.isfinite(float(v))
            ):
                param_data[k]["x"].append(float(v))
                param_data[k]["y"].append(float(r["score"]))

    print(f"PARAM_SENSITIVITY {fam}")
    print("param,n,min,max,corr_with_score")
    if not param_data:
        print("NONE,0,nan,nan,nan")
        return

    def corr(xs, ys):
        n = len(xs)
        if n < 2:
            return float("nan")
        mx = sum(xs) / n
        my = sum(ys) / n
        vx = sum((x - mx) ** 2 for x in xs)
        vy = sum((y - my) ** 2 for y in ys)
        if vx < 1e-14 or vy < 1e-14:
            return float("nan")
        cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        return cov / math.sqrt(vx * vy)

    for k in sorted(param_data):
        xs = param_data[k]["x"]
        ys = param_data[k]["y"]
        mn, mx = min(xs), max(xs)
        c = corr(xs, ys)
        cstr = f"{c:.6g}" if math.isfinite(c) else "nan"
        print(f"{k},{len(xs)},{mn:.6g},{mx:.6g},{cstr}")
        if abs(mx - mn) < 1e-12:
            print(f"WARNING_PARAM_NO_VARIANCE {k}")


def print_cpp_config(fam, result):
    params = clean_params(result.get("params", {}))
    print(
        f"\n=== BEST_CONFIG {fam} candidate={result['candidate']} tier={result['tier']} rows={result['n_rows']} ==="
    )
    print(
        f"RMS_SUMMARY mean_3d={result['mean_rms3d']:.6g} mean_z={result['mean_z_rms']:.6g} "
        f"mean_roll_pitch={result['mean_roll_pitch_rms']:.6g} mean_yaw={result['mean_yaw_rms']:.6g}"
    )
    print("// C++ snippet:")
    print("struct TunedParams {")
    for k in sorted(params):
        v = params[k]
        if k in REJECTED_KEYS or k.startswith(REJECTED_PREFIXES):
            continue
        if isinstance(v, float):
            if math.isfinite(v):
                print(f"  static constexpr float {k} = {v:.9g}f;")
        elif isinstance(v, int):
            print(f"  static constexpr int {k} = {v};")
    print("};")


def eval_candidates(fam, cand_list, tier, collect, seed, stage):
    rows = []
    total = len(cand_list)
    started = time.time()
    best_score = float("inf")
    best_candidate = "none"

    print(
        f"STAGE_START family={fam} stage={stage} tier={tier} candidates={total} seed={seed}"
    )

    for i, (cid, raw_p) in enumerate(cand_list, 1):
        p = clean_params(raw_p)
        ok, reason = validate_candidate(p)
        if not ok:
            rr = {
                "family": fam,
                "stage": stage,
                "candidate": cid,
                "seed": seed,
                "tier": tier,
                "wave_dataset": "rejected_before_run",
                "returncode": -1,
                "quality_gate_pass": False,
                "fail_reason": "rejected_before_run",
                "rejected_before_run": 1,
                "reject_reason": reason,
                "score": float("inf"),
                "_params": dict(p),
                "input_params": dict(p),
            }
            rr.update(p)
            rows.append(rr)

            elapsed = time.time() - started
            eta = (elapsed / i) * (total - i) if i > 0 else 0.0
            print(
                f"PROGRESS family={fam} stage={stage} tier={tier} i={i} total={total} candidate={cid} seed={seed} "
                f"returncode=-1 rows=1 pass_rows=0 fail_rows=1 score=inf best={best_candidate} "
                f"best_score={best_score if math.isfinite(best_score) else float('inf'):.6g} "
                f"elapsed={format_duration(elapsed)} eta={format_duration(eta)}"
            )
            continue

        run_rows = run_candidate(fam, cid, p, tier, seed, collect, stage)

        # Score current batch (already-run rows + this candidate) for live progress.
        tmp = rows + run_rows
        score_rows(tmp)

        # keep candidate rows with updated scores
        run_rows = tmp[len(rows) :]
        rows.extend(run_rows)

        cand_rows = [r for r in run_rows if r.get("candidate") == cid]
        finite_scores = [r["score"] for r in cand_rows if math.isfinite(r.get("score", float("inf")))]
        cand_score = min(finite_scores) if finite_scores else float("inf")
        if cand_score < best_score:
            best_score = cand_score
            best_candidate = cid

        pass_rows = sum(1 for r in run_rows if r.get("quality_gate_pass"))
        fail_rows = len(run_rows) - pass_rows
        ret = cand_rows[0]["returncode"] if cand_rows else -1
        elapsed = time.time() - started
        eta = (elapsed / i) * (total - i) if i > 0 else 0.0

        print(
            f"PROGRESS family={fam} stage={stage} tier={tier} i={i} total={total} candidate={cid} seed={seed} "
            f"returncode={ret} rows={len(run_rows)} pass_rows={pass_rows} fail_rows={fail_rows} "
            f"score={cand_score:.6g} best={best_candidate} best_score={best_score:.6g} "
            f"elapsed={format_duration(elapsed)} eta={format_duration(eta)}"
        )

    return rows


def sample_local_perturbations(base_params, ranges, rng, n):
    out = []
    keys = [k for k in ranges if k in base_params]
    for i in range(n):
        p = dict(base_params)
        for k in keys:
            lo, hi = ranges[k]
            v = float(p[k])
            if lo > 0.0 and hi > 0.0:
                spread = 0.18 * (math.log(hi) - math.log(lo))
                nv = math.exp(math.log(v) + rng.uniform(-spread, spread))
            else:
                spread = 0.18 * (hi - lo)
                nv = v + rng.uniform(-spread, spread)
            p[k] = min(max(nv, lo), hi)
        if "SF_BOOT_GRAV_MIN_SEC" in p and "SF_BOOT_GRAV_HOLD_SEC" in p:
            apply_timeout_constraint(p, rng)
        out.append((f"local_{i:04d}", p))
    return out


def write_csv(path, rows):
    keys = set()
    for r in rows:
        keys.update(r.keys())
    fields = sorted(keys)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def warn_identical_to_baseline(agg_rows):
    base = next((x for x in agg_rows if x["candidate"] == "baseline" and x["valid"]), None)
    if not base:
        return
    nonbase = [x for x in agg_rows if x["candidate"] != "baseline" and x["valid"]]
    if nonbase and all(abs(x["mean_score"] - base["mean_score"]) < 1e-12 for x in nonbase):
        print("WARNING_ALL_SCORES_IDENTICAL_TO_BASELINE")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["gravity", "ou", "full"], default="gravity")
    ap.add_argument("--family", choices=["OU_II", "OU_III", "both"], default="both")
    ap.add_argument("--samples", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tier", choices=["quick", "all", "final"], default="quick")
    ap.add_argument("--top-k", type=int, default=8)
    ap.add_argument("--collect-all-gates", action="store_true", default=False)
    ap.add_argument("--out-csv", default="")
    a = ap.parse_args()

    print("BUILD_START family=OU_II")
    subprocess.run(["make", "-C", str(ROOT / "tests/kalman_ou_ii"), "build"], check=True)
    print("BUILD_DONE family=OU_II")

    print("BUILD_START family=OU_III")
    subprocess.run(["make", "-C", str(ROOT / "tests/kalman_ou_iii"), "build"], check=True)
    print("BUILD_DONE family=OU_III")

    fams = ["OU_II", "OU_III"] if a.family == "both" else [a.family]
    all_rows = []

    for fam in fams:
        rng = random.Random(a.seed + (0 if fam == "OU_II" else 1000000))

        # Stage B: baseline + minimal + practical + gravity MC
        grav_mc = [
            (f"grav_mc_{i:04d}", p)
            for i, p in enumerate(sample_params(GRAVITY_RANGES[fam], a.samples, rng, gravity=True), 1)
        ]
        stage_b = [("baseline", {})] + minimal_candidates() + explicit_practical_candidates() + grav_mc
        rows_b = eval_candidates(
            fam,
            stage_b,
            tier=a.tier,
            collect=(a.collect_all_gates or a.tier != "final"),
            seed=a.seed,
            stage="B",
        )
        score_rows(rows_b)
        all_rows.extend(rows_b)

        agg_b_all = aggregate_candidates(rows_b, fam, a.tier)
        print_stage_summary("B", agg_b_all)
        print_param_sensitivity(rows_b, fam)
        warn_identical_to_baseline(agg_b_all)

        agg_b = sorted([x for x in agg_b_all if x["valid"]], key=rank_key)
        nonbase_gravity = [x for x in agg_b if x["candidate"] != "baseline"]
        gravity_winners = nonbase_gravity[: max(1, a.top_k)]

        if a.mode in ("ou", "full"):
            gravity_winners = [{"candidate": "baseline_gravity", "params": {}, "valid": True}] + gravity_winners

        final_candidates = gravity_winners

        if a.mode in ("ou", "full"):
            # Stage C: baseline + OU combinations
            ou_space = dict(OU_COMMON)
            if fam == "OU_II":
                ou_space.update(OU_II_ONLY)
            sampled_ou = sample_params(ou_space, a.samples, rng)

            stage_c = []
            for gi, g in enumerate(gravity_winners, 1):
                basep = clean_params(g["params"])
                for oi, ou in enumerate(sampled_ou, 1):
                    stage_c.append((f"comb_g{gi:02d}_ou{oi:04d}", {**basep, **ou}))

            rows_c = eval_candidates(
                fam,
                [("baseline", {})] + stage_c,
                tier=a.tier,
                collect=(a.collect_all_gates or a.tier != "final"),
                seed=a.seed,
                stage="C",
            )
            score_rows(rows_c)
            all_rows.extend(rows_c)

            agg_c_all = aggregate_candidates(rows_c, fam, a.tier)
            print_stage_summary("C", agg_c_all)
            print_param_sensitivity(rows_c, fam)
            warn_identical_to_baseline(agg_c_all)

            agg_c = sorted([x for x in agg_c_all if x["valid"]], key=rank_key)
            top_c = agg_c[: max(1, a.top_k)]

            # Stage D: baseline + local refinement
            local = []
            local_ranges = dict(GRAVITY_RANGES[fam])
            local_ranges.update(ou_space)
            per_top = max(1, a.samples // max(1, a.top_k))
            for idx, c in enumerate(top_c, 1):
                local.extend(
                    [
                        (f"refine_t{idx:02d}_{cid}", p)
                        for cid, p in sample_local_perturbations(
                            clean_params(c["params"]), local_ranges, rng, per_top
                        )
                    ]
                )

            rows_d = eval_candidates(
                fam,
                [("baseline", {})] + local,
                tier=a.tier,
                collect=(a.collect_all_gates or a.tier != "final"),
                seed=a.seed,
                stage="D",
            )
            score_rows(rows_d)
            all_rows.extend(rows_d)

            agg_d_all = aggregate_candidates(rows_d, fam, a.tier)
            print_stage_summary("D", agg_d_all)
            print_param_sensitivity(rows_d, fam)
            warn_identical_to_baseline(agg_d_all)

            agg_d = sorted([x for x in agg_d_all if x["valid"]], key=rank_key)
            final_candidates = (agg_d[: max(1, a.top_k)] or top_c or gravity_winners)

        # Stage E: multi-seed baseline+final candidates per-seed
        validate_tier = a.tier if a.tier != "quick" else "final"

        final_candidates = [c for c in final_candidates if c.get("candidate") != "baseline"]

        if a.mode == "gravity":
            final_candidates = [c for c in final_candidates if c.get("params")]
            if not final_candidates:
                print("NO_NON_BASELINE_CANDIDATES_SURVIVED")

        if not final_candidates:
            print("NO_NON_BASELINE_CANDIDATES_SURVIVED")

        rows_e = []
        for s in (a.seed, a.seed + 1, a.seed + 2):
            seed_batch = [("baseline", {})]
            for i, c in enumerate(final_candidates, 1):
                seed_batch.append((f"final_{i:02d}", clean_params(c["params"])))

            seed_rows = eval_candidates(
                fam,
                seed_batch,
                tier=validate_tier,
                collect=True,
                seed=s,
                stage="E",
            )
            score_rows(seed_rows)
            rows_e.extend(seed_rows)

        all_rows.extend(rows_e)

        agg_e_all = aggregate_candidates(rows_e, fam, validate_tier)
        print_stage_summary("E", agg_e_all)
        print_param_sensitivity(rows_e, fam)
        warn_identical_to_baseline(agg_e_all)

        agg_e = sorted(
            [x for x in agg_e_all if x["valid"] and x["candidate"] != "baseline"],
            key=rank_key,
        )

        base_final = next((x for x in agg_e_all if x["candidate"] == "baseline" and x["valid"]), None)
        if base_final:
            print(f"BASELINE_SCORE {base_final['mean_score']:.6g}")

        if not agg_e:
            print(f"\n=== BEST_CONFIG {fam} ===")
            if base_final:
                print("NO_NON_BASELINE_CANDIDATES_SURVIVED")
            else:
                print("NO_VALID_CANDIDATES")
            continue

        best = agg_e[0]
        print(f"BEST_SCORE {best['mean_score']:.6g}")

        if base_final:
            ratio = best["mean_score"] / base_final["mean_score"] if base_final["mean_score"] > 0 else float("inf")
            print(f"SCORE_RATIO {ratio:.6g}")
            if ratio > 0.99:
                print("NO_MEANINGFUL_IMPROVEMENT")

        # ratio summaries
        best_rows = [r for r in rows_e if r.get("candidate") == best["candidate"]]
        rms3d_ratios = [r.get("rms3d_ratio", float("nan")) for r in best_rows if math.isfinite(float(r.get("rms3d_ratio", float("nan"))))]
        if rms3d_ratios:
            print(f"MEAN_RMS3D_RATIO {sum(rms3d_ratios) / len(rms3d_ratios):.6g}")
            print(f"MAX_RMS3D_RATIO {max(rms3d_ratios):.6g}")

        byds = defaultdict(list)
        for r in best_rows:
            rr = r.get("rms3d_ratio", float("nan"))
            if math.isfinite(float(rr)):
                byds[r["wave_dataset"]].append(rr)

        print("PER_DATASET_RATIOS")
        for ds, vals in sorted(byds.items()):
            print(f"{ds}: {sum(vals)/len(vals):.6g}")

        print(f"BEST_PARAMS {clean_params(best['params'])}")
        print_cpp_config(fam, best)

    if a.out_csv:
        write_csv(a.out_csv, all_rows)


if __name__ == "__main__":
    main()
