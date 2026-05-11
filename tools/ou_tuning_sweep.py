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

PARAM_RE = re.compile(r"^(SF_|OU_)")

# "disp", "balanced", or "yaw"
SCORE_MODE = "yaw"

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

# kind: "log", "linear", or "int"
BOOT_SPECS = {
    "SF_RACC_WARMUP_STD": (0.35, 0.65, "log"),
    "SF_ONLINE_TUNE_WARMUP_SEC": (4.0, 10.0, "log"),
    "SF_BOOT_TILT_ACC_TAU": (1.2, 2.8, "log"),
    "SF_BOOT_GRAV_SLOW_TAU": (4.5, 12.0, "log"),
    "SF_BOOT_GRAV_ALIGN_MAX_SIN": (0.08, 0.16, "linear"),
    "SF_BOOT_GRAV_HOLD_SEC": (0.8, 2.0, "log"),
    "SF_BOOT_GRAV_MIN_SEC": (2.0, 7.0, "linear"),
    "SF_BOOT_GRAV_NORM_FRAC": (0.18, 0.65, "linear"),
}

# OU_II common ranges kept close to the original productive OU_II sweep.
# Env names are generic OU_*.
OU_COMMON_SPECS = {
    "OU_TAU_COEFF": (1.15, 1.50, "log"),
    "OU_SIGMA_COEFF": (1.20, 1.80, "log"),
    "OU_ADAPT_TAU_SEC": (1.0, 2.2, "log"),
    "OU_ADAPT_EVERY_SECS": (0.05, 0.15, "log"),
    "OU_ACC_BIAS_INIT_STD": (0.02, 0.25, "log"),
}

# OU_III also uses generic OU_* env names, but with its own productive search ranges.
OU_III_COMMON_SPECS = {
    "OU_TAU_COEFF": (0.75, 1.75, "log"),
    "OU_SIGMA_COEFF": (0.75, 2.10, "log"),
    "OU_ADAPT_TAU_SEC": (1.0, 18.0, "log"),
    "OU_ADAPT_EVERY_SECS": (0.025, 0.250, "log"),
    "OU_ACC_BIAS_INIT_STD": (0.02, 0.35, "log"),
}

# OU_II was already fine. Keep its existing family-specific knobs.
OU_II_EXTRA_SPECS = {
    "OU_P_FACTOR": (1.20, 1.80, "log"),
    "OU_R_P0_XY_FACTOR": (0.20, 0.45, "log"),
    "OU_R_P0_COEFF": (1.20, 2.50, "log"),
    "OU_R_V0_COEFF": (1.20, 2.50, "log"),
}

# OU_III-specific analogues missing from the sweep.
OU_III_EXTRA_SPECS = {
    "OU_III_R_S_XY_FACTOR": (0.20, 0.65, "log"),
    "OU_III_R_S_COEFF": (0.70, 3.20, "log"),
}

# Optional. Off by default because this can hide OU problems by tuning mag behavior.
MAG_SPECS = {
    "SF_MAG_DELAY_SEC": (0.5, 10.0, "linear"),
    "SF_MAG_GRAV_ALIGN_MAX_SIN": (0.05, 0.31, "linear"),
    "SF_MAG_GRAV_ALIGN_HOLD_SEC": (0.15, 3.0, "log"),
    "SF_MAG_GRAV_ALIGN_LPF_TAU": (0.4, 8.0, "log"),
    "SF_MAG_TILT_FALLBACK_SEC": (1.0, 12.0, "log"),
    "SF_MAG_EXTREME_GYRO_DPS": (80.0, 480.0, "log"),
    "SF_MAG_MIN_SAMPLES": (4, 60, "int"),
}


def log(msg):
    print(msg, flush=True)


def format_duration(sec):
    sec = max(0, int(sec))
    h, rem = divmod(sec, 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def family_subdir(fam):
    return "kalman_ou_ii" if fam == "OU_II" else "kalman_ou_iii"


def family_bin(fam):
    return "./kalman_ou_ii-sim" if fam == "OU_II" else "./kalman_ou_iii-sim"


def build_family(fam):
    subdir = family_subdir(fam)
    log(f"BUILD_START family={fam}")
    subprocess.run(["make", "-C", str(ROOT / "tests" / subdir), "build"], check=True)
    log(f"BUILD_DONE family={fam}")


def make_space(fam, mode, tune_mag):
    specs = {}

    if mode in ("gravity", "full"):
        specs.update(BOOT_SPECS)

    if mode in ("ou", "full"):
        if fam == "OU_II":
            specs.update(OU_COMMON_SPECS)
            specs.update(OU_II_EXTRA_SPECS)
        else:
            specs.update(OU_III_COMMON_SPECS)
            specs.update(OU_III_EXTRA_SPECS)

    if tune_mag:
        specs.update(MAG_SPECS)

    return specs


def allowed_keys(fam, mode, tune_mag):
    keys = set(make_space(fam, mode, tune_mag).keys())
    keys.add("SF_BOOT_GRAV_TIMEOUT_SEC")
    return keys


def clean_params(p, fam=None, mode="full", tune_mag=False):
    if fam is None:
        allowed = None
    else:
        allowed = allowed_keys(fam, mode, tune_mag)

    out = {}
    for k, v in p.items():
        if not PARAM_RE.match(k):
            continue
        if allowed is not None and k not in allowed:
            continue
        if isinstance(v, float) and not math.isfinite(v):
            continue
        out[k] = v
    return out


def value_from_unit(spec, u):
    lo, hi, kind = spec
    u = min(max(u, 0.0), 1.0)

    if kind == "log":
        return math.exp(math.log(lo) + u * (math.log(hi) - math.log(lo)))

    if kind == "int":
        return int(round(lo + u * (hi - lo)))

    return lo + u * (hi - lo)


def lhs_values(rng, n):
    vals = [(i + rng.random()) / n for i in range(n)]
    rng.shuffle(vals)
    return vals


def sample_lhs(space, n, rng, require_timeout=False):
    if n <= 0:
        return []

    cols = {k: lhs_values(rng, n) for k in space}
    out = []

    for i in range(n):
        p = {k: value_from_unit(space[k], cols[k][i]) for k in space}
        if require_timeout and not apply_timeout_constraint(p, rng):
            continue
        out.append(p)

    attempts = 0
    max_attempts = max(1000, 60 * n)

    while len(out) < n and attempts < max_attempts:
        attempts += 1
        p = {k: value_from_unit(spec, rng.random()) for k, spec in space.items()}
        if require_timeout and not apply_timeout_constraint(p, rng):
            continue
        out.append(p)

    if len(out) < n:
        log(f"WARNING_SAMPLE_UNDERFILLED requested={n} produced={len(out)}")

    return out[:n]


def apply_timeout_constraint(p, rng):
    if "SF_BOOT_GRAV_MIN_SEC" not in p or "SF_BOOT_GRAV_HOLD_SEC" not in p:
        return True

    min_sec = float(p["SF_BOOT_GRAV_MIN_SEC"])
    hold_sec = float(p["SF_BOOT_GRAV_HOLD_SEC"])

    margin_max = min(MAX_MARGIN_SEC, MAX_BOOT_TIMEOUT_SEC - min_sec - hold_sec)
    if margin_max < MIN_MARGIN_SEC:
        return False

    margin = MIN_MARGIN_SEC + rng.random() * (margin_max - MIN_MARGIN_SEC)
    p["SF_BOOT_GRAV_TIMEOUT_SEC"] = min_sec + hold_sec + margin
    return p["SF_BOOT_GRAV_TIMEOUT_SEC"] <= MAX_BOOT_TIMEOUT_SEC


def finite_float(x, default=float("inf")):
    try:
        v = float(x)
    except Exception:
        return default
    return v if math.isfinite(v) else default


def coerce_candidate_to_space(p, space):
    out = {}

    for k, v in p.items():
        if k == "SF_BOOT_GRAV_TIMEOUT_SEC":
            out[k] = v
            continue

        if k not in space:
            continue

        lo, hi, kind = space[k]
        fv = finite_float(v)

        if not math.isfinite(fv):
            continue

        fv = min(max(fv, lo), hi)

        if kind == "int":
            out[k] = int(round(fv))
        else:
            out[k] = fv

    return out


def validate_candidate(p, space=None):
    bad = []

    if space is not None:
        for k, v in p.items():
            if k == "SF_BOOT_GRAV_TIMEOUT_SEC":
                continue

            if k not in space:
                bad.append(f"unknown_key:{k}")
                continue

            lo, hi, kind = space[k]
            fv = float(v)

            if fv < lo or fv > hi:
                bad.append(f"out_of_range:{k}={fv} not in [{lo},{hi}]")

            if kind == "int" and abs(fv - round(fv)) > 1e-9:
                bad.append(f"not_int:{k}={fv}")

    timeout = p.get("SF_BOOT_GRAV_TIMEOUT_SEC")
    if timeout is not None:
        timeout = float(timeout)

        if timeout > MAX_BOOT_TIMEOUT_SEC:
            bad.append(f"SF_BOOT_GRAV_TIMEOUT_SEC>{MAX_BOOT_TIMEOUT_SEC}")

        if "SF_BOOT_GRAV_MIN_SEC" in p and "SF_BOOT_GRAV_HOLD_SEC" in p:
            min_sec = float(p["SF_BOOT_GRAV_MIN_SEC"])
            hold_sec = float(p["SF_BOOT_GRAV_HOLD_SEC"])

            if timeout <= min_sec + hold_sec + MIN_MARGIN_SEC:
                bad.append("SF_BOOT_GRAV_TIMEOUT_SEC<=min+hold+margin")

    return len(bad) == 0, ";".join(bad)


def candidate_values_in_space(p, space):
    ok, _ = validate_candidate(p, space)
    return ok


def add_candidate_if_valid(out, cid, p, space):
    p = coerce_candidate_to_space(p, space)
    if candidate_values_in_space(p, space):
        out.append((cid, p))


def find_space_key(space, *names):
    for name in names:
        if name in space:
            return name
    return None


def probe_candidates(space, fam):
    out = []

    for k, spec in sorted(space.items()):
        out.append((f"probe_low_{k}", {k: value_from_unit(spec, 0.08)}))
        out.append((f"probe_mid_{k}", {k: value_from_unit(spec, 0.50)}))
        out.append((f"probe_high_{k}", {k: value_from_unit(spec, 0.92)}))

    if "SF_BOOT_GRAV_MIN_SEC" in space:
        base = {
            "SF_BOOT_TILT_ACC_TAU": 2.0,
            "SF_BOOT_GRAV_SLOW_TAU": 5.0,
            "SF_BOOT_GRAV_ALIGN_MAX_SIN": 0.14,
            "SF_BOOT_GRAV_HOLD_SEC": 1.0,
            "SF_BOOT_GRAV_MIN_SEC": 6.0,
            "SF_BOOT_GRAV_NORM_FRAC": 0.35,
            "SF_ONLINE_TUNE_WARMUP_SEC": 8.0,
            "SF_RACC_WARMUP_STD": 0.5,
        }
        for i, timeout in enumerate((8.5, 10.0, 12.0, 15.0), 1):
            p = dict(base)
            p["SF_BOOT_GRAV_TIMEOUT_SEC"] = timeout
            add_candidate_if_valid(out, f"practical_boot_{i:02d}", p, space)

    tau_key = find_space_key(space, "OU_TAU_COEFF")
    sigma_key = find_space_key(space, "OU_SIGMA_COEFF")
    bias_std_key = find_space_key(space, "OU_ACC_BIAS_INIT_STD")
    adapt_tau_key = find_space_key(space, "OU_ADAPT_TAU_SEC")
    adapt_every_key = find_space_key(space, "OU_ADAPT_EVERY_SECS")
    floor_key = find_space_key(space, "OU_ACC_NOISE_FLOOR_SIGMA")
    freq_key = find_space_key(space, "OU_FREQ_INPUT_CUTOFF_HZ")

    rs_xy_key = find_space_key(space, "OU_R_P0_XY_FACTOR", "OU_III_R_S_XY_FACTOR")
    rs_coeff_key = find_space_key(space, "OU_R_P0_COEFF", "OU_III_R_S_COEFF")
    rv0_coeff_key = find_space_key(space, "OU_R_V0_COEFF", "OU_III_R_V0_COEFF")

    if tau_key and sigma_key:
        practical = [
            {tau_key: 0.60, sigma_key: 0.60},
            {tau_key: 0.60, sigma_key: 1.60},
            {tau_key: 1.00, sigma_key: 0.75},
            {tau_key: 1.00, sigma_key: 1.60},
            {tau_key: 1.50, sigma_key: 0.60},
            {tau_key: 1.50, sigma_key: 1.60},
        ]

        if bias_std_key and adapt_tau_key and adapt_every_key:
            practical.extend(
                [
                    {
                        tau_key: 0.70,
                        sigma_key: 1.30,
                        bias_std_key: 0.25,
                        adapt_tau_key: 1.25,
                        adapt_every_key: 0.025,
                    },
                    {
                        tau_key: 1.50,
                        sigma_key: 0.75,
                        bias_std_key: 0.30,
                        adapt_tau_key: 8.0,
                        adapt_every_key: 0.20,
                    },
                    {
                        tau_key: 0.60,
                        sigma_key: 2.10,
                        bias_std_key: 0.35,
                    },
                    {
                        tau_key: 1.20,
                        sigma_key: 0.75,
                        bias_std_key: 0.35,
                    },
                ]
            )

        if floor_key and freq_key:
            practical.extend(
                [
                    {
                        tau_key: 1.00,
                        sigma_key: 1.25,
                        floor_key: 0.050,
                        freq_key: 0.35,
                    },
                    {
                        tau_key: 1.05,
                        sigma_key: 1.35,
                        floor_key: 0.065,
                        freq_key: 0.42,
                    },
                    {
                        tau_key: 1.15,
                        sigma_key: 1.55,
                        floor_key: 0.080,
                        freq_key: 0.55,
                    },
                ]
            )

        if rs_xy_key and rs_coeff_key and rv0_coeff_key:
            practical.extend(
                [
                    {
                        rs_xy_key: 0.25,
                        rs_coeff_key: 1.00,
                        rv0_coeff_key: 1.00,
                    },
                    {
                        rs_xy_key: 0.35,
                        rs_coeff_key: 1.40,
                        rv0_coeff_key: 1.30,
                    },
                    {
                        rs_xy_key: 0.50,
                        rs_coeff_key: 2.00,
                        rv0_coeff_key: 1.70,
                    },
                    {
                        tau_key: 1.05,
                        sigma_key: 1.35,
                        rs_xy_key: 0.35,
                        rs_coeff_key: 1.50,
                        rv0_coeff_key: 1.50,
                    },
                ]
            )

        for i, p in enumerate(practical, 1):
            add_candidate_if_valid(out, f"practical_ou_{i:02d}", p, space)

    if (
        fam == "OU_III"
        and tau_key
        and sigma_key
        and floor_key
        and adapt_tau_key
        and adapt_every_key
        and freq_key
        and bias_std_key
        and rs_xy_key
        and rs_coeff_key
        and rv0_coeff_key
    ):
        productive_ou3 = [
            {
                tau_key: 1.05,
                sigma_key: 1.33,
                floor_key: 0.063,
                adapt_tau_key: 11.8,
                adapt_every_key: 0.146,
                freq_key: 0.423,
                bias_std_key: 0.18,
                rs_xy_key: 0.35,
                rs_coeff_key: 1.50,
                rv0_coeff_key: 1.50,
            },
            {
                tau_key: 1.00,
                sigma_key: 1.25,
                floor_key: 0.055,
                adapt_tau_key: 10.0,
                adapt_every_key: 0.120,
                freq_key: 0.400,
                bias_std_key: 0.12,
                rs_xy_key: 0.30,
                rs_coeff_key: 1.25,
                rv0_coeff_key: 1.25,
            },
            {
                tau_key: 1.10,
                sigma_key: 1.45,
                floor_key: 0.070,
                adapt_tau_key: 14.0,
                adapt_every_key: 0.180,
                freq_key: 0.500,
                bias_std_key: 0.22,
                rs_xy_key: 0.45,
                rs_coeff_key: 1.80,
                rv0_coeff_key: 1.60,
            },
            {
                tau_key: 0.95,
                sigma_key: 1.20,
                floor_key: 0.050,
                adapt_tau_key: 8.0,
                adapt_every_key: 0.100,
                freq_key: 0.350,
                bias_std_key: 0.08,
                rs_xy_key: 0.25,
                rs_coeff_key: 1.10,
                rv0_coeff_key: 1.10,
            },
        ]

        for i, p in enumerate(productive_ou3, 1):
            add_candidate_if_valid(out, f"productive_ou3_{i:02d}", p, space)

    return out


def dedupe_candidates(cands):
    seen = set()
    out = []

    for cid, p in cands:
        key = tuple(
            sorted(
                (
                    k,
                    round(float(v), 9) if isinstance(v, (float, int)) else v,
                )
                for k, v in p.items()
            )
        )
        if key in seen:
            continue
        seen.add(key)
        out.append((cid, p))

    return out


def unique_by_params(cands):
    seen = set()
    out = []

    for c in cands:
        params = c.get("params", {})
        key = tuple(
            sorted(
                (
                    k,
                    round(float(v), 9) if isinstance(v, (float, int)) else v,
                )
                for k, v in params.items()
            )
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(c)

    return out


def parse(text):
    ds = RX["ds"].findall(text)
    ang = [m.groups() for m in RX["ang"].finditer(text)]
    xyz = [m.groups() for m in RX["xyz"].finditer(text)]
    r3d = [m.group(1) for m in RX["r3d"].finditer(text)]
    accb = [m.groups() for m in RX["accb"].finditer(text)]
    gates = [m.groups() for m in RX["gate"].finditer(text)]

    n = max(len(ds), len(ang), len(xyz), len(r3d), len(accb), len(gates), 1)
    out = []

    def get(arr, i, j):
        return float(arr[i][j]) if i < len(arr) else float("nan")

    for i in range(n):
        gate_pass = (gates[i][0] == "1") if i < len(gates) else False
        reason = gates[i][1].strip() if i < len(gates) else "missing_gate_status"

        out.append(
            {
                "wave_dataset": ds[i] if i < len(ds) else f"dataset_{i}",
                "roll_rms": get(ang, i, 0),
                "pitch_rms": get(ang, i, 1),
                "yaw_rms": get(ang, i, 2),
                "x_rms": get(xyz, i, 0),
                "y_rms": get(xyz, i, 1),
                "z_rms": get(xyz, i, 2),
                "rms_3d": float(r3d[i]) if i < len(r3d) else float("nan"),
                "acc_bias_rms_3d": get(accb, i, 3),
                "quality_gate_pass": gate_pass,
                "fail_reason": reason,
            }
        )

    return out


def synthetic_failure_row(fam, cid, p, tier, seed, stage, reason, returncode):
    row = {
        "family": fam,
        "candidate": cid,
        "seed": seed,
        "tier": tier,
        "stage": stage,
        "wave_dataset": "run_failed",
        "roll_rms": float("nan"),
        "pitch_rms": float("nan"),
        "yaw_rms": float("nan"),
        "x_rms": float("nan"),
        "y_rms": float("nan"),
        "z_rms": float("nan"),
        "rms_3d": float("nan"),
        "acc_bias_rms_3d": float("nan"),
        "quality_gate_pass": False,
        "fail_reason": reason,
        "returncode": returncode,
        "roll_pitch_rms_norm": float("nan"),
        "xy_rms": float("nan"),
        "_params": dict(p),
        "input_params": dict(p),
        "rejected_before_run": 0,
        "reject_reason": "",
        "score": float("inf"),
        "score_conservative": float("inf"),
        "score_disp": float("inf"),
        "score_yaw": float("inf"),
    }
    row.update(p)
    return row


def run_candidate(
    fam,
    cid,
    p,
    tier,
    seed,
    collect,
    stage,
    run_timeout_sec,
    mode,
    tune_mag,
):
    tdir = ROOT / "tests" / family_subdir(fam)
    bin_name = family_bin(fam)

    input_params = clean_params(
        p,
        fam=fam,
        mode=mode,
        tune_mag=tune_mag,
    )

    env = os.environ.copy()
    env["W3D_SEED"] = str(seed)
    env["W3D_TIER"] = tier

    if collect:
        env["W3D_COLLECT_ALL_GATES"] = "1"

    for k, v in input_params.items():
        if isinstance(v, float):
            env[k] = f"{v:.9g}"
        else:
            env[k] = str(v)

    try:
        # IMPORTANT: no --nomag. Magnetometer remains enabled.
        pr = subprocess.run(
            [bin_name],
            cwd=tdir,
            text=True,
            capture_output=True,
            env=env,
            timeout=run_timeout_sec if run_timeout_sec and run_timeout_sec > 0 else None,
        )
    except subprocess.TimeoutExpired:
        return [
            synthetic_failure_row(
                fam,
                cid,
                input_params,
                tier,
                seed,
                stage,
                f"subprocess_timeout_{run_timeout_sec}s",
                -2,
            )
        ]

    text = pr.stdout + "\n" + pr.stderr
    rows = parse(text)

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

        for k, v in input_params.items():
            r[k] = v

    return rows


def safe_ratio(num, den):
    n = finite_float(num, 1e-9)
    d = finite_float(den, 1e-9)
    if abs(d) < 1e-9:
        d = 1e-9
    if abs(n) < 1e-9:
        n = 1e-9
    return n / d


def percentile(vals, p):
    vals = [finite_float(v) for v in vals if math.isfinite(finite_float(v))]
    if not vals:
        return float("inf")

    vals.sort()
    pos = (len(vals) - 1) * p
    lo = int(pos)
    hi = min(len(vals) - 1, lo + 1)
    frac = pos - lo
    return vals[lo] * (1.0 - frac) + vals[hi] * frac


def score_rows(rows):
    by = defaultdict(list)

    for r in rows:
        if not r.get("rejected_before_run"):
            by[(r["family"], r["candidate"], r["seed"], r["tier"], r["stage"])].append(r)

    base = {
        (r["family"], r["wave_dataset"], r["seed"], r["tier"], r["stage"]): r
        for r in rows
        if r.get("candidate") == "baseline" and not r.get("rejected_before_run")
    }

    for _, items in by.items():
        rp = []
        z = []
        r3 = []
        yaw = []
        acc = []
        any_fail = False

        for r in items:
            b = base.get((r["family"], r["wave_dataset"], r["seed"], r["tier"], r["stage"]))

            if r["candidate"] != "baseline" and b is None:
                r["score"] = float("inf")
                r["score_conservative"] = float("inf")
                r["score_disp"] = float("inf")
                r["score_yaw"] = float("inf")
                r["quality_gate_pass"] = False
                fr = (r.get("fail_reason") or "").strip(";")
                r["fail_reason"] = (fr + ";missing_matching_baseline").strip(";")
                any_fail = True
                continue

            if b is None:
                b = r

            r["roll_pitch_ratio"] = safe_ratio(r.get("roll_pitch_rms_norm"), b.get("roll_pitch_rms_norm"))
            r["z_ratio"] = safe_ratio(r.get("z_rms"), b.get("z_rms"))
            r["rms3d_ratio"] = safe_ratio(r.get("rms_3d"), b.get("rms_3d"))
            r["yaw_ratio"] = safe_ratio(r.get("yaw_rms"), b.get("yaw_rms"))
            r["acc_bias_ratio"] = safe_ratio(r.get("acc_bias_rms_3d"), b.get("acc_bias_rms_3d"))

            rp.append(r["roll_pitch_ratio"])
            z.append(r["z_ratio"])
            r3.append(r["rms3d_ratio"])
            yaw.append(r["yaw_ratio"])
            acc.append(r["acc_bias_ratio"])

            if r.get("returncode") != 0 or not r.get("quality_gate_pass", False):
                any_fail = True

        if not rp or not z or not r3 or not yaw or not acc:
            conservative_score = float("inf")
            disp_score = float("inf")
            yaw_score = float("inf")
        else:
            conservative_score = (
                3.00 * percentile(r3, 0.75)
                + 2.00 * max(r3)
                + 2.50 * percentile(z, 0.75)
                + 1.50 * max(z)
                + 1.75 * percentile(rp, 0.75)
                + 1.25 * max(rp)
                + 0.70 * percentile(yaw, 0.75)
                + 0.40 * max(yaw)
                + 0.30 * percentile(acc, 0.75)
                + 0.25 * max(acc)
                + (100.0 if any_fail else 0.0)
            )

            disp_score = (
                8.00 * percentile(r3, 0.75)
                + 4.00 * percentile(z, 0.75)
                + 1.00 * percentile(rp, 0.75)
                + 0.50 * percentile(yaw, 0.75)
                + 0.10 * percentile(acc, 0.75)
                + 0.50 * max(r3)
                + 0.25 * max(z)
                + 0.20 * max(rp)
                + 0.10 * max(acc)
                + (100.0 if any_fail else 0.0)
            )

            yaw_score = (
                6.00 * percentile(yaw, 0.75)
                + 4.00 * max(yaw)
                + 2.00 * percentile(rp, 0.75)
                + 1.25 * max(rp)
                + 1.50 * percentile(r3, 0.75)
                + 0.75 * percentile(z, 0.75)
                + 0.25 * percentile(acc, 0.75)
                + 0.20 * max(acc)
                + (100.0 if any_fail else 0.0)
            )

        for r in items:
            r["score_conservative"] = conservative_score
            r["score_disp"] = disp_score
            r["score_yaw"] = yaw_score

            if SCORE_MODE == "balanced":
                r["score"] = conservative_score
            elif SCORE_MODE == "yaw":
                r["score"] = yaw_score
            else:
                r["score"] = disp_score


def aggregate_candidates(rows, fam, tier, stage=None):
    score_rows(rows)

    by = defaultdict(list)
    for r in rows:
        if r.get("family") != fam or r.get("tier") != tier:
            continue
        if stage is not None and r.get("stage") != stage:
            continue
        by[r["candidate"]].append(r)

    out = []

    for cand, items in by.items():
        valid = all(
            rr.get("returncode") == 0
            and rr.get("quality_gate_pass")
            and not rr.get("rejected_before_run")
            and math.isfinite(finite_float(rr.get("score")))
            and math.isfinite(finite_float(rr.get("rms_3d")))
            for rr in items
        )

        params = clean_params(items[0].get("_params", items[0].get("input_params", {})))

        obj = {
            "family": fam,
            "candidate": cand,
            "tier": tier,
            "stage": items[0].get("stage", "unknown"),
            "params": params,
            "valid": valid,
            "n_rows": len(items),
        }

        if valid:
            obj.update(
                {
                    "mean_score": sum(finite_float(r["score"]) for r in items) / len(items),
                    "max_score": max(finite_float(r["score"]) for r in items),
                    "mean_score_conservative": sum(finite_float(r.get("score_conservative")) for r in items) / len(items),
                    "mean_score_disp": sum(finite_float(r.get("score_disp")) for r in items) / len(items),
                    "mean_score_yaw": sum(finite_float(r.get("score_yaw")) for r in items) / len(items),
                    "mean_rms3d": sum(finite_float(r["rms_3d"]) for r in items) / len(items),
                    "mean_z_rms": sum(finite_float(r["z_rms"]) for r in items) / len(items),
                    "mean_roll_pitch_rms": sum(finite_float(r["roll_pitch_rms_norm"]) for r in items) / len(items),
                    "mean_yaw_rms": sum(finite_float(r["yaw_rms"]) for r in items) / len(items),
                    "mean_acc_bias_rms3d": sum(finite_float(r["acc_bias_rms_3d"]) for r in items) / len(items),
                    "mean_rms3d_ratio": sum(finite_float(r.get("rms3d_ratio")) for r in items) / len(items),
                    "max_rms3d_ratio": max(finite_float(r.get("rms3d_ratio")) for r in items),
                    "mean_z_ratio": sum(finite_float(r.get("z_ratio")) for r in items) / len(items),
                    "max_z_ratio": max(finite_float(r.get("z_ratio")) for r in items),
                    "mean_acc_bias_ratio": sum(finite_float(r.get("acc_bias_ratio")) for r in items) / len(items),
                    "max_acc_bias_ratio": max(finite_float(r.get("acc_bias_ratio")) for r in items),
                    "max_roll_pitch_ratio": max(finite_float(r.get("roll_pitch_ratio")) for r in items),
                    "max_yaw_ratio": max(finite_float(r.get("yaw_ratio")) for r in items),
                }
            )
        else:
            obj.update(
                {
                    "mean_score": float("inf"),
                    "max_score": float("inf"),
                    "mean_score_conservative": float("inf"),
                    "mean_score_disp": float("inf"),
                    "mean_score_yaw": float("inf"),
                    "mean_rms3d": float("inf"),
                    "mean_z_rms": float("inf"),
                    "mean_roll_pitch_rms": float("inf"),
                    "mean_yaw_rms": float("inf"),
                    "mean_acc_bias_rms3d": float("inf"),
                    "mean_rms3d_ratio": float("inf"),
                    "max_rms3d_ratio": float("inf"),
                    "mean_z_ratio": float("inf"),
                    "max_z_ratio": float("inf"),
                    "mean_acc_bias_ratio": float("inf"),
                    "max_acc_bias_ratio": float("inf"),
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
        a["mean_acc_bias_rms3d"],
        a["mean_roll_pitch_rms"],
        a["mean_yaw_rms"],
    )


def useful_nonbaseline(agg):
    return [
        x
        for x in agg
        if x.get("valid")
        and x.get("candidate") != "baseline"
        and clean_params(x.get("params", {}))
    ]


def build_diverse_pool(agg_a, agg_b, top_k):
    cands = useful_nonbaseline(agg_a) + useful_nonbaseline(agg_b)
    valid = [c for c in cands if c.get("valid")]

    if not valid:
        return []

    picks = []

    picks.extend(sorted(valid, key=rank_key)[:top_k])
    picks.extend(sorted(valid, key=lambda c: (c["mean_rms3d"], c["max_rms3d_ratio"]))[:top_k])
    picks.extend(sorted(valid, key=lambda c: (c["mean_z_rms"], c["max_z_ratio"]))[:top_k])
    picks.extend(sorted(valid, key=lambda c: (c["mean_yaw_rms"], c["max_yaw_ratio"]))[:max(2, top_k // 2)])

    improvers = [
        c
        for c in valid
        if c.get("mean_rms3d_ratio", float("inf")) < 0.985
        or c.get("mean_z_ratio", float("inf")) < 0.985
    ]
    picks.extend(sorted(improvers, key=lambda c: (c["mean_rms3d_ratio"], c["mean_z_ratio"]))[:top_k])

    return unique_by_params(picks)[: max(top_k, 16)]


def select_final_candidates(pool, top_k):
    if not pool:
        return []

    picks = []
    picks.extend(sorted(pool, key=rank_key)[:top_k])
    picks.extend(sorted(pool, key=lambda c: (c["mean_rms3d"], c["mean_z_rms"]))[:top_k])
    picks.extend(sorted(pool, key=lambda c: (c["mean_rms3d_ratio"], c["mean_z_ratio"]))[:top_k])
    picks.extend(sorted(pool, key=lambda c: (c["mean_yaw_rms"], c["max_yaw_ratio"]))[:max(2, top_k // 2)])

    return unique_by_params(picks)[: max(top_k, 12)]


def eval_candidates(
    fam,
    candidates,
    tier,
    collect,
    seed,
    stage,
    run_timeout_sec,
    mode,
    tune_mag,
    space,
):
    rows = []
    total = len(candidates)
    started = time.time()
    best_score = float("inf")
    best_candidate = "none"

    log(f"STAGE_START family={fam} stage={stage} tier={tier} candidates={total} seed={seed}")

    for i, (cid, raw_params) in enumerate(candidates, 1):
        p = clean_params(
            raw_params,
            fam=fam,
            mode=mode,
            tune_mag=tune_mag,
        )

        ok, reason = validate_candidate(p, space)

        if not ok:
            run_rows = [
                synthetic_failure_row(
                    fam,
                    cid,
                    p,
                    tier,
                    seed,
                    stage,
                    f"rejected_before_run:{reason}",
                    -1,
                )
            ]
            run_rows[0]["rejected_before_run"] = 1
            run_rows[0]["reject_reason"] = reason
        else:
            run_rows = run_candidate(
                fam=fam,
                cid=cid,
                p=p,
                tier=tier,
                seed=seed,
                collect=collect,
                stage=stage,
                run_timeout_sec=run_timeout_sec,
                mode=mode,
                tune_mag=tune_mag,
            )

        rows.extend(run_rows)
        score_rows(rows)

        cand_rows = [r for r in rows if r.get("candidate") == cid and r.get("stage") == stage]
        finite_scores = [finite_float(r.get("score")) for r in cand_rows if math.isfinite(finite_float(r.get("score")))]
        cand_score = min(finite_scores) if finite_scores else float("inf")

        if cand_score < best_score:
            best_score = cand_score
            best_candidate = cid

        pass_rows = sum(1 for r in run_rows if r.get("quality_gate_pass"))
        fail_rows = len(run_rows) - pass_rows
        ret = run_rows[0].get("returncode", -1) if run_rows else -1

        elapsed = time.time() - started
        eta = (elapsed / i) * (total - i) if i else 0.0

        score_str = f"{cand_score:.6g}" if math.isfinite(cand_score) else "inf"
        best_str = f"{best_score:.6g}" if math.isfinite(best_score) else "inf"

        extra = ""
        if not ok:
            extra = f" reject_reason={reason}"

        log(
            f"PROGRESS family={fam} stage={stage} tier={tier} i={i} total={total} "
            f"candidate={cid} seed={seed} returncode={ret} rows={len(run_rows)} "
            f"pass_rows={pass_rows} fail_rows={fail_rows} score={score_str} "
            f"best={best_candidate} best_score={best_str} "
            f"elapsed={format_duration(elapsed)} eta={format_duration(eta)}"
            f"{extra}"
        )

    score_rows(rows)
    return rows


def print_stage_summary(stage, agg):
    valid = [x for x in agg if x["valid"]]
    log(f"STAGE_{stage}_VALID_COUNT {len(valid)}")

    top = sorted(valid, key=rank_key)[:8]
    log(f"STAGE_{stage}_TOP {','.join(t['candidate'] for t in top) if top else 'none'}")

    for t in top[:5]:
        log(
            f"STAGE_{stage}_TOP_CAND candidate={t['candidate']} "
            f"score={t['mean_score']:.6g} "
            f"score_disp={t.get('mean_score_disp', float('nan')):.6g} "
            f"score_yaw={t.get('mean_score_yaw', float('nan')):.6g} "
            f"score_old={t.get('mean_score_conservative', float('nan')):.6g} "
            f"rms3d={t['mean_rms3d']:.6g} "
            f"z={t['mean_z_rms']:.6g} rp={t['mean_roll_pitch_rms']:.6g} "
            f"yaw={t['mean_yaw_rms']:.6g} acc_bias={t['mean_acc_bias_rms3d']:.6g} "
            f"r3ratio={t['mean_rms3d_ratio']:.6g}/{t['max_rms3d_ratio']:.6g} "
            f"zratio={t['mean_z_ratio']:.6g}/{t['max_z_ratio']:.6g} "
            f"accratio={t['mean_acc_bias_ratio']:.6g}/{t['max_acc_bias_ratio']:.6g} "
            f"params={t['params']}"
        )


def print_param_sensitivity(rows, fam):
    valid = [
        r
        for r in rows
        if r.get("family") == fam
        and r.get("candidate") != "baseline"
        and r.get("returncode") == 0
        and r.get("quality_gate_pass")
        and math.isfinite(finite_float(r.get("score")))
    ]

    param_data = defaultdict(lambda: {"x": [], "y": []})

    for r in valid:
        for k, v in r.get("input_params", {}).items():
            if not PARAM_RE.match(k):
                continue
            fv = finite_float(v)
            sc = finite_float(r.get("score"))
            if math.isfinite(fv) and math.isfinite(sc):
                param_data[k]["x"].append(fv)
                param_data[k]["y"].append(sc)

    log(f"PARAM_SENSITIVITY {fam}")
    log("param,n,min,max,corr_with_score")

    if not param_data:
        log("NONE,0,nan,nan,nan")
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
        c = corr(xs, ys)
        cstr = f"{c:.6g}" if math.isfinite(c) else "nan"
        log(f"{k},{len(xs)},{min(xs):.6g},{max(xs):.6g},{cstr}")

        if abs(max(xs) - min(xs)) < 1e-12:
            log(f"WARNING_PARAM_NO_VARIANCE {k}")


def sample_local_perturbations(base_params, space, rng, n):
    base_params = {k: v for k, v in base_params.items() if k in space}
    if not base_params:
        return []

    out = []

    for i in range(n):
        p = dict(base_params)

        for k, v in list(p.items()):
            lo, hi, kind = space[k]

            if kind == "int":
                width = max(1, int(round(0.20 * (hi - lo))))
                nv = int(round(float(v))) + rng.randint(-width, width)
                p[k] = int(min(max(nv, lo), hi))
                continue

            v = float(v)

            if kind == "log" and lo > 0.0 and hi > 0.0 and v > 0.0:
                spread = 0.22 * (math.log(hi) - math.log(lo))
                nv = math.exp(math.log(v) + rng.uniform(-spread, spread))
            else:
                spread = 0.22 * (hi - lo)
                nv = v + rng.uniform(-spread, spread)

            p[k] = min(max(nv, lo), hi)

        if "SF_BOOT_GRAV_MIN_SEC" in p and "SF_BOOT_GRAV_HOLD_SEC" in p:
            if not apply_timeout_constraint(p, rng):
                continue

        out.append((f"local_{i:04d}", p))

    return out


def print_cpp_config(fam, result):
    params = clean_params(result.get("params", {}))

    log(
        f"\n=== BEST_CONFIG {fam} candidate={result['candidate']} "
        f"tier={result['tier']} rows={result['n_rows']} ==="
    )
    log(
        f"RMS_SUMMARY mean_3d={result['mean_rms3d']:.6g} "
        f"mean_z={result['mean_z_rms']:.6g} "
        f"mean_roll_pitch={result['mean_roll_pitch_rms']:.6g} "
        f"mean_yaw={result['mean_yaw_rms']:.6g} "
        f"mean_acc_bias={result['mean_acc_bias_rms3d']:.6g}"
    )
    log(
        f"RATIO_SUMMARY mean_3d_ratio={result['mean_rms3d_ratio']:.6g} "
        f"max_3d_ratio={result['max_rms3d_ratio']:.6g} "
        f"mean_z_ratio={result['mean_z_ratio']:.6g} "
        f"max_z_ratio={result['max_z_ratio']:.6g} "
        f"mean_acc_bias_ratio={result['mean_acc_bias_ratio']:.6g} "
        f"max_acc_bias_ratio={result['max_acc_bias_ratio']:.6g}"
    )

    log("// C++ snippet:")
    log("struct TunedParams {")
    for k in sorted(params):
        v = params[k]
        if isinstance(v, float):
            if math.isfinite(v):
                log(f"  static constexpr float {k} = {v:.9g}f;")
        elif isinstance(v, int):
            log(f"  static constexpr int {k} = {v};")
        else:
            log(f'  static constexpr auto {k} = "{v}";')
    log("};")

    log("// Environment form:")
    for k in sorted(params):
        v = params[k]
        if isinstance(v, float):
            log(f"export {k}={v:.9g}")
        else:
            log(f"export {k}={v}")


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fields = set()
    for r in rows:
        fields.update(r.keys())

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=sorted(fields))
        writer.writeheader()
        writer.writerows(rows)

    log(f"WROTE_CSV {path}")


def final_report(rows_e, fam, tier):
    agg = aggregate_candidates(rows_e, fam, tier, stage="E")
    print_stage_summary("E", agg)
    print_param_sensitivity(rows_e, fam)

    base = next((x for x in agg if x["candidate"] == "baseline" and x["valid"]), None)
    nonbase = sorted(useful_nonbaseline(agg), key=rank_key)

    if base:
        log(f"BASELINE_SCORE {base['mean_score']:.6g}")
        log(f"BASELINE_SCORE_OLD {base.get('mean_score_conservative', float('nan')):.6g}")
        log(
            f"BASELINE_RMS mean_3d={base['mean_rms3d']:.6g} "
            f"mean_z={base['mean_z_rms']:.6g} "
            f"mean_rp={base['mean_roll_pitch_rms']:.6g} "
            f"mean_yaw={base['mean_yaw_rms']:.6g} "
            f"mean_acc_bias={base['mean_acc_bias_rms3d']:.6g}"
        )

    if not nonbase:
        log(f"\n=== BEST_CONFIG {fam} ===")
        log("NO_NON_BASELINE_CANDIDATES_SURVIVED" if base else "NO_VALID_CANDIDATES")
        return None

    best = nonbase[0]
    log(f"BEST_SCORE {best['mean_score']:.6g}")
    log(f"BEST_SCORE_OLD {best.get('mean_score_conservative', float('nan')):.6g}")

    if base and base["mean_score"] > 0:
        ratio = best["mean_score"] / base["mean_score"]
        log(f"SCORE_RATIO {ratio:.6g}")
        if ratio > 0.995:
            log("NO_MEANINGFUL_SCORE_IMPROVEMENT")

    best_rows = [
        r
        for r in rows_e
        if r.get("candidate") == best["candidate"] and r.get("stage") == "E"
    ]

    by_ds = defaultdict(list)
    for r in best_rows:
        rr = finite_float(r.get("rms3d_ratio"))
        zr = finite_float(r.get("z_ratio"))
        ar = finite_float(r.get("acc_bias_ratio"))
        if math.isfinite(rr):
            by_ds[r["wave_dataset"]].append((rr, zr, ar))

    log("PER_DATASET_RATIOS")
    for ds, vals in sorted(by_ds.items()):
        r3avg = sum(v[0] for v in vals) / len(vals)
        zavg = sum(v[1] for v in vals) / len(vals)
        aavg = sum(v[2] for v in vals) / len(vals)
        log(f"{ds}: rms3d_ratio={r3avg:.6g} z_ratio={zavg:.6g} acc_bias_ratio={aavg:.6g}")

    log(f"BEST_PARAMS {best['params']}")
    print_cpp_config(fam, best)
    return best


def main():
    global SCORE_MODE

    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["gravity", "ou", "full"], default="full")
    ap.add_argument("--family", choices=["OU_II", "OU_III", "both"], default="both")
    ap.add_argument("--samples", type=int, default=48)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tier", choices=["quick", "all", "final"], default="quick")
    ap.add_argument("--final-tier", choices=["quick", "all", "final"], default="")
    ap.add_argument("--top-k", type=int, default=8)
    ap.add_argument("--collect-all-gates", action="store_true", default=False)
    ap.add_argument("--out-csv", default="")
    ap.add_argument("--run-timeout-sec", type=float, default=0.0)
    ap.add_argument("--no-build", action="store_true", default=False)
    ap.add_argument("--tune-mag", action="store_true", default=False)
    ap.add_argument(
        "--score-mode",
        choices=["disp", "balanced", "yaw"],
        default="disp",
        help="Score objective. disp is displacement-focused, balanced is conservative, yaw prioritizes yaw stability.",
    )

    args = ap.parse_args()
    SCORE_MODE = args.score_mode

    fams = ["OU_II", "OU_III"] if args.family == "both" else [args.family]
    all_rows = []

    if not args.no_build:
        for fam in fams:
            build_family(fam)

    for fam in fams:
        log(
            f"\n===== FAMILY_START {fam} mode={args.mode} samples={args.samples} "
            f"seed={args.seed} tune_mag={args.tune_mag} "
            f"score_mode={args.score_mode} ====="
        )

        rng = random.Random(args.seed + (0 if fam == "OU_II" else 1000000))
        space = make_space(fam, args.mode, args.tune_mag)
        require_timeout = args.mode in ("gravity", "full")

        log(f"SEARCH_SPACE family={fam} keys={','.join(sorted(space.keys()))}")

        stage_a = dedupe_candidates([("baseline", {})] + probe_candidates(space, fam))
        rows_a = eval_candidates(
            fam=fam,
            candidates=stage_a,
            tier=args.tier,
            collect=(args.collect_all_gates or args.tier != "final"),
            seed=args.seed,
            stage="A",
            run_timeout_sec=args.run_timeout_sec,
            mode=args.mode,
            tune_mag=args.tune_mag,
            space=space,
        )
        all_rows.extend(rows_a)

        agg_a = aggregate_candidates(rows_a, fam, args.tier, stage="A")
        print_stage_summary("A", agg_a)
        print_param_sensitivity(rows_a, fam)

        mc = [
            (f"mc_{i:04d}", p)
            for i, p in enumerate(sample_lhs(space, args.samples, rng, require_timeout=require_timeout), 1)
        ]

        rows_b = eval_candidates(
            fam=fam,
            candidates=[("baseline", {})] + mc,
            tier=args.tier,
            collect=(args.collect_all_gates or args.tier != "final"),
            seed=args.seed,
            stage="B",
            run_timeout_sec=args.run_timeout_sec,
            mode=args.mode,
            tune_mag=args.tune_mag,
            space=space,
        )
        all_rows.extend(rows_b)

        agg_b = aggregate_candidates(rows_b, fam, args.tier, stage="B")
        print_stage_summary("B", agg_b)
        print_param_sensitivity(rows_b, fam)

        pool = build_diverse_pool(agg_a, agg_b, args.top_k)
        top_pool = pool

        if not top_pool:
            log("NO_STAGE_A_OR_B_NONBASELINE_SURVIVORS")
            continue

        local = []
        per_top = max(2, args.samples // max(1, len(top_pool)))

        for idx, cand in enumerate(top_pool, 1):
            samples = sample_local_perturbations(cand["params"], space, rng, per_top)
            for cid, p in samples:
                local.append((f"refine_t{idx:02d}_{cid}", p))

        local = dedupe_candidates(local)

        if local:
            rows_c = eval_candidates(
                fam=fam,
                candidates=[("baseline", {})] + local,
                tier=args.tier,
                collect=(args.collect_all_gates or args.tier != "final"),
                seed=args.seed,
                stage="C",
                run_timeout_sec=args.run_timeout_sec,
                mode=args.mode,
                tune_mag=args.tune_mag,
                space=space,
            )
            all_rows.extend(rows_c)

            agg_c = aggregate_candidates(rows_c, fam, args.tier, stage="C")
            print_stage_summary("C", agg_c)
            print_param_sensitivity(rows_c, fam)

            pool = unique_by_params(pool + useful_nonbaseline(agg_c))
        else:
            log("STAGE_C_SKIPPED no_local_candidates")

        final_candidates = select_final_candidates(pool, args.top_k)

        final_tier = args.final_tier or ("final" if args.tier == "quick" else args.tier)
        rows_e = []

        for s in (args.seed, args.seed + 1, args.seed + 2):
            batch = [("baseline", {})]
            for i, cand in enumerate(final_candidates, 1):
                batch.append((f"final_{i:02d}", cand["params"]))

            seed_rows = eval_candidates(
                fam=fam,
                candidates=batch,
                tier=final_tier,
                collect=True,
                seed=s,
                stage="E",
                run_timeout_sec=args.run_timeout_sec,
                mode=args.mode,
                tune_mag=args.tune_mag,
                space=space,
            )
            rows_e.extend(seed_rows)

        all_rows.extend(rows_e)
        final_report(rows_e, fam, final_tier)

        log(f"===== FAMILY_DONE {fam} =====\n")

    if args.out_csv:
        write_csv(args.out_csv, all_rows)


if __name__ == "__main__":
    main()
