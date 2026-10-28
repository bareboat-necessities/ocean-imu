#!/usr/bin/env python3
"""Tables, paired statistics and figures of the TFG-CARMA3 study (EXPERIMENTAL).

Reads reports/results/tfg_carma3/{raw_runs,followup_runs,dataset_runs,timing_runs}.csv,
the traces and math_checks.txt written by tools/tfg_carma3_benchmark.py and
tests/kalman_carma/tfg_carma3-test, and writes summary.md, paired_stats.json and
the SVG figures next to them.  Called by `python3 tools/tfg_carma3_benchmark.py report`.

Paired statistics, candidate against the reference (CARMA3 Q-MEKF, A41,
unless stated): change = 1 - mean(candidate)/mean(reference) over the same
seeds (positive = candidate lower), 95 % bootstrap interval over seeds
(10 000 resamples, fixed seed), win rate, two-sided Wilcoxon signed-rank p
with Holm correction inside each stated family.  "Equivalent" means the whole
95 % interval lies inside the +-1 % margin fixed before any run.
"""
from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "results" / "tfg_carma3"
J29 = "jonswap_fp0.29_hs0.8_g3.3"
PRIMARY = [J29, "pm_fp0.25_hs1", "mixed_swell_wind", "changing_fp"]
EQUIV_MARGIN = 1.0   # percent

LABEL = {
    "OU3": "OU-III (shipping)", "A41": "CARMA3 Q-MEKF", "T41": "TFG-CARMA3",
    "T41R0": "TFG-CARMA3, no reset transport", "T41RA": "TFG-CARMA3, attitude-only reset",
    "T41NW": "TFG-CARMA3, invariant-coordinate writes",
    "A52": "CARMA3 Q-MEKF, CARMA(5,2)", "T52": "TFG-CARMA3, CARMA(5,2)", "T52R0": "TFG-CARMA3 (5,2), no reset",
}
# Validated categorical order (dataviz reference palette, light mode).
COL = {"OU3": "#2a78d6", "A41": "#eb6834", "T41": "#1baf7a", "T41R0": "#eda100", "T41RA": "#e87ba4",
       "T41NW": "#008300", "A52": "#6250d6", "T52": "#e34948", "T52R0": "#8a8a86", "truth": "#0b0b0b"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

SC_WORDS = {"rel0": "following", "rel45": "stern-quartering", "rel90": "beam", "rel135": "bow-quartering",
            "rel180": "head", "U2.5": "underway 2.5 m/s", "anchor": "anchored", "drift": "drifting",
            "maneuver": "manoeuvring", "maggap": "300 s magnetometer gap", "maggap900": "900 s magnetometer gap",
            "T1800": "1800 s", "T3600": "3600 s", "maggap2400": "magnetometer lost at 1200 s", "jit": "+-25 % dt jitter", "jit1": "+-1 % dt jitter, 1 us clock",
            "turns": "90 deg turns in 6 s", "gbdrift": "gyro-bias drift", "temp": "temperature drift",
            "calbias": "calibrated bias"}
SC_BASE = {J29: "JONSWAP 0.29 Hz, 0.8 m", "pm_fp0.25_hs1": "PM 0.25 Hz, 1 m",
           "mixed_swell_wind": "mixed swell + wind", "changing_fp": "changing f_p", "bimodal": "bimodal",
           "broadband": "broadband", "abrupt_change": "abrupt change", "jonswap_fp0.20_hs1_g3.3": "JONSWAP 0.20 Hz, 1 m",
           "pm_fp0.15_hs2.5": "PM 0.15 Hz, 2.5 m", "stationary": "still IMU", "small_amplitude": "5 cm sea",
           "jonswap_fp0.40_hs0.02_g3.3": "2 cm sea (noise floor)", "large_rotation": "large roll/pitch + yaw",
           "tilted": "15/5 deg heel/trim", "bias_drift": "accel-bias drift"}


def sc_label(sc):
    base, *mods = sc.split("+")
    return ", ".join([SC_BASE.get(base, base)] + [SC_WORDS.get(m, m) for m in mods])


def fnum(r, k):
    try:
        v = float(r[k])
    except (KeyError, ValueError, TypeError):
        return math.nan
    return v


def read(path):
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def by_seed(rows, sc, flt, key):
    out = {}
    for r in rows:
        if r["scenario"] == sc and r["filter"] == flt and r.get("failed", "0") == "0":
            v = fnum(r, key)
            if math.isfinite(v):
                out[int(r["seed"])] = v
    return out


def bootstrap(a, b, n=10000, seed=20261028, stat="ratio"):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(a), size=(n, len(a)))
    if stat == "ratio":
        s = 100.0 * (1.0 - b[idx].mean(axis=1) / a[idx].mean(axis=1))
    else:
        s = (b[idx] - a[idx]).mean(axis=1)
    return [float(np.percentile(s, 2.5)), float(np.percentile(s, 97.5))]


def wilcoxon_p(a, b):
    from scipy.stats import wilcoxon
    d = b - a
    if len(d) < 6 or not np.any(d != 0):
        return math.nan
    return float(wilcoxon(a, b).pvalue)


def paired(rows, sc, ref, cand, key, stat="ratio"):
    ra, rb = by_seed(rows, sc, ref, key), by_seed(rows, sc, cand, key)
    seeds = sorted(set(ra) & set(rb))
    if len(seeds) < 2:
        return None
    a = np.array([ra[s] for s in seeds])
    b = np.array([rb[s] for s in seeds])
    res = {"n": len(seeds), "ref_mean": float(a.mean()), "cand_mean": float(b.mean()),
           "win_rate": float(np.mean(b < a)), "wilcoxon_p": wilcoxon_p(a, b)}
    if stat == "ratio":
        if a.mean() == 0:   # e.g. a recovery time that is zero for every seed
            res.update(change_pct=math.nan, ci_pct=[math.nan, math.nan], equivalent=False)
            return res
        res["change_pct"] = float(100.0 * (1.0 - b.mean() / a.mean()))
        res["ci_pct"] = bootstrap(a, b)
        lo, hi = res["ci_pct"]
        res["equivalent"] = bool(-EQUIV_MARGIN <= lo and hi <= EQUIV_MARGIN)
    else:
        res["diff"] = float((b - a).mean())
        res["ci"] = bootstrap(a, b, stat="diff")
    return res


def holm(results):
    """Holm-Bonferroni over the dicts in `results` (in place: p_holm)."""
    ps = [(r["wilcoxon_p"], i) for i, r in enumerate(results) if r and math.isfinite(r["wilcoxon_p"])]
    ps.sort()
    m = len(ps)
    running = 0.0
    for k, (p, i) in enumerate(ps):
        running = max(running, min(1.0, (m - k) * p))
        results[i]["p_holm"] = running
    for r in results:
        if r is not None and "p_holm" not in r:
            r["p_holm"] = math.nan


def pct(x):
    return f"{x:+.2f}" if abs(x) >= 0.1 else f"{x:+.3f}"


def fmt_change(ps):
    """Change with its interval; `equiv.` (interval inside +-1 %) and the
    Holm-corrected Wilcoxon p when it is below 0.05 are reported independently:
    a difference can be statistically significant and practically negligible."""
    if not ps or not math.isfinite(ps["change_pct"]):
        return "--"
    lo, hi = ps["ci_pct"]
    s = f"{pct(ps['change_pct'])} % [{pct(lo)}, {pct(hi)}]"
    tags = []
    if ps.get("equivalent"):
        tags.append("equiv.")
    p = ps.get("p_holm", ps["wilcoxon_p"])
    if math.isfinite(p) and p < 0.05:
        tags.append(f"p={p:.1g}")
    return s + (" " + ", ".join(tags) if tags else "")


def mean_of(rows, sc, flt, key):
    v = list(by_seed(rows, sc, flt, key).values())
    return float(np.mean(v)) if v else math.nan


def g(v, d=4):
    return "--" if not math.isfinite(v) else f"{v:.{d}g}"


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

CORE = [("rms_z", "heave RMS (m)"), ("rms_x", "surge RMS (m)"), ("rms_y", "sway RMS (m)"),
        ("rms_h", "horizontal RMS (m)"), ("rms_3d", "3-D RMS (m)"), ("rmsv_z", "vertical velocity RMS (m/s)"),
        ("rmsv_h", "horizontal velocity RMS (m/s)"), ("lf_x", "surge 60 s mean (m)"), ("lf_z", "heave 60 s mean (m)"),
        ("peak_z", "peak heave error (m)"), ("peak_x", "peak surge error (m)"),
        ("roll_rms_deg", "roll RMS (deg)"), ("pitch_rms_deg", "pitch RMS (deg)"), ("yaw_rms_deg", "yaw RMS (deg)"),
        ("yaw_w_max_deg", "max heading excursion (deg)"), ("bg_rms_dph", "gyro-bias RMS (deg/h)"),
        ("ba_rms_mg", "accel-bias RMS (mg)"), ("bg_f60_dph", "gyro-bias, last 60 s (deg/h)"),
        ("ba_f60_mg", "accel-bias, last 60 s (mg)")]
STARTUP = [("conv_z", "heave convergence (s)"), ("conv_x", "surge convergence (s)"), ("su_rms_z", "startup heave RMS (m)"),
           ("su_rms_x", "startup surge RMS (m)"), ("su_peak_z", "startup peak heave (m)"),
           ("su_tilt_rms_deg", "startup tilt RMS (deg)"), ("su_yaw_rms_deg", "startup yaw RMS (deg)")]
CONSIST = [("nees_pos3", "position NEES (3)"), ("nees_pv6", "p,v NEES (6)"), ("att_nees", "attitude NEES (3)"),
           ("yaw_nees", "yaw NEES (1)"), ("nees9_phys", "joint NEES, common coords (9)"),
           ("nees9_native", "joint NEES, own coords (9)"), ("att_nees_in95", "attitude NEES in 95 % (0.95)"),
           ("yaw_nees_in95", "yaw NEES in 95 % (0.95)"), ("cov95_x", "surge 95 % coverage"),
           ("cov95_z", "heave 95 % coverage"), ("nis", "accel NIS mean (3)"), ("nis_in95", "accel NIS in 95 % band"),
           ("cond_max", "max corr. condition number")]


def table_primary(rows, stats):
    lines = ["## 1. Principal seas (100 paired seeds, 1200 s, scored from 300 s)", "",
             "Change of TFG-CARMA3 against CARMA3 Q-MEKF: 1 - mean(TFG)/mean(MEKF), positive = TFG lower, "
             "95 % bootstrap interval; `equiv.` = interval inside +-1 %; p = Holm-corrected Wilcoxon "
             "(family: the 4 seas x the metrics of this table).", ""]
    fam = []
    for sc in PRIMARY:
        for k, _ in CORE:
            ps = paired(rows, sc, "A41", "T41", k)
            stats.setdefault("primary", {}).setdefault(sc, {})[k] = ps
            fam.append(ps)
    holm(fam)
    for sc in PRIMARY:
        n = len(by_seed(rows, sc, "T41", "rms_z"))
        lines += [f"### {sc_label(sc)} (n = {n})", "",
                  "| metric | OU-III | CARMA3 Q-MEKF | TFG-CARMA3 | TFG vs Q-MEKF | TFG win rate |", "|---|---|---|---|---|---|"]
        for k, lab in CORE:
            ps = stats["primary"][sc][k]
            lines.append(f"| {lab} | {g(mean_of(rows, sc, 'OU3', k))} | {g(mean_of(rows, sc, 'A41', k))} | "
                         f"{g(mean_of(rows, sc, 'T41', k))} | {fmt_change(ps)} | "
                         f"{(100 * ps['win_rate']):.0f} % |" if ps else f"| {lab} | -- | -- | -- | -- | -- |")
        lines.append("")
    return lines


def all_scenarios(rows):
    seen = []
    for r in rows:
        if r["scenario"] not in seen:
            seen.append(r["scenario"])
    return seen


def table_matrix(rows, stats):
    lines = ["## 2. Every condition", "",
             "Mean over seeds; change of TFG-CARMA3 against CARMA3 Q-MEKF as in section 1 (Holm family: every "
             "condition x the four metrics).", "",
             "| condition | n | horizontal RMS, m: OU-III / Q-MEKF / TFG | TFG vs Q-MEKF | heave RMS, m: OU-III / Q-MEKF / TFG | TFG vs Q-MEKF | 3-D RMS change | yaw RMS, deg: Q-MEKF / TFG | yaw change |",
             "|---|---|---|---|---|---|---|---|---|"]
    fam = []
    scs = all_scenarios(rows)
    for sc in scs:
        for k in ("rms_h", "rms_z", "rms_3d", "yaw_rms_deg"):
            ps = paired(rows, sc, "A41", "T41", k)
            stats.setdefault("matrix", {}).setdefault(sc, {})[k] = ps
            fam.append(ps)
    holm(fam)
    for sc in scs:
        m = stats["matrix"][sc]
        n = m["rms_z"]["n"] if m["rms_z"] else 0
        h = "/".join(g(mean_of(rows, sc, f, "rms_h"), 3) for f in ("OU3", "A41", "T41"))
        z = "/".join(g(mean_of(rows, sc, f, "rms_z"), 3) for f in ("OU3", "A41", "T41"))
        y = "/".join(g(mean_of(rows, sc, f, "yaw_rms_deg"), 3) for f in ("A41", "T41"))
        lines.append(f"| {sc_label(sc)} | {n} | {h} | {fmt_change(m['rms_h'])} | {z} | {fmt_change(m['rms_z'])} | "
                     f"{fmt_change(m['rms_3d'])} | {y} | {fmt_change(m['yaw_rms_deg'])} |")
    lines.append("")
    return lines


def table_matrix_att(rows, stats):
    lines = ["## 2b. Attitude, biases and attitude consistency in every condition", "",
             "Tilt = |(roll, pitch)| error of `Log(Rhat R^T)`.  Change of TFG-CARMA3 against CARMA3 Q-MEKF as in "
             "section 1 (Holm family: every condition x the three changes); attitude NEES (3 dof, ideal 3) as "
             "means.", "",
             "| condition | tilt RMS, deg: Q-MEKF / TFG | tilt change | accel-bias RMS change | gyro-bias RMS change | "
             "attitude NEES: Q-MEKF / TFG / TFG no reset |",
             "|---|---|---|---|---|---|"]
    fam = []
    scs = all_scenarios(rows)
    for sc in scs:
        for k in ("tilt_rms_deg", "ba_rms_mg", "bg_rms_dph"):
            ps = paired(rows, sc, "A41", "T41", k)
            stats.setdefault("matrix_att", {}).setdefault(sc, {})[k] = ps
            fam.append(ps)
    holm(fam)
    for sc in scs:
        m = stats["matrix_att"][sc]
        t = "/".join(g(mean_of(rows, sc, f, "tilt_rms_deg"), 3) for f in ("A41", "T41"))
        n = "/".join(g(mean_of(rows, sc, f, "att_nees"), 3) for f in ("A41", "T41", "T41R0"))
        lines.append(f"| {sc_label(sc)} | {t} | {fmt_change(m['tilt_rms_deg'])} | {fmt_change(m['ba_rms_mg'])} | "
                     f"{fmt_change(m['bg_rms_dph'])} | {n} |")
    lines.append("")
    return lines


MECH_KEYS = [("yaw_gap_rms_deg", "yaw RMS in the gap (deg)"), ("yaw_gap_nees", "yaw NEES in the gap (1)"),
             ("pred_yaw_gap_end_deg", "predicted yaw std at gap end (deg)"),
             ("yaw_post_rms_deg", "yaw RMS, 60 s after the gap (deg)"), ("yaw_nees", "yaw NEES, scored window (1)"),
             ("pitch_rms_deg", "pitch RMS (deg)"), ("tilt_rms_deg", "tilt RMS (deg)"),
             ("att_nees", "attitude NEES (3)"), ("ba_rms_mg", "accel-bias RMS (mg)"),
             ("rms_h", "horizontal RMS (m)"), ("rms_z", "heave RMS (m)")]


def table_mechanism(mech, stats):
    if not mech:
        return []
    filters = ["OU3", "A41", "T41", "T41R0", "T41RA"]
    lines = ["## 11. Post hoc mechanism tests (30 paired seeds each)", "",
             "Designed after section 1 showed the pair separating in the mixed sea and every filter's attitude "
             "NEES far above 3; not part of the grid fixed before the runs.  Mixed sea with a 300 s magnetometer "
             "gap: without the magnetometer only the accelerometer can inform heading, which is unobservable "
             "from it, so heading information gained in the gap is false observability.  JONSWAP with a "
             "calibrated accelerometer (+-1 mg turn-on bias instead of up to 8 mg; the simulated bias random walk "
             "is unchanged): the attitude consistency when the simulated bias comes closer to the filters' "
             "accelerometer-bias prior (0.4 mg).  Changes against "
             "CARMA3 Q-MEKF as in section 1 (Holm family: the ratio rows of this section); NEES rows as "
             "differences.", ""]
    out, fam = {}, []
    for sc in all_scenarios(mech):
        for k, _ in MECH_KEYS:
            stat = "diff" if "nees" in k else "ratio"
            for c in ("T41", "T41R0"):
                ps = paired(mech, sc, "A41", c, k, stat=stat)
                out.setdefault(sc, {}).setdefault(c, {})[k] = ps
                if stat == "ratio":
                    fam.append(ps)
    holm(fam)
    stats["mechanism"] = out
    for sc in all_scenarios(mech):
        lines += [f"### {sc_label(sc)} (n = {len(by_seed(mech, sc, 'T41', 'rms_z'))})", "",
                  "| metric | " + " | ".join(LABEL[f] for f in filters) + " | TFG vs Q-MEKF | TFG no reset vs Q-MEKF |",
                  "|---|" + "---|" * (len(filters) + 2)]
        for k, lab in MECH_KEYS:
            if not any(math.isfinite(mean_of(mech, sc, f, k)) for f in filters):
                continue   # a gap metric in a run without a gap
            cells = [g(mean_of(mech, sc, f, k)) for f in filters]
            for c in ("T41", "T41R0"):
                ps = out[sc][c][k]
                if ps and "change_pct" in ps:
                    cells.append(fmt_change(ps))
                elif ps:
                    cells.append(f"{ps['diff']:+.3g} [{ps['ci'][0]:+.3g}, {ps['ci'][1]:+.3g}]")
                else:
                    cells.append("--")
            lines.append(f"| {lab} | " + " | ".join(cells) + " |")
        lines.append("")
    return lines


ABL = ["T41", "T41R0", "T41RA", "T41NW"]
GAPS = [f"{J29}+maggap", f"{J29}+T1800+maggap900", f"{J29}+T3600+maggap2400"]
ABL_SC = PRIMARY + ["broadband", "bimodal"] + GAPS + ["stationary", f"{J29}+turns", "large_rotation"]


def table_ablations(rows, stats):
    lines = ["## 3. Ablations against CARMA3 Q-MEKF", "",
             "`T41R0` no reset transport, `T41RA` J_l(-phi) on the attitude block only, `T41NW` CARMA3's covariance "
             "writes applied directly in the invariant coordinates.", ""]
    for key, lab in (("rms_h", "horizontal RMS"), ("rms_z", "heave RMS"), ("yaw_rms_deg", "yaw RMS"),
                     ("tilt_rms_deg", "tilt RMS"), ("att_nees", "attitude NEES (mean)"), ("yaw_nees", "yaw NEES (mean)")):
        stat = "diff" if "nees" in key else "ratio"
        lines += [f"### {lab}", "", "| condition | Q-MEKF | " + " | ".join(LABEL[a] for a in ABL) + " |",
                  "|---|---|" + "---|" * len(ABL)]
        for sc in ABL_SC:
            if not by_seed(rows, sc, "A41", key):
                continue
            cells = []
            for a in ABL:
                ps = paired(rows, sc, "A41", a, key, stat=stat)
                stats.setdefault("ablation", {}).setdefault(sc, {}).setdefault(a, {})[key] = ps
                if not ps:
                    cells.append("--")
                elif stat == "ratio":
                    cells.append(f"{g(ps['cand_mean'], 4)} ({fmt_change(ps)})")
                else:
                    lo, hi = ps["ci"]
                    cells.append(f"{g(ps['cand_mean'], 3)} ({ps['diff']:+.3g} [{lo:+.3g}, {hi:+.3g}])")
            lines.append(f"| {sc_label(sc)} | {g(mean_of(rows, sc, 'A41', key), 4)} | " + " | ".join(cells) + " |")
        lines.append("")
    return lines


def table_values(rows, title, metrics, scs, filters, note=""):
    lines = [title, ""] + ([note, ""] if note else [])
    for sc in scs:
        if not by_seed(rows, sc, "A41", metrics[0][0]) and not by_seed(rows, sc, "A41", "rms_z"):
            continue
        lines += [f"### {sc_label(sc)}", "", "| metric | " + " | ".join(LABEL[f] for f in filters) + " |",
                  "|---|" + "---|" * len(filters)]
        for k, lab in metrics:
            lines.append(f"| {lab} | " + " | ".join(g(mean_of(rows, sc, f, k)) for f in filters) + " |")
        lines.append("")
    return lines


def table_maggap(rows, stats):
    keys = [("yaw_gap_rms_deg", "yaw RMS in the gap (deg)"), ("yaw_gap_nees", "yaw NEES in the gap (1)"),
            ("pred_yaw_gap_end_deg", "predicted yaw std at gap end (deg)"),
            ("yaw_post_rms_deg", "yaw RMS, 60 s after the gap (deg)"), ("yaw_recover_s", "yaw recovery time (s)"),
            ("rms_h", "horizontal RMS (m)"), ("mag_forwarded", "magnetometer samples forwarded"),
            ("mag_applied", "magnetometer updates applied")]
    lines = table_values(rows, "## 5. Magnetometer outage and reacquisition", keys,
                         GAPS, ["OU3", "A41", "T41", "T41R0", "T41RA", "T41NW"],
                         "Gap from T/3; recovery = first time after the gap from which |yaw error| stays within "
                         "max(0.5 deg, 2 x pre-gap RMS) for 30 s.  Rejected = forwarded - applied.")
    for sc in GAPS:
        for k in ("yaw_gap_rms_deg", "yaw_post_rms_deg"):
            stats.setdefault("maggap", {}).setdefault(sc, {})[k] = paired(rows, sc, "A41", "T41", k)
        stats["maggap"][sc]["yaw_gap_nees"] = paired(rows, sc, "A41", "T41", "yaw_gap_nees", stat="diff")
    return lines


def table_followup(rows, stats):
    if not rows:
        return []
    import carma3_analysis as ca
    cost = {"A52": ca.cost_model(24, [5, 5, 5], 9)["total_madds_per_sample"],
            "T52": tfg_cost_model(5, "full")["total_madds_per_sample"],
            "T52R0": tfg_cost_model(5, "none")["total_madds_per_sample"]}
    stats["cost_model_followup"] = cost
    lines = ["## 7. Follow-up: CARMA(5,2) on all three axes (30 paired seeds)", "",
             f"24 error states.  Multiply-adds per sample: CARMA3 (5,2) {cost['A52']}, TFG (5,2) {cost['T52']}, "
             f"TFG (5,2) without reset transport {cost['T52R0']}.", "",
             "| condition | metric | CARMA3 (4,1) | TFG (4,1) | CARMA3 (5,2) | TFG (5,2) | TFG (5,2) vs Q-MEKF (5,2) | TFG (5,2), no reset |",
             "|---|---|---|---|---|---|---|---|"]
    for sc in all_scenarios(rows):
        for k, lab in (("rms_h", "horizontal RMS (m)"), ("rms_z", "heave RMS (m)"), ("tilt_rms_deg", "tilt RMS (deg)"),
                       ("yaw_rms_deg", "yaw RMS (deg)"), ("att_nees", "attitude NEES")):
            ps = paired(rows, sc, "A52", "T52", k, stat="diff" if "nees" in k else "ratio")
            stats.setdefault("followup", {}).setdefault(sc, {})[k] = ps
            ch = fmt_change(ps) if ps and "change_pct" in ps else (
                f"{ps['diff']:+.3g} [{ps['ci'][0]:+.3g}, {ps['ci'][1]:+.3g}]" if ps else "--")
            lines.append(f"| {sc_label(sc)} | {lab} | " + " | ".join(g(mean_of(rows, sc, f, k))
                         for f in ("A41", "T41", "A52", "T52")) + f" | {ch} | {g(mean_of(rows, sc, 'T52R0', k))} |")
    lines.append("")
    return lines


def table_dataset(rows):
    if not rows:
        return []
    fams = ["OU3", "A41", "T41", "T41R0"]
    keys = [("disp_z_pct_hs", "heave RMS (% H_s)"), ("disp_x_pct_hs", "surge RMS (% H_s)"),
            ("disp_y_pct_hs", "sway RMS (% H_s)"), ("yaw_rms_deg", "yaw RMS (deg)"),
            ("roll_rms_deg", "roll RMS (deg)"), ("pitch_rms_deg", "pitch RMS (deg)")]
    keys = [k for k in keys if k[0] in rows[0]]
    lines = ["## 8. Independent vessel-RAO dataset (8 records x 10 seeds)", "",
             "| record | " + " | ".join(f"{lab}: " + " / ".join(fams) for _, lab in keys) + " |",
             "|---|" + "---|" * len(keys)]
    inputs = sorted({r["input"] for r in rows})
    for inp in inputs:
        cells = []
        for k, _ in keys:
            vals = []
            for f in fams:
                v = [fnum(r, k) for r in rows if r["input"] == inp and r["family"] == f]
                v = [x for x in v if math.isfinite(x)]
                vals.append(g(float(np.mean(v)), 4) if v else "--")
            cells.append(" / ".join(vals))
        lines.append(f"| {inp.replace('wave_data_', '').replace('.csv', '')} | " + " | ".join(cells) + " |")
    lines.append("")
    return lines


def tfg_cost_model(nc=4, reset="full"):
    """Multiply-adds of one prediction plus one accelerometer update of
    Kalman3D_Wave_TFG_CARMA3, loop by loop, in the convention of
    carma3_analysis.cost_model (one multiply-add per product; the 3 x 3
    attitude/bias algebra of a step counted, Exp/J_l as one 3 x 3 product)."""
    nx = 9 + 3 * nc
    nblk = nx // 3
    pairs = nx * (nx + 1) // 2
    column = 9 + 9 + nc * 3 * (3 + nc) + 12              # Phi . column: phi, beta_g, the world block, beta_a
    predict = (2 * nx * column                            # Phi P Phi^T, two passes over full P
               + 2 * nc * 27 + nc * nc * (27 + 9)         # Qd world blocks [w_i]x Q_pp/Q_pb, [w_i]x Q_pp [w_j]x^T + Q_A Sigma_s
               + 5 * (27 + 27 + 54) + 3 * 54              # GL-5 nodes (J_l, C(s), C Q C^T), Q_pb, Q_bb, Q_aa
               + 5 * 27 + 3 * nc * nc)                    # step blocks (Exp, R1, E, F), mean W Phi_A^T
    h_nnz = 8                                             # [g]x (2), -I on a (3), -I on beta_a (3)
    correct = (3 * 9 + nx * 9 + nx * 3                    # S^-1 (three solves), K, a = K r
               + nx * 9 + pairs * 6                       # G = K S - P H^T; Joseph, 3 columns
               + 4 * 27 + nx * 3)                         # Exp_G(-a) o Xhat: Exp, J_l, R W, J_l rho, biases
    if reset == "full":
        correct += 2 * (nx * 9 + (nblk - 1) * nx * 18) + (nblk - 1) * 3 * 27 + 27
    elif reset == "attitude":
        correct += 2 * nblk * nx * 9 + 27                 # the implementation applies the identity blocks too
    acc = nx * h_nnz + 3 * h_nnz + correct
    mag = nx * 9 + 27 + correct
    return {"states": nx, "predict_madds": predict, "acc_update_madds": acc, "mag_update_madds": mag,
            "total_madds_per_sample": predict + acc}


def cost_models():
    import carma3_analysis as ca
    return {"OU3": ca.cost_model(21, [12], 9, "ou3"), "A41": ca.cost_model(21, [4, 4, 4], 9),
            "T41": tfg_cost_model(4, "full"), "T41R0": tfg_cost_model(4, "none"),
            "T41RA": tfg_cost_model(4, "attitude")}


def table_timing(rows, stats):
    cost = cost_models()
    stats["cost_model"] = cost
    lines = ["## 9. Cost", "",
             "Multiply-adds: one prediction plus one accelerometer update, analytic loop counts of each "
             "implementation (`cost_model` in tools/carma3_analysis.py for OU-III and CARMA3, `tfg_cost_model` "
             "here).  Host time: x86-64, -O3, one process on an otherwise idle host, mean update after Live over "
             "three JONSWAP seeds; it ranks the filters and is not an ESP32-S3 time.", "",
             "| filter | prediction | accel update | total multiply-adds | host update, fixed dt (us) | "
             "+-1 % on a 1 us clock (us) | +-25 % jitter (us) | wrapper + filter object, host (B) |",
             "|---|---|---|---|---|---|---|---|"]
    for f in ("OU3", "A41", "T41", "T41R0", "T41RA"):
        c = cost[f]
        cells = [str(c["predict_madds"]), str(c["acc_update_madds"]), str(c["total_madds_per_sample"])]
        for clk in ("fixed", "jit1", "jit"):
            v = [fnum(r, "runtime_ns") for r in rows if r["filter"] == f and r["clock"] == clk]
            cells.append(f"{np.mean(v) / 1000:.2f}" if v else "--")
        mem = [r["memory_bytes"] for r in rows if r["filter"] == f]
        cells.append(mem[0] if mem else "--")
        lines.append(f"| {LABEL[f]} | " + " | ".join(cells) + " |")
    lines += ["", f"TFG-CARMA3 magnetometer update: {cost['T41']['mag_update_madds']} multiply-adds.", ""]
    size = {r["variant"]: r for r in read(OUT / "esp32s3_size.csv")}
    link = {r["variant"]: r for r in read(OUT / "esp32s3_link.csv")}
    chain = {}
    for r in read(OUT / "esp32s3_callgraph.csv"):
        chain.setdefault(r["variant"], {})[r["root"]] = r
    if size:
        yes = {"True": "yes", "False": "no"}
        lines += ["### ESP32-S3 (Espressif GCC esp-14.2.0_20241119, -Os; build time only, no on-target timing)", "",
                  "Object sizes from `tests/kalman_carma/arduino/carma3_size_probe/build.sh`; link (with "
                  "`link_main.cpp` against the toolchain's newlib/libstdc++, default linker script) and static "
                  "call-graph chains from `callgraph.sh`.  Chain bytes sum the frames GCC reports along the deepest "
                  "direct-call path; library functions outside the translation unit count as 0.", "",
                  "| filter | code, object (B) | sizeof filter (B) | largest frame (B) | linked text (B) | "
                  "undefined refs | deepest chain from the sample entry (B) | deepest chain: time update / accel "
                  "update (B) | filter methods reach heap allocation |", "|---|---|---|---|---|---|---|---|---|"]
        for f in ("OU3", "A41", "T41"):
            s, lk, c = size.get(f), link.get(f, {}), chain.get(f, {})
            if not s:
                continue
            heap = any(c.get(k, {}).get("reaches_heap_allocation") == "True"
                       for k in ("time_update", "measurement_update_acc_only", "measurement_update_mag_only"))
            lines.append(f"| {LABEL[f]} | {s['code_bytes']} | {s['sizeof_filter']} | {s['max_stack_frame']} | "
                         f"{lk.get('text_bytes', '--')} | {lk.get('undefined_references', '--')} | "
                         f"{c.get('sample_path', {}).get('worst_stack_chain_bytes', '--')} | "
                         f"{c.get('time_update', {}).get('worst_stack_chain_bytes', '--')} / "
                         f"{c.get('measurement_update_acc_only', {}).get('worst_stack_chain_bytes', '--')} | "
                         f"{yes[str(heap)] if c else '--'} |")
        lines.append("")
    return lines


def table_math():
    p = OUT / "math_checks.txt"
    if not p.exists():
        return []
    lines = ["## 10. Mathematical verification (tests/kalman_carma/tfg_carma3-test)", "",
             "| check | measured | limit | result |", "|---|---|---|---|"]
    for line in p.read_text().splitlines():
        m = re.match(r"CHECK (.*?)\s+value=(\S+) limit=(\S+) (\w+)", line)
        if m:
            name = m.group(1).strip().replace("|", "\\|")   # a literal | would split the cell
            lines.append(f"| {name} | {m.group(2)} | {m.group(3)} | {m.group(4)} |")
        elif "checks passed" in line:
            lines += ["", line]
    lines.append("")
    return lines


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "svg.hashsalt": "tfg-carma3", "svg.fonttype": "none", "font.size": 9, "axes.titlesize": 10,
        "axes.labelsize": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK, "xtick.color": INK2,
        "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "grid.linestyle": "-",
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "lines.linewidth": 1.4,
        "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round", "figure.dpi": 100,
    })
    return plt


def save(fig, name):
    fig.savefig(OUT / name, format="svg", bbox_inches="tight", metadata={"Date": None})
    import matplotlib.pyplot as plt
    plt.close(fig)


def load_trace(sc, seed=1):
    p = OUT / "traces" / f"{sc.replace('+', '_')}_s{seed}.csv"
    if not p.exists():
        return None
    with p.open() as f:
        rd = csv.reader(f)
        hdr = next(rd)
        data = np.array([[float(x) if x not in ("", "nan", "-nan") else np.nan for x in row] for row in rd])
    return {h: data[:, i] for i, h in enumerate(hdr)}


TRACE_FLT = ["OU3", "A41", "T41"]


def thin(t, y, n=2500):
    """Min/max-preserving decimation for long time-series plots: each bucket
    of samples contributes its extremes, so the drawn envelope is the full-rate
    one with a fraction of the points (spectra use the full-rate data)."""
    if len(t) <= 2 * n:
        return t, y
    k = len(t) // n
    m = (len(t) // k) * k
    tb, yb = t[:m].reshape(-1, k), y[:m].reshape(-1, k)
    finite = np.isfinite(yb)
    lo = np.argmin(np.where(finite, yb, np.inf), axis=1)
    hi = np.argmax(np.where(finite, yb, -np.inf), axis=1)
    i1, i2 = np.minimum(lo, hi), np.maximum(lo, hi)
    r = np.arange(len(tb))
    return (np.column_stack([tb[r, i1], tb[r, i2]]).ravel(),
            np.column_stack([yb[r, i1], yb[r, i2]]).ravel())


def stride(t, y, n=2500):
    k = max(1, len(t) // n)
    return t[::k], y[::k]


def top_legend(fig, ax, title, ncol):
    """Title and a legend row above the first panel, clear of the data."""
    h, lab = ax.get_legend_handles_labels()
    fig.suptitle(title, y=1.0, fontsize=10)
    fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.5, 0.975), ncol=ncol, fontsize=8)
    fig.subplots_adjust(top=0.9)


def fig_timeseries(plt, tr, name, title, t0=600, t1=660):
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 6.0), sharex=True)
    m = (tr["t"] >= t0) & (tr["t"] <= t1)
    for ax, (k, lab) in zip(axes, (("z", "heave (m, NED down)"), ("x", "surge (m, north)"), ("y", "sway (m, east)"))):
        ax.plot(tr["t"][m], tr[f"{k}_true"][m], color=COL["truth"], lw=1.8, label="truth (wave coordinate)")
        for f in TRACE_FLT:
            ax.plot(tr["t"][m], tr[f"{k}_{f}"][m], color=COL[f], lw=1.2, label=LABEL[f])
        ax.set_ylabel(lab)
    axes[-1].set_xlabel("time (s)")
    top_legend(fig, axes[0], title, 4)
    save(fig, name)


def fig_errors(plt, tr, name, title):
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 6.0), sharex=True)
    for ax, (k, lab) in zip(axes, (("z", "heave error (m)"), ("x", "surge error (m)"), ("y", "sway error (m)"))):
        for f in TRACE_FLT:
            ax.plot(*thin(tr["t"], tr[f"e{k}_{f}"]), color=COL[f], lw=0.9, label=LABEL[f])
        ax.set_ylabel(lab)
    axes[-1].set_xlabel("time (s)")
    top_legend(fig, axes[0], title, 3)
    save(fig, name)


def fig_attitude(plt, tr, name, title):
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 6.0), sharex=True)
    for ax, (k, lab) in zip(axes, (("droll", "roll error (deg)"), ("dpitch", "pitch error (deg)"),
                                   ("dyaw", "heading error (deg)"))):
        for f in TRACE_FLT:
            ax.plot(*thin(tr["t"], tr[f"{k}_{f}"]), color=COL[f], lw=0.9, label=LABEL[f])
        ax.set_ylabel(lab)
    axes[-1].set_xlabel("time (s)")
    top_legend(fig, axes[0], title, 3)
    save(fig, name)


def welch(x, fs, nseg=1024):
    x = x[np.isfinite(x)]
    x = x - x.mean()
    w = np.hanning(nseg)
    step = nseg // 2
    segs = [x[i:i + nseg] * w for i in range(0, len(x) - nseg + 1, step)]
    P = np.mean([np.abs(np.fft.rfft(s)) ** 2 for s in segs], axis=0) * 2.0 / (fs * np.sum(w ** 2))
    return np.fft.rfftfreq(nseg, 1.0 / fs), P


def fig_psd(plt, tr, name, title):
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.2), sharey=True)
    m = tr["t"] >= 300
    fs = 1.0 / np.median(np.diff(tr["t"]))
    for ax, (k, lab) in zip(axes, (("z", "heave"), ("x", "surge"), ("y", "sway"))):
        for f in TRACE_FLT:
            fr, P = welch(tr[f"e{k}_{f}"][m], fs)
            ax.loglog(fr[1:], P[1:], color=COL[f], lw=1.2, label=LABEL[f])
        ax.set_title(f"{lab} error PSD")
        ax.set_xlabel("frequency (Hz)")
        ax.set_xlim(1e-3, 2.0)
    axes[0].set_ylabel("PSD (m^2/Hz)")
    axes[0].legend(loc="lower left", fontsize=7.5)
    fig.suptitle(title, y=1.02, fontsize=10)
    save(fig, name)


def fig_cov_vs_error(plt, rows, tr):
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8))
    ax = axes[0]
    m = tr["t"] >= 300
    for f in ("A41", "T41"):
        ax.plot(*thin(tr["t"][m], np.abs(tr[f"ez_{f}"][m])), color=COL[f], lw=0.5, alpha=0.5)
        ax.plot(*stride(tr["t"][m], 1.96 * tr[f"sz_{f}"][m]), color=COL[f], lw=1.8, label=f"{LABEL[f]}: 1.96 sigma")
    ax.set_title("heave |error| and predicted 95 % bound, JONSWAP seed 1")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("m")
    ax.legend(fontsize=7.5, loc="upper right", frameon=True, facecolor="white", edgecolor=GRID, framealpha=0.95)
    ax = axes[1]
    lim = 0
    for f in ("OU3", "A41", "T41"):
        for k, mk in (("x", "o"), ("z", "s")):
            xs, ys = [], []
            for sc in PRIMARY:
                pr, rm = by_seed(rows, sc, f, f"predstd_{k}"), by_seed(rows, sc, f, f"rms_{k}")
                s = sorted(set(pr) & set(rm))
                xs += [pr[i] for i in s]
                ys += [rm[i] for i in s]
            if xs:
                ax.scatter(xs, ys, s=10, color=COL[f], marker=mk, linewidths=0, alpha=0.7,
                           label=f"{LABEL[f]}, {'surge' if k == 'x' else 'heave'}")
                lim = max(lim, max(xs), max(ys))
    ax.plot([0, lim], [0, lim], color=INK2, lw=1)
    ax.set_xlabel("predicted RMS std (m)")
    ax.set_ylabel("empirical RMS error (m)")
    ax.set_title("predicted vs empirical, every principal-sea seed")
    ax.legend(fontsize=6.5, loc="lower right", ncol=1)
    save(fig, "covariance_vs_error.svg")


def fig_consistency(plt, rows):
    keys = [("nees_pos3", "position NEES (ideal 3)", 3), ("att_nees", "attitude NEES (ideal 3)", 3),
            ("yaw_nees", "yaw NEES (ideal 1)", 1), ("nis", "accelerometer NIS (ideal 3)", 3)]
    from matplotlib.ticker import FuncFormatter, NullFormatter
    flts = ["OU3", "A41", "T41", "T41R0"]
    regular = [sc for sc in PRIMARY if sc != "mixed_swell_wind"]
    fig, axes = plt.subplots(1, 4, figsize=(11.0, 3.4))
    for ax, (k, lab, ideal) in zip(axes, keys):
        data = [np.concatenate([list(by_seed(rows, sc, f, k).values()) for sc in regular]) for f in flts]
        bp = ax.boxplot(data, widths=0.5, patch_artist=True, showfliers=False)
        mixed = [float(np.median(list(by_seed(rows, "mixed_swell_wind", f, k).values()))) for f in flts]
        ax.scatter(range(1, len(flts) + 1), mixed, marker="D", s=26, color=INK, zorder=4)
        for patch, f in zip(bp["boxes"], flts):
            patch.set_facecolor(COL[f])
            patch.set_alpha(0.35)
            patch.set_edgecolor(COL[f])
        for el in ("whiskers", "caps", "medians"):
            for line in bp[el]:
                line.set_color(INK2)
        ax.axhline(ideal, color=INK, lw=1)
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_xticks(range(1, len(flts) + 1), [f for f in flts], rotation=0, fontsize=8)
        ax.set_title(lab, fontsize=9)
    axes[0].set_ylabel("per-run mean")
    fig.text(0.5, -0.06, "Boxes: JONSWAP, PM and changing f_p (300 runs per filter; whiskers 1.5 IQR); diamonds: "
             "mixed-sea median.  OU3 = OU-III, A41 = CARMA3 Q-MEKF, T41 = TFG-CARMA3, T41R0 = TFG, no reset "
             "transport", ha="center", fontsize=8, color=INK2)
    save(fig, "consistency_nees_nis.svg")


def fig_startup(plt, tr, name, title):
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.6), sharex=True)
    live = {f: tr["t"][np.isfinite(tr[f"ez_{f}"])][0] for f in TRACE_FLT}
    for ax, (keys, lab) in zip(axes, ((("z",), "heave"), (("x", "y"), "horizontal"))):
        for f in TRACE_FLT:
            e2 = sum(np.nan_to_num(tr[f"e{k}_{f}"]) ** 2 for k in keys)
            fs = 1.0 / np.median(np.diff(tr["t"]))
            L = int(10 * fs)
            roll = np.sqrt(np.convolve(e2, np.ones(L) / L, mode="same"))
            m = (tr["t"] >= live[f]) & (tr["t"] <= live[f] + 240)
            ax.plot(tr["t"][m] - live[f], roll[m], color=COL[f], lw=1.4, label=LABEL[f])
        ax.set_ylabel(f"{lab} 10 s RMS (m)")
    axes[0].set_title(title)
    axes[-1].set_xlabel("time after Live (s)")
    axes[0].legend(fontsize=8)
    save(fig, name)


def fig_distribution(plt, rows):
    flts = ["OU3", "A41", "T41", "T41R0"]
    fig, axes = plt.subplots(2, 4, figsize=(11.0, 5.2))
    for j, sc in enumerate(PRIMARY):
        for i, (k, lab) in enumerate((("rms_h", "horizontal RMS (m)"), ("rms_z", "heave RMS (m)"))):
            ax = axes[i, j]
            data = [list(by_seed(rows, sc, f, k).values()) for f in flts]
            for x, (f, d) in enumerate(zip(flts, data)):
                jitter = np.linspace(-0.15, 0.15, len(d)) if d else []
                ax.scatter(np.full(len(d), x + 1) + jitter, d, s=6, color=COL[f], linewidths=0, alpha=0.7)
                if d:
                    ax.plot([x + 0.7, x + 1.3], [np.median(d)] * 2, color=INK, lw=1.2)
            ax.set_xticks(range(1, len(flts) + 1), flts, fontsize=8)
            if i == 0:
                ax.set_title(sc_label(sc), fontsize=9)
            if j == 0:
                ax.set_ylabel(lab)
    fig.text(0.5, -0.02, "one dot per seed; bar = median.  OU3 = OU-III, A41 = CARMA3 Q-MEKF, T41 = TFG-CARMA3, "
             "T41R0 = TFG without reset transport", ha="center", fontsize=8, color=INK2)
    save(fig, "rms_distribution.svg")


def fig_conditions(plt, stats):
    m = stats.get("matrix", {})
    ma = stats.get("matrix_att", {})
    scs = list(m)
    keys = [("rms_h", "horizontal RMS"), ("rms_z", "heave RMS"), ("yaw_rms_deg", "yaw RMS"), ("tilt_rms_deg", "tilt RMS")]
    fig, axes = plt.subplots(1, 4, figsize=(13.6, 0.28 * len(scs) + 1.4), sharey=True)
    y = np.arange(len(scs))
    for ax, (k, lab) in zip(axes, keys):
        for i, sc in enumerate(scs):
            ps = (ma if k == "tilt_rms_deg" else m)[sc][k]
            if not ps:
                continue
            lo, hi = ps["ci_pct"]
            ax.plot([lo, hi], [i, i], color=COL["T41"], lw=2)
            ax.scatter([ps["change_pct"]], [i], s=22, color=COL["T41"], edgecolors="white", linewidths=1, zorder=3)
        ax.axvline(0, color=INK, lw=1)
        ax.axvspan(-EQUIV_MARGIN, EQUIV_MARGIN, color=GRID, alpha=0.6, lw=0)
        ax.set_title(f"{lab}: TFG vs Q-MEKF (%)", fontsize=9)
    axes[0].set_yticks(y, [sc_label(s) for s in scs], fontsize=7.5)
    axes[0].invert_yaxis()
    fig.supxlabel("1 - mean(TFG)/mean(Q-MEKF), % (> 0: TFG lower).  Dot = paired change, bar = 95 % bootstrap "
                  "interval, shaded = +-1 % equivalence margin", fontsize=8, color=INK2)
    save(fig, "paired_change_by_condition.svg")


def fig_error_level(plt, rows, truth):
    """Paired attitude change per condition against how large CARMA3's own
    state error is (left; the velocity error, the recorded state closest to
    the acceleration state) and against the true horizontal wave
    acceleration (right): first-order equivalence predicts differences that
    grow with the state corrections, not with the motion itself (section 8)."""
    acc = {}
    for r in truth:
        a = [float(x) for x in r["accel_wave_rms"].split(",")]
        acc.setdefault(r["scenario"], []).append(math.hypot(a[0], a[1]))
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.4), sharey=True)
    marks = {"mixed_swell_wind": ("mixed", (-34, -3)), "mixed_swell_wind+T3600": ("mixed, 3600 s", (-72, -3)),
             f"{J29}+maneuver": ("manoeuvring", (-28, -15)), J29: ("JONSWAP", (-40, 26)),
             f"{J29}+T3600+maggap2400": ("2400 s without magnetometer", (-170, -34)),
             "stationary": ("still IMU", (-20, 30))}
    for f, lab in (("T41", LABEL["T41"]), ("T41R0", LABEL["T41R0"])):
        xs_err, xs_acc, ys, names = [], [], [], []
        for sc in all_scenarios(rows):
            ps = paired(rows, sc, "A41", f, "tilt_rms_deg")
            err = math.hypot(mean_of(rows, sc, "A41", "rmsv_h"), mean_of(rows, sc, "A41", "rmsv_z"))
            if not ps or not math.isfinite(err):
                continue
            xs_err.append(err)
            xs_acc.append(float(np.mean(acc[sc])) if sc in acc else math.nan)
            ys.append(ps["change_pct"])
            names.append(sc)
        for ax, xs in ((axes[0], xs_err), (axes[1], xs_acc)):
            ax.scatter(xs, ys, s=30, color=COL[f], edgecolors="white", linewidths=1.2, zorder=3, label=lab)
        if f == "T41":
            for x, y, sc in zip(xs_err, ys, names):
                if sc in marks:
                    axes[0].annotate(marks[sc][0], (x, y), textcoords="offset points", xytext=marks[sc][1],
                                     fontsize=7.5, color=INK2,
                                     arrowprops={"arrowstyle": "-", "color": INK2, "lw": 0.6, "shrinkB": 4})
    for ax in axes:
        ax.axhline(0, color=INK, lw=1)
        ax.axhspan(-EQUIV_MARGIN, EQUIV_MARGIN, color=GRID, alpha=0.6, lw=0)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("CARMA3 Q-MEKF 3-D velocity RMS error of the condition (m/s, log)")
    axes[1].set_xlabel("true horizontal wave acceleration RMS (m/s^2)")
    axes[0].set_ylabel("tilt RMS change vs Q-MEKF (%, > 0: TFG lower)")
    axes[0].legend(fontsize=8, loc="upper left")
    fig.suptitle("Where the two error geometries separate: one point per condition (paired seeds)", fontsize=10)
    save(fig, "mechanism_error_level.svg")


def fig_cost(plt, timing):
    if not timing:
        return
    flts = ["OU3", "A41", "T41", "T41R0", "T41RA"]
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.0))
    cost = cost_models()
    madds = [cost[f]["total_madds_per_sample"] for f in flts]
    t = [np.mean([fnum(r, "runtime_ns") for r in timing if r["filter"] == f and r["clock"] == "fixed"]) / 1000 for f in flts]
    mem = [np.mean([fnum(r, "memory_bytes") for r in timing if r["filter"] == f]) for f in flts]
    for ax, vals, lab, fmt in ((axes[0], madds, "multiply-adds per sample (analytic)", "{:.0f}"),
                               (axes[1], t, "host update time (us), fixed dt", "{:.2f}"),
                               (axes[2], mem, "wrapper + filter object (bytes, host)", "{:.0f}")):
        ax.barh(range(len(flts)), vals, height=0.5, color=[COL[f] for f in flts])
        for i, v in enumerate(vals):
            ax.text(v, i, " " + fmt.format(v), va="center", fontsize=8, color=INK)
        ax.set_yticks(range(len(flts)), [LABEL[f] for f in flts], fontsize=8)
        ax.invert_yaxis()
        ax.set_title(lab, fontsize=9)
        ax.set_xlim(0, max(vals) * 1.3)
    for ax in axes[1:]:
        ax.set_yticklabels([])
    save(fig, "cost.svg")


def fig_maggap(plt, tr, name, title, gap0, gap1):
    flts = ["A41", "T41", "T41R0"]
    fig, axes = plt.subplots(len(flts), 1, figsize=(7.2, 6.2), sharex=True, sharey=True)
    for ax, f in zip(axes, flts):
        ts, s = stride(tr["t"], tr[f"syaw_{f}"])
        ax.fill_between(ts, -2 * s, 2 * s, color=COL[f], alpha=0.3, lw=0, label="predicted +-2 sigma")
        ax.plot(*thin(tr["t"], tr[f"dyaw_{f}"]), color=COL[f], lw=1.0, label="heading error")
        ax.axvspan(gap0, gap1, color=GRID, alpha=0.5, lw=0)
        ax.set_ylabel("deg")
        ax.set_title(LABEL[f], fontsize=9, loc="left")
        ax.legend(fontsize=7.5, loc="upper left", ncol=2)
    live = np.isfinite(tr[f"dyaw_{flts[0]}"]) & (tr["t"] > 200)
    span = max(3.0, 1.2 * max(np.nanmax(np.abs(tr[f"dyaw_{f}"][live])) for f in flts))
    axes[0].set_ylim(-span, span)
    axes[-1].set_xlabel("time (s); shaded = no magnetometer")
    fig.suptitle(title, y=1.0, fontsize=10)
    save(fig, name)


def figures(rows, stats, timing, truth):
    plt = mpl()
    tr = load_trace(J29)
    if tr is not None:
        fig_timeseries(plt, tr, "timeseries_xyz_jonswap.svg", "JONSWAP 0.29 Hz, 0.8 m, seed 1: truth and estimates")
        fig_errors(plt, tr, "errors_xyz_jonswap.svg", "JONSWAP 0.29 Hz, 0.8 m, seed 1: translational error")
        fig_attitude(plt, tr, "attitude_errors_jonswap.svg", "JONSWAP 0.29 Hz, 0.8 m, seed 1: attitude error")
        fig_psd(plt, tr, "error_psd_jonswap.svg", "JONSWAP 0.29 Hz, 0.8 m, seed 1, t >= 300 s")
        fig_startup(plt, tr, "startup_jonswap.svg", "startup in waves, JONSWAP 0.29 Hz, 0.8 m, seed 1")
        fig_cov_vs_error(plt, rows, tr)
    for sc, stem, ttl in (("stationary", "stationary", "still IMU, seed 1"),
                          ("mixed_swell_wind", "mixed", "mixed swell + wind, seed 1")):
        t2 = load_trace(sc)
        if t2 is not None:
            fig_errors(plt, t2, f"errors_xyz_{stem}.svg", f"{ttl}: translational error")
            fig_psd(plt, t2, f"error_psd_{stem}.svg", f"{ttl}, t >= 300 s")
            fig_attitude(plt, t2, f"attitude_errors_{stem}.svg", f"{ttl}: attitude error")
    for sc, gap0, gap1, stem in ((f"{J29}+maggap", 400.0, 700.0, "maggap300"),
                                 (f"{J29}+T1800+maggap900", 600.0, 1500.0, "maggap900"),
                                 (f"{J29}+T3600+maggap2400", 1200.0, 3600.0, "maggap2400")):
        t3 = load_trace(sc)
        if t3 is not None:
            fig_maggap(plt, t3, f"heading_{stem}.svg", f"{sc_label(sc)}, seed 1", gap0, gap1)
    t4 = load_trace(f"{J29}+turns")
    if t4 is not None:
        fig_attitude(plt, t4, "attitude_errors_turns.svg", "90 deg turns in 6 s every 120 s, seed 1: attitude error")
    fig_consistency(plt, rows)
    fig_distribution(plt, rows)
    fig_conditions(plt, stats)
    fig_error_level(plt, rows, truth)
    fig_cost(plt, timing)
    t5 = load_trace("mixed_swell_wind+maggap")
    if t5 is not None:
        fig_maggap(plt, t5, "heading_mixed_maggap300.svg", "mixed swell + wind, 300 s magnetometer gap, seed 1",
                   400.0, 700.0)


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, list):
        return [clean(v) for v in o]
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


def main():
    rows = read(OUT / "raw_runs.csv")
    if not rows:
        raise SystemExit("no raw runs: python3 tools/tfg_carma3_benchmark.py run")
    follow = read(OUT / "followup_runs.csv")
    dataset = read(OUT / "dataset_runs.csv")
    timing = read(OUT / "timing_runs.csv")
    mech = read(OUT / "mechanism_runs.csv")
    truth = read(OUT / "truth_levels.csv")
    stats = {}
    failed = [(r["scenario"], r["seed"], r["filter"]) for r in rows + follow if r.get("failed", "0") != "0"]
    lines = ["# TFG-CARMA3 benchmark summary (generated)", "",
             "Generated by `python3 tools/tfg_carma3_benchmark.py report` from the CSV files in this directory; "
             "see docs/tfg-carma3-research.md for the design and interpretation.  Simulation only.", "",
             f"Runs: {len({(r['scenario'], r['seed']) for r in rows})} paired synthetic runs "
             f"({len(rows)} filter runs), {len({(r['scenario'], r['seed']) for r in follow})} follow-up runs; "
             f"failed filter runs: {len(failed)} {failed[:10] if failed else ''}.", ""]
    lines += table_primary(rows, stats)
    lines += table_matrix(rows, stats)
    lines += table_matrix_att(rows, stats)
    lines += table_ablations(rows, stats)
    lines += table_values(rows, "## 4. Consistency (mean over seeds; ideal value in brackets)", CONSIST,
                          PRIMARY + GAPS + ["stationary", f"{J29}+turns"],
                          ["OU3", "A41", "T41", "T41R0", "T41RA", "T41NW"])
    lines += table_maggap(rows, stats)
    lines += table_values(rows, "## 6. Startup (first 120 s after Live) and rapid rotation", STARTUP +
                          [("turn_peak_deg", "peak attitude error after a turn (deg)"),
                           ("turn_recover_s", "attitude recovery after a turn (s)")],
                          PRIMARY + [f"{J29}+turns", "large_rotation"], ["OU3", "A41", "T41", "T41R0"])
    for sc in PRIMARY + [f"{J29}+turns"]:
        for k, _ in STARTUP + [("turn_peak_deg", ""), ("turn_recover_s", "")]:
            stats.setdefault("startup", {}).setdefault(sc, {})[k] = paired(rows, sc, "A41", "T41", k)
    lines += table_followup(follow, stats)
    lines += table_dataset(dataset)
    lines += table_timing(timing, stats)
    lines += table_math()
    lines += table_mechanism(mech, stats)
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")
    (OUT / "paired_stats.json").write_text(json.dumps(clean(stats), indent=1, sort_keys=True) + "\n")
    figures(rows, stats, timing, truth)
    print(f"report: {OUT.relative_to(ROOT)}/summary.md, paired_stats.json and figures")


if __name__ == "__main__":
    main()
