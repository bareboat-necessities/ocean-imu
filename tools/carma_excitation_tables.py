"""Validation tables of the CARMA wave-excitation study (tools/carma_excitation_report.py)."""
from __future__ import annotations

import re
from collections import defaultdict

import numpy as np

from carma_excitation_report import LABEL, aligned, fmt_imp, fnum, paired, rms_h

FULL = ["OU3", "A41@A", "A41@B0", "A41@B", "A41@C", "A41@Cnf", "A41@Cng", "T41@A", "T41@B", "T41@C"]
CORE = ["OU3", "A41@A", "A41@B0", "A41@B", "A41@C", "T41@A", "T41@C"]


def hs_of(sc):
    m = re.search(r"_hs([0-9.]+)", sc)
    return float(m.group(1)) if m else 0.0


def fp_of(sc):
    m = re.search(r"_fp([0-9.]+)", sc)
    return float(m.group(1)) if m else np.nan


def peak_h(r):
    return max(fnum(r, "peak_x"), fnum(r, "peak_y"))


def p95_h(r):
    return max(fnum(r, "p95_x"), fnum(r, "p95_y"))


def table(r, rows, header):
    r.w("| " + " | ".join(header) + " |")
    r.w("|" + "---|" * len(header))
    for row in rows:
        r.w("| " + " | ".join(row) + " |")
    r.w()


def mean_mm(v):
    v = np.asarray(v, float)
    return f"{1e3 * np.nanmean(v):.2f}" if np.isfinite(v).any() else "-"


def e1(r):
    r.w("## Experiment 1: still water (H_s = 0), 2400 s records")
    r.w()
    conds = sorted(r.cond("E1"), key=lambda k: k[1])
    if not conds:
        return
    metrics = [("heave RMS", lambda x: fnum(x, "rms_z")), ("horizontal RMS", rms_h),
               ("heave peak", lambda x: fnum(x, "peak_z")), ("horizontal peak", peak_h),
               ("heave p95", lambda x: fnum(x, "p95_z")), ("startup heave RMS", lambda x: fnum(x, "su_rms_z"))]
    for (e, sc) in conds:
        d = r.g[(e, sc)]
        fl = [f for f in FULL if f in d]
        rows = []
        for name, key in metrics:
            seeds, v = aligned(d, fl, key)
            rows.append([name + " (mm)"] + [mean_mm(v[f]) for f in fl])
            for f in fl:
                if f not in ("A41@A", "T41@A"):
                    base = "T41@A" if f.startswith("T41@A") else "A41@A"
                    ps = paired(v[f], v[base])
                    if ps:
                        r.stats.setdefault("E1", {}).setdefault(sc, {}).setdefault(name, {})[f] = ps
        seeds, v = aligned(d, fl, lambda x: fnum(x, "exc_p_frac50"))
        rows.append(["false wave detection (time with p >= 0.5)"] +
                    [("-" if f == "OU3" else f"{100 * np.nanmean(v[f]):.2f} %") for f in fl])
        seeds, v = aligned(d, fl, lambda x: fnum(x, "exc_sig_med"))
        rows.append(["CARMA sigma, median (m/s^2)"] + [("-" if f == "OU3" else f"{np.nanmean(v[f]):.2g}") for f in fl])
        seeds, v = aligned(d, fl, lambda x: fnum(x, "exc_fz_med"))
        rows.append(["CARMA frequency, median (Hz)"] + [("-" if f == "OU3" else f"{np.nanmean(v[f]):.3f}") for f in fl])
        seeds, v = aligned(d, fl, lambda x: fnum(x, "nees_pv6"))
        rows.append(["NEES p,v (6 dof)"] + [f"{np.nanmean(v[f]):.3g}" for f in fl])
        seeds, v = aligned(d, fl, lambda x: fnum(x, "tilt_rms_deg"))
        rows.append(["tilt RMS (deg)"] + [f"{np.nanmean(v[f]):.3f}" for f in fl])
        seeds, v = aligned(d, fl, lambda x: fnum(x, "ba_rms_mg"))
        rows.append(["accel bias error RMS (mg)"] + [f"{np.nanmean(v[f]):.2f}" for f in fl])
        r.w(f"### {sc} ({len(seeds)} seeds)")
        r.w()
        table(r, rows, ["metric"] + [LABEL.get(f, f) for f in fl])
    r.improvement_table("E1", ["A41@B0", "A41@B", "A41@C", "A41@Cnf", "A41@Cng"], "A41@A", "heave RMS",
                        lambda x: fnum(x, "rms_z"), "Paired change of the still-water heave RMS against the legacy CARMA3")
    r.improvement_table("E1", ["A41@B0", "A41@B", "A41@C"], "A41@A", "horizontal RMS", rms_h,
                        "Paired change of the still-water horizontal RMS against the legacy CARMA3")
    r.improvement_table("E1", ["T41@B", "T41@C"], "T41@A", "heave RMS", lambda x: fnum(x, "rms_z"),
                        "TFG-CARMA3: paired change of the still-water heave RMS against the legacy TFG-CARMA3")
    r.improvement_table("E1", ["T41@B", "T41@C"], "T41@A", "horizontal RMS", rms_h,
                        "TFG-CARMA3: paired change of the still-water horizontal RMS against the legacy TFG-CARMA3")


def small_waves(r, exp, title, filters, note=None):
    r.w(f"## {title}")
    r.w()
    if note:
        r.w(note)
        r.w()
    by = defaultdict(list)
    for (e, sc) in r.cond(exp):
        key = re.sub(r"_hs[0-9.]+", "", sc)
        by[key].append(sc)
    for fam in sorted(by):
        scs = sorted(by[fam], key=hs_of)
        rows = []
        for sc in scs:
            d = r.g[(exp, sc)]
            fl = [f for f in filters if f in d]
            seeds, vz = aligned(d, fl, lambda x: fnum(x, "rms_z"))
            _, vt = aligned(d, fl, lambda x: fnum(x, "std_z"))
            _, vh = aligned(d, fl, rms_h)
            _, vp = aligned(d, fl, lambda x: fnum(x, "exc_p_frac50"))
            _, va = aligned(d, fl, lambda x: fnum(x, "exc_pz_rms"))
            _, vaz = aligned(d, fl, lambda x: fnum(x, "true_az_rms"))
            truth = np.nanmean(vt[fl[0]])
            hs = hs_of(sc)
            cells = [f"{hs:g}", str(len(seeds)), f"{1e3 * truth:.2f}",
                     f"{1e3 * np.nanmean(vaz[fl[0]]):.1f}" if np.isfinite(vaz[fl[0]]).any() else "-"]
            for f in fl:
                cells.append(f"{1e3 * np.mean(vz[f]):.2f} ({100 * np.mean(vz[f]) / hs:.0f} %)")
            for f in ("A41@C", "T41@C"):
                if f in fl:
                    cells.append(f"{np.nanmean(va[f]) / truth:.2f} / {100 * (1 - np.nanmean(vp[f])):.1f} %")
            rows.append(cells)
            for f in fl:
                if f in ("A41@A", "T41@A", "OU3"):
                    continue
                base = "T41@A" if f.startswith("T41@A") else "A41@A"
                for name, v in (("heave RMS", vz), ("horizontal RMS", vh)):
                    ps = paired(v[f], v[base])
                    if ps:
                        r.stats.setdefault(exp, {}).setdefault(sc, {}).setdefault(name, {})[f] = ps
        fl = [f for f in filters if f in r.g[(exp, scs[0])]]
        hdr = ["H_s (m)", "n", "true heave RMS (mm)", "true vert. accel RMS (mm/s^2)"] + \
              [f"{LABEL.get(f, f)} heave RMS err (mm, % H_s)" for f in fl] + \
              [f"{LABEL.get(f, f)}: estimate/true RMS / missed-wave time" for f in ("A41@C", "T41@C") if f in fl]
        r.w(f"### {fam.replace('_g3.3', '').replace('_g5', ' (gamma 5)')}")
        r.w()
        table(r, rows, hdr)
    r.improvement_table(exp, ["A41@B0", "A41@B", "A41@C"], "A41@A", "heave RMS", lambda x: fnum(x, "rms_z"),
                        f"{title}: paired change of the heave RMS against the legacy CARMA3")
    r.improvement_table(exp, ["A41@B0", "A41@B", "A41@C"], "A41@A", "horizontal RMS", rms_h,
                        f"{title}: paired change of the horizontal RMS against the legacy CARMA3")
    r.improvement_table(exp, ["T41@C"], "T41@A", "heave RMS", lambda x: fnum(x, "rms_z"),
                        f"{title}: TFG-CARMA3, paired change of the heave RMS against the legacy TFG-CARMA3")


def e4(r):
    r.w("## Experiment 4: wave onset and disappearance")
    r.w()
    if not r.trans:
        return
    r.w("Per transition (start t_s = t_c - w/2): detection delay (first p >= 0.5 for an onset, first p < 0.5 after "
        "waves stop), time for the CARMA sigma to cover 90 % of its change, heave RMS error in [t_s, t_s + 60 s] and "
        "[t_s + 60, t_s + 180 s], horizontal RMS in the first 120 s.  Medians for times, means for errors (mm), "
        "over seeds and over the transitions of the same kind in the record.")
    r.w()
    groups = defaultdict(lambda: defaultdict(list))
    for t in r.trans:
        groups[(t["scenario"], t["kind"], t["level_from"], t["level_to"])][t["filter"]].append(t)
    rows = []
    for (sc, kind, lf, lt) in sorted(groups):
        d = groups[(sc, kind, lf, lt)]
        for f in ["A41@A", "A41@B", "A41@C", "A41@Cnf", "T41@A", "T41@C"]:
            if f not in d:
                continue
            x = d[f]
            def med(k):
                v = np.array([fnum(t, k) for t in x])
                v = v[np.isfinite(v)] if k != "detect_s" else v
                return np.median(v) if v.size else np.nan
            def mean(k):
                v = np.array([fnum(t, k) for t in x])
                return np.nanmean(v)
            det = med("detect_s") if f != "A41@A" and f != "T41@A" else med("detect_s")
            rows.append([sc.replace("jonswap_fp0.29_", "J29 ").replace("_g3.3", ""), f"{kind} {float(lf):g}->{float(lt):g}",
                         LABEL.get(f, f), str(len(x)),
                         f"{det:.1f}" if np.isfinite(det) else "never",
                         f"{med('sigma_settle_s'):.1f}" if np.isfinite(med('sigma_settle_s')) else "-",
                         f"{1e3 * mean('rms_z_0_60'):.1f}", f"{1e3 * mean('rms_z_60_180'):.1f}",
                         f"{1e3 * mean('true_z_0_60'):.1f}", f"{1e3 * mean('peak_z_0_120'):.0f}",
                         f"{1e3 * mean('rms_x_0_120'):.1f}"])
        # Paired C vs A on the early window.
        for c, b in (("A41@C", "A41@A"), ("A41@B", "A41@A"), ("T41@C", "T41@A")):
            if c in d and b in d:
                ka = {(t["seed"], t["step"]): fnum(t, "rms_z_0_60") for t in d[c]}
                kb = {(t["seed"], t["step"]): fnum(t, "rms_z_0_60") for t in d[b]}
                keys = sorted(set(ka) & set(kb))
                ps = paired([ka[k] for k in keys], [kb[k] for k in keys])
                if ps:
                    r.stats.setdefault("E4", {}).setdefault(f"{sc} {kind} {lf}->{lt}", {})[f"{c}_vs_{b}_rms_z_0_60"] = ps
    table(r, rows, ["record", "transition", "filter", "n", "detection (s)", "sigma settles (s)", "heave 0-60 s",
                    "heave 60-180 s", "true heave 0-60 s", "heave peak 0-120 s", "horizontal 0-120 s"])
    r.w("Paired change of the heave RMS error in the first 60 s after each transition:")
    r.w()
    rows = []
    for k, v in sorted(r.stats.get("E4", {}).items()):
        for name, ps in v.items():
            rows.append([k, name.replace("_rms_z_0_60", ""), str(ps["n"]), fmt_imp(ps)])
    table(r, rows, ["transition", "pair", "n", "change"])


def e5(r):
    r.w("## Experiment 5: sensor-noise uncertainty and disturbances")
    r.w()
    r.w("`/nv0.5`: the evidence model's noise variance is half the true one (underestimated by 50 %); `/nv2`: twice "
        "(overestimated by 100 %).  cnoise: 4 mg first-order Gauss-Markov noise (2 s), not in the model; nvar: white "
        "noise variance cycling 1x-2x; engine<rpm>: inboard diesel; temp: IMU temperature 25-45 C with a "
        "mis-modelled bias temperature coefficient.")
    r.w()
    fl = ["A41@A", "A41@B", "A41@B/nv0.5", "A41@B/nv2", "A41@C", "A41@C/nv0.5", "A41@C/nv2", "T41@C/nv0.5", "T41@C/nv2"]
    rows = []
    for (e, sc) in sorted(r.cond("E5"), key=lambda k: k[1]):
        d = r.g[(e, sc)]
        f2 = [f for f in fl if f in d]
        seeds, vz = aligned(d, f2, lambda x: fnum(x, "rms_z"))
        _, vh = aligned(d, f2, rms_h)
        _, vp = aligned(d, f2, lambda x: fnum(x, "exc_p_frac50"))
        cells = [sc.replace("jonswap_fp0.29_", "J29 ").replace("_g3.3", ""), str(len(seeds))]
        for f in fl:
            if f in d:
                cells.append(f"{1e3 * np.mean(vz[f]):.1f} / {1e3 * np.mean(vh[f]):.1f} / {100 * np.nanmean(vp[f]):.0f} %")
            else:
                cells.append("-")
        rows.append(cells)
        for f in f2:
            if f == "A41@A":
                continue
            ps = paired(vz[f], vz["A41@A"])
            if ps:
                r.stats.setdefault("E5", {}).setdefault(sc, {}).setdefault("heave RMS", {})[f] = ps
    table(r, rows, ["condition", "n"] + [f"{f} heave / horiz. (mm) / p >= 0.5" for f in fl])


def e6(r):
    r.w("## Experiment 6: non-wave vessel motion")
    r.w()
    fl = ["OU3", "A41@A", "A41@B", "A41@C", "T41@A", "T41@C"]
    rows = []
    for (e, sc) in sorted(r.cond("E6"), key=lambda k: k[1]):
        d = r.g[(e, sc)]
        f2 = [f for f in fl if f in d]
        seeds, vz = aligned(d, f2, lambda x: fnum(x, "rms_z"))
        _, vh = aligned(d, f2, rms_h)
        _, vt = aligned(d, f2, lambda x: fnum(x, "tilt_rms_deg"))
        _, vy = aligned(d, f2, lambda x: fnum(x, "yaw_w_rms_deg"))
        _, vb = aligned(d, f2, lambda x: fnum(x, "ba_rms_mg"))
        for f in f2:
            rows.append([sc.replace("jonswap_fp0.29_", "J29 ").replace("_g3.3", ""), LABEL.get(f, f), str(len(seeds)),
                         f"{1e3 * np.mean(vz[f]):.2f}", f"{1e3 * np.mean(vh[f]):.1f}", f"{np.nanmean(vt[f]):.3f}",
                         f"{np.nanmean(vy[f]):.2f}", f"{np.nanmean(vb[f]):.2f}"])
            if f not in ("A41@A", "T41@A", "OU3"):
                base = "T41@A" if f.startswith("T41@A") else "A41@A"
                for name, v in (("heave RMS", vz), ("horizontal RMS", vh), ("tilt RMS", vt), ("yaw RMS", vy)):
                    ps = paired(v[f], v[base])
                    if ps:
                        r.stats.setdefault("E6", {}).setdefault(sc, {}).setdefault(name, {})[f] = ps
    table(r, rows, ["condition", "filter", "n", "heave RMS (mm)", "horizontal RMS (mm)", "tilt RMS (deg)",
                    "yaw RMS (deg)", "accel bias error (mg)"])
    r.improvement_table("E6", ["A41@B", "A41@C"], "A41@A", "horizontal RMS", rms_h,
                        "Non-wave motion: paired change of the horizontal RMS against the legacy CARMA3")
    r.improvement_table("E6", ["A41@B", "A41@C"], "A41@A", "tilt RMS", lambda x: fnum(x, "tilt_rms_deg"),
                        "Non-wave motion: paired change of the tilt RMS against the legacy CARMA3")


def e7(r):
    r.w("## Experiment 7: realistic and mixed seas")
    r.w()
    fl = CORE
    rows = []
    for (e, sc) in sorted(r.cond("E7"), key=lambda k: k[1]):
        d = r.g[(e, sc)]
        f2 = [f for f in fl if f in d]
        seeds, vz = aligned(d, f2, lambda x: fnum(x, "rms_z"))
        _, vh = aligned(d, f2, rms_h)
        _, vn = aligned(d, f2, lambda x: fnum(x, "nees_pv6"))
        _, vs = aligned(d, f2, lambda x: fnum(x, "nis"))
        hs = fnum(next(iter(d[f2[0]].values())), "hs_incident")
        for f in f2:
            rows.append([sc, LABEL.get(f, f), str(len(seeds)), f"{100 * np.mean(vz[f]) / hs:.2f}",
                         f"{100 * np.mean(vh[f]) / hs:.2f}", f"{np.nanmean(vn[f]):.2f}", f"{np.nanmean(vs[f]):.3f}"])
    table(r, rows, ["sea", "filter", "n", "heave RMS (% H_s)", "horizontal RMS (% H_s)", "NEES p,v (6 dof)",
                    "accel NIS (3 dof)"])
    r.improvement_table("E7", ["A41@B", "A41@C"], "A41@A", "heave RMS", lambda x: fnum(x, "rms_z"),
                        "Mixed and realistic seas: paired change of the heave RMS against the legacy CARMA3")
    r.improvement_table("E7", ["A41@B", "A41@C"], "A41@A", "horizontal RMS", rms_h,
                        "Mixed and realistic seas: paired change of the horizontal RMS against the legacy CARMA3")
    r.improvement_table("E7", ["T41@C"], "T41@A", "heave RMS", lambda x: fnum(x, "rms_z"),
                        "TFG-CARMA3: paired change of the heave RMS against the legacy TFG-CARMA3")
    r.improvement_table("E7", ["T41@C"], "T41@A", "horizontal RMS", rms_h,
                        "TFG-CARMA3: paired change of the horizontal RMS against the legacy TFG-CARMA3")


def e8(r):
    r.w("## Experiment 8: CARMA(3,0) and CARMA(5,2)")
    r.w()
    fl = ["A30@A", "A30@C", "A52@A", "A52@C", "T52@A", "T52@C"]
    rows = []
    for (e, sc) in sorted(r.cond("E8"), key=lambda k: k[1]):
        d = r.g[(e, sc)]
        f2 = [f for f in fl if f in d]
        seeds, vz = aligned(d, f2, lambda x: fnum(x, "rms_z"))
        _, vh = aligned(d, f2, rms_h)
        rows.append([sc, str(len(seeds))] + [f"{1e3 * np.mean(vz[f]):.1f} / {1e3 * np.mean(vh[f]):.1f}" for f in f2])
        for c, b in (("A30@C", "A30@A"), ("A52@C", "A52@A"), ("T52@C", "T52@A")):
            if c in vz and b in vz:
                for name, v in (("heave RMS", vz), ("horizontal RMS", vh)):
                    ps = paired(v[c], v[b])
                    if ps:
                        r.stats.setdefault("E8", {}).setdefault(sc, {}).setdefault(name, {})[c] = ps
    table(r, rows, ["condition", "n"] + [f"{f} heave / horiz. (mm)" for f in fl])


def cost(r):
    r.w("## Host cost")
    r.w()
    acc = defaultdict(lambda: [[], []])
    for (e, sc), d in r.g.items():
        if e not in ("E7", "E2"):
            continue
        for f, s in d.items():
            for x in s.values():
                acc[f][0].append(fnum(x, "runtime_ns"))
                acc[f][1].append(fnum(x, "memory_bytes"))
    rows = [[LABEL.get(f, f), f"{np.nanmedian(v[0]) / 1e3:.2f}", f"{np.nanmedian(v[1]):.0f}"]
            for f, v in sorted(acc.items()) if f in CORE + ["A41@Cnf", "A41@Cng", "T41@B"]]
    r.w("Median per-sample update time on the benchmark host (x86-64, GCC -O3 -march=native, several filters per "
        "process, four processes in parallel: a relative figure only) and object size.")
    r.w()
    table(r, rows, ["filter", "update (us)", "wrapper + filter (bytes)"])


def validation(r):
    e1(r)
    small_waves(r, "E2", "Experiment 2: very small waves", CORE + ["A41@Cnf", "A41@Cng", "T41@B"])
    small_waves(r, "E3", "Experiment 3: long-period, low-acceleration swell", CORE,
                note="JONSWAP gamma 5; peak periods 8-25 s.")
    e4(r)
    e5(r)
    e6(r)
    e7(r)
    e8(r)
    cost(r)
