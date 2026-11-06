#!/usr/bin/env python3
"""Tables, paired statistics and figures of the CARMA wave-excitation study.

Reads reports/results/carma_excitation/{raw_runs,tuning_runs,transitions,psd_calm}.csv
and traces/, writes summary.md, paired_stats.json and figures/*.svg.  Paired
statistics follow tools/carma3_report.py: improvement = 1 - mean RMS(candidate) /
mean RMS(baseline) over the same seeds with a 95 % bootstrap interval (10 000
resamples, fixed seed), win rate and Wilcoxon signed-rank p.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

LABEL = {
    "OU3": "OU-III (shipping)", "A41@A": "A: CARMA3 legacy", "A41@B0": "B0: floor removed",
    "A41@B": "B: debiased", "A41@C": "C: evidence", "A41@Cnf": "C without frequency hold",
    "A41@Cng": "C without evidence weight", "T41@A": "A: TFG legacy", "T41@B": "B: TFG debiased",
    "T41@C": "C: TFG evidence",
}
COLOR = {"OU3": "#7f7f7f", "A41@A": "#d62728", "A41@B0": "#ff9896", "A41@B": "#ff7f0e", "A41@C": "#1f77b4",
         "T41@A": "#9467bd", "T41@C": "#17becf", "A41@Cnf": "#2ca02c", "A41@Cng": "#bcbd22", "T41@B": "#e377c2"}


def fnum(r, k):
    try:
        return float(r[k])
    except (KeyError, ValueError, TypeError):
        return np.nan


def load(path):
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def bootstrap_ci(x, stat=np.mean, n=10000, seed=20261025):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), size=(n, len(x)))
    s = stat(x[idx], axis=1)
    return (float(np.percentile(s, 2.5)), float(np.percentile(s, 97.5)))


def paired(a, b):
    """Candidate a against baseline b on the same seeds (lower is better)."""
    from scipy.stats import wilcoxon
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if len(a) < 2 or not b.mean() > 0:
        return None
    d = a - b
    rng = np.random.default_rng(20261025)
    idx = rng.integers(0, len(a), size=(10000, len(a)))
    rel_bs = 1.0 - a[idx].mean(axis=1) / b[idx].mean(axis=1)
    try:
        p = float(wilcoxon(a, b).pvalue) if len(a) >= 6 and np.any(d != 0) else np.nan
    except ValueError:
        p = np.nan
    return {"n": int(len(a)), "mean_a": float(a.mean()), "mean_b": float(b.mean()),
            "improvement": float(1.0 - a.mean() / b.mean()),
            "improvement_ci": (float(np.percentile(rel_bs, 2.5)), float(np.percentile(rel_bs, 97.5))),
            "win_rate": float(np.mean(d < 0)), "wilcoxon_p": p}


def group(rows):
    """{(exp, scenario): {filter: {seed: row}}}"""
    g = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        g[(r.get("exp", "T"), r["scenario"])][r["filter"]][int(r["seed"])] = r
    return g


def aligned(d, filters, key):
    seeds = sorted(set.intersection(*[set(d[f]) for f in filters if f in d])) if all(f in d for f in filters) else []
    return seeds, {f: np.array([key(d[f][s]) for s in seeds]) for f in filters}


def rms_h(r):
    return math.hypot(fnum(r, "rms_x"), fnum(r, "rms_y"))


def fmt_imp(ps):
    if ps is None:
        return "-"
    lo, hi = ps["improvement_ci"]
    w = "" if ps["win_rate"] in (0.0, 1.0) else f", win {100 * ps['win_rate']:.0f} %"
    if ps["win_rate"] == 1.0:
        w = ", all seeds"
    elif ps["win_rate"] == 0.0:
        w = ", no seed"
    return f"{100 * ps['improvement']:+.1f} % [{100 * lo:+.1f}, {100 * hi:+.1f}]{w}"


def mm(x):
    return f"{1e3 * x:.2f}" if np.isfinite(x) else "-"


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

class Report:
    def __init__(self, out):
        self.out = out
        self.md = []
        self.stats = {}
        raw = load(out / "raw_runs.csv")
        self.g = group(raw)
        self.tune = group(load(out / "tuning_runs.csv"))
        self.trans = load(out / "transitions.csv")

    def w(self, s=""):
        self.md.append(s)

    def cond(self, exp):
        return [(e, sc) for (e, sc) in self.g if e == exp]

    def metric_table(self, exp, filters, metrics, title, baseline="A41@A", extra=None, scen_label=None):
        self.w(f"### {title}")
        self.w()
        hdr = ["Condition", "n"] + [f"{LABEL.get(f, f)}" for f in filters]
        for _mname, _key, _fmt in metrics:
            pass
        rows = []
        for (e, sc) in sorted(self.cond(exp), key=lambda k: k[1]):
            d = self.g[(e, sc)]
            fl = [f for f in filters if f in d]
            if not fl:
                continue
            for mname, key, fmt in metrics:
                seeds, v = aligned(d, fl, key)
                if not seeds:
                    continue
                cells = []
                for f in filters:
                    if f not in v:
                        cells.append("-")
                        continue
                    cells.append(fmt(np.nanmean(v[f])))
                rows.append([f"{scen_label(sc) if scen_label else sc} ({mname})", str(len(seeds))] + cells)
                if baseline in v:
                    for f in fl:
                        if f == baseline:
                            continue
                        ps = paired(v[f], v[baseline])
                        if ps:
                            self.stats.setdefault(exp, {}).setdefault(sc, {}).setdefault(mname, {})[f] = ps
        self.w("| " + " | ".join(hdr) + " |")
        self.w("|" + "---|" * len(hdr))
        for r in rows:
            self.w("| " + " | ".join(r) + " |")
        self.w()

    def improvement_table(self, exp, cands, baseline, metric, key, title, scen_label=None):
        self.w(f"### {title}")
        self.w()
        self.w("| Condition | n | " + " | ".join(f"{LABEL.get(c, c)} vs {LABEL.get(baseline, baseline)}" for c in cands) + " |")
        self.w("|" + "---|" * (2 + len(cands)))
        for (e, sc) in sorted(self.cond(exp), key=lambda k: k[1]):
            d = self.g[(e, sc)]
            if baseline not in d:
                continue
            cells = []
            n = 0
            for c in cands:
                if c not in d:
                    cells.append("-")
                    continue
                seeds, v = aligned(d, [c, baseline], key)
                ps = paired(v[c], v[baseline])
                n = len(seeds)
                cells.append(fmt_imp(ps))
                if ps:
                    self.stats.setdefault(exp, {}).setdefault(sc, {}).setdefault(metric, {})[f"{c}_vs_{baseline}"] = ps
            self.w(f"| {scen_label(sc) if scen_label else sc} | {n} | " + " | ".join(cells) + " |")
        self.w()

    # -- Tuning ------------------------------------------------------------
    def tuning(self):
        if not self.tune:
            return
        self.w("## Tuning (seeds 101-110, not used for validation)")
        self.w()
        self.w("Mean heave / horizontal RMS error (mm) of the evidence mode and its parameter variants on the tuning "
               "conditions.  Defaults: u = 0.5, z = (2, 4), zf = (8, 16), k = (2, 2.5).  The last column is the "
               "geometric mean over conditions of each variant's (3-D RMS + 1 mm) relative to the default's.")
        self.w()
        filters = None
        tab = []
        ratios = defaultdict(list)
        for (_e, sc), d in sorted(self.tune.items(), key=lambda k: k[0][1]):
            if filters is None:
                filters = sorted(d.keys(), key=lambda f: (f != "A41@A", f != "A41@B", f != "A41@C", f))
            seeds, vz = aligned(d, filters, lambda r: fnum(r, "rms_z"))
            _, vh = aligned(d, filters, rms_h)
            tab.append((sc, len(seeds), {f: (np.mean(vz[f]), np.mean(vh[f])) for f in filters}))
            # 1 mm is added to both sides so still-water rows, where the
            # default reads 0.0 mm, do not dominate the ratio.
            ref = np.mean(np.hypot(vz["A41@C"], vh["A41@C"])) + 1e-3
            for f in filters:
                ratios[f].append((np.mean(np.hypot(vz[f], vh[f])) + 1e-3) / ref)
        self.w("| Condition | n | " + " | ".join(filters) + " |")
        self.w("|" + "---|" * (2 + len(filters)))
        for sc, n, v in tab:
            self.w(f"| {sc} | {n} | " + " | ".join(f"{1e3 * v[f][0]:.1f} / {1e3 * v[f][1]:.1f}" for f in filters) + " |")
        self.w("| geometric mean 3-D ratio | | " + " | ".join(f"{np.exp(np.mean(np.log(ratios[f]))):.3f}" for f in filters) + " |")
        self.w()
        self.stats["tuning_geomean_ratio"] = {f: float(np.exp(np.mean(np.log(ratios[f])))) for f in filters}

    def write(self):
        (self.out / "summary.md").write_text("\n".join(self.md) + "\n")
        def clean(o):
            if isinstance(o, dict):
                return {k: clean(v) for k, v in o.items()}
            if isinstance(o, (list, tuple)):
                return [clean(v) for v in o]
            if isinstance(o, float) and not math.isfinite(o):
                return None
            return o
        (self.out / "paired_stats.json").write_text(json.dumps(clean(self.stats), indent=1, sort_keys=True) + "\n")


def mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    spec = importlib.util.spec_from_file_location("svg_determinism", ROOT / "plots" / "svg_determinism.py")
    svgd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(svgd)
    svgd.configure_svg(matplotlib)
    return plt, svgd


def style(ax):
    ax.grid(True, which="both", color="#e6e6e3", lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def load_trace(path):
    with path.open() as f:
        hdr = f.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return {h: data[:, i] for i, h in enumerate(hdr)}


def main(out):
    out = Path(out)
    r = Report(out)
    r.w("# CARMA wave-excitation study: generated results")
    r.w()
    r.w("Generated by `python3 tools/carma_excitation_study.py report` from the raw CSVs in this directory; "
        "docs/carma-wave-excitation.md is the analysis.  RMS errors in mm over t >= 300 s unless stated; "
        "improvements are paired against the legacy CARMA3 (A41) or, for TFG, the legacy TFG-CARMA3 (T41).")
    r.w()
    r.tuning()
    import carma_excitation_tables as tables
    tables.validation(r)
    r.write()
    plt, svgd = mpl()
    import carma_excitation_figures as figs
    figs.all_figures(r, plt, svgd, out)


if __name__ == "__main__":
    main(ROOT / "reports" / "results" / "carma_excitation")
