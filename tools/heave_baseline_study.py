#!/usr/bin/env python3
"""Paired conventional heave baselines against the shipped estimators.

Replays every record of the OU-III article protocol once per shipped estimator
(OU-III, OU-II, TFG, adaptive PII, TVG-NLO) through tests/hpdi/heave_export-*,
which run each simulator's own code. The OU-III export also carries the
measurement-only inputs the baselines receive: the Mahony-proxy vertical
acceleration a_up (Eq. 73) and the canonical period estimate Tz_hat (Eq. 80).
The baselines (tools/heave_baselines.py: causal HPDI, offline FDDI) are tuned
on a development split and frozen before the evaluation split is scored; every
method is scored by the harness's own scorer, on identical records and seeds.

Splits
  eval  the ten predeclared seed triplets of tools/ou_validation.py on the
        eight stationary RAO records, the controlled crossfade and the
        bidirectional low-high-low record; secondary cells reuse the triplets.
  dev   the eight pinned records under the comparator-retuning sensor and
        initialization draws (reports/results/comparator_rao_retuning/
        selection.json): default, 11, 23, 61001, 62003. These draws and the
        pinned (not phase-randomized) wave histories are disjoint from eval.

Usage
  python3 tools/heave_baseline_study.py run [--work DIR] [--jobs N]
  python3 tools/heave_baseline_study.py analyze [--work DIR]
`run` replays and exports (DIR/records), verifies, tunes and evaluates;
`analyze` rebuilds the tables, statistics and plot from DIR.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
import heave_baselines as hb  # noqa: E402
import ou_roundtrip_transition as rt  # noqa: E402
import ou_validation as ov  # noqa: E402

ROOT = TOOLS.parent
OUT = ROOT / "reports" / "results" / "heave_baselines"
HPDI_TESTS = ROOT / "tests" / "hpdi"
HPDI_BIN = HPDI_TESTS / "hpdi_replay"

FAMILIES = ("OU_III", "OU_II", "TFG", "PII", "TVG_NLO")
LABEL = {"OU_III": "OU-III", "OU_II": "OU-II", "TFG": "TFG", "PII": "PII", "TVG_NLO": "TVG-NLO"}
EXPORT = {f: HPDI_TESTS / f"heave_export-{s}" for f, s in
          zip(FAMILIES, ("ou3", "ou2", "tfg", "pii", "nlo"))}
SHIPPED = {
    "OU_III": ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim",
    "OU_II": ROOT / "tests" / "kalman_ou_ii" / "kalman_ou_ii-sim",
    "TFG": ROOT / "tests" / "kalman_tfg" / "kalman_tfg-sim",
}
# Adaptive settings exactly as tools/ou_validation.py (OU-II/OU-III) and
# tools/tfg_comparison.py (TFG) pass them. PII and TVG-NLO take none.
FAMILY_ENV = {
    "OU_III": {"W3D_TUNING_MODE": "adaptive", "W3D_AW_COV_SYNC": "periodic"},
    "OU_II": {"W3D_TUNING_MODE": "adaptive", "W3D_AW_COV_SYNC": "periodic"},
    "TFG": {"TFG_TUNING": "adaptive", "TFG_AW_COV_SYNC": "1"},
    "PII": {},
    "TVG_NLO": {},
}
ENV_PREFIXES = ("W3D_", "TFG_", "SF_", "OU_II_", "OU_III_", "OU_SIGMA_", "PII_", "NLO_")

DURATION_S = 1200.0
WINDOW_S = 900.0
DEV_DRAWS = ("default", "11", "23", "61001", "62003")
BASELINES = ("HPDI-classic", "HPDI-compensated", "FDDI")
REFERENCES = ("OU-III", "OU-II")
METHOD_ORDER = ("OU-III", "TFG", "OU-II", "PII", "TVG-NLO", "HPDI-classic", "HPDI-compensated",
                "FDDI", "HPDI-classic float32", "HPDI-compensated float32")
STATS_SEED = 20260317
RESAMPLES = 10_000

# Crossfade of tools/ou_validation.py: 1.5 m / 5.7 s -> 4.0 m / 11.4 s, C2
# quintic over 540-660 s of the 1200 s record; scored at the final Hs = 4 m.
CROSSFADE = (ov.TRANSITION_START_FRACTION * DURATION_S, ov.TRANSITION_END_FRACTION * DURATION_S)
RAMP = (585.0, 615.0)  # tools/ou_robustness_core.py "rapid"
ENGINE_CRUISE = {"W3D_ENGINE_RPM": "2400", "W3D_ENGINE_LEVEL_MPS2": "0.60",
                 "W3D_ENGINE_BANDWIDTH_HZ": "80"}  # tools/engine_noise_degradation.py nominal cell


@dataclass(frozen=True)
class Unit:
    name: str
    cohort: str          # stationary, crossfade, lowhighlow, low_motion, ramp, default_draw, noise_free, engine, dev
    split: str           # eval or dev
    family: str          # JONSWAP, PM-Stokes, Crossfade, LowHighLow
    scenario: str        # ou_validation scenario slug where one exists
    seed: str            # pairing key
    kind: str            # stationary, crossfade, ramp, lowhighlow, low_motion, pinned
    source: str
    filename: str
    imu_seed: str = ""
    init_seed: str = ""
    wave_seed: int = 0
    no_noise: bool = False
    env: tuple[tuple[str, str], ...] = ()
    segments: tuple[tuple[str, float, float], ...] = ()


def crossfade_segments(t0: float, t1: float) -> tuple[tuple[str, float, float], ...]:
    """The four intervals tools/ou_validation.py scores inside the final window."""
    start = DURATION_S - WINDOW_S
    recover = min(DURATION_S, t1 + (t1 - t0))
    spans = (("start", start, t0), ("blend", t0, t1), ("recover", t1, recover), ("end", recover, DURATION_S))
    return tuple((n, max(a, start), min(b, DURATION_S)) for n, a, b in spans if min(b, DURATION_S) > max(a, start))


def spectrum_family(path: Path) -> str:
    return "PM-Stokes" if "pmstokes" in path.name else "JONSWAP"


def height_of(name: str) -> float:
    return float(re.search(r"_H([0-9.]+)_", name).group(1))


def build_units(data_dir: Path, secondary: bool) -> list[Unit]:
    triplets = ov.broadcast_seed_triplets(ov.DEFAULT_FULL_WAVE_SEEDS, ov.DEFAULT_FULL_IMU_SEEDS,
                                          ov.DEFAULT_FULL_INIT_SEEDS)
    stationary = sorted(data_dir.glob("wave_data_jonswap_*.csv")) + sorted(data_dir.glob("wave_data_pmstokes_*.csv"))
    if len(stationary) != 8:
        raise SystemExit(f"expected the eight pinned RAO records in {data_dir}")
    nominal = ov.find_default_input(data_dir, "1.500", "50.710")
    low = ov.find_default_input(data_dir, "0.270", "14.047")
    units: list[Unit] = []
    for seed in triplets:
        key = f"{seed.wave_phase_seed}-{seed.imu_noise_seed}-{seed.initialization_seed}"
        common = dict(split="eval", seed=key, imu_seed=str(seed.imu_noise_seed),
                      init_seed=str(seed.initialization_seed), wave_seed=seed.wave_phase_seed)
        for path in stationary:
            units.append(Unit(name=f"eval_{ov._scenario_slug(path)}_{key}", cohort="stationary",
                              family=spectrum_family(path), scenario=ov._scenario_slug(path), kind="stationary",
                              source=str(path), filename=path.name, **common))
        units.append(Unit(name=f"eval_crossfade_{key}", cohort="crossfade", family="Crossfade",
                          scenario="nonstationary_H1_5_to_H4_0_Tp5_7_to_11_4", kind="crossfade",
                          source=str(nominal), filename="wave_data_jonswap_H4.000_L202.839_A-30.00_P120.00.csv",
                          segments=crossfade_segments(*CROSSFADE), **common))
        units.append(Unit(name=f"eval_lowhighlow_{key}", cohort="lowhighlow", family="LowHighLow",
                          scenario="roundtrip_H1_5_H4_0_H1_5", kind="lowhighlow", source=str(nominal),
                          filename=rt.ROUNDTRIP_WAVE_NAME, segments=rt.roundtrip_segments(), **common))
        if secondary:
            units.append(Unit(name=f"eval_lowmotion_{key}", cohort="low_motion", family="JONSWAP",
                              scenario="low_motion_H0_050", kind="low_motion", source=str(low),
                              filename=low.name.replace("_H0.270_", "_H0.050_"), **common))
            units.append(Unit(name=f"eval_ramp_{key}", cohort="ramp", family="Crossfade",
                              scenario="rapid_H1_5_to_H4_0", kind="ramp", source=str(nominal),
                              filename="wave_data_jonswap_H4.000_L202.839_A-30.00_P120.00.csv",
                              segments=crossfade_segments(*RAMP), **common))
    for draw in DEV_DRAWS:
        for path in stationary:
            seeds = {} if draw == "default" else {"imu_seed": draw, "init_seed": draw}
            units.append(Unit(name=f"dev_{ov._scenario_slug(path)}_{draw}", cohort="dev", split="dev",
                              family=spectrum_family(path), scenario=ov._scenario_slug(path), seed=f"draw-{draw}",
                              kind="pinned", source=str(path), filename=path.name, **seeds))
    for path in stationary:
        units.append(Unit(name=f"default_{ov._scenario_slug(path)}", cohort="default_draw", split="eval",
                          family=spectrum_family(path), scenario=ov._scenario_slug(path), seed="default",
                          kind="pinned", source=str(path), filename=path.name))
        if secondary:
            units.append(Unit(name=f"noisefree_{ov._scenario_slug(path)}", cohort="noise_free", split="eval",
                              family=spectrum_family(path), scenario=ov._scenario_slug(path), seed="no-noise",
                              kind="pinned", source=str(path), filename=path.name, no_noise=True))
            units.append(Unit(name=f"engine_{ov._scenario_slug(path)}", cohort="engine", split="eval",
                              family=spectrum_family(path), scenario=ov._scenario_slug(path), seed="default",
                              kind="pinned", source=str(path), filename=path.name,
                              env=tuple(sorted(ENGINE_CRUISE.items()))))
    return units


def generate(unit: Unit) -> tuple[list[str], np.ndarray]:
    source = Path(unit.source)
    columns, data = ov.read_wave_csv(source, DURATION_S)
    if unit.kind == "stationary":
        return columns, ov.phase_randomize_wave(columns, data, unit.wave_seed)
    if unit.kind == "low_motion":
        randomized = ov.phase_randomize_wave(columns, data, unit.wave_seed)
        return columns, ov.scale_wave_motion(columns, randomized, 0.05 / 0.27)
    end_path = ov.find_default_input(source.parent, "8.500", "202.839")
    _, end = ov.read_wave_csv(end_path, DURATION_S)
    scale = 4.0 / height_of(end_path.name)
    if unit.kind == "lowhighlow":
        return columns, rt.make_roundtrip_wave(ov, columns, data, end, unit.wave_seed, scale)
    t0, t1 = CROSSFADE if unit.kind == "crossfade" else RAMP
    return columns, ov.make_nonstationary_wave(columns, data, end, unit.wave_seed, scale, t0, t1)


def unit_env(unit: Unit, family: str) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith(ENV_PREFIXES)}
    env.update(W3D_WRITE_TIMESERIES="0", W3D_COLLECT_ALL_GATES="1",
               W3D_VALIDATION_WINDOW_SEC=f"{WINDOW_S:.9g}",
               W3D_VALIDATION_SEGMENTS=",".join(f"{n}:{a:.9g}:{b:.9g}" for n, a, b in unit.segments))
    if unit.imu_seed:
        env["W3D_IMU_SEED"] = unit.imu_seed
    if unit.init_seed:
        env["W3D_INIT_SEED"] = unit.init_seed
    env.update(FAMILY_ENV[family])
    env.update(dict(unit.env))
    return env


def metric_lines(stdout: str) -> list[str]:
    return [line for line in stdout.splitlines() if line.startswith(("VALIDATION_METRICS ", "VALIDATION_SEGMENT "))]


def run_unit(unit: Unit, work: Path, check_shipped: bool) -> dict[str, Any]:
    record = work / "records" / f"{unit.name}.csv"
    meta_path = work / "records" / f"{unit.name}.json"
    if record.exists() and meta_path.exists():
        return json.loads(meta_path.read_text())
    scratch = work / "scratch" / unit.name
    scratch.mkdir(parents=True, exist_ok=True)
    wave = scratch / unit.filename
    if unit.kind == "pinned":
        # The pinned file itself, byte for byte, as the deterministic tables
        # and the comparator retuning replayed it.
        shutil.copyfile(unit.source, wave)
        columns, data = ov.read_wave_csv(wave)
    else:
        columns, data = generate(unit)
        ov.write_wave_csv(wave, columns, data)
    # The harness steps every estimator at dt = 1/200 s and never reads the
    # record's time column (which the pinned records print with six
    # significant digits, so it repeats after 1000 s). The baselines get the
    # same sample clock.
    times = np.arange(data.shape[0]) / 200.0
    del data
    meta: dict[str, Any] = {"unit": asdict(unit), "families": {}}
    columns_text: dict[str, list[list[str]]] = {}
    for family in FAMILIES:
        out = scratch / f"{family}.csv"
        args = ["--input", str(wave)] + (["--no-noise"] if unit.no_noise else [])
        env = unit_env(unit, family)
        done = subprocess.run([str(EXPORT[family]), "--out", str(out), *args], cwd=EXPORT[family].parent, env=env,
                              text=True, capture_output=True)
        if done.returncode not in (0, 1):
            raise RuntimeError(f"{unit.name} {family}: exit {done.returncode}\n{done.stderr[-2000:]}")
        lines = metric_lines(done.stdout)
        hs = float(re.search(r"HEAVE_EXPORT family=\S+ samples=\d+ hs=(\S+)", done.stdout).group(1))
        shipped_identical = None
        if check_shipped and family in SHIPPED:
            ref = subprocess.run([str(SHIPPED[family]), *args], cwd=SHIPPED[family].parent, env=env, text=True,
                                 capture_output=True)
            shipped_identical = metric_lines(ref.stdout) == lines and bool(lines)
            if not shipped_identical:
                raise RuntimeError(f"{unit.name} {family}: export metrics differ from the shipped simulator")
        metrics = ov.parse_validation_metrics(done.stdout)
        meta["families"][family] = {"hs": hs, "metrics": metrics, "shipped_identical": shipped_identical,
                                    "metric_lines_sha256": hashlib.sha256("\n".join(lines).encode()).hexdigest()}
        with open(out) as fh:
            columns_text[family] = [line.rstrip("\n").split(",") for line in fh]
    # Merge: t from the record, a_up/Tz_hat and z_ref from the OU-III export.
    ou3 = columns_text["OU_III"]
    head = ou3[0]
    z_ref = [row[head.index("z_ref")] for row in ou3]
    for family, rows in columns_text.items():
        if [row[0] for row in rows] != z_ref:
            raise RuntimeError(f"{unit.name}: z_ref of {family} differs from OU-III's")
    if len(times) != len(ou3) - 1:
        raise RuntimeError(f"{unit.name}: record has {len(times)} samples, export {len(ou3) - 1}")
    tmp = record.with_suffix(".tmp")
    tmp.parent.mkdir(parents=True, exist_ok=True)
    with open(tmp, "w") as fh:
        header = ["t", "a_up", "Tz_hat", "z_ref"] + [f"{p}_{f}" for f in FAMILIES for p in ("err", "zhat")]
        fh.write(",".join(header) + "\n")
        ia, it = head.index("a_up"), head.index("Tz_hat")
        for k in range(1, len(ou3)):
            cells = [f"{times[k - 1]:.6f}", ou3[k][ia], ou3[k][it], z_ref[k]]
            for family in FAMILIES:
                cells += columns_text[family][k][1:3]
            fh.write(",".join(cells) + "\n")
    tmp.replace(record)
    meta["record"] = str(record.relative_to(work))
    meta["hs"] = meta["families"]["OU_III"]["hs"]
    meta_path.write_text(json.dumps(meta, indent=1, sort_keys=True))
    shutil.rmtree(scratch)
    return meta


# --------------------------------------------------------------------------- verification


def rescore_record(args) -> dict[str, Any]:
    path, hs, segments = args
    with open(path) as fh:
        header = fh.readline().strip().split(",")
    names = [h for h in header if h.startswith("err_")]
    use = [header.index(h) for h in ("a_up", "z_ref")] + [header.index(h) for h in names]
    data = np.loadtxt(path, delimiter=",", skiprows=1, usecols=use)
    a_up, z_ref = data[:, 0], data[:, 1]
    out: dict[str, Any] = {}
    for i, name in enumerate(names):
        err = data[:, 2 + i].astype(np.float32)
        k = hb.harness_count(err.size, WINDOW_S)
        out[name[4:]] = {"": hb.pct_hs(hb.harness_rms(err[-k:]), hs)}
        for seg, t0, t1 in segments:
            part = err[min(err.size, hb.harness_index(t0)):min(err.size, hb.harness_index(t1))]
            out[name[4:]][seg] = hb.pct_hs(hb.harness_rms(part), hs)
    # Sign convention: a_up must track the second derivative of z_ref.
    acc_ref = np.gradient(np.gradient(z_ref, 0.005), 0.005)
    tail = slice(-hb.harness_count(z_ref.size, WINDOW_S), None)
    out["_sign_corr"] = float(np.corrcoef(a_up[tail], acc_ref[tail])[0, 1])
    return out


def verify(work: Path, metas: list[dict[str, Any]], jobs: int) -> dict[str, Any]:
    tasks = [(work / m["record"], m["hs"], [tuple(s) for s in m["unit"]["segments"]]) for m in metas]
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        scores = list(ex.map(rescore_record, tasks, chunksize=4))
    worst = 0.0
    checked = 0
    sign = []
    for meta, score in zip(metas, scores):
        sign.append(score.pop("_sign_corr"))
        for family, by_segment in score.items():
            metrics = meta["families"][family]["metrics"]
            for seg, value in by_segment.items():
                key = "disp_z_pct_hs" if not seg else f"seg_{seg}_disp_z_pct_hs"
                harness = float(metrics[key])
                rel = abs(value / harness - 1.0)
                worst = max(worst, rel)
                checked += 1
                if np.float32(value) != np.float32(harness):
                    raise RuntimeError(f"{meta['unit']['name']} {family} {seg or 'window'}: "
                                       f"rescored {value!r} != harness {harness!r}")
    return {"scores_checked": checked, "max_relative_difference": worst,
            "exact_float32_equality": True, "a_up_vs_d2z_ref_correlation_min": float(min(sign)),
            "z_ref_sign": "up-positive (no --flip-z)" if min(sign) > 0.5 else "CHECK"}


def committed_comparison(metas: list[dict[str, Any]], reproduced_raw: Path | None) -> dict[str, Any]:
    """Relative deviation of the replayed metrics from committed evidence rows."""
    report: dict[str, Any] = {}

    def compare(name: str, pairs: list[tuple[float, float]]) -> None:
        rel = [abs(a / b - 1.0) for a, b in pairs if b]
        report[name] = {"n": len(rel), "max_relative_difference": max(rel) if rel else None,
                        "identical": sum(1 for a, b in pairs if np.float32(a) == np.float32(b))}

    by_key = {}
    for m in metas:
        u = m["unit"]
        if u["cohort"] in ("stationary", "crossfade"):
            by_key[(u["scenario"], u["seed"])] = m
    for label, path in (("ou_validation committed", ROOT / "reports/results/ou_validation/ou_validation_raw.csv"),
                        ("ou_validation reproduced on this machine", reproduced_raw)):
        if path is None or not Path(path).exists():
            continue
        for family in ("OU_III", "OU_II"):
            pairs = []
            for r in csv.DictReader(open(path)):
                if r["family"] != family or r["mode"] != "Adaptive":
                    continue
                m = by_key.get((r["scenario"], f"{r['wave_phase_seed']}-{r['imu_noise_seed']}-{r['initialization_seed']}"))
                if m is None:
                    continue
                metrics = m["families"][family]["metrics"]
                for key in [k for k in r if k == "disp_z_pct_hs" or (k.startswith("seg_") and k.endswith("disp_z_pct_hs"))]:
                    if r[key] not in ("", None) and key in metrics:
                        pairs.append((float(metrics[key]), float(r[key])))
            compare(f"{label} {family}", pairs)
    pairs = []
    for r in csv.DictReader(open(ROOT / "reports/results/tfg_comparison/raw_runs.csv")):
        m = by_key.get((r["scenario"], f"{r['wave_phase_seed']}-{r['imu_noise_seed']}-{r['initialization_seed']}"))
        if m is not None:
            pairs.append((float(m["families"][r["family"]]["metrics"]["disp_z_pct_hs"]), float(r["disp_z_pct_hs"])))
    compare("tfg_comparison committed OU_III+TFG", pairs)
    robustness = ROOT / "reports/results/ou_robustness/ou_robustness_raw.csv"
    pairs = []
    cohort = {"low_motion": ("low_motion", "Hs0.05"), "ramp": ("transition_rate", "rapid")}
    lookup = {(cohort[m["unit"]["cohort"]], m["unit"]["seed"]): m for m in metas if m["unit"]["cohort"] in cohort}
    for r in csv.DictReader(open(robustness)):
        if r.get("mode") != "Adaptive":
            continue
        key = ((r["experiment"], r["case"]), f"{r['wave_phase_seed']}-{r['imu_noise_seed']}-{r['initialization_seed']}")
        if key in lookup:
            pairs.append((float(lookup[key]["families"]["OU_III"]["metrics"]["disp_z_pct_hs"]), float(r["disp_z_pct_hs"])))
    compare("ou_robustness committed OU_III (Hs 0.05, 30 s ramp)", pairs)
    pairs = []
    rows = {r["segment"]: r for r in csv.DictReader(open(ROOT / "reports/results/ou_rs_law/ou_rs_roundtrip_scores.csv"))}
    for m in metas:
        if m["unit"]["cohort"] == "lowhighlow" and m["unit"]["seed"] == "11-101-1009":
            metrics = m["families"]["OU_III"]["metrics"]
            for seg, r in rows.items():
                pairs.append((float(metrics[f"seg_{seg}_disp_z_pct_hs"]), float(r["disp_z_pct_hs"])))
    compare("ou_rs_law committed low-high-low OU_III (seed 11-101-1009)", pairs)
    return report


def leakage_audit() -> dict[str, Any]:
    """Where reference motion is read by the baseline code paths."""
    hits = {}
    for rel in ("src/hpdi/HeaveHPDI.h", "tools/hpdi_replay.cpp", "tools/heave_baselines.py"):
        text = (ROOT / rel).read_text().splitlines()
        hits[rel] = [f"{i + 1}: {line.strip()}" for i, line in enumerate(text)
                     if re.search(r"\b(z_ref|rec\.z|r\.z|disp_|\.z\[)|\bz_col\b", line)]
    return hits


# --------------------------------------------------------------------------- baselines


def write_manifest(work: Path, metas: list[dict[str, Any]]) -> Path:
    path = work / "manifest.csv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["path", "family", "hs", "split", "seed", "cohort", "scenario", "unit"])
        for m in metas:
            u = m["unit"]
            w.writerow([m["record"], u["family"], f"{m['hs']:.9g}", u["split"], u["seed"], u["cohort"],
                        u["scenario"], u["name"]])
    return path


def baseline_record_scores(args) -> dict[str, Any]:
    """Frozen baselines on one record: whole window and segments, harness scorer."""
    path, hs, segments, frozen, float32 = args
    rec = hb.load_record(Path(path), hb.Columns())
    out: dict[str, Any] = {}
    estimates: dict[str, np.ndarray] = {}
    fcfg = frozen["FDDI"]
    estimates["FDDI"] = hb.FDDI(rec.t, rec.a, frozen["protocol"]["fddi_pad_s"]).displacement(
        hb.fddi_cutoff(fcfg, rec, WINDOW_S), fcfg["taper"])
    for name in ("HPDI-classic", "HPDI-compensated"):
        cfg = frozen[name]
        for suffix, flag in (("", []), (" float32", ["--float"])):
            if suffix and not float32:
                continue
            cmd = [str(HPDI_BIN), "--in", str(path), "--dump", "-", "--n", str(cfg["n"]), "--m", str(cfg["m"]),
                   "--mode", cfg["mode"], "--value", repr(cfg["value"]),
                   "--smoothing-periods", repr(frozen["protocol"]["hpdi_smoothing_periods"]),
                   "--refresh", repr(frozen["protocol"]["hpdi_refresh_s"]),
                   "--fallback-period", repr(frozen["protocol"]["hpdi_fallback_period_s"]), *flag]
            res = subprocess.run(cmd, capture_output=True, text=True, check=True)
            estimates[name + suffix] = np.loadtxt(res.stdout.splitlines()[1:], delimiter=",", usecols=1)
    for method, z_hat in estimates.items():
        k = hb.harness_count(rec.z.size, WINDOW_S)
        rms = hb.harness_rms(hb.float_error(z_hat[-k:], rec.z[-k:]))
        out[method] = {"": (hb.pct_hs(rms, hs), rms)}
        for seg, t0, t1 in segments:
            r = hb.segment_rms(z_hat, rec.z, t0, t1)
            out[method][seg] = (hb.pct_hs(r, hs), r)
    return out


def run_baselines(work: Path, metas: list[dict[str, Any]], jobs: int, float32: bool) -> list[dict[str, Any]]:
    frozen = json.loads((work / "frozen.json").read_text())
    evals = [m for m in metas if m["unit"]["split"] == "eval"]
    tasks = [(str(work / m["record"]), m["hs"], [tuple(s) for s in m["unit"]["segments"]], frozen, float32)
             for m in evals]
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        results = list(ex.map(baseline_record_scores, tasks, chunksize=2))
    rows = []
    for m, res in zip(evals, results):
        for method, by_seg in res.items():
            for seg, (pct, rms) in by_seg.items():
                rows.append(long_row(m, method, seg, pct, rms))
    return rows


def long_row(meta: dict[str, Any], method: str, segment: str, pct: float, rms: float) -> dict[str, Any]:
    u = meta["unit"]
    return {"method": method, "family": u["family"], "hs": f"{meta['hs']:.9g}", "seed": u["seed"],
            "cohort": u["cohort"], "split": u["split"], "segment": segment, "pct_hs": f"{pct:.9g}",
            "rms_m": f"{rms:.9g}", "scenario": u["scenario"], "record": meta["record"]}


def existing_rows(metas: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for m in metas:
        for family in FAMILIES:
            metrics = m["families"][family]["metrics"]
            rows.append(long_row(m, LABEL[family], "", float(metrics["disp_z_pct_hs"]),
                                 float(metrics["disp_z_rms_m"])))
            for seg, _, _ in m["unit"]["segments"]:
                rows.append(long_row(m, LABEL[family], seg, float(metrics[f"seg_{seg}_disp_z_pct_hs"]),
                                     float(metrics[f"seg_{seg}_disp_z_rms_m"])))
    return rows


# --------------------------------------------------------------------------- statistics


def seed_means(rows, method: str, cohort: str, family: str, segment: str = "") -> dict[str, float]:
    groups: dict[str, list[float]] = {}
    for r in rows:
        if r["method"] == method and r["cohort"] == cohort and r["family"] == family and r["segment"] == segment:
            groups.setdefault(r["seed"], []).append(float(r["pct_hs"]))
    counts = {len(v) for v in groups.values()}
    if len(counts) > 1:
        raise ValueError(f"unbalanced records for {method} {cohort} {family}")
    return {k: float(np.mean(v)) for k, v in groups.items()}


def contrast(rows, left: str, right: str, cohort: str, family: str, segment: str = "") -> dict[str, Any]:
    a = seed_means(rows, left, cohort, family, segment)
    b = seed_means(rows, right, cohort, family, segment)
    keys = sorted(set(a) & set(b))
    if set(a) != set(b) or len(keys) < 2:
        raise ValueError(f"contrast {left} - {right} on {cohort}/{family} is not paired")
    d = np.array([a[k] - b[k] for k in keys])
    rng = np.random.default_rng(STATS_SEED)
    low, high = ov._bootstrap_mean_ci(d, RESAMPLES, rng)
    inference = ov.paired_inference(d, np.random.default_rng(STATS_SEED))
    return {"left": left, "right": right, "cohort": cohort, "family": family, "segment": segment,
            "n_seeds": len(keys), "left_mean": float(np.mean([a[k] for k in keys])),
            "right_mean": float(np.mean([b[k] for k in keys])), "mean_difference": float(np.mean(d)),
            "bootstrap_ci95_low": low, "bootstrap_ci95_high": high,
            "t_ci95_low": inference["t_ci95_low"], "t_ci95_high": inference["t_ci95_high"],
            "t_p_value": inference["t_p_value"], "sign_flip_p_value": inference["randomization_p_value"],
            "sign_flip_exact": inference["randomization_exact"], "negative": inference["sign_negative"],
            "positive": inference["sign_positive"]}


def article_primary_check(metas: list[dict[str, Any]]) -> dict[str, Any]:
    """The article's own aggregate routine on the replayed OU-III/OU-II rows."""
    rows = []
    for m in metas:
        u = m["unit"]
        if u["cohort"] != "stationary":
            continue
        w, i, s = (int(x) for x in u["seed"].split("-"))
        for family in ("OU_II", "OU_III"):
            rows.append({"scenario": u["scenario"], "mode": "Adaptive", "family": family, "wave_phase_seed": w,
                         "imu_noise_seed": i, "initialization_seed": s,
                         "disp_z_pct_hs": m["families"][family]["metrics"]["disp_z_pct_hs"]})
    out = {}
    for spectrum in ("jonswap", "pmstokes"):
        agg = ov.stationary_normalized_aggregate(rows, RESAMPLES, STATS_SEED, spectrum)
        d = agg["OU_III_minus_OU_II"]
        out[spectrum] = {"OU_III": [agg["OU_III"]["mean"], agg["OU_III"]["std"]],
                         "OU_II": [agg["OU_II"]["mean"], agg["OU_II"]["std"]],
                         "difference": d["mean_paired_difference"],
                         "bootstrap_ci95": [d["bootstrap_ci95_low"], d["bootstrap_ci95_high"]],
                         "sign_flip_p": d["randomization_p_value"]}
    return out


# --------------------------------------------------------------------------- reporting


def fmt(x: float, nd: int = 2) -> str:
    return f"{x:.{nd}f}"


def per_sea_table(rows, methods, cohort: str, families=("JONSWAP", "PM-Stokes"), segment: str = "") -> list[str]:
    lines = ["| Spectrum | Hs [m] | " + " | ".join(methods) + " |",
             "|---|---:|" + "---:|" * len(methods)]
    for family in families:
        heights = sorted({float(r["hs"]) for r in rows if r["cohort"] == cohort and r["family"] == family})
        for hs in heights:
            cells = []
            for method in methods:
                v = [float(r["pct_hs"]) for r in rows if r["method"] == method and r["cohort"] == cohort
                     and r["family"] == family and float(r["hs"]) == hs and r["segment"] == segment]
                cells.append("–" if not v else (fmt(v[0]) if len(v) == 1 else
                                                f"{np.mean(v):.2f} ± {np.std(v, ddof=1):.2f}"))
            lines.append(f"| {family} | {hs:.2f} | " + " | ".join(cells) + " |")
    return lines


def contrast_table(cs: list[dict[str, Any]]) -> list[str]:
    lines = ["| Contrast | Ensemble | n | Left | Right | Difference [pp] | Bootstrap 95% CI | t 95% CI | sign-flip p | t p |",
             "|---|---|---:|---:|---:|---:|---|---|---:|---:|"]
    for c in cs:
        ens = c["family"] + (f" / {c['segment']}" if c["segment"] else "") + (
            "" if c["cohort"] in ("stationary",) else f" ({c['cohort']})")
        lines.append(f"| {c['left']} − {c['right']} | {ens} | {c['n_seeds']} | {c['left_mean']:.3f} | "
                     f"{c['right_mean']:.3f} | {c['mean_difference']:+.3f} | "
                     f"[{c['bootstrap_ci95_low']:+.3f}, {c['bootstrap_ci95_high']:+.3f}] | "
                     f"[{c['t_ci95_low']:+.3f}, {c['t_ci95_high']:+.3f}] | {c['sign_flip_p_value']:.4f} | "
                     f"{c['t_p_value']:.2g} |")
    return lines


def write_plot(rows, path: Path) -> None:
    """Fig. 6 of the article (ov.write_metric_plot), with the baselines added."""
    summary = []
    rng = np.random.default_rng(STATS_SEED)
    methods = [m for m in METHOD_ORDER if "float32" not in m]
    for family, tag in (("JONSWAP", "jonswap"), ("PM-Stokes", "pmstokes")):
        for hs in sorted({float(r["hs"]) for r in rows if r["cohort"] == "stationary" and r["family"] == family}):
            h = f"{hs:.3f}".replace(".", "_")
            scenario = f"stationary_{tag}_H{h}_"
            for method in methods:
                v = np.array([float(r["pct_hs"]) for r in rows if r["method"] == method and r["cohort"] == "stationary"
                              and r["family"] == family and float(r["hs"]) == hs and r["segment"] == ""])
                lo, hi = ov._bootstrap_mean_ci(v, RESAMPLES, rng)
                summary.append({"metric": "disp_z_pct_hs", "mode": "", "scenario": scenario,
                                "family": method + (" (offline)" if method == "FDDI" else ""),
                                "mean": float(v.mean()), "bootstrap_ci95_low": lo, "bootstrap_ci95_high": hi})
    ov.write_metric_plot(path, summary, "disp_z_pct_hs", r"Vertical RMS error [% $H_s$]", modes=("",))


def write_long_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["method", "family", "hs", "seed", "pct_hs", "rms_m", "cohort", "split",
                                           "segment", "scenario", "record"], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def analyze(work: Path, out: Path) -> dict[str, Any]:
    rows = list(csv.DictReader(open(out / "heave_baselines_long.csv")))
    contrasts = []
    for cohort, family, segment in (("stationary", "JONSWAP", ""), ("stationary", "PM-Stokes", ""),
                                    ("crossfade", "Crossfade", ""), ("lowhighlow", "LowHighLow", "")):
        for left in BASELINES + ("TFG", "PII", "TVG-NLO"):
            for right in REFERENCES:
                if left in ("TFG", "PII", "TVG-NLO") and right == "OU-II":
                    continue
                contrasts.append(contrast(rows, left, right, cohort, family, segment))
        if cohort == "stationary" and family == "JONSWAP":
            contrasts.append(contrast(rows, "OU-III", "OU-II", cohort, family, segment))
            contrasts.append(contrast(rows, "HPDI-compensated", "HPDI-classic", cohort, family, segment))
            contrasts.append(contrast(rows, "FDDI", "HPDI-compensated", cohort, family, segment))
    for cohort, family in (("low_motion", "JONSWAP"), ("ramp", "Crossfade")):
        if any(r["cohort"] == cohort for r in rows):
            for left in BASELINES:
                contrasts.append(contrast(rows, left, "OU-III", cohort, family))
    (out / "contrasts.json").write_text(json.dumps(contrasts, indent=1))
    return {"contrasts": contrasts, "rows": rows}


def segment_table(rows, methods, cohort: str, family: str, segments: Sequence[str]) -> list[str]:
    lines = ["| Interval | " + " | ".join(methods) + " |", "|---|" + "---:|" * len(methods)]
    for seg in ("",) + tuple(segments):
        cells = []
        for method in methods:
            v = [float(r["pct_hs"]) for r in rows if r["method"] == method and r["cohort"] == cohort
                 and r["family"] == family and r["segment"] == seg]
            cells.append(f"{np.mean(v):.2f} ± {np.std(v, ddof=1):.2f}" if len(v) > 1 else
                         (f"{v[0]:.2f}" if v else "–"))
        lines.append(f"| {seg or 'final 900 s'} | " + " | ".join(cells) + " |")
    return lines


def deterministic_table(rows, methods, cohort: str) -> list[str]:
    """Table X layout: one record per sea, Z RMS in metres and in % Hs."""
    lines = ["| Spectrum | Hs | " + " | ".join(f"{m} [m]" for m in methods) + " | "
             + " | ".join(f"{m} [%Hs]" for m in methods) + " |",
             "|---|---:|" + "---:|" * (2 * len(methods))]
    for family in ("JONSWAP", "PM-Stokes"):
        for hs in sorted({float(r["hs"]) for r in rows if r["cohort"] == cohort and r["family"] == family}):
            sel = {r["method"]: r for r in rows if r["cohort"] == cohort and r["family"] == family
                   and float(r["hs"]) == hs and r["segment"] == ""}
            m_cells = [f"{float(sel[m]['rms_m']):.3f}" if m in sel else "–" for m in methods]
            p_cells = [f"{float(sel[m]['pct_hs']):.2f}" if m in sel else "–" for m in methods]
            lines.append(f"| {family} | {hs:.2f} | " + " | ".join(m_cells) + " | " + " | ".join(p_cells) + " |")
    return lines


def find(cs, left, right, cohort, family, segment=""):
    for c in cs:
        if (c["left"], c["right"], c["cohort"], c["family"], c["segment"]) == (left, right, cohort, family, segment):
            return c
    return None


def describe_cfg(cfg: dict) -> str:
    unit = "Hz" if cfg["mode"] == "fixed" else "/T_z"
    if "n" in cfg:
        return f"n={cfg['n']}, m={cfg['m']}, {cfg['mode']} cutoff {cfg['value']:.4g} {unit}"
    return f"{cfg['mode']} f1 = {cfg['value']:.4g} {unit}, taper {cfg['taper']}"


def write_tex(out: Path, result: dict[str, Any]) -> None:
    """Rows for docs/heave_baselines.tex-part: ten-seed paired vertical endpoint."""
    rows, cs = result["rows"], result["contrasts"]
    frozen = json.loads((out / "frozen.json").read_text())
    label = {"OU-III": r"OU--III", "OU-II": r"OU--II", "TFG": "TFG", "PII": "Adaptive PII", "TVG-NLO": r"TVG--NLO",
             "HPDI-classic": rf"HPDI (classic, $n={frozen['HPDI-classic']['n']}$)",
             "HPDI-compensated": rf"HPDI (compensated, $n={frozen['HPDI-compensated']['n']},"
                                 rf"m={frozen['HPDI-compensated']['m']}$)",
             "FDDI": r"FDDI (offline, non-causal)"}

    def cell(method: str, family: str) -> str:
        v = np.array(list(seed_means(rows, method, "stationary", family).values()))
        return rf"{v.mean():.2f}$\pm${v.std(ddof=1):.2f}"

    def diff(method: str, family: str) -> str:
        if method == "OU-III":
            return "--"
        c = find(cs, method, "OU-III", "stationary", family)
        return (rf"{c['mean_difference']:+.2f} [{c['bootstrap_ci95_low']:+.2f}, {c['bootstrap_ci95_high']:+.2f}]"
                .replace("+", "$+$").replace("-", "$-$"))

    lines = [r"\begin{table}[t]", r"\centering",
             r"\caption{Ten-seed vertical-displacement RMS error over the final \SI{900}{s} (\% of $H_s$), "
             r"seed-level means over the four stationary seas of each spectrum (mean$\pm$sd, $n=10$), and the "
             r"paired difference from OU--III with its percentile-bootstrap 95\% interval. Baseline parameters "
             r"were frozen on development draws. FDDI is non-causal and is an offline reference.}",
             r"\label{tab:conventional-baselines}", r"\footnotesize", r"\setlength{\tabcolsep}{2.5pt}",
             r"\begin{tabular}{@{}lrrrr@{}}", r"\toprule",
             r"& \multicolumn{2}{c}{JONSWAP} & \multicolumn{2}{c}{PM--Stokes} \\",
             r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}",
             r"Method & \% $H_s$ & $-$OU--III [pp] & \% $H_s$ & $-$OU--III [pp] \\", r"\midrule"]
    for method in ("OU-III", "TFG", "OU-II", "PII", "TVG-NLO", "HPDI-classic", "HPDI-compensated", "FDDI"):
        lines.append(f"{label[method]} & {cell(method, 'JONSWAP')} & {diff(method, 'JONSWAP')} & "
                     f"{cell(method, 'PM-Stokes')} & {diff(method, 'PM-Stokes')} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    (out / "heave_baselines_table.tex").write_text("\n".join(lines) + "\n")


def write_summary(work: Path, out: Path, result: dict[str, Any]) -> None:
    rows, cs = result["rows"], result["contrasts"]
    checks = json.loads((out / "checks.json").read_text())
    frozen = json.loads((out / "frozen.json").read_text())
    main_methods = ["OU-III", "TFG", "OU-II", "PII", "TVG-NLO", "HPDI-classic", "HPDI-compensated", "FDDI"]
    L: list[str] = ["# Conventional heave baselines: HPDI and FDDI against the shipped estimators", "",
                    "Generated by `tools/heave_baseline_study.py`. FDDI is non-causal: every table labels it an "
                    "**offline reference**, not a real-time competitor. All methods are scored on the vertical "
                    "endpoint only, with the simulation harness's own scorer, on identical records and seeds.", ""]
    ap = checks["article_primary"]
    L += ["## Reproduction", "",
          "The article's own aggregate routine (`ou_validation.stationary_normalized_aggregate`) applied to "
          "the replayed OU-III/OU-II rows:", "",
          "| Ensemble | OU-III %Hs | OU-II %Hs | OU-III − OU-II [pp] | bootstrap 95% CI | exact sign-flip p |",
          "|---|---|---|---:|---|---:|"]
    for spectrum, label in (("jonswap", "JONSWAP (primary)"), ("pmstokes", "PM-Stokes")):
        a = ap[spectrum]
        L.append(f"| {label} | {a['OU_III'][0]:.2f} ± {a['OU_III'][1]:.2f} | {a['OU_II'][0]:.2f} ± {a['OU_II'][1]:.2f} | "
                 f"{a['difference']:+.3f} | [{a['bootstrap_ci95'][0]:+.3f}, {a['bootstrap_ci95'][1]:+.3f}] | "
                 f"{a['sign_flip_p']:.4f} |")
    L += ["", "## Verification", "",
          f"- Rescoring: {checks['rescore']['scores_checked']} exported per-sample error series (all five "
          "estimators, every record, whole window and every segment) rescored by `tools/heave_baselines.py` "
          f"reproduce the harness's printed %Hs with exact float32 equality (max relative difference "
          f"{checks['rescore']['max_relative_difference']:.1e}).",
          "- Export on vs off: for every record, the OU-III, OU-II and TFG export builds print VALIDATION_METRICS "
          "and VALIDATION_SEGMENT lines byte-identical to the shipped simulators' "
          f"(all identical: {checks['shipped_metric_lines_identical']}).",
          f"- Reference sign: min correlation of a_up with d²z_ref/dt² over the scored window is "
          f"{checks['rescore']['a_up_vs_d2z_ref_correlation_min']:.3f}; z_ref is "
          f"{checks['rescore']['z_ref_sign']}.",
          "- Against committed evidence (same seeds, same records):", ""]
    L += ["| Evidence | n | identical | max relative difference |", "|---|---:|---:|---:|"]
    for name, c in checks["committed"].items():
        L.append(f"| {name} | {c['n']} | {c['identical']} | "
                 f"{'–' if c['max_relative_difference'] is None else format(c['max_relative_difference'], '.2e')} |")
    L += ["", "Reference motion in the baseline code paths (the only reads are the scoring and selection code):", ""]
    for rel, hits in checks["leakage_audit"].items():
        L.append(f"- `{rel}`: " + ("; ".join(f"`{h}`" for h in hits) if hits else "no reads"))
    proto = frozen["protocol"]
    L += ["", "## Frozen parameters (dev split)", "",
          f"Selected on {len(proto['dev_records'])} dev records (8 pinned records × draws "
          f"{', '.join(DEV_DRAWS)}), objective: mean %Hs. Frozen before the eval split was scored.", "",
          "| Method | Configuration | dev mean %Hs | grid edge |", "|---|---|---:|---|"]
    for name in ("HPDI-classic", "HPDI-compensated", "FDDI", "HPDI-best"):
        cfg = frozen[name]
        L.append(f"| {name}{' (offline reference)' if name == 'FDDI' else ''} | {describe_cfg(cfg)} | "
                 f"{cfg['dev_pct_hs']:.3f} | {'yes' if cfg.get('grid_edge') else 'no'} |")
    L += ["", "## Primary-style paired contrasts", "",
          "Seed-level mean over the four stationary JONSWAP seas; PM-Stokes reported separately, never pooled. "
          "Bootstrap: 10,000 percentile resamples, PCG64 seed 20260317; exact sign-flip (1,024 patterns); "
          "paired Student-t. Differences are left minus right; negative favours the left method.", ""]
    L += contrast_table([c for c in cs if c["cohort"] == "stationary"])
    L += ["", "## Per-sea vertical RMS (mean ± sd over ten seed triplets, % Hs)", "",
          "FDDI = offline reference (non-causal).", ""]
    L += per_sea_table(rows, main_methods, "stationary")
    seg_x = [s for s, _, _ in crossfade_segments(*CROSSFADE)]
    L += ["", "## Controlled crossfade (1.5 m / 5.7 s → 4.0 m / 11.4 s, 540–660 s quintic)", "",
          "Scored at the final Hs = 4 m, as in the article; intervals as `tools/ou_validation.py`.", ""]
    L += segment_table(rows, main_methods, "crossfade", "Crossfade", seg_x)
    L += [""] + contrast_table([c for c in cs if c["cohort"] == "crossfade"])
    seg_r = [s for s, _, _ in rt.roundtrip_segments()]
    L += ["", "## Bidirectional low–high–low (article interval partition, % of Hs = 4 m)", ""]
    L += segment_table(rows, main_methods, "lowhighlow", "LowHighLow", seg_r)
    L += [""] + contrast_table([c for c in cs if c["cohort"] == "lowhighlow"])
    L += ["", "## Secondary cells", ""]
    if any(r["cohort"] == "low_motion" for r in rows):
        L += ["### Hs = 0.05 m low motion (ten triplets)", ""]
        L += segment_table(rows, main_methods, "low_motion", "JONSWAP", ())
        L += [""] + contrast_table([c for c in cs if c["cohort"] == "low_motion"])
        L += ["", "### 30 s ramp (585–615 s, ten triplets)", ""]
        L += segment_table(rows, main_methods, "ramp", "Crossfade", [s for s, _, _ in crossfade_segments(*RAMP)])
        L += [""] + contrast_table([c for c in cs if c["cohort"] == "ramp"])
        L += ["", "### Noise-free mismatch floor (pinned records, `--no-noise`)", ""]
        L += deterministic_table(rows, main_methods, "noise_free")
        L += ["", "### Nominal-cruise engine vibration (2400 rpm, 0.60 m/s², 80 Hz; deployed guard enabled)", ""]
        L += deterministic_table(rows, main_methods, "engine")
    fl = ["HPDI-classic", "HPDI-classic float32", "HPDI-compensated", "HPDI-compensated float32"]
    L += ["", "### HPDI replayed in float32, as deployed (ten triplets, % Hs)", ""]
    L += per_sea_table(rows, fl, "stationary")
    L += ["", "## Deterministic default draw (article Table X layout)", "",
          "The pinned records with the simulators' default sensor draw. This draw is one of the dev draws, so "
          "the baseline entries here are in-sample (as the comparators' were during their retuning).", ""]
    L += deterministic_table(rows, main_methods, "default_draw")
    (out / "summary.md").write_text("\n".join(L) + "\n")


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("command", choices=("run", "analyze"))
    p.add_argument("--work", type=Path, default=ROOT / "runs" / "heave_baselines")
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--data-dir", type=Path, default=ROOT / "plots" / "kalman_ou_ii")
    p.add_argument("--jobs", type=int, default=os.cpu_count() or 2)
    p.add_argument("--no-secondary", action="store_true", help="skip the secondary cells")
    p.add_argument("--no-shipped-check", action="store_true",
                   help="skip the byte comparison against the shipped simulators")
    p.add_argument("--reproduced-raw", type=Path, help="ou_validation_raw.csv of a local full replay")
    p.add_argument("--skip-build", action="store_true")
    a = p.parse_args(argv)
    a.work = a.work.resolve()
    a.out = a.out.resolve()
    a.work.mkdir(parents=True, exist_ok=True)
    a.out.mkdir(parents=True, exist_ok=True)
    units = build_units(a.data_dir.resolve(), not a.no_secondary)

    if a.command == "run":
        if not a.skip_build:
            subprocess.run(["make", "-C", str(HPDI_TESTS), "build", f"-j{a.jobs}"], check=True)
            for binary in SHIPPED.values():
                subprocess.run(["make", "-C", str(binary.parent), binary.name, f"-j{a.jobs}"], check=True)
        with ThreadPoolExecutor(max_workers=a.jobs) as pool:
            futures = [pool.submit(run_unit, u, a.work, not a.no_shipped_check) for u in units]
            metas = []
            for i, f in enumerate(futures, 1):
                metas.append(f.result())
                print(f"[{i}/{len(units)}] {metas[-1]['unit']['name']}", flush=True)
        manifest = write_manifest(a.work, metas)
        checks = {"rescore": verify(a.work, metas, a.jobs),
                  "shipped_metric_lines_identical": all(
                      m["families"][f]["shipped_identical"] for m in metas for f in SHIPPED)
                  if not a.no_shipped_check else None,
                  "committed": committed_comparison(metas, a.reproduced_raw),
                  "article_primary": article_primary_check(metas),
                  "leakage_audit": leakage_audit()}
        (a.work / "checks.json").write_text(json.dumps(checks, indent=1))
        print(json.dumps({k: v for k, v in checks.items() if k != "leakage_audit"}, indent=1))
        tune = subprocess.run([sys.executable, str(TOOLS / "heave_baselines.py"), "tune", "--manifest", str(manifest),
                               "--out", str(a.work / "frozen.json"), "--hpdi-bin", str(HPDI_BIN),
                               "--jobs", str(a.jobs)], capture_output=True, text=True, check=True)
        print(tune.stdout)
        (a.work / "tune.log").write_text(tune.stdout)
        if "grid edge" in tune.stdout:
            raise SystemExit("an optimum is on a grid edge: widen the grid in tools/heave_baselines.py and rerun")
        evaluate = subprocess.run([sys.executable, str(TOOLS / "heave_baselines.py"), "evaluate", "--manifest",
                                   str(manifest), "--params", str(a.work / "frozen.json"), "--out",
                                   str(a.work / "evaluate.csv"), "--hpdi-bin", str(HPDI_BIN), "--jobs", str(a.jobs)],
                                  capture_output=True, text=True, check=True)
        (a.work / "evaluate.log").write_text(evaluate.stdout)
        rows = existing_rows(metas) + run_baselines(a.work, metas, a.jobs, float32=True)
        # The CLI evaluate path and this driver must score the frozen baselines identically.
        cli = {(r["method"], r["path"]): float(r["pct_hs"]) for r in csv.DictReader(open(a.work / "evaluate.csv"))}
        mismatches = [r for r in rows if r["segment"] == "" and (r["method"], str(a.work / r["record"])) in cli
                      and np.float32(cli[(r["method"], str(a.work / r["record"]))]) != np.float32(float(r["pct_hs"]))]
        if mismatches:
            raise SystemExit(f"{len(mismatches)} baseline scores differ between evaluate and the driver")
        write_long_csv(a.out / "heave_baselines_long.csv", rows)
        shutil.copy(a.work / "frozen.json", a.out / "frozen.json")
        (a.out / "checks.json").write_text(json.dumps(checks, indent=1))
    result = analyze(a.work, a.out)
    write_plot(result["rows"], a.out / "heave_baselines_vertical.svg")
    write_summary(a.work, a.out, result)
    write_tex(a.out, result)
    print(f"{len(result['contrasts'])} contrasts -> {a.out / 'contrasts.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
