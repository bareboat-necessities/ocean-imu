#!/usr/bin/env python3
"""Temporary CI harness for the sea-period EMA experiment.

The branch stays based on current main.  This script patches the runner checkout
only, exposing two dimensionless OU-III knobs without changing any existing
sea-scaled coefficients:

  OU_III_ADAPT_TAU_PERIODS    tau/sigma applied EMA horizon / T_sea
  OU_III_TUNER_FREQ_PERIODS   tuner frequency EMA horizon / T_sea

It also removes the 0.1 s candidate-commit decimation: every valid physical
measurement updates the EMA and the resulting candidate is staged for the next
IMU sample.  The 0.1 s cadence remains only on posterior a_w covariance sync.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OU3 = ROOT / "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
TUNER = ROOT / "src/tuner/SeaStateAutoTuner.h"
SIM = ROOT / "tests/kalman_ou_iii/kalman_ou_iii-sim.cpp"
OUTDIR = ROOT / "reports/results/ou_sea_period_ema"


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{label}: expected one anchor, found {n}")
    return text.replace(old, new, 1)


def patch_tuner() -> None:
    text = TUNER.read_text()

    text = replace_once(
        text,
        """        last_dt_freq = -1.0f;\n        alpha_freq = 0.0f;\n        tau_var_sec = 0.0f;\n""",
        """        alpha_freq = 0.0f;\n        tau_freq_applied_sec = tau_freq;\n        tau_var_sec = 0.0f;\n""",
        "tuner reset",
    )

    text = replace_once(
        text,
        """        // Update alpha for frequency smoothing (fixed horizon in seconds)\n        updateAlphaFreq(dt_s);\n""",
        """        // The frequency EMA is updated on every physical sample.  When the\n        // period-scaled mode is enabled, only its horizon changes with the\n        // externally measured sea period; no sample is decimated or skipped.\n        updateAlphaFreq(dt_s, f_input_hz);\n""",
        "tuner updateAlphaFreq call",
    )

    text = replace_once(
        text,
        """    inline void setTauFreq(float t) {\n        tau_freq = std::max(1e-3f, t);\n        last_dt_freq = -1.0f;  // force recompute on next update\n    }\n\n""",
        """    // Fixed-seconds compatibility mode.  Calling this explicitly disables\n    // sea-period scaling for the frequency EMA.\n    inline void setTauFreq(float t) {\n        tau_freq = std::max(1e-3f, t);\n        freq_periods = 0.0f;\n        tau_freq_applied_sec = tau_freq;\n    }\n\n    // Frequency-EMA horizon in periods of the externally supplied wave-band\n    // frequency.  The input is measurement-only in the OU-III wrapper.\n    inline void setFrequencySmoothingPeriods(float k) {\n        if (std::isfinite(k) && k > 0.0f) freq_periods = k;\n    }\n    inline float getFrequencySmoothingPeriods() const { return freq_periods; }\n    inline float getFrequencySmoothingHorizonSec() const { return tau_freq_applied_sec; }\n\n""",
        "tuner setTauFreq",
    )

    text = replace_once(
        text,
        """    float tau_var_sec = 0.0f;       // last horizon used, for inspection\n    float tau_freq  = 1.0f;   // seconds\n    float last_dt_freq  = -1.0f;\n    float alpha_freq = 0.0f;\n\n    DebiasedEMA A_mean, A_sq, A_var;\n    DebiasedEMA Freq_smoothed;\n\n    inline void updateAlphaFreq(float dt_s) {\n        if (dt_s == last_dt_freq) return;\n        last_dt_freq = dt_s;\n        alpha_freq = 1.0f - std::exp(-dt_s / tau_freq);\n    }\n""",
        """    float tau_var_sec = 0.0f;       // last variance horizon used, for inspection\n    float tau_freq = 1.0f;             // fixed-seconds compatibility horizon\n    float freq_periods = 0.0f;         // >0: frequency EMA horizon / T_eff\n    float tau_freq_applied_sec = 1.0f; // last frequency-EMA horizon used\n    float alpha_freq = 0.0f;\n\n    DebiasedEMA A_mean, A_sq, A_var;\n    DebiasedEMA Freq_smoothed;\n\n    inline void updateAlphaFreq(float dt_s, float f_input_hz) {\n        float horizon = tau_freq;\n        if (freq_periods > 0.0f) {\n            float f = f_input_hz;\n            if (!std::isfinite(f)) f = f_min_hz;\n            f = std::max(f_min_hz, std::min(f_max_hz, f));\n            horizon = freq_periods / f;\n        }\n        tau_freq_applied_sec = std::max(1e-3f, horizon);\n        alpha_freq = 1.0f - std::exp(-dt_s / tau_freq_applied_sec);\n    }\n""",
        "tuner private frequency EMA",
    )

    TUNER.write_text(text)


def patch_ou3() -> None:
    text = OU3.read_text()

    text = replace_once(
        text,
        """constexpr float ADAPT_TAU_SEC              = 1.8f;\nconstexpr float ADAPT_EVERY_SECS           = 0.1f;\n""",
        """// Legacy fixed-second horizon is retained for ablation/backward compatibility.\n// A positive period multiplier makes the common tau/sigma EMA self-similar:\n//     tau_ema = ADAPT_TAU_PERIODS * T_sea.\n// Zero selects the legacy ADAPT_TAU_SEC path.  The experiment overrides these\n// at runtime; the final measured value is committed only after confirmation.\nconstexpr float ADAPT_TAU_SEC              = 1.8f;\nconstexpr float ADAPT_TAU_PERIODS          = 0.0f;\nconstexpr float TUNER_FREQ_EMA_PERIODS     = 0.0f;\n// This cadence belongs only to posterior a_w covariance maintenance.  Online\n// tuner candidates are staged after every valid physical measurement.\nconstexpr float ADAPT_EVERY_SECS           = 0.1f;\n""",
        "OU3 constants",
    )

    text = replace_once(
        text,
        """    void setAdaptationTimeConstants(float tau_sec) {\n        if (std::isfinite(tau_sec) && tau_sec > 0.0f)   adapt_tau_sec_   = tau_sec;\n     }\n\n""",
        """    // Fixed-seconds compatibility mode.  Explicitly selecting it disables\n    // the period-scaled common tau/sigma EMA.\n    void setAdaptationTimeConstants(float tau_sec) {\n        if (std::isfinite(tau_sec) && tau_sec > 0.0f) {\n            adapt_tau_sec_ = tau_sec;\n            adapt_tau_periods_ = 0.0f;\n        }\n    }\n\n    // Common tau/sigma EMA horizon as a fraction/multiple of measured T_sea.\n    // This changes only averaging memory; tau_coeff_, sigma_coeff_, the sigma\n    // variance K-period horizon, and the already swept r_S horizon are untouched.\n    void setAdaptationPeriods(float periods) {\n        if (std::isfinite(periods) && periods > 0.0f) {\n            adapt_tau_periods_ = periods;\n        }\n    }\n\n    // Frequency EMA inside SeaStateAutoTuner, likewise expressed in T_sea.\n    void setTunerFreqSmoothingPeriods(float periods) {\n        if (std::isfinite(periods) && periods > 0.0f) {\n            tuner_freq_ema_periods_ = periods;\n            tuner_.setFrequencySmoothingPeriods(periods);\n        }\n    }\n\n    // Fixed-seconds frequency-EMA compatibility mode for controlled ablations.\n    void setTunerFreqSmoothingTimeConstant(float tau_sec) {\n        if (std::isfinite(tau_sec) && tau_sec > 0.0f) {\n            tuner_freq_ema_periods_ = 0.0f;\n            tuner_.setTauFreq(tau_sec);\n        }\n    }\n\n""",
        "OU3 adaptation setters",
    )

    text = replace_once(
        text,
        """        adapt_mekf(dt, tau_target_, sigma_target_, RS_target_);\n    }\n\n    void adapt_mekf(float dt, float tau_t, float sigma_t, float RS_t) {\n        const float alpha = 1.0f - std::exp(-dt / adapt_tau_sec_);\n\n        const float RS_sec = seastate::common::adaptiveSmoothingHorizonSec(\n            adapt_RS_mult_, tau_t, RS_t, tune_.RS_applied, adapt_RS_slew_log_, dt);\n        const float alpha_RS = 1.0f - std::exp(-dt / RS_sec);\n\n        tune_.tau_applied   += alpha    * (tau_t   - tune_.tau_applied);\n        tune_.sigma_applied += alpha    * (sigma_t - tune_.sigma_applied);\n        tune_.RS_applied    += alpha_RS * (RS_t    - tune_.RS_applied);\n\n        if (time_ - last_adapt_time_sec_ > adapt_every_secs_) {\n            // y_k may change the smoothed candidate here, but the MEKF keeps\n            // the schedule that was active before y_k arrived. Commit this\n            // candidate at the beginning of updateTime(k+1).\n            online_tune_apply_pending_ = true;\n            last_adapt_time_sec_ = time_;\n        }\n    }\n""",
        """        const float sea_period_sec = 1.0f / f_tune;\n        adapt_mekf(dt, tau_target_, sigma_target_, RS_target_, sea_period_sec);\n    }\n\n    void adapt_mekf(float dt, float tau_t, float sigma_t, float RS_t,\n                    float sea_period_sec) {\n        float adapt_sec = adapt_tau_sec_;\n        if (adapt_tau_periods_ > 0.0f &&\n            std::isfinite(sea_period_sec) && sea_period_sec > 0.0f) {\n            adapt_sec = std::max(dt, adapt_tau_periods_ * sea_period_sec);\n        }\n        const float alpha = 1.0f - std::exp(-dt / adapt_sec);\n\n        const float RS_sec = seastate::common::adaptiveSmoothingHorizonSec(\n            adapt_RS_mult_, tau_t, RS_t, tune_.RS_applied, adapt_RS_slew_log_, dt);\n        const float alpha_RS = 1.0f - std::exp(-dt / RS_sec);\n\n        // Every valid physical sample contributes to every EMA.  The candidate\n        // is staged immediately and becomes active at the *next* IMU boundary,\n        // preserving predictability without decimating measurements.\n        tune_.tau_applied   += alpha    * (tau_t   - tune_.tau_applied);\n        tune_.sigma_applied += alpha    * (sigma_t - tune_.sigma_applied);\n        tune_.RS_applied    += alpha_RS * (RS_t    - tune_.RS_applied);\n        online_tune_apply_pending_ = true;\n    }\n""",
        "OU3 adapt_mekf",
    )

    text = replace_once(
        text,
        """    float adapt_tau_sec_          = ADAPT_TAU_SEC;\n    float adapt_RS_mult_          = ADAPT_RS_MULT;\n""",
        """    float adapt_tau_sec_          = ADAPT_TAU_SEC;\n    float adapt_tau_periods_      = ADAPT_TAU_PERIODS;\n    float tuner_freq_ema_periods_ = TUNER_FREQ_EMA_PERIODS;\n    float adapt_RS_mult_          = ADAPT_RS_MULT;\n""",
        "OU3 members",
    )

    OU3.write_text(text)


def patch_sim() -> None:
    text = SIM.read_text()
    text = replace_once(
        text,
        """            if (env_float(\"OU_III_ADAPT_TAU_SEC\", v)) {\n                filter.setAdaptationTimeConstants(v);\n            }\n\n            // Smoothing horizon of the r_S EMA, in units of tau_target.\n""",
        """            if (env_float(\"OU_III_ADAPT_TAU_SEC\", v)) {\n                filter.setAdaptationTimeConstants(v);\n            }\n\n            // Common tau/sigma EMA horizon in measured sea periods.\n            if (env_float(\"OU_ADAPT_TAU_PERIODS\", v)) {\n                filter.setAdaptationPeriods(v);\n            }\n            if (env_float(\"OU_III_ADAPT_TAU_PERIODS\", v)) {\n                filter.setAdaptationPeriods(v);\n            }\n\n            // Frequency EMA inside the tuner, also in measured sea periods.\n            if (env_float(\"OU_TUNER_FREQ_PERIODS\", v)) {\n                filter.setTunerFreqSmoothingPeriods(v);\n            }\n            if (env_float(\"OU_III_TUNER_FREQ_PERIODS\", v)) {\n                filter.setTunerFreqSmoothingPeriods(v);\n            }\n            if (env_float(\"OU_TUNER_FREQ_TAU_SEC\", v)) {\n                filter.setTunerFreqSmoothingTimeConstant(v);\n            }\n            if (env_float(\"OU_III_TUNER_FREQ_TAU_SEC\", v)) {\n                filter.setTunerFreqSmoothingTimeConstant(v);\n            }\n\n            // Smoothing horizon of the r_S EMA, in units of tau_target.\n""",
        "sim env knobs",
    )
    SIM.write_text(text)


def patch_sources() -> None:
    patch_tuner()
    patch_ou3()
    patch_sim()
    print("PATCHED sea-period EMA experiment sources", flush=True)


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def point_args(adapt_values: list[float], freq_values: list[float]) -> list[str]:
    args: list[str] = []
    for a, f in itertools.product(adapt_values, freq_values):
        label = f"a={a:g},f={f:g}"
        spec = (f"{label}:OU_III_ADAPT_TAU_PERIODS={a:g}:"
                f"OU_III_TUNER_FREQ_PERIODS={f:g}")
        args += ["--point", spec]
    return args


def run_coarse() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    # Covers horizons from roughly 0.2 s to 5 s over the calibrated T_z range,
    # bracketing the legacy 1.0 s/1.8 s constants on both short and long seas.
    adapt = [0.05, 0.10, 0.20, 0.35, 0.60]
    freq = [0.05, 0.10, 0.20, 0.35, 0.60]
    cmd = [
        sys.executable, "tools/ou_ema_adapt_study.py", "transition",
        "--family", "OU_III",
        "--transition-stage", "confirm",
        "--transition-pair", "large,small",
        "--transition-dir", "up,down",
        "--transition-seeds", "11,12",
        "--seeds", "default",
        "--jobs", str(max(1, min(4, os.cpu_count() or 2))),
        "--out", str(OUTDIR / "coarse.csv"),
    ] + point_args(adapt, freq)
    run(cmd)


def print_ranking() -> None:
    path = OUTDIR / "coarse_summary.csv"
    if not path.exists():
        raise SystemExit(f"missing {path}")
    rows = list(csv.DictReader(path.open(newline="")))
    by_label: dict[str, dict[str, dict[str, str]]] = {}
    for r in rows:
        by_label.setdefault(r["label"], {})[r["segment"]] = r

    ranked = []
    for label, segs in by_label.items():
        if label == "baseline":
            continue
        if not all(k in segs for k in ("blend", "recover", "window")):
            continue
        b, rec, win = segs["blend"], segs["recover"], segs["window"]
        regressions = sum(int(float(x.get("gate_regressions", "0") or 0))
                          for x in (b, rec, win))
        objective = (0.45 * float(b["score"]) +
                     0.35 * float(rec["score"]) +
                     0.20 * float(win["score"]) +
                     100.0 * regressions)
        ranked.append((objective, label, b, rec, win))
    ranked.sort()

    print("SEA_EMA_COARSE_RANKING", flush=True)
    for objective, label, b, rec, win in ranked[:10]:
        print(
            f"RANK label={label} objective={objective:.6f} "
            f"blend_3d={float(b['rms_3d_gm']):.6f} blend_z={float(b['z_rms_gm']):.6f} "
            f"blend_3d_max={float(b['rms_3d_max']):.6f} "
            f"recover_3d={float(rec['rms_3d_gm']):.6f} recover_z={float(rec['z_rms_gm']):.6f} "
            f"window_3d={float(win['rms_3d_gm']):.6f} window_z={float(win['z_rms_gm']):.6f} "
            f"regressions={int(float(b['gate_regressions'])) + int(float(rec['gate_regressions'])) + int(float(win['gate_regressions']))}",
            flush=True,
        )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=["patch", "coarse", "rank"])
    ns = ap.parse_args()
    if ns.command == "patch":
        patch_sources()
    elif ns.command == "coarse":
        run_coarse()
    else:
        print_ranking()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
