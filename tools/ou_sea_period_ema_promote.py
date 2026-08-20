#!/usr/bin/env python3
"""Promote the measured OU-III sea-period EMA law into production sources.

This is a one-time branch helper.  It starts from current main sources, applies
the already validated experiment transformation, restores the deployed 0.1 s
parameter-activation cadence, then expresses the measured winner in the paper's
sea-time convention

    T_sea = T_z / 2,
    tau_adapt = 0.40 T_sea,
    tau_freq  = 0.10 T_sea.

Those are exactly the 0.20 T_z / 0.05 T_z settings measured in the CI sweep.
No existing sea-scaled coefficient (c_tau, c_sigma, sigma K-periods, r_S EMA,
pseudo cadence, or regularizer law) is changed.
"""

from __future__ import annotations

import runpy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OU3 = ROOT / "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
TUNER = ROOT / "src/tuner/SeaStateAutoTuner.h"
SIM = ROOT / "tests/kalman_ou_iii/kalman_ou_iii-sim.cpp"
TEST = ROOT / "tests/kalman_ou_iii/tuner_schedule-test.cpp"
DOC = ROOT / "docs/ou-ema-adaptation-tuning.md"
PAPER = ROOT / "doc/kalman_ou_iii/w3d-fus-methods.tex-part"


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{label}: expected one anchor, found {n}")
    return text.replace(old, new, 1)


# Reuse the exact source transformation that was compiled and swept, then put
# back the deployed activation cadence so the production change is only EMA
# memory scaling.
exp = runpy.run_path(str(ROOT / "tools/ou_sea_period_ema_experiment.py"))
exp["patch_sources"]()
runpy.run_path(str(ROOT / "tools/ou_sea_period_ema_preserve_cadence.py"))

# ---------------------------------------------------------------------------
# Generic tuner: retain fixed-seconds default for OU-II/TFG; expose an opt-in
# sea-time mode.  The experiment used T_z directly, so convert to T_sea=T_z/2.
# ---------------------------------------------------------------------------
text = TUNER.read_text()
text = text.replace("setFrequencySmoothingPeriods", "setFrequencySmoothingSeaPeriods")
text = text.replace("getFrequencySmoothingPeriods", "getFrequencySmoothingSeaPeriods")
text = text.replace("freq_periods", "freq_sea_periods")
text = text.replace(
    "Frequency-EMA horizon in periods of the externally supplied wave-band\n"
    "    // frequency.  The input is measurement-only in the OU-III wrapper.\n",
    "Frequency-EMA horizon in sea-time units T_sea=T_z/2 of the externally\n"
    "    // supplied wave-band frequency. The input is measurement-only in OU-III.\n",
)
text = text.replace(
    "float freq_sea_periods = 0.0f;         // >0: frequency EMA horizon / T_eff",
    "float freq_sea_periods = 0.0f;         // >0: frequency EMA horizon / T_sea",
)
text = replace_once(
    text,
    "horizon = freq_sea_periods / f;",
    "horizon = freq_sea_periods * (0.5f / f);",
    "tuner T_sea conversion",
)
TUNER.write_text(text)

# ---------------------------------------------------------------------------
# OU-III: measured defaults in the same T_sea convention used by c_tau.
# ---------------------------------------------------------------------------
text = OU3.read_text()
text = text.replace("ADAPT_TAU_PERIODS", "ADAPT_TAU_SEA_PERIODS")
text = text.replace("TUNER_FREQ_EMA_PERIODS", "TUNER_FREQ_EMA_SEA_PERIODS")
text = text.replace("adapt_tau_periods_", "adapt_tau_sea_periods_")
text = text.replace("tuner_freq_ema_periods_", "tuner_freq_ema_sea_periods_")
text = text.replace("setAdaptationPeriods", "setAdaptationSeaPeriods")
text = text.replace("setTunerFreqSmoothingPeriods", "setTunerFreqSmoothingSeaPeriods")
text = text.replace("setFrequencySmoothingPeriods", "setFrequencySmoothingSeaPeriods")
text = replace_once(
    text,
    "constexpr float ADAPT_TAU_SEA_PERIODS          = 0.0f;\n"
    "constexpr float TUNER_FREQ_EMA_SEA_PERIODS     = 0.0f;",
    "constexpr float ADAPT_TAU_SEA_PERIODS          = 0.40f;\n"
    "constexpr float TUNER_FREQ_EMA_SEA_PERIODS     = 0.10f;",
    "OU-III measured defaults",
)
text = text.replace(
    "A positive period multiplier makes the common tau/sigma EMA self-similar:\n"
    "//     tau_ema = ADAPT_TAU_SEA_PERIODS * T_sea.\n"
    "// Zero selects the legacy ADAPT_TAU_SEC path.  The experiment overrides these\n"
    "// at runtime; the final measured value is committed only after confirmation.\n",
    "A positive sea-period multiplier makes the common tau/sigma EMA self-similar:\n"
    "//     tau_ema = ADAPT_TAU_SEA_PERIODS * T_sea,  T_sea = T_z/2.\n"
    "// The measured default 0.40 is the 0.20*T_z winner from the paired sweep.\n"
    "// Zero selects the legacy ADAPT_TAU_SEC path for controlled ablations.\n",
)
text = text.replace(
    "Common tau/sigma EMA horizon as a fraction/multiple of measured T_sea.",
    "Common tau/sigma EMA horizon as a fraction/multiple of measured T_sea=T_z/2.",
)
text = text.replace(
    "Frequency EMA inside SeaStateAutoTuner, likewise expressed in T_sea.",
    "Frequency EMA inside SeaStateAutoTuner, likewise expressed in T_sea=T_z/2.",
)
text = replace_once(
    text,
    "        freq_stillness_.setTargetFreqHz(min_freq_hz_);\n"
    "        startup_stage_   = StartupStage::Cold;",
    "        freq_stillness_.setTargetFreqHz(min_freq_hz_);\n"
    "        // OU-III opts into the measured sea-scaled tuner-frequency EMA.\n"
    "        // SeaStateAutoTuner itself remains fixed-seconds by default so OU-II/TFG\n"
    "        // are unchanged by this OU-III-only calibration.\n"
    "        tuner_.setFrequencySmoothingSeaPeriods(TUNER_FREQ_EMA_SEA_PERIODS);\n"
    "        startup_stage_   = StartupStage::Cold;",
    "OU-III constructor tuner mode",
)
text = replace_once(
    text,
    "        const float sea_period_sec = 1.0f / f_tune;\n"
    "        adapt_mekf(dt, tau_target_, sigma_target_, RS_target_, sea_period_sec);",
    "        const float sea_time_sec = 0.5f / f_tune;  // T_sea = T_z/2\n"
    "        adapt_mekf(dt, tau_target_, sigma_target_, RS_target_, sea_time_sec);",
    "OU-III T_sea target",
)
text = text.replace("float sea_period_sec)", "float sea_time_sec)")
text = text.replace("std::isfinite(sea_period_sec) && sea_period_sec > 0.0f", "std::isfinite(sea_time_sec) && sea_time_sec > 0.0f")
text = text.replace("adapt_tau_sea_periods_ * sea_period_sec", "adapt_tau_sea_periods_ * sea_time_sec")
OU3.write_text(text)

# ---------------------------------------------------------------------------
# Simulator environment names: explicit SEA_PERIODS avoids the T_z/T_sea
# ambiguity that the temporary experiment names had.
# ---------------------------------------------------------------------------
text = SIM.read_text()
text = text.replace("OU_ADAPT_TAU_PERIODS", "OU_ADAPT_TAU_SEA_PERIODS")
text = text.replace("OU_III_ADAPT_TAU_PERIODS", "OU_III_ADAPT_TAU_SEA_PERIODS")
text = text.replace("OU_TUNER_FREQ_PERIODS", "OU_TUNER_FREQ_SEA_PERIODS")
text = text.replace("OU_III_TUNER_FREQ_PERIODS", "OU_III_TUNER_FREQ_SEA_PERIODS")
text = text.replace("setAdaptationPeriods", "setAdaptationSeaPeriods")
text = text.replace("setTunerFreqSmoothingPeriods", "setTunerFreqSmoothingSeaPeriods")
text = text.replace("Common tau/sigma EMA horizon in measured sea periods.", "Common tau/sigma EMA horizon in measured sea-time units T_sea=T_z/2.")
text = text.replace("Frequency EMA inside the tuner, also in measured sea periods.", "Frequency EMA inside the tuner, also in measured sea-time units T_sea=T_z/2.")
SIM.write_text(text)

# ---------------------------------------------------------------------------
# Regression contract: every sample moves EMA state even if the 0.1 s activation
# gate does not fire; the next activation still remains predictable.
# ---------------------------------------------------------------------------
text = TEST.read_text()
text = replace_once(
    text,
    "    const float active_tau_before = f.mekf_->tau_aw;\n",
    "    if (!near(f.adapt_tau_sea_periods_, 0.40f) ||\n"
    "        !near(f.tuner_.getFrequencySmoothingSeaPeriods(), 0.10f)) {\n"
    "        std::cerr << \"FAIL: OU-III sea-scaled EMA defaults are not 0.40/0.10 of T_sea\\n\";\n"
    "        return 1;\n"
    "    }\n\n"
    "    // Physical measurements are never cadence-decimated.  Even inside the\n"
    "    // 0.1 s parameter-activation interval, a valid sample must move the EMA\n"
    "    // candidate while leaving the active MEKF schedule untouched.\n"
    "    {\n"
    "        Filter s(false);\n"
    "        s.time_ = 2.0;\n"
    "        s.last_adapt_time_sec_ = 2.0;\n"
    "        s.online_tune_apply_pending_ = false;\n"
    "        const float before = s.tune_.tau_applied;\n"
    "        constexpr float T_sea = 2.5f;\n"
    "        constexpr float dt_ema = 0.1f;\n"
    "        const float expected_alpha = 1.0f - std::exp(-dt_ema / (0.40f * T_sea));\n"
    "        const float expected = before + expected_alpha * (3.0f - before);\n"
    "        s.adapt_mekf(dt_ema, 3.0f, 0.8f, 1.2f, T_sea);\n"
    "        if (!near(s.tune_.tau_applied, expected, 1e-5f)) {\n"
    "            std::cerr << \"FAIL: sea-scaled EMA did not consume the physical sample\\n\";\n"
    "            return 1;\n"
    "        }\n"
    "        if (s.online_tune_apply_pending_) {\n"
    "            std::cerr << \"FAIL: EMA sample incorrectly bypassed activation cadence\\n\";\n"
    "            return 1;\n"
    "        }\n"
    "    }\n\n"
    "    const float active_tau_before = f.mekf_->tau_aw;\n",
    "test sea EMA contract",
)
text = replace_once(
    text,
    "    f.adapt_mekf(0.1f, 3.0f, 0.8f, 1.2f);",
    "    f.adapt_mekf(0.1f, 3.0f, 0.8f, 1.2f, 2.5f);",
    "test adapt_mekf signature",
)
TEST.write_text(text)

# ---------------------------------------------------------------------------
# Engineering study note: append a concise, reproducible third pass rather than
# rewriting the earlier historical measurements.
# ---------------------------------------------------------------------------
text = DOC.read_text()
marker = "### The other two averaging horizons\n"
insert = r'''### Third pass: make the remaining fixed-second EMAs self-similar

The earlier fixed-seconds sweep above answered only whether `1.8 s` should be
shorter.  It did not test the stronger physical hypothesis that the remaining
adaptation memories should be expressed in the same sea time scale as the OU
operating point.  A joint sweep therefore replaced only the two fixed-second
EMAs while leaving `c_tau`, `c_sigma`, `OU_SIGMA_VAR_K_PERIODS`,
`ADAPT_RS_MULT`, the pseudo-update cadence, and the SpectralMSE law unchanged.
Every physical sample still updates the tuner and the EMA state; the existing
0.1 s timer only gates when the already-smoothed candidate becomes active.

The sweep used the same two transition endpoint pairs in both directions.  A
coarse joint grid was followed by confirmation over 3 wave-phase seeds x 3 IMU
seed triplets (36 transition realizations per candidate) plus the eight
stationary JONSWAP/PM-Stokes seas at the same three IMU seeds (24 realizations).
The measured winner in the paper's convention `T_sea=T_z/2` is

```
tau_adapt = 0.40 T_sea   (= 0.20 T_z)
tau_freq  = 0.10 T_sea   (= 0.05 T_z)
```

Ratios below are to the previous fixed `1.8 s / 1.0 s` baseline:

| ensemble / interval | 3D RMS | Z RMS | roll/pitch | yaw | accel bias |
| --- | ---: | ---: | ---: | ---: | ---: |
| transition blend | 0.99493 | 0.99298 | 0.99863 | 0.99964 | 0.99933 |
| transition run-on | 0.99778 | 0.99658 | 0.99862 | 0.99852 | 0.99905 |
| transition whole 900 s | 0.99866 | 0.99659 | 0.99723 | 0.99951 | 0.99735 |
| stationary 900 s | 0.99913 | 0.99614 | 0.99937 | 0.99960 | 0.99941 |

There were no new quality-gate regressions.  Paired bootstrap intervals on the
primary endpoints exclude unity: stationary 3D RMS is 0.99913
[0.99898, 0.99928] and stationary Z RMS 0.99614 [0.99512, 0.99714]; transition
whole-window 3D RMS is 0.99866 [0.99843, 0.99890] and Z RMS 0.99659
[0.99586, 0.99731].

A deliberately aggressive fixed-seconds control (`0.5 s / 0.25 s`) wins some
transition 3D averages, but has a worse stationary vertical tail (maximum Z
ratio 1.00133 versus 0.99984 for the sea-scaled winner).  A fixed control
matched near the nominal 5 s zero-crossing period (`1.0 s / 0.25 s`) is also
slightly better in mean 3D RMS, but the sea-scaled law is better in Z RMS and in
the composite/tail score on the transition and stationary ensembles.  The
chosen law is therefore not merely "make the EMA faster"; it buys the expected
cross-sea normalization while preserving the stationary tail.

Reproduce the comparison with `tools/ou_ema_adapt_study.py` using the simulator
environment knobs `OU_III_ADAPT_TAU_SEA_PERIODS` and
`OU_III_TUNER_FREQ_SEA_PERIODS`.  Fixed-seconds controls remain available as
`OU_III_ADAPT_TAU_SEC` and `OU_III_TUNER_FREQ_TAU_SEC`.

'''
if marker not in text:
    raise RuntimeError("docs insertion marker missing")
text = text.replace(marker, insert + marker, 1)
DOC.write_text(text)

# ---------------------------------------------------------------------------
# Manuscript implementation text: same normalized sea-time convention.
# ---------------------------------------------------------------------------
text = PAPER.read_text()
text = replace_once(
    text,
    "Once $\\widehat T_z$ is valid, a fixed-horizon EMA forms\n"
    "\\begin{equation}\n"
    "f_{\\mathrm{tune}}\n"
    "=\\clip\\!\\left(\\overline{1/\\widehat{T}_z},\n"
    " f_{\\mathrm{tune,min}},f_{\\mathrm{tune,max}}\\right).\n"
    "\\label{eq:tuning-frequency}\n"
    "\\end{equation}\n"
    "Before then, a bounded wave-band prior is used.  The applied smoothing time is\n"
    "\\SI{1}{s}.  This frequency sets the OU time scale, acceleration-statistic band\n",
    "Once $\\widehat T_z$ is valid, the reciprocal period is smoothed on a\n"
    "sea-scaled horizon,\n"
    "\\begin{equation}\n"
    "\\tau_f=c_f T_{\\mathrm{sea}},\\qquad\n"
    "T_{\\mathrm{sea}}=\\widehat T_z/2,\\qquad\n"
    "\\alpha_f=1-e^{-\\Delta t/\\tau_f},\n"
    "\\label{eq:tuning-frequency-horizon}\n"
    "\\end{equation}\n"
    "and $f_{\\mathrm{tune}}$ is the bounded EMA of $1/\\widehat T_z$.  The\n"
    "measured value is $c_f=0.10$.  Before a valid period is available, a bounded\n"
    "wave-band prior is used.  This frequency sets the OU time scale, acceleration-statistic band\n",
    "paper frequency EMA",
)
text = replace_once(
    text,
    "with $\\tau_{\\mathrm{adapt}}=\\SI{1.8}{s}$.  The regularizer uses the slower\n",
    "with a sea-scaled horizon\n"
    "\\begin{equation}\n"
    "\\tau_{\\mathrm{adapt}}=c_{\\mathrm{adapt}}T_{\\mathrm{sea}},\\qquad\n"
    "c_{\\mathrm{adapt}}=0.40.\n"
    "\\label{eq:adapt-ema-sea-time}\n"
    "\\end{equation}\n"
    "The regularizer uses the slower\n",
    "paper common EMA",
)
text = replace_once(
    text,
    "where $m_{r_S}=3$.  Candidates are committed every \\SI{0.1}{s} and staged by\n"
    "one IMU sample, so parameters used at step $k$ are measurable with respect to\n"
    "data through $k-1$.\n",
    "where $m_{r_S}=3$.  Every valid physical sample contributes to the tuner and\n"
    "to the EMAs above.  The resulting candidate is activated every \\SI{0.1}{s}\n"
    "and staged by one IMU sample; this activation cadence does not decimate the\n"
    "measurement stream.  Parameters used at step $k$ are therefore measurable\n"
    "with respect to data through $k-1$.\n",
    "paper activation cadence",
)
PAPER.write_text(text)

print("PROMOTED measured OU-III sea-period EMA defaults (0.40, 0.10 of T_sea)")
