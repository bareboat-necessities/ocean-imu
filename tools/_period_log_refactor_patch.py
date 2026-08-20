#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def replace_once(rel, old, new):
    path = ROOT / rel
    text = path.read_text(encoding="utf-8")
    if new in text:
        return False
    if old not in text:
        raise SystemExit(f"pattern not found in {rel}: {old[:120]!r}")
    text = text.replace(old, new, 1)
    path.write_text(text, encoding="utf-8")
    return True


def replace_all(rel, old, new):
    path = ROOT / rel
    text = path.read_text(encoding="utf-8")
    if old not in text:
        return False
    text = text.replace(old, new)
    path.write_text(text, encoding="utf-8")
    return True


changed = False

# Expose the two statistically meaningful wave-period horizons for controlled
# sweeps and preserve their configuration across tracking resets.
for rel in (
    "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
    "src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h",
    "src/kalman_tfg/SeaStateFusionFilter_TFG.h",
):
    changed |= replace_once(
        rel,
        "    inline float getWavePeriodSec() const noexcept { return wave_period_.getPeriodSec(); }\n"
        "    inline bool wavePeriodReady() const noexcept { return wave_period_.isReady(); }\n",
        "    inline float getWavePeriodSec() const noexcept { return wave_period_.getPeriodSec(); }\n"
        "    inline float getRawWavePeriodSec() const noexcept { return wave_period_.getRawPeriodSec(); }\n"
        "    inline bool wavePeriodReady() const noexcept { return wave_period_.isReady(); }\n"
        "    void setWavePeriodMomentHorizonPeriods(float periods) {\n"
        "        wave_period_.setMomentHorizonPeriods(periods);\n"
        "    }\n"
        "    void setWavePeriodLogSmoothingPeriods(float periods) {\n"
        "        wave_period_.setLogSmoothingPeriods(periods);\n"
        "    }\n"
        "    float getWavePeriodMomentHorizonPeriods() const noexcept {\n"
        "        return wave_period_.getMomentHorizonPeriods();\n"
        "    }\n"
        "    float getWavePeriodLogSmoothingPeriods() const noexcept {\n"
        "        return wave_period_.getLogSmoothingPeriods();\n"
        "    }\n"
    )
    changed |= replace_all(rel, "wave_period_ = WavePeriodEstimator{};", "wave_period_.reset();")

# Add simulator sweep controls. Generic names are shared by all families;
# family-specific names win when both are present.
for family, rel in (
    ("OU_III", "tests/kalman_ou_iii/kalman_ou_iii-sim.cpp"),
    ("OU_II", "tests/kalman_ou_ii/kalman_ou_ii-sim.cpp"),
):
    marker = (
        "            // Ablate the wave-period estimator's input away from the\n"
        "            // complementary-levelled default."
    )
    insertion = (
        "            // Wave-period statistical horizons.  The moment horizon is\n"
        "            // measured in T_z units; the second is the canonical log(T_z)\n"
        "            // smoothing horizon, also in T_z units.\n"
        "            if (env_float(\"OU_WAVE_PERIOD_MOMENT_K_PERIODS\", v)) {\n"
        "                filter.setWavePeriodMomentHorizonPeriods(v);\n"
        "            }\n"
        f"            if (env_float(\"{family}_WAVE_PERIOD_MOMENT_K_PERIODS\", v)) {{\n"
        "                filter.setWavePeriodMomentHorizonPeriods(v);\n"
        "            }\n"
        "            if (env_float(\"OU_WAVE_PERIOD_LOG_SMOOTH_PERIODS\", v)) {\n"
        "                filter.setWavePeriodLogSmoothingPeriods(v);\n"
        "            }\n"
        f"            if (env_float(\"{family}_WAVE_PERIOD_LOG_SMOOTH_PERIODS\", v)) {{\n"
        "                filter.setWavePeriodLogSmoothingPeriods(v);\n"
        "            }\n\n"
        + marker
    )
    changed |= replace_once(rel, marker, insertion)

# TFG uses the same shared estimator but a slightly different simulator block;
# add generic controls before its wave-period-input ablation if present.
tfg = ROOT / "tests/kalman_tfg/kalman_tfg-sim.cpp"
if tfg.exists():
    text = tfg.read_text(encoding="utf-8")
    marker = "            if (const char* src = std::getenv(\"W3D_WAVE_PERIOD_INPUT\")) {"
    if marker in text and "OU_WAVE_PERIOD_MOMENT_K_PERIODS" not in text:
        ins = (
            "            if (env_float(\"OU_WAVE_PERIOD_MOMENT_K_PERIODS\", v)) {\n"
            "                filter.setWavePeriodMomentHorizonPeriods(v);\n"
            "            }\n"
            "            if (env_float(\"TFG_WAVE_PERIOD_MOMENT_K_PERIODS\", v)) {\n"
            "                filter.setWavePeriodMomentHorizonPeriods(v);\n"
            "            }\n"
            "            if (env_float(\"OU_WAVE_PERIOD_LOG_SMOOTH_PERIODS\", v)) {\n"
            "                filter.setWavePeriodLogSmoothingPeriods(v);\n"
            "            }\n"
            "            if (env_float(\"TFG_WAVE_PERIOD_LOG_SMOOTH_PERIODS\", v)) {\n"
            "                filter.setWavePeriodLogSmoothingPeriods(v);\n"
            "            }\n\n"
        )
        tfg.write_text(text.replace(marker, ins + marker, 1), encoding="utf-8")
        changed = True

# The OU-III simulator comment still called the cubic arm deployed even after
# SpectralMSE became the default.
changed |= replace_all(
    "tests/kalman_ou_iii/kalman_ou_iii-sim.cpp",
    "            // 0 = Cubic (deployed), 1 = StrongRiccati, 2 = PosteriorRiccati,\n"
    "            // 3 = SpectralMSE (bias-variance).",
    "            // 0 = Cubic (legacy fallback), 1 = StrongRiccati, 2 = PosteriorRiccati,\n"
    "            // 3 = SpectralMSE (deployed default).",
)

# The schedule tests should no longer assert the removed second frequency EMA.
changed |= replace_once(
    "tests/kalman_ou_iii/tuner_schedule-test.cpp",
    "    if (!near(f.adapt_tau_sea_periods_, 0.40f) ||\n"
    "        !near(f.tuner_.getFrequencySmoothingSeaPeriods(), 0.10f)) {\n"
    "        std::cerr << \"FAIL: OU-III sea-scaled EMA defaults are not 0.40/0.10 of T_sea\\n\";\n"
    "        return 1;\n"
    "    }\n",
    "    if (!near(f.adapt_tau_sea_periods_, 0.40f)) {\n"
    "        std::cerr << \"FAIL: OU-III parameter-slew default is not 0.40 of T_sea\\n\";\n"
    "        return 1;\n"
    "    }\n"
    "    SeaStateAutoTuner direct_freq(2.0f);\n"
    "    direct_freq.setFrequencySmoothingSeaPeriods(0.10f); // compatibility no-op\n"
    "    direct_freq.update(0.005f, 0.0f, 0.25f);\n"
    "    if (!near(direct_freq.getFrequencyHz(), 0.25f) ||\n"
    "        !near(direct_freq.getFrequencySmoothingHorizonSec(), 0.0f)) {\n"
    "        std::cerr << \"FAIL: tuner reintroduced a second wave-frequency smoother\\n\";\n"
    "        return 1;\n"
    "    }\n",
)

changed |= replace_once(
    "tests/kalman_ou_ii/tuner_schedule-test.cpp",
    "    if (!near(f.getAdaptationSeaPeriods(), 0.40f) ||\n"
    "        !near(f.getTunerFreqSmoothingSeaPeriods(), 0.10f)) {\n"
    "        std::cerr << \"FAIL: OU-II sea-scaled EMA defaults are not 0.40/0.10 of T_sea\\n\";\n"
    "        return 1;\n"
    "    }\n",
    "    if (!near(f.getAdaptationSeaPeriods(), 0.40f)) {\n"
    "        std::cerr << \"FAIL: OU-II parameter-slew default is not 0.40 of T_sea\\n\";\n"
    "        return 1;\n"
    "    }\n",
)

# Rewrite the universal-tuner guard check to assert direct canonical frequency
# consumption rather than the removed frequency EMA horizon.
ou2 = ROOT / "tests/kalman_ou_ii/tuner_schedule-test.cpp"
text = ou2.read_text(encoding="utf-8")
old = '''    // The shared tuner must use the same guards.  At absurdly short/long input
    // periods the deployed 0.10*T_sea frequency EMA becomes 0.05/0.60 s, and
    // the K=2 variance horizon becomes 2/24 s.  These values are far outside
    // the eight calibrated seas but finite and deliberately bounded.
    {
        SeaStateAutoTuner fast(2.0f, 1.0f);
        fast.setFrequencySmoothingSeaPeriods(0.10f);
        fast.update(0.005f, 0.0f, 100.0f);
        if (!near(fast.getFrequencySmoothingHorizonSec(), 0.05f) ||
            !near(fast.getVarianceHorizonSec(), 2.0f)) {
            std::cerr << "FAIL: short-period tuner horizons escaped universal clamps\\n";
            return 1;
        }

        SeaStateAutoTuner slow(2.0f, 1.0f);
        slow.setFrequencySmoothingSeaPeriods(0.10f);
        slow.update(0.005f, 0.0f, 0.001f);
        if (!near(slow.getFrequencySmoothingHorizonSec(), 0.60f) ||
            !near(slow.getVarianceHorizonSec(), 24.0f)) {
            std::cerr << "FAIL: long-period tuner horizons escaped universal clamps\\n";
            return 1;
        }
    }
'''
new = '''    // The shared tuner consumes the canonical wave frequency directly.  Only
    // the variance horizon is dynamically bounded here: K=2 gives 2 s at the
    // short-time guard and 24 s at the long-time guard.
    {
        SeaStateAutoTuner fast(2.0f);
        fast.update(0.005f, 0.0f, 100.0f);
        if (!near(fast.getFrequencyHz(), 5.0f) ||
            !near(fast.getVarianceHorizonSec(), 2.0f)) {
            std::cerr << "FAIL: short-period tuner guard/direct-frequency contract failed\\n";
            return 1;
        }

        SeaStateAutoTuner slow(2.0f);
        slow.update(0.005f, 0.0f, 0.001f);
        if (!near(slow.getFrequencyHz(), 0.05f) ||
            !near(slow.getVarianceHorizonSec(), 24.0f)) {
            std::cerr << "FAIL: long-period tuner guard/direct-frequency contract failed\\n";
            return 1;
        }
    }
'''
if old in text:
    ou2.write_text(text.replace(old, new, 1), encoding="utf-8")
    changed = True
elif new not in text:
    raise SystemExit("OU-II shared tuner guard block not found")

print("period/log refactor patch changed files" if changed else "period/log refactor patch already applied")
