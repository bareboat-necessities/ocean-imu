#define EIGEN_NON_ARDUINO

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

#include "tuner/AccelVibrationGuard.h"

using seastate::tuner::AccelVibrationGuard;

namespace {

bool check(bool condition, const char* message) {
    if (!condition) std::cerr << "FAIL: " << message << '\n';
    return condition;
}

constexpr float FS = 200.0f;
constexpr float DT = 1.0f / FS;
constexpr float TWO_PI = 6.28318530717958647692f;
constexpr float G = 9.80665f;

// The deployed per-axis accelerometer white noise, matching what the shared
// simulation harness injects: 1.51e-3 g.
constexpr float ACC_WHITE_SIGMA = 1.51e-3f * G;

// The engine lines the study records at the nominal cruise condition, as
// recorded amplitudes in m/s^2 at their aliased frequencies.
struct Line { float hz; float amp; };
constexpr Line CRUISE_LINES[] = {
    {60.0f, 0.346f}, {80.0f, 0.151f}, {20.0f, 0.057f}, {40.0f, 0.026f},
    {20.0f, 0.080f}, {40.0f, 0.158f}, {80.0f, 0.068f}, {15.4f, 0.182f},
    {46.2f, 0.222f}, {92.3f, 0.056f},
};

// Runs a guard over a synthetic record and reports what it measured and did.
struct Outcome {
    float removed_rms = 0.0f;
    float engagement = 0.0f;
    float wave_amplitude_out = 0.0f;
    bool ever_moved = false;
};

Outcome run(AccelVibrationGuard& guard, bool with_engine, float wave_hz,
            float wave_amp, float seconds)
{
    std::mt19937 rng(4242u);
    std::normal_distribution<float> white(0.0f, ACC_WHITE_SIGMA);

    const int n = static_cast<int>(FS * seconds);
    const int score_from = n / 2;
    Outcome out;
    double sum_sq = 0.0;
    int count = 0;

    for (int k = 0; k < n; ++k) {
        const float t = static_cast<float>(k) * DT;
        // A gravity vector with a wave-band oscillation on the vertical axis.
        Eigen::Vector3f acc(white(rng), white(rng),
                            G + wave_amp * std::sin(TWO_PI * wave_hz * t) + white(rng));
        if (with_engine) {
            for (const auto& line : CRUISE_LINES) {
                acc.z() += line.amp * std::sin(TWO_PI * line.hz * t);
                acc.x() += 0.75f * line.amp * std::sin(TWO_PI * line.hz * t + 1.0f);
                acc.y() += 0.55f * line.amp * std::sin(TWO_PI * line.hz * t + 2.0f);
            }
        }
        const Eigen::Vector3f before = acc;
        const Eigen::Vector3f after = guard.step(acc, DT);
        if (after != before) out.ever_moved = true;

        if (k >= score_from) {
            const float centered = after.z() - G;
            sum_sq += static_cast<double>(centered) * static_cast<double>(centered);
            ++count;
        }
    }
    out.removed_rms = guard.removedRms();
    out.engagement = guard.engagement();
    out.wave_amplitude_out =
        static_cast<float>(std::sqrt(2.0 * sum_sq / static_cast<double>(count)));
    return out;
}

bool test_disabled_guard_is_bit_transparent() {
    AccelVibrationGuard guard;  // no cutoff set
    bool ok = check(!guard.enabled(), "a guard with no cutoff reports disabled");

    std::mt19937 rng(7u);
    std::normal_distribution<float> white(0.0f, 1.0f);
    for (int k = 0; k < 2000; ++k) {
        const Eigen::Vector3f acc(white(rng), white(rng), G + white(rng));
        ok &= check(guard.step(acc, DT) == acc, "disabled guard returns its input");
    }
    ok &= check(guard.groupDelaySec() == 0.0f, "disabled guard adds no delay");
    return ok;
}

bool test_clean_input_never_engages_the_guard() {
    // The whole point of the engagement gate: a quiet installation must keep
    // the unconditioned measurement path, bit for bit.
    AccelVibrationGuard guard;
    guard.setCutoffHz(14.0f);

    const Outcome clean = run(guard, /*with_engine=*/false, 0.25f, 1.5f, 300.0f);
    bool ok = check(clean.removed_rms < AccelVibrationGuard::kEngageLoDefault,
                    "clean out-of-band level sits below the engagement floor");
    ok &= check(clean.engagement == 0.0f, "clean input leaves the guard dormant");
    ok &= check(!clean.ever_moved, "dormant guard never alters a sample");
    std::cout << "  clean out-of-band RMS = " << clean.removed_rms
              << " m/s^2 (floor " << AccelVibrationGuard::kEngageLoDefault << ")\n";
    return ok;
}

bool test_detector_floor_does_not_follow_the_motion() {
    // The property the whole default rests on: the detector reads the sensor,
    // not the sea, so a big sea never looks like machinery and a near-still one
    // does not lower the bar.  Swept over a 100:1 range of wave amplitude.
    bool ok = true;
    float lowest = 1e9f, highest = 0.0f;
    for (const float amplitude : {0.05f, 0.5f, 1.5f, 5.0f}) {
        AccelVibrationGuard guard;
        guard.setCutoffHz(14.0f);
        const Outcome out = run(guard, /*with_engine=*/false, 0.25f, amplitude, 300.0f);
        ok &= check(out.engagement == 0.0f, "no wave amplitude engages the guard");
        ok &= check(!out.ever_moved, "no wave amplitude alters a sample");
        lowest = std::min(lowest, out.removed_rms);
        highest = std::max(highest, out.removed_rms);
    }
    ok &= check(highest < 1.15f * lowest,
                "detector floor is within 15% across a 100:1 amplitude range");
    std::cout << "  detector floor over 100:1 wave amplitude = " << lowest
              << " to " << highest << " m/s^2\n";
    return ok;
}

bool test_engine_vibration_engages_and_is_removed() {
    AccelVibrationGuard guard;
    guard.setCutoffHz(14.0f);

    const Outcome engine = run(guard, /*with_engine=*/true, 0.25f, 1.5f, 300.0f);
    bool ok = check(engine.removed_rms > AccelVibrationGuard::kEngageHiDefault,
                    "engine out-of-band level clears the engagement ceiling");
    ok &= check(engine.engagement == 1.0f, "engine vibration fully engages the guard");

    // What survives on the vertical axis should be the wave, not the engine.
    // Unguarded, the engine lines alone carry ~0.34 m/s^2 RMS there.
    ok &= check(std::fabs(engine.wave_amplitude_out - 1.5f) < 0.25f,
                "the wave-band amplitude survives the guard");
    std::cout << "  engine out-of-band RMS = " << engine.removed_rms
              << " m/s^2, recovered wave amplitude = "
              << engine.wave_amplitude_out << " m/s^2\n";
    return ok;
}

// Feeds the guard for `seconds` and returns the time in seconds at which the
// predicate first held, or -1 if it never did.
template <typename Predicate>
float time_until(AccelVibrationGuard& guard, bool with_engine, float seconds,
                 Predicate held)
{
    std::mt19937 rng(99u);
    std::normal_distribution<float> white(0.0f, ACC_WHITE_SIGMA);
    const int n = static_cast<int>(FS * seconds);
    for (int k = 0; k < n; ++k) {
        const float t = static_cast<float>(k) * DT;
        Eigen::Vector3f acc(white(rng), white(rng),
                            G + 1.5f * std::sin(TWO_PI * 0.25f * t) + white(rng));
        if (with_engine) {
            for (const auto& line : CRUISE_LINES) {
                acc.z() += line.amp * std::sin(TWO_PI * line.hz * t);
                acc.x() += 0.75f * line.amp * std::sin(TWO_PI * line.hz * t + 1.0f);
                acc.y() += 0.55f * line.amp * std::sin(TWO_PI * line.hz * t + 2.0f);
            }
        }
        guard.step(acc, DT);
        if (held(guard)) return static_cast<float>(k) * DT;
    }
    return -1.0f;
}

bool test_guard_releases_when_the_engine_stops() {
    // Engaging is only half of it: an estimator that keeps paying the guard's
    // group delay, and keeps telling the MEKF to distrust a quiet
    // accelerometer, long after the engine is shut down would be worse than
    // one that never engaged.  Both signals have to come back down.
    AccelVibrationGuard guard;
    guard.setCutoffHz(14.0f);

    time_until(guard, /*with_engine=*/true, 300.0f,
               [](const AccelVibrationGuard&) { return false; });
    bool ok = check(guard.engagement() == 1.0f, "fully engaged before shutdown");
    ok &= check(guard.excessRms() > 0.0f, "excess is driving before shutdown");

    AccelVibrationGuard excess_guard = guard;
    const float excess_release = time_until(
        excess_guard, /*with_engine=*/false, 600.0f,
        [](const AccelVibrationGuard& g) { return g.excessRms() == 0.0f; });

    AccelVibrationGuard engage_guard = guard;
    const float engage_release = time_until(
        engage_guard, /*with_engine=*/false, 600.0f,
        [](const AccelVibrationGuard& g) { return g.engagement() == 0.0f; });

    ok &= check(excess_release >= 0.0f, "the covariance drive returns to zero");
    ok &= check(engage_release >= 0.0f, "the guard fully disengages");
    // Both must release in well under the time a passage spends under power,
    // and the covariance must let go no later than the conditioning does --
    // releasing R first is the safe order, since it stops understating the
    // measurement while the low pass is still settling.
    ok &= check(excess_release < 60.0f, "covariance drive releases inside a minute");
    ok &= check(engage_release < 180.0f, "guard disengages inside three minutes");
    ok &= check(excess_release <= engage_release,
                "the covariance releases no later than the conditioning");
    std::cout << "  release after shutdown: covariance drive " << excess_release
              << " s, conditioning " << engage_release << " s\n";
    return ok;
}

bool test_separation_between_clean_and_engine_is_wide() {
    AccelVibrationGuard clean_guard;
    clean_guard.setCutoffHz(14.0f);
    const float clean = run(clean_guard, false, 0.25f, 1.5f, 300.0f).removed_rms;

    AccelVibrationGuard engine_guard;
    engine_guard.setCutoffHz(14.0f);
    const float engine = run(engine_guard, true, 0.25f, 1.5f, 300.0f).removed_rms;

    // The thresholds are only defensible if the two populations are far apart.
    const bool ok = check(engine > 5.0f * clean,
                          "engine and clean out-of-band levels differ by >5x");
    std::cout << "  separation = " << (engine / std::max(1e-6f, clean)) << "x\n";
    return ok;
}

bool test_group_delay_follows_the_configuration() {
    AccelVibrationGuard guard;
    guard.setCutoffHz(14.0f);
    guard.setPoles(2);
    guard.setEngagement(0.0f, 0.0f, 0.1f);  // lo <= 0 engages unconditionally

    run(guard, /*with_engine=*/false, 0.25f, 1.5f, 60.0f);
    const float expected = 2.0f / (TWO_PI * 14.0f);
    bool ok = check(std::fabs(guard.groupDelaySec() - expected) < 1e-4f,
                    "fully engaged group delay is poles / (2 pi fc)");
    ok &= check(guard.engagement() == 1.0f,
                "a zero lower threshold engages the guard unconditionally");
    return ok;
}

bool test_poles_are_clamped_to_the_supported_range() {
    AccelVibrationGuard guard;
    guard.setPoles(0);
    bool ok = check(guard.poles() == 1, "pole count clamps up to one");
    guard.setPoles(99);
    ok &= check(guard.poles() == AccelVibrationGuard::kMaxPoles,
                "pole count clamps down to the maximum");
    return ok;
}

bool test_reset_restores_the_dormant_state() {
    AccelVibrationGuard guard;
    guard.setCutoffHz(14.0f);
    run(guard, /*with_engine=*/true, 0.25f, 1.5f, 120.0f);
    bool ok = check(guard.engagement() > 0.0f, "guard engaged before reset");

    guard.reset();
    ok &= check(guard.engagement() == 0.0f, "reset clears engagement");
    ok &= check(guard.removedRms() == 0.0f, "reset clears the measurement");

    const Eigen::Vector3f acc(0.0f, 0.0f, G);
    ok &= check(guard.step(acc, DT) == acc, "first sample after reset is passed through");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_disabled_guard_is_bit_transparent();
    ok &= test_clean_input_never_engages_the_guard();
    ok &= test_detector_floor_does_not_follow_the_motion();
    ok &= test_engine_vibration_engages_and_is_removed();
    ok &= test_guard_releases_when_the_engine_stops();
    ok &= test_separation_between_clean_and_engine_is_wide();
    ok &= test_group_delay_follows_the_configuration();
    ok &= test_poles_are_clamped_to_the_supported_range();
    ok &= test_reset_restores_the_dormant_state();

    if (!ok) {
        std::cerr << "accel_vibration_guard-test FAILED\n";
        return 1;
    }
    std::cout << "accel_vibration_guard-test OK\n";
    return 0;
}
