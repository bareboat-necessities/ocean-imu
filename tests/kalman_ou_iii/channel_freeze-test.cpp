// Guards the OU-III adaptation-channel ablation.
//
// The deployed tuning law couples the integral pseudo-measurement scale to the
// OU parameters through r_S = clip(1.2 sigma_aw tau^3).  Because of that
// coupling, an ablation that freezes tau and sigma_aw with setFixedTuning()
// also freezes r_S, and cannot say which channel carries the benefit.
// setChannelFreeze() exists to break exactly that coupling, so the properties
// checked here are the ones the ablation's validity rests on.
#define EIGEN_NON_ARDUINO

#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

#include <cmath>
#include <iostream>

namespace {

using Filter = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;

bool check(bool condition, const char* message) {
    if (!condition) std::cerr << "FAIL: " << message << '\n';
    return condition;
}

constexpr float NOMINAL_TAU = 1.56533504f;
constexpr float NOMINAL_SIGMA = 1.26339269f;
constexpr float NOMINAL_RS = 5.81489031f;

bool test_rejects_ambiguous_requests() {
    Filter filter;
    bool ok = true;

    // Freezing both channels is setFixedTuning(); freezing neither is the
    // deployed filter.  Neither may be reachable through this entry point, or
    // the four ablation modes stop being distinct.
    ok &= check(
        !filter.setChannelFreeze(true, NOMINAL_TAU, NOMINAL_SIGMA, true, NOMINAL_RS),
        "freezing both channels must be rejected");
    ok &= check(
        !filter.setChannelFreeze(false, NOMINAL_TAU, NOMINAL_SIGMA, false, NOMINAL_RS),
        "freezing neither channel must be rejected");

    ok &= check(!filter.setChannelFreeze(true, -1.0f, NOMINAL_SIGMA, false, 0.0f),
                "non-positive tau must be rejected");
    ok &= check(!filter.setChannelFreeze(true, NOMINAL_TAU, NAN, false, 0.0f),
                "non-finite sigma_aw must be rejected");
    ok &= check(!filter.setChannelFreeze(false, 0.0f, 0.0f, true, -3.0f),
                "non-positive r_S must be rejected");

    ok &= check(!filter.frozenOUChannel() && !filter.frozenRSChannel(),
                "a rejected request must leave both channels adapting");
    return ok;
}

bool test_freezes_exactly_one_channel() {
    bool ok = true;

    Filter rs_only;
    ok &= check(
        rs_only.setChannelFreeze(true, NOMINAL_TAU, NOMINAL_SIGMA, false, 0.0f),
        "freezing the OU channel alone must be accepted");
    ok &= check(rs_only.frozenOUChannel() && !rs_only.frozenRSChannel(),
                "r_S-only mode must hold the OU channel and free r_S");
    ok &= check(std::abs(rs_only.getTauApplied() - NOMINAL_TAU) <= 1e-6f,
                "frozen tau must be applied immediately");
    ok &= check(std::abs(rs_only.getSigmaApplied() - NOMINAL_SIGMA) <= 1e-6f,
                "frozen sigma_aw must be applied immediately");

    Filter ou_only;
    ok &= check(ou_only.setChannelFreeze(false, 0.0f, 0.0f, true, NOMINAL_RS),
                "freezing the r_S channel alone must be accepted");
    ok &= check(!ou_only.frozenOUChannel() && ou_only.frozenRSChannel(),
                "OU-only mode must hold r_S and free the OU channel");
    ok &= check(std::abs(ou_only.getRSApplied() - NOMINAL_RS) <= 1e-6f,
                "frozen r_S must be applied immediately");

    // The tuner must keep running in both cases.  If setChannelFreeze
    // disabled it the way setFixedTuning does, the free channel would stop
    // tracking the sea and the ablation would compare two frozen points.
    ok &= check(rs_only.tunerEnabled() && ou_only.tunerEnabled(),
                "partial freezes must leave the online tuner running");
    return ok;
}

bool test_fixed_tuning_clears_the_freezes() {
    Filter filter;
    bool ok = true;
    ok &= check(
        filter.setChannelFreeze(true, NOMINAL_TAU, NOMINAL_SIGMA, false, 0.0f),
        "setup freeze must be accepted");
    ok &= check(
        filter.setFixedTuning(NOMINAL_TAU, NOMINAL_SIGMA, NOMINAL_RS),
        "setFixedTuning must be accepted");
    // Otherwise a fixed-mode run that followed a partial-mode call would keep
    // a stale channel freeze and silently score a different configuration.
    ok &= check(!filter.frozenOUChannel() && !filter.frozenRSChannel(),
                "setFixedTuning must clear both channel freezes");
    return ok;
}

}  // namespace

int main() {
    if (!test_rejects_ambiguous_requests()) return 1;
    if (!test_freezes_exactly_one_channel()) return 1;
    if (!test_fixed_tuning_clears_the_freezes()) return 1;
    std::cout << "channel_freeze-test: OK\n";
    return 0;
}
