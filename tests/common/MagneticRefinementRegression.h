#pragma once

#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <iostream>

// Physical gravity/field observations share one constant attitude. A field
// rotated by 90 degrees for an entire refinement window must not become an
// absolute heading assignment after the Kalman core has rejected that field.
template<class Wrapper>
bool magneticRefinementRegression() {
    using V = Eigen::Vector3f;
    const Eigen::Quaternionf truth = Eigen::AngleAxisf(.4f, V::UnitZ()) *
        Eigen::AngleAxisf(.1f, V::UnitY()) * Eigen::AngleAxisf(.15f, V::UnitX());
    const V acc = truth.conjugate() * V(0, 0, -9.80665f);
    const V clean = truth.conjugate() * V(20, 0, 43);
    const V shock = truth.conjugate() * V(0, -20, 43);
    typename Wrapper::Config cfg;
    cfg.mag_continuous_hard_iron = false;
    cfg.mag_refine_start_sec = 90;
    cfg.mag_refine_window_sec = 10;
    cfg.online_tune_warmup_sec = 10;
    cfg.proxy_startup_timeout_sec = 60;
    Wrapper f;
    f.begin(cfg);
    float max_error = 0;
    float max_shock_step = 0;
    float initial_error = 0;
    bool rejected = true;
    bool refined_during_shock = false;
    for (int i = 0; i < 36000; ++i) {
        const float t = float(i) * .005f;
        const bool inject = t >= 89 && t < 112;
        const auto before = f.attitudeQuat();
        f.update(.005f, V::Zero(), acc);
        if (i % 5 == 0) {
            f.updateMag(inject ? shock : clean);
            if (inject) rejected = rejected && !f.raw().mekf().lastMagDiag().accepted;
        }
        if (i == 17800) {
            initial_error = before.angularDistance(truth);
            rejected = rejected && f.isLive() &&
                f.raw().mekf().magneticHealth().accepted_count > 0;
        }
        if (inject) max_shock_step = std::max(max_shock_step,
            before.angularDistance(f.attitudeQuat()));
        if (t >= 89) {
            max_error = std::max(max_error, f.attitudeQuat().angularDistance(truth));
            if (t < 112 && f.hasRefinedMagReference()) refined_during_shock = true;
        }
    }
    // The existing coarse north lock need not yet match the eventual fine
    // reference. Bound disturbance-time changes, then require truth recovery.
    const bool ok = rejected && !refined_during_shock && max_shock_step < .01f &&
        f.hasRefinedMagReference() && f.attitudeQuat().angularDistance(truth) < .01f;
    std::cout << "  shared reference-window shock: max error=" << max_error
              << " rad, initial=" << initial_error << ", shock step=" << max_shock_step
              << ", final=" << f.attitudeQuat().angularDistance(truth)
              << ", recovered=" << f.hasRefinedMagReference() << '\n';
    if (!ok) std::cerr << "FAIL: shared reference refinement bypassed magnetic rejection/recovery\n";
    return ok;
}
