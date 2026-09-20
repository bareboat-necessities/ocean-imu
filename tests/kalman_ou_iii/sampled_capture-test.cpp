// Regression for the indistinguishable physical histories proved in
// docs/ou3-sampled-capture-obstruction.md. This finite source replay checks
// the literal wrapper path; the analytic certificate supplies all-time bounds.
#define EIGEN_NON_ARDUINO
#include <algorithm>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

const float g_std = 9.80665f;
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;

static int fail(const char* message) {
    std::cerr << "FAIL: " << message << '\n';
    return 1;
}

int main() {
    Fusion::Config cfg;
    cfg.sigma_a.setConstant(0.12f);
    cfg.sigma_g.setConstant(0.00135f);
    cfg.sigma_m.setConstant(0.8f);
    cfg.mag_delay_sec = 0.0f;
    cfg.mag_init_min_mag_norm = 5.0f;
    Fusion filter;
    filter.begin(cfg);
    int live = -1, refined = -1, active = -1, accepted = 0;
    for (int k = 1; k <= 45000; ++k) {
        // Both nonzero continuous motions have these exact samples at k/200.
        filter.update(0.005f, Eigen::Vector3f::Zero(),
                      Eigen::Vector3f(0.0f, 0.0f, -g_std));
        if (k % 8 == 0) {
            filter.updateMag(Eigen::Vector3f(75.0f, 0.0f, 0.0f));
            if (filter.raw().mekf().lastMagDiag().accepted) ++accepted;
            if (active >= 0 && !filter.raw().mekf().lastMagDiag().accepted)
                return fail("stationary A21 magnetic correction was not applied");
        }
        if (live < 0 && filter.isLive()) live = k;
        if (refined < 0 && filter.hasRefinedMagReference()) refined = k;
        if (active < 0 && filter.raw().mekf().acc_bias_updates_enabled()) active = k;
        const auto q = filter.raw().mekf().quaternion_boat();
        if (!q.coeffs().allFinite() || q.vec().norm() > 1e-6f)
            return fail("common stationary sample record did not retain level attitude");
        const auto p = filter.raw().mekf().covariance_full();
        // The heading/axial-bias pair is genuinely invariant in this record.
        for (int i : {2, 5}) {
            for (int j = 0; j < 21; ++j) {
                if (j != 2 && j != 5 && std::abs(p(i,j)) > 1e-15f)
                    return fail("heading pair coupled to another axis group");
            }
        }
    }
    if (live < 0 || refined < live || active < refined || accepted < 500)
        return fail("literal capture/refinement/release path did not complete");
    // True down is (0,+/-4/5,3/5) in the body. The estimated down is (0,0,1).
    const double physical_tilt_rad = std::acos(3.0/5.0);
    if (!(physical_tilt_rad > 0.10471975511965977))
        return fail("constructed truth unexpectedly lies in the declared tilt region");
    std::cout << "sampled capture obstruction source replay PASS: live=" << live
              << " refined=" << refined << " active=" << active
              << " accepted=" << accepted << " true_tilt_rad="
              << physical_tilt_rad << '\n';
    return 0;
}
