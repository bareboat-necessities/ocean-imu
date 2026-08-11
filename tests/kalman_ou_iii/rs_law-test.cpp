// Pins the integral-regularizer adaptation-law family used by the
// amplitude-exponent ablation in the OU-III paper.
//
// The paper's ablation argument rests on three structural claims about
//     r_S = sqrt(2 r_a) tau^3 / (sqrt(T_S) kappa^3) * (sigma/sigma_ref)^p:
//   1. p = 1 reproduces the deployed Cubic schedule exactly, so the family
//      passes through the shipped filter rather than near it;
//   2. the ratio between family members is independent of tau, so the
//      ablation cannot be confounded by an overall regularizer gain change;
//   3. the PosteriorRiccati transition law reduces to the StrongRiccati
//      asymptote when zeta >> 1 and to the Cubic law when zeta << 1.
// If any of these regress the published comparison stops meaning what it says.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std = 9.80665f;
using Filter = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;

static bool rel_near(float a, float b, float rel) {
    return std::fabs(a - b) <= rel * std::max(1e-12f, std::max(std::fabs(a), std::fabs(b)));
}

static int fail(const char* msg) {
    std::cerr << "FAIL: " << msg << "\n";
    return 1;
}

int main() {
    Filter f(false);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.25f));

    // The deployed default must remain the Cubic schedule: this file adds an
    // ablation surface, it does not change what ships.
    if (f.getRSLaw() != RSAdaptationLaw::Cubic)
        return fail("default OU-III r_S law is no longer Cubic");
    if (!rel_near(f.getRSSigmaExponent(), 0.0f, 1e-6f) &&
        f.getRSSigmaExponent() != 0.0f)
        return fail("default amplitude exponent changed");

    const float r_a = f.getRSAccelNoiseDensity();
    const float kappa = f.getRSPoleKappa();
    const float k3 = kappa * kappa * kappa;
    const float sigma_ref = std::sqrt(2.0f * r_a) /
        (k3 * f.R_S_coeff_ * std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S));

    // kappa is calibrated so sigma_ref is the sigma_aw of the Hs = 1.5 m
    // nominal calibration sea (Table "Frozen operating points").
    if (!rel_near(sigma_ref, 0.724f, 2e-3f)) {
        std::cerr << "  sigma_ref=" << sigma_ref << " expected ~0.724\n";
        return fail("kappa no longer anchors sigma_ref at the nominal sea");
    }

    const float taus[] = {0.6f, 1.25f, 2.179f, 3.59f, 8.0f};
    const float sigmas[] = {0.2f, 0.353f, 0.724f, 1.124f, 2.5f};

    for (float tau : taus) {
        const float TS = f.pseudo_update_period_for_(tau);
        const float cadence = std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S / TS);
        float ratio_at_first_sigma = 0.0f;
        for (int i = 0; i < 5; ++i) {
            const float sg = sigmas[i];

            // (1) p = 1 must reproduce the deployed Cubic filter input, which
            //     is the clamped base times the cadence renormalization.
            f.rs_law_ = RSAdaptationLaw::Cubic;
            const float cubic_eff = f.rs_target_from_law_(tau, sg) * cadence;
            f.rs_law_ = RSAdaptationLaw::StrongRiccati;
            f.rs_sigma_exponent_ = 1.0f;
            const float p1 = f.rs_target_from_law_(tau, sg);
            if (!rel_near(p1, cubic_eff, 1e-4f)) {
                std::cerr << "  tau=" << tau << " sigma=" << sg
                          << " p1=" << p1 << " cubic_eff=" << cubic_eff << "\n";
                return fail("p=1 does not reproduce the deployed Cubic schedule");
            }

            // (2) the p=0 / p=1 ratio must not depend on tau, otherwise the
            //     ablation confounds amplitude shape with overall gain.
            f.rs_sigma_exponent_ = 0.0f;
            const float p0 = f.rs_target_from_law_(tau, sg);
            const float ratio = p0 / cubic_eff;
            if (!rel_near(ratio, sigma_ref / sg, 1e-3f)) {
                std::cerr << "  tau=" << tau << " sigma=" << sg
                          << " ratio=" << ratio << " expected=" << sigma_ref / sg << "\n";
                return fail("strong/cubic ratio is not sigma_ref/sigma");
            }
            if (i == 0) ratio_at_first_sigma = ratio * sigmas[0];
            else if (!rel_near(ratio * sg, ratio_at_first_sigma, 1e-3f))
                return fail("strong/cubic ratio acquired a sigma-dependent factor");
        }
    }

    // (3) transition-law asymptotics.  zeta = 2 sigma^2 tau / r_a.
    f.rs_sigma_exponent_ = 0.0f;
    {
        const float tau = 2.179f;
        // Deep strong branch: must match the StrongRiccati asymptote.
        const float sg_strong = std::sqrt(1e6f * r_a / (2.0f * tau));
        f.rs_law_ = RSAdaptationLaw::PosteriorRiccati;
        const float post_s = f.rs_target_from_law_(tau, sg_strong);
        f.rs_law_ = RSAdaptationLaw::StrongRiccati;
        const float str_s = f.rs_target_from_law_(tau, sg_strong);
        if (!rel_near(post_s, str_s, 2e-3f))
            return fail("PosteriorRiccati does not reduce to StrongRiccati at large zeta");

        // Deep weak branch: q_eff -> 2 sigma^2 tau, so with T_S = c_T tau the
        // pole-preserving schedule becomes r_S ~ sigma * tau^3.  That is the
        // *pre-cadence-change* cubic exponent, not the current filter input,
        // which the self-similar cadence moved to sigma * tau^(5/2).  Assert
        // the exponent rather than equality with the deployed schedule.
        f.rs_law_ = RSAdaptationLaw::PosteriorRiccati;
        float weak_coeff = 0.0f;
        for (float t : {1.5f, 3.0f, 6.0f}) {
            // Hold zeta fixed and deep in the weak branch as tau varies.
            const float sg_weak = std::sqrt(1e-6f * r_a / (2.0f * t));
            const float post_w = f.rs_target_from_law_(t, sg_weak);
            const float coeff = post_w / (sg_weak * t * t * t);
            if (weak_coeff == 0.0f) weak_coeff = coeff;
            else if (!rel_near(coeff, weak_coeff, 5e-3f)) {
                std::cerr << "  tau=" << t << " coeff=" << coeff
                          << " expected=" << weak_coeff << "\n";
                return fail("PosteriorRiccati weak branch is not proportional to sigma*tau^3");
            }
        }
        // ...and that this differs from the deployed tau^(5/2) input by
        // exactly one power of sqrt(tau), which is the whole point of
        // distinguishing the two branches.
        f.rs_law_ = RSAdaptationLaw::Cubic;
        float cubic_coeff = 0.0f;
        for (float t : {1.5f, 3.0f, 6.0f}) {
            const float eff = f.rs_target_from_law_(t, 1.0f) *
                std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S / f.pseudo_update_period_for_(t));
            const float coeff = eff / std::pow(t, 2.5f);
            if (cubic_coeff == 0.0f) cubic_coeff = coeff;
            else if (!rel_near(coeff, cubic_coeff, 5e-3f)) {
                std::cerr << "  tau=" << t << " coeff=" << coeff
                          << " expected=" << cubic_coeff << "\n";
                return fail("deployed Cubic input is not proportional to sigma*tau^(5/2)");
            }
        }
    }

    // The deployed envelope really is the strong branch: the smallest
    // calibrated operating point must still have zeta >> 1.
    f.tune_.tau_applied = 1.247f;
    f.tune_.sigma_applied = 0.353f;
    const float zeta_min = f.getAccelInformationRatio();
    if (!(zeta_min > 1.0e4f)) {
        std::cerr << "  zeta at the Hs=0.27 m point = " << zeta_min << "\n";
        return fail("smallest calibrated sea is no longer deep in the strong branch");
    }

    std::cout << "OU-III r_S law family passed (sigma_ref=" << sigma_ref
              << " m/s^2, kappa=" << kappa
              << ", zeta[Hs=0.27 m]=" << zeta_min << ")\n";
    return 0;
}
