// Pins the integral-regularizer adaptation-law family used by the
// amplitude-exponent ablation in the OU-III paper.
//
// The paper's ablation argument rests on four structural claims about
//     r_S = sqrt(2 r_a) tau^3 / (sqrt(T_S) kappa^3) * (sigma/sigma_ref)^p:
//   1. p = 1 reproduces the sigma-tilted schedule c_R sigma_aw tau^3 exactly,
//      so the family passes through that arm rather than near it;
//   2. p = 0 reproduces the deployed Cubic schedule c_R tau^3 up to the
//      constant gain sigma_ref, so the deployed law is the family's
//      sigma-independent member and not a fourth unrelated schedule;
//   3. the ratio between family members is independent of tau, so the
//      ablation cannot be confounded by an overall regularizer gain change;
//   4. the PosteriorRiccati transition law reduces to the StrongRiccati
//      asymptote when zeta >> 1 and to the sigma tau^3 law when zeta << 1.
// If any of these regress the published comparison stops meaning what it says.
//
// The deployed Cubic base carries no sigma_aw factor.  Its acceleration scale
// is the accelerometer noise, r_S,base = C_R sqrt(R_a) tau^3, Eq.
// (adapt-rs-base).  That is the arm the amplitude ablation used to call p = 0,
// so the roles of p = 0 and p = 1 relative to the shipped filter are the
// reverse of what they were before the multiplier was removed.  Two further
// properties of that base are pinned below: it is proportional to sqrt(R_a),
// and at the analytical C_R = sqrt(2h/T_S,0)/kappa^3 the cadence-normalized
// base is identical to the StrongRiccati pole-placement law.
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

    // The deployed default is the bias-variance law.  The three pole-placement
    // laws below it remain the ablation family and are checked as such.
    if (f.getRSLaw() != RSAdaptationLaw::SpectralMSE)
        return fail("default OU-III r_S law is no longer SpectralMSE");

    // The SpectralMSE law must carry the exponents the balance derives:
    //     r_S ~ q_eff^(1/14) sigma_a,B^(6/7) tau^(24/7) / sqrt(T_S),
    // i.e. tau^(41/14) at the filter input away from cadence clamps.  Measure
    // them rather than trusting the expression, and check that dividing
    // c_sigma back out really makes r_S independent of the OU-prior
    // coefficient, which is what lets C_J and c_sigma be tuned separately.
    {
        Filter g(false);
        g.initialize(Eigen::Vector3f::Constant(0.0148f),
                     Eigen::Vector3f::Constant(0.00157f),
                     Eigen::Vector3f::Constant(0.25f));
        g.rs_law_ = RSAdaptationLaw::SpectralMSE;
        const float tau = 2.179f, sg = 0.724f, d = 1.02f;
        const float r0 = g.rs_spectral_mse_target_(tau, sg);
        const float p_a = std::log(g.rs_spectral_mse_target_(tau, sg * d) / r0)
                        / std::log(d);
        if (!rel_near(p_a, 6.0f / 7.0f, 2e-3f)) {
            std::cerr << "  amplitude exponent " << p_a << " expected 6/7\n";
            return fail("SpectralMSE amplitude exponent is not 6/7");
        }
        // tau exponent at the filter input, away from the cadence clamps.
        const float p_t = std::log(g.rs_spectral_mse_target_(tau * d, sg) / r0)
                        / std::log(d);
        if (!rel_near(p_t, 41.0f / 14.0f, 3e-3f)) {
            std::cerr << "  tau exponent " << p_t << " expected 41/14\n";
            return fail("SpectralMSE tau exponent is not 41/14");
        }
        // c_sigma invariance: sigma_aw and sigma_coeff scale together, so the
        // physical sigma_a,B and hence r_S must not move.
        const float before = g.rs_spectral_mse_target_(tau, sg);
        g.setSigmaCoeff(2.0f * g.sigma_coeff_);
        const float after = g.rs_spectral_mse_target_(tau, 2.0f * sg);
        if (!rel_near(after, before, 1e-4f))
            return fail("SpectralMSE r_S is not invariant to c_sigma");
    }

    // Cubic remains a supported deployment configuration, not just an ablation
    // surface: it is the low-cost option for embedded targets.  Pin that it is
    // selectable and that its schedule is still the tau^(5/2) one, so a future
    // change to the default cannot silently strand it.
    {
        Filter mcu(false);
        mcu.initialize(Eigen::Vector3f::Constant(0.0148f),
                       Eigen::Vector3f::Constant(0.00157f),
                       Eigen::Vector3f::Constant(0.25f));
        mcu.setRSLaw(RSAdaptationLaw::Cubic);
        if (mcu.getRSLaw() != RSAdaptationLaw::Cubic)
            return fail("Cubic is no longer selectable as the low-cost law");
        float coeff = 0.0f;
        for (float t : {1.5f, 3.0f, 6.0f}) {
            const float cad = std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S /
                                        mcu.pseudo_update_period_for_(t));
            const float c = mcu.rs_target_from_law_(t, 0.724f) * cad
                          / std::pow(t, 2.5f);
            if (coeff == 0.0f) coeff = c;
            else if (!rel_near(c, coeff, 5e-3f))
                return fail("Cubic is no longer the tau^(5/2) schedule");
        }
    }

    f.rs_law_ = RSAdaptationLaw::Cubic;
    if (!rel_near(f.getRSSigmaExponent(), 0.0f, 1e-6f) &&
        f.getRSSigmaExponent() != 0.0f)
        return fail("default amplitude exponent changed");

    const float r_a = f.getRSAccelNoiseDensity();
    const float kappa = f.getRSPoleKappa();
    const float k3 = kappa * kappa * kappa;
    const float sqrt_Ra = std::sqrt(r_a / FREQ_SMOOTHER_DT);
    const float sigma_ref = std::sqrt(2.0f * r_a) /
        (k3 * f.R_S_coeff_ * sqrt_Ra * std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S));

    // The analytical reference must be what Eq. (adapt-cR-kappa) says it is,
    // C_R = sqrt(2h/T_S,0)/kappa^3, evaluated at the default kappa.
    {
        const float C_R_analytic =
            std::sqrt(2.0f * FREQ_SMOOTHER_DT / PSEUDO_UPDATE_PERIOD_NOMINAL_S) / k3;
        if (!rel_near(C_R_analytic, R_S_COEFF_ANALYTICAL_REFERENCE, 1e-3f)) {
            std::cerr << "  sqrt(2h/T_S0)/kappa^3=" << C_R_analytic
                      << " constant=" << R_S_COEFF_ANALYTICAL_REFERENCE << "\n";
            return fail("R_S_COEFF_ANALYTICAL_REFERENCE is not the pole-placement value");
        }
    }

    // The Cubic base takes its acceleration scale from the sensor: quadrupling
    // the accelerometer noise variance must scale r_S,base by two, with no
    // other change.  This is the property that distinguishes C_R sqrt(R_a)
    // tau^3 from a bare constant times tau^3, which is numerically identical
    // only at the one sensor the coefficient happened to be fitted on.
    {
        const float base_1 = f.rs_target_from_law_(2.179f, 0.724f);
        f.setRSAccelNoiseDensity(4.0f * r_a);
        const float base_4 = f.rs_target_from_law_(2.179f, 0.724f);
        f.setRSAccelNoiseDensity(r_a);
        if (!rel_near(base_4, 2.0f * base_1, 1e-4f)) {
            std::cerr << "  r_a x4 gave " << base_4 << ", expected " << 2.0f * base_1 << "\n";
            return fail("Cubic base is not proportional to sqrt(R_a)");
        }
    }

    // sigma_ref is the amplitude anchor of the tilt family: it is where the
    // p = 0 and p = 1 members cross, and it scales as 1/C_R.  It is no longer
    // pinned to the sigma_aw of a particular sea.  The structural claims below
    // hold for any anchor, so pin it only as positive and finite -- but pin
    // separately that at the analytical C_R it is exactly 1, which is the
    // statement that the cadence-normalized cubic base and the pole-placement
    // law are the same schedule.
    if (!(std::isfinite(sigma_ref) && sigma_ref > 0.0f))
        return fail("tilt anchor sigma_ref is not a positive finite amplitude");
    {
        Filter probe(false);
        probe.initialize(Eigen::Vector3f::Constant(0.0148f),
                         Eigen::Vector3f::Constant(0.00157f),
                         Eigen::Vector3f::Constant(0.25f));
        probe.setRSCoeff(R_S_COEFF_ANALYTICAL_REFERENCE);
        const float ref_anchor = std::sqrt(2.0f * r_a) /
            (k3 * probe.R_S_coeff_ * sqrt_Ra *
             std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S));
        if (!rel_near(ref_anchor, 1.0f, 2e-3f)) {
            std::cerr << "  sigma_ref at the analytical C_R = " << ref_anchor << "\n";
            return fail("analytical C_R does not make the cubic base equal the pole law");
        }
        // ...and say the same thing directly at the filter input: base times
        // the cadence renormalization must equal the StrongRiccati target.
        for (float t : {1.25f, 2.179f, 4.0f}) {
            const float cad = std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S /
                                        probe.pseudo_update_period_for_(t));
            probe.rs_law_ = RSAdaptationLaw::Cubic;
            const float cubic_in = probe.rs_target_from_law_(t, 0.724f) * cad;
            probe.rs_law_ = RSAdaptationLaw::StrongRiccati;
            probe.rs_sigma_exponent_ = 0.0f;
            const float strong_in = probe.rs_target_from_law_(t, 0.724f);
            if (!rel_near(cubic_in, strong_in, 2e-3f)) {
                std::cerr << "  tau=" << t << " cubic=" << cubic_in
                          << " strong=" << strong_in << "\n";
                return fail("analytical C_R does not reproduce StrongRiccati at the input");
            }
        }
    }

    const float taus[] = {0.6f, 1.25f, 2.179f, 3.59f, 8.0f};
    const float sigmas[] = {0.2f, 0.353f, 0.724f, 1.124f, 2.5f};

    for (float tau : taus) {
        const float TS = f.pseudo_update_period_for_(tau);
        const float cadence = std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S / TS);
        float ratio_at_first_sigma = 0.0f;
        for (int i = 0; i < 5; ++i) {
            const float sg = sigmas[i];

            // The deployed Cubic filter input: the base times the cadence
            // renormalization.  It no longer carries sigma_aw.
            f.rs_law_ = RSAdaptationLaw::Cubic;
            const float cubic_eff = f.rs_target_from_law_(tau, sg) * cadence;
            // The sigma-tilted arm the ablation compares against is that same
            // input multiplied by sigma_aw, i.e. the pre-change deployed law.
            const float tilted_eff = cubic_eff * sg;

            // (1) p = 1 must reproduce the sigma-tilted schedule exactly.
            f.rs_law_ = RSAdaptationLaw::StrongRiccati;
            f.rs_sigma_exponent_ = 1.0f;
            const float p1 = f.rs_target_from_law_(tau, sg);
            if (!rel_near(p1, tilted_eff, 1e-4f)) {
                std::cerr << "  tau=" << tau << " sigma=" << sg
                          << " p1=" << p1 << " tilted_eff=" << tilted_eff << "\n";
                return fail("p=1 does not reproduce the sigma-tilted schedule");
            }

            // (2) p = 0 must reproduce the deployed Cubic input up to the
            //     constant gain sigma_ref, with no residual sigma or tau
            //     dependence: the deployed law is the family's p = 0 member.
            f.rs_sigma_exponent_ = 0.0f;
            const float p0 = f.rs_target_from_law_(tau, sg);
            if (!rel_near(p0, sigma_ref * cubic_eff, 1e-3f)) {
                std::cerr << "  tau=" << tau << " sigma=" << sg
                          << " p0=" << p0 << " expected=" << sigma_ref * cubic_eff << "\n";
                return fail("p=0 is not the deployed Cubic schedule up to sigma_ref");
            }

            // (3) the p=0 / p=1 ratio must not depend on tau, otherwise the
            //     ablation confounds amplitude shape with overall gain.
            const float ratio = p0 / tilted_eff;
            if (!rel_near(ratio, sigma_ref / sg, 1e-3f)) {
                std::cerr << "  tau=" << tau << " sigma=" << sg
                          << " ratio=" << ratio << " expected=" << sigma_ref / sg << "\n";
                return fail("strong/tilted ratio is not sigma_ref/sigma");
            }
            if (i == 0) ratio_at_first_sigma = ratio * sigmas[0];
            else if (!rel_near(ratio * sg, ratio_at_first_sigma, 1e-3f))
                return fail("strong/tilted ratio acquired a sigma-dependent factor");
        }
    }

    // (4) transition-law asymptotics.  zeta = 2 sigma^2 tau / r_a.
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
        // distinguishing the two branches.  The deployed input must also be
        // flat in sigma_aw, which is the change the retune introduced: the
        // base is c_R tau^3 and the cadence renormalization carries no
        // amplitude either.
        f.rs_law_ = RSAdaptationLaw::Cubic;
        float cubic_coeff = 0.0f;
        for (float t : {1.5f, 3.0f, 6.0f}) {
            const float cadence =
                std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S / f.pseudo_update_period_for_(t));
            const float eff = f.rs_target_from_law_(t, 1.0f) * cadence;
            const float coeff = eff / std::pow(t, 2.5f);
            if (cubic_coeff == 0.0f) cubic_coeff = coeff;
            else if (!rel_near(coeff, cubic_coeff, 5e-3f)) {
                std::cerr << "  tau=" << t << " coeff=" << coeff
                          << " expected=" << cubic_coeff << "\n";
                return fail("deployed Cubic input is not proportional to tau^(5/2)");
            }
            for (float sg : {0.05f, 0.353f, 1.5f, 4.0f}) {
                const float eff_sg = f.rs_target_from_law_(t, sg) * cadence;
                if (!rel_near(eff_sg, eff, 1e-6f)) {
                    std::cerr << "  tau=" << t << " sigma=" << sg
                              << " eff=" << eff_sg << " expected=" << eff << "\n";
                    return fail("deployed Cubic input still depends on sigma_aw");
                }
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
