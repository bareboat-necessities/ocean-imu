// Pins the dual pseudo-measurement adaptation-law family of OU-II.
//
// The deployed law is the joint physical-MSE schedule of
// doc/kalman_ou_ii/ou2-dual-regularization-mse.tex,
//     r_p = C_P q_eff^(1/10) sigma_a,B^(4/5) tau^(12/5) / sqrt(T_S),
//     r_v = C_V q_eff^(1/10) sigma_a,B^(4/5) tau^(7/5)  / sqrt(T_S),
// and the note's claims about it are structural rather than numerical.  Five
// of them are pinned here, because the manuscript's comparison stops meaning
// what it says if any regresses:
//   1. the amplitude exponent is 4/5, not the empirical 1;
//   2. the period exponents are 12/5 and 7/5 at the base, i.e. 19/10 and 9/10
//      at the filter input away from the cadence clamps -- both channels
//      moving by the same 2/5 above the empirical law, which is the note's
//      "common shape correction" claim;
//   3. Corollary (optimal pseudo-channel ratio) holds exactly in the
//      implementation: r_p/r_v = (C_P/C_V) tau, with no tau^0 residue.  This
//      is the half of the prediction that does not depend on q_eff or on the
//      cadence, and it is also what lets the velocity channel be free;
//   4. r_p is invariant to c_sigma, because the law divides it back out to
//      recover the physical band RMS.  C_P and c_sigma stay separately
//      identifiable only while that holds;
//   5. the Empirical law stays selectable and stays sigma_aw tau^2 /
//      sigma_aw tau at the base, so the ablation arm is the schedule the
//      calibrated coefficients were fitted for.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#define private public
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
#undef private

using Filter = SeaStateFusionFilter_OU_II<TrackerType::KALMANF>;

static bool rel_near(float a, float b, float rel) {
    return std::fabs(a - b) <= rel * std::max(1e-12f, std::max(std::fabs(a), std::fabs(b)));
}

static int fail(const char* msg) {
    std::cerr << "FAIL: " << msg << "\n";
    return 1;
}

static void make(Filter& f) {
    f.initialize(Eigen::Vector3f::Constant(0.0148f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.25f));
}

// Nominal OU-II operating point of the H_s = 1.5 m JONSWAP reference record.
static constexpr float TAU_REF   = 2.1815f;
static constexpr float SIGMA_REF = 0.6805f;

int main() {
    Filter f(false);
    make(f);

    if (f.getPseudoLaw() != PseudoAdaptationLaw::PhysicalMSE)
        return fail("default OU-II pseudo law is no longer PhysicalMSE");

    // (1) and (2): measure the exponents rather than trusting the expression.
    // The cadence is unclamped at tau = 2.18 s, so the filter-input period
    // exponents are the base ones less the 1/2 the sqrt(T_S) contributes.
    {
        const float d = 1.02f;
        float rp0 = 0.0f, rv0 = 0.0f, rp1 = 0.0f, rv1 = 0.0f;
        f.pseudo_targets_from_law_(TAU_REF, SIGMA_REF, rp0, rv0);

        f.pseudo_targets_from_law_(TAU_REF, SIGMA_REF * d, rp1, rv1);
        const float a_p = std::log(rp1 / rp0) / std::log(d);
        const float a_v = std::log(rv1 / rv0) / std::log(d);
        if (!rel_near(a_p, 0.8f, 2e-3f) || !rel_near(a_v, 0.8f, 2e-3f)) {
            std::cerr << "  amplitude exponents " << a_p << ", " << a_v
                      << " expected 4/5 on both channels\n";
            return fail("PhysicalMSE amplitude exponent is not 4/5");
        }

        f.pseudo_targets_from_law_(TAU_REF * d, SIGMA_REF, rp1, rv1);
        const float t_p = std::log(rp1 / rp0) / std::log(d);
        const float t_v = std::log(rv1 / rv0) / std::log(d);
        if (!rel_near(t_p, 1.9f, 2e-3f) || !rel_near(t_v, 0.9f, 2e-3f)) {
            std::cerr << "  filter-input period exponents " << t_p << ", " << t_v
                      << " expected 19/10 and 9/10\n";
            return fail("PhysicalMSE period exponents are not 19/10 and 9/10");
        }
        // The note's central structural claim: both channels sit the same 2/5
        // above the empirical tau^(3/2) and tau^(1/2), so the theory asks for
        // a common shape correction and not a new relative law.
        if (!rel_near(t_p - 1.5f, 0.4f, 5e-3f) ||
            !rel_near(t_v - 0.5f, 0.4f, 5e-3f)) {
            return fail("PhysicalMSE does not sit a common tau^(2/5) above Empirical");
        }
    }

    // (3): the channel ratio is exactly (C_P/C_V) tau at every operating point,
    // including inside the cadence clamps, where T_S stops tracking tau.  The
    // small-sea end of the envelope is below the 5 ms floor at tau = 0.3 s.
    for (float tau : {0.3f, 1.1f, TAU_REF, 4.2f, 11.0f}) {
        for (float sigma : {0.05f, SIGMA_REF, 3.0f}) {
            float rp = 0.0f, rv = 0.0f;
            f.pseudo_targets_from_law_(tau, sigma, rp, rv);
            if (!rel_near(rp / rv, f.getPseudoMseRatio() * tau, 1e-5f)) {
                std::cerr << "  tau=" << tau << " sigma=" << sigma
                          << " ratio " << rp / rv << " expected "
                          << f.getPseudoMseRatio() * tau << "\n";
                return fail("PhysicalMSE channel ratio is not (C_P/C_V) tau");
            }
        }
    }

    // (4): c_sigma moves the OU prior, never the pseudo strength.
    {
        Filter g(false);
        make(g);
        float rp_a = 0.0f, rv_a = 0.0f, rp_b = 0.0f, rv_b = 0.0f;
        g.setSigmaCoeff(0.85f);
        g.pseudo_targets_from_law_(TAU_REF, SIGMA_REF, rp_a, rv_a);
        g.setSigmaCoeff(1.70f);
        // The physical band RMS is what the law depends on, so doubling
        // c_sigma at a doubled sigma_aw must leave both channels alone.
        g.pseudo_targets_from_law_(TAU_REF, SIGMA_REF * 2.0f, rp_b, rv_b);
        if (!rel_near(rp_a, rp_b, 1e-5f) || !rel_near(rv_a, rv_b, 1e-5f))
            return fail("PhysicalMSE is not invariant to c_sigma at fixed sigma_a,B");
    }

    // (5): the Empirical arm is still the calibrated schedule, exactly.
    {
        Filter g(false);
        make(g);
        g.setPseudoLaw(PseudoAdaptationLaw::Empirical);
        if (g.getPseudoLaw() != PseudoAdaptationLaw::Empirical)
            return fail("Empirical law is no longer selectable");
        float rp = 0.0f, rv = 0.0f;
        g.pseudo_targets_from_law_(TAU_REF, SIGMA_REF, rp, rv);
        if (!rel_near(rp, 0.65f * SIGMA_REF * TAU_REF * TAU_REF, 1e-6f) ||
            !rel_near(rv, 1.30f * SIGMA_REF * TAU_REF, 1e-6f))
            return fail("Empirical base pair is not c_p sigma tau^2 / c_v sigma tau");
        // Empirical states a base at the reference cadence, so it is the law
        // that gets the sqrt(T_0/T_S) renormalization; PhysicalMSE must not.
        g.tune_.tau_applied = TAU_REF;
        g.apply_pseudo_update_cadence_();
        const float emp_scale = g.pseudo_update_information_rate_scale_();
        g.setPseudoLaw(PseudoAdaptationLaw::PhysicalMSE);
        const float mse_scale = g.pseudo_update_information_rate_scale_();
        if (!rel_near(emp_scale,
                      std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S
                                / g.getPseudoUpdatePeriodSec()), 1e-5f))
            return fail("Empirical law lost its cadence renormalization");
        if (mse_scale != 1.0f)
            return fail("PhysicalMSE double-counts the pseudo-update cadence");
    }

    // The analytical constants are the deployed ones; a re-fit has to move the
    // header, the ou_validation mirror and this line together.
    if (!rel_near(f.getPseudoMseCoeff(), 0.1116f, 1e-5f) ||
        !rel_near(f.getPseudoMseRatio(), 0.4611f, 1e-5f))
        return fail("deployed PhysicalMSE constants moved without updating the test");

    std::cout << "OU-II dual pseudo-measurement adaptation-law family passed\n";
    return 0;
}
