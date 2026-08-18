// Is the lower clamp on either drift-band regularizer a guard, or is it the
// thing actually setting the operating point?
//
// OU-III's r_S floor was 0.4 m*s while its own schedule asked for 0.24 m*s at
// the calibrated H_s = 0.27 m point, so the floor -- not the law -- set the
// regularizer in every low-motion sea. A full sweep of the three tuner
// multipliers left both low-motion records constant to three decimals, because
// nothing the multipliers could say survived the clamp. Dropping the floor to
// 0.15 recovered 8-10% on those records.
//
// OU-II's laws are of the same shape -- the deployed PhysicalMSE law asks for
// r_p ~ sigma_a,B^(4/5) tau^(12/5)/sqrt(T_S) against OU-III's tau^(24/7) -- so
// the same failure is available to it and has to be checked rather than
// assumed. This test does the checking against the deployed law at the
// operating points the eight scored records actually produce, so it also fails
// if a future re-fit, or a switch of deployed law, walks the schedule down into
// the clamp.
//
// It is asked of the *deployed* law rather than of a hard-coded schedule,
// because a floor is only clear of the law that ships. PhysicalMSE is tighter
// at the small-sea end than the Empirical law it replaced -- that is what a
// distortion penalty scaling with wave energy does -- and moving to it is what
// took the r_p0 floor from 0.05 m to 0.03 m.
//
// Also pinned: under the Empirical law the clamp is applied to the *base*
// tuner value and the cadence renormalization is deliberately applied after
// it, so the filter input is allowed below the base floor once T_S exceeds the
// nominal 15 ms. That ordering is what stops the tau-scaled cadence from
// reintroducing exactly the clipping the floor discussion is about.
// PhysicalMSE reaches the same place by a different route: it evaluates at the
// realized T_S and returns the filter input itself, so there is no second
// factor to re-clamp and the clamp already applies where it matters. Both
// halves are checked.
#define EIGEN_NON_ARDUINO

#include <cmath>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#define private public
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
#undef private

const float g_std = 9.80665f;

namespace {

using Filter = SeaStateFusionFilter_OU_II<TrackerType::KALMANF>;

bool check(bool condition, const char* message) {
    if (!condition) std::cerr << "FAIL: " << message << '\n';
    return condition;
}

// (tau_target, sigma_target) reported by kalman_ou_ii-sim at the end of each
// of the eight scored records, deployed configuration. The smallest-motion
// entries are the ones that decide whether a floor binds.
struct OperatingPoint {
    const char* record;
    float tau;
    float sigma;
};

// Re-measured on the filter that ships.  Both channels are exogenous of the
// regularizer coefficients by construction -- the tuner reads the levelled
// acceleration and the wave-band period, never an estimator state -- and the
// (r_p0, r_v0) re-fit reproduces all sixteen numbers bit-for-bit, which is that
// exogeneity measured rather than assumed.  What had moved them since they were
// last written down is the attitude and magnetic work, which reaches the tuner
// through its levelled input.
constexpr OperatingPoint CALIBRATED[] = {
    {"jonswap  H0.27", 1.28034f, 0.339120f},
    {"jonswap  H1.50", 2.18073f, 0.688910f},
    {"jonswap  H4.00", 3.59057f, 1.06115f},
    {"jonswap  H8.50", 4.22184f, 1.34176f},
    {"pmstokes H0.27", 1.17259f, 0.364507f},
    {"pmstokes H1.50", 2.04482f, 0.746216f},
    {"pmstokes H4.00", 3.28991f, 1.06387f},
    {"pmstokes H8.50", 4.10027f, 1.41327f},
};

// The near-still stress case OU-III used to justify moving its floor. Scaled
// from the smallest calibrated sea by significant wave height at a comparable
// period, which is what makes it the hardest case for a lower clamp.
constexpr float STILL_HS = 0.05f;
constexpr float SMALLEST_CALIBRATED_HS = 0.27f;

// The two pseudo standard deviations the deployed law puts at the filter
// input for a given operating point, including each law's cadence contract.
float law_input_p(Filter& f, float tau, float sigma) {
    float r_p = 0.0f, r_v = 0.0f;
    f.pseudo_targets_from_law_(tau, sigma, r_p, r_v);
    if (f.getPseudoLaw() != PseudoAdaptationLaw::Empirical) return r_p;
    return r_p * std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S
                           / f.pseudo_update_period_for_(tau));
}

float law_input_v(Filter& f, float tau, float sigma) {
    float r_p = 0.0f, r_v = 0.0f;
    f.pseudo_targets_from_law_(tau, sigma, r_p, r_v);
    if (f.getPseudoLaw() != PseudoAdaptationLaw::Empirical) return r_v;
    return r_v * std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S
                           / f.pseudo_update_period_for_(tau));
}

// The floors have to stay clear of the schedule across the calibrated
// envelope, or they, and not the law, are setting the regularizer.
bool test_floors_are_clear_of_the_schedule() {
    bool ok = true;

    Filter f(false);

    float worst_p0 = 1e30f;
    float worst_v0 = 1e30f;
    const char* worst_p0_rec = "";
    const char* worst_v0_rec = "";
    float worst_tau = 0.0f;
    float worst_sigma = 0.0f;

    // Ask the deployed law, and compare at the filter input -- the quantity
    // the clamps actually see once each law's cadence contract is applied.
    for (const auto& op : CALIBRATED) {
        const float r_p0 = law_input_p(f, op.tau, op.sigma);
        const float r_v0 = law_input_v(f, op.tau, op.sigma);
        if (r_p0 < worst_p0) {
            worst_p0 = r_p0; worst_p0_rec = op.record;
            worst_tau = op.tau; worst_sigma = op.sigma;
        }
        if (r_v0 < worst_v0) { worst_v0 = r_v0; worst_v0_rec = op.record; }
    }

    ok &= check(worst_p0 > 2.0f * MIN_R_p0_std,
                "the r_p0 floor must stay clear of the calibrated schedule");
    ok &= check(worst_v0 > 2.0f * MIN_R_v0_std,
                "the r_v0 floor must stay clear of the calibrated schedule");

    std::cout << "  smallest calibrated demand: r_p0 = " << worst_p0
              << " m (" << worst_p0_rec << ") against a " << MIN_R_p0_std
              << " floor, " << (worst_p0 / MIN_R_p0_std) << "x clear\n";
    std::cout << "  smallest calibrated demand: r_v0 = " << worst_v0
              << " m/s (" << worst_v0_rec << ") against a " << MIN_R_v0_std
              << " floor, " << (worst_v0 / MIN_R_v0_std) << "x clear\n";

    // The near-still case is where a floor first starts to bite, so it is
    // checked separately and with a smaller margin.  Scale the *input* to the
    // law rather than its output: the amplitude exponent is 4/5 under the
    // deployed law and 1 under Empirical, so scaling the output would assume
    // the schedule that a change of law is allowed to move.
    const float still_scale = STILL_HS / SMALLEST_CALIBRATED_HS;
    const float still_p0 = law_input_p(f, worst_tau, worst_sigma * still_scale);
    const float still_v0 = law_input_v(f, worst_tau, worst_sigma * still_scale);

    ok &= check(still_p0 > MIN_R_p0_std,
                "the r_p0 floor must not clip the near-still stress case");
    ok &= check(still_v0 > MIN_R_v0_std,
                "the r_v0 floor must not clip the near-still stress case");

    std::cout << "  near-still (Hs=" << STILL_HS << " m) demand: r_p0 = "
              << still_p0 << " m, r_v0 = " << still_v0 << " m/s\n";

    return ok;
}

// The clamp applies to the base tuner value; the information-rate
// renormalization is applied afterwards and is deliberately not re-clamped.
// Re-clamping it would put the clipping straight back at the small-sea end,
// where T_S > 15 ms drives the filter input below the base floor on purpose.
bool test_cadence_renormalization_is_not_reclamped() {
    bool ok = true;

    Filter f(false);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.25f));
    f.initialize_from_acc(Eigen::Vector3f(0.0f, 0.0f, -g_std));
    f.startup_stage_ = Filter::StartupStage::Live;
    f.mekf_->set_linear_block_enabled(true);
    // The clamp-then-renormalize ordering is a property of the Empirical law's
    // cadence contract, so select it explicitly rather than depending on which
    // law happens to ship.
    f.setPseudoLaw(PseudoAdaptationLaw::Empirical);

    // Sit the base value exactly on the floor and stretch the cadence past
    // nominal, which is what a long-period low sea does.
    f.tune_.tau_applied = 4.0f;
    f.tune_.R_p0_std_applied = MIN_R_p0_std;
    f.tune_.R_v0_std_applied = MIN_R_v0_std;
    f.apply_ou_tune_(false);
    f.apply_R_p0_tune_();
    f.apply_R_v0_tune_();

    const float period = f.getPseudoUpdatePeriodSec();
    ok &= check(period > PSEUDO_UPDATE_PERIOD_NOMINAL_S,
                "tau = 4 s must stretch the cadence past the nominal 15 ms");

    const float scale = std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S / period);
    const float applied_p0 = std::sqrt(f.mekf_->R_p0(2, 2));
    const float applied_v0 = std::sqrt(f.mekf_->R_v0(2, 2));

    ok &= check(applied_p0 < MIN_R_p0_std,
                "the filter input must be allowed below the base r_p0 floor");
    ok &= check(applied_v0 < MIN_R_v0_std,
                "the filter input must be allowed below the base r_v0 floor");
    ok &= check(std::fabs(applied_p0 - MIN_R_p0_std * scale) < 1e-6f,
                "the r_p0 filter input must be exactly the renormalized base");
    ok &= check(std::fabs(applied_v0 - MIN_R_v0_std * scale) < 1e-6f,
                "the r_v0 filter input must be exactly the renormalized base");

    std::cout << "  Empirical: at tau = 4 s, T_S = " << period
              << " s, base floor " << MIN_R_p0_std << " -> filter input "
              << applied_p0 << " m\n";

    // PhysicalMSE has no second factor to re-clamp: the law already contains
    // the realized T_S, so its clamped value *is* the filter input.  Applying
    // the Empirical renormalization on top would divide the regularizer by
    // sqrt(T_S/T_0) a second time, which at tau = 4 s is a 90 % error.
    f.setPseudoLaw(PseudoAdaptationLaw::PhysicalMSE);
    f.apply_R_p0_tune_();
    f.apply_R_v0_tune_();
    ok &= check(std::fabs(std::sqrt(f.mekf_->R_p0(2, 2)) - MIN_R_p0_std) < 1e-6f,
                "the PhysicalMSE r_p0 clamp must land at the filter input");
    ok &= check(std::fabs(std::sqrt(f.mekf_->R_v0(2, 2)) - MIN_R_v0_std) < 1e-6f,
                "the PhysicalMSE r_v0 clamp must land at the filter input");

    return ok;
}

}  // namespace

int main() {
    bool ok = true;
    ok &= test_floors_are_clear_of_the_schedule();
    ok &= test_cadence_renormalization_is_not_reclamped();

    if (!ok) {
        std::cerr << "regularizer floor checks FAILED\n";
        return 1;
    }
    std::cout << "all regularizer floor checks passed\n";
    return 0;
}
