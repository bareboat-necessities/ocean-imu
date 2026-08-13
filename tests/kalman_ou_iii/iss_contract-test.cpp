// Pins the implementation-side hypotheses of the OU-III local ISS theorem
// (doc/kalman_ou_iii/w3d-iss-stability.tex-part and its hypothesis audit).
//
// Environmental excitation and the absence of resets cannot be made true by a
// unit test and remain explicit Live-interval assumptions.  What *can* regress
// mechanically is checked here: the 21-state proof configuration, the finite
// residual accelerometer-bias OU contraction, the proof-compatible covariance
// policy, the wrapper clamps that witness the bounded exogenous schedule, and
// the pseudo-update cadence contract underlying the translational
// observability lemma.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std = 9.80665f;

using Filter = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;

static bool near(float a, float b, float rel = 2e-5f) {
    return std::fabs(a - b) <= rel * std::max(1.0f, std::max(std::fabs(a), std::fabs(b)));
}

static int fail(const char* msg) {
    std::cerr << "FAIL: " << msg << "\n";
    return 1;
}

// A_3(t) = 1/2 \int_0^t (t-s)^2 E(s) ds with E(s) = exp(-s/tau), evaluated by
// composite Simpson quadrature.  The closed form is a difference of terms of
// order tau^3, so for t << tau it loses most of its significant digits; the
// quadrature form mirrors the definition used in the lemma and stays accurate.
static double A3_of(double t, double tau) {
    if (t <= 0.0) return 0.0;
    const int n = 20000;                 // even
    const double dt = t / double(n);
    double acc = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double s = double(i) * dt;
        const double f = (t - s) * (t - s) * std::exp(-s / tau);
        const double w = (i == 0 || i == n) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
        acc += w * f;
    }
    return 0.5 * acc * dt / 3.0;
}

// Observation matrix of the S=0 pseudo-measurement over four update instants,
// in state order [v, p, S, a]:  row_j = [t_j^2/2, t_j, 1, A_3(t_j)].
static double detOS(const double t[4], double tau) {
    Eigen::Matrix4d O;
    for (int j = 0; j < 4; ++j) {
        O(j, 0) = 0.5 * t[j] * t[j];
        O(j, 1) = t[j];
        O(j, 2) = 1.0;
        O(j, 3) = A3_of(t[j], tau);
    }
    return O.determinant();
}

int main() {
    Filter f(true);
    f.initialize(Eigen::Vector3f::Constant(0.0148f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.25f));
    if (!f.mekf_) return fail("OU-III core was not created");

    // ---- Proof configuration -------------------------------------------
    // 21 states: attitude(3), gyro bias(3), v/p/S/a_w(12), residual accel
    // bias(3).  The 21-state detectability lemma is stated for exactly this.
    if (f.mekf_->covariance_full().rows() != 21 ||
        f.mekf_->covariance_full().cols() != 21) {
        return fail("default OU-III state dimension is no longer 21");
    }

    // A finite tau_b is what makes the residual-bias block UES; the random-walk
    // limit is explicitly outside the theorem.
    if (!near(f.mekf_->get_acc_bias_time_constant(), 5000.0f))
        return fail("default residual accel-bias OU time constant changed");
    if (!(f.mekf_->get_acc_bias_time_constant() < std::numeric_limits<float>::max()))
        return fail("residual accel-bias block degenerated to a random walk");

    // Only the prediction-time PSD inflation path is covered by the proof.
    if (f.mekf_->legacy_aw_covariance_replacement())
        return fail("legacy raw a_w covariance replacement became default");
    if (!f.periodicAwCovarianceSync())
        return fail("default PSD a_w covariance synchronization was disabled");

    // ---- Bounded exogenous schedule -------------------------------------
    if (!near(MIN_TAU_S, 0.02f) || !near(MAX_TAU_S, 12.0f) ||
        !near(MAX_SIGMA_A, 6.0f) ||
        !near(MIN_R_S, 0.15f) || !near(MAX_R_S, 400.0f) ||
        !near(ACC_NOISE_FLOOR_SIGMA_DEFAULT, 0.12f)) {
        return fail("proof-relevant OU-III wrapper clamps changed");
    }
    if (!near(f.S_factor_, 1.87f))
        return fail("default OU-III acceleration anisotropy no longer matches proof audit");
    if (!near(f.R_S_xy_factor_, 1.0f))
        return fail("default OU-III integral regularizer is no longer isotropic");
    // The deployed pair is (S_sigma, rho_xy) = (1.87, 1), which is *not* the
    // equal-normalized-pole point of the similarity law: sigma_S ~ sigma_aw
    // tau^3 puts that at rho_xy = S_sigma.  docs/ou-iii-anisotropy-consistency.md
    // measures both and keeps rho_xy = 1, so the setter must be able to express
    // rho_xy > 1 for the ablation to mean anything.  It was clamped to 1, which
    // silently turned that configuration back into the deployed one.
    {
        Filter probe;
        probe.setRSXYFactor(1.87f);
        if (!near(probe.R_S_xy_factor_, 1.87f))
            return fail("rho_xy > 1 is not expressible; the similarity-law ablation is a no-op");
        probe.setRSXYFactor(9.0f);
        if (!near(probe.R_S_xy_factor_, 4.0f))
            return fail("rho_xy upper bound changed");
        probe.setRSXYFactor(-1.0f);
        if (!near(probe.R_S_xy_factor_, 0.0f))
            return fail("rho_xy lower bound changed");
    }

    // ---- Pseudo-update cadence contract ---------------------------------
    // The audit quotes T_S = clip(c_T * tau_applied, T_min, T_max).  Pin the
    // constants that turn the bounded-gap hypothesis into a checked statement.
    if (!near(PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT, 0.005f) ||
        !near(PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT, 0.25f) ||
        !near(PSEUDO_UPDATE_TAU_RATIO_DEFAULT, 0.015f / 1.1f)) {
        return fail("deployed OU-III pseudo-update cadence constants changed");
    }
    // Over the clamped tau envelope the configured period stays inside the
    // interval quoted by the audit.
    const float TS_at_tau_min =
        std::min(std::max(PSEUDO_UPDATE_TAU_RATIO_DEFAULT * MIN_TAU_S,
                          PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT),
                 PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT);
    const float TS_at_tau_max =
        std::min(std::max(PSEUDO_UPDATE_TAU_RATIO_DEFAULT * MAX_TAU_S,
                          PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT),
                 PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT);
    if (!near(TS_at_tau_min, 0.005f) || TS_at_tau_max > 0.165f ||
        TS_at_tau_max < 0.160f) {
        return fail("configured pseudo-update period left the audited envelope");
    }

    // ---- Realized cadence, observed from the running scheduler -----------
    f.startup_stage_ = Filter::StartupStage::Live;
    f.mekf_->set_linear_block_enabled(true);
    if (!f.mekf_->linear_block_enabled())
        return fail("Live state did not enable the OU-III linear block");

    const float h = 0.005f;
    const float Tp = 7.0f;
    const float w = 2.0f * float(M_PI) / Tp;

    std::vector<double> fire_times;
    double t_now = 0.0;
    float prev_elapsed = f.mekf_->pseudo_update_elapsed_s_;
    for (int k = 0; k < 60000; ++k) {
        const Eigen::Vector3f gyro = Eigen::Vector3f::Zero();
        const Eigen::Vector3f acc(0.0f, 0.0f, -g_std + 1.2f * std::sin(w * float(t_now)));
        f.updateTime(h, gyro, acc);
        t_now += double(h);
        const float elapsed = f.mekf_->pseudo_update_elapsed_s_;
        // The elapsed accumulator only decreases when an update is applied.
        if (elapsed < prev_elapsed && t_now > 200.0) fire_times.push_back(t_now);
        prev_elapsed = elapsed;
    }
    if (fire_times.size() < 64)
        return fail("pseudo-update scheduler did not run");

    const float tau_applied = f.getTauApplied();
    const float TS_cfg = f.getPseudoUpdatePeriodSec();
    if (!(tau_applied >= MIN_TAU_S && tau_applied <= MAX_TAU_S))
        return fail("applied tau left its clamp interval");
    // T_S/tau is the similarity group the deployed cadence holds fixed.
    if (!near(TS_cfg / tau_applied, PSEUDO_UPDATE_TAU_RATIO_DEFAULT, 1e-3f))
        return fail("realized T_S/tau is not the configured self-similar ratio");

    // Bounded-gap witness of Lemma "Uniform observability of one OU-III axis":
    // a gap is never shorter than one IMU step, never longer than one
    // configured period plus one step.
    double gap_min = 1e30, gap_max = 0.0;
    for (size_t i = 1; i < fire_times.size(); ++i) {
        const double g = fire_times[i] - fire_times[i - 1];
        gap_min = std::min(gap_min, g);
        gap_max = std::max(gap_max, g);
    }
    const double h_d = double(h);
    if (gap_min < h_d - 1e-9)
        return fail("realized pseudo-update gap fell below one IMU step");
    if (gap_max > double(PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT) + h_d + 1e-9)
        return fail("realized pseudo-update gap exceeded T_S_max + h");

    // The gaps are quantized to the IMU grid and are NOT equal: this is the
    // regression guard for the paper's bounded-gap (not exact-cadence)
    // formulation.  If a future change restores an exactly periodic cadence
    // this fires, and the proof text should be revisited rather than silently
    // left stronger than needed.
    if (gap_max - gap_min < 0.5 * h_d)
        return fail("realized cadence became (near) exactly periodic; revisit the ISS cadence text");

    // ---- Determinant witness for the translational UCO lemma ------------
    // Take four consecutive realized instants and check the lemma's bound
    //      |det O_S| >= T_-^6 exp(-3 T_+ / tau_min)
    // holds on the unequally spaced grid the scheduler actually produces.
    for (size_t i = 0; i + 3 < fire_times.size(); i += 7) {
        double t[4];
        for (int j = 0; j < 4; ++j) t[j] = fire_times[i + j] - fire_times[i];
        double lo = 1e30, hi = 0.0;
        for (int j = 1; j < 4; ++j) {
            const double g = t[j] - t[j - 1];
            lo = std::min(lo, g);
            hi = std::max(hi, g);
        }
        const double bound = std::pow(lo, 6.0) * std::exp(-3.0 * hi / double(tau_applied));
        const double det = std::fabs(detOS(t, double(tau_applied)));
        if (!(det > 0.0) || !(det >= bound)) {
            std::cerr << "  det=" << det << " bound=" << bound
                      << " gaps=[" << lo << "," << hi << "]\n";
            return fail("OU-III four-point observability determinant violated the lemma bound");
        }
    }

    // ---- Effective r_S bounds after cadence renormalization --------------
    // apply_RS_tune_ scales the clamped base by sqrt(T_S0 / T_S) and does not
    // re-clamp, so the audited enclosure is [0.121, 693] m*s rather than the
    // base interval [0.15, 400].
    const float scale_lo = std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S / TS_at_tau_max);
    const float scale_hi = std::sqrt(PSEUDO_UPDATE_PERIOD_NOMINAL_S / TS_at_tau_min);
    if (!(scale_lo > 0.30f && scale_lo < 0.31f) || !near(scale_hi, std::sqrt(3.0f), 1e-3f))
        return fail("cadence information-rate factor left its audited range");
    const float rs_eff_lo = MIN_R_S * scale_lo;
    const float rs_eff_hi = MAX_R_S * scale_hi;
    if (!(rs_eff_lo > 0.044f && rs_eff_lo < 0.046f) || !(rs_eff_hi > 690.0f && rs_eff_hi < 694.0f))
        return fail("effective r_S enclosure no longer matches the ISS audit");
    const float rs_live = std::sqrt(f.mekf_->R_S(2, 2));
    if (!(rs_live >= rs_eff_lo && rs_live <= rs_eff_hi) || !std::isfinite(rs_live))
        return fail("live r_S left the audited effective enclosure");

    // ---- Residual accelerometer-bias OU contraction ----------------------
    f.mekf_->set_linear_block_enabled(false);
    f.mekf_->set_acc_bias_updates_enabled(true);
    f.mekf_->set_acc_bias_time_constant(2.0f);
    const Eigen::Vector3f b0(0.20f, -0.10f, 0.05f);
    f.mekf_->set_initial_acc_bias(b0);
    f.mekf_->time_update(Eigen::Vector3f::Zero(), 0.10f);
    const Eigen::Vector3f expected = b0 * std::exp(-0.10f / 2.0f);
    if (!f.mekf_->get_acc_bias().isApprox(expected, 2e-6f))
        return fail("active residual accel-bias OU block is not contractive");

    // In the certified post-unlock Live regime the residual-bias OU state must
    // be propagating, not held by the startup freeze.
    f.accel_bias_locked_ = false;
    f.startup_stage_ = Filter::StartupStage::TunerWarm;
    f.enterLive_();
    if (!f.mekf_->acc_bias_updates_enabled_)
        return fail("post-unlock Live state did not enable accel-bias OU propagation");
    if (!f.mekf_->linear_block_enabled())
        return fail("default Live state did not re-enable the OU-III linear block");

    std::cout << "OU-III ISS proof contract passed"
              << " (tau=" << tau_applied
              << " s, T_S=" << TS_cfg
              << " s, realized gaps [" << gap_min << "," << gap_max << "] s)\n";
    return 0;
}
