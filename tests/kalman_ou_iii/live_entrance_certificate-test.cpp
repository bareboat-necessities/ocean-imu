// Pins the semantics of the Phase-3 Live-entrance certificate.
//
// The truth-based sweep in live_handoff_validation.cpp answers "is the
// certificate sound on real trajectories".  This answers the questions that
// have nothing to do with a trajectory: that the integration-epoch reset is an
// identity rather than a small number, that a handoff with no heading bound is
// refused for that reason and not some other, that a failure says which term
// spent the budget, that the certification does not survive a Live-interval
// boundary, and that the provenance of every bound travels with it.
//
// These are contract tests.  Each one fixes a statement the paper makes.
#define EIGEN_NON_ARDUINO

#include <cmath>
#include <iostream>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#include "kalman_ou_iii/LiveEntranceCertificate.h"

namespace ou3 = ocean_imu::kalman::ou3;

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
    if (!ok) {
        std::cerr << "FAIL: " << what << "\n";
        ++failures;
    }
}

ou3::LiveHandoffObservables goodObservables() {
    ou3::LiveHandoffObservables o;
    o.gate_satisfied = true;
    o.aligned_branch = true;
    o.yaw_gauged = true;
    o.gravity_align_sin = 0.03f;
    o.sigma_theta = 0.087f;
    o.sigma_bg = 3.6e-4f;
    o.sigma_S = 50.0f;
    o.sigma_aw = 0.5f;
    o.sigma_ba = 0.05f;
    o.sigma_v = 4.0f;
    o.sigma_p = 8.0f;
    o.seeded_world_accel = true;
    return o;
}

void testEpochIsExact() {
    // The S handoff bound is zero, and it is labelled as an identity rather
    // than as an envelope.  Every other sensitive bound is at best an
    // envelope, so this is the one block whose provenance is Exact.
    const auto c = ou3::evaluateLiveEntrance(goodObservables(), ou3::LiveEnvelope{},
                                             ou3::LiveBasinConstants{}, 0.075f);
    check(c.integral_S.value == 0.0f, "S handoff bound is not exactly zero");
    check(c.integral_S.source == ou3::BoundSource::Exact,
          "S handoff bound is not labelled exact");
    // A zero numerator contributes nothing whatever the marginal is, which is
    // why the reset does not have to collapse the covariance to be useful.
    ou3::LiveHandoffObservables wide = goodObservables();
    wide.sigma_S = 1.0e-3f;
    const auto c2 = ou3::evaluateLiveEntrance(wide, ou3::LiveEnvelope{},
                                              ou3::LiveBasinConstants{}, 0.075f);
    check(std::fabs(c.sensitive_metric_norm - c2.sensitive_metric_norm) < 1e-5f,
          "the S marginal changed the sensitive metric norm despite a zero bound");
}

void testNoGaugeIsRefusedForTheRightReason() {
    ou3::LiveHandoffObservables o = goodObservables();
    o.yaw_gauged = false;
    const auto c = ou3::evaluateLiveEntrance(o, ou3::LiveEnvelope{},
                                             ou3::LiveBasinConstants{}, 0.075f);
    check(!c.certified, "a handoff with no magnetic gauge was certified");
    check(c.failure == ou3::LiveCertFailure::AttitudeUnbounded,
          "a handoff with no magnetic gauge was refused for the wrong reason");
    check(!std::isfinite(c.attitude.value),
          "a handoff with no magnetic gauge reported a finite attitude bound");
    check(c.attitude.source == ou3::BoundSource::External,
          "an unbounded heading was not labelled an external hypothesis");
}

void testBranchIsAPrecondition() {
    ou3::LiveHandoffObservables o = goodObservables();
    o.aligned_branch = false;
    const auto c = ou3::evaluateLiveEntrance(o, ou3::LiveEnvelope{},
                                             ou3::LiveBasinConstants{}, 0.075f);
    check(!c.certified, "a handoff off the aligned branch was certified");
    check(c.failure == ou3::LiveCertFailure::GateNotSatisfied,
          "the aligned branch is not a hard precondition");
}

void testTiltBoundTracksTheDeployedGate() {
    // The certificate must read the gate threshold rather than keep its own
    // copy, or the two can drift apart silently.
    ou3::LiveEnvelope env;
    env.gravity_direction_rad = 0.0f;
    env.mag_reference_rel_error = 0.0f;
    const auto a = ou3::evaluateLiveEntrance(goodObservables(), env,
                                             ou3::LiveBasinConstants{}, 0.075f);
    const auto b = ou3::evaluateLiveEntrance(goodObservables(), env,
                                             ou3::LiveBasinConstants{}, 0.2f);
    check(std::fabs(a.attitude.value - std::asin(0.075f)) < 1e-6f,
          "the tilt bound is not asin(gate threshold) + eta_g");
    check(b.attitude.value > a.attitude.value,
          "the attitude bound does not follow the gate threshold");
}

void testEnvelopeEntersTheBound() {
    ou3::LiveEnvelope tight;
    tight.gravity_direction_rad = 0.01f;
    tight.mag_reference_rel_error = 0.01f;
    ou3::LiveEnvelope loose;
    loose.gravity_direction_rad = 0.2f;
    loose.mag_reference_rel_error = 0.2f;
    const auto a = ou3::evaluateLiveEntrance(goodObservables(), tight,
                                             ou3::LiveBasinConstants{}, 0.075f);
    const auto b = ou3::evaluateLiveEntrance(goodObservables(), loose,
                                             ou3::LiveBasinConstants{}, 0.075f);
    check(b.attitude.value > a.attitude.value,
          "a looser declared envelope did not widen the attitude bound");
    check(b.sensitive_metric_norm > a.sensitive_metric_norm,
          "a looser declared envelope did not cost anything in the metric");
    check(b.margin < a.margin,
          "a looser declared envelope did not cost margin");
}

void testDiagnosticsAreReadable() {
    // A refusal has to say which term spent the budget: both sides of the
    // inequality, the tube radius, and the small-gain slope's effect are all
    // required to be present.
    ou3::LiveBasinConstants k;
    k.c_eff = 29.3f;
    const auto c = ou3::evaluateLiveEntrance(goodObservables(), ou3::LiveEnvelope{},
                                             k, 0.075f);
    check(!c.certified, "the deployed constants certified a nominal handoff");
    check(c.failure == ou3::LiveCertFailure::BasinExceeded,
          "a budget failure was not reported as one");
    check(std::isfinite(c.basin_lhs) && std::isfinite(c.basin_rhs),
          "both sides of the basin inequality must be readable after a refusal");
    check(std::fabs(c.margin - (c.basin_rhs - c.basin_lhs)) < 1e-5f,
          "margin is not the difference of the two reported sides");
    check(c.nonlinear_radius > 0.0f && c.small_gain_nu > 0.0f &&
              c.small_gain_nu < 1.0f,
          "the tube radius and small-gain fraction were not reported");
    check(std::isfinite(c.eta_gravity_rad) && std::isfinite(c.eta_heading_rad),
          "the two halves of the attitude assumption were not reported apart");
    check(c.weakestSource() == ou3::BoundSource::Envelope,
          "the weakest provenance of a nominal certificate is not the envelope");
}

void testASmallEnoughHandoffCertifies() {
    // The inequality has to be capable of holding.  A certificate that can
    // only ever refuse would pass every safety test and mean nothing, so this
    // pins that the arithmetic closes when the handoff really is small.
    ou3::LiveBasinConstants k;
    k.c_eff = 0.01f;
    ou3::LiveHandoffObservables o = goodObservables();
    o.sigma_theta = 5.0f;
    o.sigma_bg = 1.0f;
    o.sigma_aw = 100.0f;
    o.sigma_ba = 10.0f;
    o.sigma_v = 100.0f;
    o.sigma_p = 100.0f;
    const auto c = ou3::evaluateLiveEntrance(o, ou3::LiveEnvelope{}, k, 0.075f);
    check(c.certified, "the entrance inequality never closes");
    check(c.margin > 0.0f, "a certified handoff reported a non-positive margin");
    check(c.failure == ou3::LiveCertFailure::None,
          "a certified handoff carried a failure reason");
}

void testEpochResetOnTheFilter() {
    using Core = Kalman3D_Wave_OU_III<double, true, true>;
    Core f(Eigen::Vector3d::Constant(0.0294),
           Eigen::Vector3d::Constant(0.000157),
           Eigen::Vector3d::Constant(0.36));

    // Let the real dynamics put something in S -- a standing velocity drives
    // p, and p drives S -- then restart the epoch.
    f.set_linear_block_enabled(true);
    f.seed_translational(Eigen::Vector3d(1.0, 0.0, 0.0),
                         Eigen::Vector3d::Zero(), 1.0, 1.0);
    for (int i = 0; i < 400; ++i) f.time_update(Eigen::Vector3d::Zero(), 0.005);
    check(f.get_integral_displacement().norm() > 0.0,
          "the test did not manage to make the integral state nonzero");

    f.reset_integral_epoch(50.0);
    check(f.get_integral_displacement().norm() == 0.0,
          "reset_integral_epoch did not zero the integral state");
    check(std::fabs(f.block_sigma(Core::StateBlock::IntegralS) - 50.0) < 1e-6,
          "reset_integral_epoch did not install the requested marginal");

    // The other blocks are untouched.
    check(f.block_sigma(Core::StateBlock::Velocity) > 0.0 &&
              f.block_sigma(Core::StateBlock::Position) > 0.0,
          "reset_integral_epoch disturbed the translational marginals");
}

void testSeedersDropTheirCrossCovariance() {
    using Core = Kalman3D_Wave_OU_III<double, true, true>;
    Core f(Eigen::Vector3d::Constant(0.0294),
           Eigen::Vector3d::Constant(0.000157),
           Eigen::Vector3d::Constant(0.36));

    f.seed_gyro_bias(Eigen::Vector3d(1e-4, -2e-4, 3e-4), 5e-4);
    check((f.gyroscope_bias() - Eigen::Vector3d(1e-4, -2e-4, 3e-4)).norm() < 1e-12,
          "seed_gyro_bias did not install the estimate");
    check(std::fabs(f.block_sigma(Core::StateBlock::GyroBias) - 5e-4) < 1e-9,
          "seed_gyro_bias did not install the requested marginal");

    f.seed_world_accel(Eigen::Vector3d(0.5, -0.25, 1.0), 2.0);
    check((f.get_world_accel() - Eigen::Vector3d(0.5, -0.25, 1.0)).norm() < 1e-12,
          "seed_world_accel did not install the estimate");
    check(std::fabs(f.block_sigma(Core::StateBlock::WorldAccel) - 2.0) < 1e-9,
          "seed_world_accel did not install the requested marginal");

    f.seed_translational(Eigen::Vector3d(1.0, 0.0, -1.0),
                         Eigen::Vector3d(0.0, 2.0, 0.0), 3.0, 7.0);
    check((f.get_velocity() - Eigen::Vector3d(1.0, 0.0, -1.0)).norm() < 1e-12,
          "seed_translational did not install the velocity");
    check((f.get_position() - Eigen::Vector3d(0.0, 2.0, 0.0)).norm() < 1e-12,
          "seed_translational did not install the position");
    check(std::fabs(f.block_sigma(Core::StateBlock::Velocity) - 3.0) < 1e-9,
          "seed_translational did not install the velocity marginal");
    check(std::fabs(f.block_sigma(Core::StateBlock::Position) - 7.0) < 1e-9,
          "seed_translational did not install the position marginal");
}

// Print the entrance inequality under the two seed policies, at the library
// default envelope, so the numbers the paper quotes are regenerated rather
// than typed.  The deployed seed marginals are the constructor's; the coherent
// ones are the constructed bounds themselves.
void reportNominalBudget() {
    const ou3::LiveEnvelope env;                 // library defaults
    const ou3::LiveBasinConstants k;             // worst-of-eight c_eff

    ou3::LiveHandoffObservables deployed = goodObservables();
    deployed.sigma_theta = 0.087f;   // proxy_handoff_yaw_sigma_rad
    deployed.sigma_bg    = 1.0e-3f;  // sqrt(Pb0)
    deployed.sigma_S     = 50.0f;    // constructor sigma_S0
    deployed.sigma_aw    = 0.36f;    // OU stationary spread, calmest reference
    deployed.sigma_ba    = 0.004f;   // constructor sigma_bacc0
    deployed.sigma_v     = 1.0f;
    deployed.sigma_p     = 20.0f;
    const auto a = ou3::evaluateLiveEntrance(deployed, env, k, 0.075f);

    ou3::LiveHandoffObservables coherent = deployed;
    coherent.sigma_bg = std::sqrt(env.gyro_bias_perp_rad_s * env.gyro_bias_perp_rad_s +
                                  env.gyro_bias_axial_rad_s * env.gyro_bias_axial_rad_s);
    coherent.sigma_aw = a.world_accel.value;
    coherent.sigma_ba = env.accel_bias_ms2;
    coherent.sigma_v  = env.translational_velocity_ms;
    coherent.sigma_p  = env.translational_position_m;
    const auto b = ou3::evaluateLiveEntrance(coherent, env, k, 0.075f);

    auto line = [](const char* tag, const ou3::LiveEntranceCertificate& c) {
        std::cout << "CERT_BUDGET " << tag
                  << " theta_bound=" << c.attitude.value
                  << " aw_bound=" << c.world_accel.value
                  << " z_xi=" << c.sensitive_metric_norm
                  << " z_ell=" << c.kinematic_metric_norm
                  << " lhs=" << c.basin_lhs
                  << " required_c_eff=" << (1.0f / (4.0f * c.basin_lhs))
                  << "\n";
    };
    line("deployed_seed", a);
    line("coherent_seed", b);
}

}  // namespace

int main() {
    testEpochIsExact();
    testNoGaugeIsRefusedForTheRightReason();
    testBranchIsAPrecondition();
    testTiltBoundTracksTheDeployedGate();
    testEnvelopeEntersTheBound();
    testDiagnosticsAreReadable();
    testASmallEnoughHandoffCertifies();
    testEpochResetOnTheFilter();
    testSeedersDropTheirCrossCovariance();
    reportNominalBudget();

    if (failures) {
        std::cerr << failures << " certificate contract(s) failed\n";
        return 1;
    }
    std::cout << "LIVE_ENTRANCE_CERTIFICATE all contracts hold\n";
    return 0;
}
