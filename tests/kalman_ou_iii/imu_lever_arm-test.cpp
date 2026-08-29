#define EIGEN_NON_ARDUINO

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

#include "util/W3dSimCommon.h"

namespace {

bool check(bool condition, const char* message) {
    if (!condition) std::cerr << "FAIL: " << message << '\n';
    return condition;
}

constexpr float DT = 1.0f / 200.0f;
constexpr float TWO_PI = 6.28318530717958647692f;

// A rigid body rocking about one axis, so both rotational terms have a closed
// form the test can compare against.
struct Rocking {
    float amplitude_radps = 0.25f;
    float freq_hz = 0.4f;

    Vector3f omega(float t) const {
        return Vector3f(0.0f, 0.0f,
                        amplitude_radps * std::sin(TWO_PI * freq_hz * t));
    }
    Vector3f alpha(float t) const {
        return Vector3f(0.0f, 0.0f,
                        amplitude_radps * TWO_PI * freq_hz *
                            std::cos(TWO_PI * freq_hz * t));
    }
};

float rms(const std::vector<float>& values) {
    if (values.empty()) return 0.0f;
    double acc = 0.0;
    for (float v : values) acc += double(v) * double(v);
    return float(std::sqrt(acc / double(values.size())));
}

// --- the closed-form kinematics -------------------------------------------

bool test_centripetal_term_points_inward() {
    const Vector3f a = w3d_lever_acceleration(
        Vector3f(0.0f, 0.0f, 2.0f), Vector3f::Zero(), Vector3f(0.1f, 0.0f, 0.0f));
    return check(std::fabs(a.x() + 0.4f) < 1e-6f &&
                     std::fabs(a.y()) < 1e-6f && std::fabs(a.z()) < 1e-6f,
                 "centripetal term must be -omega^2 r along the arm");
}

bool test_tangential_term_is_perpendicular() {
    const Vector3f a = w3d_lever_acceleration(
        Vector3f::Zero(), Vector3f(0.0f, 0.0f, 3.0f), Vector3f(0.1f, 0.0f, 0.0f));
    return check(std::fabs(a.x()) < 1e-6f && std::fabs(a.y() - 0.3f) < 1e-6f &&
                     std::fabs(a.z()) < 1e-6f,
                 "tangential term must be alpha x r");
}

bool test_no_offset_installs_nothing() {
    W3dLeverArmConfig cfg;
    cfg.model = W3dLeverArmConfig::Model::Exact;
    W3dLeverArm arm(cfg, DT);
    const Vector3f acc(1.0f, -2.0f, 9.7f);
    const Vector3f out = arm.install(acc, Vector3f(0.1f, 0.2f, 0.3f));
    return check(!arm.installs() && (out - acc).norm() == 0.0f,
                 "a zero offset must leave the specific force untouched");
}

bool test_an_arm_along_the_rate_installs_nothing() {
    // Both terms vanish when r is parallel to omega, which is why a purely
    // vertical installation offset is benign under pure yaw.
    W3dLeverArmConfig cfg;
    cfg.offset_body_zu = Vector3f(0.0f, 0.0f, 0.30f);
    W3dLeverArm arm(cfg, DT);
    Rocking body;
    const Vector3f acc_cg(0.0f, 0.0f, 9.81f);
    float worst = 0.0f;
    for (int i = 0; i < 2000; ++i) {
        const float t = float(i) * DT;
        worst = std::max(worst, (arm.install(acc_cg, body.omega(t)) - acc_cg).norm());
    }
    return check(arm.installs() && worst < 1e-6f,
                 "an arm parallel to the rate must install nothing");
}

// --- the derivative --------------------------------------------------------

bool test_backward_difference_tracks_a_known_derivative() {
    Rocking body;
    W3dRateDerivative derivative(DT);
    float worst = 0.0f;
    for (int i = 0; i < 4000; ++i) {
        const float t = float(i) * DT;
        const Vector3f got = derivative.update(body.omega(t));
        if (i < 4) continue;  // second-order start-up
        worst = std::max(worst, (got - body.alpha(t)).norm());
    }
    // The scale to beat is the derivative itself, ~0.63 rad/s^2 here.
    return check(worst < 1e-3f,
                 "causal backward difference must track a smooth truth rate");
}

// --- the two stages --------------------------------------------------------

bool test_exact_model_cancels_the_installation() {
    Rocking body;
    W3dLeverArmConfig cfg;
    // Perpendicular to the rocking axis, so both rotational terms are live.
    cfg.offset_body_zu = Vector3f(0.30f, 0.0f, 0.0f);
    cfg.model = W3dLeverArmConfig::Model::Exact;
    W3dLeverArm arm(cfg, DT);

    const Vector3f acc_cg(0.2f, -0.1f, 9.81f);
    float worst = 0.0f;
    for (int i = 0; i < 4000; ++i) {
        const float t = float(i) * DT;
        const Vector3f omega = body.omega(t);
        const Vector3f installed = arm.install(acc_cg, omega);
        // The oracle runs on the same sample, after sensor corruption would
        // have happened.  Corruption is additive, so pass the term through.
        const Vector3f fused = arm.compensate(installed, omega);
        worst = std::max(worst, (fused - acc_cg).norm());
    }
    return check(worst < 1e-5f && arm.residual_rms_mps2() < 1e-5f,
                 "the oracle model must return the CG specific force");
}

bool test_unmodeled_arm_leaves_the_whole_term() {
    Rocking body;
    W3dLeverArmConfig cfg;
    cfg.offset_body_zu = Vector3f(0.30f, 0.0f, 0.0f);
    W3dLeverArm arm(cfg, DT);

    const Vector3f acc_cg(0.0f, 0.0f, 9.81f);
    std::vector<float> injected;
    for (int i = 0; i < 4000; ++i) {
        const float t = float(i) * DT;
        const Vector3f omega = body.omega(t);
        const Vector3f installed = arm.install(acc_cg, omega);
        const Vector3f fused = arm.compensate(installed, omega);
        injected.push_back((fused - acc_cg).norm());
    }
    const float measured = rms(injected);
    return check(!arm.compensates() && measured > 0.05f &&
                     std::fabs(measured - arm.installed_rms_mps2()) < 1e-3f &&
                     std::fabs(arm.residual_rms_mps2() -
                               arm.installed_rms_mps2()) < 1e-6f,
                 "without a model the residual must equal the installed term");
}

bool test_measured_gyro_model_removes_most_of_the_term() {
    Rocking body;
    W3dLeverArmConfig cfg;
    cfg.offset_body_zu = Vector3f(0.30f, 0.0f, 0.0f);
    cfg.model = W3dLeverArmConfig::Model::MeasuredGyro;
    W3dLeverArm arm(cfg, DT);

    // The deployed harness injects 0.00157 rad/s per-sample gyro white noise
    // and a fixed bias; the model has to work through both.
    std::mt19937 rng(20260829u);
    std::normal_distribution<float> white(0.0f, 0.00157f);
    const Vector3f gyro_bias(4.0e-4f, -3.0e-4f, 6.0e-4f);

    const Vector3f acc_cg(0.0f, 0.0f, 9.81f);
    for (int i = 0; i < 40000; ++i) {
        const float t = float(i) * DT;
        const Vector3f omega = body.omega(t);
        const Vector3f installed = arm.install(acc_cg, omega);
        const Vector3f gyro_meas = omega + gyro_bias +
            Vector3f(white(rng), white(rng), white(rng));
        arm.compensate(installed, gyro_meas);
    }
    const float installed_rms = arm.installed_rms_mps2();
    const float residual_rms = arm.residual_rms_mps2();
    return check(installed_rms > 0.01f && residual_rms < 0.25f * installed_rms,
                 "the gyro-derived model must remove most of the term");
}

bool test_a_too_narrow_derivative_band_is_worse_than_none() {
    // The design point matters.  A 0.2 Hz corner phase-shifts a term whose
    // amplitude is already right, so the "correction" adds error instead of
    // removing it.  The study charts this; the test pins that it is real.
    Rocking body;
    W3dLeverArmConfig cfg;
    cfg.offset_body_zu = Vector3f(0.30f, 0.0f, 0.0f);
    cfg.model = W3dLeverArmConfig::Model::MeasuredGyro;
    cfg.derivative_cutoff_hz = 0.2f;
    W3dLeverArm arm(cfg, DT);

    const Vector3f acc_cg(0.0f, 0.0f, 9.81f);
    for (int i = 0; i < 40000; ++i) {
        const float t = float(i) * DT;
        const Vector3f omega = body.omega(t);
        arm.compensate(arm.install(acc_cg, omega), omega);
    }
    return check(arm.residual_rms_mps2() > arm.installed_rms_mps2(),
                 "an over-narrow derivative band must be a net loss");
}

// --- the environment interface ---------------------------------------------

bool test_env_defaults_to_no_lever_arm() {
    unsetenv("W3D_IMU_LEVER_ARM_M");
    unsetenv("W3D_IMU_LEVER_ARM_MODEL");
    unsetenv("W3D_IMU_LEVER_ARM_CUTOFF_HZ");
    return check(!w3d_lever_arm_config_from_env().has_value() &&
                     w3d_lever_arm_from_env(DT) == nullptr,
                 "an unset environment must leave the historical path alone");
}

bool test_env_parses_offset_and_model() {
    setenv("W3D_IMU_LEVER_ARM_M", "0.1,-0.2,0.3", 1);
    setenv("W3D_IMU_LEVER_ARM_MODEL", "gyro", 1);
    setenv("W3D_IMU_LEVER_ARM_CUTOFF_HZ", "8", 1);
    const auto cfg = w3d_lever_arm_config_from_env();
    const bool ok = cfg.has_value() &&
        std::fabs(cfg->offset_body_zu.x() - 0.1f) < 1e-6f &&
        std::fabs(cfg->offset_body_zu.y() + 0.2f) < 1e-6f &&
        std::fabs(cfg->offset_body_zu.z() - 0.3f) < 1e-6f &&
        cfg->model == W3dLeverArmConfig::Model::MeasuredGyro &&
        std::fabs(cfg->derivative_cutoff_hz - 8.0f) < 1e-6f;
    unsetenv("W3D_IMU_LEVER_ARM_M");
    unsetenv("W3D_IMU_LEVER_ARM_MODEL");
    unsetenv("W3D_IMU_LEVER_ARM_CUTOFF_HZ");
    return check(ok, "the environment must carry the offset and the model");
}

bool test_env_rejects_a_malformed_offset() {
    setenv("W3D_IMU_LEVER_ARM_M", "0.1,0.2", 1);
    bool threw = false;
    try {
        (void)w3d_lever_arm_config_from_env();
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    unsetenv("W3D_IMU_LEVER_ARM_M");
    return check(threw, "a two-component offset must be rejected");
}

bool test_env_rejects_an_unknown_model() {
    setenv("W3D_IMU_LEVER_ARM_M", "0,0,0.3", 1);
    setenv("W3D_IMU_LEVER_ARM_MODEL", "magic", 1);
    bool threw = false;
    try {
        (void)w3d_lever_arm_config_from_env();
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    unsetenv("W3D_IMU_LEVER_ARM_M");
    unsetenv("W3D_IMU_LEVER_ARM_MODEL");
    return check(threw, "an unknown model name must be rejected");
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_centripetal_term_points_inward();
    ok &= test_tangential_term_is_perpendicular();
    ok &= test_no_offset_installs_nothing();
    ok &= test_an_arm_along_the_rate_installs_nothing();
    ok &= test_backward_difference_tracks_a_known_derivative();
    ok &= test_exact_model_cancels_the_installation();
    ok &= test_unmodeled_arm_leaves_the_whole_term();
    ok &= test_measured_gyro_model_removes_most_of_the_term();
    ok &= test_a_too_narrow_derivative_band_is_worse_than_none();
    ok &= test_env_defaults_to_no_lever_arm();
    ok &= test_env_parses_offset_and_model();
    ok &= test_env_rejects_a_malformed_offset();
    ok &= test_env_rejects_an_unknown_model();

    if (!ok) {
        std::cerr << "imu_lever_arm-test FAILED\n";
        return 1;
    }
    std::cout << "imu_lever_arm-test OK\n";
    return 0;
}
