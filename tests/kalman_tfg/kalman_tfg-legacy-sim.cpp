/*
    Comparison-only simulator for TFG's embedded-friendly legacy r_S law.

    It deliberately uses the SAME current TFG front end, startup, magnetic
    acquisition, hard-iron estimator and sea-state statistics as the default
    simulator.  The only changed estimator choice is

        SpectralMSE -> LegacyCubic

    so comparing this executable with kalman_tfg-sim isolates the new 24/7
    spectral law from the front-end parity changes.  Comparing this executable
    with the pre-PR main results prices the parity changes themselves.
*/

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#define EIGEN_NON_ARDUINO

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "util/W3dSimCommon.h"
#include "kalman_tfg/SeaStateFusionFilter_TFG.h"

using Eigen::Quaternionf;
using Eigen::Vector3f;

namespace {

bool env_float(const char* name, float& out) {
    if (const char* s = std::getenv(name)) {
        out = static_cast<float>(std::atof(s));
        return true;
    }
    return false;
}

class FusionAdapter_TFG_Legacy final : public IW3dFusionAdapter {
public:
    using Fusion = ocean_imu::tfg::SeaStateFusionFilter_TFG<>;

    FusionAdapter_TFG_Legacy(bool with_mag,
                             const Vector3f& sigma_a_init,
                             const Vector3f& sigma_g,
                             const Vector3f& sigma_m) {
        Fusion::Config cfg;
        cfg.with_mag = with_mag;
        cfg.sigma_a = sigma_a_init;
        cfg.gyro_noise_density = sigma_g.x();
        cfg.sigma_m = sigma_m;
        cfg.mag_delay_sec = ocean_imu::tfg::MAG_DELAY_SEC;
        cfg.freeze_acc_bias_until_live = true;
        cfg.Racc_warmup_std = 0.5f;

        if (float v = 0.0f; env_float("SF_SIGMA_A_SCALE", v)) cfg.sigma_a *= v;
        if (float v = 0.0f; env_float("SF_SIGMA_G_SCALE", v)) cfg.gyro_noise_density *= v;
        if (float v = 0.0f; env_float("SF_SIGMA_M_SCALE", v)) cfg.sigma_m *= v;
        if (float v = 0.0f; env_float("SF_GYRO_BIAS_RW_VAR", v)) cfg.gyro_bias_rw_var = v;
        if (float v = 0.0f; env_float("SF_PQ0", v)) cfg.initial_covariance = v;
        if (float v = 0.0f; env_float("SF_MAG_DELAY_SEC", v)) cfg.mag_delay_sec = v;

        fusion_.begin(cfg);
        fusion_.setEmbeddedFriendlyLegacyRSLaw(true);
    }

    void updateMag(const Vector3f& mag_body_ned) override {
        fusion_.updateMag(mag_body_ned);
    }

    void update(float dt,
                const Vector3f& gyr_meas_ned,
                const Vector3f& acc_meas_ned,
                float temperature_c) override {
        fusion_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);
    }

    FilterSnapshot snapshot() const override {
        FilterSnapshot s;
        s.disp_est_zu = ned_to_zu(fusion_.get_position());
        s.vel_est_zu  = ned_to_zu(fusion_.get_velocity());
        s.acc_est_zu  = ned_to_zu(fusion_.get_world_accel());

        const Quaternionf q_bw_ned = fusion_.quaternion().normalized();
        float roll_deg = 0.0f, pitch_deg = 0.0f, yaw_deg = 0.0f;
        quat_to_euler_nautical(q_bw_ned, roll_deg, pitch_deg, yaw_deg);
        s.euler_nautical_deg = Vector3f(roll_deg, pitch_deg, wrapDeg(yaw_deg));

        s.acc_bias_est_ned  = fusion_.mekf().get_acc_bias();
        s.gyro_bias_est_ned = fusion_.mekf().gyroscope_bias();
        s.tau_target     = fusion_.getTauTarget();
        s.sigma_target   = fusion_.getSigmaTarget();
        s.tuning_target  = fusion_.getRSTarget();
        s.tau_applied    = fusion_.getTauApplied();
        s.sigma_applied  = fusion_.getSigmaApplied();
        s.tuning_applied = fusion_.getRSApplied();
        s.wave_period_sec = fusion_.getWavePeriodSec();
        s.accel_variance  = fusion_.getAccelVariance();
        return s;
    }

private:
    Fusion fusion_{};
};

void process_one(const std::string& filename,
                 float dt,
                 bool with_mag,
                 const W3dRandomSeeds& seeds,
                 bool write_timeseries,
                 float validation_window_sec) {
    constexpr float MAG_ODR_HZ = 25.0f;
    auto result = process_wave_file_for_tracker<FusionAdapter_TFG_Legacy>(
        filename, dt, with_mag, true, MAG_ODR_HZ,
        "_fusion_tfg_legacy", "_fusion_tfg_legacy_nomag", seeds, write_timeseries);
    if (!result) return;

    if (validation_window_sec > 0.0f) {
        print_validation_metrics(*result, dt, validation_window_sec, "TFG_LEGACY");
    }

    // These are the pre-PR deterministic TFG bars.  This executable is an
    // ablation/comparison instrument, so callers normally leave
    // W3D_COLLECT_ALL_GATES unset; the bars still make regressions visible in
    // the log without converting expected comparison differences into CI
    // failures.
    static constexpr W3dFailureLimits kOldBars{
        .err_limit_percent_z_jonswap   = 4.812f,
        .err_limit_percent_z_pmstokes  = 4.71f,
        .err_limit_yaw_deg             = 1.352f,
        .err_limit_percent_3d_jonswap  = 20.43f,
        .err_limit_percent_3d_pmstokes = 19.64f,
        .acc_z_bias_percent            = 5.12f,
        .bias_3d_percent               = 164.9f,
    };
    static constexpr W3dSummaryLabels kLabels{ .target = "RS_target",
                                               .applied = "RS_applied" };
    print_summary_and_fail_if_needed(*result, dt, kOldBars, kLabels);
}

}  // namespace

int main(int argc, char* argv[]) {
    const float dt = 1.0f / 200.0f;
    bool with_mag = true;
    std::vector<std::string> requested_files;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--nomag") {
            with_mag = false;
        } else if (arg == "--input") {
            if (++i >= argc) {
                std::cerr << "ERROR: --input requires a CSV path\n";
                return 2;
            }
            requested_files.emplace_back(argv[i]);
        } else if (arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [--nomag] [--input PATH]...\n";
            return 0;
        } else {
            std::cerr << "ERROR: unknown argument: " << arg << "\n";
            return 2;
        }
    }

    W3dRandomSeeds seeds;
    try {
        seeds = w3d_random_seeds_from_env();
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }

    bool write_timeseries = true;
    if (const char* value = std::getenv("W3D_WRITE_TIMESERIES")) {
        write_timeseries = std::string(value) != "0";
    }
    float validation_window_sec = 0.0f;
    env_float("W3D_VALIDATION_WINDOW_SEC", validation_window_sec);

    std::cout << "TFG legacy comparison simulation starting, law=LegacyCubic"
              << ", with_mag=" << (with_mag ? "true" : "false") << "\n";

    const auto files = requested_files.empty()
        ? collect_wave_data_files(".")
        : requested_files;

    try {
        for (const auto& fname : files) {
            process_one(fname, dt, with_mag, seeds, write_timeseries,
                        validation_window_sec);
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }

    if (std::getenv("W3D_COLLECT_ALL_GATES") && w3d_any_quality_gate_failed()) return 1;
    return 0;
}
