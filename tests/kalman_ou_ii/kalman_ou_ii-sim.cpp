#include <filesystem>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>

/*
   Copyright (c) 2025-2026  Mikhail Grushinskiy
*/

#define EIGEN_NON_ARDUINO

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "util/W3dSimCommon.h"
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"

using Eigen::Vector3f;
using Eigen::Matrix3f;
using Eigen::Quaternionf;

bool add_noise = true;
bool attitude_only = false;

namespace {

bool env_float(const char* name, float& out)
{
    if (const char* s = std::getenv(name)) {
        out = static_cast<float>(std::atof(s));
        return true;
    }
    return false;
}

bool env_int(const char* name, int& out)
{
    if (const char* s = std::getenv(name)) {
        out = std::atoi(s);
        return true;
    }
    return false;
}

} // namespace

class FusionAdapter_OU_II final : public IW3dFusionAdapter {
public:
    FusionAdapter_OU_II(bool with_mag,
                        const Vector3f& sigma_a_init,
                        const Vector3f& sigma_g,
                        const Vector3f& sigma_m)
        : with_mag_(with_mag)
    {
        cfg_.with_mag = with_mag;
        cfg_.sigma_a = sigma_a_init;
        cfg_.sigma_g = sigma_g;
        cfg_.sigma_m = sigma_m;
        cfg_.mag_delay_sec = MAG_DELAY_SEC;
        cfg_.freeze_acc_bias_until_live = true;
        cfg_.Racc_warmup_std = 0.5f;
        apply_env_overrides();

        fusion_.begin(cfg_);
        auto& filter = fusion_.raw();

        if (attitude_only) {
            filter.enableLinearBlock(false);
            filter.mekf().set_initial_acc_bias(Vector3f::Zero());
            filter.mekf().set_initial_acc_bias_std(0.0f);
            filter.mekf().set_Q_bacc_rw(Vector3f::Zero());
            filter.mekf().set_Racc_std(Vector3f::Constant(0.4f));
        } else {
            filter.enableLinearBlock(true);
            filter.enableTuner(true);
            filter.enableClamp(true);

            float v = 0.0f;

            // Generic OU_* names are accepted for compatibility.
            // OU_II_* names are applied afterward and win if both are set.

            if (env_float("OU_P_FACTOR", v)) {
                filter.setPFactor(v);
            }
            if (env_float("OU_II_P_FACTOR", v)) {
                filter.setPFactor(v);
            }

            if (env_float("OU_R_P0_XY_FACTOR", v)) {
                filter.setR_p0_XYFactor(v);
            }
            if (env_float("OU_II_R_P0_XY_FACTOR", v)) {
                filter.setR_p0_XYFactor(v);
            }

            if (env_float("OU_TAU_COEFF", v)) {
                filter.setTauCoeff(v);
            }
            if (env_float("OU_II_TAU_COEFF", v)) {
                filter.setTauCoeff(v);
            }

            if (env_float("OU_SIGMA_COEFF", v)) {
                filter.setSigmaCoeff(v);
            }
            if (env_float("OU_II_SIGMA_COEFF", v)) {
                filter.setSigmaCoeff(v);
            }

            if (env_float("OU_R_P0_COEFF", v)) {
                filter.setR_p0_Coeff(v);
            }
            if (env_float("OU_II_R_P0_COEFF", v)) {
                filter.setR_p0_Coeff(v);
            }

            if (env_float("OU_R_V0_COEFF", v)) {
                filter.setR_v0_Coeff(v);
            }
            if (env_float("OU_II_R_V0_COEFF", v)) {
                filter.setR_v0_Coeff(v);
            }

            if (env_float("OU_ACC_NOISE_FLOOR_SIGMA", v)) {
                filter.setAccNoiseFloorSigma(v);
            }
            if (env_float("OU_II_ACC_NOISE_FLOOR_SIGMA", v)) {
                filter.setAccNoiseFloorSigma(v);
            }

            if (env_float("OU_ADAPT_TAU_SEC", v)) {
                filter.setAdaptationTimeConstants(v);
            }
            if (env_float("OU_II_ADAPT_TAU_SEC", v)) {
                filter.setAdaptationTimeConstants(v);
            }

            if (env_float("OU_ADAPT_EVERY_SECS", v)) {
                filter.setAdaptationUpdatePeriod(v);
            }
            if (env_float("OU_II_ADAPT_EVERY_SECS", v)) {
                filter.setAdaptationUpdatePeriod(v);
            }

            if (env_float("OU_FREQ_INPUT_CUTOFF_HZ", v)) {
                filter.setFreqInputCutoffHz(v);
            }
            if (env_float("OU_II_FREQ_INPUT_CUTOFF_HZ", v)) {
                filter.setFreqInputCutoffHz(v);
            }

            if (env_float("OU_ACC_BIAS_INIT_STD", v)) {
                filter.mekf().set_initial_acc_bias_std(v);
            }
            if (env_float("OU_II_ACC_BIAS_INIT_STD", v)) {
                filter.mekf().set_initial_acc_bias_std(v);
            }

            Vector3f b = filter.mekf().get_acc_bias();
            bool bias_changed = false;

            if (env_float("OU_ACC_BIAS_INIT_X", v)) {
                b.x() = v;
                bias_changed = true;
            }
            if (env_float("OU_II_ACC_BIAS_INIT_X", v)) {
                b.x() = v;
                bias_changed = true;
            }

            if (env_float("OU_ACC_BIAS_INIT_Y", v)) {
                b.y() = v;
                bias_changed = true;
            }
            if (env_float("OU_II_ACC_BIAS_INIT_Y", v)) {
                b.y() = v;
                bias_changed = true;
            }

            if (env_float("OU_ACC_BIAS_INIT_Z", v)) {
                b.z() = v;
                bias_changed = true;
            }
            if (env_float("OU_II_ACC_BIAS_INIT_Z", v)) {
                b.z() = v;
                bias_changed = true;
            }

            if (bias_changed) {
                filter.mekf().set_initial_acc_bias(b);
            }
        }
    }

    void apply_env_overrides() {
        float vf = 0.0f;
        int vi = 0;

        if (env_float("SF_MAG_DELAY_SEC", vf)) cfg_.mag_delay_sec = vf;
        if (env_float("SF_MAG_GRAV_ALIGN_MAX_SIN", vf)) cfg_.mag_gravity_align_max_sin = vf;
        if (env_float("SF_MAG_GRAV_ALIGN_HOLD_SEC", vf)) cfg_.mag_gravity_align_hold_sec = vf;
        if (env_float("SF_MAG_GRAV_ALIGN_LPF_TAU", vf)) cfg_.mag_gravity_align_lpf_tau = vf;
        if (env_float("SF_MAG_TILT_FALLBACK_SEC", vf)) cfg_.mag_tilt_fallback_sec = vf;
        if (env_float("SF_MAG_EXTREME_GYRO_DPS", vf)) cfg_.mag_extreme_gyro_dps = vf;
        if (env_float("SF_MAG_INIT_MIN_MAG_NORM", vf)) cfg_.mag_init_min_mag_norm = vf;
        if (env_int("SF_MAG_MIN_SAMPLES", vi)) cfg_.mag_min_samples = vi;

        if (env_float("SF_RACC_WARMUP_STD", vf)) cfg_.Racc_warmup_std = vf;
        if (env_float("SF_ONLINE_TUNE_WARMUP_SEC", vf)) cfg_.online_tune_warmup_sec = vf;

        if (env_float("SF_BOOT_TILT_ACC_TAU", vf)) cfg_.bootstrap_tilt_obs_acc_tau_sec = vf;
        if (env_float("SF_BOOT_GRAV_SLOW_TAU", vf)) cfg_.bootstrap_gravity_slow_tau_sec = vf;
        if (env_float("SF_BOOT_GRAV_ALIGN_MAX_SIN", vf)) cfg_.bootstrap_gravity_align_max_sin = vf;
        if (env_float("SF_BOOT_GRAV_HOLD_SEC", vf)) cfg_.bootstrap_gravity_hold_sec = vf;
        if (env_float("SF_BOOT_GRAV_MIN_SEC", vf)) cfg_.bootstrap_gravity_min_sec = vf;
        if (env_float("SF_BOOT_GRAV_TIMEOUT_SEC", vf)) cfg_.bootstrap_gravity_timeout_sec = vf;
        if (env_float("SF_BOOT_GRAV_NORM_FRAC", vf)) cfg_.bootstrap_gravity_norm_frac = vf;
    }

    void updateMag(const Vector3f& mag_body_ned) override {
        fusion_.updateMag(mag_body_ned);
    }

    void update(float dt,
                const Vector3f& gyr_meas_ned,
                const Vector3f& acc_meas_ned,
                float temperature_c) override
    {
        fusion_.update(dt, gyr_meas_ned, acc_meas_ned, temperature_c);
    }

    FilterSnapshot snapshot() const override {
        const auto& filter = fusion_.raw();
        const auto& d = filter.dir();

        FilterSnapshot s;
        s.disp_est_zu = ned_to_zu(filter.mekf().get_position());
        s.vel_est_zu  = ned_to_zu(filter.mekf().get_velocity());
        s.acc_est_zu  = ned_to_zu(filter.mekf().get_world_accel());

// Filter attitude is BODY->WORLD in NED.
//
// In IMU-only mode, after mag lock this is BODY->WORLD in the learned
// magnetic-NED frame, not true-north NED.
//
// Do not apply WMM/declination correction here. A real IMU does not know true
// north unless an external declination/location model is explicitly supplied.
const Quaternionf q_bw_ned = filter.mekf().quaternion_boat().normalized();

float roll_deg  = 0.0f;
float pitch_deg = 0.0f;
float yaw_deg   = 0.0f;
quat_to_euler_nautical(q_bw_ned, roll_deg, pitch_deg, yaw_deg);

s.euler_nautical_deg = Vector3f(roll_deg, pitch_deg, wrapDeg(yaw_deg));
       
        s.acc_bias_est_ned    = filter.mekf().get_acc_bias();
        s.gyro_bias_est_ned   = filter.mekf().gyroscope_bias();
        s.mag_bias_est_ned_uT = get_mag_bias_est_uT(filter.mekf());

        s.tau_target     = filter.getTauTarget();
        s.sigma_target   = filter.getSigmaTarget();
        s.tuning_target  = p0_s_from_sigma_tau(s.sigma_target, s.tau_target);

        s.tau_applied    = filter.getTauApplied();
        s.sigma_applied  = filter.getSigmaApplied();
        s.tuning_applied = p0_s_from_sigma_tau(s.sigma_applied, s.tau_applied);

        s.freq_hz             = filter.getFreqHz();
        s.period_sec          = filter.getPeriodSec();
        s.accel_variance      = filter.getAccelVariance();
        s.displacement_scale_m = filter.getDisplacementScale();
        s.velocity_scale_mps   = filter.getVerticalSpeedEnvelopeMps(true);

        s.direction.phase = d.getPhase();
        s.direction.direction_deg = d.getDirectionDegrees();
        s.direction.direction_deg_generator_signed = dirDegGeneratorSignedFromVec(d.getDirection());
        s.direction.uncertainty_deg = d.getDirectionUncertaintyDegrees();
        s.direction.confidence = d.getLastStableConfidence();
        s.direction.amplitude = d.getAmplitude();
        s.direction.direction_vec = d.getDirection();
        s.direction.filtered_signal = d.getFilteredSignal();

        constexpr float CONF_THRESH = 20.0f;
        constexpr float AMP_THRESH  = 0.08f;
        if (s.direction.confidence > CONF_THRESH && s.direction.amplitude > AMP_THRESH) {
            s.direction.sign = filter.getDirSignState();
            s.direction.sign_num =
                (s.direction.sign == FORWARD) ? 1 :
                (s.direction.sign == BACKWARD ? -1 : 0);
        }

        return s;
    }

private:
    bool with_mag_ = true;
    using Fusion = SeaStateFusion_OU_II<TrackerType::KALMANF>;
    mutable Fusion fusion_;
    Fusion::Config cfg_{};
};

static constexpr W3dFailureLimits FAIL_LIMITS{
    .err_limit_percent_z_jonswap   = 9.80f,
    .err_limit_percent_z_pmstokes  = 9.10f,
    .err_limit_yaw_deg             = 2.62f,
    .err_limit_percent_3d_jonswap  = 34.5f,
    .err_limit_percent_3d_pmstokes = 35.0f,
    .acc_z_bias_percent            = 12.0f,
    .bias_3d_percent               = 300.0f,
};

static constexpr W3dSummaryLabels SUMMARY_LABELS{
    .target = "p0_S_target",
    .applied = "p0_S_applied",
};

static void process_wave_file_for_tracker(const std::string& filename, float dt, bool with_mag)
{
    constexpr float MAG_ODR_HZ = 25.0f;
    auto result = process_wave_file_for_tracker<FusionAdapter_OU_II>(
        filename, dt, with_mag, add_noise, MAG_ODR_HZ,
        "_fusion_ou2", "_fusion_ou2_nomag");

    if (!result) return;
    print_summary_and_fail_if_needed(*result, dt, FAIL_LIMITS, SUMMARY_LABELS);
}

int main(int argc, char* argv[]) {
    const float dt = 1.0f / 200.0f;
    bool with_mag = true;
    add_noise = true;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--nomag") {
            with_mag = false;
        } else if (arg == "--no-noise") {
            add_noise = false;
        }
    }

    std::cout << "Simulation starting with_mag=" << (with_mag ? "true" : "false")
              << ", mag_delay=" << MAG_DELAY_SEC
              << " sec, noise=" << (add_noise ? "true" : "false")
              << "\n";

    const auto files = collect_wave_data_files(".");

    for (const auto& fname : files) {
        process_wave_file_for_tracker(fname, dt, with_mag);
    }

    if (std::getenv("W3D_COLLECT_ALL_GATES") && w3d_any_quality_gate_failed()) return 1;
    return 0;
}
