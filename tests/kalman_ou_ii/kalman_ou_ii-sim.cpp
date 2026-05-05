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
            if (const char* s = std::getenv("OU_P_FACTOR")) filter.setPFactor(std::atof(s));
            if (const char* s = std::getenv("OU_R_P0_XY_FACTOR")) filter.setR_p0_XYFactor(std::atof(s));
            if (const char* s = std::getenv("OU_TAU_COEFF")) filter.setTauCoeff(std::atof(s));
            if (const char* s = std::getenv("OU_SIGMA_COEFF")) filter.setSigmaCoeff(std::atof(s));
            if (const char* s = std::getenv("OU_R_P0_COEFF")) filter.setR_p0_Coeff(std::atof(s));
            if (const char* s = std::getenv("OU_R_V0_COEFF")) filter.setR_v0_Coeff(std::atof(s));
            if (const char* s = std::getenv("OU_ACC_NOISE_FLOOR_SIGMA")) filter.setAccNoiseFloorSigma(std::atof(s));
            if (const char* s = std::getenv("OU_ADAPT_TAU_SEC")) filter.setAdaptationTimeConstants(std::atof(s));
            if (const char* s = std::getenv("OU_ADAPT_EVERY_SECS")) filter.setAdaptationUpdatePeriod(std::atof(s));
            if (const char* s = std::getenv("OU_FREQ_INPUT_CUTOFF_HZ")) filter.setFreqInputCutoffHz(std::atof(s));
if (const char* s = std::getenv("OU_ACC_BIAS_INIT_STD")) {
    filter.mekf().set_initial_acc_bias_std(std::atof(s));
}

if (const char* s = std::getenv("OU_ACC_BIAS_INIT_X")) {
    Vector3f b = filter.mekf().get_acc_bias();
    b.x() = std::atof(s);
    filter.mekf().set_initial_acc_bias(b);
}
if (const char* s = std::getenv("OU_ACC_BIAS_INIT_Y")) {
    Vector3f b = filter.mekf().get_acc_bias();
    b.y() = std::atof(s);
    filter.mekf().set_initial_acc_bias(b);
}
if (const char* s = std::getenv("OU_ACC_BIAS_INIT_Z")) {
    Vector3f b = filter.mekf().get_acc_bias();
    b.z() = std::atof(s);
    filter.mekf().set_initial_acc_bias(b);
}         
        }
    }
    void apply_env_overrides() {
        if (const char* s = std::getenv("SF_MAG_DELAY_SEC")) cfg_.mag_delay_sec = std::atof(s);
        if (const char* s = std::getenv("SF_MAG_GRAV_ALIGN_MAX_SIN")) cfg_.mag_gravity_align_max_sin = std::atof(s);
        if (const char* s = std::getenv("SF_MAG_GRAV_ALIGN_HOLD_SEC")) cfg_.mag_gravity_align_hold_sec = std::atof(s);
        if (const char* s = std::getenv("SF_MAG_GRAV_ALIGN_LPF_TAU")) cfg_.mag_gravity_align_lpf_tau = std::atof(s);
        if (const char* s = std::getenv("SF_MAG_TILT_FALLBACK_SEC")) cfg_.mag_tilt_fallback_sec = std::atof(s);
        if (const char* s = std::getenv("SF_MAG_EXTREME_GYRO_DPS")) cfg_.mag_extreme_gyro_dps = std::atof(s);
        if (const char* s = std::getenv("SF_MAG_INIT_MIN_MAG_NORM")) cfg_.mag_init_min_mag_norm = std::atof(s);
        if (const char* s = std::getenv("SF_MAG_MIN_SAMPLES")) cfg_.mag_min_samples = std::atoi(s);
        if (const char* s = std::getenv("SF_RACC_WARMUP_STD")) cfg_.Racc_warmup_std = std::atof(s);
        if (const char* s = std::getenv("SF_ONLINE_TUNE_WARMUP_SEC")) cfg_.online_tune_warmup_sec = std::atof(s);
        if (const char* s = std::getenv("SF_BOOT_TILT_ACC_TAU")) cfg_.bootstrap_tilt_obs_acc_tau_sec = std::atof(s);
        if (const char* s = std::getenv("SF_BOOT_GRAV_SLOW_TAU")) cfg_.bootstrap_gravity_slow_tau_sec = std::atof(s);
        if (const char* s = std::getenv("SF_BOOT_GRAV_ALIGN_MAX_SIN")) cfg_.bootstrap_gravity_align_max_sin = std::atof(s);
        if (const char* s = std::getenv("SF_BOOT_GRAV_HOLD_SEC")) cfg_.bootstrap_gravity_hold_sec = std::atof(s);
        if (const char* s = std::getenv("SF_BOOT_GRAV_MIN_SEC")) cfg_.bootstrap_gravity_min_sec = std::atof(s);
        if (const char* s = std::getenv("SF_BOOT_GRAV_TIMEOUT_SEC")) cfg_.bootstrap_gravity_timeout_sec = std::atof(s);
        if (const char* s = std::getenv("SF_BOOT_GRAV_NORM_FRAC")) cfg_.bootstrap_gravity_norm_frac = std::atof(s);
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
        // Before mag lock this is just the filter's current world frame.
        // After mag lock it is BODY->WORLD in MAGNETIC-NORTH world.
        Quaternionf q_bw_ned = filter.mekf().quaternion_boat().normalized();

        // Convert magnetic-world attitude into true/world attitude only after
        // magnetic north has actually been learned/locked.
        if (with_mag_ && fusion_.hasMagNorthLock()) {
            // East-positive declination:
            // heading_true = heading_mag + declination
            // Therefore:
            //   C_bw_true = C_true<-mag * C_bw_mag
            const Quaternionf q_mag_to_true_ned =
                quat_from_euler(0.0f, 0.0f, MagSim_WMM::default_declination_deg);
            q_bw_ned = (q_mag_to_true_ned * q_bw_ned).normalized();
        }

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
    .err_limit_percent_z_jonswap   = 9.92f,
    .err_limit_percent_z_pmstokes  = 9.1f,
    .err_limit_yaw_deg             = 2.8f,
    .err_limit_percent_3d_jonswap  = 38.0f,
    .err_limit_percent_3d_pmstokes = 38.0f,
    .acc_z_bias_percent            = 21.0f,
    .bias_3d_percent               = 325.0f,
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
