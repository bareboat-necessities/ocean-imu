#include <filesystem>
#include <iostream>
#include <cstdlib>
#include <stdexcept>
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
    // OU-II's own MEKF sensor variances, as multiples of the ones the shared
    // harness hands every family.  The harness builds each as a fixed multiple
    // of the white noise it injects -- 2.8x on accel, 2.0x on gyro, 1.2x on
    // mag -- and those multiples had never been swept for any family.
    // docs/ou-ii-qmekf-variances.md is the sweep.
    //
    // sigma_g is a units correction, not a fit, and it is the same one OU-III
    // carries: the harness multiplies a *per-sample* gyro standard deviation
    // at 200 Hz, but Kalman3D_Wave_OU_II integrates this argument as a noise
    // *density* (Q_AA = Qbase * Ts), so the deployed value overstated the
    // angular random walk by sqrt(200) = 14.1x on top of its own 2x inflation
    // -- 28.3x in std, 800x in variance.  0.05 puts the argument back on the
    // injected density and keeps a sqrt(2) inflation over it, which is where
    // the measured optimum sits.  One error in two places, not two errors.
    //
    // The other two are empirical: accel to 1.4x the injected white (from
    // 2.8x) and mag to 2.4x (from 1.2x).  The mag value is worth a note,
    // because moving it alone is a *loss* here -- 2x costs pitch and 3D in the
    // one-at-a-time sweep -- and only becomes a gain once sigma_g is
    // corrected.  It was adopted from the joint round, not the axis round.
    //
    // Paired over 8 records x 5 seeds against the previous point: roll -12.5%,
    // pitch -11.6%, yaw -1.2%, vertical displacement -1.7%, 3D displacement
    // -23.8%, accelerometer bias -0.7%, gyro bias -36.6%.  tau_applied and
    // sigma_applied are bit-for-bit unchanged, so the OU schedule is untouched.
    //
    // The SF_SIGMA_*_SCALE overrides below multiply these, so a scale of 1
    // reproduces the deployed point and a re-run of the sweep re-centres on it.
    static constexpr float SIGMA_A_RESCALE = 0.5f;   // 2.8x -> 1.4x injected accel white
    static constexpr float SIGMA_G_RESCALE = 0.05f;  // 2.0x sample std -> sqrt(2)x density
    static constexpr float SIGMA_M_RESCALE = 2.0f;   // 1.2x -> 2.4x injected mag white

    FusionAdapter_OU_II(bool with_mag,
                        const Vector3f& sigma_a_init,
                        const Vector3f& sigma_g,
                        const Vector3f& sigma_m)
        : with_mag_(with_mag)
    {
        cfg_.with_mag = with_mag;
        cfg_.sigma_a = sigma_a_init * SIGMA_A_RESCALE;
        cfg_.sigma_g = sigma_g * SIGMA_G_RESCALE;
        cfg_.sigma_m = sigma_m * SIGMA_M_RESCALE;
        cfg_.mag_delay_sec = MAG_DELAY_SEC;
        apply_env_overrides();
        load_fixed_tuning();

        fusion_.begin(cfg_);
        auto& filter = fusion_.raw();

        filter.setPeriodicAwCovarianceSync(load_periodic_aw_cov_sync());

        {
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

            // X and Y are independent; there is deliberately no combined knob,
            // so a sweep that means to move both has to say so twice.
            if (env_float("OU_R_P0_X_FACTOR", v)) {
                filter.setR_p0_XFactor(v);
            }
            if (env_float("OU_II_R_P0_X_FACTOR", v)) {
                filter.setR_p0_XFactor(v);
            }
            if (env_float("OU_R_P0_Y_FACTOR", v)) {
                filter.setR_p0_YFactor(v);
            }
            if (env_float("OU_II_R_P0_Y_FACTOR", v)) {
                filter.setR_p0_YFactor(v);
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

            // Dual pseudo-measurement adaptation law.
            // 0 = Empirical (c_p sigma_aw tau^2, c_v sigma_aw tau),
            // 1 = PhysicalMSE (deployed joint displacement-MSE law).
            if (env_float("OU_II_PSEUDO_LAW", v)) {
                const int law = static_cast<int>(v);
                if (law == 0) {
                    filter.setPseudoLaw(PseudoAdaptationLaw::Empirical);
                } else {
                    filter.setPseudoLaw(PseudoAdaptationLaw::PhysicalMSE);
                }
            }
            // C_P and the channel ratio C_P/C_V of the PhysicalMSE law.
            if (env_float("OU_II_PSEUDO_MSE_COEFF", v)) {
                filter.setPseudoMseCoeff(v);
            }
            if (env_float("OU_II_PSEUDO_MSE_RATIO", v)) {
                filter.setPseudoMseRatio(v);
            }
            if (env_float("OU_II_PSEUDO_RA", v)) {
                filter.setPseudoAccelNoiseDensity(v);
            }

            // Clamps on the two drift-band regularizers.  Exposed because
            // OU-III found that its equivalent floor, not the schedule, was
            // setting the operating point in every low-motion sea, and the
            // only way to notice that is to move the floor and watch whether
            // anything responds.
            {
                float lo = MIN_R_p0_std, hi = MAX_R_p0_std;
                const bool got_lo = env_float("OU_II_R_P0_MIN", lo);
                const bool got_hi = env_float("OU_II_R_P0_MAX", hi);
                if (got_lo || got_hi) {
                    filter.setR_p0_Bounds(lo, hi);
                }
            }
            {
                float lo = MIN_R_v0_std, hi = MAX_R_v0_std;
                const bool got_lo = env_float("OU_II_R_V0_MIN", lo);
                const bool got_hi = env_float("OU_II_R_V0_MAX", hi);
                if (got_lo || got_hi) {
                    filter.setR_v0_Bounds(lo, hi);
                }
            }

            if (env_float("OU_ACC_NOISE_FLOOR_SIGMA", v)) {
                filter.setAccNoiseFloorSigma(v);
            }
            if (env_float("OU_II_ACC_NOISE_FLOOR_SIGMA", v)) {
                filter.setAccNoiseFloorSigma(v);
            }

            // Out-of-band accelerometer guard ahead of the proxy and the
            // MEKF.  Armed by default at the deployed corner; these override
            // it, and a zero cutoff removes it entirely.
            {
                float guard_hz = 0.0f;
                int guard_poles = ACC_VIBRATION_GUARD_POLES_DEFAULT;
                env_int("OU_II_ACC_GUARD_POLES", guard_poles);
                if (env_float("OU_II_ACC_GUARD_HZ", guard_hz)) {
                    filter.setAccelVibrationGuard(guard_hz, guard_poles);
                }
                float racc_gain = 0.0f;
                if (env_float("OU_II_ACC_GUARD_RACC_GAIN", racc_gain)) {
                    filter.setAccelVibrationRaccGain(racc_gain);
                }
                float engage_lo = 0.0f, engage_hi = 0.0f, engage_tau = 0.0f;
                const bool lo_set = env_float("OU_II_ACC_GUARD_ENGAGE_LO", engage_lo);
                const bool hi_set = env_float("OU_II_ACC_GUARD_ENGAGE_HI", engage_hi);
                const bool tau_set = env_float("OU_II_ACC_GUARD_ENGAGE_TAU", engage_tau);
                if (lo_set || hi_set || tau_set) {
                    filter.setAccelVibrationEngagement(
                        lo_set ? engage_lo : -1.0f,
                        hi_set ? engage_hi : -1.0f,
                        tau_set ? engage_tau : -1.0f);
                }
            }

            if (env_float("OU_ADAPT_TAU_SEC", v)) {
                filter.setAdaptationTimeConstants(v);
            }
            if (env_float("OU_II_ADAPT_TAU_SEC", v)) {
                filter.setAdaptationTimeConstants(v);
            }

            // Smoothing horizons of the two drift-correction EMAs, in units
            // of tau_target.
            if (env_float("OU_ADAPT_R_P0_MULT", v)) {
                filter.setR_p0_AdaptMult(v);
            }
            if (env_float("OU_II_ADAPT_R_P0_MULT", v)) {
                filter.setR_p0_AdaptMult(v);
            }

            if (env_float("OU_ADAPT_R_V0_MULT", v)) {
                filter.setR_v0_AdaptMult(v);
            }
            if (env_float("OU_II_ADAPT_R_V0_MULT", v)) {
                filter.setR_v0_AdaptMult(v);
            }

            // Discrepancy threshold that shortens both horizons when the sea
            // state actually moves.  0 keeps the plain proportional horizon.
            if (env_float("OU_ADAPT_R_SLEW_LOG", v)) {
                filter.setR_AdaptSlewLog(v);
            }
            if (env_float("OU_II_ADAPT_R_SLEW_LOG", v)) {
                filter.setR_AdaptSlewLog(v);
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

            // Accelerometer-bias random walk.  The bias competes with the OU
            // acceleration for the low-frequency content, and the wave-band
            // operating point moves the OU corner down toward it, so this is
            // the knob that prices that competition.
            if (env_float("OU_II_ACC_BIAS_RW", v)) {
                filter.mekf().set_Q_bacc_rw(Eigen::Vector3f::Constant(v));
            }

            // Knobs that no longer exist.  tau and the sigma band are
            // wave-band quantities at every instant of the run, and every
            // consumer of a vertical acceleration reads the private Mahony
            // observer.  A stale sweep script must fail here rather than
            // silently reporting the deployed configuration as an ablation.
            for (const char* removed : {"W3D_TUNER_FREQ_SOURCE",
                                        "W3D_TUNING_BAND",
                                        "W3D_FREQ_TRACKER_INPUT"}) {
                if (std::getenv(removed) != nullptr) {
                    throw std::runtime_error(
                        std::string(removed) +
                        " was removed: the tuning frequency is always the wave"
                        " band and the frequency tracker always runs on the"
                        " complementary observer");
                }
            }

            // Wave-band prior used until the period estimator has a value.
            if (env_float("OU_TUNE_FREQ_PRIOR_HZ", v)) {
                filter.setTuneFreqPriorHz(v);
            }

            // sigma_a averaging horizon, in periods of the tuning frequency,
            // and its absolute clamps in seconds.
            if (env_float("OU_SIGMA_VAR_K_PERIODS", v)) {
                filter.setSigmaVarianceKPeriods(v);
            }
            {
                float lo = 0.3f, hi = 60.0f;
                const bool got_lo = env_float("OU_SIGMA_VAR_HORIZON_MIN_S", lo);
                const bool got_hi = env_float("OU_SIGMA_VAR_HORIZON_MAX_S", hi);
                if (got_lo || got_hi) {
                    filter.setSigmaVarianceHorizonBounds(lo, hi);
                }
            }

            // Ablate the wave-period estimator's input away from the
            // complementary-levelled default.  "leveled" restores the older
            // behaviour, which levels with the attitude solution and so closes
            // the tuner coupling.
            if (const char* src = std::getenv("W3D_WAVE_PERIOD_INPUT")) {
                const std::string value = src;
                if (value == "complementary") {
                    filter.setWavePeriodInput(
                        WavePeriodInputSource::Complementary);
                } else if (value == "leveled") {
                    filter.setWavePeriodInput(WavePeriodInputSource::Leveled);
                } else {
                    throw std::runtime_error(
                        "W3D_WAVE_PERIOD_INPUT must be leveled or "
                        "complementary");
                }
            }

            // Gains of that private observer, so the correction corner can be
            // swept against the wave band it must stay below.  It also solves
            // the startup attitude, which is why two_ki now defaults nonzero.
            {
                float two_kp = STARTUP_PROXY_TWO_KP_DEFAULT;
                float two_ki = STARTUP_PROXY_TWO_KI_DEFAULT;
                const bool kp = env_float("W3D_WAVE_PERIOD_MAHONY_KP", two_kp);
                const bool ki = env_float("W3D_WAVE_PERIOD_MAHONY_KI", two_ki);
                if (kp || ki) {
                    filter.setWavePeriodComplementaryGains(two_kp, two_ki);
                }
            }

            // Ablate the tau-scaled pseudo-measurement cadence back to the
            // historical fixed 15 ms one, so the information-rate
            // renormalization can be priced rather than assumed.
            if (const char* raw = std::getenv("W3D_PSEUDO_CADENCE")) {
                const std::string value = raw;
                if (value == "fixed") {
                    filter.setTauScaledPseudoUpdateCadence(false);
                } else if (value == "tau_scaled") {
                    filter.setTauScaledPseudoUpdateCadence(true);
                } else {
                    throw std::runtime_error(
                        "W3D_PSEUDO_CADENCE must be tau_scaled or fixed");
                }
            }
            if (env_float("OU_II_PSEUDO_TAU_RATIO", v)) {
                filter.setPseudoUpdateTauRatio(v);
            }
            {
                float lo = PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT;
                float hi = PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT;
                const bool got_lo = env_float("OU_II_PSEUDO_PERIOD_MIN_S", lo);
                const bool got_hi = env_float("OU_II_PSEUDO_PERIOD_MAX_S", hi);
                if (got_lo || got_hi) {
                    filter.setPseudoUpdatePeriodBounds(lo, hi);
                }
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
        if (env_float("SF_MAG_GRAV_ALIGN_LPF_TAU", vf)) cfg_.mag_gravity_align_world_tau_sec = vf;
        if (env_float("SF_MAG_GRAV_ALIGN_WARMUP_SEC", vf)) cfg_.mag_gravity_align_world_warmup_sec = vf;
        if (env_float("SF_MAG_TILT_FALLBACK_SEC", vf)) cfg_.mag_tilt_fallback_sec = vf;
        if (env_float("SF_MAG_EXTREME_GYRO_DPS", vf)) cfg_.mag_extreme_gyro_dps = vf;
        if (env_float("SF_MAG_INIT_MIN_MAG_NORM", vf)) cfg_.mag_init_min_mag_norm = vf;
        if (env_int("SF_MAG_MIN_SAMPLES", vi)) cfg_.mag_min_samples = vi;
        if (env_float("SF_MAG_MIN_WINDOW_SEC", vf)) cfg_.mag_min_window_sec = vf;

        if (env_float("SF_ONLINE_TUNE_WARMUP_SEC", vf)) cfg_.online_tune_warmup_sec = vf;

        // The MEKF variances the Kalman3D_Wave_OU_II constructor takes; see
        // the SIGMA_*_RESCALE block above and docs/ou-ii-qmekf-variances.md.
        // The three sensor sigmas are swept as scale factors on the deployed
        // point, so a scale of 1 leaves it in place.
        if (env_float("SF_SIGMA_A_SCALE", vf)) cfg_.sigma_a *= vf;
        if (env_float("SF_SIGMA_G_SCALE", vf)) cfg_.sigma_g *= vf;
        if (env_float("SF_SIGMA_M_SCALE", vf)) cfg_.sigma_m *= vf;

        if (env_float("SF_PQ0", vf)) cfg_.Pq0 = vf;
        if (env_float("SF_PB0", vf)) cfg_.Pb0 = vf;
        if (env_float("SF_GYRO_BIAS_RW_VAR", vf)) cfg_.b0 = vf;
        if (env_float("SF_RP0_NOISE_VAR", vf)) cfg_.R_p0_noise = vf;
        if (env_float("SF_RV0_NOISE_VAR", vf)) cfg_.R_v0_noise = vf;

        if (env_float("SF_PROXY_START_MIN_SEC", vf)) cfg_.proxy_startup_min_sec = vf;
        if (env_float("SF_PROXY_START_TIMEOUT_SEC", vf)) cfg_.proxy_startup_timeout_sec = vf;
        if (env_int("SF_ACC_BIAS_UNLOCK_MAG_UPDATES", vi)) cfg_.acc_bias_unlock_mag_updates = vi;
        if (env_float("SF_PROXY_MAG_SETTLE_SEC", vf)) cfg_.proxy_mag_settle_sec = vf;
        if (env_float("SF_MAG_REFINE_START_SEC", vf)) cfg_.mag_refine_start_sec = vf;
        if (env_float("SF_MAG_REFINE_WINDOW_SEC", vf)) cfg_.mag_refine_window_sec = vf;
        if (const char* r = std::getenv("SF_MAG_REFINE")) cfg_.mag_refine_enabled = (std::string(r) != "0");
        // Continuous exogenous hard-iron estimation.  Same names as OU-III's,
        // so a paired study can set one environment and run both families.
        if (const char* h = std::getenv("SF_MAG_CONT_HI")) cfg_.mag_continuous_hard_iron = (std::string(h) != "0");
        if (env_float("SF_MAG_HI_MEMORY_SEC", vf)) cfg_.mag_hi_memory_sec = vf;
        if (env_float("SF_MAG_HI_RIDGE", vf)) cfg_.mag_hi_model_ridge = vf;
        if (env_float("SF_MAG_HI_RIDGE_REL", vf)) cfg_.mag_hi_model_ridge_relative = vf;
        if (env_float("SF_MAG_HI_MIN_INFO", vf)) cfg_.mag_hi_min_information = vf;
        if (env_float("SF_MAG_HI_MIN_WEIGHT", vf)) cfg_.mag_hi_min_effective_weight = vf;
        if (env_float("SF_MAG_HI_MAX_RESID", vf)) cfg_.mag_hi_max_residual_rms_uT = vf;
        if (env_float("SF_MAG_HI_FRACTION", vf)) cfg_.mag_hi_apply_fraction = vf;
        if (env_float("SF_MAG_HI_SLEW_TAU", vf)) cfg_.mag_hi_slew_tau_sec = vf;
        if (env_float("SF_PROXY_TILT_SIGMA", vf)) cfg_.proxy_handoff_tilt_sigma_rad = vf;
        if (env_float("SF_PROXY_YAW_SIGMA", vf)) cfg_.proxy_handoff_yaw_sigma_rad = vf;
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

        // Startup timing marks, to stderr so stdout parsing is untouched.
        if (!reported_lock_ && fusion_.hasMagNorthLock()) {
            reported_lock_ = true;
            std::cerr << "STARTUP first_heading_s=" << fusion_.magNorthLockTimeSec() << "\n";
        }
        if (!reported_live_ && fusion_.isLive()) {
            reported_live_ = true;
            std::cerr << "STARTUP live_s=" << fusion_.liveTimeSec() << "\n";
        }
        if (!reported_refine_ && fusion_.hasRefinedMagReference()) {
            reported_refine_ = true;
            std::cerr << "STARTUP mag_refined_s=" << fusion_.magRefineTimeSec() << "\n";
        }

        auto& filter = fusion_.raw();
        if (fixed_tuning_ && !fixed_tuning_applied_ && filter.isAdaptiveLive()) {
            if (!filter.setFixedTuning(
                    fixed_tau_s_, fixed_sigma_a_, fixed_R_p0_std_, fixed_R_v0_std_))
            {
                throw std::runtime_error("invalid fixed OU-II tuning point");
            }
            fixed_tuning_applied_ = true;
        }
    }

    // Vibration-guard telemetry, so a replay can say whether the guard engaged
    // and how much out-of-band accelerometer content it was seeing.  The shared
    // runner owns the adapter, so end-of-record is the destructor; silent
    // unless a cutoff was configured, which keeps an unguarded run unchanged.
    ~FusionAdapter_OU_II() override {
        const auto& filter = fusion_.raw();
        if (!(filter.accelVibrationGuardCutoffHz() > 0.0f)) return;
        std::cout << "ACC_GUARD cutoff_hz=" << filter.accelVibrationGuardCutoffHz()
                  << " poles=" << filter.accelVibrationGuardPoles()
                  << " engagement=" << filter.accelVibrationGuardEngagement()
                  << " out_of_band_rms_mps2=" << filter.accelVibrationRms()
                  << " delay_sec=" << filter.accelVibrationGuardDelaySec()
                  << " racc_gain=" << filter.accelVibrationRaccGain()
                  << " racc_std_mps2=" << filter.accelVibrationRaccStd().x()
                  << "\n";
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

        // The MEKF has no magnetometer-bias state, so get_mag_bias_est_uT()
        // returns zero here and the harness's "Bias error RMS (mag)" was
        // reporting the injected offset itself rather than anything the filter
        // had done about it.  The wrapper's continuous hard-iron correction is
        // the estimate of that offset -- it is subtracted from every sample the
        // MEKF sees -- so it is what belongs in this slot, and reporting it
        // makes the correction measurable in uT instead of only through yaw.
        s.mag_bias_est_ned_uT = get_mag_bias_est_uT(filter.mekf()) +
                                fusion_.magHardIronBodyUT();

        s.tau_target     = filter.getTauTarget();
        s.sigma_target   = filter.getSigmaTarget();
        s.tuning_target  = p0_s_from_sigma_tau(s.sigma_target, s.tau_target);

        s.tau_applied    = filter.getTauApplied();
        s.sigma_applied  = filter.getSigmaApplied();
        s.tuning_applied = p0_s_from_sigma_tau(s.sigma_applied, s.tau_applied);

        s.freq_hz             = filter.getFreqHz();
        s.wave_period_sec     = filter.getWavePeriodSec();
        s.period_sec          = filter.getPeriodSec();
        s.accel_variance      = filter.getAccelVariance();
        s.displacement_scale_m = filter.getDisplacementScale();
        s.velocity_scale_mps   = filter.getVerticalSpeedEnvelopeMps(true);

        s.direction.phase = d.getPhase();
        s.direction.direction_deg = d.getAxisDegrees();
        s.direction.apparent_to_deg = filter.getApparentWaveDirectionToDeg();
        s.direction.apparent_from_deg = filter.getApparentWaveDirectionFromDeg();
        s.direction.sense_coherence = filter.getDirSenseCoherence();
        s.direction.direction_deg_generator_signed = dirDegGeneratorSignedFromVec(d.getAxis());
        s.direction.uncertainty_deg = d.getAxisUncertaintyDegrees();
        s.direction.confidence = d.getLastStableConfidence();
        s.direction.amplitude = d.getAmplitude();
        s.direction.direction_vec = d.getAxis();
        s.direction.filtered_signal = d.getFilteredSignal();

        constexpr float CONF_THRESH = 20.0f;
        constexpr float AMP_THRESH  = 0.08f;
        if (s.direction.confidence > CONF_THRESH && s.direction.amplitude > AMP_THRESH) {
            s.direction.sign = filter.getDirSignState();
            s.direction.sign_num =
                (s.direction.sign == FORWARD) ? 1 :
                (s.direction.sign == BACKWARD ? -1 : 0);

            // Physical directed propagation vector.  The class above is
            // relative to the estimator's own axis representative; this is not.
            const float travel_x = filter.dir_sign().getDirectedX();
            const float travel_y = filter.dir_sign().getDirectedY();
            if (std::isfinite(travel_x) && std::isfinite(travel_y)) {
                s.direction.travel_vec_boat = Eigen::Vector2f(travel_x, travel_y);
            }
        }

        return s;
    }

private:
    // Selects the a_w covariance-synchronization policy under test.
    // "periodic" (default, and the deployed policy) re-aligns the
    // latent-acceleration marginal with its stationary prior once per
    // adaptation period; "reconfigure" restricts that to discrete
    // reconfiguration events and is the matched ablation.
    static bool load_periodic_aw_cov_sync()
    {
        const char* raw = std::getenv("W3D_AW_COV_SYNC");
        const std::string policy = (raw && *raw) ? raw : "periodic";
        if (policy == "reconfigure") return false;
        if (policy == "periodic") return true;
        throw std::runtime_error(
            "W3D_AW_COV_SYNC must be reconfigure or periodic");
    }

    void load_fixed_tuning()
    {
        const char* raw_mode = std::getenv("W3D_TUNING_MODE");
        const std::string mode = raw_mode ? raw_mode : "adaptive";
        if (mode == "adaptive") return;
        if (!mode.starts_with("fixed")) {
            throw std::runtime_error(
                "W3D_TUNING_MODE must be adaptive or start with fixed");
        }

        fixed_tuning_ =
            env_float("W3D_FIXED_TAU_S", fixed_tau_s_) &&
            env_float("W3D_FIXED_SIGMA_A", fixed_sigma_a_) &&
            env_float("W3D_FIXED_R_P0_STD", fixed_R_p0_std_) &&
            env_float("W3D_FIXED_R_V0_STD", fixed_R_v0_std_);
        if (!fixed_tuning_ ||
            !(std::isfinite(fixed_tau_s_) && fixed_tau_s_ > 0.0f &&
              std::isfinite(fixed_sigma_a_) && fixed_sigma_a_ > 0.0f &&
              std::isfinite(fixed_R_p0_std_) && fixed_R_p0_std_ > 0.0f &&
              std::isfinite(fixed_R_v0_std_) && fixed_R_v0_std_ > 0.0f))
        {
            throw std::runtime_error(
                "fixed OU-II mode requires positive W3D_FIXED_TAU_S, "
                "W3D_FIXED_SIGMA_A, W3D_FIXED_R_P0_STD, and "
                "W3D_FIXED_R_V0_STD");
        }
    }

    bool with_mag_ = true;
    bool fixed_tuning_ = false;
    bool fixed_tuning_applied_ = false;
    float fixed_tau_s_ = NAN;
    float fixed_sigma_a_ = NAN;
    float fixed_R_p0_std_ = NAN;
    float fixed_R_v0_std_ = NAN;
    bool reported_lock_ = false;
    bool reported_live_ = false;
    bool reported_refine_ = false;
    using Fusion = SeaStateFusion_OU_II<TrackerType::KALMANF>;
    mutable Fusion fusion_;
    Fusion::Config cfg_{};
};

// Regression sentinels for the deterministic single-realization protocol, not
// targets.  Each is the worst value the current filter produces across the
// scored records plus about half a percent, rounded up in the last digit the
// channel is quoted in -- a tenth for the percentage channels, a hundredth for
// yaw, which at about two degrees would otherwise be handed several times the
// margin the rule asks for.
//
// That margin is deliberately small because the metrics are deterministic, and
// how deterministic is measured rather than assumed: rebuilding at
// -march=x86-64 instead of the host's native cascadelake moves the gated
// numbers here by at most 4.4e-4 relative (accelerometer bias 3D, jonswap
// H8.5), and yaw by 2.4e-4.  The 6e-6 this comment used to claim still holds
// for the simpler observers -- NLO and the PII observer are within 1.5e-6 --
// but not for a filter carrying this many matrix solves.  Half a percent
// leaves better than a factor of ten on every gate below, which is the check
// to redo before cutting one finer.
//
// Setting one below what the filter currently achieves makes it fail every run
// rather than catching a regression.
//
// Re-derived for the 900 s scoring window: a sentinel fitted to the previous
// 60 s window is not a sentinel for this one, it is just a number the filter
// passes by a wide margin.
//
// bias_3d_percent re-derived again when the r_p0 and r_v0 smoothing horizons
// were shortened from 5 to 3 wave-period-halves.  The binding record moves from
// pmstokes H0.27 to jonswap H1.5, where the horizontal accelerometer bias is
// already unobservable -- the error exceeds the true bias with either horizon
// -- and the displacement error on that record is unchanged.  See
// docs/ou-ema-adaptation-tuning.md.
//
// bias_3d_percent re-derived once more when the frequency tracker moved onto
// the complementary-levelled vertical acceleration.  It is the only sentinel
// that moved: 81.75 -> 84.97 on the same binding record, jonswap H1.5.  Worth
// being explicit that this is a re-derivation and not a relaxation, because
// raising a sentinel to admit one's own change is exactly how these stop
// meaning anything.
//
// The quantity did not get worse; the realization moved.  Replaying that record
// under six seeds gives a mean of 74.2% on the old input and 74.1% on the new
// one, and the per-seed spread reaches 89% either way -- the distribution is
// unchanged and this sentinel samples one point of it.  The record is also the
// one where, as noted above, the horizontal accelerometer bias is unobservable
// to begin with: an error already larger than the true bias is not measuring
// estimation quality, which is why a 3 pp move in it costs no displacement
// accuracy (3D RMS 20.77% -> 20.78%, vertical unchanged to four figures).
// Re-derived once more for the OU-III parity change: the period-scaled sigma
// band, the tau-scaled pseudo-update cadence and the Mahony-proxy startup
// policy.  Five of the seven move.  Three tighten and two loosen, and the two
// that loosen need saying out loud, because raising a sentinel to admit one's
// own change is exactly how these stop meaning anything.
//
// Both loosened limits are realization moves, not quality regressions, and the
// paired ensemble is what establishes that.  Five IMU noise seeds across the
// eight scored records (n = 40 paired records per metric), deployed
// configuration against the pre-change filter:
//
//     3D RMS % of max |disp|    18.997 -> 19.042    +0.24% +/- 0.26%   n.s.
//     accel-bias 3D % of true   82.94  -> 70.14     -15.4% +/- 7.2%    better
//
// So the aggregate accelerometer-bias error improves by 15% while this
// single-realization sentinel gets 8% worse, and the 3D displacement error
// does not move at all.  The binding record for bias_3d_percent is again one
// where the horizontal accelerometer bias is unobservable -- the error exceeds
// the true bias either way -- which is the same property the previous
// re-derivation of this sentinel noted, and the reason a several-percent move
// in it costs no displacement accuracy.
//
// The three that tighten are where the improvement is deterministic enough to
// show up in one realization as well as in the ensemble.
//
// Yaw re-cut to hundredths, 2.2 -> 2.18.  The filter did not move; the tenth
// was worth 1.8 percent of slack against a rule that asks for half of one.
// The six percentage-channel gates are already inside a point of the rule at
// the precision they are quoted in and keep their values.
//
// Re-derived once more for the continuous hard-iron correction, which this
// family now carries on the same defaults as OU-III.  Four of the seven move,
// and they move the same way and for the same reasons they did there.
//
// Yaw halves -- 1.887 to 0.813 deg mean over the eight records, worst 2.161 to
// 1.089 -- because the standing heading error was never a tracking error.  It
// was the hard-iron offset absorbed into the world reference at acquisition,
// which is a gauge, and correcting the stream while moving the reference by a
// delta walks the filter off it.  That gate has to come down with the error or
// it stops being a sentinel: 2.18 -> 1.10, which lands within three hundredths
// of OU-III's own 1.07.
//
// Three go up, and saying so is the point of writing this down, because
// raising a sentinel to admit one's own change is exactly how these stop
// meaning anything:
//
//   acc Z bias %      5.33 -> 5.41   (jonswap H8.5)
//   acc 3D bias %    91.66 -> 93.90  (jonswap H4.0)
//   3D % PM-Stokes   21.02 -> 21.19  (pmstokes H8.5)
//
// The correction walks the heading onto the corrected field during the run,
// and the horizontal accelerometer bias absorbs part of that motion.  It is
// the least observable quantity scored here -- an error above 90 percent of
// the true bias means the error exceeds the thing being estimated, under every
// configuration this family has shipped -- so a two-point move in it is not a
// measurable loss of bias accuracy.  What it is not allowed to be is hidden,
// hence the numbers above.  The displacement channels are flat: vertical mean
// 6.454 -> 6.461 percent of Hs, 3D mean 18.90 -> 19.03 percent, and pitch
// improves (0.289 -> 0.255 deg) for the same reason yaw does.
//
// SF_MAG_CONT_HI=0 is the matched ablation and reproduces the pre-correction
// filter to within 2.6e-4 relative, which is the noise a rebuild of this
// family produces anyway; see docs/quality-gate-regauge.md.  It exceeds the
// yaw gate, as it should -- these are fitted to the filter that ships.  Score
// it with W3D_COLLECT_ALL_GATES.
// Then cut to the rule, which the tenth-quantum was not delivering on the two
// single-digit percentage channels.  A tenth is 1.5 percent of a 6.8 and 1.9
// percent of a 5.4, so rounding a half-percent margin up to one hands back
// three times what the rule asks for; a hundredth costs nothing to write and
// gives it back.  Yaw goes to a thousandth for the same reason -- 1.0949
// rounded to 1.10 is a 0.96 percent bar on a 0.5 percent rule.
//
//   Z %Hs PM-Stokes   6.9  -> 6.85    worst 6.8061, margin 1.38% -> 0.65%
//   yaw deg           1.10 -> 1.095   worst 1.0895, margin 0.96% -> 0.50%
//   acc Z bias %      5.5  -> 5.44    worst 5.4059, margin 1.74% -> 0.63%
//
// The other four were already inside the rule at the precision they are quoted
// in and are left alone.  These three are checked against the drift measured
// for this family rather than against the rule alone: the binding records move
// by 9.4e-5, 1.7e-4 and 3.0e-4 relative between a native and an -march=x86-64
// build, so the smallest of the three margins is still 21 times the spread it
// has to survive.  docs/quality-gate-regauge.md carries that measurement and
// the command that redoes it, which is the check to repeat before cutting any
// of these finer.
//
// Then tightened twice more, with the filter standing still for both -- the
// same two passes OU-III took, since the two families are gated by one rule
// and one script.
//
// First the quantum.  Cutting a tenth to a hundredth, above, fixed the two
// channels where a tenth was worth over a percent, but a fixed absolute step
// still cannot deliver one margin across values from 1 to 94: a hundredth is
// 0.15 percent of a 6.8 and 0.01 percent of a 94.  Quoting every channel to
// four significant figures instead -- so the quantum is a thousandth of the
// value everywhere -- puts all seven between 0.50 and 0.54 percent:
//
//   Z %Hs JONSWAP    6.9  -> 6.899    worst 6.8644   0.52% -> 0.50%
//   Z %Hs PM-Stokes  6.85 -> 6.841    worst 6.8062   0.64% -> 0.51%
//   acc Z bias %     5.44 -> 5.435    worst 5.4073   0.61% -> 0.51%
//   bias 3D %        94.4 -> 94.37    worst 93.8979  0.53% -> 0.50%
//
// Then roll and pitch, which this simulator has measured every run and gated
// never.  What this family has taken on over its last several changes -- the
// OU-III parity work, the Mahony-proxy startup policy, and the continuous
// hard-iron correction, which improved pitch from 0.289 to 0.255 deg mean --
// is largely attitude work, and yaw was the only attitude channel carrying a
// sentinel.  Roll and pitch run 0.2511 to 0.4753 and 0.1801 to 0.3620 deg
// across the eight records, well clear of the bars fitted to them.
//
// They are gated like yaw, on the magnetometer-on protocol only, and for the
// same reason: without it worst-case roll goes to 0.5382 and worst-case pitch
// to 0.7113 deg, both past the bars below.  That is the largest IMU-only
// attitude penalty of the two OU families -- OU-III loses 18 percent of
// worst-case pitch and gains on roll -- and it is a fact about the filter, not
// about the gates, which is why the gates decline to score it.
//
// All nine limits are what tools/ou_regauge_gates.py --family ou_ii prints for
// the filter that ships.  They are checked against this family's own build
// drift rather than against the rule alone, and this family is the noisier of
// the two: the binding records move by 8.0e-6 to 5.6e-4 relative between a
// native and an -march=x86-64 build, and both builds pass all nine.  The
// thinnest margin-to-drift ratio in the set is pitch at 9.3x -- the tightest
// anywhere in the five families, and the first bar to re-measure rather than
// re-cut if a rebuild ever breaches it.  docs/quality-gate-regauge.md carries
// that measurement and the command that redoes it.
//
// Re-derived for the (r_p0, r_v0) coefficient re-fit,
// docs/ou-ii-pseudo-variance-tuning.md.  Every one of the nine still passed on
// its previous value, so this pass is the rule following a filter that moved,
// not a breach being papered over -- but pitch had come down to 0.0001 deg of
// margin against a channel whose measured -march rebuild drift is about 2e-4
// deg, so leaving the set alone would have shipped a bar that a rebuild
// decides, and not the filter.  Reverting both coefficients through
// OU_R_P0_COEFF/OU_R_V0_COEFF puts all nine back at the rule to the digit,
// which is the check that this set moved for the re-fit and for nothing else.
//
// Five tighten and four loosen:
//
//   Z %Hs JONSWAP    6.899 -> 6.865   worst 6.8644 -> 6.8300 (jonswap H0.27)
//   Z %Hs PM-Stokes  6.841 -> 6.848   worst 6.8062 -> 6.8139 (pmstokes H0.27)
//   yaw deg          1.095 -> 1.089   worst 1.0895 -> 1.0833 (jonswap H1.5)
//   roll deg        0.4778 -> 0.4792  worst 0.4753 -> 0.4768 (jonswap H4.0)
//   pitch deg       0.3639 -> 0.3657  worst 0.3620 -> 0.3638 (jonswap H8.5)
//   3D % JONSWAP      21.1 -> 20.92   worst 20.9867 -> 20.8140 (jonswap H1.5)
//   3D % PM-Stokes    21.3 -> 21.03   worst 21.1935 -> 20.9203 (pmstokes H8.5)
//   acc Z bias %     5.435 -> 5.324   worst 5.4073 -> 5.2969 (jonswap H8.5)
//   bias 3D %        94.37 -> 94.47   worst 93.8979 -> 93.9911 (jonswap H4.0)
//
// Both displacement gates come down, by 0.9 and 1.3 percent, which is where a
// change to the translational regularizer is supposed to show up.
//
// Of the four that loosen, three are single-realization moves against an
// ensemble that goes the other way.  Pooled over the eight records at six IMU
// seed triplets the re-fit reads 0.9958 vertical, 0.9796 pitch, 1.0003 roll and
// 1.0006 accelerometer-bias 3D, so pitch improves 2 percent while its
// deterministic worst record moves up 0.0018 deg, and roll is flat within its
// own realization noise.  The vertical one is not noise: it is the small-sea
// end of the position-coefficient trade, and holding the per-record mean
// vertical error where it is -- seven of eight records improve and the eighth
// is 1.0003 -- is what bounded R_p0_coeff at 0.65.  The bias one is the least
// observable quantity scored here, with an error above 90 percent of the true
// bias on the binding record under every configuration this family has shipped.
// Then all nine at once, and all nine downward, when the MEKF sensor
// variances were swept for the first time (docs/ou-ii-qmekf-variances.md).
// The gyro term of that sweep is a units correction -- the harness was
// handing a per-sample standard deviation to an argument the filter
// integrates as a noise density -- so the filter that ships now is a
// materially better one rather than the same one re-drawn, and no channel
// moved the wrong way.  OU-III took the identical correction in the same
// round; the two families share the harness that supplied the argument.
//
//   Z %Hs JONSWAP    6.865  -> 6.776    worst 6.8300 -> 6.7420
//   Z %Hs PM-Stokes  6.848  -> 6.803    worst 6.8139 -> 6.7688
//   yaw deg          1.089  -> 1.074    worst 1.0833 -> 1.0681
//   roll deg         0.4792 -> 0.4352   worst 0.4768 -> 0.4330
//   pitch deg        0.3657 -> 0.279    worst 0.3638 -> 0.2775
//   3D % JONSWAP     20.92  -> 16.54    worst 20.8140 -> 16.4569
//   3D % PM-Stokes   21.03  -> 17.67    worst 20.9203 -> 17.5728
//   acc Z bias %     5.324  -> 4.802    worst 5.2969 -> 4.7776
//   acc 3D bias %    94.47  -> 92.35    worst 93.9911 -> 91.8873
//
// Everything here is what tools/ou_regauge_gates.py prints for the filter as
// it now stands, cut to the same rule as every line above it.
//
// Then all nine again for the deployed physical-MSE pseudo-measurement law
// (docs/ou-ii-dual-mse-adaptation.md).  This is the first regauge in a while
// where bars move in both directions, and the split is the change's own shape
// rather than realization noise:
//
//   Z %Hs JONSWAP    6.776  -> 6.707    worst 6.7420 -> 6.6736 (jonswap H0.27)
//   Z %Hs PM-Stokes  6.803  -> 6.649    worst 6.7688 -> 6.6152 (pmstokes H0.27)
//   yaw deg          1.074  -> 1.073    worst 1.0681 -> 1.0668 (jonswap H0.27)
//   roll deg         0.4352 -> 0.4357   worst 0.4330 -> 0.4335 (jonswap H4.0)
//   pitch deg        0.279  -> 0.2833   worst 0.2775 -> 0.2819 (jonswap H8.5)
//   3D % JONSWAP     16.54  -> 16.84    worst 16.4569 -> 16.7559 (jonswap H8.5)
//   3D % PM-Stokes   17.67  -> 18.27    worst 17.5728 -> 18.1773 (pmstokes H8.5)
//   acc Z bias %     4.802  -> 4.791    worst 4.7776 -> 4.7668 (jonswap H8.5)
//   bias 3D %        92.35  -> 92.33    worst 91.8873 -> 91.8692 (jonswap H4.0)
//
// Both vertical bars come down, by 1.0 and 2.3 percent, and the binding record
// for each is the smallest sea -- which is where the derived law differs most
// from the empirical one it replaced and where the paired multi-seed
// comparison puts its whole gain (-0.063 and -0.104 %Hs at Hs = 0.27 m against
// ties inside their intervals at Hs = 8.5 m).  Yaw and both bias bars come down
// slightly with them.
//
// Three bars loosen, and unlike the usual single-realization moves these are
// real.  Both 3D bars and the pitch bar are bound by the *largest* sea, which
// is the end where the derived law regularizes more loosely than the empirical
// one, and the paired comparison sees the same thing at the same place:
// +0.053 m and +0.026 m of 3D RMS on the two Hs = 8.5 m scenarios, with all
// four small and medium seas improving.  That is the reduced model's stated
// limitation -- one scalar residual-acceleration intensity for all three axes,
// when only the vertical carries gravity leakage -- and it is the accepted cost
// of a law whose primary endpoint improves by 0.0254 +/- 0.0173 %Hs.  Pitch
// moves 1.6 percent on a channel the paired comparison puts at
// +0.0004 +/- 0.0004 deg, i.e. detectable and not meaningful.
//
// Everything here is what tools/ou_regauge_gates.py prints for the filter as
// it now stands, cut to the same rule as every line above it.
//
// Then for the continuous hard-iron re-tune, which cut the estimator's
// absolute ridge floor from 4e-3 to 5e-4.  OU-II takes that calibration
// unchanged from OU-III, for the reason the shared table has always given: a
// sweep that moved one family and not the other would be comparing two
// calibrations rather than two filters.  See docs/continuous-mag-hard-iron.md.
//
// Per record, old ridge -> new: 1.0668 -> 0.6437, 1.0397 -> 0.6533,
// 0.8316 -> 1.0779, 0.3907 -> 0.4910, 0.6896 -> 0.7121, 0.9773 -> 0.5648,
// 0.4867 -> 0.7194, 0.5278 -> 0.4287.  Mean yaw 0.7513 -> 0.6614, five of
// eight records improve.
//
// The yaw bar goes *up*, by one percent, and it is the one number here that
// has to be argued rather than reported.  The worst record changes identity --
// jonswap H0.27 at 1.0668 was the binding one and is now 0.6437, while jonswap
// H4.0 goes 0.8316 -> 1.0779 and takes its place.  What the re-tune does is
// hand back most of the fit on the poorly excited records and a little more on
// the well excited ones; H4.0 is a well excited record on a draw where the
// extra fit is partly aliased distortion, so it pays.  Pooled over five
// magnetometer-calibration draws the same change takes OU-III's yaw down 11
// percent with its own worst record down 10, so this is one record on one
// draw, not the correction over-applying in general.  The bar is re-cut to it
// because a sentinel is fitted to what the filter produces.
//
// Everything else moves in the third digit or better and is re-cut to the
// rule.  All nine are what tools/ou_regauge_gates.py prints for the filter
// that ships.
// Then all nine again when the startup gravity gate moved into the world
// frame (docs/ou-startup-gravity-gate.md).  The gate is what certifies the
// tilt the magnetic reference is framed in and the MEKF is seeded with, and
// in waves it was measuring the sea rather than the levelling error, so it
// only closed by luck: time to a live filter ran 22 s in the calm records and
// 76 to 150 s in the big ones, the worst of them by timeout.  It now closes on
// quality in 22 to 33 s everywhere.
//
// Five records are untouched to four significant figures -- the two calm ones
// are bit-identical, because there the old gate already closed on the first
// quiet stretch -- and the four big-sea records re-mix:
//
//   jonswap  H4.0   roll 0.4345 -> 0.3876, acc 3D bias 91.83 -> 83.77,
//                   yaw 1.0678 -> 1.0419, gyro 3D bias 15.37 -> 14.47
//   pmstokes H8.5   roll 0.2774 -> 0.1667, pitch 0.1772 -> 0.1276,
//                   acc 3D bias 62.42 -> 34.49, gyro 3D bias 19.01 -> 13.70
//   jonswap  H8.5   roll 0.3665 -> 0.1933, acc 3D bias 81.71 -> 64.19,
//                   but pitch 0.2804 -> 0.3280 and 3D 16.7113 -> 17.0920
//
// The two bars that go up are jonswap H8.5's, and its pitch is the one number
// here that is a systematic move rather than a re-draw: paired over six IMU
// seeds it is +12% on every one of them (ratios 1.087 to 1.155).  It is not
// re-cut to hide a regression, it is re-cut because that record is a genuine
// loser in a trade the other three big-sea records win by more -- pooled over
// all eight records and six seeds, pitch is 1.0112 [0.9900, 1.0327], i.e. no
// effect, while vertical error is 0.9966 [0.9940, 0.9992] and gyro bias
// 0.9628 [0.9374, 0.9889], both better at 95%.  Seven of the nine bars come
// down, one of them (the accelerometer bias aggregate) by 9%.
//
// Re-cut for R_p0_x_factor = R_p0_y_factor 1 -> 0.72, the horizontal
// position-regularizer retune that came with the per-axis split; see
// docs/ou-horizontal-anisotropy-per-axis-split.md.  Every bar holds or comes
// down, and the two the retune aims at come down hard: 3D JONSWAP
// 17.18 -> 15.54 and PM-Stokes 17.94 -> 16.62, 9.6% and 7.4% tighter.  Nothing
// here is a loosening.
static constexpr W3dFailureLimits FAIL_LIMITS{
    .err_limit_percent_z_jonswap   = 6.674f,  // was 6.672,  worst 6.6408 (jonswap H0.27)
    .err_limit_percent_z_pmstokes  = 6.605f,  // was 6.605,  worst 6.5717 (pmstokes H0.27)
    .err_limit_yaw_deg             = 1.041f,  // was 1.048,  worst 1.0357 (jonswap H4.0)
    .err_limit_roll_deg            = 0.3866f, // was 0.3896, worst 0.3846 (jonswap H4.0)
    .err_limit_pitch_deg           = 0.3237f, // was 0.3296, worst 0.3221 (jonswap H8.5)
    .err_limit_percent_3d_jonswap  = 15.54f,  // was 17.18,  worst 15.4571 (jonswap H8.5)
    .err_limit_percent_3d_pmstokes = 16.62f,  // was 17.94,  worst 16.5347 (pmstokes H8.5)
    .acc_z_bias_percent            = 4.663f,  // was 4.67,   worst 4.6391 (jonswap H8.5)
    .bias_3d_percent               = 83.97f,  // was 84.2,   worst 83.5458 (jonswap H4.0, accel)
};

static constexpr W3dSummaryLabels SUMMARY_LABELS{
    .target = "p0_S_target",
    .applied = "p0_S_applied",
};

static void process_wave_file_for_tracker(const std::string& filename,
                                          float dt,
                                          bool with_mag,
                                          const W3dRandomSeeds& seeds,
                                          bool write_timeseries,
                                          float validation_window_sec)
{
    constexpr float MAG_ODR_HZ = 25.0f;
    auto result = process_wave_file_for_tracker<FusionAdapter_OU_II>(
        filename, dt, with_mag, add_noise, MAG_ODR_HZ,
        "_fusion_ou2", "_fusion_ou2_nomag", seeds, write_timeseries);

    if (!result) return;
    if (validation_window_sec > 0.0f) {
        print_validation_metrics(*result, dt, validation_window_sec, "OU_II");
    }
    print_summary_and_fail_if_needed(*result, dt, FAIL_LIMITS, SUMMARY_LABELS);
}

int main(int argc, char* argv[]) {
    const float dt = 1.0f / 200.0f;
    bool with_mag = true;
    add_noise = true;
    std::vector<std::string> requested_files;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--nomag") {
            with_mag = false;
        } else if (arg == "--no-noise") {
            add_noise = false;
        } else if (arg == "--input") {
            if (++i >= argc) {
                std::cerr << "ERROR: --input requires a CSV path\n";
                return 2;
            }
            requested_files.emplace_back(argv[i]);
        } else if (arg == "--help") {
            std::cout << "Usage: " << argv[0]
                      << " [--nomag] [--no-noise] [--input PATH]...\n";
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

    std::cout << "Simulation starting with_mag=" << (with_mag ? "true" : "false")
              << ", mag_delay=" << MAG_DELAY_SEC
              << " sec, noise=" << (add_noise ? "true" : "false")
              << "\n";
    std::cout << "RANDOM_SEEDS accel_noise=" << seeds.accel_noise
              << " gyro_noise=" << seeds.gyro_noise
              << " mag_noise=" << seeds.mag_noise
              << " accel_initialization=" << seeds.accel_initialization
              << " gyro_initialization=" << seeds.gyro_initialization
              << " mag_initialization=" << seeds.mag_initialization << "\n";

    const auto files = requested_files.empty()
        ? collect_wave_data_files(".")
        : requested_files;

    try {
        for (const auto& fname : files) {
            process_wave_file_for_tracker(
                fname, dt, with_mag, seeds, write_timeseries,
                validation_window_sec);
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }

    if (std::getenv("W3D_COLLECT_ALL_GATES") && w3d_any_quality_gate_failed()) return 1;
    return 0;
}
