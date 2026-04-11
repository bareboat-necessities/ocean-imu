#include <cmath>
#include <type_traits>

class FusionAdapterAdaptivePIIMahony final : public IW3dFusionAdapter {
public:
    using HeaveFilter = marine_obs::AdaptiveVerticalPIIMahony<float, true, TrackerType::PLL>;

    FusionAdapterAdaptivePIIMahony(bool with_mag,
                                   const Vector3f& sigma_a_init,
                                   const Vector3f& sigma_g,
                                   const Vector3f& sigma_m,
                                   const Vector3f& mag_world_a)
        : with_mag_(with_mag),
          filter_(make_config_(with_mag, sigma_a_init, sigma_g, sigma_m, mag_world_a))
    {
    }

    void updateMag(const Vector3f& mag_body_ned) override {
        // Runner supplies body-frame NED at mag ODR.
        // Cache only the latest fresh mag sample.
        last_mag_body_ned_ = mag_body_ned;
        have_mag_ = true;
        mag_fresh_ = true;
    }

    void update(float dt,
                const Vector3f& gyr_meas_ned,
                const Vector3f& acc_meas_ned,
                float temperature_c) override
    {
        (void)temperature_c;

        // Feed the Mahony wrapper in the sim's native body Z-up nautical frame.
        // body NED (North, East, Down) -> body Z-up / ENU-like nautical (East, North, Up)
        const Vector3f gyr_body_zu = ned_to_zu(gyr_meas_ned);
        const Vector3f acc_body_zu = ned_to_zu(acc_meas_ned);

        // Consume magnetometer only once per fresh mag tick.
        // This avoids replaying a stale mag sample at 200 Hz.
        if (with_mag_ && have_mag_ && mag_fresh_) {
            const Vector3f mag_body_zu = ned_to_zu(last_mag_body_ned_);

            filter_.updateIMUMag(
                gyr_body_zu.x(), gyr_body_zu.y(), gyr_body_zu.z(),
                acc_body_zu.x(), acc_body_zu.y(), acc_body_zu.z(),
                mag_body_zu.x(), mag_body_zu.y(), mag_body_zu.z(),
                dt
            );
            mag_fresh_ = false;
        } else {
            filter_.updateIMU(
                gyr_body_zu.x(), gyr_body_zu.y(), gyr_body_zu.z(),
                acc_body_zu.x(), acc_body_zu.y(), acc_body_zu.z(),
                dt
            );
        }
    }

    FilterSnapshot snapshot() const override {
        FilterSnapshot s;

        const auto hs = filter_.snapshot();
        const auto& obs = hs.core.observer;

        // This filter is vertical-only in translation.
        s.disp_est_zu = Vector3f(0.0f, 0.0f, filter_.displacement());
        s.vel_est_zu  = Vector3f(0.0f, 0.0f, filter_.velocity());
        s.acc_est_zu  = Vector3f(0.0f, 0.0f, filter_.accelFiltered());

        // Wrapper exposes world->body quaternion in Z-up world.
        // Because we now feed body Z-up directly, this is directly comparable
        // to sim truth q_wb_zu.
        Quaternionf q_wb_zu(
            hs.q_world_to_body.w,
            hs.q_world_to_body.x,
            hs.q_world_to_body.y,
            hs.q_world_to_body.z
        );

        const bool finite =
            std::isfinite(q_wb_zu.w()) && std::isfinite(q_wb_zu.x()) &&
            std::isfinite(q_wb_zu.y()) && std::isfinite(q_wb_zu.z());

        if (!finite || q_wb_zu.norm() < 1e-6f) {
            q_wb_zu = Quaternionf::Identity();
        } else {
            q_wb_zu.normalize();
        }

        float roll_deg  = 0.0f;
        float pitch_deg = 0.0f;
        float yaw_deg   = 0.0f;
        quat_wb_zu_to_euler_nautical(q_wb_zu, roll_deg, pitch_deg, yaw_deg);

        s.euler_nautical_deg = Vector3f(roll_deg, pitch_deg, yaw_deg);

        s.acc_bias_est_ned    = Vector3f::Zero();
        s.gyro_bias_est_ned   = Vector3f::Zero();
        s.mag_bias_est_ned_uT = Vector3f::Zero();

        const float r_active      = obs.r;
        const float tau_a_active  = obs.tau_a;
        const float tau_d_active  = obs.tau_d;
        const float kb_active     = obs.kb;

        const float sigma_raw     = hs.core.accel_sigma;
        const float sigma_used    = obs.sigma_a_filt;

        const float f_raw_hz      = hs.core.accel_freq_hz;
        const float f_used_hz     = (obs.f_disp_filt_hz > 1e-6f) ? obs.f_disp_filt_hz : f_raw_hz;
        const float omega_used    = (f_used_hz > 1e-6f) ? (2.0f * float(M_PI) * f_used_hz) : NAN;

        s.tau_target     = tau_d_active;
        s.sigma_target   = sigma_raw;
        s.tuning_target  = kb_active;

        s.tau_applied    = tau_a_active;
        s.sigma_applied  = sigma_used;
        s.tuning_applied = r_active;

        s.freq_hz = f_used_hz;
        s.period_sec = (f_used_hz > 1e-6f) ? (1.0f / f_used_hz) : NAN;
        s.accel_variance = hs.core.accel_var;

        if (std::isfinite(omega_used) && omega_used > 1e-6f &&
            std::isfinite(sigma_used) && sigma_used >= 0.0f) {
            s.displacement_scale_m = sigma_used / (omega_used * omega_used);
            s.velocity_scale_mps   = sigma_used / omega_used;
        } else {
            s.displacement_scale_m = NAN;
            s.velocity_scale_mps   = NAN;
        }

        s.direction.phase = NAN;
        s.direction.direction_deg = NAN;
        s.direction.direction_deg_generator_signed = NAN;
        s.direction.uncertainty_deg = NAN;
        s.direction.confidence = NAN;
        s.direction.amplitude = NAN;
        s.direction.direction_vec = Eigen::Vector2f::Zero();
        s.direction.filtered_signal = Eigen::Vector2f::Zero();
        s.direction.sign = UNCERTAIN;
        s.direction.sign_num = 0;

        return s;
    }

private:
    static HeaveFilter::Config make_config_(bool with_mag,
                                            const Vector3f& sigma_a_init,
                                            const Vector3f& sigma_g,
                                            const Vector3f& sigma_m,
                                            const Vector3f& mag_world_a)
    {
        (void)sigma_a_init;
        (void)sigma_g;
        (void)sigma_m;
        (void)mag_world_a;

        HeaveFilter::Config cfg{};

        cfg.core.observer.r          = 0.150f;
        cfg.core.observer.tau_a      = 0.68f;
        cfg.core.observer.tau_d      = 49.0f;
        cfg.core.observer.kb         = 2.5e-5f;
        cfg.core.observer.lambda_b   = 3.0e-3f;
        cfg.core.observer.bias_limit = 0.12f;

        cfg.core.observer.a_f_limit = 50.0f;
        cfg.core.observer.v_limit   = 50.0f;
        cfg.core.observer.p_limit   = 20.0f;
        cfg.core.observer.S_limit   = 200.0f;
        cfg.core.observer.d_limit   = 20.0f;

        cfg.core.adaptation.enabled = true;
        cfg.core.adaptation.min_confidence = 0.22f;

        cfg.core.adaptation.f_disp_ref_hz    = 0.12f;
        cfg.core.adaptation.sigma_a_ref      = 0.95f;
        cfg.core.adaptation.input_smooth_tau = 4.5f;
        cfg.core.adaptation.param_smooth_tau = 7.5f;

        cfg.core.adaptation.r_freq_exp   = 0.28f;
        cfg.core.adaptation.r_sigma_exp  = 0.02f;

        cfg.core.adaptation.tau_a_freq_exp   = -0.40f;
        cfg.core.adaptation.tau_a_sigma_exp  = -0.03f;

        cfg.core.adaptation.tau_d_freq_exp   = -0.03f;
        cfg.core.adaptation.tau_d_sigma_exp  = -0.01f;

        cfg.core.adaptation.kb_freq_exp   = 0.02f;
        cfg.core.adaptation.kb_sigma_exp  = 0.08f;

        cfg.core.adaptation.r_min = 0.145f;
        cfg.core.adaptation.r_max = 0.225f;

        cfg.core.adaptation.tau_a_min = 0.50f;
        cfg.core.adaptation.tau_a_max = 0.90f;

        cfg.core.adaptation.tau_d_min = 44.0f;
        cfg.core.adaptation.tau_d_max = 58.0f;

        cfg.core.adaptation.kb_min = 5e-6f;
        cfg.core.adaptation.kb_max = 6e-5f;

        cfg.core.auto_schedule_from_accel_freq = true;
        cfg.core.auto_schedule_period_s = 0.50f;
        cfg.core.force_enable_adaptation_when_auto_schedule = true;
        cfg.core.fallback_confidence_floor = 0.52f;
        cfg.core.fallback_confidence_when_locked = 0.82f;
        cfg.core.coarse_schedule_blend = 0.48f;
        cfg.core.coarse_schedule_confidence_floor = 0.62f;

        cfg.core.accel_freq_tracker =
            marine_obs::detail::make_default_tracker_config<
                std::remove_cvref_t<decltype(cfg.core.accel_freq_tracker)>, float>();

        cfg.mahony_twoKp = 0.45f;
        cfg.mahony_twoKi = 0.015f;
        cfg.gravity_mps2 = g_std;
        cfg.use_mag = with_mag;

        cfg.adapt_mahony_gains = true;
        cfg.mahony_twoKp_calm  = 0.90f;
        cfg.mahony_twoKp_rough = 0.35f;
        cfg.mahony_twoKi_calm  = 0.025f;
        cfg.mahony_twoKi_rough = 0.010f;
        cfg.mahony_sigma_ref = 0.18f;
        cfg.mahony_norm_err_ref = 0.08f;
        cfg.mahony_gain_smooth_tau_s = 2.0f;
        cfg.mahony_acc_trust_min = 0.05f;

        return cfg;
    }

private:
    bool with_mag_ = true;
    bool have_mag_ = false;
    bool mag_fresh_ = false;

    Vector3f last_mag_body_ned_ = Vector3f::Zero();
    HeaveFilter filter_;
};
