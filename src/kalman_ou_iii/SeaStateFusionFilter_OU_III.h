template<TrackerType trackerT>
class SeaStateFusion_OU_III {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    struct Config {
        bool with_mag = true;

        float mag_delay_sec          = MAG_DELAY_SEC;
        float online_tune_warmup_sec = 10.0f;

        bool  freeze_acc_bias_until_live = true;
        float Racc_warmup_std = 1.2f;

        Eigen::Vector3f sigma_a = Eigen::Vector3f(0.2f, 0.2f, 0.2f);
        Eigen::Vector3f sigma_g = Eigen::Vector3f(0.01f, 0.01f, 0.01f);
        Eigen::Vector3f sigma_m = Eigen::Vector3f(0.3f, 0.3f, 0.3f);

        // Mag-start gate.
        float mag_gravity_align_max_sin   = 0.070f;
        float mag_gravity_align_hold_sec  = 2.0f;
        float mag_gravity_align_lpf_tau   = 1.0f;
        float mag_tilt_fallback_sec       = 30.0f;
        float mag_extreme_gyro_dps        = 45.0f;
        float mag_init_min_mag_norm       = 1e-3f;

        // Real-device mag acquisition.
        //
        // The important rule:
        //
        //   MagAutoTuner must receive tilt from accel/gravity only.
        //   Never feed it tilt extracted from the MEKF quaternion, because that
        //   quaternion contains arbitrary unobservable IMU yaw before mag lock.
        //
        int   mag_min_samples              = 295;
        float mag_min_window_sec           = 0.0f;
        float mag_max_window_sec           = 0.0f;
        float mag_sample_dt_sec            = 1.0f / 200.0f;

        bool  mag_enable_quality_weighting = false;
        float mag_min_effective_weight     = 0.0f;
        float mag_acc_norm_rel_soft        = 0.22f;
        float mag_gyro_soft_dps            = 45.0f;

        // Bootstrap tilt observer for dynamic motion in waves.
        float bootstrap_tilt_obs_acc_tau_sec  = 2.15f;
        float bootstrap_gravity_slow_tau_sec  = 6.0f;
        float bootstrap_gravity_align_max_sin = 0.070f;
        float bootstrap_gravity_hold_sec      = 2.0f;
        float bootstrap_gravity_min_sec       = 6.87f;
        float bootstrap_gravity_timeout_sec   = 15.0f;
        float bootstrap_gravity_norm_frac     = 0.22f;

        bool enable_displacement_detrend = false;
        bool use_custom_displacement_detrend_cfg = false;
        AdaptiveWaveDetrender3D::Config displacement_detrend_cfg{};
    };

    void begin(const Config& cfg) {
        cfg_ = cfg;

        begun_ = true;
        stage_ = Stage::Uninitialized;
        t_ = 0.0f;

        gravity_gate_acc_lpf_.reset();
        mag_gravity_good_sec_ = 0.0f;
        mag_init_eligible_t0_ = NAN;
        last_mag_sample_t_ = NAN;

        mag_ref_set_ = false;

        last_mag_tilt_frame_yaw_rad_ = NAN;
        last_mag_startup_yaw_correction_rad_ = NAN;

        MagAutoTuner::Config mag_cfg;
        mag_cfg.mag_norm_min = cfg_.mag_init_min_mag_norm;
        mag_cfg.min_samples  = cfg_.mag_min_samples;

        mag_cfg.min_window_sec = cfg_.mag_min_window_sec;
        mag_cfg.max_window_sec = cfg_.mag_max_window_sec;
        mag_cfg.sample_dt_sec  = cfg_.mag_sample_dt_sec;

        mag_cfg.gravity_ref = g_std;
        mag_cfg.enable_quality_weighting = cfg_.mag_enable_quality_weighting;
        mag_cfg.min_effective_weight     = cfg_.mag_min_effective_weight;
        mag_cfg.acc_norm_rel_soft        = cfg_.mag_acc_norm_rel_soft;
        mag_cfg.gyro_soft_dps            = cfg_.mag_gyro_soft_dps;

        mag_auto_tuner_.setConfig(mag_cfg);

        resetTiltInit_();

        last_acc_body_ned_.setZero();
        last_gyro_body_ned_.setZero();
        have_last_imu_ = false;

        impl_.setWithMag(cfg_.with_mag);
        impl_.setFreezeAccBiasUntilLive(cfg_.freeze_acc_bias_until_live);
        impl_.setWarmupRaccStd(cfg_.Racc_warmup_std);
        impl_.setMagDelaySec(0.0f); // outer wrapper owns mag delay
        impl_.setOnlineTuneWarmupSec(cfg_.online_tune_warmup_sec);

        impl_.initialize(cfg_.sigma_a, cfg_.sigma_g, cfg_.sigma_m);
        last_impl_startup_stage_ = impl_.getStartupStage();

        impl_.setNominalRaccStd(cfg_.sigma_a);

        displacement_up_m_.setZero();
        displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};

        if (cfg_.enable_displacement_detrend) {
            if (cfg_.use_custom_displacement_detrend_cfg) {
                displacement_detrender_.setConfig(cfg_.displacement_detrend_cfg);
            } else {
                displacement_detrender_.setConfig(
                    seastate::common::defaultDisplacementDetrenderConfig<AdaptiveWaveDetrender3D::Config>(FREQ_GUESS));
            }

            displacement_detrender_.reset(0.0f, 0.0f, 0.0f);
        }
    }

    void update(float dt,
                const Eigen::Vector3f& gyro_body_ned,
                const Eigen::Vector3f& acc_body_ned,
                float tempC = 35.0f)
    {
        if (!begun_) return;
        if (!(dt > 0.0f) || !std::isfinite(dt)) return;

        t_ += dt;

        if (stage_ == Stage::Uninitialized) {
            const bool tilt_ready = seastate::common::runStartupGravityInit(
                gyro_body_ned,
                acc_body_ned,
                dt,
                t_,
                g_std,
                cfg_.bootstrap_tilt_obs_acc_tau_sec,
                cfg_.bootstrap_gravity_slow_tau_sec,
                cfg_.bootstrap_gravity_align_max_sin,
                cfg_.bootstrap_gravity_hold_sec,
                cfg_.bootstrap_gravity_min_sec,
                cfg_.bootstrap_gravity_timeout_sec,
                cfg_.bootstrap_gravity_norm_frac,
                bootstrap_tilt_obs_,
                bootstrap_gravity_slow_lpf_,
                bootstrap_gravity_good_sec_,
                [this](const Eigen::Vector3f& acc_init) {
                    impl_.initialize_from_acc(acc_init);
                });

            if (tilt_ready) {
                stage_ = Stage::Warming;
            }
        }

        last_acc_body_ned_  = acc_body_ned;
        last_gyro_body_ned_ = gyro_body_ned;
        have_last_imu_      = true;

        if (stage_ != Stage::Uninitialized) {
            impl_.updateTime(dt, gyro_body_ned, acc_body_ned, tempC);

            const Eigen::Vector3f acc_gate_lp =
                gravity_gate_acc_lpf_.step(
                    acc_body_ned,
                    dt,
                    cfg_.mag_gravity_align_lpf_tau);

            const float align_sin =
                seastate::common::gravityAlignResidualSin(
                    impl_.mekf().quaternion_boat(),
                    acc_gate_lp);

            const float gyro_dps =
                gyro_body_ned.norm() * 57.295779513f;

            const bool extreme_motion =
                !std::isfinite(gyro_dps) ||
                (gyro_dps > cfg_.mag_extreme_gyro_dps);

            const bool gravity_good_now =
                std::isfinite(align_sin) &&
                (align_sin <= cfg_.mag_gravity_align_max_sin) &&
                !extreme_motion;

            if (gravity_good_now) {
                mag_gravity_good_sec_ += dt;
                if (mag_gravity_good_sec_ > 10.0f) {
                    mag_gravity_good_sec_ = 10.0f;
                }
            } else {
                mag_gravity_good_sec_ =
                    std::max(0.0f, mag_gravity_good_sec_ - 2.0f * dt);
            }

            const Eigen::Vector3f pos_ned_m = impl_.mekf().get_position();

            displacement_up_m_ =
                Eigen::Vector3f(
                    pos_ned_m.x(),
                    pos_ned_m.y(),
                    -pos_ned_m.z());

            if (cfg_.enable_displacement_detrend) {
                const float wave_hz = impl_.getFreqHz();

                const bool ext_freq_valid =
                    isLive() &&
                    std::isfinite(wave_hz) &&
                    (wave_hz >= displacement_detrender_.config().min_wave_freq_hz) &&
                    (wave_hz <= displacement_detrender_.config().max_wave_freq_hz);

                displacement_det_out_ =
                    displacement_detrender_.update(
                        displacement_up_m_,
                        dt,
                        wave_hz,
                        ext_freq_valid);
            } else {
                displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};
                displacement_det_out_.input = displacement_up_m_;
                displacement_det_out_.baseline_slow = Eigen::Vector3f::Zero();
                displacement_det_out_.wave_raw = displacement_up_m_;
                displacement_det_out_.wave_clean = displacement_up_m_;
            }
        }

        const auto cur_stage = impl_.getStartupStage();

        if (cur_stage != last_impl_startup_stage_) {
            if (cur_stage == SeaStateFusionFilter_OU_III<trackerT>::StartupStage::Cold) {
                mag_ref_set_ = false;
                mag_auto_tuner_.reset();

                gravity_gate_acc_lpf_.reset();
                mag_gravity_good_sec_ = 0.0f;
                mag_init_eligible_t0_ = NAN;
                last_mag_sample_t_ = NAN;

                last_mag_tilt_frame_yaw_rad_ = NAN;
                last_mag_startup_yaw_correction_rad_ = NAN;

                if (stage_ != Stage::Live) {
                    stage_ = Stage::Warming;

                    displacement_up_m_.setZero();
                    displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};

                    if (cfg_.enable_displacement_detrend) {
                        displacement_detrender_.reset(0.0f, 0.0f, 0.0f);
                    }
                }
            }

            last_impl_startup_stage_ = cur_stage;
        }

        if (stage_ == Stage::Warming && impl_.isAdaptiveLive()) {
            stage_ = Stage::Live;
        }
    }

    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!begun_ || !cfg_.with_mag) return;
        if (stage_ == Stage::Uninitialized) return;
        if (t_ < cfg_.mag_delay_sec) return;

        if (!std::isfinite(mag_init_eligible_t0_)) {
            mag_init_eligible_t0_ = t_;
        }

        const bool gravity_trusted =
            (mag_gravity_good_sec_ >= cfg_.mag_gravity_align_hold_sec);

        const bool fallback_ok =
            ((t_ - mag_init_eligible_t0_) >= cfg_.mag_tilt_fallback_sec);

        if (!mag_ref_set_) {
            if (!gravity_trusted && !fallback_ok) {
                return;
            }

            if (have_last_imu_) {
                const float dt_mag =
                    (std::isfinite(last_mag_sample_t_) && t_ > last_mag_sample_t_)
                        ? (t_ - last_mag_sample_t_)
                        : cfg_.mag_sample_dt_sec;

                last_mag_sample_t_ = t_;

                Eigen::Quaternionf q_tilt_bw;
                if (!magTiltQuatFromGravityOnly_(q_tilt_bw)) {
                    return;
                }

                if (mag_auto_tuner_.addSampleWithTiltQuatDt(
                        dt_mag,
                        q_tilt_bw,
                        last_acc_body_ned_,
                        last_gyro_body_ned_,
                        mag_body_ned))
                {
                    Eigen::Vector3f mag_world_ref_uT;

                    if (mag_auto_tuner_.getMagWorldRef(mag_world_ref_uT) &&
                        mag_world_ref_uT.allFinite() &&
                        mag_world_ref_uT.norm() > cfg_.mag_init_min_mag_norm)
                    {
                        // This reference was learned in a gravity-only tilt frame.
                        // No MEKF yaw was used to build it.
                        impl_.mekf().set_mag_world_ref(mag_world_ref_uT);

                        const float mag_tilt_yaw_rad =
                            mag_auto_tuner_.getYawGaugeCorrectionRad();

                        if (std::isfinite(mag_tilt_yaw_rad)) {
                            // One-time yaw-gauge lock.
                            //
                            // q_tilt_bw is gravity-only BODY->WORLD.
                            // mag_tilt_yaw_rad is the heading of averaged magnetic
                            // north in that same gravity-only frame.
                            //
                            // Therefore Rz(-mag_tilt_yaw_rad) * q_tilt_bw makes the
                            // learned magnetic north align with +X.
                            const Eigen::Quaternionf q_yaw(
                                Eigen::AngleAxisf(
                                    wrapPi_(-mag_tilt_yaw_rad),
                                    Eigen::Vector3f::UnitZ()));

                            Eigen::Quaternionf q_new = q_yaw * q_tilt_bw;
                            q_new.normalize();

                            if (q_new.coeffs().allFinite()) {
                                impl_.mekf().set_quaternion_boat(q_new);

                                last_mag_tilt_frame_yaw_rad_ =
                                    wrapPi_(mag_tilt_yaw_rad);

                                last_mag_startup_yaw_correction_rad_ =
                                    wrapPi_(-mag_tilt_yaw_rad);
                            }
                        }

                        mag_ref_set_ = true;
                    }
                }
            }
        }

        if (mag_ref_set_) {
            impl_.updateMag(mag_body_ned);
        }
    }

    bool hasMagNorthLock() const noexcept {
        return mag_ref_set_;
    }

    bool isLive() const {
        return stage_ == Stage::Live;
    }

    float freqHz() const {
        return impl_.getFreqHz();
    }

    float waveDirectionDeg() const {
        return impl_.getWaveDirectionDeg();
    }

    Eigen::Vector3f eulerNauticalDeg() const {
        return impl_.getEulerNautical();
    }

    const Eigen::Vector3f& displacementUpMeters() const {
        return displacement_up_m_;
    }

    const AdaptiveWaveDetrender3D::Output& displacementDetrend() const {
        return displacement_det_out_;
    }

    SeaStateFusionFilter_OU_III<trackerT>& raw() {
        return impl_;
    }

    const SeaStateFusionFilter_OU_III<trackerT>& raw() const {
        return impl_;
    }

    int magAcceptedCount() const noexcept {
        return mag_auto_tuner_.acceptedCount();
    }

    int magRejectedCount() const noexcept {
        return mag_auto_tuner_.rejectedCount();
    }

    float magAcceptedWindowSec() const noexcept {
        return mag_auto_tuner_.acceptedWindowSec();
    }

    float magEffectiveWeight() const noexcept {
        return mag_auto_tuner_.effectiveWeight();
    }

    float magTiltFrameYawDeg() const noexcept {
        return std::isfinite(last_mag_tilt_frame_yaw_rad_)
            ? last_mag_tilt_frame_yaw_rad_ * 57.29577951308232f
            : NAN;
    }

    float magStartupYawCorrectionDeg() const noexcept {
        return std::isfinite(last_mag_startup_yaw_correction_rad_)
            ? last_mag_startup_yaw_correction_rad_ * 57.29577951308232f
            : NAN;
    }

private:
    enum class Stage {
        Uninitialized,
        Warming,
        Live
    };

    struct Vec3LPF {
        Eigen::Vector3f state = Eigen::Vector3f::Zero();
        bool initialized = false;

        void reset() {
            state.setZero();
            initialized = false;
        }

        Eigen::Vector3f step(const Eigen::Vector3f& x,
                             float dt,
                             float tau_sec)
        {
            if (!x.allFinite()) return state;

            const float tau = std::max(1.0e-3f, tau_sec);
            const float alpha = 1.0f - std::exp(-dt / tau);

            if (!initialized) {
                state = x;
                initialized = true;
                return state;
            }

            state += alpha * (x - state);
            return state;
        }
    };

    using StartupTiltObserver = seastate::common::StartupTiltObserver;

    void resetTiltInit_() {
        bootstrap_tilt_obs_.reset();
        bootstrap_gravity_slow_lpf_.reset();
        bootstrap_gravity_good_sec_ = 0.0f;
    }

    static float wrapPi_(float a)
    {
        constexpr float PI_F = 3.14159265358979323846f;
        constexpr float TWO_PI_F = 2.0f * PI_F;

        if (!std::isfinite(a)) return NAN;

        while (a > PI_F) {
            a -= TWO_PI_F;
        }

        while (a <= -PI_F) {
            a += TWO_PI_F;
        }

        return a;
    }

    static Eigen::Quaternionf quatFromUnitVectorToUnitVector_(
        const Eigen::Vector3f& from_unit,
        const Eigen::Vector3f& to_unit)
    {
        float d = from_unit.dot(to_unit);
        d = std::max(-1.0f, std::min(1.0f, d));

        Eigen::Vector3f axis = from_unit.cross(to_unit);
        float axis_n = axis.norm();

        if (axis_n < 1.0e-6f) {
            if (d > 0.0f) {
                return Eigen::Quaternionf::Identity();
            }

            Eigen::Vector3f fallback =
                std::fabs(from_unit.x()) < 0.8f
                    ? Eigen::Vector3f::UnitX()
                    : Eigen::Vector3f::UnitY();

            axis = from_unit.cross(fallback);
            axis_n = axis.norm();

            if (!(axis_n > 1.0e-6f) || !axis.allFinite()) {
                return Eigen::Quaternionf::Identity();
            }

            axis /= axis_n;

            Eigen::Quaternionf q(
                Eigen::AngleAxisf(
                    3.14159265358979323846f,
                    axis));

            q.normalize();

            if (!q.coeffs().allFinite()) {
                return Eigen::Quaternionf::Identity();
            }

            return q;
        }

        axis /= axis_n;

        Eigen::Quaternionf q(
            Eigen::AngleAxisf(
                std::acos(d),
                axis));

        q.normalize();

        if (!q.coeffs().allFinite()) {
            return Eigen::Quaternionf::Identity();
        }

        return q;
    }

    static bool tiltOnlyQuatFromAcc_(
        const Eigen::Vector3f& acc_body_ned,
        Eigen::Quaternionf& q_tilt_bw_out)
    {
        q_tilt_bw_out = Eigen::Quaternionf::Identity();

        if (!acc_body_ned.allFinite()) {
            return false;
        }

        const float an = acc_body_ned.norm();

        if (!(an > 1.0e-6f) || !std::isfinite(an)) {
            return false;
        }

        // Your convention:
        //
        //   at rest: acc.z() ≈ -g
        //
        // Therefore gravity-down direction in BODY is:
        //
        //   down_body = -acc / |acc|
        //
        const Eigen::Vector3f down_body =
            (-acc_body_ned / an).normalized();

        const Eigen::Vector3f down_world =
            Eigen::Vector3f::UnitZ();

        q_tilt_bw_out =
            quatFromUnitVectorToUnitVector_(
                down_body,
                down_world);

        q_tilt_bw_out.normalize();

        return q_tilt_bw_out.coeffs().allFinite();
    }

    bool magTiltQuatFromGravityOnly_(
        Eigen::Quaternionf& q_tilt_bw_out) const
    {
        q_tilt_bw_out = Eigen::Quaternionf::Identity();

        Eigen::Vector3f acc_for_tilt = last_acc_body_ned_;

        // Prefer the same slow/LPF accel used by the mag-start gravity gate.
        // This is still gravity-only; it does not contain MEKF yaw.
        if (gravity_gate_acc_lpf_.initialized &&
            gravity_gate_acc_lpf_.state.allFinite() &&
            gravity_gate_acc_lpf_.state.norm() > 1.0e-6f)
        {
            acc_for_tilt = gravity_gate_acc_lpf_.state;
        }

        return tiltOnlyQuatFromAcc_(acc_for_tilt, q_tilt_bw_out);
    }

private:
    Config cfg_{};
    SeaStateFusionFilter_OU_III<trackerT> impl_{false};

    bool begun_ = false;

    Stage stage_ = Stage::Uninitialized;
    float t_ = 0.0f;

    typename SeaStateFusionFilter_OU_III<trackerT>::StartupStage last_impl_startup_stage_ =
        SeaStateFusionFilter_OU_III<trackerT>::StartupStage::Cold;

    Eigen::Vector3f last_acc_body_ned_  = Eigen::Vector3f::Zero();
    Eigen::Vector3f last_gyro_body_ned_ = Eigen::Vector3f::Zero();
    bool have_last_imu_ = false;

    bool mag_ref_set_ = false;
    MagAutoTuner mag_auto_tuner_{};

    float last_mag_sample_t_ = NAN;

    float last_mag_tilt_frame_yaw_rad_ = NAN;
    float last_mag_startup_yaw_correction_rad_ = NAN;

    AdaptiveWaveDetrender3D displacement_detrender_{};
    AdaptiveWaveDetrender3D::Output displacement_det_out_{};
    Eigen::Vector3f displacement_up_m_ = Eigen::Vector3f::Zero();

    Vec3LPF gravity_gate_acc_lpf_{};
    float   mag_gravity_good_sec_ = 0.0f;
    float   mag_init_eligible_t0_ = NAN;

    StartupTiltObserver bootstrap_tilt_obs_{};
    Vec3LPF             bootstrap_gravity_slow_lpf_{};
    float               bootstrap_gravity_good_sec_ = 0.0f;
};
