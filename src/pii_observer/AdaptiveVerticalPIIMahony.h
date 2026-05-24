#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy
*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

#include "ahrs/Mahony_AHRS.h"
#include "pii_observer/AdaptiveVerticalPII.h"

namespace marine_obs {

template<typename T = float,
         bool WithBias = true,
         TrackerType TT = TrackerType::PLL>
class AdaptiveVerticalPIIMahony {
    static_assert(std::is_floating_point<T>::value,
                  "AdaptiveVerticalPIIMahony<T>: T must be a floating-point type.");

public:
    using Core = AdaptiveVerticalPII<T, WithBias, TT>;
    using CoreConfig = typename Core::Config;
    using CoreSnapshot = typename Core::Snapshot;

    struct QuaternionT {
        T w = T(1);
        T x = T(0);
        T y = T(0);
        T z = T(0);
    };

    struct Config {
        CoreConfig core = [] {
            CoreConfig cfg{};
            cfg.observer.r = T(0.150);
            cfg.observer.tau_a = T(0.68);
            cfg.observer.tau_d = T(49.0);
            cfg.observer.kb = T(2.5e-5);
            cfg.observer.lambda_b = T(3.0e-3);
            cfg.observer.bias_limit = T(0.12);
            cfg.observer.a_f_limit = T(50.0);
            cfg.observer.v_limit = T(50.0);
            cfg.observer.p_limit = T(20.0);
            cfg.observer.S_limit = T(200.0);
            cfg.observer.d_limit = T(20.0);
            cfg.adaptation.enabled = true;
            cfg.adaptation.min_confidence = T(0.22);
            cfg.adaptation.f_disp_ref_hz = T(0.12);
            cfg.adaptation.sigma_a_ref = T(0.95);
            cfg.adaptation.input_smooth_tau = T(4.5);
            cfg.adaptation.param_smooth_tau = T(7.5);
            cfg.adaptation.r_freq_exp = T(0.28);
            cfg.adaptation.r_sigma_exp = T(0.02);
            cfg.adaptation.tau_a_freq_exp = T(-0.40);
            cfg.adaptation.tau_a_sigma_exp = T(-0.03);
            cfg.adaptation.tau_d_freq_exp = T(-0.03);
            cfg.adaptation.tau_d_sigma_exp = T(-0.01);
            cfg.adaptation.kb_freq_exp = T(0.02);
            cfg.adaptation.kb_sigma_exp = T(0.08);
            cfg.adaptation.r_min = T(0.145);
            cfg.adaptation.r_max = T(0.225);
            cfg.adaptation.tau_a_min = T(0.50);
            cfg.adaptation.tau_a_max = T(0.90);
            cfg.adaptation.tau_d_min = T(44.0);
            cfg.adaptation.tau_d_max = T(58.0);
            cfg.adaptation.kb_min = T(5e-6);
            cfg.adaptation.kb_max = T(6e-5);
            cfg.auto_schedule_from_accel_freq = true;
            cfg.auto_schedule_period_s = T(0.50);
            cfg.force_enable_adaptation_when_auto_schedule = true;
            cfg.fallback_confidence_floor = T(0.52);
            cfg.fallback_confidence_when_locked = T(0.82);
            cfg.coarse_schedule_blend = T(0.48);
            cfg.coarse_schedule_confidence_floor = T(0.62);
            cfg.accel_freq_tracker =
                detail::make_default_tracker_config<typename Core::AccelFreqTrackerConfig, T>();
            return cfg;
        }();

        T mahony_twoKp = static_cast<T>(1.40);
        T mahony_twoKi = static_cast<T>(0.060);

        T gravity_mps2 = static_cast<T>(9.80665);

        bool use_mag = true;

        bool adapt_mahony_gains = true;

        T mahony_twoKp_calm  = static_cast<T>(0.90);
        T mahony_twoKp_rough = static_cast<T>(0.35);

        T mahony_twoKi_calm  = static_cast<T>(0.025);
        T mahony_twoKi_rough = static_cast<T>(0.010);

        T mahony_sigma_ref     = static_cast<T>(0.18);
        T mahony_norm_err_ref  = static_cast<T>(0.08);
        T mahony_innov_ref     = static_cast<T>(0.12);

        T mahony_gain_smooth_tau_s = static_cast<T>(2.0);

        T mahony_acc_trust_min = static_cast<T>(0.05);
    };

    struct Snapshot {
        CoreSnapshot core{};

        QuaternionT q_world_to_body{};
        QuaternionT q_body_to_world{};

        T roll_deg  = T(0);
        T pitch_deg = T(0);
        T yaw_deg   = T(0);

        T ax_world = T(0);
        T ay_world = T(0);
        T az_world = T(0);

        T vertical_world_accel_up = T(0);

        T mahony_twoKp = static_cast<T>(1.40);
        T mahony_twoKi = static_cast<T>(0.060);
        T gravity_mps2 = static_cast<T>(9.80665);

        T mahony_accel_norm_err = T(0);
        T mahony_accel_innov_norm = T(0);
        T mahony_accel_trust = T(1);
    };

public:
    explicit AdaptiveVerticalPIIMahony(const Config& cfg = Config())
        : core_(cfg.core) {
        configure(cfg);
        reset();
    }

    void configure(const Config& cfg) {
        cfg_ = sanitizeConfig_(cfg);

        core_.configure(cfg_.core);

        gravity_mps2_ = cfg_.gravity_mps2;

        mahony_twoKp_active_ = cfg_.mahony_twoKp;
        mahony_twoKi_active_ = cfg_.mahony_twoKi;

        mahony_AHRS_init(&mahony_, mahony_twoKp_active_, mahony_twoKi_active_);

        last_roll_deg_ = last_pitch_deg_ = last_yaw_deg_ = T(0);
        ax_world_ = ay_world_ = az_world_ = T(0);
        vertical_world_accel_up_ = T(0);

        last_accel_norm_err_ = T(0);
        last_accel_innov_norm_ = T(0);
        last_accel_trust_ = T(1);
    }

    void reset(T p0 = T(0),
               T v0 = T(0),
               T a_f0 = T(0),
               T S0 = T(0),
               T d0 = T(0),
               T b0 = T(0))
    {
        mahony_ = Mahony_AHRS<T>{};

        mahony_twoKp_active_ = cfg_.mahony_twoKp;
        mahony_twoKi_active_ = cfg_.mahony_twoKi;
        mahony_AHRS_init(&mahony_, mahony_twoKp_active_, mahony_twoKi_active_);

        last_roll_deg_ = T(0);
        last_pitch_deg_ = T(0);
        last_yaw_deg_ = T(0);

        ax_world_ = ay_world_ = az_world_ = T(0);
        vertical_world_accel_up_ = T(0);

        last_accel_norm_err_ = T(0);
        last_accel_innov_norm_ = T(0);
        last_accel_trust_ = T(1);

        core_.reset(p0, v0, a_f0, S0, d0, b0);
    }

    void setMahonyGains(T twoKp, T twoKi) {
        if (std::isfinite(twoKp)) {
            mahony_.twoKp = twoKp;
            mahony_twoKp_active_ = twoKp;
            cfg_.mahony_twoKp = twoKp;
        }
        if (std::isfinite(twoKi)) {
            mahony_.twoKi = twoKi;
            mahony_twoKi_active_ = twoKi;
            cfg_.mahony_twoKi = twoKi;
        }
    }

    T updateIMU(T gx_rad_s, T gy_rad_s, T gz_rad_s,
                T ax_mps2, T ay_mps2, T az_mps2,
                T dt_s)
    {
        if (!(std::isfinite(dt_s) && dt_s > T(0))) {
            return core_.displacement();
        }

        adaptMahonyGains_(ax_mps2, ay_mps2, az_mps2, dt_s);

        T pitch_deg = T(0);
        T roll_deg  = T(0);
        T yaw_deg   = T(0);

        mahony_AHRS_update(&mahony_,
                           gx_rad_s, gy_rad_s, gz_rad_s,
                           ax_mps2, ay_mps2, az_mps2,
                           &pitch_deg, &roll_deg, &yaw_deg,
                           dt_s);

        last_pitch_deg_ = pitch_deg;
        last_roll_deg_  = roll_deg;
        last_yaw_deg_   = yaw_deg;

        computeWorldAccelAndVertical_(ax_mps2, ay_mps2, az_mps2);

        return core_.update(vertical_world_accel_up_, dt_s);
    }

    T updateIMUMag(T gx_rad_s, T gy_rad_s, T gz_rad_s,
                   T ax_mps2, T ay_mps2, T az_mps2,
                   T mx, T my, T mz,
                   T dt_s)
    {
        if (!(std::isfinite(dt_s) && dt_s > T(0))) {
            return core_.displacement();
        }

        adaptMahonyGains_(ax_mps2, ay_mps2, az_mps2, dt_s);

        T pitch_deg = T(0);
        T roll_deg  = T(0);
        T yaw_deg   = T(0);

        mahony_AHRS_update_mag(&mahony_,
                               gx_rad_s, gy_rad_s, gz_rad_s,
                               ax_mps2, ay_mps2, az_mps2,
                               mx, my, mz,
                               &pitch_deg, &roll_deg, &yaw_deg,
                               dt_s);

        last_pitch_deg_ = pitch_deg;
        last_roll_deg_  = roll_deg;
        last_yaw_deg_   = yaw_deg;

        computeWorldAccelAndVertical_(ax_mps2, ay_mps2, az_mps2);

        return core_.update(vertical_world_accel_up_, dt_s);
    }

    void updateAdaptationFromDisplacementFrequency(T f_disp_hz,
                                                   T dt_est,
                                                   T confidence = std::numeric_limits<T>::quiet_NaN())
    {
        core_.updateAdaptationFromDisplacementFrequency(f_disp_hz, dt_est, confidence);
    }

    void updateAdaptationExternal(T f_disp_hz,
                                  T sigma_a,
                                  T dt_est,
                                  T confidence = std::numeric_limits<T>::quiet_NaN())
    {
        core_.updateAdaptationExternal(f_disp_hz, sigma_a, dt_est, confidence);
    }

    void updateAdaptationFromAccelFrequencyProxy(T dt_est) {
        core_.updateAdaptationFromAccelFrequencyProxy(dt_est);
    }

    void setAutoScheduleFromAccelFreq(bool on) {
        core_.setAutoScheduleFromAccelFreq(on);
    }

    void setAutoSchedulePeriod(T period_s) {
        core_.setAutoSchedulePeriod(period_s);
    }

    Core& core() { return core_; }
    const Core& core() const { return core_; }

    Mahony_AHRS<T>& mahonyState() { return mahony_; }
    const Mahony_AHRS<T>& mahonyState() const { return mahony_; }

    QuaternionT quaternionWorldToBody() const {
        QuaternionT q;
        q.w = static_cast<T>(mahony_.q0);
        q.x = static_cast<T>(mahony_.q1);
        q.y = static_cast<T>(mahony_.q2);
        q.z = static_cast<T>(mahony_.q3);
        return q;
    }

    QuaternionT quaternionBodyToWorld() const {
        QuaternionT q = quaternionWorldToBody();
        q.x = -q.x;
        q.y = -q.y;
        q.z = -q.z;
        return q;
    }

    T rollDeg() const  { return last_roll_deg_;  }
    T pitchDeg() const { return last_pitch_deg_; }
    T yawDeg() const   { return last_yaw_deg_;   }

    T worldAccelX() const { return ax_world_; }
    T worldAccelY() const { return ay_world_; }
    T worldAccelZUpSpecificForce() const { return az_world_; }

    T verticalWorldAccelUp() const { return vertical_world_accel_up_; }

    T displacement() const { return core_.displacement(); }
    T velocity() const { return core_.velocity(); }
    T accelFiltered() const { return core_.accelFiltered(); }
    T accelSigma() const { return core_.accelSigma(); }
    T accelFrequencyHz() const { return core_.accelFrequencyHz(); }

    T displacementScale(bool smoothed = true) const {
        return core_.displacementScale(smoothed);
    }

    T verticalSpeedEnvelope(bool smoothed = true) const {
        return core_.verticalSpeedEnvelope(smoothed);
    }

    bool envelopeReady() const {
        return core_.envelopeReady();
    }

    T envelopeFrequencyHz() const {
        return core_.envelopeFrequencyHz();
    }

    T envelopeAccelVariance() const {
        return core_.envelopeAccelVariance();
    }

    T envelopeSigmaTarget() const {
        return core_.envelopeSigmaTarget();
    }

    T envelopeSigmaApplied() const {
        return core_.envelopeSigmaApplied();
    }

    T envelopeTauTarget() const {
        return core_.envelopeTauTarget();
    }

    T envelopeTauApplied() const {
        return core_.envelopeTauApplied();
    }

    T mahonyTwoKpActive() const { return mahony_twoKp_active_; }
    T mahonyTwoKiActive() const { return mahony_twoKi_active_; }
    T mahonyAccelNormErr() const { return last_accel_norm_err_; }
    T mahonyAccelInnovNorm() const { return last_accel_innov_norm_; }
    T mahonyAccelTrust() const { return last_accel_trust_; }

    Snapshot snapshot() const {
        Snapshot s;
        s.core = core_.snapshot();
        s.q_world_to_body = quaternionWorldToBody();
        s.q_body_to_world = quaternionBodyToWorld();
        s.roll_deg = last_roll_deg_;
        s.pitch_deg = last_pitch_deg_;
        s.yaw_deg = last_yaw_deg_;
        s.ax_world = ax_world_;
        s.ay_world = ay_world_;
        s.az_world = az_world_;
        s.vertical_world_accel_up = vertical_world_accel_up_;
        s.mahony_twoKp = mahony_twoKp_active_;
        s.mahony_twoKi = mahony_twoKi_active_;
        s.gravity_mps2 = gravity_mps2_;
        s.mahony_accel_norm_err = last_accel_norm_err_;
        s.mahony_accel_innov_norm = last_accel_innov_norm_;
        s.mahony_accel_trust = last_accel_trust_;
        return s;
    }

private:
    static constexpr T eps_() {
        return T(1e-9);
    }

    static T clamp01_(T x) {
        return std::clamp(x, T(0), T(1));
    }

    static T finite_or_default_(T x, T def) {
        return std::isfinite(x) ? x : def;
    }

    static T one_pole_alpha_(T dt, T tau) {
        if (!(std::isfinite(dt) && dt > T(0))) return T(0);
        if (!(std::isfinite(tau) && tau > T(0))) return T(1);
        const T a = dt / (tau + dt);
        return std::clamp(a, T(0), T(1));
    }

    static T safe_rsqrt_or_zero_(T x) {
        if (!(std::isfinite(x) && x > eps_())) return T(0);
        return T(1) / std::sqrt(x);
    }

    static Config sanitizeConfig_(Config cfg) {
        if (!(std::isfinite(cfg.gravity_mps2) && cfg.gravity_mps2 > T(0))) {
            cfg.gravity_mps2 = static_cast<T>(9.80665);
        }

        if (!std::isfinite(cfg.mahony_twoKp)) cfg.mahony_twoKp = static_cast<T>(1.40);
        if (!std::isfinite(cfg.mahony_twoKi)) cfg.mahony_twoKi = static_cast<T>(0.060);

        if (!std::isfinite(cfg.mahony_twoKp_calm))  cfg.mahony_twoKp_calm  = cfg.mahony_twoKp;
        if (!std::isfinite(cfg.mahony_twoKp_rough)) cfg.mahony_twoKp_rough = cfg.mahony_twoKp;
        if (!std::isfinite(cfg.mahony_twoKi_calm))  cfg.mahony_twoKi_calm  = cfg.mahony_twoKi;
        if (!std::isfinite(cfg.mahony_twoKi_rough)) cfg.mahony_twoKi_rough = cfg.mahony_twoKi;

        cfg.mahony_twoKp_calm  = std::max(cfg.mahony_twoKp_calm,  T(0));
        cfg.mahony_twoKp_rough = std::max(cfg.mahony_twoKp_rough, T(0));
        cfg.mahony_twoKi_calm  = std::max(cfg.mahony_twoKi_calm,  T(0));
        cfg.mahony_twoKi_rough = std::max(cfg.mahony_twoKi_rough, T(0));

        cfg.mahony_sigma_ref = std::max(finite_or_default_(cfg.mahony_sigma_ref, T(0.18)), eps_());
        cfg.mahony_norm_err_ref = std::max(finite_or_default_(cfg.mahony_norm_err_ref, T(0.08)), eps_());
        cfg.mahony_innov_ref = std::max(finite_or_default_(cfg.mahony_innov_ref, T(0.12)), eps_());
        cfg.mahony_gain_smooth_tau_s =
            std::max(finite_or_default_(cfg.mahony_gain_smooth_tau_s, T(2.0)), eps_());
        cfg.mahony_acc_trust_min =
            clamp01_(finite_or_default_(cfg.mahony_acc_trust_min, T(0.05)));

        return cfg;
    }

    static void rotateBodyToWorldZUp_(const Mahony_AHRS<T>& m,
                                      T vx, T vy, T vz,
                                      T& ox, T& oy, T& oz)
    {
        const T qw = static_cast<T>(m.q0);
        const T qx = static_cast<T>(m.q1);
        const T qy = static_cast<T>(m.q2);
        const T qz = static_cast<T>(m.q3);

        const T tx = T(2) * (qy * vz - qz * vy);
        const T ty = T(2) * (qz * vx - qx * vz);
        const T tz = T(2) * (qx * vy - qy * vx);

        ox = vx + qw * tx + (qy * tz - qz * ty);
        oy = vy + qw * ty + (qz * tx - qx * tz);
        oz = vz + qw * tz + (qx * ty - qy * tx);
    }

    void computeWorldAccelAndVertical_(T ax_body, T ay_body, T az_body) {
        rotateBodyToWorldZUp_(mahony_, ax_body, ay_body, az_body,
                              ax_world_, ay_world_, az_world_);

        vertical_world_accel_up_ = az_world_ - gravity_mps2_;

        if (!std::isfinite(vertical_world_accel_up_)) {
            vertical_world_accel_up_ = T(0);
        }
    }

    void computeMahonyAccelMetrics_(T ax_body, T ay_body, T az_body,
                                    T& norm_err, T& innov_norm) const
    {
        norm_err = T(0);
        innov_norm = T(0);

        const T a2 = ax_body * ax_body + ay_body * ay_body + az_body * az_body;
        const T inv_a = safe_rsqrt_or_zero_(a2);
        if (!(inv_a > T(0))) {
            return;
        }

        const T ax_n = ax_body * inv_a;
        const T ay_n = ay_body * inv_a;
        const T az_n = az_body * inv_a;

        const T a_norm = T(1) / inv_a;
        norm_err = std::abs(a_norm - gravity_mps2_) / std::max(gravity_mps2_, eps_());

        const T q0 = mahony_.q0;
        const T q1 = mahony_.q1;
        const T q2 = mahony_.q2;
        const T q3 = mahony_.q3;

        const T gx_p = T(2) * (q1 * q3 - q0 * q2);
        const T gy_p = T(2) * (q0 * q1 + q2 * q3);
        const T gz_p = (q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3);

        const T ex = ay_n * gz_p - az_n * gy_p;
        const T ey = az_n * gx_p - ax_n * gz_p;
        const T ez = ax_n * gy_p - ay_n * gx_p;
        innov_norm = std::sqrt(ex * ex + ey * ey + ez * ez);

        if (!std::isfinite(norm_err)) norm_err = T(0);
        if (!std::isfinite(innov_norm)) innov_norm = T(0);
    }

    void adaptMahonyGains_(T ax_body, T ay_body, T az_body, T dt_s)
    {
        if (!cfg_.adapt_mahony_gains) {
            mahony_.init(mahony_twoKp_active_, mahony_twoKi_active_);
            return;
        }

        T norm_err = T(0);
        T innov_norm = T(0);
        computeMahonyAccelMetrics_(ax_body, ay_body, az_body, norm_err, innov_norm);

        last_accel_norm_err_ = norm_err;
        last_accel_innov_norm_ = innov_norm;

        T sigma_a = core_.accelSigma();
        if (!(std::isfinite(sigma_a) && sigma_a > T(0))) {
            sigma_a = cfg_.mahony_sigma_ref;
        }

        const T n = norm_err / std::max(cfg_.mahony_norm_err_ref, eps_());
        const T u = innov_norm / std::max(cfg_.mahony_innov_ref, eps_());
        const T s = sigma_a / std::max(cfg_.mahony_sigma_ref, eps_());

        const T rho_n = T(1) / (T(1) + n * n);
        const T rho_e = T(1) / (T(1) + u * u);
        const T rho_s = T(1) / (T(1) + s * s);

        T trust = rho_n * rho_e * rho_s;
        if (!std::isfinite(trust)) {
            trust = cfg_.mahony_acc_trust_min;
        }
        trust = std::clamp(trust, cfg_.mahony_acc_trust_min, T(1));
        last_accel_trust_ = trust;

        const T kp_cmd =
            cfg_.mahony_twoKp_rough +
            trust * (cfg_.mahony_twoKp_calm - cfg_.mahony_twoKp_rough);

        const T trust_i = trust * trust;
        const T ki_cmd =
            cfg_.mahony_twoKi_rough +
            trust_i * (cfg_.mahony_twoKi_calm - cfg_.mahony_twoKi_rough);

        const T alpha = one_pole_alpha_(dt_s, cfg_.mahony_gain_smooth_tau_s);

        mahony_twoKp_active_ += alpha * (kp_cmd - mahony_twoKp_active_);
        mahony_twoKi_active_ += alpha * (ki_cmd - mahony_twoKi_active_);

        if (!(std::isfinite(mahony_twoKp_active_) && mahony_twoKp_active_ >= T(0))) {
            mahony_twoKp_active_ = cfg_.mahony_twoKp;
        }
        if (!(std::isfinite(mahony_twoKi_active_) && mahony_twoKi_active_ >= T(0))) {
            mahony_twoKi_active_ = cfg_.mahony_twoKi;
        }

        mahony_.init(mahony_twoKp_active_, mahony_twoKi_active_);
    }

private:
    Config cfg_{};

    Core core_;
    Mahony_AHRS<T> mahony_{};

    T gravity_mps2_ = static_cast<T>(9.80665);

    T mahony_twoKp_active_ = static_cast<T>(Mahony_AHRS<T>::twoKpDef);
    T mahony_twoKi_active_ = static_cast<T>(Mahony_AHRS<T>::twoKiDef);

    T last_roll_deg_  = T(0);
    T last_pitch_deg_ = T(0);
    T last_yaw_deg_   = T(0);

    T ax_world_ = T(0);
    T ay_world_ = T(0);
    T az_world_ = T(0);

    T vertical_world_accel_up_ = T(0);

    T last_accel_norm_err_ = T(0);
    T last_accel_innov_norm_ = T(0);
    T last_accel_trust_ = T(1);
};

} // namespace marine_obs
