#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Time-varying-gain nonlinear observer, Arduino/Eigen version.

  Based on Bryne/Fossen/Johansen:
  "Nonlinear Observer with Time-Varying Gains for Inertial Navigation
   Aided by Satellite Reference Systems in Dynamic Positioning"

  Compile-time options:
    WithGNSS = true:
      horizontal GNSS/PosRef correction p_xy - p_hat_xy is used.
    WithGNSS = false:
      horizontal GNSS correction is compiled out logically.
      Horizontal TMO states are held at zero because they are unobservable.
      Vertical VVR still works.

    Mag = None:
      no yaw reference. Yaw and z-gyro bias are unobservable.
    Mag = Compass:
      uses compass heading psi_c like the paper.
    Mag = Magnetometer:
      uses calibrated 3D body-frame magnetometer, projected horizontally,
      as a yaw-only magnetic reference.

  Frames:
    State/output frame: NED, +Z down.
    IMU input frame: BODY with NED/marine signs:
      x forward, y starboard, z down.
    q_nb maps BODY -> NED:
      v_n = R_nb * v_b

  IMU input:
    gyro_b_rad_s: BODY angular rate, rad/s.
    specific_force_b_mps2: BODY specific force, m/s^2.
      Level and still in NED convention is roughly [0, 0, -g].
*/

#if defined(ARDUINO)
  #include <Arduino.h>
  #if __has_include(<ArduinoEigenDense.h>)
    #include <ArduinoEigenDense.h>
  #else
    #include <Eigen/Dense>
    #include <Eigen/Geometry>
  #endif
#else
  #include <cmath>
  #include <cstdint>
  #include <Eigen/Dense>
  #include <Eigen/Geometry>
#endif

#include <cstdint>

enum class NloMagType : std::uint8_t {
    None,
    Compass,
    Magnetometer
};

template <bool WithGNSS = true,
          NloMagType Mag = NloMagType::Compass,
          typename R = float>
class TimeVaryingGainNLO {
public:
    using Vec2 = Eigen::Matrix<R, 2, 1>;
    using Vec3 = Eigen::Matrix<R, 3, 1>;
    using Mat2 = Eigen::Matrix<R, 2, 2>;
    using Mat3 = Eigen::Matrix<R, 3, 3>;
    using Quat = Eigen::Quaternion<R>;

    static constexpr bool kWithGNSS = WithGNSS;
    static constexpr NloMagType kMagType = Mag;

    struct GnssXY {
        Vec2 position_n_m = Vec2::Zero();
        R rms_xy_m = R(1);
        bool valid = false;
    };

    struct Aux {
        R compass_heading_rad = R(0);
        bool compass_valid = true;

        Vec3 mag_b = Vec3::Zero();
        bool mag_valid = false;

        GnssXY gnss;
    };

    struct Config {
        R gravity_mps2 = R(9.80665);

        R gyro_bias_limit_rad_s = R(0.20);
        R max_specific_force_mps2 = R(30.0);

        /*
          No-GNSS roll/pitch attitude leak.

          This is used only for quaternion / gyro-bias correction when
          WithGNSS=false. It is not used for xi_dot / VVR dynamics.

          Keep this weak. Too much attitude leak improves roll/pitch but
          damages heave in steep waves.
        */        
        R no_gnss_attitude_sigma_scale = R(0.03);
        R no_gnss_attitude_sigma_limit_rad_s = R(0.0020);

        bool use_time_varying_attitude_gains = true;
        R attitude_gain_tau_s = R(25);
        R attitude_gain_switch_s = R(100);

        R k1_initial = R(20);
        R k2_initial = R(20);
        R kI_initial = R(1);

        R k1_nominal = R(0.55);
        R k2_nominal = R(1.0);
        R kI_nominal = R(0.01);

        R K_p0z_p0z = R(5.4295);
        R K_pz_p0z  = R(2.2396);
        R K_vz_p0z  = R(0.4454);
        R K_xiz_p0z = R(0.0354);

        R K_pp_scalar  = R(0.9513);
        R K_vp_scalar  = R(0.3275);
        R K_xip_scalar = R(0.0354);

        R theta = R(1);

        bool use_time_varying_tmo_gain = true;
        R vartheta0 = R(0.5);
        R vartheta1_a = R(2.0);
        R vartheta1_b = R(1.5);
        R vartheta1_without_gnss = R(0.0);
        R gnss_rms_lpf_tau_s = R(125);

        R vartheta2_tau_s = R(25);
        R vartheta2_switch_s = R(100);

        R p0z_highpass_tau_s = R(600);

        bool use_triad_style_force_injection = true;

        R magnetic_declination_rad = R(0);

        R expected_mag_norm = R(0);
        R max_mag_relative_norm_error = R(0.35);
    };

    explicit TimeVaryingGainNLO(const Config& cfg = Config{})
        : cfg_(cfg)
    {
        reset();
    }

    void reset() {
        t_s_ = R(0);

        q_nb_.setIdentity();
        gyro_bias_b_.setZero();
        accel_bias_b_.setZero();

        p_n_.setZero();
        v_n_.setZero();
        xi_n_.setZero();
        fhat_n_.setZero();

        p0z_hat_ = R(0);
        p0z_hp_state_ = R(0);

        k1_ = cfg_.use_time_varying_attitude_gains ? cfg_.k1_initial : cfg_.k1_nominal;
        k2_ = cfg_.use_time_varying_attitude_gains ? cfg_.k2_initial : cfg_.k2_nominal;
        kI_ = cfg_.use_time_varying_attitude_gains ? cfg_.kI_initial : cfg_.kI_nominal;

        gnss_rms_lpf_ = R(1);
        vartheta2_ = R(1);
        vartheta_ = R(1);

        sigma_b_.setZero();

        lockUnobservableStates_();
    }

    bool initializeFromAccel(const Vec3& specific_force_b_mps2,
                             R yaw_seed_rad = R(0)) {
        return initializeFromAccelCompass(specific_force_b_mps2, yaw_seed_rad);
    }

    bool initializeFromAccelCompass(const Vec3& specific_force_b_mps2,
                                    R heading_rad) {
        const R n = specific_force_b_mps2.norm();
        if (!isFinite(n) || n < R(1e-6)) {
            return false;
        }

        const Vec3 down_b = (-specific_force_b_mps2 / n).normalized();

        Vec3 north_b;
        north_b << std::cos(heading_rad), -std::sin(heading_rad), R(0);

        const bool ok = setAttitudeFromDownAndHorizontalReference(
            down_b,
            north_b,
            nedNorth()
        );

        lockUnobservableStates_();
        return ok;
    }

    bool initializeFromAccelMag(const Vec3& specific_force_b_mps2,
                                const Vec3& mag_b_calibrated) {
        const R fn = specific_force_b_mps2.norm();
        const R mn = mag_b_calibrated.norm();

        if (!isFinite(fn) || !isFinite(mn) || fn < R(1e-6) || mn < R(1e-6)) {
            return false;
        }

        const Vec3 down_b = (-specific_force_b_mps2 / fn).normalized();

        Vec3 m_b = mag_b_calibrated / mn;
        Vec3 mag_h_b = m_b - down_b * down_b.dot(m_b);

        const bool ok = setAttitudeFromDownAndHorizontalReference(
            down_b,
            mag_h_b,
            magneticReferenceNED()
        );

        lockUnobservableStates_();
        return ok;
    }

    void setQuaternionBodyToNED(const Quat& q_nb) {
        q_nb_ = q_nb;
        q_nb_.normalize();
        lockUnobservableStates_();
    }

    void setPositionNED(const Vec3& p_n_m) {
        p_n_ = p_n_m;
        if constexpr (!WithGNSS) {
            p_n_.x() = R(0);
            p_n_.y() = R(0);
        }
    }

    void setVelocityNED(const Vec3& v_n_mps) {
        v_n_ = v_n_mps;
        if constexpr (!WithGNSS) {
            v_n_.x() = R(0);
            v_n_.y() = R(0);
        }
    }

    void setGyroBiasBody(const Vec3& bias_rad_s) {
        gyro_bias_b_ = clampNorm(bias_rad_s, cfg_.gyro_bias_limit_rad_s);
        if constexpr (Mag == NloMagType::None) {
            gyro_bias_b_.z() = R(0);
        }
    }

    void setAccelBiasBody(const Vec3& bias_mps2) {
        accel_bias_b_ = bias_mps2;
    }

    void setMagneticDeclinationRad(R declination_rad) {
        cfg_.magnetic_declination_rad = declination_rad;
    }

    void update(R dt,
                const Vec3& gyro_b_rad_s,
                const Vec3& specific_force_b_mps2,
                const Aux& aux = Aux{}) {
        if (!(dt > R(0)) || !isFinite(dt)) {
            return;
        }

        t_s_ += dt;

        lockUnobservableStates_();

        const Vec3 f_b = specific_force_b_mps2 - accel_bias_b_;

        updateGainSchedules(dt, aux);

        const Mat3 R_nb_before = q_nb_.toRotationMatrix();
        fhat_n_ = R_nb_before * f_b + xi_n_;

        if constexpr (!WithGNSS) {
            xi_n_.x() = R(0);
            xi_n_.y() = R(0);
            fhat_n_ = R_nb_before * f_b + xi_n_;
        }

        /*
          Split the two sigma roles:

            sigma_tmo_b:
              Original paper-style force feedback. This is used for xi_dot /
              TMO/VVR correction dynamics. It must remain the full-force
              version, otherwise the vertical heave loop can run away.

            sigma_b_:
              Public/attitude sigma used for quaternion correction and gyro
              bias learning. In no-GNSS mode, the original full-force feedback
              self-cancels because fhat_n ~= R*f_b. Use a small vertical-only
              force-reference leak for roll/pitch, but do not feed that into
              xi_dot.
        */
        const Vec3 sigma_tmo_b =
            computeSigmaBody(f_b, fhat_n_, aux, R_nb_before);

        sigma_b_ = sigma_tmo_b;

        if constexpr (!WithGNSS) {
            sigma_b_ =
                computeNoGnssAttitudeSigmaBody(f_b, fhat_n_, aux, R_nb_before);
        }

        Vec3 bias_dot =
            projectBallDerivative(
                gyro_bias_b_,
                -kI_ * sigma_b_,
                cfg_.gyro_bias_limit_rad_s
            );

        if constexpr (Mag == NloMagType::None) {
            bias_dot.z() = R(0);
            gyro_bias_b_.z() = R(0);
        }

        gyro_bias_b_ += dt * bias_dot;
        gyro_bias_b_ = clampNorm(gyro_bias_b_, cfg_.gyro_bias_limit_rad_s);

        if constexpr (Mag == NloMagType::None) {
            gyro_bias_b_.z() = R(0);
        }

        const Vec3 omega_corr_b = gyro_b_rad_s - gyro_bias_b_ + sigma_b_;
        integrateQuaternionRight(dt, omega_corr_b);

        const Mat3 R_nb = q_nb_.toRotationMatrix();
        fhat_n_ = R_nb * f_b + xi_n_;

        const R raw_p0z_err = R(0) - p0z_hat_;
        const R p0z_err = highpassP0zInnovation(dt, raw_p0z_err);

        Vec2 pxy_err = Vec2::Zero();
        const bool have_gnss = hasUsableGNSS(aux);
        if (have_gnss) {
            pxy_err = aux.gnss.position_n_m - p_n_.template head<2>();
        }

        const R th = cfg_.theta;
        const R th2 = th * th;
        const R th3 = th2 * th;
        const R th4 = th2 * th2;

        const R s1 = vartheta_ * th;
        const R s2 = vartheta_ * th2;
        const R s3 = vartheta_ * th3;
        const R s4 = vartheta_ * th4;

        const Mat2 Kpp  = cfg_.K_pp_scalar  * Mat2::Identity();
        const Mat2 Kvp  = cfg_.K_vp_scalar  * Mat2::Identity();
        const Mat2 Kxip = cfg_.K_xip_scalar * Mat2::Identity();

        const R p0z_dot =
            p_n_.z() + s1 * cfg_.K_p0z_p0z * p0z_err;

        Vec3 p_dot = v_n_;
        if (have_gnss) {
            p_dot.template head<2>() += s2 * (Kpp * pxy_err);
        }
        p_dot.z() += s2 * cfg_.K_pz_p0z * p0z_err;

        Vec3 gravity_n;
        gravity_n << R(0), R(0), cfg_.gravity_mps2;

        Vec3 v_dot = fhat_n_ + gravity_n;
        if (have_gnss) {
            v_dot.template head<2>() += s3 * (Kvp * pxy_err);
        }
        v_dot.z() += s3 * cfg_.K_vz_p0z * p0z_err;

        /*
          Use sigma_tmo_b here, not sigma_b_.

          sigma_b_ may be a small no-GNSS vertical-reference attitude leak.
          Feeding that into xi_dot destabilizes the VVR/heave loop.
        */
        Vec3 xi_dot = -R_nb * (skew(sigma_tmo_b) * f_b);

        if (have_gnss) {
            xi_dot.template head<2>() += s4 * (Kxip * pxy_err);
        }
        xi_dot.z() += s4 * cfg_.K_xiz_p0z * p0z_err;

        if constexpr (!WithGNSS) {
            p_dot.x() = R(0);
            p_dot.y() = R(0);

            v_dot.x() = R(0);
            v_dot.y() = R(0);

            xi_dot.x() = R(0);
            xi_dot.y() = R(0);

            p_n_.x() = R(0);
            p_n_.y() = R(0);
            v_n_.x() = R(0);
            v_n_.y() = R(0);
            xi_n_.x() = R(0);
            xi_n_.y() = R(0);
        }

        p0z_hat_ += dt * p0z_dot;
        p_n_     += dt * p_dot;
        v_n_     += dt * v_dot;
        xi_n_    += dt * xi_dot;

        lockUnobservableStates_();

        fhat_n_ = q_nb_.toRotationMatrix() * f_b + xi_n_;
    }

    const Config& config() const { return cfg_; }
    Config& config() { return cfg_; }

    const Quat& quaternionBodyToNED() const { return q_nb_; }
    Mat3 rotationBodyToNED() const { return q_nb_.toRotationMatrix(); }

    const Vec3& positionNED() const { return p_n_; }
    const Vec3& velocityNED() const { return v_n_; }
    const Vec3& specificForceNED() const { return fhat_n_; }
    const Vec3& xiNED() const { return xi_n_; }

    const Vec3& gyroBiasBody() const { return gyro_bias_b_; }
    const Vec3& sigmaBody() const { return sigma_b_; }

    R integratedVerticalPositionState() const { return p0z_hat_; }

    R gainK1() const { return k1_; }
    R gainK2() const { return k2_; }
    R gainKI() const { return kI_; }
    R gainVartheta() const { return vartheta_; }
    R timeSeconds() const { return t_s_; }

    Vec3 eulerRad() const {
        const Mat3 Rnb = q_nb_.toRotationMatrix();

        Vec3 e;
        e.x() = std::atan2(Rnb(2, 1), Rnb(2, 2));
        e.y() = -std::asin(clamp(Rnb(2, 0), R(-1), R(1)));
        e.z() = std::atan2(Rnb(1, 0), Rnb(0, 0));
        return e;
    }

private:
    Config cfg_;

    R t_s_ = R(0);

    Quat q_nb_ = Quat::Identity();

    Vec3 gyro_bias_b_ = Vec3::Zero();
    Vec3 accel_bias_b_ = Vec3::Zero();

    Vec3 p_n_ = Vec3::Zero();
    Vec3 v_n_ = Vec3::Zero();
    Vec3 xi_n_ = Vec3::Zero();
    Vec3 fhat_n_ = Vec3::Zero();

    R p0z_hat_ = R(0);
    R p0z_hp_state_ = R(0);

    R k1_ = R(0.55);
    R k2_ = R(1.0);
    R kI_ = R(0.01);

    R gnss_rms_lpf_ = R(1);
    R vartheta2_ = R(1);
    R vartheta_ = R(1);

    Vec3 sigma_b_ = Vec3::Zero();

    void lockUnobservableStates_() {
        if constexpr (!WithGNSS) {
            p_n_.x() = R(0);
            p_n_.y() = R(0);
            v_n_.x() = R(0);
            v_n_.y() = R(0);
            xi_n_.x() = R(0);
            xi_n_.y() = R(0);
        }

        if constexpr (Mag == NloMagType::None) {
            gyro_bias_b_.z() = R(0);
        }
    }

    static bool isFinite(R x) {
        using std::isfinite;
        return isfinite(x);
    }

    static R clamp(R x, R lo, R hi) {
        return (x < lo) ? lo : ((x > hi) ? hi : x);
    }

    static Vec3 nedNorth() {
        Vec3 v;
        v << R(1), R(0), R(0);
        return v;
    }

    static Vec3 nedDown() {
        Vec3 v;
        v << R(0), R(0), R(1);
        return v;
    }

    Vec3 magneticReferenceNED() const {
        Vec3 v;
        v << std::cos(cfg_.magnetic_declination_rad),
             std::sin(cfg_.magnetic_declination_rad),
             R(0);
        return v;
    }

    static Mat3 skew(const Vec3& x) {
        Mat3 S;
        S << R(0),   -x.z(),  x.y(),
             x.z(),  R(0),   -x.x(),
            -x.y(),  x.x(),  R(0);
        return S;
    }

    static Vec3 clampNorm(const Vec3& v, R max_norm) {
        const R n = v.norm();
        if (!(max_norm > R(0)) || n <= max_norm || n <= R(1e-12)) {
            return v;
        }
        return v * (max_norm / n);
    }

    static Vec3 normalizeSafe(const Vec3& v, const Vec3& fallback) {
        const R n = v.norm();
        if (n > R(1e-9) && isFinite(n)) {
            return v / n;
        }
        return fallback;
    }

    static Vec3 projectBallDerivative(const Vec3& x, const Vec3& u, R radius) {
        if (!(radius > R(0))) {
            return Vec3::Zero();
        }

        const R n = x.norm();
        if (n < radius || n <= R(1e-12)) {
            return u;
        }

        const Vec3 nhat = x / n;
        const R outward = nhat.dot(u);
        if (outward <= R(0)) {
            return u;
        }

        return u - outward * nhat;
    }

    bool setAttitudeFromDownAndHorizontalReference(const Vec3& down_b_raw,
                                                   const Vec3& h_b_raw,
                                                   const Vec3& h_n_raw) {
        Vec3 down_b = normalizeSafe(down_b_raw, Vec3(R(0), R(0), R(1)));

        Vec3 h_b = h_b_raw - down_b * down_b.dot(h_b_raw);
        if (h_b.norm() < R(1e-6)) {
            return false;
        }
        h_b.normalize();

        const Vec3 down_n = nedDown();

        Vec3 h_n = h_n_raw - down_n * down_n.dot(h_n_raw);
        if (h_n.norm() < R(1e-6)) {
            return false;
        }
        h_n.normalize();

        Vec3 east_b = down_b.cross(h_b);
        if (east_b.norm() < R(1e-6)) {
            return false;
        }
        east_b.normalize();

        Vec3 third_b = h_b.cross(east_b);
        third_b.normalize();

        Vec3 east_n = down_n.cross(h_n);
        if (east_n.norm() < R(1e-6)) {
            return false;
        }
        east_n.normalize();

        Vec3 third_n = h_n.cross(east_n);
        third_n.normalize();

        Mat3 T_b;
        T_b.col(0) = h_b;
        T_b.col(1) = east_b;
        T_b.col(2) = third_b;

        Mat3 T_n;
        T_n.col(0) = h_n;
        T_n.col(1) = east_n;
        T_n.col(2) = third_n;

        const Mat3 R_bn = T_b * T_n.transpose();
        const Mat3 R_nb = R_bn.transpose();

        q_nb_ = Quat(R_nb);
        q_nb_.normalize();

        return true;
    }

    void integrateQuaternionRight(R dt, const Vec3& omega_b) {
        const R omega_norm = omega_b.norm();
        const R angle = omega_norm * dt;

        Quat dq;
        if (angle < R(1e-8) || omega_norm < R(1e-12)) {
            const Vec3 half = R(0.5) * dt * omega_b;
            dq = Quat(R(1), half.x(), half.y(), half.z());
        } else {
            const Vec3 axis = omega_b / omega_norm;
            const R half_angle = R(0.5) * angle;
            const R s = std::sin(half_angle);
            dq = Quat(std::cos(half_angle),
                      s * axis.x(),
                      s * axis.y(),
                      s * axis.z());
        }

        q_nb_ = q_nb_ * dq;
        q_nb_.normalize();
    }

    bool hasUsableGNSS(const Aux& aux) const {
        if constexpr (!WithGNSS) {
            (void)aux;
            return false;
        } else {
            return aux.gnss.valid;
        }
    }

    void updateGainSchedules(R dt, const Aux& aux) {
        if (cfg_.use_time_varying_attitude_gains) {
            Vec3 target;
            if (t_s_ <= cfg_.attitude_gain_switch_s) {
                target << cfg_.k1_initial, cfg_.k2_initial, cfg_.kI_initial;
            } else {
                target << cfg_.k1_nominal, cfg_.k2_nominal, cfg_.kI_nominal;
            }

            const R a = clamp(dt / cfg_.attitude_gain_tau_s, R(0), R(1));
            k1_ += a * (target.x() - k1_);
            k2_ += a * (target.y() - k2_);
            kI_ += a * (target.z() - kI_);
        } else {
            k1_ = cfg_.k1_nominal;
            k2_ = cfg_.k2_nominal;
            kI_ = cfg_.kI_nominal;
        }

        if (!cfg_.use_time_varying_tmo_gain) {
            vartheta_ = R(1);
            return;
        }

        R vartheta1 = cfg_.vartheta1_without_gnss;

        if constexpr (WithGNSS) {
            if (aux.gnss.valid && aux.gnss.rms_xy_m > R(0)) {
                const R a = clamp(dt / cfg_.gnss_rms_lpf_tau_s, R(0), R(1));
                gnss_rms_lpf_ += a * (aux.gnss.rms_xy_m - gnss_rms_lpf_);
                vartheta1 =
                    cfg_.vartheta1_b * std::exp(-cfg_.vartheta1_a * gnss_rms_lpf_);
            }
        } else {
            (void)aux;
        }

        const R target2 = (t_s_ <= cfg_.vartheta2_switch_s) ? R(1) : R(0);
        const R a2 = clamp(dt / cfg_.vartheta2_tau_s, R(0), R(1));
        vartheta2_ += a2 * (target2 - vartheta2_);

        vartheta_ = cfg_.vartheta0 + vartheta1 + vartheta2_;

        if (vartheta_ < R(1e-3)) {
            vartheta_ = R(1e-3);
        }
    }

    R highpassP0zInnovation(R dt, R raw) {
        const R T = cfg_.p0z_highpass_tau_s;
        if (!(T > R(0))) {
            return raw;
        }

        p0z_hp_state_ += dt * (raw - p0z_hp_state_ / T);
        return raw - p0z_hp_state_ / T;
    }

    bool makeCompassReference(const Aux& aux,
                              Vec3& c_b,
                              Vec3& c_n) const {
        if constexpr (Mag != NloMagType::Compass) {
            (void)aux;
            (void)c_b;
            (void)c_n;
            return false;
        } else {
            if (!aux.compass_valid) {
                return false;
            }

            c_b << std::cos(aux.compass_heading_rad),
                  -std::sin(aux.compass_heading_rad),
                   R(0);

            c_b = normalizeSafe(c_b, Vec3(R(1), R(0), R(0)));
            c_n = nedNorth();
            return true;
        }
    }

    bool makeMagReference(const Aux& aux,
                          const Mat3& R_nb,
                          Vec3& c_b,
                          Vec3& c_n) const {
        if constexpr (Mag != NloMagType::Magnetometer) {
            (void)aux;
            (void)R_nb;
            (void)c_b;
            (void)c_n;
            return false;
        } else {
            if (!aux.mag_valid) {
                return false;
            }

            const R mn = aux.mag_b.norm();
            if (!isFinite(mn) || mn < R(1e-9)) {
                return false;
            }

            if (cfg_.expected_mag_norm > R(0)) {
                const R rel =
                    std::abs(mn - cfg_.expected_mag_norm) / cfg_.expected_mag_norm;
                if (rel > cfg_.max_mag_relative_norm_error) {
                    return false;
                }
            }

            Vec3 m_b = aux.mag_b / mn;

            const Vec3 down_b = R_nb.transpose() * nedDown();

            Vec3 mh_b = m_b - down_b * down_b.dot(m_b);
            if (mh_b.norm() < R(1e-6)) {
                return false;
            }

            c_b = mh_b.normalized();
            c_n = magneticReferenceNED();
            return true;
        }
    }

    bool makeYawReference(const Aux& aux,
                          const Mat3& R_nb,
                          Vec3& c_b,
                          Vec3& c_n) const {
        if constexpr (Mag == NloMagType::None) {
            (void)aux;
            (void)R_nb;
            (void)c_b;
            (void)c_n;
            return false;
        } else if constexpr (Mag == NloMagType::Compass) {
            return makeCompassReference(aux, c_b, c_n);
        } else if constexpr (Mag == NloMagType::Magnetometer) {
            return makeMagReference(aux, R_nb, c_b, c_n);
        } else {
            return false;
        }
    }

    Vec3 computeNoGnssAttitudeSigmaBody(const Vec3& f_b,
                                        const Vec3& fhat_n,
                                        const Aux& aux,
                                        const Mat3& R_nb) const {
        Vec3 fhat_att_n = fhat_n;
        fhat_att_n.x() = R(0);
        fhat_att_n.y() = R(0);

        const Vec3 fhat_sat_n =
            clampNorm(fhat_att_n, cfg_.max_specific_force_mps2);

        Vec3 sigma_force = Vec3::Zero();

        if (!cfg_.use_triad_style_force_injection) {
            sigma_force = k1_ * f_b.cross(R_nb.transpose() * fhat_sat_n);
        } else {
            const Vec3 f_b_u =
                normalizeSafe(f_b, Vec3(R(0), R(0), R(-1)));

            const Vec3 f_n_u =
                normalizeSafe(fhat_sat_n, Vec3(R(0), R(0), R(-1)));

            sigma_force =
                k1_ * f_b_u.cross(R_nb.transpose() * f_n_u);
        }

        /*
          Deliberately weak leak path:
            - improves no-GNSS roll/pitch self-cancellation
            - does not drive xi_dot / vertical TMO/VVR
            - must stay small enough not to corrupt heave in steep waves
        */
        sigma_force *= cfg_.no_gnss_attitude_sigma_scale;
        sigma_force = clampNorm(
            sigma_force,
            cfg_.no_gnss_attitude_sigma_limit_rad_s
        );

        if constexpr (Mag == NloMagType::None) {
            sigma_force.z() = R(0);
        }

        Vec3 sigma_yaw = Vec3::Zero();

        Vec3 c_b;
        Vec3 c_n;
        if (makeYawReference(aux, R_nb, c_b, c_n)) {
            if (!cfg_.use_triad_style_force_injection) {
                sigma_yaw = k2_ * c_b.cross(R_nb.transpose() * c_n);
            } else {
                const Vec3 f_b_u =
                    normalizeSafe(f_b, Vec3(R(0), R(0), R(-1)));

                const Vec3 f_n_u =
                    normalizeSafe(fhat_sat_n, Vec3(R(0), R(0), R(-1)));

                Vec3 b2 = f_b_u.cross(c_b);
                Vec3 n2 = f_n_u.cross(c_n);

                b2 = normalizeSafe(b2, Vec3(R(0), R(1), R(0)));
                n2 = normalizeSafe(n2, Vec3(R(0), R(1), R(0)));

                sigma_yaw = k2_ * b2.cross(R_nb.transpose() * n2);
            }
        }

        return sigma_force + sigma_yaw;
    }

    Vec3 computeSigmaBody(const Vec3& f_b,
                          const Vec3& fhat_n,
                          const Aux& aux,
                          const Mat3& R_nb) const {
        const Vec3 fhat_sat_n =
            clampNorm(fhat_n, cfg_.max_specific_force_mps2);

        if (!cfg_.use_triad_style_force_injection) {
            Vec3 sigma = k1_ * f_b.cross(R_nb.transpose() * fhat_sat_n);

            Vec3 c_b;
            Vec3 c_n;
            if (makeYawReference(aux, R_nb, c_b, c_n)) {
                sigma += k2_ * c_b.cross(R_nb.transpose() * c_n);
            }

            return sigma;
        }

        const Vec3 f_b_u =
            normalizeSafe(f_b, Vec3(R(0), R(0), R(-1)));

        const Vec3 f_n_u =
            normalizeSafe(fhat_sat_n, Vec3(R(0), R(0), R(-1)));

        Vec3 sigma =
            k1_ * f_b_u.cross(R_nb.transpose() * f_n_u);

        Vec3 c_b;
        Vec3 c_n;
        if (!makeYawReference(aux, R_nb, c_b, c_n)) {
            return sigma;
        }

        Vec3 b2 = f_b_u.cross(c_b);
        Vec3 n2 = f_n_u.cross(c_n);

        b2 = normalizeSafe(b2, Vec3(R(0), R(1), R(0)));
        n2 = normalizeSafe(n2, Vec3(R(0), R(1), R(0)));

        sigma += k2_ * b2.cross(R_nb.transpose() * n2);

        return sigma;
    }
};
