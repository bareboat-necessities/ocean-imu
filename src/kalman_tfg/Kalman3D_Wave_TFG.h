#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    Right-invariant two-frame Lie-group error-state EKF for marine wave state.

    R = R_bw (body -> world), eta = Xhat o X^-1, and corrections are applied
    on the left. The world state is [v p S a_w]; gyro bias is a body-frame
    random walk and the residual accelerometer bias is the same slow OU process
    used by current OU-III.

    This is intentionally described as a right-invariant Lie-group error-state
    EKF with invariant-style measurements. The complete marine process is not
    claimed to be a group-affine/log-linear InEKF: the world/bias couplings
    below still depend on the estimated world state.

    Tangent layout matches OU-III:

        xi = [ phi | beta_g | rho_v rho_p rho_S rho_a | beta_a ]
               0      3        6     9     12    15      18

    Level 2 (the default) uses beta_b = R (bhat-b). Level 1 keeps additive
    body-frame bias errors as an ablation. The two are related during
    propagation by the time-varying coordinate map beta = R delta_b.

    The world OU chain is discretized with IntegratedOUChain. The deterministic
    gyro-bias/world coupling and the gyro/gyro-bias process covariance use a
    full structured Gauss--Legendre discretization of the invariant error
    dynamics. This transports each noise source through attitude and every
    world column instead of applying zero-rate or first-order leakage formulas.

    After every finite measurement correction a=K r, covariance is transported
    into the post-injection tangent coordinates by the full retraction Jacobian

        J_reset = d/d(delta) Log_G( Exp_G(-a) Exp_G(a+delta) ) | delta=0,

    rather than by OU-III's attitude-only MEKF reset.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>

#include "kalman_ou_common/KalmanOUCoreMath.h"
#include "lie/SO3Jacobians.h"
#include "lie/TwoFrameGroup.h"

namespace ocean_imu::kalman {

namespace tfg_detail {

using ocean_imu::kalman::ou_detail::safe_inv_tau;

template<typename T>
inline void integral_transition_axis(T tau, T h, Eigen::Matrix<T,4,4>& Psi) {
    const T inv_tau = safe_inv_tau(tau);
    const T x = h * inv_tau;
    const T tau2 = tau * tau;
    const T tau3 = tau2 * tau;
    const T tau4 = tau3 * tau;

    Psi.setZero();
    Psi(0,0) = h;
    Psi(1,0) = T(0.5) * h * h;
    Psi(2,0) = h * h * h / T(6);
    Psi(1,1) = h;
    Psi(2,1) = T(0.5) * h * h;
    Psi(2,2) = h;

    if (std::abs(x) < T(1e-2)) {
        const T x2 = x*x, x3 = x2*x, x4 = x3*x, x5 = x4*x;
        const T x6 = x5*x, x7 = x6*x, x8 = x7*x;
        Psi(3,3) = tau  * (x - x2/T(2) + x3/T(6) - x4/T(24) + x5/T(120) - x6/T(720));
        Psi(0,3) = tau2 * (x2/T(2) - x3/T(6) + x4/T(24) - x5/T(120) + x6/T(720)
                           - x7/T(5040));
        Psi(1,3) = tau3 * (x3/T(6) - x4/T(24) + x5/T(120) - x6/T(720)
                           + x7/T(5040) - x8/T(40320));
        Psi(2,3) = tau4 * (x4/T(24) - x5/T(120) + x6/T(720) - x7/T(5040)
                           + x8/T(40320));
    } else {
        const T em1 = std::expm1(-x);
        const T x2 = x*x, x3 = x2*x;
        Psi(3,3) = -tau  * em1;
        Psi(0,3) =  tau2 * (em1 + x);
        Psi(1,3) =  tau3 * (x2/T(2) - em1 - x);
        Psi(2,3) =  tau4 * (x3/T(6) - x2/T(2) + em1 + x);
    }
}

} // namespace tfg_detail

template <typename T = float,
          bool with_gyro_bias = true,
          bool with_accel_bias = true,
          bool two_frame_bias = true>
class Kalman3D_Wave_TFG {
    static_assert(with_gyro_bias,
                  "the two-frame tangent layout pins the gyro bias at offset 3");

  public:
    using Group   = ocean_imu::lie::TwoFrameGroup<T, 4, with_accel_bias ? 2 : 1>;
    using Scalar  = T;
    using Vector3 = Eigen::Matrix<T,3,1>;
    using Matrix3 = Eigen::Matrix<T,3,3>;

    static constexpr int NX = Group::NX;
    static constexpr int OFF_PHI = Group::OFF_PHI;
    static constexpr int OFF_BG  = Group::OFF_BG;
    static constexpr int OFF_V   = Group::world_offset(0);
    static constexpr int OFF_P   = Group::world_offset(1);
    static constexpr int OFF_S   = Group::world_offset(2);
    static constexpr int OFF_AW  = Group::world_offset(3);
    static constexpr int OFF_BA  = with_accel_bias ? Group::OFF_BA : -1;

    using MatrixNX = Eigen::Matrix<T,NX,NX>;
    using Tangent  = Eigen::Matrix<T,NX,1>;

    static constexpr T STD_GRAVITY = T(9.80665);

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    explicit Kalman3D_Wave_TFG(T gyro_noise_density_rad_sqrt_s = T(0.005),
                               T gravity_magnitude = STD_GRAVITY)
        : gravity_magnitude_(gravity_magnitude) {
        set_gyro_noise_density_rad_sqrt_s(gyro_noise_density_rad_sqrt_s);
        P_.setIdentity();
        P_ *= T(1e-4);
    }

    // The seed is an argument so it can be swept; the default is the value
    // this filter has always used.  See docs/tfg-qmekf-variances.md.
    void initialize_identity(T initial_covariance = T(1e-4)) {
        X_ = Group::Identity();
        P_.setIdentity();
        P_ *= std::max(T(0), initial_covariance);
    }

    void initialize_from_truth(const Eigen::Quaternion<T>& q_bw,
                               const Vector3& v, const Vector3& p,
                               const Vector3& S, const Vector3& a_w,
                               const Vector3& b_g = Vector3::Zero(),
                               const Vector3& b_a = Vector3::Zero()) {
        Eigen::Quaternion<T> q = q_bw;
        q.normalize();
        X_.R = q.toRotationMatrix();
        X_.X.col(0) = v;
        X_.X.col(1) = p;
        X_.X.col(2) = S;
        X_.X.col(3) = a_w;
        X_.B.col(0) = b_g;
        if constexpr (with_accel_bias) X_.B.col(1) = b_a;
    }

    [[nodiscard]] Eigen::Quaternion<T> quaternion() const {
        Eigen::Quaternion<T> q(X_.R);
        q.normalize();
        return q;
    }
    [[nodiscard]] Matrix3 R_bw() const { return X_.R; }
    [[nodiscard]] Matrix3 R_wb() const { return X_.R.transpose(); }
    [[nodiscard]] Vector3 get_velocity() const              { return X_.X.col(0); }
    [[nodiscard]] Vector3 get_position() const              { return X_.X.col(1); }
    [[nodiscard]] Vector3 get_integral_displacement() const { return X_.X.col(2); }
    [[nodiscard]] Vector3 get_world_accel() const           { return X_.X.col(3); }
    [[nodiscard]] Vector3 gyroscope_bias() const            { return X_.B.col(0); }
    [[nodiscard]] Vector3 get_acc_bias() const {
        if constexpr (with_accel_bias) return X_.B.col(1);
        return Vector3::Zero();
    }
    [[nodiscard]] const MatrixNX& covariance_full() const { return P_; }
    [[nodiscard]] MatrixNX& covariance_full() { return P_; }
    [[nodiscard]] const Group& state() const { return X_; }
    [[nodiscard]] Group& state() { return X_; }

    bool initialize_from_acc(const Vector3& acc_body) {
        if (!acc_body.allFinite()) return false;
        const T n = acc_body.norm();
        if (!(n > T(1e-6))) return false;

        const Vector3 up_body = acc_body / n;
        const Vector3 up_world(T(0), T(0), T(-1));
        const Vector3 v = up_body.cross(up_world);
        const T c = up_body.dot(up_world);
        const T s2 = v.squaredNorm();
        Matrix3 R;
        if (s2 < T(1e-18)) {
            if (c > T(0)) {
                R = Matrix3::Identity();
            } else {
                Vector3 axis = up_body.cross(Vector3::UnitX());
                if (axis.norm() < T(1e-6)) axis = up_body.cross(Vector3::UnitY());
                axis.normalize();
                R = ocean_imu::lie::Exp<T>(Vector3(axis * T(M_PI)));
            }
        } else {
            const Matrix3 V = ocean_imu::lie::skew<T>(v);
            R = Matrix3::Identity() + V + V*V*((T(1)-c)/s2);
        }
        X_.R = R;
        reorthonormalize_();
        return true;
    }

    void set_linear_block_enabled(bool on) { linear_block_enabled_ = on; }
    [[nodiscard]] bool linear_block_enabled() const { return linear_block_enabled_; }

    void reset_aw_covariance_to_stationary() {
        P_.template block<3,3>(OFF_AW, OFF_AW) = Sigma_aw_;
        for (int i = 0; i < NX; ++i) {
            if (i >= OFF_AW && i < OFF_AW + 3) continue;
            for (int k = 0; k < 3; ++k) {
                P_(i, OFF_AW+k) = T(0);
                P_(OFF_AW+k, i) = T(0);
            }
        }
    }

    // Re-align the posterior a_w marginal with the stationary prior without
    // throwing away what the filter has learned. Unlike
    // reset_aw_covariance_to_stationary() above, which replaces the block and
    // decorrelates it, this only ever adds the PSD part of the shortfall, and
    // it does so at the next time update rather than in the middle of a
    // measurement sequence. OU-III runs exactly this at the adaptation cadence
    // so that a stationary scale which has moved since the last discrete
    // reconfiguration still reaches the marginal.
    void synchronize_aw_covariance_to_stationary() {
        aw_cov_floor_target_ = T(0.5) * (Sigma_aw_ + Sigma_aw_.transpose());
        aw_cov_floor_pending_ = true;
    }

    // Rewrite the heading of the attitude, keeping the estimated tilt, the
    // world columns and the biases untouched.
    //
    // This is a state correction and not a gauge change: the world frame is
    // already the magnetic one fixed at first lock, and what the second-stage
    // magnetic acquisition finds is that the *heading in it* was wrong. Turning
    // this into apply_world_yaw_gauge() would rotate v, p and S along with R
    // and corrupt the displacement history the correction has nothing to say
    // about. Covariance is left alone, matching OU-III's set_quaternion_boat().
    void set_attitude_yaw_absolute(T yaw_abs_rad) {
        if (!std::isfinite(yaw_abs_rad)) return;
        const Matrix3 R = X_.R;
        const T yaw = std::atan2(R(1,0), R(0,0));
        if (!std::isfinite(yaw)) return;
        X_.R = (ocean_imu::lie::Exp<T>(Vector3(T(0), T(0), yaw_abs_rad - yaw)) * R).eval();
        reorthonormalize_();
    }

    void set_aw_time_constant(T tau) { tau_aw_ = std::max(T(1e-3), tau); }
    [[nodiscard]] T aw_time_constant() const { return tau_aw_; }

    void set_aw_stationary_std(const Vector3& sigma) {
        Sigma_aw_ = sigma.cwiseAbs2().asDiagonal();
    }
    void set_aw_stationary_cov_full(const Matrix3& Sigma) {
        Sigma_aw_ = T(0.5) * (Sigma + Sigma.transpose());
    }
    [[nodiscard]] const Matrix3& aw_stationary_cov() const { return Sigma_aw_; }

    void set_gyro_noise_density_rad_sqrt_s(T density) {
        const T d = std::max(T(0), density);
        Q_gyro_ = Matrix3::Identity() * (d*d);
    }
    void set_Q_bgyro_rw(const Vector3& var_per_s) {
        Q_bg_ = var_per_s.cwiseAbs().asDiagonal();
    }

    // Same API and semantics as current OU-III: the argument is a continuous
    // driving-noise standard deviation, despite the historical method name.
    void set_Q_bacc_rw(const Vector3& driving_std_per_sqrt_s) {
        if constexpr (with_accel_bias)
            Q_ba_ = driving_std_per_sqrt_s.array().square().matrix().asDiagonal();
    }
    void set_acc_bias_time_constant(T tau_seconds) {
        if constexpr (with_accel_bias) tau_ba_ = std::max(T(1e-3), tau_seconds);
    }
    [[nodiscard]] T get_acc_bias_time_constant() const noexcept { return tau_ba_; }
    void set_acc_bias_ou_stationary_std(const Vector3& stationary_std, T tau_seconds) {
        if constexpr (with_accel_bias) {
            set_acc_bias_time_constant(tau_seconds);
            const Vector3 s = stationary_std.cwiseAbs();
            Q_ba_ = (T(2)/tau_ba_) * s.array().square().matrix().asDiagonal();
        }
    }
    void set_initial_acc_bias_std(T s) {
        if constexpr (with_accel_bias) {
            sigma_ba0_ = std::max(T(0), s);
            P_.template block<3,3>(OFF_BA, OFF_BA) =
                Matrix3::Identity() * sigma_ba0_ * sigma_ba0_;
        }
    }
    void set_initial_acc_bias(const Vector3& b0) {
        if constexpr (with_accel_bias) X_.B.col(1) = b0;
    }
    [[nodiscard]] Vector3 get_acc_bias_at_temperature(T tempC) const {
        return accel_bias_at_(tempC);
    }

    void set_initial_linear_uncertainty(T std_v, T std_p, T std_S, T std_aw) {
        P_.template block<12,12>(OFF_V, OFF_V).setZero();
        const T d[4] = {std_v*std_v, std_p*std_p, std_S*std_S, std_aw*std_aw};
        for (int c = 0; c < 4; ++c)
            for (int k = 0; k < 3; ++k) P_(OFF_V + 3*c + k, OFF_V + 3*c + k) = d[c];
    }

    void propagate_mean(const Vector3& gyro_meas, T Ts) {
        if (!(Ts > T(0)) || !std::isfinite(Ts)) return;
        const Vector3 omega = gyro_meas - X_.B.col(0);
        last_omega_ = omega;

        X_.R = X_.R * ocean_imu::lie::Exp<T>(Vector3(omega * Ts));
        reorthonormalize_();

        if (linear_block_enabled_) {
            Eigen::Matrix<T,4,4> Phi_axis;
            ou_detail::IntegratedOUChain<T,3>::transition(tau_aw_, Ts, Phi_axis);
            X_.X = (X_.X * Phi_axis.transpose()).eval();
        }

        if constexpr (with_accel_bias) {
            if (acc_bias_updates_enabled_) {
                X_.B.col(1) *= std::exp(-Ts / tau_ba_);
            }
        }
    }

    void build_transition(T Ts, const Vector3& gyro_meas, MatrixNX& Phi) const {
        Phi.setIdentity();
        if (!(Ts > T(0)) || !std::isfinite(Ts)) return;

        const Matrix3 R0 = X_.R;
        const Vector3 omega = gyro_meas - X_.B.col(0);
        const Matrix3 R_end = R0 * ocean_imu::lie::Exp<T>(Vector3(omega * Ts));

        // beta_g is world-referred at Level 2 and additive body-referred at
        // Level 1.  Fh = int_0^h R(t) dt maps a physical body-bias error into
        // the right-invariant attitude error exactly for constant measured
        // angular rate over the sample.
        const Matrix3 Fh =
            R0 * (Ts * ocean_imu::lie::left_jacobian<T>(Vector3(omega * Ts)));
        const Matrix3 bg0_to_body =
            two_frame_bias ? Matrix3(R0.transpose()) : Matrix3::Identity();
        Phi.template block<3,3>(OFF_PHI, OFF_BG) = -Fh * bg0_to_body;
        if constexpr (two_frame_bias)
            Phi.template block<3,3>(OFF_BG, OFF_BG) = R_end * R0.transpose();

        if constexpr (with_accel_bias) {
            const T alpha_ba = acc_bias_updates_enabled_ ? std::exp(-Ts/tau_ba_) : T(1);
            if constexpr (two_frame_bias)
                Phi.template block<3,3>(OFF_BA, OFF_BA) =
                    alpha_ba * R_end * R0.transpose();
            else
                Phi.template block<3,3>(OFF_BA, OFF_BA) =
                    alpha_ba * Matrix3::Identity();
        }

        // Exact OU-chain state transition for the active linear block; when
        // it is gated off the physical world columns are frozen and the world
        // transition is identity.
        Eigen::Matrix<T,4,4> Phi_axis;
        world_transition_(Ts, Phi_axis);
        Phi.template block<12,12>(OFF_V, OFF_V).setZero();
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (Phi_axis(i,j) == T(0)) continue;
                for (int k = 0; k < 3; ++k)
                    Phi(OFF_V + 3*i + k, OFF_V + 3*j + k) = Phi_axis(i,j);
            }
        }

        // rho_dot = A rho - [x(t)]x R(t) delta_b_g.  Integrate the actual
        // time-varying nominal x(t) and R(t), rather than replacing them by
        // their start-of-step values.  D maps a physical body bias to the
        // final world-error columns; Level 2 converts beta_g(0) back to that
        // physical body-bias coordinate with R0^T.
        Eigen::Matrix<T,12,3> D;
        integrate_world_bias_impulse_(T(0), Ts, omega, D);
        Phi.template block<12,3>(OFF_V, OFF_BG) = -D * bg0_to_body;
    }

    void build_process_noise(T Ts, const Vector3& gyro_meas, MatrixNX& Qd) const {
        Qd.setZero();
        if (!(Ts > T(0)) || !std::isfinite(Ts)) return;

        const Matrix3 R0 = X_.R;
        const Vector3 omega = gyro_meas - X_.B.col(0);
        const Matrix3 R_end = R0 * ocean_imu::lie::Exp<T>(Vector3(omega * Ts));
        const Matrix3 Fh =
            R0 * (Ts * ocean_imu::lie::left_jacobian<T>(Vector3(omega * Ts)));

        // Five-point Gauss--Legendre quadrature of the *complete structured
        // impulse response* of the stated continuous invariant error model.
        // For gyro white noise an impulse at s enters phi through R(s) and all
        // world columns through [x(s)]x R(s), followed by the remaining OU
        // chain.  A gyro-bias RW impulse persists from s to h, therefore also
        // drives phi and every world column for the rest of the sample.
        static constexpr double nodes[5] = {
            -0.9061798459386640, -0.5384693101056831, 0.0,
             0.5384693101056831,  0.9061798459386640
        };
        static constexpr double weights[5] = {
            0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
            0.4786286704993665, 0.2369268850561891
        };

        for (int q = 0; q < 5; ++q) {
            const T s = T(0.5) * Ts * (T(1) + T(nodes[q]));
            const T wq = T(0.5) * Ts * T(weights[q]);
            const Matrix3 R_s = R0 * ocean_imu::lie::Exp<T>(Vector3(omega * s));

            Eigen::Matrix<T,4,4> Phi_s, Phi_rem;
            world_transition_(s, Phi_s);
            world_transition_(Ts - s, Phi_rem);
            const Eigen::Matrix<T,3,4> X_s = X_.X * Phi_s.transpose();

            Eigen::Matrix<T,NX,3> Lg = Eigen::Matrix<T,NX,3>::Zero();
            Lg.template block<3,3>(OFF_PHI, 0) = R_s;
            for (int i = 0; i < 4; ++i) {
                Matrix3 M = Matrix3::Zero();
                for (int j = 0; j < 4; ++j) {
                    if (Phi_rem(i,j) == T(0)) continue;
                    M.noalias() += Phi_rem(i,j) *
                        (ocean_imu::lie::skew<T>(Vector3(X_s.col(j))) * R_s);
                }
                Lg.template block<3,3>(OFF_V + 3*i, 0) = M;
            }
            Qd.noalias() += wq * (Lg * Q_gyro_ * Lg.transpose());

            Eigen::Matrix<T,NX,3> Lb = Eigen::Matrix<T,NX,3>::Zero();
            const Matrix3 Fs =
                R0 * (s * ocean_imu::lie::left_jacobian<T>(Vector3(omega * s)));
            const Matrix3 C = Fh - Fs;  // int_s^h R(t) dt
            Lb.template block<3,3>(OFF_PHI, 0) = -C;
            Lb.template block<3,3>(OFF_BG, 0) =
                two_frame_bias ? R_end : Matrix3::Identity();

            Eigen::Matrix<T,12,3> D;
            integrate_world_bias_impulse_(s, Ts, omega, D);
            Lb.template block<12,3>(OFF_V, 0) = -D;
            Qd.noalias() += wq * (Lb * Q_bg_ * Lb.transpose());
        }

        // The OU acceleration driving noise is already available in exact
        // discrete form for the linear chain and is independent of gyro and
        // gyro-bias driving noise.
        if (linear_block_enabled_) add_ou_process_noise_(Ts, Qd);

        if constexpr (with_accel_bias) {
            if (acc_bias_updates_enabled_) {
                // A body-frame OU noise impulse decays in body coordinates;
                // its final Level-2 coordinate is simply R_end times that
                // decayed body increment, so the scalar OU integral is exact.
                const T qint = T(0.5) * tau_ba_ * (-std::expm1(T(-2)*Ts/tau_ba_));
                Qd.template block<3,3>(OFF_BA, OFF_BA) +=
                    (two_frame_bias ? Matrix3(R_end * Q_ba_ * R_end.transpose()) : Q_ba_) * qint;
            }
        }

        Qd = T(0.5) * (Qd + Qd.transpose()).eval();
        ou_detail::project_psd_ou_iii<T,NX>(Qd);
    }

    [[nodiscard]] Tangent error_from(const Kalman3D_Wave_TFG& reference) const {
        Tangent xi;
        xi.setZero();
        const Vector3 phi = ocean_imu::lie::Log<T>(Matrix3(X_.R * reference.X_.R.transpose()));
        const Matrix3 dR  = ocean_imu::lie::Exp<T>(phi);
        const Matrix3 Jli = ocean_imu::lie::inv_left_jacobian<T>(phi);
        xi.template segment<3>(OFF_PHI) = phi;
        for (int c = 0; c < 4; ++c)
            xi.template segment<3>(OFF_V + 3*c) =
                Jli * (X_.X.col(c) - dR * reference.X_.X.col(c));

        if constexpr (two_frame_bias) {
            const Matrix3 Jri = ocean_imu::lie::inv_left_jacobian<T>(Vector3(-phi));
            const Matrix3 JriR = Jri * reference.X_.R;
            xi.template segment<3>(OFF_BG) = JriR * (X_.B.col(0) - reference.X_.B.col(0));
            if constexpr (with_accel_bias)
                xi.template segment<3>(OFF_BA) = JriR * (X_.B.col(1) - reference.X_.B.col(1));
        } else {
            xi.template segment<3>(OFF_BG) = X_.B.col(0) - reference.X_.B.col(0);
            if constexpr (with_accel_bias)
                xi.template segment<3>(OFF_BA) = X_.B.col(1) - reference.X_.B.col(1);
        }
        return xi;
    }

    void time_update(const Vector3& gyro_meas, T Ts) {
        if (!(Ts > T(0)) || !std::isfinite(Ts)) return;
        MatrixNX Phi, Qd;
        build_transition(Ts, gyro_meas, Phi);
        build_process_noise(Ts, gyro_meas, Qd);
        propagate_mean(gyro_meas, Ts);
        P_ = (Phi * P_ * Phi.transpose()).eval();
        P_ += Qd;
        apply_pending_aw_covariance_inflation_();
        P_ = T(0.5) * (P_ + P_.transpose()).eval();
    }

    void inject(const Tangent& xi) {
        const Vector3 phi = xi.template segment<3>(OFF_PHI);
        const Matrix3 E_R = ocean_imu::lie::Exp<T>(phi);
        const Matrix3 Jl  = ocean_imu::lie::left_jacobian<T>(phi);
        const Matrix3 R_pre = X_.R;

        Eigen::Matrix<T,3,4> rho;
        for (int c = 0; c < 4; ++c) rho.col(c) = xi.template segment<3>(OFF_V + 3*c);
        X_.X = (E_R * X_.X + Jl * rho).eval();
        X_.R = (E_R * X_.R).eval();
        reorthonormalize_();

        if constexpr (two_frame_bias) {
            const Matrix3 Jr = ocean_imu::lie::left_jacobian<T>(Vector3(-phi));
            const Matrix3 RtJr = R_pre.transpose() * Jr;
            X_.B.col(0) += RtJr * xi.template segment<3>(OFF_BG);
            if constexpr (with_accel_bias)
                X_.B.col(1) += RtJr * xi.template segment<3>(OFF_BA);
        } else {
            X_.B.col(0) += xi.template segment<3>(OFF_BG);
            if constexpr (with_accel_bias) X_.B.col(1) += xi.template segment<3>(OFF_BA);
        }
    }

    // Full finite-retraction reset Jacobian. For Level 2 this is the right
    // Jacobian of the complete TwoFrameGroup. Level-1 additive bias blocks stay
    // identity while attitude/world blocks retain the group retraction map.
    void build_reset_jacobian(const Tangent& correction, MatrixNX& J) const {
        J.setIdentity();
        const Vector3 phi = correction.template segment<3>(OFF_PHI);
        const Matrix3 Jr = ocean_imu::lie::left_jacobian<T>(Vector3(-phi));
        J.template block<3,3>(OFF_PHI, OFF_PHI) = Jr;
        for (int c = 0; c < 4; ++c)
            set_reset_vector_block_(J, OFF_V + 3*c, phi,
                                    correction.template segment<3>(OFF_V + 3*c), Jr);
        if constexpr (two_frame_bias) {
            set_reset_vector_block_(J, OFF_BG, phi,
                                    correction.template segment<3>(OFF_BG), Jr);
            if constexpr (with_accel_bias)
                set_reset_vector_block_(J, OFF_BA, phi,
                                        correction.template segment<3>(OFF_BA), Jr);
        }
    }

    // A magnetic first lock is a world-frame gauge choice, not a heading
    // measurement update. Rotate every world-referred nominal quantity and the
    // complete tangent covariance coherently. Body biases remain body-frame
    // states; their Level-2 tangent coordinates rotate because beta=R delta_b.
    void apply_world_yaw_gauge(T yaw_delta) {
        if (!std::isfinite(yaw_delta)) return;
        const Matrix3 Q = ocean_imu::lie::Exp<T>(Vector3(T(0), T(0), yaw_delta));
        X_.R = (Q * X_.R).eval();
        X_.X = (Q * X_.X).eval();
        reorthonormalize_();

        MatrixNX G = MatrixNX::Identity();
        G.template block<3,3>(OFF_PHI, OFF_PHI) = Q;
        for (int c = 0; c < 4; ++c)
            G.template block<3,3>(OFF_V + 3*c, OFF_V + 3*c) = Q;
        if constexpr (two_frame_bias) {
            G.template block<3,3>(OFF_BG, OFF_BG) = Q;
            if constexpr (with_accel_bias) G.template block<3,3>(OFF_BA, OFF_BA) = Q;
        }
        P_ = (G * P_ * G.transpose()).eval();
        P_ = T(0.5) * (P_ + P_.transpose()).eval();

        Sigma_aw_ = (Q * Sigma_aw_ * Q.transpose()).eval();
        R_S_ = (Q * R_S_ * Q.transpose()).eval();
        if (have_mag_reference_) B_w_ = Q * B_w_;
    }

    struct MeasDiag3 {
        Vector3 r{Vector3::Zero()};
        Matrix3 S{Matrix3::Zero()};
        T nis{T(0)};
        bool accepted{false};
    };

    [[nodiscard]] const MeasDiag3& lastAccDiag() const { return last_acc_diag_; }
    [[nodiscard]] const MeasDiag3& lastMagDiag() const { return last_mag_diag_; }
    [[nodiscard]] const MeasDiag3& lastIntegralDiag() const { return last_S_diag_; }

    void set_Racc_std(const Vector3& std_body) { Racc_ = std_body.cwiseAbs2().asDiagonal(); }
    void set_Racc_matrix(const Matrix3& R_body) { Racc_ = T(0.5)*(R_body + R_body.transpose()); }
    void set_Rmag_std(const Vector3& std_body)  { Rmag_ = std_body.cwiseAbs2().asDiagonal(); }
    void set_Rmag(T variance) { Rmag_ = Matrix3::Identity() * std::abs(variance); }
    void set_RS_noise(const Vector3& sigma_S) {
        R_S_ = sigma_S.array().square().matrix().asDiagonal();
    }
    void set_RS_noise(T sigma_S) { set_RS_noise(Vector3::Constant(std::abs(sigma_S))); }
    void set_RS_noise_matrix(const Matrix3& R_S) { R_S_ = T(0.5)*(R_S + R_S.transpose()); }

    void set_magnetic_reference_world(const Vector3& B_w) {
        if (B_w.allFinite() && B_w.norm() > T(1e-9)) {
            B_w_ = B_w;
            have_mag_reference_ = true;
        }
    }
    [[nodiscard]] const Vector3& magnetic_reference_world() const { return B_w_; }
    [[nodiscard]] bool has_magnetic_reference() const { return have_mag_reference_; }

    void set_accel_bias_temp_coeff(const Vector3& k_a) { k_a_ = k_a; }
    void set_acc_bias_updates_enabled(bool on) {
        if (acc_bias_updates_enabled_ == on) return;
        if constexpr (with_accel_bias) {
            if (!on) {
                // A frozen bias must not be changed indirectly by another
                // measurement through stale cross-covariance.
                for (int i = 0; i < NX; ++i) {
                    if (i >= OFF_BA && i < OFF_BA+3) continue;
                    for (int k = 0; k < 3; ++k) {
                        P_(i, OFF_BA+k) = T(0);
                        P_(OFF_BA+k, i) = T(0);
                    }
                }
            } else {
                Matrix3 Pba = T(0.5) *
                    (P_.template block<3,3>(OFF_BA, OFF_BA) +
                     P_.template block<3,3>(OFF_BA, OFF_BA).transpose()).eval();
                const T floor = sigma_ba0_ * sigma_ba0_;
                for (int i = 0; i < 3; ++i) Pba(i,i) = std::max(Pba(i,i), floor);
                P_.template block<3,3>(OFF_BA, OFF_BA) = Pba;
            }
        }
        acc_bias_updates_enabled_ = on;
    }
    [[nodiscard]] bool acc_bias_updates_enabled() const { return acc_bias_updates_enabled_; }

    void set_meas_gate_nis(T threshold) { nis_gate_ = std::max(T(0), threshold); }
    [[nodiscard]] Vector3 gravity_world() const { return Vector3(T(0), T(0), gravity_magnitude_); }

    void accel_residual(const Vector3& acc_meas_body, T tempC,
                        Vector3& r, Eigen::Matrix<T,3,NX>& H, Matrix3& Rw) const {
        const Vector3 ba = accel_bias_at_(tempC);
        r = X_.R * (acc_meas_body - ba) - (X_.X.col(3) - gravity_world());
        H.setZero();
        H.template block<3,3>(0, OFF_PHI) = ocean_imu::lie::skew<T>(gravity_world());
        if (linear_block_enabled_) H.template block<3,3>(0, OFF_AW) = -Matrix3::Identity();
        if constexpr (with_accel_bias) {
            if (acc_bias_updates_enabled_)
                H.template block<3,3>(0, OFF_BA) =
                    two_frame_bias ? Matrix3(-Matrix3::Identity()) : Matrix3(-X_.R);
        }
        Rw = X_.R * Racc_ * X_.R.transpose();
        if (!linear_block_enabled_) Rw += Sigma_aw_;
    }

    void mag_residual(const Vector3& mag_meas_body,
                      Vector3& r, Eigen::Matrix<T,3,NX>& H, Matrix3& Rw) const {
        r = X_.R * mag_meas_body - B_w_;
        H.setZero();
        H.template block<3,3>(0, OFF_PHI) = -ocean_imu::lie::skew<T>(B_w_);
        Rw = X_.R * Rmag_ * X_.R.transpose();
    }

    void integral_residual(Vector3& r, Eigen::Matrix<T,3,NX>& H, Matrix3& Rw) const {
        r = -X_.X.col(2);
        H.setZero();
        H.template block<3,3>(0, OFF_PHI) = ocean_imu::lie::skew<T>(Vector3(X_.X.col(2)));
        H.template block<3,3>(0, OFF_S) = -Matrix3::Identity();
        Rw = R_S_;
    }

    bool measurement_update_acc_only(const Vector3& acc_meas_body, T tempC = T(35)) {
        if (!acc_meas_body.allFinite()) { last_acc_diag_ = MeasDiag3{}; return false; }
        Vector3 r; Eigen::Matrix<T,3,NX> H; Matrix3 Rw;
        accel_residual(acc_meas_body, tempC, r, H, Rw);
        return apply_update3_(r, H, Rw, last_acc_diag_);
    }

    bool measurement_update_mag_only(const Vector3& mag_meas_body) {
        if (!have_mag_reference_ || !mag_meas_body.allFinite() ||
            !(mag_meas_body.norm() > T(1e-9))) {
            last_mag_diag_ = MeasDiag3{};
            return false;
        }
        Vector3 r; Eigen::Matrix<T,3,NX> H; Matrix3 Rw;
        mag_residual(mag_meas_body, r, H, Rw);
        return apply_update3_(r, H, Rw, last_mag_diag_);
    }

    bool applyIntegralZeroPseudoMeas() {
        if (!linear_block_enabled_) { last_S_diag_ = MeasDiag3{}; return false; }
        Vector3 r; Eigen::Matrix<T,3,NX> H; Matrix3 Rw;
        integral_residual(r, H, Rw);
        return apply_update3_(r, H, Rw, last_S_diag_);
    }

    bool apply_update3(const Vector3& r, const Eigen::Matrix<T,3,NX>& H,
                       const Matrix3& Rw, MeasDiag3& diag) {
        return apply_update3_(r, H, Rw, diag);
    }

  private:
    void world_transition_(T dt, Eigen::Matrix<T,4,4>& Phi) const {
        Phi.setIdentity();
        if (linear_block_enabled_)
            ou_detail::IntegratedOUChain<T,3>::transition(tau_aw_, dt, Phi);
    }

    // D(s,h) = int_s^h Phi_A(h-t) [x(t)]x R(t) dt.  It is the physical
    // body-gyro-bias impulse map into [rho_v rho_p rho_S rho_a] at h.  Keeping
    // it in body-bias coordinates makes the same function valid for Level 1
    // and Level 2; only the endpoint bias coordinate differs.
    void integrate_world_bias_impulse_(T s, T h, const Vector3& omega,
                                       Eigen::Matrix<T,12,3>& D) const {
        D.setZero();
        if (!(h > s)) return;
        static constexpr double nodes[5] = {
            -0.9061798459386640, -0.5384693101056831, 0.0,
             0.5384693101056831,  0.9061798459386640
        };
        static constexpr double weights[5] = {
            0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
            0.4786286704993665, 0.2369268850561891
        };

        const T span = h - s;
        for (int q = 0; q < 5; ++q) {
            const T t = s + T(0.5) * span * (T(1) + T(nodes[q]));
            const T wq = T(0.5) * span * T(weights[q]);
            Eigen::Matrix<T,4,4> Phi_t, Phi_rem;
            world_transition_(t, Phi_t);
            world_transition_(h - t, Phi_rem);
            const Eigen::Matrix<T,3,4> X_t = X_.X * Phi_t.transpose();
            const Matrix3 R_t = X_.R * ocean_imu::lie::Exp<T>(Vector3(omega * t));

            for (int i = 0; i < 4; ++i) {
                Matrix3 M = Matrix3::Zero();
                for (int j = 0; j < 4; ++j) {
                    if (Phi_rem(i,j) == T(0)) continue;
                    M.noalias() += Phi_rem(i,j) *
                        (ocean_imu::lie::skew<T>(Vector3(X_t.col(j))) * R_t);
                }
                D.template block<3,3>(3*i, 0).noalias() += wq * M;
            }
        }
    }

    [[nodiscard]] Vector3 accel_bias_at_(T tempC) const {
        if constexpr (with_accel_bias)
            return X_.B.col(1) + k_a_ * (tempC-kTempRefC);
        (void)tempC;
        return Vector3::Zero();
    }

    // Off-diagonal vector/attitude block of the full group right Jacobian.
    // For ad_xi = [Omega 0; Z Omega], J_r = sum (-ad)^n/(n+1)!.
    static Matrix3 reset_cross_block_(const Vector3& phi, const Vector3& z) {
        const Matrix3 Omega = ocean_imu::lie::skew<T>(phi);
        const Matrix3 Z = ocean_imu::lie::skew<T>(z);
        Matrix3 D = Matrix3::Identity();
        Matrix3 L = Matrix3::Zero();
        Matrix3 Q = Matrix3::Zero();
        T coeff = T(1);
        const T tol = T(8) * std::numeric_limits<T>::epsilon();
        for (int n = 1; n <= 32; ++n) {
            const Matrix3 Lnext = L * Omega + D * Z;
            const Matrix3 Dnext = D * Omega;
            L = Lnext;
            D = Dnext;
            coeff *= -T(1) / T(n+1);
            const Matrix3 term = coeff * L;
            Q += term;
            if (n >= 5 && term.cwiseAbs().maxCoeff() <=
                              tol * (T(1) + Q.cwiseAbs().maxCoeff())) break;
        }
        return Q;
    }

    static void set_reset_vector_block_(MatrixNX& J, int off,
                                        const Vector3& phi, const Vector3& z,
                                        const Matrix3& Jr) {
        J.template block<3,3>(off, off) = Jr;
        J.template block<3,3>(off, OFF_PHI) = reset_cross_block_(phi, z);
    }

    bool apply_update3_(const Vector3& r, const Eigen::Matrix<T,3,NX>& H,
                        const Matrix3& Rw, MeasDiag3& diag) {
        diag = MeasDiag3{};
        diag.r = r;
        if (!r.allFinite() || !H.allFinite() || !Rw.allFinite()) return false;

        const Eigen::Matrix<T,NX,3> PHt = P_ * H.transpose();
        Matrix3 S = H * PHt + Rw;
        S = T(0.5) * (S + S.transpose()).eval();
        diag.S = S;

        Eigen::LDLT<Matrix3> ldlt(S);
        if (ldlt.info() != Eigen::Success) return false;
        const Vector3 Sinv_r = ldlt.solve(r);
        if (!Sinv_r.allFinite()) return false;
        diag.nis = r.dot(Sinv_r);
        if (nis_gate_ > T(0) && !(diag.nis <= nis_gate_)) return false;

        const Eigen::Matrix<T,NX,3> K = ldlt.solve(PHt.transpose()).transpose();
        if (!K.allFinite()) return false;

        const Tangent correction = K * r;
        const MatrixNX IKH = MatrixNX::Identity() - K * H;
        const MatrixNX Pj =
            (IKH * P_ * IKH.transpose() + K * Rw * K.transpose()).eval();

        MatrixNX Jreset;
        build_reset_jacobian(correction, Jreset);
        inject(Tangent(-correction));
        P_ = (Jreset * Pj * Jreset.transpose()).eval();
        P_ = T(0.5) * (P_ + P_.transpose()).eval();

        diag.accepted = true;
        return true;
    }

    void reorthonormalize_() {
        Eigen::Quaternion<T> q(X_.R);
        q.normalize();
        X_.R = q.toRotationMatrix();
    }

    // Add only the positive-semidefinite part of the shortfall against the
    // stationary target, so the sync can never remove information the
    // measurements put there and can never break the PSD of P.
    void apply_pending_aw_covariance_inflation_() {
        if (!aw_cov_floor_pending_ || !linear_block_enabled_) return;
        aw_cov_floor_pending_ = false;

        Matrix3 P_aw = P_.template block<3,3>(OFF_AW, OFF_AW);
        P_aw = T(0.5) * (P_aw + P_aw.transpose());
        Matrix3 Delta = aw_cov_floor_target_ - P_aw;
        Delta = T(0.5) * (Delta + Delta.transpose());

        Eigen::SelfAdjointEigenSolver<Matrix3> es(Delta);
        if (es.info() != Eigen::Success) return;
        Vector3 evals = es.eigenvalues();
        for (int i = 0; i < 3; ++i) evals(i) = std::max(T(0), evals(i));
        Delta = es.eigenvectors() * evals.asDiagonal() * es.eigenvectors().transpose();
        Delta = T(0.5) * (Delta + Delta.transpose());

        P_.template block<3,3>(OFF_AW, OFF_AW) += Delta;
    }

    void add_ou_process_noise_(T Ts, MatrixNX& Qd) const {
        const bool correlated = !ou_detail::is_isotropic3<T>(Sigma_aw_) &&
                                !is_diagonal_(Sigma_aw_);
        if (correlated) {
            Eigen::Matrix<T,4,4> Q_unit;
            ou_detail::IntegratedOUChain<T,3>::process_covariance(tau_aw_, Ts, T(1), Q_unit);
            Q_unit = T(0.5) * (Q_unit + Q_unit.transpose()).eval();
            const Matrix3 Sig = T(0.5) * (Sigma_aw_ + Sigma_aw_.transpose());
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    Qd.template block<3,3>(OFF_V + 3*i, OFF_V + 3*j).noalias()
                        += Sig * Q_unit(i,j);
        } else {
            for (int axis = 0; axis < 3; ++axis) {
                Eigen::Matrix<T,4,4> Q_axis;
                ou_detail::IntegratedOUChain<T,3>::process_covariance(
                    tau_aw_, Ts, Sigma_aw_(axis,axis), Q_axis);
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j)
                        Qd(OFF_V + 3*i + axis, OFF_V + 3*j + axis) += Q_axis(i,j);
            }
        }
    }

    static bool is_diagonal_(const Matrix3& M) {
        Matrix3 off = M;
        off.diagonal().setZero();
        return off.cwiseAbs().maxCoeff() <=
               T(1e-12) * std::max(T(1), M.cwiseAbs().maxCoeff());
    }

    T gravity_magnitude_;
    Group X_{};
    MatrixNX P_{MatrixNX::Identity()};

    T tau_aw_{T(1.0)};
    Matrix3 Sigma_aw_{Matrix3::Identity()};
    Matrix3 aw_cov_floor_target_{Matrix3::Identity()};
    bool aw_cov_floor_pending_{false};
    Matrix3 Q_gyro_{Matrix3::Identity() * T(2.5e-5)};
    Matrix3 Q_bg_{Matrix3::Identity() * T(1e-10)};
    Matrix3 Q_ba_{Matrix3::Identity() * T(2.5e-7)};
    T tau_ba_{T(5000.0)};
    T sigma_ba0_{T(0.004)};

    Matrix3 Racc_{Matrix3::Identity() * T(0.16)};
    Matrix3 Rmag_{Matrix3::Identity() * T(0.01)};
    Matrix3 R_S_{Matrix3::Identity() * T(1.5)};

    Vector3 B_w_{Vector3::UnitX()};
    bool have_mag_reference_{false};

    Vector3 k_a_{Vector3::Constant(T(0.002))};
    bool acc_bias_updates_enabled_{true};
    T nis_gate_{T(0)};
    bool linear_block_enabled_{true};

    static constexpr T kTempRefC = T(35);

    MeasDiag3 last_acc_diag_{};
    MeasDiag3 last_mag_diag_{};
    MeasDiag3 last_S_diag_{};
    Vector3 last_omega_{Vector3::Zero()};
};

} // namespace ocean_imu::kalman
