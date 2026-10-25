#pragma once

/*
  Copyright (c) 2025-2026 Mikhail Grushinskiy

  Based on: https://github.com/thomaspasser/q-mekf
  MIT License, Copyright (c) 2023 Thomas Passer

  Kalman3D_Wave_CARMA -- EXPERIMENTAL research filter, not a shipping
  estimator.  Kalman3D_Wave_OU_III remains the deployed and authoritative
  implementation; this class exists to measure whether a CARMA vertical prior
  improves heave estimation against it.  See docs/carma-heave-filter.md.

  State (NED, z down):
     dtheta (3)          attitude error, left-multiplicative WORLD->BODY
     b_g (3)             gyro bias (optional)
     horizontal (8)      [v_x, v_y, p_x, p_y, S_x, S_y, a_x, a_y]: the OU-III
                         integrated OU chain and S = 0 pseudo-measurement,
                         per horizontal axis, unchanged
     vertical (N [+1])   CARMA(N, N-3) prior in [h, v_z, a_z, latent ...]
                         coordinates (CarmaVerticalModel.h), with h = p_z;
                         optionally S_z = int h dt (with_vertical_integral)
     b_a (3)             residual accelerometer bias OU (optional)

  For CARMA(5,2) without S_z this is 6 + 8 + 5 + 3 = 22 states, against the
  21 of OU-III.

  Retained from OU-III, line for line in their mathematics: the attitude and
  gyro-bias prediction (exact rotation transition, structured Qd), the
  accelerometer and magnetometer measurement models and Jacobians, the
  error-state reset, the residual accelerometer-bias OU with its measurement
  freeze, the horizontal OU chain and its S pseudo-measurement, the queued
  PSD a_w covariance floor (horizontal axes only), Joseph covariance updates
  and the input guards.  Not carried over: wind heel, IMU lever arm and the
  commanded-parameter random-walk filter, none of which the deployed wrapper
  enables.

  The vertical acceleration enters the accelerometer model exactly where
  OU-III's a_w,z does,
      f_b = R_wb (a_w - g) + b_a,   a_w = [a_x, a_y, a_z],  a_z = x_vertical[2],
  so the Jacobian with respect to the vertical block is R_wb e_z e_a^T.
  The vertical block has no zero pseudo-measurement unless
  with_vertical_integral is set; its stationary prior bounds the displacement
  instead.  With with_vertical_integral the explicit primitive S_z' = h is
  propagated exactly together with the CARMA states and S_z = 0 is applied at
  the same cadence as the horizontal S update, with the variance law chosen by
  the wrapper.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

#include "kalman_ou_common/KalmanOUCoreMath.h"
#include "kalman_carma/CarmaVerticalModel.h"

namespace ocean_imu::kalman {

template <typename T = float, int NCARMA = 5, bool with_vertical_integral = false,
          bool with_gyro_bias = true, bool with_accel_bias = true>
class Kalman3D_Wave_CARMA {
public:
    static constexpr int BASE_N = with_gyro_bias ? 6 : 3;
    static constexpr int NH = 8;
    static constexpr int NZ = NCARMA + (with_vertical_integral ? 1 : 0);
    static constexpr int NL = NH + NZ;
    static constexpr int NX = BASE_N + NL + (with_accel_bias ? 3 : 0);

    static constexpr int OFF_L  = BASE_N;
    static constexpr int OFF_VH = BASE_N + 0;   // v_x, v_y
    static constexpr int OFF_PH = BASE_N + 2;   // p_x, p_y
    static constexpr int OFF_SH = BASE_N + 4;   // S_x, S_y
    static constexpr int OFF_AH = BASE_N + 6;   // a_x, a_y
    static constexpr int OFF_Z  = BASE_N + 8;   // h, v_z, a_z, latent...
    static constexpr int OFF_SZ = with_vertical_integral ? (OFF_Z + NCARMA) : -1;
    static constexpr int OFF_BA = with_accel_bias ? (BASE_N + NL) : -1;

    using Vector3 = Eigen::Matrix<T,3,1>;
    using Matrix3 = Eigen::Matrix<T,3,3>;
    using MatrixBaseN = Eigen::Matrix<T,BASE_N,BASE_N>;
    using MatrixNX = Eigen::Matrix<T,NX,NX>;
    using MatrixNX3 = Eigen::Matrix<T,NX,3>;
    using MatrixLL = Eigen::Matrix<T,NL,NL>;
    using MatrixZZ = Eigen::Matrix<T,NZ,NZ>;
    using MatrixC = Eigen::Matrix<T,NCARMA,NCARMA>;
    using VerticalModel = carma::CarmaVerticalModel<T,NCARMA>;

    static constexpr T STD_GRAVITY = T(9.80665);
    static constexpr T tempC_ref = T(35.0);

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    struct MeasDiag3 {
        Vector3 r = Vector3::Zero();
        Matrix3 S = Matrix3::Zero();
        T nis = std::numeric_limits<T>::quiet_NaN();
        bool accepted = false;
    };

    const MeasDiag3& lastAccDiag() const noexcept { return last_acc_diag_; }
    const MeasDiag3& lastMagDiag() const noexcept { return last_mag_diag_; }

    Kalman3D_Wave_CARMA(Vector3 const& sigma_a, Vector3 const& gyro_noise_density_rad_sqrt_s,
                        Vector3 const& sigma_m,
                        T Pq0 = T(5e-4), T Pb0 = T(1e-6), T b0 = T(1e-11), T R_S_noise_var = T(1.5),
                        T gravity_magnitude = T(STD_GRAVITY),
                        const carma::CarmaShape& shape = carma::CarmaDefaultShape<NCARMA>::value);

    bool initialize_from_acc_mag(Vector3 const& acc, Vector3 const& mag);
    bool initialize_from_acc(Vector3 const& acc);
    bool initialize_from_acc_preserve_yaw(Vector3 const& acc);
    bool initialize_from_attitude(Eigen::Quaternion<T> const& q_bw, T tilt_sigma_rad, T yaw_sigma_rad);
    static Eigen::Quaternion<T> quaternion_from_acc(Vector3 const& acc);

    void set_mag_world_ref(const Vector3& B_world) { v2ref = B_world; }

    void time_update(Vector3 const& gyr, T Ts);
    void measurement_update_acc_only(Vector3 const& acc, T tempC = tempC_ref);
    void measurement_update_mag_only(Vector3 const& mag);

    // Horizontal S_x, S_y = 0 (and S_z = 0 with with_vertical_integral).
    void applyIntegralZeroPseudoMeas();

    void measurement_update_position_pseudo(const Vector3& p_meas, const Vector3& sigma_meas);
    void measurement_update_velocity_pseudo(const Vector3& v_meas, const Vector3& sigma_meas);
    void measurement_update_vert_velocity_pseudo(T vz_meas, T sigma_meas);

    [[nodiscard]] Eigen::Quaternion<T> quaternion() const { return qref.conjugate(); }
    [[nodiscard]] Eigen::Quaternion<T> quaternion_boat() const { return quaternion(); }
    [[nodiscard]] MatrixBaseN covariance_base() const { return Pext.template topLeftCorner<BASE_N,BASE_N>(); }
    [[nodiscard]] const MatrixNX& covariance_full() const { return Pext; }

    void set_quaternion_boat(const Eigen::Quaternion<T>& q_bw) {
        project_gyro_bias_();
        Eigen::Quaternion<T> q = q_bw;
        if (!(q.norm() > T(1e-8))) return;
        q.normalize();
        qref = q.conjugate();
        qref.normalize();
        xext.template segment<3>(0).setZero();
    }

    [[nodiscard]] Vector3 gyroscope_bias() const {
        if constexpr (with_gyro_bias) return xext.template segment<3>(3);
        else return Vector3::Zero();
    }
    [[nodiscard]] Vector3 get_acc_bias() const {
        if constexpr (with_accel_bias) return xext.template segment<3>(OFF_BA);
        else return Vector3::Zero();
    }

    bool acc_bias_updates_enabled() const { return acc_bias_updates_enabled_; }
    void set_acc_bias_updates_enabled(bool en);

    [[nodiscard]] Vector3 get_velocity() const {
        return Vector3(xext(OFF_VH), xext(OFF_VH + 1), xext(OFF_Z + 1));
    }
    [[nodiscard]] Vector3 get_position() const {
        return Vector3(xext(OFF_PH), xext(OFF_PH + 1), xext(OFF_Z));
    }
    [[nodiscard]] Vector3 get_world_accel() const {
        return Vector3(xext(OFF_AH), xext(OFF_AH + 1), xext(OFF_Z + 2));
    }
    [[nodiscard]] Vector3 get_integral_displacement() const {
        if constexpr (with_vertical_integral)
            return Vector3(xext(OFF_SH), xext(OFF_SH + 1), xext(OFF_SZ));
        else
            return Vector3(xext(OFF_SH), xext(OFF_SH + 1), std::numeric_limits<T>::quiet_NaN());
    }
    // Vertical CARMA block mean [h, v_z, a_z, latent...] and covariance.
    [[nodiscard]] Eigen::Matrix<T,NZ,1> vertical_state() const { return xext.template segment<NZ>(OFF_Z); }
    [[nodiscard]] MatrixZZ vertical_covariance() const { return Pext.template block<NZ,NZ>(OFF_Z, OFF_Z); }
    // [p_z, v_z] covariance, the block scored by vertical NEES.
    [[nodiscard]] Eigen::Matrix<T,2,2> vertical_pv_covariance() const {
        return Pext.template block<2,2>(OFF_Z, OFF_Z);
    }

    // Horizontal OU chain: active parameters.
    [[nodiscard]] T get_aw_time_constant() const { return tau_aw; }
    [[nodiscard]] Vector3 get_aw_stationary_std() const {
        // z reports the CARMA stationary acceleration std, for symmetry with
        // OU-III's getter.
        return Vector3(std::sqrt(std::max(T(0), Sigma_aw_h(0))),
                       std::sqrt(std::max(T(0), Sigma_aw_h(1))),
                       vmodel_.accel_std());
    }
    [[nodiscard]] Vector3 get_RS_noise_std() const {
        return R_S.diagonal().cwiseMax(T(0)).cwiseSqrt();
    }

    void set_aw_time_constant(T tau_seconds) { tau_aw = std::max(T(1e-3), tau_seconds); }
    void set_pseudo_update_period_s(T period_s) {
        if (!(period_s > T(0)) || !std::isfinite(period_s)) return;
        const T new_period = std::max(T(1e-4), period_s);
        pseudo_update_elapsed_s_ =
            ou_detail::retarget_period_elapsed_progress_preserving(pseudo_update_elapsed_s_, new_period);
        pseudo_update_period_s_ = new_period;
    }
    [[nodiscard]] T get_pseudo_update_period_s() const { return pseudo_update_period_s_; }

    // Horizontal OU stationary std; the z component is ignored (the vertical
    // prior is configured by set_vertical_carma()).
    void set_aw_stationary_std(const Vector3& std_aw) {
        Sigma_aw_h(0) = std_aw.x() * std_aw.x();
        Sigma_aw_h(1) = std_aw.y() * std_aw.y();
    }

    // Horizontal a_w marginal reset plus vertical reset to the CARMA
    // stationary covariance; called by the shared initialization and at Live.
    void reset_aw_covariance_to_stationary() {
        for (int k = 0; k < 2; ++k) {
            const int a = OFF_AH + k;
            for (int i = 0; i < NX; ++i) {
                if (i == a) continue;
                Pext(a, i) = T(0);
                Pext(i, a) = T(0);
            }
            Pext(a, a) = Sigma_aw_h(k);
        }
        aw_covariance_floor_pending_ = false;
        reset_vertical_covariance_to_stationary();
        symmetrize_Pext_();
    }

    // Queued PSD a_w floor at the next prediction, horizontal axes; see the
    // OU-III counterpart.  The vertical CARMA block has no such operation.
    void synchronize_aw_covariance_to_stationary() {
        aw_covariance_floor_target_ = Sigma_aw_h;
        aw_covariance_floor_pending_ = true;
    }

    // Configure the vertical prior: kappa multiplies every frequency of the
    // shape, accel_std is the stationary std of a_z [m/s^2].  Parameters are
    // replaced; the state and covariance are kept, so the change reaches the
    // estimate through F and Q, as an OU-III tau or sigma change does.
    bool set_vertical_carma(T kappa, T accel_std) { return vmodel_.configure(kappa, accel_std); }
    [[nodiscard]] const VerticalModel& vertical_model() const noexcept { return vmodel_; }

    void reset_vertical_covariance_to_stationary() {
        MatrixC P0;
        if (!vmodel_.stationary_covariance(P0)) return;
        for (int i = 0; i < NX; ++i) {
            for (int j = 0; j < NZ; ++j) {
                if (i >= OFF_Z && i < OFF_Z + NZ) continue;
                Pext(i, OFF_Z + j) = T(0);
                Pext(OFF_Z + j, i) = T(0);
            }
        }
        Pext.template block<NZ,NZ>(OFF_Z, OFF_Z).setZero();
        Pext.template block<NCARMA,NCARMA>(OFF_Z, OFF_Z) = P0;
        if constexpr (with_vertical_integral) {
            Pext(OFF_SZ, OFF_SZ) = sigma_S0_ * sigma_S0_;
        }
        symmetrize_Pext_();
    }

    // R_S per axis as a standard deviation [m s]: x, y horizontal; z used by
    // the S_z variant only.
    void set_RS_noise(const Vector3& sigma_S) {
        R_S = sigma_S.array().square().matrix().asDiagonal();
    }
    void set_Racc_std(const Vector3& sigma_acc) {
        Racc = sigma_acc.array().square().matrix().asDiagonal();
    }
    void set_Rmag(const Vector3& sigma_mag) {
        Rmag = sigma_mag.array().square().matrix().asDiagonal();
    }

    void set_initial_acc_bias_std(T s) {
        if constexpr (with_accel_bias) {
            sigma_bacc0_ = std::max(T(0), s);
            Pext.template block<3,3>(OFF_BA, OFF_BA) = Matrix3::Identity() * sigma_bacc0_ * sigma_bacc0_;
        }
    }
    void set_Q_bacc_rw(const Vector3& driving_std_per_sqrt_s) {
        if constexpr (with_accel_bias)
            Q_bacc_ = driving_std_per_sqrt_s.array().square().matrix().asDiagonal();
    }
    void set_acc_bias_time_constant(T tau_seconds) {
        if constexpr (with_accel_bias) tau_bacc_ = std::max(T(1e-3), tau_seconds);
    }
    [[nodiscard]] T get_acc_bias_time_constant() const noexcept { return tau_bacc_; }
    [[nodiscard]] Vector3 get_acc_bias_at_temperature(T tempC) const {
        if constexpr (with_accel_bias) return xext.template segment<3>(OFF_BA) + k_a_ * (tempC - tempC_ref);
        return Vector3::Zero();
    }
    [[nodiscard]] Vector3 get_acc_bias_body_at_temperature(T tempC) const {
        return get_acc_bias_at_temperature(tempC);
    }
    void set_initial_acc_bias(const Vector3& b0) {
        if constexpr (with_accel_bias) xext.template segment<3>(OFF_BA) = b0;
    }
    void set_accel_bias_temp_coeff(const Vector3& ka_per_degC) { k_a_ = ka_per_degC; }
    void set_accel_bias_limit(T limit_mps2) { acc_bias_limit_ = limit_mps2; }
    [[nodiscard]] T accel_bias_limit() const { return acc_bias_limit_; }

    void set_gyro_noise_density_rad_sqrt_s(const Vector3& density) {
        const Vector3 d = density.cwiseAbs();
        Qbase.template topLeftCorner<3,3>() = d.array().square().matrix().asDiagonal();
    }
    void set_exact_att_bias_Qd(bool on) { use_exact_att_bias_Qd_ = on; }

private:
    const T gravity_magnitude_ = T(STD_GRAVITY);

    Eigen::Quaternion<T> qref;
    Vector3 v2ref = Vector3::UnitX();

    Eigen::Matrix<T,NX,1> xext;
    MatrixNX Pext;

    Vector3 last_gyr_bias_corrected{};

    T sigma_bacc0_ = T(0.004);
    Matrix3 Q_bacc_ = Matrix3::Identity() * T(2.5e-7);
    T tau_bacc_ = T(5000.0);
    Vector3 k_a_ = Vector3::Constant(T(0.002));
    T acc_bias_limit_ = T(0.4);

    Matrix3 Rmag;
    MatrixBaseN Qbase;
    Matrix3 Racc;
    Matrix3 R_S;

    // Horizontal OU parameters.
    T tau_aw = T(2.1);
    Eigen::Matrix<T,2,1> Sigma_aw_h = Eigen::Matrix<T,2,1>::Constant(T(2.2 * 2.2));
    bool aw_covariance_floor_pending_ = false;
    Eigen::Matrix<T,2,1> aw_covariance_floor_target_ = Eigen::Matrix<T,2,1>::Zero();

    // Vertical prior and the cached exact discretization for the last step.
    VerticalModel vmodel_;
    MatrixZZ F_Z_cache_ = MatrixZZ::Identity();
    MatrixZZ Q_Z_cache_ = MatrixZZ::Zero();
    T cache_dt_ = T(-1);
    uint32_t cache_generation_ = 0xffffffffu;

    T sigma_S0_ = T(50.0);
    T pseudo_update_elapsed_s_ = T(0);
    T pseudo_update_period_s_ = T(0.015);

    bool acc_bias_updates_enabled_ = true;
    bool use_exact_att_bias_Qd_ = true;

    // Scratch (avoids large stack frames in the updates)
    MatrixBaseN F_AA_scratch_;
    MatrixBaseN Q_AA_scratch_;
    Eigen::Matrix<T,4,4> Phi_h_ = Eigen::Matrix<T,4,4>::Identity();
    Eigen::Matrix<T,4,4> Qd_h_[2];
    MatrixLL tmpLL_scratch_;
    MatrixBaseN tmpAA_scratch_;
    Eigen::Matrix<T,BASE_N,NL> tmpAL_scratch_;
    MatrixNX3 PCt_scratch_;
    MatrixNX3 K_scratch_;
    Matrix3 S_scratch_;

    MeasDiag3 last_acc_diag_;
    MeasDiag3 last_mag_diag_;

    void symmetrize_Pext_() {
        for (int i = 0; i < NX; ++i)
            for (int j = i + 1; j < NX; ++j) {
                const T v = T(0.5) * (Pext(i,j) + Pext(j,i));
                Pext(i,j) = v;
                Pext(j,i) = v;
            }
    }

    void apply_pending_aw_covariance_inflation_() {
        if (!aw_covariance_floor_pending_) return;
        aw_covariance_floor_pending_ = false;
        Eigen::Matrix<T,2,2> P_aw;
        P_aw << Pext(OFF_AH, OFF_AH), Pext(OFF_AH, OFF_AH + 1),
                Pext(OFF_AH + 1, OFF_AH), Pext(OFF_AH + 1, OFF_AH + 1);
        P_aw = T(0.5) * (P_aw + P_aw.transpose()).eval();
        // Delta = Pi_+(target - P_aw), the positive part in closed form for
        // a symmetric 2 x 2 [[a, b], [b, c]].
        const T a = aw_covariance_floor_target_(0) - P_aw(0,0);
        const T c = aw_covariance_floor_target_(1) - P_aw(1,1);
        const T b = -P_aw(0,1);
        const T mean = T(0.5) * (a + c);
        const T rad = std::sqrt(T(0.25) * (a - c) * (a - c) + b * b);
        const T l_max = mean + rad, l_min = mean - rad;
        T d00 = a, d01 = b, d11 = c;
        if (!(l_max > T(0))) {
            d00 = d01 = d11 = T(0);
        } else if (l_min < T(0)) {
            // l_max u u^T with u the unit eigenvector of l_max.
            T ux = b, uy = l_max - a;
            if (std::abs(ux) + std::abs(uy) < std::abs(l_max - c) + std::abs(b)) {
                ux = l_max - c;
                uy = b;
            }
            const T n2 = ux * ux + uy * uy;
            if (n2 > T(0)) {
                d00 = l_max * ux * ux / n2;
                d01 = l_max * ux * uy / n2;
                d11 = l_max * uy * uy / n2;
            } else {
                d00 = std::max(a, T(0));
                d01 = T(0);
                d11 = std::max(c, T(0));
            }
        }
        if (!(std::isfinite(d00) && std::isfinite(d01) && std::isfinite(d11))) return;
        Pext(OFF_AH, OFF_AH)         += d00;
        Pext(OFF_AH, OFF_AH + 1)     += d01;
        Pext(OFF_AH + 1, OFF_AH)      = Pext(OFF_AH, OFF_AH + 1);
        Pext(OFF_AH + 1, OFF_AH + 1) += d11;
    }

    Matrix3 R_wb() const { return qref.toRotationMatrix(); }
    Matrix3 R_bw() const { return qref.toRotationMatrix().transpose(); }

    void set_accel_only_attitude_covariance_(T tilt_sigma_rad = T(0.035), T yaw_sigma_rad = T(1.5708)) {
        Vector3 u = R_wb() * Vector3::UnitZ();
        const T n = u.norm();
        if (!(n > T(1e-8)) || !u.allFinite()) {
            Pext.template block<3,3>(0,0) = Matrix3::Identity() * (yaw_sigma_rad * yaw_sigma_rad);
            return;
        }
        u /= n;
        const Matrix3 Py = u * u.transpose();
        Matrix3 Patt = (tilt_sigma_rad * tilt_sigma_rad) * (Matrix3::Identity() - Py)
                     + (yaw_sigma_rad * yaw_sigma_rad) * Py;
        Pext.template block<3,3>(0,0) = T(0.5) * (Patt + Patt.transpose());
    }

    void project_gyro_bias_() {
        if constexpr (with_gyro_bias) {
            auto bias = xext.template segment<3>(3);
            ou_detail::project_gyro_bias<T>(bias);
        }
    }

    static Matrix3 skew_(const Vector3& v) {
        Matrix3 M;
        M << 0, -v(2), v(1),
             v(2), 0, -v(0),
            -v(1), v(0), 0;
        return M;
    }

    void applyQuaternionCorrectionFromErrorState();
    void project_acc_bias_();

    bool safe_ldlt3_(Matrix3& S, Eigen::LDLT<Matrix3>& ldlt, T noise_scale) const {
        ldlt.compute(S);
        if (ldlt.info() == Eigen::Success) return true;
        const T bump = std::max(std::numeric_limits<T>::epsilon(), T(1e-6) * (noise_scale + T(1)));
        S.diagonal().array() += bump;
        ldlt.compute(S);
        return ldlt.info() == Eigen::Success;
    }

    T nis3_(const Eigen::LDLT<Matrix3>& ldlt, const Vector3& r) const {
        const Vector3 x = ldlt.solve(r);
        if (!x.allFinite()) return std::numeric_limits<T>::quiet_NaN();
        const T v = r.dot(x);
        return std::isfinite(v) ? v : std::numeric_limits<T>::quiet_NaN();
    }

    void gain3_(const MatrixNX3& PCt, const Eigen::LDLT<Matrix3>& ldlt, MatrixNX3& K) const {
        for (int i = 0; i < NX; ++i) {
            const Vector3 rhs(PCt(i,0), PCt(i,1), PCt(i,2));
            const Vector3 sol = ldlt.solve(rhs);
            K(i,0) = sol(0); K(i,1) = sol(1); K(i,2) = sol(2);
        }
    }

    // Joseph update P <- P - K PCt^T - PCt K^T + K S K^T over the first M
    // measurement columns of K and PCt (M = 1, 2, 3).
    template<int M>
    void joseph_update_(const MatrixNX3& K, const Matrix3& S, const MatrixNX3& PCt) {
        for (int i = 0; i < NX; ++i) {
            for (int j = i; j < NX; ++j) {
                T kcp_ij = T(0), kcp_ji = T(0), ksk = T(0);
                for (int l = 0; l < M; ++l) {
                    kcp_ij += K(i,l) * PCt(j,l);
                    kcp_ji += K(j,l) * PCt(i,l);
                }
                for (int a = 0; a < M; ++a)
                    for (int b = 0; b < M; ++b) ksk += K(i,a) * S(a,b) * K(j,b);
                Pext(i,j) += -(kcp_ij + kcp_ji) + ksk;
                if (j != i) Pext(j,i) = Pext(i,j);
            }
        }
        symmetrize_Pext_();
    }

    // Linear update on M state indices: z = x[idx] + n, n ~ N(0, R).
    template<int M>
    void select_update_(const int (&idx)[M], const Vector3& r, const Matrix3& R, bool freeze_ba);

    void zero_AL_cross_cov_once_() {
        Pext.template block<BASE_N,NL>(0, OFF_L).setZero();
        Pext.template block<NL,BASE_N>(OFF_L, 0).setZero();
    }

    void freeze_acc_bias_rows_(MatrixNX3& M) const {
        if constexpr (with_accel_bias) M.template block<3,3>(OFF_BA, 0).setZero();
    }

    void refresh_vertical_discretization_(T Ts);

    // out = F_LL in.
    void apply_F_LL_(const T (&in)[NL], T (&out)[NL]) const {
        for (int ax = 0; ax < 2; ++ax)
            for (int i = 0; i < 4; ++i) {
                T acc = T(0);
                for (int j = 0; j < 4; ++j) acc += Phi_h_(i,j) * in[2 * j + ax];
                out[2 * i + ax] = acc;
            }
        for (int i = 0; i < NZ; ++i) {
            T acc = T(0);
            for (int j = 0; j < NZ; ++j) acc += F_Z_cache_(i,j) * in[NH + j];
            out[NH + i] = acc;
        }
    }
};

// ------------------------------------------------------------------------
// Implementation
// ------------------------------------------------------------------------

template <typename T, int NC, bool VI, bool GB, bool AB>
Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::Kalman3D_Wave_CARMA(
    Vector3 const& sigma_a, Vector3 const& gyro_noise_density_rad_sqrt_s, Vector3 const& sigma_m,
    T Pq0, T Pb0, T b0, T R_S_noise_var, T gravity_magnitude, const carma::CarmaShape& shape)
  : gravity_magnitude_(gravity_magnitude),
    Rmag(sigma_m.array().square().matrix().asDiagonal()),
    Racc(sigma_a.array().square().matrix().asDiagonal()),
    vmodel_(shape)
{
    qref.setIdentity();
    Qbase.setZero();
    if constexpr (GB) {
        Qbase.template topLeftCorner<3,3>() = gyro_noise_density_rad_sqrt_s.array().square().matrix().asDiagonal();
        Qbase.template bottomRightCorner<3,3>() = Matrix3::Identity() * b0;
    } else {
        Qbase = gyro_noise_density_rad_sqrt_s.array().square().matrix().asDiagonal();
    }
    R_S = Matrix3::Identity() * R_S_noise_var;

    xext.setZero();
    Pext.setZero();
    Pext.template topLeftCorner<3,3>() = Matrix3::Identity() * Pq0;
    if constexpr (GB) Pext.template block<3,3>(3,3) = Matrix3::Identity() * Pb0;

    // Horizontal chain: the OU-III seeds.
    for (int k = 0; k < 2; ++k) {
        Pext(OFF_VH + k, OFF_VH + k) = T(1.0);
        Pext(OFF_PH + k, OFF_PH + k) = T(20.0 * 20.0);
        Pext(OFF_SH + k, OFF_SH + k) = sigma_S0_ * sigma_S0_;
        Pext(OFF_AH + k, OFF_AH + k) = Sigma_aw_h(k);
    }
    if constexpr (AB) {
        Pext.template block<3,3>(OFF_BA, OFF_BA) = Matrix3::Identity() * sigma_bacc0_ * sigma_bacc0_;
    }
    // Vertical: stationary prior at the default operating point until the
    // wrapper configures it.
    vmodel_.configure(T(1), T(1));
    reset_vertical_covariance_to_stationary();
}

template <typename T, int NC, bool VI, bool GB, bool AB>
void Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::set_acc_bias_updates_enabled(bool en) {
    if (acc_bias_updates_enabled_ == en) return;
    if constexpr (AB) {
        if (!en) {
            Pext.template block<3,BASE_N>(OFF_BA, 0).setZero();
            Pext.template block<BASE_N,3>(0, OFF_BA).setZero();
            Pext.template block<3,NL>(OFF_BA, OFF_L).setZero();
            Pext.template block<NL,3>(OFF_L, OFF_BA).setZero();
        } else {
            const T target = sigma_bacc0_ * sigma_bacc0_;
            for (int i = 0; i < 3; ++i)
                Pext(OFF_BA + i, OFF_BA + i) = std::max(Pext(OFF_BA + i, OFF_BA + i), target);
        }
    }
    acc_bias_updates_enabled_ = en;
}

template <typename T, int NC, bool VI, bool GB, bool AB>
bool Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::initialize_from_acc_mag(Vector3 const& acc, Vector3 const& mag) {
    if (!acc.allFinite() || !mag.allFinite()) return false;
    const T anorm = acc.norm();
    if (!(anorm >= T(1e-8))) return false;
    if (!(mag.norm() >= T(1e-8))) return false;
    const Vector3 z_world = -acc / anorm;
    Vector3 mag_h = mag - (mag.dot(z_world)) * z_world;
    const T mh = mag_h.norm();
    if (!(mh >= T(1e-8))) return false;
    mag_h /= mh;
    project_gyro_bias_();
    const Vector3 y_world = z_world.cross(mag_h).normalized();
    Matrix3 R_wb_new;
    R_wb_new.col(0) = mag_h;
    R_wb_new.col(1) = y_world;
    R_wb_new.col(2) = z_world;
    qref = Eigen::Quaternion<T>(R_wb_new);
    qref.normalize();
    v2ref = R_bw() * mag;
    xext.template head<3>().setZero();
    Pext.template block<3,3>(0,0) = Matrix3::Identity() * T(5e-4);
    if constexpr (GB) {
        Pext.template block<3,3>(0,3).setZero();
        Pext.template block<3,3>(3,0).setZero();
        for (int i = 0; i < 3; ++i) Pext(3 + i, 3 + i) = std::max(Pext(3 + i, 3 + i), T(1e-6));
    }
    zero_AL_cross_cov_once_();
    if constexpr (AB) {
        Pext.template block<3,3>(0, OFF_BA).setZero();
        Pext.template block<3,3>(OFF_BA, 0).setZero();
    }
    symmetrize_Pext_();
    return true;
}

template <typename T, int NC, bool VI, bool GB, bool AB>
Eigen::Quaternion<T> Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::quaternion_from_acc(Vector3 const& acc) {
    const Vector3 an = acc.normalized();
    const Vector3 zb = Vector3::UnitZ();
    const Vector3 target = -an;
    T cos_theta = zb.dot(target);
    Vector3 axis = zb.cross(target);
    const T norm_axis = axis.norm();
    if (norm_axis < T(1e-8)) {
        return cos_theta > 0 ? Eigen::Quaternion<T>::Identity() : Eigen::Quaternion<T>(0, 1, 0, 0);
    }
    axis /= norm_axis;
    cos_theta = std::max(T(-1), std::min(T(1), cos_theta));
    Eigen::Quaternion<T> q(Eigen::AngleAxis<T>(std::acos(cos_theta), axis));
    q.normalize();
    return q;
}

template <typename T, int NC, bool VI, bool GB, bool AB>
bool Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::initialize_from_acc(Vector3 const& acc) {
    const T anorm = acc.norm();
    if (!acc.allFinite() || !(anorm >= T(1e-8))) return false;
    project_gyro_bias_();
    qref = quaternion_from_acc(acc / anorm);
    qref.normalize();
    xext.template head<3>().setZero();
    set_accel_only_attitude_covariance_();
    if constexpr (GB) {
        Pext.template block<3,3>(0,3).setZero();
        Pext.template block<3,3>(3,0).setZero();
    }
    zero_AL_cross_cov_once_();
    if constexpr (AB) {
        Pext.template block<3,3>(0, OFF_BA).setZero();
        Pext.template block<3,3>(OFF_BA, 0).setZero();
    }
    symmetrize_Pext_();
    return true;
}

template <typename T, int NC, bool VI, bool GB, bool AB>
bool Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::initialize_from_attitude(
    Eigen::Quaternion<T> const& q_bw, T tilt_sigma_rad, T yaw_sigma_rad)
{
    if (!q_bw.coeffs().allFinite() || !(q_bw.norm() > T(1e-8))) return false;
    set_quaternion_boat(q_bw);
    const T tilt = (std::isfinite(tilt_sigma_rad) && tilt_sigma_rad > T(0)) ? tilt_sigma_rad : T(0.035);
    const T yaw = (std::isfinite(yaw_sigma_rad) && yaw_sigma_rad > T(0)) ? yaw_sigma_rad : T(1.5708);
    set_accel_only_attitude_covariance_(tilt, yaw);
    if constexpr (GB) {
        Pext.template block<3,3>(0,3).setZero();
        Pext.template block<3,3>(3,0).setZero();
    }
    zero_AL_cross_cov_once_();
    if constexpr (AB) {
        Pext.template block<3,3>(0, OFF_BA).setZero();
        Pext.template block<3,3>(OFF_BA, 0).setZero();
    }
    symmetrize_Pext_();
    return true;
}

template <typename T, int NC, bool VI, bool GB, bool AB>
bool Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::initialize_from_acc_preserve_yaw(Vector3 const& acc_body) {
    const Eigen::Quaternion<T> q_old = quaternion_boat();
    if (!q_old.coeffs().allFinite()) return initialize_from_acc(acc_body);
    const T yaw_old = std::atan2(T(2) * (q_old.w() * q_old.z() + q_old.x() * q_old.y()),
                                 T(1) - T(2) * (q_old.y() * q_old.y() + q_old.z() * q_old.z()));
    if (!initialize_from_acc(acc_body)) return false;
    Eigen::Quaternion<T> qt = quaternion_boat();
    if (!(qt.norm() > T(1e-8))) return true;
    qt.normalize();
    const T sinp = std::max(T(-1), std::min(T(1), T(2) * (qt.w() * qt.y() - qt.z() * qt.x())));
    const T pitch = std::asin(sinp);
    const T roll = std::atan2(T(2) * (qt.w() * qt.x() + qt.y() * qt.z()),
                              T(1) - T(2) * (qt.x() * qt.x() + qt.y() * qt.y()));
    Eigen::Quaternion<T> q_new =
        Eigen::Quaternion<T>(Eigen::AngleAxis<T>(yaw_old, Vector3::UnitZ())) *
        Eigen::Quaternion<T>(Eigen::AngleAxis<T>(pitch, Vector3::UnitY())) *
        Eigen::Quaternion<T>(Eigen::AngleAxis<T>(roll, Vector3::UnitX()));
    q_new.normalize();
    set_quaternion_boat(q_new);
    return true;
}

template <typename T, int NC, bool VI, bool GB, bool AB>
void Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::refresh_vertical_discretization_(T Ts) {
    if (Ts == cache_dt_ && vmodel_.generation() == cache_generation_) return;
    if constexpr (VI) {
        MatrixZZ A = MatrixZZ::Zero();
        MatrixZZ W = MatrixZZ::Zero();
        A.template topLeftCorner<NC,NC>() = vmodel_.A();
        A(NC, 0) = T(1);    // S_z' = h
        W.template topLeftCorner<NC,NC>() = vmodel_.W();
        carma::discretize_lti<T,NZ>(A, W, Ts, F_Z_cache_, Q_Z_cache_);
    } else {
        carma::discretize_lti<T,NZ>(vmodel_.A(), vmodel_.W(), Ts, F_Z_cache_, Q_Z_cache_);
    }
    cache_dt_ = Ts;
    cache_generation_ = vmodel_.generation();
}

template <typename T, int NC, bool VI, bool GB, bool AB>
void Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::time_update(Vector3 const& gyr, T Ts) {
    project_gyro_bias_();
    if (!(Ts > T(0)) || !std::isfinite(Ts)) return;

    Vector3 gyro_bias = Vector3::Zero();
    if constexpr (GB) gyro_bias = xext.template segment<3>(3);
    last_gyr_bias_corrected = gyr - gyro_bias;

    const Eigen::Quaternion<T> dq_wb = ou_detail::quat_from_delta_theta((-last_gyr_bias_corrected * Ts).eval());
    qref = dq_wb * qref;
    qref.normalize();

    // Attitude(+gyro bias) block: the OU-III transition and structured Qd.
    const Vector3 w = last_gyr_bias_corrected;
    const T omega = w.norm();
    const T theta = omega * Ts;
    MatrixBaseN& F_AA = F_AA_scratch_;
    F_AA.setIdentity();
    const Matrix3 I3 = Matrix3::Identity();
    if (theta < T(1e-5)) {
        const Matrix3 Wx = skew_(w);
        F_AA.template topLeftCorner<3,3>() = I3 - Wx * Ts + (Wx * Wx) * (Ts * Ts / T(2));
    } else {
        const Matrix3 Wn = skew_(w / (omega + std::numeric_limits<T>::epsilon()));
        const T s = std::sin(theta), c = std::cos(theta);
        F_AA.template topLeftCorner<3,3>() = I3 - s * Wn + (T(1) - c) * (Wn * Wn);
    }
    if constexpr (GB) {
        Matrix3 Rstep, Bstep;
        ou_detail::rot_and_B_from_wt(w, Ts, Rstep, Bstep);
        F_AA.template topLeftCorner<3,3>() = Rstep;
        F_AA.template block<3,3>(0,3) = Bstep;
    }

    MatrixBaseN& Q_AA = Q_AA_scratch_;
    Q_AA.setZero();
    if (!use_exact_att_bias_Qd_) {
        Q_AA = Qbase * Ts;
    } else {
        const Matrix3 Qg = Qbase.template topLeftCorner<3,3>();
        Matrix3 Qbg = Matrix3::Zero();
        if constexpr (GB) Qbg = Qbase.template bottomRightCorner<3,3>();
        Matrix3 I_R;
        if (ou_detail::is_isotropic3(Qg)) I_R = Matrix3::Identity() * (Qg(0,0) * Ts);
        else I_R = ou_detail::simpson_R_Q_RT(w, Ts, Qg);
        Matrix3 I_BB = Matrix3::Zero();
        if constexpr (GB) I_BB = ou_detail::simpson_B_Q_BT(w, Ts, Qbg);
        Q_AA.template topLeftCorner<3,3>() = I_R + I_BB;
        if constexpr (GB) {
            Matrix3 IB;
            ou_detail::integral_B_ds(w, Ts, IB);
            const Matrix3 Qtb = IB * Qbg;
            Q_AA.template topRightCorner<3,3>() = Qtb;
            Q_AA.template bottomLeftCorner<3,3>() = Qtb.transpose();
            Q_AA.template bottomRightCorner<3,3>() = Qbg * Ts;
        }
        Q_AA = T(0.5) * (Q_AA + Q_AA.transpose()).eval();
        ou_detail::project_psd_ou_iii<T,BASE_N>(Q_AA, T(1e-12));
    }
    {
        MatrixBaseN& tmpAA = tmpAA_scratch_;
        tmpAA.noalias() = F_AA * Pext.template block<BASE_N,BASE_N>(0,0);
        Pext.template block<BASE_N,BASE_N>(0,0).noalias() = tmpAA * F_AA.transpose();
        Pext.template block<BASE_N,BASE_N>(0,0) += Q_AA;
    }

    // Linear block.  F_LL is block diagonal -- one OU-III chain per horizontal
    // axis on the interleaved indices {ax, 2+ax, 4+ax, 6+ax} and the vertical
    // block on [NH, NL) -- so it is applied through apply_F_LL_() instead of
    // being formed densely.  Besides the cost (2 NL (32 + NZ^2) instead of
    // 2 NL^3 multiply-adds), the dense 11 x 11 loop of the CARMA(3,0) build
    // faulted under GCC 13 -O3 on x86-64 (an aligned vector load emitted for
    // a misaligned column; clean under clang, ASan and UBSan); the
    // structured form avoids that pattern.
    {
        const T tau = std::max(T(1e-6), tau_aw);
        ou_detail::IntegratedOUChain<T,3>::transition(tau, Ts, Phi_h_);
        for (int ax = 0; ax < 2; ++ax)
            ou_detail::IntegratedOUChain<T,3>::process_covariance(tau, Ts, Sigma_aw_h(ax), Qd_h_[ax]);
    }
    refresh_vertical_discretization_(Ts);

    T in[NL], out[NL];

    // Mean: x_L <- F_LL x_L.
    for (int k = 0; k < NL; ++k) in[k] = xext(OFF_L + k);
    apply_F_LL_(in, out);
    for (int k = 0; k < NL; ++k) xext(OFF_L + k) = out[k];

    // P_LL <- F (F P_LL)^T + Q: columns of P_LL, then rows of the product.
    MatrixLL& FP = tmpLL_scratch_;
    for (int j = 0; j < NL; ++j) {
        for (int k = 0; k < NL; ++k) in[k] = Pext(OFF_L + k, OFF_L + j);
        apply_F_LL_(in, out);
        for (int k = 0; k < NL; ++k) FP(k, j) = out[k];
    }
    for (int i = 0; i < NL; ++i) {
        for (int k = 0; k < NL; ++k) in[k] = FP(i, k);
        apply_F_LL_(in, out);
        for (int j = 0; j < NL; ++j) Pext(OFF_L + i, OFF_L + j) = out[j];
    }
    for (int ax = 0; ax < 2; ++ax)
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                Pext(OFF_L + 2 * a + ax, OFF_L + 2 * b + ax) += Qd_h_[ax](a, b);
    for (int a = 0; a < NZ; ++a)
        for (int b = 0; b < NZ; ++b)
            Pext(OFF_Z + a, OFF_Z + b) += Q_Z_cache_(a, b);

    // P_AL <- F_AA P_AL F_LL^T
    Eigen::Matrix<T,BASE_N,NL>& tmpAL = tmpAL_scratch_;
    for (int i = 0; i < BASE_N; ++i)
        for (int j = 0; j < NL; ++j) {
            T acc = T(0);
            for (int k = 0; k < BASE_N; ++k) acc += F_AA(i,k) * Pext(k, OFF_L + j);
            tmpAL(i,j) = acc;
        }
    for (int i = 0; i < BASE_N; ++i) {
        for (int k = 0; k < NL; ++k) in[k] = tmpAL(i, k);
        apply_F_LL_(in, out);
        for (int j = 0; j < NL; ++j) {
            Pext(i, OFF_L + j) = out[j];
            Pext(OFF_L + j, i) = out[j];
        }
    }

    // Residual accelerometer-bias OU (as OU-III).
    if constexpr (AB) {
        Matrix3 P_BB = Pext.template block<3,3>(OFF_BA, OFF_BA);
        const T tau_b = std::max(T(1e-3), tau_bacc_);
        const T phi_b = acc_bias_updates_enabled_ ? std::exp(-Ts / tau_b) : T(1);
        if (acc_bias_updates_enabled_) {
            xext.template segment<3>(OFF_BA) *= phi_b;
            P_BB *= phi_b * phi_b;
            const T qd_scale = -T(0.5) * tau_b * std::expm1(-T(2) * Ts / tau_b);
            P_BB.noalias() += Q_bacc_ * qd_scale;
        }
        Pext.template block<3,3>(OFF_BA, OFF_BA) = P_BB;

        Eigen::Matrix<T,BASE_N,3> tmpAB;
        for (int i = 0; i < BASE_N; ++i)
            for (int j = 0; j < 3; ++j) {
                T acc = T(0);
                for (int k = 0; k < BASE_N; ++k) acc += F_AA(i,k) * Pext(k, OFF_BA + j);
                tmpAB(i,j) = acc * phi_b;
            }
        Pext.template block<BASE_N,3>(0, OFF_BA) = tmpAB;
        Pext.template block<3,BASE_N>(OFF_BA, 0) = tmpAB.transpose();

        // P_LB <- F_LL P_LB phi_b
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < NL; ++k) in[k] = Pext(OFF_L + k, OFF_BA + j);
            apply_F_LL_(in, out);
            for (int k = 0; k < NL; ++k) {
                Pext(OFF_L + k, OFF_BA + j) = out[k] * phi_b;
                Pext(OFF_BA + j, OFF_L + k) = out[k] * phi_b;
            }
        }
    }

    apply_pending_aw_covariance_inflation_();
    symmetrize_Pext_();

    if (ou_detail::periodic_update_due(Ts, pseudo_update_period_s_, pseudo_update_elapsed_s_)) {
        applyIntegralZeroPseudoMeas();
    }
}

template <typename T, int NC, bool VI, bool GB, bool AB>
void Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::measurement_update_acc_only(Vector3 const& acc_meas, T tempC) {
    last_acc_diag_ = MeasDiag3{};
    if (!acc_meas.allFinite() || !std::isfinite(tempC)) return;
    [[maybe_unused]] const bool use_ba = (AB && acc_bias_updates_enabled_);

    const Vector3 g_world(0, 0, gravity_magnitude_);
    Vector3 ba_term = Vector3::Zero();
    if constexpr (AB) ba_term = xext.template segment<3>(OFF_BA) + k_a_ * (tempC - tempC_ref);

    const Matrix3 Rwb = R_wb();
    const Vector3 aw(xext(OFF_AH), xext(OFF_AH + 1), xext(OFF_Z + 2));
    const Vector3 f_cog_b = Rwb * (aw - g_world);
    const Vector3 r = acc_meas - (f_cog_b + ba_term);
    last_acc_diag_.r = r;

    // Sparse C: attitude (J_att = -[f]x), a_x, a_y, a_z (columns of R_wb),
    // b_a (identity).
    const Matrix3 J_att = -skew_(f_cog_b);
    constexpr int NCOL = 6 + (AB ? 3 : 0);
    int col[NCOL];
    Vector3 cvec[NCOL];
    for (int k = 0; k < 3; ++k) { col[k] = k; cvec[k] = J_att.col(k); }
    col[3] = OFF_AH;     cvec[3] = Rwb.col(0);
    col[4] = OFF_AH + 1; cvec[4] = Rwb.col(1);
    col[5] = OFF_Z + 2;  cvec[5] = Rwb.col(2);
    if constexpr (AB) {
        for (int k = 0; k < 3; ++k) { col[6 + k] = OFF_BA + k; cvec[6 + k] = Vector3::Unit(k); }
    }

    // Full P C^T (bias column always: S marginalizes the bias uncertainty
    // whether or not it is being learned).
    MatrixNX3& PCt = PCt_scratch_;
    PCt.setZero();
    for (int i = 0; i < NX; ++i) {
        Vector3 acc = Vector3::Zero();
        for (int c = 0; c < NCOL; ++c) acc += Pext(i, col[c]) * cvec[c];
        PCt(i,0) = acc(0); PCt(i,1) = acc(1); PCt(i,2) = acc(2);
    }
    Matrix3& S_mat = S_scratch_;
    S_mat = Racc;
    for (int c = 0; c < NCOL; ++c) {
        const Vector3 pc(PCt(col[c],0), PCt(col[c],1), PCt(col[c],2));
        S_mat.noalias() += cvec[c] * pc.transpose();
    }
    S_mat = T(0.5) * (S_mat + S_mat.transpose()).eval();

    if constexpr (AB) {
        if (!use_ba) {
            // OU-III: the gain does not route through the frozen bias column.
            for (int i = 0; i < NX; ++i)
                for (int k = 0; k < 3; ++k) PCt(i,k) -= Pext(i, OFF_BA + k);
            freeze_acc_bias_rows_(PCt);
        }
    }

    Eigen::LDLT<Matrix3> ldlt;
    if (!safe_ldlt3_(S_mat, ldlt, Racc.norm())) {
        last_acc_diag_.S = S_mat;
        return;
    }
    last_acc_diag_.S = S_mat;
    last_acc_diag_.nis = nis3_(ldlt, r);
    MatrixNX3& K = K_scratch_;
    gain3_(PCt, ldlt, K);
    if constexpr (AB) {
        if (!use_ba) freeze_acc_bias_rows_(K);
    }
    xext.noalias() += K * r;
    joseph_update_<3>(K, S_mat, PCt);
    applyQuaternionCorrectionFromErrorState();
    last_acc_diag_.accepted = true;
}

template <typename T, int NC, bool VI, bool GB, bool AB>
void Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::measurement_update_mag_only(const Vector3& mag_meas) {
    last_mag_diag_ = MeasDiag3{};
    if (!mag_meas.allFinite()) return;
    if (!(mag_meas.norm() > T(1e-6))) return;
    const Vector3 v2hat = R_wb() * v2ref;
    const Vector3 r = mag_meas - v2hat;
    last_mag_diag_.r = r;
    const Matrix3 J_att = -skew_(v2hat);
    Matrix3& S_mat = S_scratch_;
    S_mat = Rmag;
    S_mat.noalias() += J_att * Pext.template block<3,3>(0,0) * J_att.transpose();
    MatrixNX3& PCt = PCt_scratch_;
    PCt.noalias() = Pext.template block<NX,3>(0,0) * J_att.transpose();
    if constexpr (AB) {
        if (!acc_bias_updates_enabled_) freeze_acc_bias_rows_(PCt);
    }
    Eigen::LDLT<Matrix3> ldlt;
    if (!safe_ldlt3_(S_mat, ldlt, Rmag.norm())) {
        last_mag_diag_.S = S_mat;
        return;
    }
    last_mag_diag_.S = S_mat;
    last_mag_diag_.nis = nis3_(ldlt, r);
    MatrixNX3& K = K_scratch_;
    gain3_(PCt, ldlt, K);
    if constexpr (AB) {
        if (!acc_bias_updates_enabled_) freeze_acc_bias_rows_(K);
    }
    xext.noalias() += K * r;
    joseph_update_<3>(K, S_mat, PCt);
    applyQuaternionCorrectionFromErrorState();
    last_mag_diag_.accepted = true;
}

template <typename T, int NC, bool VI, bool GB, bool AB>
template <int M>
void Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::select_update_(const int (&idx)[M], const Vector3& r,
                                                        const Matrix3& R, bool freeze_ba) {
    for (int k = 0; k < M; ++k) if (!std::isfinite(r(k))) return;
    Matrix3& S_mat = S_scratch_;
    S_mat.setIdentity();
    for (int a = 0; a < M; ++a)
        for (int b = 0; b < M; ++b) S_mat(a,b) = Pext(idx[a], idx[b]) + R(a,b);
    MatrixNX3& PCt = PCt_scratch_;
    PCt.setZero();
    for (int i = 0; i < NX; ++i)
        for (int k = 0; k < M; ++k) PCt(i,k) = Pext(i, idx[k]);
    if (freeze_ba) freeze_acc_bias_rows_(PCt);
    Eigen::LDLT<Matrix3> ldlt;
    if (!safe_ldlt3_(S_mat, ldlt, R.norm())) return;
    MatrixNX3& K = K_scratch_;
    gain3_(PCt, ldlt, K);
    if (freeze_ba) freeze_acc_bias_rows_(K);
    for (int i = 0; i < NX; ++i) {
        T s = T(0);
        for (int k = 0; k < M; ++k) s += K(i,k) * r(k);
        xext(i) += s;
    }
    joseph_update_<M>(K, S_mat, PCt);
    applyQuaternionCorrectionFromErrorState();
}

template <typename T, int NC, bool VI, bool GB, bool AB>
void Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::applyIntegralZeroPseudoMeas() {
    const bool freeze = AB && !acc_bias_updates_enabled_;
    if constexpr (VI) {
        const int idx[3] = {OFF_SH, OFF_SH + 1, OFF_SZ};
        const Vector3 r(-xext(OFF_SH), -xext(OFF_SH + 1), -xext(OFF_SZ));
        select_update_<3>(idx, r, R_S, freeze);
    } else {
        // M = 2 uses the leading 2x2 of the 3x3 scratch; the third diagonal is
        // the identity filler and never enters K (PCt column 2 is zero).
        const int idx[2] = {OFF_SH, OFF_SH + 1};
        const Vector3 r(-xext(OFF_SH), -xext(OFF_SH + 1), T(0));
        select_update_<2>(idx, r, R_S, freeze);
    }
}

template <typename T, int NC, bool VI, bool GB, bool AB>
void Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::measurement_update_position_pseudo(
    const Vector3& p_meas, const Vector3& sigma_meas) {
    const int idx[3] = {OFF_PH, OFF_PH + 1, OFF_Z};
    const Vector3 r = p_meas - get_position();
    const Vector3 s = sigma_meas.cwiseMax(T(0));
    select_update_<3>(idx, r, Matrix3(s.array().square().matrix().asDiagonal()), false);
}

template <typename T, int NC, bool VI, bool GB, bool AB>
void Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::measurement_update_velocity_pseudo(
    const Vector3& v_meas, const Vector3& sigma_meas) {
    const int idx[3] = {OFF_VH, OFF_VH + 1, OFF_Z + 1};
    const Vector3 r = v_meas - get_velocity();
    const Vector3 s = sigma_meas.cwiseMax(T(0));
    select_update_<3>(idx, r, Matrix3(s.array().square().matrix().asDiagonal()), false);
}

template <typename T, int NC, bool VI, bool GB, bool AB>
void Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::measurement_update_vert_velocity_pseudo(T vz_meas, T sigma_meas) {
    if (!std::isfinite(vz_meas)) return;
    const int idx[1] = {OFF_Z + 1};
    const T s = std::max(T(0), sigma_meas);
    Matrix3 R = Matrix3::Identity();
    R(0,0) = s * s;
    select_update_<1>(idx, Vector3(vz_meas - xext(OFF_Z + 1), T(0), T(0)), R, false);
}

template <typename T, int NC, bool VI, bool GB, bool AB>
void Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::applyQuaternionCorrectionFromErrorState() {
    project_gyro_bias_();
    const Vector3 dtheta = xext.template segment<3>(0);
    if (!dtheta.allFinite()) {
        xext.template head<3>().setZero();
        return;
    }
    qref = ou_detail::quat_from_delta_theta(dtheta) * qref;
    qref.normalize();
    ou_detail::apply_left_error_reset<T, NX>(Pext, dtheta);
    xext.template head<3>().setZero();
    project_acc_bias_();
}

template <typename T, int NC, bool VI, bool GB, bool AB>
void Kalman3D_Wave_CARMA<T,NC,VI,GB,AB>::project_acc_bias_() {
    if constexpr (AB) {
        if (!(acc_bias_limit_ > T(0))) return;
        auto b = xext.template segment<3>(OFF_BA);
        if (!b.allFinite()) { b.setZero(); return; }
        const T n = b.norm();
        if (n > acc_bias_limit_) b *= (acc_bias_limit_ / n);
    }
}

}  // namespace ocean_imu::kalman
