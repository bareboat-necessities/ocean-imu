#pragma once

/*
  Copyright (c) 2025-2026 Mikhail Grushinskiy

  Based on: https://github.com/thomaspasser/q-mekf
  MIT License, Copyright (c) 2023 Thomas Passer

  Kalman3D_Wave_CARMA3 -- EXPERIMENTAL three-axis research filter.  Not a
  shipping estimator: Kalman3D_Wave_OU_III remains the deployed, authoritative
  implementation, and nothing in the deployed configuration includes this
  header; see docs/carma-three-axis.md.

  Every translational axis has a CARMA prior in the library coordinates
  [p, v, a, latent...] of CarmaVerticalModel.h:

     horizontal    CARMA(NHC, NHC-3) per axis, NHC = 3, 4, 5, optionally
                   augmented with S = int p dt (SH) and its zero
                   pseudo-measurement;
     vertical      CARMA(NZC, NZC-3) [h, v, a, latent...] (+ S_z with SZ).

  State (NED, z down): [dtheta(3), b_g(3), x block, y block, z block, b_a(3)].

  Horizontal CARMA prior.  Both horizontal axes share one transition matrix
  A_h (one shape, one frequency scale: the same waves drive both axes) and
  have the joint driving-noise intensity Sigma_h (x) W_1, W_1 the unit
  acceleration-variance intensity of A_h and Sigma_h = [[s_x^2, r s_x s_y],
  [r s_x s_y, s_y^2]] the 2 x 2 directional acceleration covariance.  Because
  A is shared, the exact discretization is F = I_2 (x) F_1 and
  Q = Sigma_h (x) Q_1 (Kronecker), so one exact discretization of A_h serves
  both axes and a correlated (r != 0) prior costs nothing extra.  r = 0 is
  independent axes; the coupled variant uses r != 0.

  Retained from OU-III line for line in its mathematics: attitude and
  gyro-bias prediction, accelerometer and magnetometer models and Jacobians,
  error-state reset, residual accelerometer-bias OU and its freeze, Joseph
  updates, input guards.  The optional IMU lever arm
  (set_imu_lever_arm_body, off by default) is OU-III's measurement model with
  the band-limited rate derivative of kalman_common/ImuLeverArm.h.
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

#include "kalman_common/CovarianceProjection.h"
#include "kalman_common/ImuLeverArm.h"
#include "kalman_common/MekfAttitudeMath.h"
#include "kalman_common/PeriodicSchedule.h"
#include "kalman_carma/CarmaVerticalModel.h"

namespace ocean_imu::kalman {

template <typename T = float, int NHC = 4, bool SH = false, int NZC = 4, bool SZ = false,
          bool with_gyro_bias = true, bool with_accel_bias = true>
class Kalman3D_Wave_CARMA3 {
    static_assert(NHC >= 3 && NHC <= 5, "horizontal CARMA order 3..5");
    static_assert(NZC >= 3 && NZC <= 5, "vertical CARMA order 3..5");
public:
    static constexpr int BASE_N = with_gyro_bias ? 6 : 3;
    static constexpr int NHA = NHC + (SH ? 1 : 0);   // states per horizontal axis
    static constexpr int NH = 2 * NHA;
    static constexpr int NZ = NZC + (SZ ? 1 : 0);
    static constexpr int NL = NH + NZ;
    static constexpr int NX = BASE_N + NL + (with_accel_bias ? 3 : 0);

    static constexpr int OFF_L = BASE_N;
    static constexpr int OFF_Z = BASE_N + NH;
    static constexpr int OFF_SZ = SZ ? (OFF_Z + NZC) : -1;
    static constexpr int OFF_BA = with_accel_bias ? (BASE_N + NL) : -1;

    // Per-axis indices of the horizontal block.
    static constexpr int idx_p(int ax) { return BASE_N + ax * NHA + 0; }
    static constexpr int idx_v(int ax) { return BASE_N + ax * NHA + 1; }
    static constexpr int idx_a(int ax) { return BASE_N + ax * NHA + 2; }
    static constexpr int idx_S(int ax) { return SH ? BASE_N + ax * NHA + NHC : -1; }
    static constexpr int pos_index(int ax) { return ax < 2 ? idx_p(ax) : OFF_Z; }
    static constexpr int vel_index(int ax) { return ax < 2 ? idx_v(ax) : OFF_Z + 1; }

    using Vector3 = Eigen::Matrix<T,3,1>;
    using Matrix3 = Eigen::Matrix<T,3,3>;
    using MatrixBaseN = Eigen::Matrix<T,BASE_N,BASE_N>;
    using MatrixNX = Eigen::Matrix<T,NX,NX>;
    using MatrixNX3 = Eigen::Matrix<T,NX,3>;
    using MatrixLL = Eigen::Matrix<T,NL,NL>;
    using MatrixZZ = Eigen::Matrix<T,NZ,NZ>;
    using MatrixCZ = Eigen::Matrix<T,NZC,NZC>;
    using MatrixHA = Eigen::Matrix<T,NHA,NHA>;
    using MatrixCH = Eigen::Matrix<T,NHC,NHC>;
    using VerticalModel = carma::CarmaVerticalModel<T,NZC>;
    using HorizontalModel = carma::CarmaVerticalModel<T,NHC>;

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

    Kalman3D_Wave_CARMA3(Vector3 const& sigma_a, Vector3 const& gyro_noise_density_rad_sqrt_s,
                         Vector3 const& sigma_m,
                         T Pq0 = T(5e-4), T Pb0 = T(1e-6), T b0 = T(1e-11), T R_S_noise_var = T(1.5),
                         T gravity_magnitude = T(STD_GRAVITY),
                         const carma::CarmaShape& vshape = carma::CarmaDefaultShape<NZC>::value,
                         const carma::CarmaShape& hshape = carma::CarmaDefaultShape<NHC>::value);

    bool initialize_from_acc_mag(Vector3 const& acc, Vector3 const& mag);
    bool initialize_from_acc(Vector3 const& acc);
    bool initialize_from_acc_preserve_yaw(Vector3 const& acc);
    bool initialize_from_attitude(Eigen::Quaternion<T> const& q_bw, T tilt_sigma_rad, T yaw_sigma_rad);
    static Eigen::Quaternion<T> quaternion_from_acc(Vector3 const& acc);

    void set_mag_world_ref(const Vector3& B_world) { v2ref = B_world; }

    void time_update(Vector3 const& gyr, T Ts);
    void measurement_update_acc_only(Vector3 const& acc, T tempC = tempC_ref);
    void measurement_update_mag_only(Vector3 const& mag);
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
        return Vector3(xext(idx_v(0)), xext(idx_v(1)), xext(OFF_Z + 1));
    }
    [[nodiscard]] Vector3 get_position() const {
        return Vector3(xext(idx_p(0)), xext(idx_p(1)), xext(OFF_Z));
    }
    [[nodiscard]] Vector3 get_world_accel() const {
        return Vector3(xext(idx_a(0)), xext(idx_a(1)), xext(OFF_Z + 2));
    }
    [[nodiscard]] Vector3 get_integral_displacement() const {
        const T nan = std::numeric_limits<T>::quiet_NaN();
        T sx = nan, sy = nan, sz = nan;
        if constexpr (SH) { sx = xext(idx_S(0)); sy = xext(idx_S(1)); }
        if constexpr (SZ) sz = xext(OFF_SZ);
        return Vector3(sx, sy, sz);
    }
    // [p, v] covariance of one axis (0 = x, 1 = y, 2 = z, NED).
    [[nodiscard]] Eigen::Matrix<T,2,2> pv_covariance(int axis) const {
        const int ip = pos_index(axis), iv = vel_index(axis);
        Eigen::Matrix<T,2,2> M;
        M << Pext(ip, ip), Pext(ip, iv), Pext(iv, ip), Pext(iv, iv);
        return M;
    }
    [[nodiscard]] Eigen::Matrix<T,2,2> vertical_pv_covariance() const { return pv_covariance(2); }
    [[nodiscard]] Matrix3 position_covariance() const {
        Matrix3 M;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) M(i,j) = Pext(pos_index(i), pos_index(j));
        return M;
    }
    [[nodiscard]] Eigen::Matrix<T,NZ,1> vertical_state() const { return xext.template segment<NZ>(OFF_Z); }

    // Committed operating-point time constant (reported by the wrapper; the
    // CARMA priors are scaled through kappa, not through tau).
    [[nodiscard]] T get_aw_time_constant() const { return tau_aw; }
    [[nodiscard]] Vector3 get_aw_stationary_std() const {
        return Vector3(std::sqrt(std::max(T(0), Sigma_h_(0,0))), std::sqrt(std::max(T(0), Sigma_h_(1,1))),
                       vmodel_.accel_std());
    }
    [[nodiscard]] Vector3 get_RS_noise_std() const { return R_S.diagonal().cwiseMax(T(0)).cwiseSqrt(); }

    void set_aw_time_constant(T tau_seconds) { tau_aw = std::max(T(1e-3), tau_seconds); }
    void set_pseudo_update_period_s(T period_s) {
        if (!(period_s > T(0)) || !std::isfinite(period_s)) return;
        const T new_period = std::max(T(1e-4), period_s);
        pseudo_update_elapsed_s_ =
            filter_detail::retarget_period_elapsed_progress_preserving(pseudo_update_elapsed_s_, new_period);
        pseudo_update_period_s_ = new_period;
    }
    [[nodiscard]] T get_pseudo_update_period_s() const { return pseudo_update_period_s_; }

    // Stationary prior covariance of every CARMA block (the name is the
    // shared wrapper interface).
    void reset_aw_covariance_to_stationary() {
        reset_horizontal_covariance_to_stationary();
        reset_vertical_covariance_to_stationary();
        symmetrize_Pext_();
    }

    bool set_vertical_carma(T kappa, T accel_std) { return vmodel_.configure(kappa, accel_std); }
    [[nodiscard]] const VerticalModel& vertical_model() const noexcept { return vmodel_; }

    // Horizontal CARMA prior: frequency scale kappa (both axes) and the
    // directional acceleration covariance [[sx^2, r sx sy], [r sx sy, sy^2]].
    bool set_horizontal_carma(T kappa, T sx, T sy, T r = T(0)) {
        if (!(sx > T(0)) || !(sy > T(0)) || !std::isfinite(sx) || !std::isfinite(sy)) return false;
        if (!(std::abs(r) <= T(1))) return false;
        if (!hmodel_.configure(kappa, T(1))) return false;   // unit-variance A_h, W_1
        Sigma_h_ << sx * sx, r * sx * sy, r * sx * sy, sy * sy;
        return true;
    }
    bool set_horizontal_shape(const carma::CarmaShape& shape) { return hmodel_.set_shape(shape); }
    bool set_vertical_shape(const carma::CarmaShape& shape) { return vmodel_.set_shape(shape); }
    [[nodiscard]] const HorizontalModel& horizontal_model() const noexcept { return hmodel_; }
    [[nodiscard]] Eigen::Matrix<T,2,2> horizontal_sigma() const { return Sigma_h_; }

    void reset_vertical_covariance_to_stationary() {
        MatrixCZ P0;
        if (!vmodel_.stationary_covariance(P0)) return;
        for (int i = 0; i < NX; ++i) {
            if (i >= OFF_Z && i < OFF_Z + NZ) continue;
            for (int j = 0; j < NZ; ++j) {
                Pext(i, OFF_Z + j) = T(0);
                Pext(OFF_Z + j, i) = T(0);
            }
        }
        Pext.template block<NZ,NZ>(OFF_Z, OFF_Z).setZero();
        Pext.template block<NZC,NZC>(OFF_Z, OFF_Z) = P0;
        if constexpr (SZ) Pext(OFF_SZ, OFF_SZ) = sigma_S0_ * sigma_S0_;
        symmetrize_Pext_();
    }

    void reset_horizontal_covariance_to_stationary() {
        MatrixCH P1;
        if (!hmodel_.stationary_covariance(P1)) return;
        for (int i = 0; i < NX; ++i) {
            if (i >= OFF_L && i < OFF_L + NH) continue;
            for (int j = 0; j < NH; ++j) {
                Pext(i, OFF_L + j) = T(0);
                Pext(OFF_L + j, i) = T(0);
            }
        }
        Pext.template block<NH,NH>(OFF_L, OFF_L).setZero();
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b)
                Pext.template block<NHC,NHC>(OFF_L + a * NHA, OFF_L + b * NHA) = Sigma_h_(a,b) * P1;
        if constexpr (SH)
            for (int a = 0; a < 2; ++a) Pext(idx_S(a), idx_S(a)) = sigma_S0_ * sigma_S0_;
        symmetrize_Pext_();
    }

    // R_S per axis as a standard deviation [m s].
    void set_RS_noise(const Vector3& sigma_S) { R_S = sigma_S.array().square().matrix().asDiagonal(); }
    void set_Racc_std(const Vector3& sigma_acc) { Racc = sigma_acc.array().square().matrix().asDiagonal(); }
    void set_Rmag(const Vector3& sigma_mag) { Rmag = sigma_mag.array().square().matrix().asDiagonal(); }

    void set_initial_acc_bias_std(T s) {
        if constexpr (with_accel_bias) {
            sigma_bacc0_ = std::max(T(0), s);
            Pext.template block<3,3>(OFF_BA, OFF_BA) = Matrix3::Identity() * sigma_bacc0_ * sigma_bacc0_;
        }
    }
    void set_Q_bacc_rw(const Vector3& driving_std_per_sqrt_s) {
        if constexpr (with_accel_bias) Q_bacc_ = driving_std_per_sqrt_s.array().square().matrix().asDiagonal();
    }
    void set_acc_bias_time_constant(T tau_seconds) {
        if constexpr (with_accel_bias) tau_bacc_ = std::max(T(1e-3), tau_seconds);
    }
    [[nodiscard]] T get_acc_bias_time_constant() const noexcept { return tau_bacc_; }
    [[nodiscard]] Vector3 get_acc_bias_at_temperature(T tempC) const {
        if constexpr (with_accel_bias) return xext.template segment<3>(OFF_BA) + k_a_ * (tempC - tempC_ref);
        return Vector3::Zero();
    }
    [[nodiscard]] Vector3 get_acc_bias_body_at_temperature(T tempC) const { return get_acc_bias_at_temperature(tempC); }
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

    // IMU installation lever arm (off-CoG correction), kalman_common/
    // ImuLeverArm.h.  r_b is the IMU position relative to the CoG in the body
    // frame (x forward, y starboard, z down), metres.  Zero, the default,
    // switches the correction off and leaves the filter exactly as it is
    // without it.
    void set_imu_lever_arm_body(const Vector3& r_b) { lever_.set_offset_body(r_b); }
    void clear_imu_lever_arm() { lever_.clear(); }
    [[nodiscard]] const Vector3& imu_lever_arm_body() const noexcept { return lever_.offset_body(); }
    // Corner of the low-pass ahead of the rate derivative (default 15 Hz).
    void set_imu_lever_arm_cutoff_hz(T hz) { lever_.set_derivative_cutoff_hz(hz); }
    [[nodiscard]] T imu_lever_arm_cutoff_hz() const noexcept { return lever_.derivative_cutoff_hz(); }
    // Ablation only: OU-III's model in place of the default, and its alpha
    // smoothing (used by that model alone).  docs/carma-imu-lever-arm.md.
    using LeverArmModel = typename filter_detail::ImuLeverArm<T>::Model;
    void set_imu_lever_arm_model(LeverArmModel m) { lever_.set_model(m); }
    void set_alpha_smoothing_tau(T tau_sec) { lever_.set_alpha_smoothing_tau(tau_sec); }

private:
    const T gravity_magnitude_ = T(STD_GRAVITY);

    Eigen::Quaternion<T> qref;
    Vector3 v2ref = Vector3::UnitX();

    Eigen::Matrix<T,NX,1> xext;
    MatrixNX Pext;

    Vector3 last_gyr_bias_corrected{};
    filter_detail::ImuLeverArm<T> lever_;

    T sigma_bacc0_ = T(0.004);
    Matrix3 Q_bacc_ = Matrix3::Identity() * T(2.5e-7);
    T tau_bacc_ = T(5000.0);
    Vector3 k_a_ = Vector3::Constant(T(0.002));
    T acc_bias_limit_ = T(0.4);

    Matrix3 Rmag;
    MatrixBaseN Qbase;
    Matrix3 Racc;
    Matrix3 R_S;

    T tau_aw = T(2.1);

    // Horizontal CARMA prior (unit acceleration variance) and its
    // discretization at the current step; Sigma_h_ scales Q and the
    // stationary covariance.
    HorizontalModel hmodel_;
    Eigen::Matrix<T,2,2> Sigma_h_ = Eigen::Matrix<T,2,2>::Identity();
    carma::StepDiscretization<T,NHA> disc_h_;

    // Vertical prior and its discretization at the current step.
    VerticalModel vmodel_;
    carma::StepDiscretization<T,NZ> disc_z_;

    T sigma_S0_ = T(50.0);
    T pseudo_update_elapsed_s_ = T(0);
    T pseudo_update_period_s_ = T(0.015);

    bool acc_bias_updates_enabled_ = true;
    bool use_exact_att_bias_Qd_ = true;

    MatrixBaseN F_AA_scratch_;
    MatrixBaseN Q_AA_scratch_;
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
            filter_detail::project_gyro_bias<T>(bias);
        }
    }

    static Matrix3 skew_(const Vector3& v) {
        Matrix3 M;
        M << 0, -v(2), v(1),
             v(2), 0, -v(0),
            -v(1), v(0), 0;
        return M;
    }

    // filter_detail::simpson_B_Q_BT evaluated with the step's B(Ts) = B1, which
    // time_update already has, and B(0) = 0: (Ts/6) (4 Bm Q Bm^T + B1 Q B1^T).
    static Matrix3 simpson_B_Q_BT_(const Vector3& w, T Ts, const Matrix3& Q, const Matrix3& B1) {
        Matrix3 Rm, Bm;
        filter_detail::rot_and_B_from_wt(w, T(0.5) * Ts, Rm, Bm);
        T BQm[3][3], BQ1[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                T am = T(0), a1 = T(0);
                for (int k = 0; k < 3; ++k) {
                    am += Bm(i,k) * Q(k,j);
                    a1 += B1(i,k) * Q(k,j);
                }
                BQm[i][j] = am;
                BQ1[i][j] = a1;
            }
        const T c = Ts / T(6);
        Matrix3 out;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                T am = T(0), a1 = T(0);
                for (int k = 0; k < 3; ++k) {
                    am += BQm[i][k] * Bm(j,k);
                    a1 += BQ1[i][k] * B1(j,k);
                }
                out(i,j) = c * (T(4) * am + a1);
            }
        return out;
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

    // K = PCt S^-1 through the three columns of S^-1 (three vector solves, not
    // one per state).
    void gain3_(const MatrixNX3& PCt, const Eigen::LDLT<Matrix3>& ldlt, MatrixNX3& K) const {
        Matrix3 Sinv;
        for (int k = 0; k < 3; ++k) {
            Vector3 e = Vector3::Zero();
            e(k) = T(1);
            const Vector3 c = ldlt.solve(e);
            Sinv(0,k) = c(0); Sinv(1,k) = c(1); Sinv(2,k) = c(2);
        }
        for (int i = 0; i < NX; ++i)
            for (int b = 0; b < 3; ++b)
                K(i,b) = PCt(i,0) * Sinv(0,b) + PCt(i,1) * Sinv(1,b) + PCt(i,2) * Sinv(2,b);
    }

    // Joseph form P - K C P - (K C P)^T + K S K^T, with C P = PCt^T for the
    // symmetric P.  With G = K S - PCt the update is K G^T - PCt K^T, the
    // same expression for any gain (frozen rows included) at 2M instead of
    // M (M + 2) products per entry.  Both triangles get the same value, so
    // Pext stays exactly symmetric.  The gain is copied to locals, which the
    // stores into Pext cannot alias.
    template<int M>
    void joseph_update_(const MatrixNX3& K, const Matrix3& S, const MatrixNX3& PCt) {
        T Kl[NX][M], G[NX][M];
        for (int i = 0; i < NX; ++i)
            for (int b = 0; b < M; ++b) {
                Kl[i][b] = K(i,b);
                T g = -PCt(i,b);
                for (int a = 0; a < M; ++a) g += K(i,a) * S(a,b);
                G[i][b] = g;
            }
        for (int i = 0; i < NX; ++i) {
            T pi[M];
            for (int l = 0; l < M; ++l) pi[l] = PCt(i,l);
            for (int j = i; j < NX; ++j) {
                T d = T(0);
                for (int l = 0; l < M; ++l) d += Kl[i][l] * G[j][l] - pi[l] * Kl[j][l];
                const T v = Pext(i,j) + d;
                Pext(i,j) = v;
                Pext(j,i) = v;
            }
        }
    }

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
    void refresh_horizontal_discretization_(T Ts);

    // With gyro bias F_AA = [[R, B], [0, I]]: rows 3..5 are unit rows.
    static constexpr bool aa_unit_row_(int i) { return with_gyro_bias && i >= 3; }

    // Pext block (rows ra.., cols rb..) <- Fa P_ab Fb^T, mirrored into (b, a).
    // A diagonal block (ra == rb, Fa == Fb) computes its upper triangle.
    // Blocks below the diagonal are never read, so the upper block triangle
    // can be propagated in place.
    template<int NA, int NB>
    void propagate_block_(int ra, int rb, const Eigen::Matrix<T,NA,NA>& Fa,
                          const Eigen::Matrix<T,NB,NB>& Fb, bool diag) {
        T X[NA][NB];
        for (int i = 0; i < NA; ++i)
            for (int j = 0; j < NB; ++j) {
                T acc = T(0);
                for (int k = 0; k < NA; ++k) acc += Fa(i,k) * Pext(ra + k, rb + j);
                X[i][j] = acc;
            }
        for (int i = 0; i < NA; ++i)
            for (int j = diag ? i : 0; j < NB; ++j) {
                T acc = T(0);
                for (int k = 0; k < NB; ++k) acc += X[i][k] * Fb(j,k);
                Pext(ra + i, rb + j) = acc;
                Pext(rb + j, ra + i) = acc;
            }
    }

    // out = F_LL in.
    void apply_F_LL_(const T (&in)[NL], T (&out)[NL]) const {
        const MatrixHA& F_H = disc_h_.Phi();
        const MatrixZZ& F_Z = disc_z_.Phi();
        for (int ax = 0; ax < 2; ++ax)
            for (int i = 0; i < NHA; ++i) {
                T acc = T(0);
                for (int j = 0; j < NHA; ++j) acc += F_H(i,j) * in[ax * NHA + j];
                out[ax * NHA + i] = acc;
            }
        for (int i = 0; i < NZ; ++i) {
            T acc = T(0);
            for (int j = 0; j < NZ; ++j) acc += F_Z(i,j) * in[NH + j];
            out[NH + i] = acc;
        }
    }
};

// ------------------------------------------------------------------------
// Implementation
// ------------------------------------------------------------------------

#define CARMA3_TPL template <typename T, int NHC, bool SH, int NZC, bool SZ, bool GB, bool AB>
#define CARMA3_CLS Kalman3D_Wave_CARMA3<T,NHC,SH,NZC,SZ,GB,AB>

CARMA3_TPL
CARMA3_CLS::Kalman3D_Wave_CARMA3(
    Vector3 const& sigma_a, Vector3 const& gyro_noise_density_rad_sqrt_s, Vector3 const& sigma_m,
    T Pq0, T Pb0, T b0, T R_S_noise_var, T gravity_magnitude,
    const carma::CarmaShape& vshape, const carma::CarmaShape& hshape)
  : gravity_magnitude_(gravity_magnitude),
    Rmag(sigma_m.array().square().matrix().asDiagonal()),
    Racc(sigma_a.array().square().matrix().asDiagonal()),
    hmodel_(hshape),
    vmodel_(vshape)
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

    set_horizontal_carma(T(1), T(1), T(1), T(0));
    reset_horizontal_covariance_to_stationary();
    if constexpr (AB) {
        Pext.template block<3,3>(OFF_BA, OFF_BA) = Matrix3::Identity() * sigma_bacc0_ * sigma_bacc0_;
    }
    vmodel_.configure(T(1), T(1));
    reset_vertical_covariance_to_stationary();
}

CARMA3_TPL
void CARMA3_CLS::set_acc_bias_updates_enabled(bool en) {
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

CARMA3_TPL
bool CARMA3_CLS::initialize_from_acc_mag(Vector3 const& acc, Vector3 const& mag) {
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

CARMA3_TPL
Eigen::Quaternion<T> CARMA3_CLS::quaternion_from_acc(Vector3 const& acc) {
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

CARMA3_TPL
bool CARMA3_CLS::initialize_from_acc(Vector3 const& acc) {
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

CARMA3_TPL
bool CARMA3_CLS::initialize_from_attitude(Eigen::Quaternion<T> const& q_bw, T tilt_sigma_rad, T yaw_sigma_rad) {
    if (!q_bw.coeffs().allFinite() || !(q_bw.norm() > T(1e-8))) return false;
    set_quaternion_boat(q_bw);
    // The wrapper hands over here after holding the prediction; the lever
    // arm's rate history must not difference across that gap.
    lever_.reset_kinematics();
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

CARMA3_TPL
bool CARMA3_CLS::initialize_from_acc_preserve_yaw(Vector3 const& acc_body) {
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

CARMA3_TPL
void CARMA3_CLS::refresh_vertical_discretization_(T Ts) {
    if (disc_z_.advance(vmodel_.generation(), Ts)) return;
    if constexpr (SZ) {
        MatrixZZ A = MatrixZZ::Zero();
        MatrixZZ W = MatrixZZ::Zero();
        A.template topLeftCorner<NZC,NZC>() = vmodel_.A();
        A(NZC, 0) = T(1);    // S_z' = h
        W.template topLeftCorner<NZC,NZC>() = vmodel_.W();
        disc_z_.refresh(A, W, vmodel_.generation(), Ts);
    } else {
        disc_z_.refresh(vmodel_.A(), vmodel_.W(), vmodel_.generation(), Ts);
    }
}

CARMA3_TPL
void CARMA3_CLS::refresh_horizontal_discretization_(T Ts) {
    if (disc_h_.advance(hmodel_.generation(), Ts)) return;
    if constexpr (SH) {
        MatrixHA A = MatrixHA::Zero();
        MatrixHA W = MatrixHA::Zero();
        A.template topLeftCorner<NHC,NHC>() = hmodel_.A();
        A(NHC, 0) = T(1);   // S' = p
        W.template topLeftCorner<NHC,NHC>() = hmodel_.W();
        disc_h_.refresh(A, W, hmodel_.generation(), Ts);
    } else {
        disc_h_.refresh(hmodel_.A(), hmodel_.W(), hmodel_.generation(), Ts);
    }
}

CARMA3_TPL
void CARMA3_CLS::time_update(Vector3 const& gyr, T Ts) {
    project_gyro_bias_();
    if (!(Ts > T(0)) || !std::isfinite(Ts)) return;

    Vector3 gyro_bias = Vector3::Zero();
    if constexpr (GB) gyro_bias = xext.template segment<3>(3);
    last_gyr_bias_corrected = gyr - gyro_bias;
    lever_.update_rate(last_gyr_bias_corrected, Ts);

    const Eigen::Quaternion<T> dq_wb = filter_detail::quat_from_delta_theta((-last_gyr_bias_corrected * Ts).eval());
    qref = dq_wb * qref;
    qref.normalize();

    const Vector3 w = last_gyr_bias_corrected;
    MatrixBaseN& F_AA = F_AA_scratch_;
    F_AA.setIdentity();
    [[maybe_unused]] Matrix3 Bstep;
    if constexpr (GB) {
        Matrix3 Rstep;
        filter_detail::rot_and_B_from_wt(w, Ts, Rstep, Bstep);
        F_AA.template topLeftCorner<3,3>() = Rstep;
        F_AA.template block<3,3>(0,3) = Bstep;
    } else {
        const T omega = w.norm();
        const T theta = omega * Ts;
        const Matrix3 I3 = Matrix3::Identity();
        if (theta < T(1e-5)) {
            const Matrix3 Wx = skew_(w);
            F_AA.template topLeftCorner<3,3>() = I3 - Wx * Ts + (Wx * Wx) * (Ts * Ts / T(2));
        } else {
            const Matrix3 Wn = skew_(w / (omega + std::numeric_limits<T>::epsilon()));
            const T s = std::sin(theta), c = std::cos(theta);
            F_AA.template topLeftCorner<3,3>() = I3 - s * Wn + (T(1) - c) * (Wn * Wn);
        }
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
        if (filter_detail::is_isotropic3(Qg)) I_R = Matrix3::Identity() * (Qg(0,0) * Ts);
        else I_R = filter_detail::simpson_R_Q_RT(w, Ts, Qg);
        Matrix3 I_BB = Matrix3::Zero();
        if constexpr (GB) I_BB = simpson_B_Q_BT_(w, Ts, Qbg, Bstep);
        Q_AA.template topLeftCorner<3,3>() = I_R + I_BB;
        if constexpr (GB) {
            Matrix3 IB;
            filter_detail::integral_B_ds(w, Ts, IB);
            const Matrix3 Qtb = IB * Qbg;
            Q_AA.template topRightCorner<3,3>() = Qtb;
            Q_AA.template bottomLeftCorner<3,3>() = Qtb.transpose();
            Q_AA.template bottomRightCorner<3,3>() = Qbg * Ts;
        }
        Q_AA = T(0.5) * (Q_AA + Q_AA.transpose()).eval();
        filter_detail::project_psd_symmetrize_large<T,BASE_N>(Q_AA, T(1e-12));
    }
    // Every block below is written to both triangles with one value, so the
    // propagated covariance is exactly symmetric.  Unit rows of F_AA are
    // applied as copies.
    {
        MatrixBaseN& tmpAA = tmpAA_scratch_;
        for (int i = 0; i < BASE_N; ++i)
            for (int j = 0; j < BASE_N; ++j) {
                if (aa_unit_row_(i)) { tmpAA(i,j) = Pext(i,j); continue; }
                T acc = T(0);
                for (int k = 0; k < BASE_N; ++k) acc += F_AA(i,k) * Pext(k,j);
                tmpAA(i,j) = acc;
            }
        for (int i = 0; i < BASE_N; ++i)
            for (int j = i; j < BASE_N; ++j) {
                T acc;
                if (aa_unit_row_(j)) {
                    acc = Q_AA(i,j) + tmpAA(i,j);
                } else {
                    acc = Q_AA(i,j);
                    for (int k = 0; k < BASE_N; ++k) acc += tmpAA(i,k) * F_AA(j,k);
                }
                Pext(i,j) = acc;
                Pext(j,i) = acc;
            }
    }

    // Linear block: block-diagonal F_LL applied structurally.
    refresh_horizontal_discretization_(Ts);
    refresh_vertical_discretization_(Ts);
    const MatrixHA& F_H = disc_h_.Phi();
    const MatrixZZ& F_Z = disc_z_.Phi();

    {
        T in[NL], out[NL];
        for (int k = 0; k < NL; ++k) in[k] = xext(OFF_L + k);
        apply_F_LL_(in, out);
        for (int k = 0; k < NL; ++k) xext(OFF_L + k) = out[k];
    }

    // P_LL <- F_LL P_LL F_LL^T over the upper block triangle (x, y, z).
    constexpr int LX = OFF_L, LY = OFF_L + NHA;
    propagate_block_<NHA,NHA>(LX, LX, F_H, F_H, true);
    propagate_block_<NHA,NHA>(LX, LY, F_H, F_H, false);
    propagate_block_<NHA,NZ>(LX, OFF_Z, F_H, F_Z, false);
    propagate_block_<NHA,NHA>(LY, LY, F_H, F_H, true);
    propagate_block_<NHA,NZ>(LY, OFF_Z, F_H, F_Z, false);
    propagate_block_<NZ,NZ>(OFF_Z, OFF_Z, F_Z, F_Z, true);
    // Q_H = Sigma_h (x) Q_1.
    const MatrixHA& Q_H1 = disc_h_.Q();
    const MatrixZZ& Q_Z = disc_z_.Q();
    for (int ax = 0; ax < 2; ++ax)
        for (int bx = 0; bx < 2; ++bx) {
            const T s = Sigma_h_(ax, bx);
            if (s == T(0)) continue;
            for (int a = 0; a < NHA; ++a)
                for (int b = 0; b < NHA; ++b)
                    Pext(OFF_L + ax * NHA + a, OFF_L + bx * NHA + b) += s * Q_H1(a, b);
        }
    for (int a = 0; a < NZ; ++a)
        for (int b = 0; b < NZ; ++b)
            Pext(OFF_Z + a, OFF_Z + b) += Q_Z(a, b);

    // P_AL <- F_AA P_AL F_LL^T.
    Eigen::Matrix<T,BASE_N,NL>& tmpAL = tmpAL_scratch_;
    for (int i = 0; i < BASE_N; ++i)
        for (int j = 0; j < NL; ++j) {
            if (aa_unit_row_(i)) { tmpAL(i,j) = Pext(i, OFF_L + j); continue; }
            T acc = T(0);
            for (int k = 0; k < BASE_N; ++k) acc += F_AA(i,k) * Pext(k, OFF_L + j);
            tmpAL(i,j) = acc;
        }
    for (int i = 0; i < BASE_N; ++i) {
        for (int ax = 0; ax < 2; ++ax)
            for (int j = 0; j < NHA; ++j) {
                T acc = T(0);
                for (int k = 0; k < NHA; ++k) acc += F_H(j,k) * tmpAL(i, ax * NHA + k);
                Pext(i, OFF_L + ax * NHA + j) = acc;
                Pext(OFF_L + ax * NHA + j, i) = acc;
            }
        for (int j = 0; j < NZ; ++j) {
            T acc = T(0);
            for (int k = 0; k < NZ; ++k) acc += F_Z(j,k) * tmpAL(i, NH + k);
            Pext(i, OFF_Z + j) = acc;
            Pext(OFF_Z + j, i) = acc;
        }
    }

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
                T acc;
                if (aa_unit_row_(i)) {
                    acc = Pext(i, OFF_BA + j);
                } else {
                    acc = T(0);
                    for (int k = 0; k < BASE_N; ++k) acc += F_AA(i,k) * Pext(k, OFF_BA + j);
                }
                tmpAB(i,j) = acc * phi_b;
            }
        Pext.template block<BASE_N,3>(0, OFF_BA) = tmpAB;
        Pext.template block<3,BASE_N>(OFF_BA, 0) = tmpAB.transpose();

        // P_LB <- F_LL P_LB phi_b.
        for (int j = 0; j < 3; ++j) {
            T in[NL], out[NL];
            for (int k = 0; k < NL; ++k) in[k] = Pext(OFF_L + k, OFF_BA + j);
            apply_F_LL_(in, out);
            for (int k = 0; k < NL; ++k) {
                Pext(OFF_L + k, OFF_BA + j) = out[k] * phi_b;
                Pext(OFF_BA + j, OFF_L + k) = out[k] * phi_b;
            }
        }
    }

    if constexpr (SH || SZ) {
        if (filter_detail::periodic_update_due(Ts, pseudo_update_period_s_, pseudo_update_elapsed_s_)) {
            applyIntegralZeroPseudoMeas();
        }
    }
}

CARMA3_TPL
void CARMA3_CLS::measurement_update_acc_only(Vector3 const& acc_meas, T tempC) {
    last_acc_diag_ = MeasDiag3{};
    if (!acc_meas.allFinite() || !std::isfinite(tempC)) return;
    [[maybe_unused]] const bool use_ba = (AB && acc_bias_updates_enabled_);

    const Vector3 g_world(0, 0, gravity_magnitude_);
    Vector3 ba_term = Vector3::Zero();
    if constexpr (AB) ba_term = xext.template segment<3>(OFF_BA) + k_a_ * (tempC - tempC_ref);

    const Matrix3 Rwb = R_wb();
    const Vector3 aw(xext(idx_a(0)), xext(idx_a(1)), xext(OFF_Z + 2));
    const Vector3 f_cog_b = Rwb * (aw - g_world);
    Vector3 f_pred = f_cog_b + ba_term;
    // Off-CoG IMU: the rotational term at the sensor (ImuLeverArm.h).  Its
    // attitude Jacobian is zero; its gyro-bias Jacobian is appended below.
    if (lever_.enabled()) f_pred += lever_.acceleration();
    const Vector3 r = acc_meas - f_pred;
    last_acc_diag_.r = r;

    const Matrix3 J_att = -skew_(f_cog_b);
    constexpr int NCOL_MAX = 6 + (AB ? 3 : 0) + (GB ? 3 : 0);
    int ncol = 6 + (AB ? 3 : 0);
    int col[NCOL_MAX];
    Vector3 cvec[NCOL_MAX];
    for (int k = 0; k < 3; ++k) { col[k] = k; cvec[k] = J_att.col(k); }
    col[3] = idx_a(0);   cvec[3] = Rwb.col(0);
    col[4] = idx_a(1);   cvec[4] = Rwb.col(1);
    col[5] = OFF_Z + 2;  cvec[5] = Rwb.col(2);
    if constexpr (AB) {
        for (int k = 0; k < 3; ++k) { col[6 + k] = OFF_BA + k; cvec[6 + k] = Vector3::Unit(k); }
    }
    if constexpr (GB) {
        if (lever_.enabled()) {
            const Matrix3 J_bg = lever_.bias_jacobian();
            for (int k = 0; k < 3; ++k) { col[ncol] = 3 + k; cvec[ncol] = J_bg.col(k); ++ncol; }
        }
    }

    MatrixNX3& PCt = PCt_scratch_;
    for (int i = 0; i < NX; ++i) {
        T pc0 = T(0), pc1 = T(0), pc2 = T(0);
        for (int c = 0; c < ncol; ++c) {
            const T v = Pext(i, col[c]);
            pc0 += v * cvec[c](0); pc1 += v * cvec[c](1); pc2 += v * cvec[c](2);
        }
        PCt(i,0) = pc0; PCt(i,1) = pc1; PCt(i,2) = pc2;
    }
    Matrix3& S_mat = S_scratch_;
    S_mat = Racc;
    for (int c = 0; c < ncol; ++c)
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) S_mat(a,b) += cvec[c](a) * PCt(col[c], b);
    for (int a = 0; a < 3; ++a)
        for (int b = a + 1; b < 3; ++b) {
            const T v = T(0.5) * (S_mat(a,b) + S_mat(b,a));
            S_mat(a,b) = v;
            S_mat(b,a) = v;
        }

    if constexpr (AB) {
        if (!use_ba) {
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

CARMA3_TPL
void CARMA3_CLS::measurement_update_mag_only(const Vector3& mag_meas) {
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

CARMA3_TPL
template <int M>
void CARMA3_CLS::select_update_(const int (&idx)[M], const Vector3& r, const Matrix3& R, bool freeze_ba) {
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

CARMA3_TPL
void CARMA3_CLS::applyIntegralZeroPseudoMeas() {
    const bool freeze = AB && !acc_bias_updates_enabled_;
    if constexpr (SH && SZ) {
        const int idx[3] = {idx_S(0), idx_S(1), OFF_SZ};
        const Vector3 r(-xext(idx[0]), -xext(idx[1]), -xext(idx[2]));
        select_update_<3>(idx, r, R_S, freeze);
    } else if constexpr (SH) {
        const int idx[2] = {idx_S(0), idx_S(1)};
        const Vector3 r(-xext(idx[0]), -xext(idx[1]), T(0));
        select_update_<2>(idx, r, R_S, freeze);
    } else if constexpr (SZ) {
        const int idx[1] = {OFF_SZ};
        Matrix3 R = Matrix3::Identity();
        R(0,0) = R_S(2,2);
        select_update_<1>(idx, Vector3(-xext(OFF_SZ), T(0), T(0)), R, freeze);
    }
}

CARMA3_TPL
void CARMA3_CLS::measurement_update_position_pseudo(const Vector3& p_meas, const Vector3& sigma_meas) {
    const int idx[3] = {idx_p(0), idx_p(1), OFF_Z};
    const Vector3 r = p_meas - get_position();
    const Vector3 s = sigma_meas.cwiseMax(T(0));
    select_update_<3>(idx, r, Matrix3(s.array().square().matrix().asDiagonal()), false);
}

CARMA3_TPL
void CARMA3_CLS::measurement_update_velocity_pseudo(const Vector3& v_meas, const Vector3& sigma_meas) {
    const int idx[3] = {idx_v(0), idx_v(1), OFF_Z + 1};
    const Vector3 r = v_meas - get_velocity();
    const Vector3 s = sigma_meas.cwiseMax(T(0));
    select_update_<3>(idx, r, Matrix3(s.array().square().matrix().asDiagonal()), false);
}

CARMA3_TPL
void CARMA3_CLS::measurement_update_vert_velocity_pseudo(T vz_meas, T sigma_meas) {
    if (!std::isfinite(vz_meas)) return;
    const int idx[1] = {OFF_Z + 1};
    const T s = std::max(T(0), sigma_meas);
    Matrix3 R = Matrix3::Identity();
    R(0,0) = s * s;
    select_update_<1>(idx, Vector3(vz_meas - xext(OFF_Z + 1), T(0), T(0)), R, false);
}

CARMA3_TPL
void CARMA3_CLS::applyQuaternionCorrectionFromErrorState() {
    project_gyro_bias_();
    const Vector3 dtheta = xext.template segment<3>(0);
    if (!dtheta.allFinite()) {
        xext.template head<3>().setZero();
        return;
    }
    qref = filter_detail::quat_from_delta_theta(dtheta) * qref;
    qref.normalize();
    filter_detail::apply_left_error_reset<T, NX>(Pext, dtheta);
    xext.template head<3>().setZero();
    project_acc_bias_();
}

CARMA3_TPL
void CARMA3_CLS::project_acc_bias_() {
    if constexpr (AB) {
        if (!(acc_bias_limit_ > T(0))) return;
        auto b = xext.template segment<3>(OFF_BA);
        if (!b.allFinite()) { b.setZero(); return; }
        const T n = b.norm();
        if (n > acc_bias_limit_) b *= (acc_bias_limit_ / n);
    }
}

#undef CARMA3_TPL
#undef CARMA3_CLS

}  // namespace ocean_imu::kalman
