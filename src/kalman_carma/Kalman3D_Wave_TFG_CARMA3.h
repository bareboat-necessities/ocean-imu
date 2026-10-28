#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Kalman3D_Wave_TFG_CARMA3 -- EXPERIMENTAL three-axis research filter.  Not a
  shipping estimator: Kalman3D_Wave_OU_III remains the deployed, authoritative
  implementation, and nothing in the deployed configuration includes this
  header; see docs/tfg-carma3-research.md.

  The motion prior of Kalman3D_Wave_CARMA3 (one CARMA(NC, NC-3) shape, frequency
  scale and exact discretization on all three world axes, no S state) on the
  right-invariant two-frame Lie group of src/lie/TwoFrameGroup.h, replacing
  CARMA3's quaternion MEKF error.

  Nominal state X = (R, W, B), TwoFrameGroup<T, NC, 2>:

      R = R_bw                       body -> world (NED, z down)
      W = [p v a l_1 .. l_{NC-3}]    3 x NC, world frame; row k is the CARMA
                                     state of axis k in the library coordinates
                                     of CarmaVerticalModel.h
      B = [b_g b_a]                  3 x 2, body frame

  Convention contract (TwoFrameGroup.h): R = R_bw, right-invariant error
  eta = Xhat X^-1 = Exp_G(xi), left correction Xhat+ = Exp_G(-xihat) Xhat.

      xi = [ phi | beta_g | rho_p rho_v rho_a rho_l.. | beta_a ]
             0     3        6     9     12    15        6 + 3 NC
      dR = Rhat R^T = Exp(phi),  rho_i = what_i - dR w_i,  beta = R (bhat - b)

  Continuous model (unit-free, gyro noise n_g, gyro-bias walk n_bg, CARMA
  driving noise xi_k per world axis with spatial covariance Sigma_s):

      R'   = R [w_m - b_g - n_g]x
      W'   = W A^T + g e_a xi^T            (x_k' = A x_k + g e_a xi_k)
      b_g' = n_bg
      b_a' = -b_a / tau_b + n_ba            (active);  b_a' = 0 (held)

  Linearized invariant error dynamics (w = w_m - bhat_g, w_w = Rhat w):

      phi'    = -beta_g + Rhat n_g
      beta_g' = [w_w]x beta_g - Rhat n_bg
      rho_i'  = sum_j A_ij rho_j - [what_i]x beta_g + [what_i]x Rhat n_g - g e_a,i xi
      beta_a' = [w_w]x beta_a - beta_a / tau_b - Rhat n_ba

  There is no attitude -> translation term: with one A on every axis the
  column mixing commutes with every rotation, so A (phi x w) = phi x (A w)
  cancels exactly.  The same identity makes the structured impulse response
  collapse: Phi_A(h-t) [W(t)]x = [W(h)]x because Phi_A(h-t) Phi_A(t) =
  Phi_A(h), so over one step with constant measured rate

      Phi_{phi,bg}  = F  = -Rhat0 h J_l(w h) Rhat0^T
      Phi_{bg,bg}   = E  = Rhat1 Rhat0^T = Exp(w_w h)
      Phi_{rho,rho} = Phi_A (x) I_3        (exact CARMA discretization)
      Phi_{rho,bg}  = [W(h)]x F
      Phi_{ba,ba}   = phi_b E,  phi_b = exp(-h/tau_b) active, 1 held

  and every gyro and gyro-bias noise column of rho is [W(h)]x times the
  matching attitude column, so the discrete process noise needs only 3 x 3
  integrals (five-point Gauss-Legendre where no closed form is used).  This
  replaces the structured quadrature of Kalman3D_Wave_TFG for this prior; the
  test checks both against a fine-step integration of the continuous model.

  The model is not group-affine (the beta_g -> rho term depends on the
  estimate), so this is a right-invariant two-frame Lie-group error-state
  EKF, not a log-linear invariant EKF.  Its accelerometer and magnetometer
  Jacobians are state independent:

      r_a = Rhat (f_m - bhat_a - k_a dT) - (ahat - g),   H_a = [ [g]x 0 .. -I(rho_a) .. -I(beta_a) ]
      r_m = Rhat m - B_w,                               H_m = [ -[B_w]x 0 .. ]

  with noise Rhat R Rhat^T.  A finite correction a = K r is injected as
  Exp_G(-a) Xhat and the covariance transported by the full group right
  Jacobian J_r(a) = d Log_G(Exp_G(-a) Exp_G(a + d)) / d d.

  Policy parity with Kalman3D_Wave_CARMA3.  The fusion wrapper is the same code
  for both filters (SeaStateFusionFilter_TFG_CARMA3.h).  Every covariance
  *write* that CARMA3 performs in its own coordinates (attitude seeds, the
  stationary CARMA prior, freezing the accelerometer bias, absolute attitude
  writes) is performed here in the CARMA3-equivalent coordinates

      e = [ dtheta_b, db_g, dw_0 .. dw_{NC-1}, db_a ],   C_true = Exp(dtheta_b) Chat,
                                                          dx = x_true - xhat
      phi = Rhat dtheta_b,  beta = -Rhat db,  rho_i = -dw_i + [what_i]x Rhat dtheta_b

  and mapped back, so both filters start every phase with the same physical
  uncertainty (CovarianceWriteMode::Physical, the default).  The bias
  projections, the frozen-bias measurement semantics and the absence of
  innovation gating are CARMA3's.  CovarianceWriteMode::Native and the
  ResetMode values exist for the controlled ablations of the research report.

  One CARMA model serves every axis: the vertical and horizontal frequency
  scales and shapes must agree (the wrapper's mapping gives them the same
  value); the axes may differ in driving amplitude and the horizontal pair may
  be correlated through Sigma_s.  A different horizontal scale is recorded and
  ignored (anisotropic dynamics would add an attitude -> translation term; see
  the report), never silently approximated.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

#include "kalman_common/MekfAttitudeMath.h"
#include "kalman_carma/CarmaVerticalModel.h"
#include "lie/SO3Jacobians.h"
#include "lie/TwoFrameGroup.h"

namespace ocean_imu::kalman {

namespace tfg_carma_detail {

/*
    Off-diagonal (vector <- attitude) block of the right Jacobian of the
    two-frame group.  Every vector block (world column or world-referred body
    bias) has ad_xi = [[ [phi]x, 0 ], [ [z]x, [phi]x ]], so

        J_r = sum_n (-ad)^n / (n+1)!  =  [[ J_l(-phi), 0 ], [ Q_r(phi, z), J_l(-phi) ]],
        Q_r(phi, z) = Q_l(-phi, -z),

    Q_l being the SE(3) left-Jacobian block (Barfoot, State Estimation for
    Robotics, eq. 7.86).  Its scalar coefficients are evaluated from their
    alternating series below |phi| = 1, where the closed forms cancel in
    float.
*/
template<typename T>
inline Eigen::Matrix<T,3,3> right_jacobian_cross(const Eigen::Matrix<T,3,1>& phi,
                                                 const Eigen::Matrix<T,3,1>& z) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    using Vector3 = Eigen::Matrix<T,3,1>;
    const Matrix3 P = lie::skew<T>(Vector3(-phi));
    const Matrix3 Z = lie::skew<T>(Vector3(-z));
    const T t2 = phi.squaredNorm();
    T c1, c2, c3;
    if (t2 < T(1)) {
        // c1 = sum (-1)^k t^2k/(2k+3)!, c2 = sum (-1)^k t^2k/(2k+4)!,
        // c3 = sum (-1)^k (k+1) t^2k/(2k+5)!; nine terms leave < 1e-19.
        static constexpr double f3[9] = {1.0/6.0, 1.0/120.0, 1.0/5040.0, 1.0/362880.0,
            1.0/39916800.0, 1.0/6227020800.0, 1.0/1307674368000.0, 1.0/355687428096000.0,
            1.0/121645100408832000.0};
        static constexpr double f4[9] = {1.0/24.0, 1.0/720.0, 1.0/40320.0, 1.0/3628800.0,
            1.0/479001600.0, 1.0/87178291200.0, 1.0/20922789888000.0, 1.0/6402373705728000.0,
            1.0/2432902008176640000.0};
        static constexpr double f5[9] = {1.0/120.0, 2.0/5040.0, 3.0/362880.0, 4.0/39916800.0,
            5.0/6227020800.0, 6.0/1307674368000.0, 7.0/355687428096000.0,
            8.0/121645100408832000.0, 9.0/51090942171709440000.0};
        c1 = T(f3[8]); c2 = T(f4[8]); c3 = T(f5[8]);
        for (int k = 7; k >= 0; --k) {
            c1 = T(f3[k]) - t2 * c1;
            c2 = T(f4[k]) - t2 * c2;
            c3 = T(f5[k]) - t2 * c3;
        }
    } else {
        const T t = std::sqrt(t2);
        const T s = std::sin(t), c = std::cos(t);
        c1 = (t - s) / (t2 * t);
        c2 = (t2 + T(2) * c - T(2)) / (T(2) * t2 * t2);
        c3 = (T(2) * t - T(3) * s + t * c) / (T(2) * t2 * t2 * t);
    }
    const Matrix3 PZ = P * Z, ZP = Z * P, PZP = PZ * P;
    return T(0.5) * Z + c1 * (PZ + ZP + PZP)
         + c2 * (P * PZ + ZP * P - T(3) * PZP)
         + c3 * (PZP * P + P * PZP);
}

// Five-point Gauss-Legendre nodes and weights on [-1, 1].
inline constexpr double kGL5Nodes[5] = {
    -0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640};
inline constexpr double kGL5Weights[5] = {
    0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665,
    0.2369268850561891};

}  // namespace tfg_carma_detail

template <typename T = float, int NC = 4>
class Kalman3D_Wave_TFG_CARMA3 {
    static_assert(NC >= 3 && NC <= 5, "CARMA order 3..5");
public:
    using Group = lie::TwoFrameGroup<T, NC, 2>;
    static constexpr int NX = Group::NX;          // 3 + 3 NC + 6
    static constexpr int NBLK = NX / 3;           // 3-vector blocks
    static constexpr int OFF_PHI = Group::OFF_PHI;
    static constexpr int OFF_BG = Group::OFF_BG;
    static constexpr int OFF_W = Group::OFF_W;
    static constexpr int OFF_BA = Group::OFF_BA;
    static constexpr int IP = 0, IV = 1, IA = 2;  // world columns
    static constexpr int OFF_P = Group::world_offset(IP);
    static constexpr int OFF_V = Group::world_offset(IV);
    static constexpr int OFF_A = Group::world_offset(IA);
    static constexpr int BLK_PHI = 0, BLK_BG = 1, BLK_W = 2, BLK_BA = 2 + NC;
    // Index of state `axis` of world column `col` (axis 0 = x, 1 = y, 2 = z).
    static constexpr int world_index(int col, int axis) { return OFF_W + 3 * col + axis; }
    static constexpr int pos_index(int axis) { return world_index(IP, axis); }
    static constexpr int vel_index(int axis) { return world_index(IV, axis); }

    using Vector3 = Eigen::Matrix<T,3,1>;
    using Matrix3 = Eigen::Matrix<T,3,3>;
    using MatrixNX = Eigen::Matrix<T,NX,NX>;
    using MatrixNX3 = Eigen::Matrix<T,NX,3>;
    using Tangent = Eigen::Matrix<T,NX,1>;
    using MatC = Eigen::Matrix<T,NC,NC>;
    using Model = carma::CarmaVerticalModel<T,NC>;

    static constexpr T STD_GRAVITY = T(9.80665);
    static constexpr T tempC_ref = T(35.0);

    // Ablations of the research report; the defaults are the filter.
    enum class ResetMode : uint8_t {
        Full,          // full two-frame right Jacobian J_r(a)
        AttitudeOnly,  // J_l(-phi) on the attitude block only (MEKF-like)
        None           // no transport
    };
    enum class CovarianceWriteMode : uint8_t {
        Physical,      // CARMA3's covariance writes in CARMA3-equivalent coordinates
        Native         // the same writes applied directly to the invariant covariance
    };

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    struct MeasDiag3 {
        Vector3 r = Vector3::Zero();   // world-frame residual
        Matrix3 S = Matrix3::Zero();
        T nis = std::numeric_limits<T>::quiet_NaN();
        bool accepted = false;
    };

    const MeasDiag3& lastAccDiag() const noexcept { return last_acc_diag_; }
    const MeasDiag3& lastMagDiag() const noexcept { return last_mag_diag_; }

    Kalman3D_Wave_TFG_CARMA3(Vector3 const& sigma_a, Vector3 const& gyro_noise_density_rad_sqrt_s,
                             Vector3 const& sigma_m,
                             T Pq0 = T(5e-4), T Pb0 = T(1e-6), T b0 = T(1e-11), T R_S_noise_var = T(1.5),
                             T gravity_magnitude = T(STD_GRAVITY),
                             const carma::CarmaShape& shape = carma::CarmaDefaultShape<NC>::value);

    bool initialize_from_acc_mag(Vector3 const& acc, Vector3 const& mag);
    bool initialize_from_acc(Vector3 const& acc);
    bool initialize_from_acc_preserve_yaw(Vector3 const& acc);
    bool initialize_from_attitude(Eigen::Quaternion<T> const& q_bw, T tilt_sigma_rad, T yaw_sigma_rad);
    static Eigen::Quaternion<T> quaternion_from_acc(Vector3 const& acc);

    void set_mag_world_ref(const Vector3& B_world) { B_w_ = B_world; }
    [[nodiscard]] const Vector3& mag_world_ref() const noexcept { return B_w_; }

    void time_update(Vector3 const& gyr, T Ts);
    void measurement_update_acc_only(Vector3 const& acc, T tempC = tempC_ref);
    void measurement_update_mag_only(Vector3 const& mag);

    [[nodiscard]] Eigen::Quaternion<T> quaternion() const {
        Eigen::Quaternion<T> q(X_.R);
        q.normalize();
        return q;
    }
    [[nodiscard]] Eigen::Quaternion<T> quaternion_boat() const { return quaternion(); }
    [[nodiscard]] const MatrixNX& covariance_full() const { return P_; }
    [[nodiscard]] MatrixNX& covariance_full() { return P_; }
    [[nodiscard]] const Group& state() const { return X_; }
    [[nodiscard]] Group& state() { return X_; }

    // Absolute attitude write keeping the world columns and biases (CARMA3's
    // set_quaternion_boat).  The covariance follows the write mode.
    void set_quaternion_boat(const Eigen::Quaternion<T>& q_bw);

    [[nodiscard]] Vector3 gyroscope_bias() const { return X_.B.col(0); }
    [[nodiscard]] Vector3 get_acc_bias() const { return X_.B.col(1); }
    bool acc_bias_updates_enabled() const { return acc_bias_updates_enabled_; }
    void set_acc_bias_updates_enabled(bool en);

    [[nodiscard]] Vector3 get_position() const { return X_.X.col(IP); }
    [[nodiscard]] Vector3 get_velocity() const { return X_.X.col(IV); }
    [[nodiscard]] Vector3 get_world_accel() const { return X_.X.col(IA); }

    // Covariance in the CARMA3-equivalent coordinates e (see the header).
    [[nodiscard]] MatrixNX covariance_physical() const {
        MatrixNX P = P_, scratch;
        map_coordinates_(P, scratch, false);
        return P;
    }
    // Cov([dp; dv]) of the additive world errors, 6 x 6 ([x y z] of p, then v).
    [[nodiscard]] Eigen::Matrix<T,6,6> physical_pv_covariance() const;
    [[nodiscard]] Matrix3 position_covariance() const {
        return physical_pv_covariance().template topLeftCorner<3,3>();
    }

    [[nodiscard]] T get_aw_time_constant() const { return tau_aw_; }
    [[nodiscard]] Vector3 get_aw_stationary_std() const {
        return Sigma_s_.diagonal().cwiseMax(T(0)).cwiseSqrt();
    }
    [[nodiscard]] Vector3 get_RS_noise_std() const { return R_S_.diagonal().cwiseMax(T(0)).cwiseSqrt(); }
    void set_aw_time_constant(T tau_seconds) { tau_aw_ = std::max(T(1e-3), tau_seconds); }
    // No S state: the wrapper's cadence is kept for interface parity only.
    void set_pseudo_update_period_s(T period_s) {
        if (!(period_s > T(0)) || !std::isfinite(period_s)) return;
        pseudo_update_period_s_ = std::max(T(1e-4), period_s);
    }
    [[nodiscard]] T get_pseudo_update_period_s() const { return pseudo_update_period_s_; }

    // Stationary CARMA prior covariance Sigma_s (x) P_inf, decorrelated from
    // everything else (CARMA3's reset_aw_covariance_to_stationary).
    void reset_aw_covariance_to_stationary();

    // Vertical amplitude and the common frequency scale.
    bool set_vertical_carma(T kappa, T accel_std);
    // Horizontal amplitudes [[sx^2, r sx sy], [r sx sy, sy^2]]; kappa must
    // equal the common scale.
    bool set_horizontal_carma(T kappa, T sx, T sy, T r = T(0));
    // One shape serves every axis.
    bool set_carma_shape(const carma::CarmaShape& shape) { return model_.set_shape(shape); }
    bool set_horizontal_shape(const carma::CarmaShape& shape) { return same_shape_(shape); }
    bool set_vertical_shape(const carma::CarmaShape& shape) { return same_shape_(shape); }
    [[nodiscard]] const Model& vertical_model() const noexcept { return model_; }
    [[nodiscard]] const Model& horizontal_model() const noexcept { return model_; }
    [[nodiscard]] Eigen::Matrix<T,2,2> horizontal_sigma() const {
        return Sigma_s_.template topLeftCorner<2,2>();
    }
    [[nodiscard]] const Matrix3& spatial_sigma() const noexcept { return Sigma_s_; }
    // Count of horizontal configurations whose frequency scale differed from
    // the common one (not representable here; see the header).
    [[nodiscard]] uint32_t rejected_anisotropic_configurations() const noexcept { return kappa_rejects_; }

    void set_RS_noise(const Vector3& sigma_S) { R_S_ = sigma_S.array().square().matrix().asDiagonal(); }
    void set_Racc_std(const Vector3& sigma_acc) { Racc_ = sigma_acc.array().square().matrix().asDiagonal(); }
    void set_Rmag(const Vector3& sigma_mag) { Rmag_ = sigma_mag.array().square().matrix().asDiagonal(); }

    void set_initial_acc_bias_std(T s);
    void set_Q_bacc_rw(const Vector3& driving_std_per_sqrt_s) {
        Q_bacc_ = driving_std_per_sqrt_s.array().square().matrix().asDiagonal();
    }
    void set_acc_bias_time_constant(T tau_seconds) { tau_bacc_ = std::max(T(1e-3), tau_seconds); }
    [[nodiscard]] T get_acc_bias_time_constant() const noexcept { return tau_bacc_; }
    [[nodiscard]] Vector3 get_acc_bias_at_temperature(T tempC) const {
        return X_.B.col(1) + k_a_ * (tempC - tempC_ref);
    }
    [[nodiscard]] Vector3 get_acc_bias_body_at_temperature(T tempC) const { return get_acc_bias_at_temperature(tempC); }
    void set_initial_acc_bias(const Vector3& b0) { X_.B.col(1) = b0; }
    void set_accel_bias_temp_coeff(const Vector3& ka_per_degC) { k_a_ = ka_per_degC; }
    void set_accel_bias_limit(T limit_mps2) { acc_bias_limit_ = limit_mps2; }
    [[nodiscard]] T accel_bias_limit() const { return acc_bias_limit_; }

    void set_gyro_noise_density_rad_sqrt_s(const Vector3& density) {
        const Vector3 d = density.cwiseAbs();
        Q_gyro_ = d.array().square().matrix().asDiagonal();
    }
    // Gyro-bias random-walk variance density per axis, (rad/s)^2/s.
    void set_Q_bgyro_rw(const Vector3& var_per_s) { Q_bg_ = var_per_s.cwiseAbs().asDiagonal(); }
    // The attitude/bias process noise is always the exact one; the flag is
    // kept for interface parity.
    void set_exact_att_bias_Qd(bool on) { use_exact_att_bias_Qd_ = on; }

    void set_reset_mode(ResetMode m) { reset_mode_ = m; }
    [[nodiscard]] ResetMode reset_mode() const { return reset_mode_; }
    void set_covariance_write_mode(CovarianceWriteMode m) { write_mode_ = m; }
    [[nodiscard]] CovarianceWriteMode covariance_write_mode() const { return write_mode_; }

    // Diagnostics: innovation-covariance factorizations that needed the
    // diagonal bump, and corrections skipped as non-finite.
    [[nodiscard]] uint32_t ldlt_bumps() const noexcept { return ldlt_bumps_; }
    [[nodiscard]] uint32_t skipped_updates() const noexcept { return skipped_updates_; }

    // ---- Building blocks, public for the mathematical tests ------------------

    // Transition of the invariant error over one step of measured rate gyr from
    // the current nominal state (before the mean propagation).
    void build_transition(const Vector3& gyr, T Ts, MatrixNX& Phi) const;
    // Discrete process noise of the same step.
    void build_process_noise(const Vector3& gyr, T Ts, MatrixNX& Qd) const;
    // J_r(a) (Full), J_l(-phi) on the attitude block (AttitudeOnly) or I.
    void build_reset_jacobian(const Tangent& a, MatrixNX& J) const;
    // Inject the estimated error a: Xhat <- Exp_G(-a) Xhat.
    void inject(const Tangent& a);
    // Right-invariant error of this filter's state relative to a reference
    // (the truth), xi = Log_G(Xhat Xref^-1).
    [[nodiscard]] Tangent error_from(const Group& reference) const { return X_.local(reference); }
    // P <- T P T^T (to_invariant) or T^-1 P T^-T: e <-> xi at the current
    // nominal state, in place through `scratch`.
    void map_coordinates_(MatrixNX& P, MatrixNX& scratch, bool to_invariant) const;

private:
    const T gravity_magnitude_ = T(STD_GRAVITY);

    Group X_{};
    MatrixNX P_;
    Vector3 B_w_ = Vector3::UnitX();
    Vector3 last_gyr_bias_corrected_ = Vector3::Zero();

    T sigma_bacc0_ = T(0.004);
    Matrix3 Q_bacc_ = Matrix3::Identity() * T(2.5e-7);
    T tau_bacc_ = T(5000.0);
    Vector3 k_a_ = Vector3::Constant(T(0.002));
    T acc_bias_limit_ = T(0.4);

    Matrix3 Q_gyro_;
    Matrix3 Q_bg_;
    Matrix3 Racc_;
    Matrix3 Rmag_;
    Matrix3 R_S_;

    T tau_aw_ = T(2.1);
    T pseudo_update_period_s_ = T(0.015);

    // One CARMA model at unit acceleration variance and its discretization;
    // Sigma_s scales the driving noise per world axis.
    Model model_;
    Matrix3 Sigma_s_ = Matrix3::Identity();
    carma::StepDiscretization<T,NC> disc_;
    uint32_t kappa_rejects_ = 0;

    bool acc_bias_updates_enabled_ = true;
    bool use_exact_att_bias_Qd_ = true;
    ResetMode reset_mode_ = ResetMode::Full;
    CovarianceWriteMode write_mode_ = CovarianceWriteMode::Physical;
    uint32_t ldlt_bumps_ = 0;
    uint32_t skipped_updates_ = 0;

    // Scratch kept off the task stack.
    MatrixNX M_scratch_;
    MatrixNX3 PHt_scratch_;
    MatrixNX3 K_scratch_;

    MeasDiag3 last_acc_diag_;
    MeasDiag3 last_mag_diag_;

    struct StepBlocks {
        Matrix3 F;              // Phi_{phi,bg}
        Matrix3 E;              // Phi_{bg,bg}
        Matrix3 R1;             // attitude at the end of the step
        Matrix3 Xs[NC];         // [w_i(h)]x
        T phi_b;
    };
    struct StepNoise {
        Matrix3 Qpp, Qpb, Qbb, Qaa;
    };

    bool same_shape_(const carma::CarmaShape& s) const {
        const carma::CarmaShape& m = model_.shape();
        return s.f1_hz == m.f1_hz && s.zeta1 == m.zeta1 && s.f2_hz == m.f2_hz && s.zeta2 == m.zeta2 &&
               s.f_lambda_hz == m.f_lambda_hz && s.fz_hz == m.fz_hz && s.zeta_z == m.zeta_z;
    }

    void reorthonormalize_() {
        Eigen::Quaternion<T> q(X_.R);
        q.normalize();
        X_.R = q.toRotationMatrix();
    }

    void symmetrize_P_() {
        for (int i = 0; i < NX; ++i)
            for (int j = i + 1; j < NX; ++j) {
                const T v = T(0.5) * (P_(i,j) + P_(j,i));
                P_(i,j) = v;
                P_(j,i) = v;
            }
    }

    void project_gyro_bias_() {
        auto bias = X_.B.col(0);
        filter_detail::project_gyro_bias<T>(bias);
    }
    void project_acc_bias_() {
        if (!(acc_bias_limit_ > T(0))) return;
        auto b = X_.B.col(1);
        if (!b.allFinite()) { b.setZero(); return; }
        const T n = b.norm();
        if (n > acc_bias_limit_) b *= (acc_bias_limit_ / n);
    }

    // Covariance writes: P_ is brought into CARMA3-equivalent coordinates at
    // the current state, edited, and mapped back at the (possibly new) state
    // (Physical), or edited as it is (Native).  In place, no NX x NX local.
    void begin_write_() {
        if (write_mode_ == CovarianceWriteMode::Physical) map_coordinates_(P_, M_scratch_, false);
    }
    void end_write_() {
        if (write_mode_ == CovarianceWriteMode::Physical) map_coordinates_(P_, M_scratch_, true);
        symmetrize_P_();
    }
    // out = T in, or T in^T when `transposed` (T^-1 when !to_invariant).
    void apply_T_rows_(const MatrixNX& in, bool transposed, bool to_invariant, MatrixNX& out) const;

    // CARMA3's body-frame accel-only attitude seed and decorrelation, in e.
    void seed_attitude_(MatrixNX& Pe, T tilt_sigma_rad, T yaw_sigma_rad) const;
    static void decorrelate_attitude_(MatrixNX& Pe);

    void refresh_discretization_(T Ts) {
        if (disc_.advance(model_.generation(), Ts)) return;
        disc_.refresh(model_.A(), model_.W(), model_.generation(), Ts);
    }

    void step_blocks_(const Vector3& w, T Ts, const Matrix3& R0, const MatC& Fc,
                      const Eigen::Matrix<T,3,NC>& W1, StepBlocks& b) const;
    void step_noise_(const Vector3& w, T Ts, const Matrix3& R0, const StepBlocks& b, StepNoise& q) const;
    // out = Phi * in (in read transposed when `transposed`).
    void apply_phi_(const MatrixNX& in, bool transposed, const StepBlocks& b, const MatC& Fc,
                    MatrixNX& out) const;
    void add_noise_(MatrixNX& P, const StepBlocks& b, const StepNoise& q, const MatC& Qc) const;

    bool safe_ldlt3_(Matrix3& S, Eigen::LDLT<Matrix3>& ldlt, T noise_scale) {
        ldlt.compute(S);
        if (ldlt.info() == Eigen::Success) return true;
        ++ldlt_bumps_;
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
    void gain3_(const MatrixNX3& PHt, const Eigen::LDLT<Matrix3>& ldlt, MatrixNX3& K) const {
        Matrix3 Sinv;
        for (int k = 0; k < 3; ++k) {
            Vector3 e = Vector3::Zero();
            e(k) = T(1);
            const Vector3 c = ldlt.solve(e);
            Sinv(0,k) = c(0); Sinv(1,k) = c(1); Sinv(2,k) = c(2);
        }
        for (int i = 0; i < NX; ++i)
            for (int b = 0; b < 3; ++b)
                K(i,b) = PHt(i,0) * Sinv(0,b) + PHt(i,1) * Sinv(1,b) + PHt(i,2) * Sinv(2,b);
    }
    // CARMA3's Joseph form P + K S K^T - K PHt^T - PHt K^T, exactly symmetric.
    void joseph_update_(const MatrixNX3& K, const Matrix3& S, const MatrixNX3& PHt) {
        T Kl[NX][3], G[NX][3];
        for (int i = 0; i < NX; ++i)
            for (int b = 0; b < 3; ++b) {
                Kl[i][b] = K(i,b);
                T g = -PHt(i,b);
                for (int a = 0; a < 3; ++a) g += K(i,a) * S(a,b);
                G[i][b] = g;
            }
        for (int i = 0; i < NX; ++i) {
            const T pi[3] = {PHt(i,0), PHt(i,1), PHt(i,2)};
            for (int j = i; j < NX; ++j) {
                T d = T(0);
                for (int l = 0; l < 3; ++l) d += Kl[i][l] * G[j][l] - pi[l] * Kl[j][l];
                const T v = P_(i,j) + d;
                P_(i,j) = v;
                P_(j,i) = v;
            }
        }
    }
    void freeze_acc_bias_rows_(MatrixNX3& M) const { M.template block<3,3>(OFF_BA, 0).setZero(); }

    // Gain, Joseph update, transport and injection shared by both sensors.
    bool correct3_(const Vector3& r, Matrix3& S, MatrixNX3& PHt, T noise_scale, bool freeze_ba,
                   MeasDiag3& diag);
    void transport_(const Tangent& a);
};

// ------------------------------------------------------------------------
// Implementation
// ------------------------------------------------------------------------

#define TFGC3_TPL template <typename T, int NC>
#define TFGC3_CLS Kalman3D_Wave_TFG_CARMA3<T,NC>

TFGC3_TPL
TFGC3_CLS::Kalman3D_Wave_TFG_CARMA3(
    Vector3 const& sigma_a, Vector3 const& gyro_noise_density_rad_sqrt_s, Vector3 const& sigma_m,
    T Pq0, T Pb0, T b0, T R_S_noise_var, T gravity_magnitude, const carma::CarmaShape& shape)
  : gravity_magnitude_(gravity_magnitude),
    Racc_(sigma_a.array().square().matrix().asDiagonal()),
    Rmag_(sigma_m.array().square().matrix().asDiagonal()),
    model_(shape)
{
    Q_gyro_ = gyro_noise_density_rad_sqrt_s.array().square().matrix().asDiagonal();
    Q_bg_ = Matrix3::Identity() * b0;
    R_S_ = Matrix3::Identity() * R_S_noise_var;

    // CARMA3's constructor covariance; at the identity state the coordinate
    // map only flips signs, so it is written directly.
    P_.setZero();
    P_.template block<3,3>(OFF_PHI, OFF_PHI) = Matrix3::Identity() * Pq0;
    P_.template block<3,3>(OFF_BG, OFF_BG) = Matrix3::Identity() * Pb0;
    P_.template block<3,3>(OFF_BA, OFF_BA) = Matrix3::Identity() * (sigma_bacc0_ * sigma_bacc0_);
    model_.configure(T(1), T(1));
    Sigma_s_.setIdentity();
    reset_aw_covariance_to_stationary();
}

TFGC3_TPL
bool TFGC3_CLS::set_vertical_carma(T kappa, T accel_std) {
    if (!(accel_std > T(0)) || !std::isfinite(accel_std)) return false;
    if (!model_.configure(kappa, T(1))) return false;
    Sigma_s_(2,2) = accel_std * accel_std;
    return true;
}

TFGC3_TPL
bool TFGC3_CLS::set_horizontal_carma(T kappa, T sx, T sy, T r) {
    if (!(sx > T(0)) || !(sy > T(0)) || !std::isfinite(sx) || !std::isfinite(sy)) return false;
    if (!(std::abs(r) <= T(1))) return false;
    if (!(kappa > T(0)) || !std::isfinite(kappa)) return false;
    Sigma_s_(0,0) = sx * sx;
    Sigma_s_(1,1) = sy * sy;
    Sigma_s_(0,1) = Sigma_s_(1,0) = r * sx * sy;
    if (std::abs(kappa - model_.kappa()) > T(1e-6) * model_.kappa()) {
        ++kappa_rejects_;
        return false;
    }
    return true;
}

TFGC3_TPL
void TFGC3_CLS::apply_T_rows_(const MatrixNX& in, bool transposed, bool to_invariant, MatrixNX& out) const {
    // T   = [phi <- dtheta: R], [bg: -R], [w_i: -I, w_i <- dtheta: [w_i]x R], [ba: -R];
    // T^-1 = [dtheta <- phi: R^T], [bg: -R^T], [w_i: -I, w_i <- phi: [w_i]x], [ba: -R^T].
    const Matrix3 A = to_invariant ? Matrix3(X_.R) : Matrix3(X_.R.transpose());
    Matrix3 Sx[NC];
    for (int i = 0; i < NC; ++i) {
        const Matrix3 S = lie::skew<T>(Vector3(X_.X.col(i)));
        Sx[i] = to_invariant ? Matrix3(S * X_.R) : S;
    }
    auto at = [&](int r, int c) { return transposed ? in(c, r) : in(r, c); };
    for (int c = 0; c < NX; ++c) {
        T ph[3], bg[3], ba[3];
        for (int k = 0; k < 3; ++k) {
            ph[k] = at(OFF_PHI + k, c);
            bg[k] = at(OFF_BG + k, c);
            ba[k] = at(OFF_BA + k, c);
        }
        for (int k = 0; k < 3; ++k) {
            out(OFF_PHI + k, c) = A(k,0) * ph[0] + A(k,1) * ph[1] + A(k,2) * ph[2];
            out(OFF_BG + k, c) = -(A(k,0) * bg[0] + A(k,1) * bg[1] + A(k,2) * bg[2]);
            out(OFF_BA + k, c) = -(A(k,0) * ba[0] + A(k,1) * ba[1] + A(k,2) * ba[2]);
        }
        for (int i = 0; i < NC; ++i)
            for (int k = 0; k < 3; ++k)
                out(world_index(i, k), c) = Sx[i](k,0) * ph[0] + Sx[i](k,1) * ph[1] + Sx[i](k,2) * ph[2]
                                          - at(world_index(i, k), c);
    }
}

TFGC3_TPL
void TFGC3_CLS::map_coordinates_(MatrixNX& P, MatrixNX& scratch, bool to_invariant) const {
    // P' = T P T^T = T (T P)^T for symmetric P: two row passes.
    apply_T_rows_(P, false, to_invariant, scratch);
    apply_T_rows_(scratch, true, to_invariant, P);
    for (int i = 0; i < NX; ++i)
        for (int j = i + 1; j < NX; ++j) {
            const T v = T(0.5) * (P(i,j) + P(j,i));
            P(i,j) = v;
            P(j,i) = v;
        }
}

TFGC3_TPL
Eigen::Matrix<T,6,6> TFGC3_CLS::physical_pv_covariance() const {
    // dp = -rho_p + [p]x phi, dv = -rho_v + [v]x phi.
    Eigen::Matrix<T,6,NX> J = Eigen::Matrix<T,6,NX>::Zero();
    J.template block<3,3>(0, OFF_PHI) = lie::skew<T>(Vector3(X_.X.col(IP)));
    J.template block<3,3>(0, OFF_P) = -Matrix3::Identity();
    J.template block<3,3>(3, OFF_PHI) = lie::skew<T>(Vector3(X_.X.col(IV)));
    J.template block<3,3>(3, OFF_V) = -Matrix3::Identity();
    Eigen::Matrix<T,6,6> C = J * P_ * J.transpose();
    return T(0.5) * (C + C.transpose());
}

TFGC3_TPL
void TFGC3_CLS::reset_aw_covariance_to_stationary() {
    MatC P1;
    if (!model_.stationary_covariance(P1)) return;
    begin_write_();
    for (int i = 0; i < NX; ++i) {
        if (i >= OFF_W && i < OFF_BA) continue;
        for (int j = OFF_W; j < OFF_BA; ++j) {
            P_(i, j) = T(0);
            P_(j, i) = T(0);
        }
    }
    for (int a = 0; a < NC; ++a)
        for (int b = 0; b < NC; ++b)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    P_(world_index(a, k), world_index(b, l)) = P1(a, b) * Sigma_s_(k, l);
    end_write_();
}

TFGC3_TPL
void TFGC3_CLS::set_initial_acc_bias_std(T s) {
    sigma_bacc0_ = std::max(T(0), s);
    // A diagonal body-frame block; Rhat D Rhat^T in the invariant frame
    // (Physical), D itself (Native).
    begin_write_();
    P_.template block<3,3>(OFF_BA, OFF_BA) = Matrix3::Identity() * (sigma_bacc0_ * sigma_bacc0_);
    end_write_();
}

TFGC3_TPL
void TFGC3_CLS::set_acc_bias_updates_enabled(bool en) {
    if (acc_bias_updates_enabled_ == en) return;
    begin_write_();
    if (!en) {
        for (int i = 0; i < NX; ++i) {
            if (i >= OFF_BA) continue;
            for (int k = 0; k < 3; ++k) {
                P_(i, OFF_BA + k) = T(0);
                P_(OFF_BA + k, i) = T(0);
            }
        }
    } else {
        const T target = sigma_bacc0_ * sigma_bacc0_;
        for (int i = 0; i < 3; ++i) P_(OFF_BA + i, OFF_BA + i) = std::max(P_(OFF_BA + i, OFF_BA + i), target);
    }
    end_write_();
    acc_bias_updates_enabled_ = en;
}

TFGC3_TPL
void TFGC3_CLS::seed_attitude_(MatrixNX& Pe, T tilt_sigma_rad, T yaw_sigma_rad) const {
    // CARMA3: body-frame block with yaw along u = R_wb e_z (the world vertical
    // in the body), decorrelated from every other state.  Native writes the
    // same physical statement in the world-referred attitude coordinates.
    const bool physical = write_mode_ == CovarianceWriteMode::Physical;
    Vector3 u = physical ? Vector3(X_.R.transpose() * Vector3::UnitZ()) : Vector3(Vector3::UnitZ());
    const T n = u.norm();
    Matrix3 Patt;
    if (!(n > T(1e-8)) || !u.allFinite()) {
        Patt = Matrix3::Identity() * (yaw_sigma_rad * yaw_sigma_rad);
    } else {
        u /= n;
        const Matrix3 Py = u * u.transpose();
        Patt = (tilt_sigma_rad * tilt_sigma_rad) * (Matrix3::Identity() - Py)
             + (yaw_sigma_rad * yaw_sigma_rad) * Py;
        Patt = T(0.5) * (Patt + Patt.transpose());
    }
    decorrelate_attitude_(Pe);
    Pe.template block<3,3>(OFF_PHI, OFF_PHI) = Patt;
}

TFGC3_TPL
void TFGC3_CLS::decorrelate_attitude_(MatrixNX& Pe) {
    // CARMA3's pattern: attitude against every other state, and the gyro bias
    // against the world block (zero_AL_cross_cov_once_ spans both rows).
    for (int i = 3; i < NX; ++i)
        for (int k = 0; k < 3; ++k) {
            Pe(i, k) = T(0);
            Pe(k, i) = T(0);
        }
    for (int i = OFF_W; i < OFF_BA; ++i)
        for (int k = OFF_BG; k < OFF_BG + 3; ++k) {
            Pe(i, k) = T(0);
            Pe(k, i) = T(0);
        }
}

TFGC3_TPL
void TFGC3_CLS::set_quaternion_boat(const Eigen::Quaternion<T>& q_bw) {
    project_gyro_bias_();
    Eigen::Quaternion<T> q = q_bw;
    if (!(q.norm() > T(1e-8))) return;
    q.normalize();
    // CARMA3 leaves its own covariance untouched: the same e covariance is
    // re-expressed at the new attitude (Physical), or xi's is kept (Native).
    begin_write_();
    X_.R = q.toRotationMatrix();
    reorthonormalize_();
    end_write_();
}

TFGC3_TPL
bool TFGC3_CLS::initialize_from_acc_mag(Vector3 const& acc, Vector3 const& mag) {
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
    begin_write_();
    Eigen::Quaternion<T> qwb(R_wb_new);
    qwb.normalize();
    X_.R = qwb.toRotationMatrix().transpose();
    reorthonormalize_();
    B_w_ = X_.R * mag;
    // CARMA3: attitude 5e-4 I, decorrelated; gyro-bias variance floored.
    decorrelate_attitude_(P_);
    P_.template block<3,3>(OFF_PHI, OFF_PHI) = Matrix3::Identity() * T(5e-4);
    for (int i = 0; i < 3; ++i) P_(OFF_BG + i, OFF_BG + i) = std::max(P_(OFF_BG + i, OFF_BG + i), T(1e-6));
    end_write_();
    return true;
}

TFGC3_TPL
Eigen::Quaternion<T> TFGC3_CLS::quaternion_from_acc(Vector3 const& acc) {
    // CARMA3's world-to-body quaternion from gravity (returned as such).
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

TFGC3_TPL
bool TFGC3_CLS::initialize_from_acc(Vector3 const& acc) {
    const T anorm = acc.norm();
    if (!acc.allFinite() || !(anorm >= T(1e-8))) return false;
    project_gyro_bias_();
    begin_write_();
    Eigen::Quaternion<T> qwb = quaternion_from_acc(acc / anorm);
    qwb.normalize();
    X_.R = qwb.toRotationMatrix().transpose();
    reorthonormalize_();
    seed_attitude_(P_, T(0.035), T(1.5708));
    end_write_();
    return true;
}

TFGC3_TPL
bool TFGC3_CLS::initialize_from_attitude(Eigen::Quaternion<T> const& q_bw, T tilt_sigma_rad, T yaw_sigma_rad) {
    if (!q_bw.coeffs().allFinite() || !(q_bw.norm() > T(1e-8))) return false;
    set_quaternion_boat(q_bw);
    const T tilt = (std::isfinite(tilt_sigma_rad) && tilt_sigma_rad > T(0)) ? tilt_sigma_rad : T(0.035);
    const T yaw = (std::isfinite(yaw_sigma_rad) && yaw_sigma_rad > T(0)) ? yaw_sigma_rad : T(1.5708);
    begin_write_();
    seed_attitude_(P_, tilt, yaw);
    end_write_();
    return true;
}

TFGC3_TPL
bool TFGC3_CLS::initialize_from_acc_preserve_yaw(Vector3 const& acc_body) {
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

TFGC3_TPL
void TFGC3_CLS::step_blocks_(const Vector3& w, T Ts, const Matrix3& R0, const MatC& Fc,
                             const Eigen::Matrix<T,3,NC>& W1, StepBlocks& b) const {
    (void)Fc;
    const Vector3 wh = w * Ts;
    const Matrix3 Rinc = lie::Exp<T>(wh);
    b.R1 = R0 * Rinc;
    b.E = b.R1 * R0.transpose();
    b.F = -(R0 * (Ts * lie::left_jacobian<T>(wh)) * R0.transpose());
    for (int i = 0; i < NC; ++i) b.Xs[i] = lie::skew<T>(Vector3(W1.col(i)));
    b.phi_b = acc_bias_updates_enabled_ ? std::exp(-Ts / std::max(T(1e-3), tau_bacc_)) : T(1);
}

TFGC3_TPL
void TFGC3_CLS::step_noise_(const Vector3& w, T Ts, const Matrix3& R0, const StepBlocks& b,
                            StepNoise& q) const {
    namespace d = tfg_carma_detail;
    const Matrix3 Jh = Ts * lie::left_jacobian<T>(Vector3(w * Ts));   // int_0^h Exp(w t) dt
    Matrix3 Ig = Matrix3::Zero();   // int R Qg R^T
    Matrix3 Ib = Matrix3::Zero();   // int C Qbg C^T,  C(s) = int_s^h R
    Matrix3 IC = Matrix3::Zero();   // int C
    const bool iso = filter_detail::is_isotropic3<T>(Q_gyro_);
    for (int k = 0; k < 5; ++k) {
        const T s = T(0.5) * Ts * (T(1) + T(d::kGL5Nodes[k]));
        const T wq = T(0.5) * Ts * T(d::kGL5Weights[k]);
        const Vector3 ws = w * s;
        const Matrix3 Cs = R0 * (Jh - s * lie::left_jacobian<T>(ws));
        if (!iso) {
            const Matrix3 Rs = R0 * lie::Exp<T>(ws);
            Ig.noalias() += wq * (Rs * Q_gyro_ * Rs.transpose());
        }
        Ib.noalias() += wq * (Cs * Q_bg_ * Cs.transpose());
        IC.noalias() += wq * Cs;
    }
    if (iso) Ig = Matrix3::Identity() * (Q_gyro_(0,0) * Ts);
    q.Qpp = Ig + Ib;
    q.Qpp = T(0.5) * (q.Qpp + q.Qpp.transpose()).eval();
    q.Qpb = -(IC * Q_bg_ * b.R1.transpose());
    q.Qbb = Ts * (b.R1 * Q_bg_ * b.R1.transpose());
    q.Qbb = T(0.5) * (q.Qbb + q.Qbb.transpose()).eval();
    if (acc_bias_updates_enabled_) {
        const T tau_b = std::max(T(1e-3), tau_bacc_);
        const T qd_scale = -T(0.5) * tau_b * std::expm1(-T(2) * Ts / tau_b);
        q.Qaa = (b.R1 * Q_bacc_ * b.R1.transpose()) * qd_scale;
        q.Qaa = T(0.5) * (q.Qaa + q.Qaa.transpose()).eval();
    } else {
        q.Qaa.setZero();
    }
}

TFGC3_TPL
void TFGC3_CLS::apply_phi_(const MatrixNX& in, bool transposed, const StepBlocks& b, const MatC& Fc,
                           MatrixNX& out) const {
    auto at = [&](int r, int c) { return transposed ? in(c, r) : in(r, c); };
    for (int c = 0; c < NX; ++c) {
        T pphi[3], pbg[3], fpb[3];
        for (int k = 0; k < 3; ++k) { pphi[k] = at(OFF_PHI + k, c); pbg[k] = at(OFF_BG + k, c); }
        for (int k = 0; k < 3; ++k) {
            fpb[k] = b.F(k,0) * pbg[0] + b.F(k,1) * pbg[1] + b.F(k,2) * pbg[2];
            out(OFF_PHI + k, c) = pphi[k] + fpb[k];
            out(OFF_BG + k, c) = b.E(k,0) * pbg[0] + b.E(k,1) * pbg[1] + b.E(k,2) * pbg[2];
        }
        for (int i = 0; i < NC; ++i) {
            const Matrix3& S = b.Xs[i];
            for (int k = 0; k < 3; ++k) {
                T acc = S(k,0) * fpb[0] + S(k,1) * fpb[1] + S(k,2) * fpb[2];
                for (int j = 0; j < NC; ++j) acc += Fc(i,j) * at(world_index(j, k), c);
                out(world_index(i, k), c) = acc;
            }
        }
        T pba[3];
        for (int k = 0; k < 3; ++k) pba[k] = at(OFF_BA + k, c);
        for (int k = 0; k < 3; ++k)
            out(OFF_BA + k, c) = b.phi_b * (b.E(k,0) * pba[0] + b.E(k,1) * pba[1] + b.E(k,2) * pba[2]);
    }
}

TFGC3_TPL
void TFGC3_CLS::add_noise_(MatrixNX& P, const StepBlocks& b, const StepNoise& q, const MatC& Qc) const {
    auto add3 = [&](int r, int c, const Matrix3& M) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) P(r + i, c + j) += M(i,j);
    };
    add3(OFF_PHI, OFF_PHI, q.Qpp);
    add3(OFF_PHI, OFF_BG, q.Qpb);
    add3(OFF_BG, OFF_PHI, Matrix3(q.Qpb.transpose()));
    add3(OFF_BG, OFF_BG, q.Qbb);
    Matrix3 GQp[NC], GQb[NC];
    for (int i = 0; i < NC; ++i) {
        GQp[i] = b.Xs[i] * q.Qpp;
        GQb[i] = b.Xs[i] * q.Qpb;
    }
    for (int i = 0; i < NC; ++i) {
        const int oi = OFF_W + 3 * i;
        add3(oi, OFF_PHI, GQp[i]);
        add3(OFF_PHI, oi, Matrix3(GQp[i].transpose()));
        add3(oi, OFF_BG, GQb[i]);
        add3(OFF_BG, oi, Matrix3(GQb[i].transpose()));
        for (int j = 0; j < NC; ++j) {
            const int oj = OFF_W + 3 * j;
            const Matrix3 M = GQp[i] * b.Xs[j].transpose() + Qc(i,j) * Sigma_s_;
            add3(oi, oj, M);
        }
    }
    add3(OFF_BA, OFF_BA, q.Qaa);
}

TFGC3_TPL
void TFGC3_CLS::build_transition(const Vector3& gyr, T Ts, MatrixNX& Phi) const {
    Phi.setIdentity();
    if (!(Ts > T(0)) || !std::isfinite(Ts)) return;
    MatC Fc, Qc;
    if (!carma::discretize_lti<T,NC>(model_.A(), model_.W(), Ts, Fc, Qc)) return;
    const Vector3 w = gyr - X_.B.col(0);
    const Eigen::Matrix<T,3,NC> W1 = X_.X * Fc.transpose();
    StepBlocks b;
    step_blocks_(w, Ts, X_.R, Fc, W1, b);
    MatrixNX I = MatrixNX::Identity();
    apply_phi_(I, false, b, Fc, Phi);
}

TFGC3_TPL
void TFGC3_CLS::build_process_noise(const Vector3& gyr, T Ts, MatrixNX& Qd) const {
    Qd.setZero();
    if (!(Ts > T(0)) || !std::isfinite(Ts)) return;
    MatC Fc, Qc;
    if (!carma::discretize_lti<T,NC>(model_.A(), model_.W(), Ts, Fc, Qc)) return;
    const Vector3 w = gyr - X_.B.col(0);
    const Eigen::Matrix<T,3,NC> W1 = X_.X * Fc.transpose();
    StepBlocks b;
    step_blocks_(w, Ts, X_.R, Fc, W1, b);
    StepNoise q;
    step_noise_(w, Ts, X_.R, b, q);
    add_noise_(Qd, b, q, Qc);
}

TFGC3_TPL
void TFGC3_CLS::time_update(Vector3 const& gyr, T Ts) {
    project_gyro_bias_();
    if (!(Ts > T(0)) || !std::isfinite(Ts)) return;

    const Vector3 w = gyr - X_.B.col(0);
    last_gyr_bias_corrected_ = w;
    const Matrix3 R0 = X_.R;

    refresh_discretization_(Ts);
    const MatC& Fc = disc_.Phi();
    const MatC& Qc = disc_.Q();

    // Mean: attitude by the bias-corrected rate, every world axis by the exact
    // CARMA transition, the accelerometer bias by its OU decay (held: frozen).
    StepBlocks b;
    const Eigen::Matrix<T,3,NC> W1 = X_.X * Fc.transpose();
    step_blocks_(w, Ts, R0, Fc, W1, b);
    X_.R = b.R1;
    reorthonormalize_();
    X_.X = W1;
    if (acc_bias_updates_enabled_) X_.B.col(1) *= b.phi_b;

    // P <- Phi P Phi^T + Qd, with the structured Phi (two passes of Phi * .).
    MatrixNX& M = M_scratch_;
    apply_phi_(P_, false, b, Fc, M);
    apply_phi_(M, true, b, Fc, P_);
    StepNoise q;
    step_noise_(w, Ts, R0, b, q);
    add_noise_(P_, b, q, Qc);
    symmetrize_P_();
}

TFGC3_TPL
void TFGC3_CLS::build_reset_jacobian(const Tangent& a, MatrixNX& J) const {
    J.setIdentity();
    if (reset_mode_ == ResetMode::None) return;
    const Vector3 phi = a.template segment<3>(OFF_PHI);
    const Matrix3 Jr = lie::left_jacobian<T>(Vector3(-phi));
    J.template block<3,3>(OFF_PHI, OFF_PHI) = Jr;
    if (reset_mode_ == ResetMode::AttitudeOnly) return;
    for (int blk = 1; blk < NBLK; ++blk) {
        const int o = 3 * blk;
        J.template block<3,3>(o, o) = Jr;
        J.template block<3,3>(o, OFF_PHI) =
            tfg_carma_detail::right_jacobian_cross<T>(phi, Vector3(a.template segment<3>(o)));
    }
}

TFGC3_TPL
void TFGC3_CLS::transport_(const Tangent& a) {
    if (reset_mode_ == ResetMode::None) return;
    // P <- J P J^T with J = blockdiag(Jr) + C (C only in the attitude column):
    // (J P)_k = Jr P_k + C_k P_phi;  P'_{.,k} = (JP)_{.,k} Jr^T + (JP)_{.,phi} C_k^T.
    const Vector3 phi = a.template segment<3>(OFF_PHI);
    const Matrix3 Jr = lie::left_jacobian<T>(Vector3(-phi));
    Matrix3 C[NBLK];
    C[0].setZero();
    const bool full = reset_mode_ == ResetMode::Full;
    for (int blk = 1; blk < NBLK; ++blk)
        C[blk] = full ? tfg_carma_detail::right_jacobian_cross<T>(phi, Vector3(a.template segment<3>(3 * blk)))
                      : Matrix3::Zero();
    auto diag_of = [&](int blk) -> Matrix3 {
        return (full || blk == 0) ? Jr : Matrix3::Identity();
    };
    MatrixNX& M = M_scratch_;
    for (int blk = 0; blk < NBLK; ++blk) {
        const Matrix3 D = diag_of(blk);
        const int o = 3 * blk;
        for (int c = 0; c < NX; ++c)
            for (int k = 0; k < 3; ++k) {
                T acc = D(k,0) * P_(o, c) + D(k,1) * P_(o + 1, c) + D(k,2) * P_(o + 2, c);
                if (blk != 0)
                    acc += C[blk](k,0) * P_(0, c) + C[blk](k,1) * P_(1, c) + C[blk](k,2) * P_(2, c);
                M(o + k, c) = acc;
            }
    }
    for (int blk = 0; blk < NBLK; ++blk) {
        const Matrix3 D = diag_of(blk);
        const int o = 3 * blk;
        for (int r = 0; r < NX; ++r)
            for (int k = 0; k < 3; ++k) {
                T acc = M(r, o) * D(k,0) + M(r, o + 1) * D(k,1) + M(r, o + 2) * D(k,2);
                if (blk != 0)
                    acc += M(r, 0) * C[blk](k,0) + M(r, 1) * C[blk](k,1) + M(r, 2) * C[blk](k,2);
                P_(r, o + k) = acc;
            }
    }
    symmetrize_P_();
}

TFGC3_TPL
void TFGC3_CLS::inject(const Tangent& a) {
    X_ = Group::Exp_G(Tangent(-a)).compose(X_);
    reorthonormalize_();
}

TFGC3_TPL
bool TFGC3_CLS::correct3_(const Vector3& r, Matrix3& S, MatrixNX3& PHt, T noise_scale, bool freeze_ba,
                          MeasDiag3& diag) {
    Eigen::LDLT<Matrix3> ldlt;
    if (!safe_ldlt3_(S, ldlt, noise_scale)) {
        diag.S = S;
        return false;
    }
    diag.S = S;
    diag.nis = nis3_(ldlt, r);
    MatrixNX3& K = K_scratch_;
    gain3_(PHt, ldlt, K);
    if (freeze_ba) freeze_acc_bias_rows_(K);
    Tangent a;
    for (int i = 0; i < NX; ++i) a(i) = K(i,0) * r(0) + K(i,1) * r(1) + K(i,2) * r(2);
    if (!a.allFinite()) {
        ++skipped_updates_;
        return false;
    }
    joseph_update_(K, S, PHt);
    transport_(a);
    inject(a);
    project_gyro_bias_();
    project_acc_bias_();
    diag.accepted = true;
    return true;
}

TFGC3_TPL
void TFGC3_CLS::measurement_update_acc_only(Vector3 const& acc_meas, T tempC) {
    last_acc_diag_ = MeasDiag3{};
    if (!acc_meas.allFinite() || !std::isfinite(tempC)) return;
    const bool use_ba = acc_bias_updates_enabled_;

    const Vector3 g_world(0, 0, gravity_magnitude_);
    const Vector3 ba_term = X_.B.col(1) + k_a_ * (tempC - tempC_ref);
    const Matrix3& R = X_.R;
    const Vector3 r = R * (acc_meas - ba_term) - (Vector3(X_.X.col(IA)) - g_world);
    last_acc_diag_.r = r;
    const Matrix3 Rw = R * Racc_ * R.transpose();

    // H = [ [g]x | 0 | 0 0 -I 0.. | -I ];  [g]x = [[0,-g,0],[g,0,0],[0,0,0]].
    const T g = gravity_magnitude_;
    MatrixNX3& PHt = PHt_scratch_;
    for (int i = 0; i < NX; ++i) {
        PHt(i,0) = -g * P_(i, OFF_PHI + 1) - P_(i, OFF_A + 0) - P_(i, OFF_BA + 0);
        PHt(i,1) =  g * P_(i, OFF_PHI + 0) - P_(i, OFF_A + 1) - P_(i, OFF_BA + 1);
        PHt(i,2) = -P_(i, OFF_A + 2) - P_(i, OFF_BA + 2);
    }
    Matrix3 S = Rw;
    for (int b = 0; b < 3; ++b) {
        S(0,b) += -g * PHt(OFF_PHI + 1, b) - PHt(OFF_A + 0, b) - PHt(OFF_BA + 0, b);
        S(1,b) +=  g * PHt(OFF_PHI + 0, b) - PHt(OFF_A + 1, b) - PHt(OFF_BA + 1, b);
        S(2,b) += -PHt(OFF_A + 2, b) - PHt(OFF_BA + 2, b);
    }
    for (int a = 0; a < 3; ++a)
        for (int b = a + 1; b < 3; ++b) {
            const T v = T(0.5) * (S(a,b) + S(b,a));
            S(a,b) = v;
            S(b,a) = v;
        }
    if (!use_ba) {
        // CARMA3's held bias: in S, not in the cross-covariance or the gain.
        for (int i = 0; i < NX; ++i)
            for (int k = 0; k < 3; ++k) PHt(i,k) += P_(i, OFF_BA + k);
        freeze_acc_bias_rows_(PHt);
    }
    correct3_(r, S, PHt, Racc_.norm(), !use_ba, last_acc_diag_);
}

TFGC3_TPL
void TFGC3_CLS::measurement_update_mag_only(const Vector3& mag_meas) {
    last_mag_diag_ = MeasDiag3{};
    if (!mag_meas.allFinite()) return;
    if (!(mag_meas.norm() > T(1e-6))) return;
    const Matrix3& R = X_.R;
    const Vector3 r = R * mag_meas - B_w_;
    last_mag_diag_.r = r;
    const Matrix3 Hphi = -lie::skew<T>(B_w_);
    const Matrix3 Rw = R * Rmag_ * R.transpose();
    MatrixNX3& PHt = PHt_scratch_;
    for (int i = 0; i < NX; ++i)
        for (int b = 0; b < 3; ++b)
            PHt(i,b) = P_(i, OFF_PHI) * Hphi(b,0) + P_(i, OFF_PHI + 1) * Hphi(b,1) + P_(i, OFF_PHI + 2) * Hphi(b,2);
    Matrix3 S = Rw;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            S(a,b) += Hphi(a,0) * PHt(OFF_PHI, b) + Hphi(a,1) * PHt(OFF_PHI + 1, b) + Hphi(a,2) * PHt(OFF_PHI + 2, b);
    const bool freeze = !acc_bias_updates_enabled_;
    if (freeze) freeze_acc_bias_rows_(PHt);
    correct3_(r, S, PHt, Rmag_.norm(), freeze, last_mag_diag_);
}

#undef TFGC3_TPL
#undef TFGC3_CLS

}  // namespace ocean_imu::kalman
