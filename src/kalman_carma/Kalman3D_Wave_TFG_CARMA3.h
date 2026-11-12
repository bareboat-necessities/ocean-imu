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
  Jacobians are state independent (without the optional lever arm below):

      r_a = Rhat (f_m - bhat_a - k_a dT) - (ahat - g),   H_a = [ [g]x 0 .. -I(rho_a) .. -I(beta_a) ]
      r_m = Rhat m - B_w,                               H_m = [ -[B_w]x 0 .. ]

  with noise Rhat R Rhat^T.  A finite correction a = K r is injected as
  Exp_G(-a) Xhat and the covariance transported by the full group right
  Jacobian J_r(a) = d Log_G(Exp_G(-a) Exp_G(a + d)) / d d.

  Optional IMU lever arm (set_imu_lever_arm_body, off by default): the
  rotational term L = alpha x r + omega x (omega x r) of an IMU at body
  offset r is removed from f_m with CARMA3's model,
  kalman_common/ImuLeverArm.h.  L depends on the gyro bias through the
  bias-corrected rate, d L / d b_g = J_L, and b_g - bhat_g = -Rhat^T beta_g to
  first order, so the accelerometer Jacobian gains

      H_a,beta_g = -Rhat J_L Rhat^T,

  the one estimate-dependent block, present only while a lever arm is set.

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
#include <type_traits>

#include "kalman_common/ImuLeverArm.h"
#include "kalman_common/MekfAttitudeMath.h"
#include "kalman_common/ConditionalCovariance.h"
#include "kalman_carma/CarmaVerticalModel.h"
#include "lie/SO3Jacobians.h"
#include "lie/TwoFrameGroup.h"

namespace ocean_imu::kalman {

namespace tfg_carma_detail {
// Empty policy storage; the default estimator retains its original model.
struct EmptyMotion {
    explicit EmptyMotion(const carma::CarmaShape& = carma::CARMA41_JONSWAP_SHAPE) {}
};


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
/*
    Q_r(phi, .) is linear in z and its matrix products collapse: with
    P = [-phi]x, Z = [-z]x, s = phi . z, t2 = |phi|^2 and
    [a]x [b]x = b a^T - (a . b) I,

        P Z + Z P = z phi^T + phi z^T - 2 s I,   P Z P = -s P,
        P P Z + Z P P = 2 s P - t2 Z,            P Z P P + P P Z P = -2 s P^2,

    so

        Q_r = -(1/2 - c2 t2) [z]x + c1 (z phi^T + phi z^T) + s N,
        N   = -2 c1 I - (2 c2 - c1) [phi]x - 2 c3 (phi phi^T - t2 I).

    The coefficients and N depend on phi alone: one reset evaluates them once
    for all of its vector blocks.  Plain arrays, as the sample path below.
*/
template<typename T>
class RightJacobianCross {
public:
    explicit RightJacobianCross(const Eigen::Matrix<T,3,1>& phi) {
        const T x = phi(0), y = phi(1), z = phi(2);
        const T t2 = x * x + y * y + z * z;
        T c1, c2, c3;
        coefficients_(t2, c1, c2, c3);
        phi_[0] = x; phi_[1] = y; phi_[2] = z;
        k0_ = T(0.5) - c2 * t2;
        k1_ = c1;
        const T k2 = T(2) * c2 - c1, k3 = T(2) * c3;
        const T d = -T(2) * c1 + k3 * t2;
        n_[0][0] = d - k3 * x * x;       n_[0][1] = k2 * z - k3 * x * y;  n_[0][2] = -k2 * y - k3 * x * z;
        n_[1][0] = -k2 * z - k3 * x * y; n_[1][1] = d - k3 * y * y;       n_[1][2] = k2 * x - k3 * y * z;
        n_[2][0] = k2 * y - k3 * x * z;  n_[2][1] = -k2 * x - k3 * y * z; n_[2][2] = d - k3 * z * z;
    }
    // out = Q_r(phi, z), row-major.
    void operator()(const T* z, T (&out)[3][3]) const {
        const T s = phi_[0] * z[0] + phi_[1] * z[1] + phi_[2] * z[2];
        const T y[3] = {k1_ * z[0], k1_ * z[1], k1_ * z[2]};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) out[i][j] = y[i] * phi_[j] + phi_[i] * y[j] + s * n_[i][j];
        const T h[3] = {k0_ * z[0], k0_ * z[1], k0_ * z[2]};
        out[0][1] += h[2]; out[0][2] -= h[1];
        out[1][0] -= h[2]; out[1][2] += h[0];
        out[2][0] += h[1]; out[2][1] -= h[0];
    }

private:
    static void coefficients_(T t2, T& c1, T& c2, T& c3);
    T phi_[3], n_[3][3], k0_, k1_;
};

template<typename T>
void RightJacobianCross<T>::coefficients_(T t2, T& c1, T& c2, T& c3) {
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
}

template<typename T>
inline Eigen::Matrix<T,3,3> right_jacobian_cross(const Eigen::Matrix<T,3,1>& phi,
                                                 const Eigen::Matrix<T,3,1>& z) {
    T q[3][3];
    RightJacobianCross<T>{phi}(z.data(), q);
    Eigen::Matrix<T,3,3> Q;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) Q(i,j) = q[i][j];
    return Q;
}

// Row-major plain copy of a 3 x 3 matrix.  The sample-path covariance kernels
// index plain arrays: a firmware build (-Os, Eigen assertions on) keeps
// Eigen's checked coefficient accessors and small products out of line.
template<typename T>
inline void rows3(const Eigen::Matrix<T,3,3>& M, T (&out)[3][3]) {
    const T* m = M.data();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) out[i][j] = m[i + 3 * j];
}

// Coefficients of J_l(phi) = I + a [phi]x + b [phi]x^2 from t2 = |phi|^2:
// lie::left_jacobian()'s own expressions and small-angle branch.
template<typename T>
inline void left_jacobian_coefficients(T t2, T& a, T& b) {
    const T t = std::sqrt(t2);
    if (t < lie::kSmallAngle<T>) {
        const T t4 = t2 * t2;
        a = T(0.5);
        a = std::fma(-t2, T(1)/T(24), a);
        a = std::fma( t4, T(1)/T(720), a);
        b = T(1)/T(6);
        b = std::fma(-t2, T(1)/T(120), b);
        b = std::fma( t4, T(1)/T(5040), b);
    } else {
        a = (T(1) - std::cos(t)) / t2;
        b = (t - std::sin(t)) / (t2 * t);
    }
}

// Coefficients of Exp(phi) = I + a [phi]x + b [phi]x^2, as lie::Exp().
template<typename T>
inline void exp_coefficients(T t2, T& a, T& b) {
    const T t = std::sqrt(t2);
    if (t < lie::kSmallAngle<T>) {
        const T t4 = t2 * t2;
        a = T(1);
        a = std::fma(-t2, T(1)/T(6), a);
        a = std::fma( t4, T(1)/T(120), a);
        b = T(0.5);
        b = std::fma(-t2, T(1)/T(24), b);
        b = std::fma( t4, T(1)/T(720), b);
    } else {
        a = std::sin(t) / t;
        b = (T(1) - std::cos(t)) / t2;
    }
}

// out = I + a [v]x + b [v]x^2, with [v]x^2 = v v^T - |v|^2 I.
template<typename T>
inline void so3_series(const T (&v)[3], T a, T b, T (&out)[3][3]) {
    const T n2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) out[i][j] = b * (v[i] * v[j]);
    const T d = T(1) - b * n2;
    out[0][0] += d; out[1][1] += d; out[2][2] += d;
    out[0][1] -= a * v[2]; out[0][2] += a * v[1];
    out[1][0] += a * v[2]; out[1][2] -= a * v[0];
    out[2][0] -= a * v[1]; out[2][1] += a * v[0];
}

// out = A B, or A B^T when `transposed`.
template<typename T>
inline void mul3(const T (&A)[3][3], const T (&B)[3][3], bool transposed, T (&out)[3][3]) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            out[i][j] = transposed ? A[i][0] * B[j][0] + A[i][1] * B[j][1] + A[i][2] * B[j][2]
                                   : A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
}

// Five-point Gauss-Legendre nodes and weights on [-1, 1].
inline constexpr double kGL5Nodes[5] = {
    -0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640};
inline constexpr double kGL5Weights[5] = {
    0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665,
    0.2369268850561891};

}  // namespace tfg_carma_detail

// NC counts world columns, not necessarily the CARMA autoregressive order.
// Motion=void retains the original independent NC-state realization. A motion
// policy provides its exact common dynamics and Hermitian covariance pair;
// SharedModesMotion has six columns but still models CARMA(4,1).
template <typename T = float, int NC = 4, class Motion = void>
class Kalman3D_Wave_TFG_CARMA3 {
    static constexpr bool kMotion = !std::is_void_v<Motion>;
    static_assert((!kMotion && NC >= 3 && NC <= 5) || (kMotion && NC == 6),
                  "Independent CARMA order 3..5, or six-state shared CARMA(4,1)");
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
    using Model = carma::CarmaVerticalModel<T, kMotion ? 4 : NC>;

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
                             const carma::CarmaShape& shape = carma::CarmaDefaultShape<Model::ORDER>::value);

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
        if constexpr (kMotion) return motion_.accel_std();
        else return Sigma_s_.diagonal().cwiseMax(T(0)).cwiseSqrt();
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
    bool set_carma_shape(const carma::CarmaShape& shape) { if constexpr (kMotion) return motion_.set_shape(shape);
        else return model_.set_shape(shape); }
    bool set_horizontal_shape(const carma::CarmaShape& shape) { if constexpr (kMotion) return motion_.set_shape(shape);
        else return same_shape_(shape); }
    bool set_vertical_shape(const carma::CarmaShape& shape) { if constexpr (kMotion) return motion_.set_shape(shape);
        else return same_shape_(shape); }
    [[nodiscard]] const Model& vertical_model() const noexcept { if constexpr (kMotion) return motion_.vertical_model(); else return model_; }
    [[nodiscard]] const Model& horizontal_model() const noexcept { if constexpr (kMotion) return motion_.horizontal_model(); else return model_; }
    [[nodiscard]] Eigen::Matrix<T,2,2> horizontal_sigma() const {
        if constexpr (kMotion) return motion_.horizontal_sigma();
        else return Sigma_s_.template topLeftCorner<2,2>();
    }
    [[nodiscard]] const Matrix3& spatial_sigma() const noexcept { if constexpr (kMotion) return motion_.excitation().re; else return Sigma_s_; }
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

    // IMU installation lever arm (off-CoG correction), CARMA3's API and model
    // (kalman_common/ImuLeverArm.h).  r_b is the IMU position relative to the
    // CoG in the body frame (x forward, y starboard, z down), metres.  Zero,
    // the default, switches the correction off and leaves the filter exactly
    // as it is without it.
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

    // Shared motion policy; not available on the independent estimator.
    auto& motion() noexcept requires(kMotion) { return motion_; }
    const auto& motion() const noexcept requires(kMotion) { return motion_; }

    // Conditional Schur shrink is independent of the error representation:
    // conditioning on attitude and biases removes rho's attitude shear.
    bool shrink_motion_vertical_coupling(T alpha) requires(kMotion) {
        int r[9], h[2*NC], z[NC];
        for (int i = 0; i < 6; ++i) r[i] = i;
        for (int i = 0; i < 3; ++i) r[6+i] = OFF_BA+i;
        for (int i = 0; i < NC; ++i) {
            h[i] = world_index(i,0); h[NC+i] = world_index(i,1); z[i] = world_index(i,2);
        }
        return filter_detail::shrink_conditional_cross<T>(P_, r, h, z, alpha);
    }

    // x_axis <- M_axis x_axis, P transported through the physical additive
    // coordinates. This includes [w'_i]x - M_axis [w_i]x in the invariant
    // attitude column when the axis maps do not commute with rotations.
    // Only 6-vectors are on the stack; the covariance scratch is a member.
    void apply_motion_recoordinatization(const MatC (&maps)[3]) requires(kMotion) {
        map_coordinates_(P_, M_scratch_, false);
        T* const p = P_.data();
        T* const xw = X_.X.data();                 // X(ax, i) = xw[ax + 3 i]
        for (int ax = 0; ax < 3; ++ax) {
            const T* const mp = maps[ax].data();   // map(i, j) = mp[i + NC j]
            T v[NC], y[NC];
            auto apply = [&] {
                for (int i = 0; i < NC; ++i) {
                    T s = T(0);
                    for (int j = 0; j < NC; ++j) s += mp[i + NC * j] * v[j];
                    y[i] = s;
                }
            };
            for (int i = 0; i < NC; ++i) v[i] = xw[ax + 3 * i];
            apply();
            for (int i = 0; i < NC; ++i) xw[ax + 3 * i] = y[i];
            for (int c = 0; c < NX; ++c) {
                for (int i = 0; i < NC; ++i) v[i] = p[world_index(i,ax) + NX * c];
                apply();
                for (int i = 0; i < NC; ++i) p[world_index(i,ax) + NX * c] = y[i];
            }
            for (int r = 0; r < NX; ++r) {
                for (int i = 0; i < NC; ++i) v[i] = p[r + NX * world_index(i,ax)];
                apply();
                for (int i = 0; i < NC; ++i) p[r + NX * world_index(i,ax)] = y[i];
            }
        }
        map_coordinates_(P_, M_scratch_, true);
    }

private:
    const T gravity_magnitude_ = T(STD_GRAVITY);

    Group X_{};
    MatrixNX P_;
    Vector3 B_w_ = Vector3::UnitX();
    Vector3 last_gyr_bias_corrected_ = Vector3::Zero();
    filter_detail::ImuLeverArm<T> lever_;

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
    [[no_unique_address]] std::conditional_t<kMotion, tfg_carma_detail::EmptyMotion, Model> model_;
    [[no_unique_address]] std::conditional_t<kMotion, Motion, tfg_carma_detail::EmptyMotion> motion_;
    Matrix3 Sigma_s_ = Matrix3::Identity();
    [[no_unique_address]] std::conditional_t<kMotion, tfg_carma_detail::EmptyMotion, carma::StepDiscretization<T,NC>> disc_;
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
        // Not "F": the Arduino core defines a function-like F() macro.
        Matrix3 Phi_pb;         // Phi_{phi,bg}
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
    if constexpr (kMotion) {
        motion_.set_heading(std::atan2(X_.R(1,0), X_.R(0,0)));
        motion_.refresh(Ts);
        recoordinate_motion_();
    } else {
        if (disc_.advance(model_.generation(), Ts)) return;
        disc_.refresh(model_.A(), model_.W(), model_.generation(), Ts);

    }
}

    void recoordinate_motion_() {
        if constexpr (kMotion) {
            MatC maps[3];
            if (motion_.take_axis_recoordinatization(maps)) apply_motion_recoordinatization(maps);
        }
    }
    const MatC& motion_F_() const {
        if constexpr (kMotion) return motion_.block_F(); else return disc_.Phi();
    }
    const MatC& motion_Q_() const {
        if constexpr (kMotion) return motion_.block_QR(); else return disc_.Q();
    }
    const MatC* motion_QI_() const {
        if constexpr (kMotion) return &motion_.block_QI(); else return nullptr;
    }
    bool discretize_motion_(T h, MatC& F, MatC& QR, MatC& QI) const {
        if constexpr (kMotion) return motion_.discretize(h, F, QR, QI);
        else return carma::discretize_lti<T,NC>(model_.A(), model_.W(), h, F, QR);
    }

    void step_blocks_(const Vector3& w, T Ts, const Matrix3& R0, const MatC& Fc,
                      const Eigen::Matrix<T,3,NC>& W1, StepBlocks& b) const;
    void step_noise_(const Vector3& w, T Ts, const Matrix3& R0, const StepBlocks& b, StepNoise& q) const;
    // out = Phi * in (in read transposed when `transposed`).
    void apply_phi_(const MatrixNX& in, bool transposed, const StepBlocks& b, const MatC& Fc,
                    MatrixNX& out) const;
    void add_noise_(MatrixNX& P, const StepBlocks& b, const StepNoise& q, const MatC& Qc, const MatC* Qi = nullptr) const;
    // P_ <- Phi P_ Phi^T + Qd, the sample path's form of apply_phi_() twice
    // and add_noise_(): M = Phi P, then only the upper block triangle of
    // M Phi^T + Qd, mirrored.
    void propagate_covariance_(const StepBlocks& b, const StepNoise& q, const MatC& Fc, const MatC& Qc,
                               const MatC* Qi);

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
        T Si[3][3];
        tfg_carma_detail::rows3<T>(Sinv, Si);
        const T* const h = PHt.data();
        T* const k = K.data();
        for (int b = 0; b < 3; ++b)
            for (int i = 0; i < NX; ++i)
                k[i + NX * b] = h[i] * Si[0][b] + h[i + NX] * Si[1][b] + h[i + 2 * NX] * Si[2][b];
    }
    // CARMA3's Joseph form P + K S K^T - K PHt^T - PHt K^T, exactly symmetric.
    void joseph_update_(const MatrixNX3& K, const Matrix3& S, const MatrixNX3& PHt) {
        T Kl[NX][3], G[NX][3], Sl[3][3];
        tfg_carma_detail::rows3<T>(S, Sl);
        const T* const k = K.data();
        const T* const h = PHt.data();
        for (int i = 0; i < NX; ++i)
            for (int b = 0; b < 3; ++b) Kl[i][b] = k[i + NX * b];
        for (int i = 0; i < NX; ++i)
            for (int b = 0; b < 3; ++b) {
                T g = -h[i + NX * b];
                for (int a = 0; a < 3; ++a) g += Kl[i][a] * Sl[a][b];
                G[i][b] = g;
            }
        T* const p = P_.data();
        for (int i = 0; i < NX; ++i) {
            const T k0 = Kl[i][0], k1 = Kl[i][1], k2 = Kl[i][2];
            const T h0 = h[i], h1 = h[i + NX], h2 = h[i + 2 * NX];
            for (int j = i; j < NX; ++j) {
                T d = k0 * G[j][0] - h0 * Kl[j][0];
                d += k1 * G[j][1] - h1 * Kl[j][1];
                d += k2 * G[j][2] - h2 * Kl[j][2];
                const T v = p[i + NX * j] + d;
                p[i + NX * j] = v;
                p[j + NX * i] = v;
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

#define TFGC3_TPL template <typename T, int NC, class Motion>
#define TFGC3_CLS Kalman3D_Wave_TFG_CARMA3<T,NC,Motion>

TFGC3_TPL
TFGC3_CLS::Kalman3D_Wave_TFG_CARMA3(
    Vector3 const& sigma_a, Vector3 const& gyro_noise_density_rad_sqrt_s, Vector3 const& sigma_m,
    T Pq0, T Pb0, T b0, T R_S_noise_var, T gravity_magnitude, const carma::CarmaShape& shape)
  : gravity_magnitude_(gravity_magnitude),
    Q_gyro_(gyro_noise_density_rad_sqrt_s.array().square().matrix().asDiagonal()),
    Q_bg_(Matrix3::Identity() * b0),
    Racc_(sigma_a.array().square().matrix().asDiagonal()),
    Rmag_(sigma_m.array().square().matrix().asDiagonal()),
    R_S_(Matrix3::Identity() * R_S_noise_var),
    model_(shape), motion_(shape)
{
    // CARMA3's constructor covariance; at the identity state the coordinate
    // map only flips signs, so it is written directly.
    P_.setZero();
    P_.template block<3,3>(OFF_PHI, OFF_PHI) = Matrix3::Identity() * Pq0;
    P_.template block<3,3>(OFF_BG, OFF_BG) = Matrix3::Identity() * Pb0;
    P_.template block<3,3>(OFF_BA, OFF_BA) = Matrix3::Identity() * (sigma_bacc0_ * sigma_bacc0_);
    if constexpr (kMotion) {
        motion_.set_horizontal(T(1), T(1), T(1), T(0));
        motion_.set_vertical(T(1), T(1));
    } else model_.configure(T(1), T(1));
    Sigma_s_.setIdentity();
    reset_aw_covariance_to_stationary();
}

TFGC3_TPL
bool TFGC3_CLS::set_vertical_carma(T kappa, T accel_std) {
    if constexpr (kMotion) {
        return motion_.set_vertical(kappa, accel_std);
    } else {
    if (!(accel_std > T(0)) || !std::isfinite(accel_std)) return false;
    if (!model_.configure(kappa, T(1))) return false;
    Sigma_s_(2,2) = accel_std * accel_std;
    return true;

    }
}

TFGC3_TPL
bool TFGC3_CLS::set_horizontal_carma(T kappa, T sx, T sy, T r) {
    if constexpr (kMotion) {
        return motion_.set_horizontal(kappa, sx, sy, r);
    } else {
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
}

TFGC3_TPL
void TFGC3_CLS::apply_T_rows_(const MatrixNX& in, bool transposed, bool to_invariant, MatrixNX& out) const {
    // T   = [phi <- dtheta: R], [bg: -R], [w_i: -I, w_i <- dtheta: [w_i]x R], [ba: -R];
    // T^-1 = [dtheta <- phi: R^T], [bg: -R^T], [w_i: -I, w_i <- phi: [w_i]x], [ba: -R^T].
    T A[3][3], Sx[NC][3][3];
    tfg_carma_detail::rows3<T>(to_invariant ? Matrix3(X_.R) : Matrix3(X_.R.transpose()), A);
    for (int i = 0; i < NC; ++i) {
        const Matrix3 S = lie::skew<T>(Vector3(X_.X.col(i)));
        tfg_carma_detail::rows3<T>(to_invariant ? Matrix3(S * X_.R) : S, Sx[i]);
    }
    const T* const ip = in.data();
    T* const op = out.data();
    auto at = [&](int r, int c) { return transposed ? ip[c + NX * r] : ip[r + NX * c]; };
    for (int c = 0; c < NX; ++c) {
        T* const oc = op + NX * c;
        T ph[3], bg[3], ba[3];
        for (int k = 0; k < 3; ++k) {
            ph[k] = at(OFF_PHI + k, c);
            bg[k] = at(OFF_BG + k, c);
            ba[k] = at(OFF_BA + k, c);
        }
        for (int k = 0; k < 3; ++k) {
            oc[OFF_PHI + k] = A[k][0] * ph[0] + A[k][1] * ph[1] + A[k][2] * ph[2];
            oc[OFF_BG + k] = -(A[k][0] * bg[0] + A[k][1] * bg[1] + A[k][2] * bg[2]);
            oc[OFF_BA + k] = -(A[k][0] * ba[0] + A[k][1] * ba[1] + A[k][2] * ba[2]);
        }
        for (int i = 0; i < NC; ++i)
            for (int k = 0; k < 3; ++k)
                oc[world_index(i, k)] = Sx[i][k][0] * ph[0] + Sx[i][k][1] * ph[1] + Sx[i][k][2] * ph[2]
                                      - at(world_index(i, k), c);
    }
}

TFGC3_TPL
void TFGC3_CLS::map_coordinates_(MatrixNX& P, MatrixNX& scratch, bool to_invariant) const {
    // P' = T P T^T = T (T P)^T for symmetric P: two row passes.
    apply_T_rows_(P, false, to_invariant, scratch);
    apply_T_rows_(scratch, true, to_invariant, P);
    T* const p = P.data();
    for (int i = 0; i < NX; ++i)
        for (int j = i + 1; j < NX; ++j) {
            const T v = T(0.5) * (p[i + NX * j] + p[j + NX * i]);
            p[i + NX * j] = v;
            p[j + NX * i] = v;
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
    if constexpr (kMotion) {
        motion_.refresh(motion_.step());
        recoordinate_motion_();
        MatC PR, Pimag;
        if (!motion_.stationary_pair(PR, Pimag)) return;
        begin_write_();
        for (int i = 0; i < NX; ++i) {
            if (i >= OFF_W && i < OFF_BA) continue;
            for (int j = OFF_W; j < OFF_BA; ++j) P_(i,j) = P_(j,i) = T(0);
        }
        for (int i = 0; i < NC; ++i)
            for (int j = 0; j < NC; ++j)
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b)
                        P_(world_index(i,a),world_index(j,b)) =
                            PR(i,j) * motion_.excitation().re(a,b) + Pimag(i,j) * motion_.excitation().im(a,b);
        end_write_();
    } else {
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
    // As in CARMA3: the wrapper hands over here after holding the prediction.
    lever_.reset_kinematics();
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
    namespace d = tfg_carma_detail;
    const T wh[3] = {w(0) * Ts, w(1) * Ts, w(2) * Ts};
    const T t2 = wh[0] * wh[0] + wh[1] * wh[1] + wh[2] * wh[2];
    T ea, eb, ja, jb, R0a[3][3], Rinc[3][3], Jh[3][3], R1[3][3], X[3][3], Y[3][3];
    d::exp_coefficients<T>(t2, ea, eb);
    d::left_jacobian_coefficients<T>(t2, ja, jb);
    d::rows3<T>(R0, R0a);
    d::so3_series<T>(wh, ea, eb, Rinc);            // Exp(w h)
    d::so3_series<T>(wh, ja, jb, Jh);              // J_l(w h)
    d::mul3<T>(R0a, Rinc, false, R1);
    d::mul3<T>(R1, R0a, true, X);                  // E = R1 R0^T
    d::mul3<T>(R0a, Jh, false, Y);
    d::mul3<T>(Y, R0a, true, Jh);                  // R0 J_l R0^T
    T* const r1 = b.R1.data();
    T* const e = b.E.data();
    T* const pb = b.Phi_pb.data();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            r1[i + 3 * j] = R1[i][j];
            e[i + 3 * j] = X[i][j];
            pb[i + 3 * j] = -(Ts * Jh[i][j]);
        }
    const T* const w1 = W1.data();
    for (int i = 0; i < NC; ++i) {
        const T x = w1[3 * i], y = w1[3 * i + 1], z = w1[3 * i + 2];
        T* const xs = b.Xs[i].data();              // [w_i]x, column-major
        xs[0] = T(0); xs[1] = z;    xs[2] = -y;
        xs[3] = -z;   xs[4] = T(0); xs[5] = x;
        xs[6] = y;    xs[7] = -x;   xs[8] = T(0);
    }
    b.phi_b = acc_bias_updates_enabled_ ? std::exp(-Ts / std::max(T(1e-3), tau_bacc_)) : T(1);
}

TFGC3_TPL
void TFGC3_CLS::step_noise_(const Vector3& w, T Ts, const Matrix3& R0, const StepBlocks& b,
                            StepNoise& q) const {
    namespace d = tfg_carma_detail;
    // With W = [w]x, s J_l(w s) = int_0^s Exp(w t) dt = s I + a s^2 W + b s^3 W^2,
    // (a, b) the J_l coefficients at |w| s.  The rotation R0 factors out of
    // every integral: C(s) = int_s^h R = R0 K(s), R(s) = R0 Exp(w s).
    T R0a[3][3], R1a[3][3], Qg[3][3], Qb[3][3];
    d::rows3<T>(R0, R0a);
    d::rows3<T>(b.R1, R1a);
    d::rows3<T>(Q_gyro_, Qg);
    d::rows3<T>(Q_bg_, Qb);
    const T wv[3] = {w(0), w(1), w(2)};
    const T n2 = wv[0] * wv[0] + wv[1] * wv[1] + wv[2] * wv[2];
    const T W[3][3] = {{T(0), -wv[2], wv[1]}, {wv[2], T(0), -wv[0]}, {-wv[1], wv[0], T(0)}};
    T W2[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) W2[i][j] = wv[i] * wv[j] - (i == j ? n2 : T(0));
    T ah, bh;
    d::left_jacobian_coefficients<T>(n2 * Ts * Ts, ah, bh);
    const T h1 = ah * Ts * Ts, h2 = bh * Ts * Ts * Ts;
    T Ig[3][3] = {}, Ib[3][3] = {}, IK[3][3] = {}, K[3][3], X[3][3], Y[3][3];
    const bool iso = filter_detail::is_isotropic3<T>(Q_gyro_);
    for (int k = 0; k < 5; ++k) {
        const T s = T(0.5) * Ts * (T(1) + T(d::kGL5Nodes[k]));
        const T wq = T(0.5) * Ts * T(d::kGL5Weights[k]);
        T as, bs;
        d::left_jacobian_coefficients<T>(n2 * s * s, as, bs);
        const T k0 = Ts - s, k1 = h1 - as * s * s, k2 = h2 - bs * s * s * s;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) K[i][j] = (i == j ? k0 : T(0)) + k1 * W[i][j] + k2 * W2[i][j];
        d::mul3<T>(K, Qb, false, X);
        d::mul3<T>(X, K, true, Y);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                Ib[i][j] += wq * Y[i][j];          // int K Qbg K^T
                IK[i][j] += wq * K[i][j];          // int K
            }
        if (!iso) {
            d::rows3<T>(lie::Exp<T>(Vector3(w * s)), K);
            d::mul3<T>(K, Qg, false, X);
            d::mul3<T>(X, K, true, Y);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) Ig[i][j] += wq * Y[i][j];   // int Exp Qg Exp^T
        }
    }
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) Ib[i][j] += Ig[i][j];
    T Qpp[3][3], Qpb[3][3], Qbb[3][3];
    d::mul3<T>(R0a, Ib, false, X);
    d::mul3<T>(X, R0a, true, Qpp);
    if (iso)
        for (int i = 0; i < 3; ++i) Qpp[i][i] += Qg[0][0] * Ts;
    d::mul3<T>(R0a, IK, false, X);
    d::mul3<T>(X, Qb, false, Y);
    d::mul3<T>(Y, R1a, true, Qpb);
    d::mul3<T>(R1a, Qb, false, X);
    d::mul3<T>(X, R1a, true, Qbb);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            q.Qpp(i,j) = T(0.5) * (Qpp[i][j] + Qpp[j][i]);
            q.Qpb(i,j) = -Qpb[i][j];
            q.Qbb(i,j) = T(0.5) * Ts * (Qbb[i][j] + Qbb[j][i]);
        }
    if (acc_bias_updates_enabled_) {
        const T tau_b = std::max(T(1e-3), tau_bacc_);
        const T qd_scale = -T(0.5) * tau_b * std::expm1(-T(2) * Ts / tau_b);
        d::rows3<T>(Q_bacc_, Qb);
        d::mul3<T>(R1a, Qb, false, X);
        d::mul3<T>(X, R1a, true, Y);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) q.Qaa(i,j) = T(0.5) * qd_scale * (Y[i][j] + Y[j][i]);
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
            fpb[k] = b.Phi_pb(k,0) * pbg[0] + b.Phi_pb(k,1) * pbg[1] + b.Phi_pb(k,2) * pbg[2];
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
void TFGC3_CLS::add_noise_(MatrixNX& P, const StepBlocks& b, const StepNoise& q, const MatC& Qc, const MatC* Qi) const {
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
            Matrix3 M = GQp[i] * b.Xs[j].transpose();
            if constexpr (kMotion)
                M += Qc(i,j) * motion_.excitation().re + (*Qi)(i,j) * motion_.excitation().im;
            else M += Qc(i,j) * Sigma_s_;
            add3(oi, oj, M);
        }
    }
    add3(OFF_BA, OFF_BA, q.Qaa);
}

TFGC3_TPL
void TFGC3_CLS::build_transition(const Vector3& gyr, T Ts, MatrixNX& Phi) const {
    Phi.setIdentity();
    if (!(Ts > T(0)) || !std::isfinite(Ts)) return;
    MatC Fc, Qc, Qi;
    if (!discretize_motion_(Ts, Fc, Qc, Qi)) return;
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
    MatC Fc, Qc, Qi;
    if (!discretize_motion_(Ts, Fc, Qc, Qi)) return;
    const Vector3 w = gyr - X_.B.col(0);
    const Eigen::Matrix<T,3,NC> W1 = X_.X * Fc.transpose();
    StepBlocks b;
    step_blocks_(w, Ts, X_.R, Fc, W1, b);
    StepNoise q;
    step_noise_(w, Ts, X_.R, b, q);
    add_noise_(Qd, b, q, Qc, kMotion ? &Qi : nullptr);
}

TFGC3_TPL
void TFGC3_CLS::time_update(Vector3 const& gyr, T Ts) {
    project_gyro_bias_();
    if (!(Ts > T(0)) || !std::isfinite(Ts) || !gyr.allFinite()) return;

    const Vector3 w = gyr - X_.B.col(0);
    last_gyr_bias_corrected_ = w;
    lever_.update_rate(w, Ts);
    const Matrix3 R0 = X_.R;

    refresh_discretization_(Ts);
    const MatC& Fc = motion_F_();
    const MatC& Qc = motion_Q_();

    // Mean: attitude by the bias-corrected rate, every world axis by the exact
    // CARMA transition, the accelerometer bias by its OU decay (held: frozen).
    StepBlocks b;
    Eigen::Matrix<T,3,NC> W1;
    {
        const T* const x = X_.X.data();            // X(ax, j) = x[ax + 3 j]
        const T* const f = Fc.data();
        T* const o = W1.data();
        for (int i = 0; i < NC; ++i)
            for (int ax = 0; ax < 3; ++ax) {
                T s = T(0);
                for (int j = 0; j < NC; ++j) s += x[ax + 3 * j] * f[i + NC * j];
                o[ax + 3 * i] = s;
            }
    }
    step_blocks_(w, Ts, R0, Fc, W1, b);
    X_.R = b.R1;
    reorthonormalize_();
    X_.X = W1;
    if (acc_bias_updates_enabled_) X_.B.col(1) *= b.phi_b;

    StepNoise q;
    step_noise_(w, Ts, R0, b, q);
    propagate_covariance_(b, q, Fc, Qc, motion_QI_());
}

TFGC3_TPL
void TFGC3_CLS::propagate_covariance_(const StepBlocks& b, const StepNoise& q, const MatC& Fc,
                                      const MatC& Qc, const MatC* Qi) {
    namespace d = tfg_carma_detail;
    // Phi by rows: phi' = phi + B bg, bg' = E bg, rho_i' = w_i x (B bg) +
    // sum_j F_ij rho_j ([w_i]x is b.Xs[i]), ba' = phi_b E ba.
    T B[3][3], E[3][3], Eb[3][3], F[NC][NC], w[NC][3];
    d::rows3<T>(b.Phi_pb, B);
    d::rows3<T>(b.E, E);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) Eb[i][j] = b.phi_b * E[i][j];
    for (int i = 0; i < NC; ++i) {
        for (int j = 0; j < NC; ++j) F[i][j] = Fc.data()[i + NC * j];
        const T* x = b.Xs[i].data();
        w[i][0] = x[2 + 3 * 1];
        w[i][1] = x[0 + 3 * 2];
        w[i][2] = x[1 + 3 * 0];
    }
    T* const p = P_.data();
    T* const m = M_scratch_.data();

    // M = Phi P, column by column.
    for (int c = 0; c < NX; ++c) {
        const T* pc = p + NX * c;
        T* mc = m + NX * c;
        const T g0 = pc[OFF_BG], g1 = pc[OFF_BG + 1], g2 = pc[OFF_BG + 2];
        const T f0 = B[0][0] * g0 + B[0][1] * g1 + B[0][2] * g2;
        const T f1 = B[1][0] * g0 + B[1][1] * g1 + B[1][2] * g2;
        const T f2 = B[2][0] * g0 + B[2][1] * g1 + B[2][2] * g2;
        mc[OFF_PHI] = pc[OFF_PHI] + f0;
        mc[OFF_PHI + 1] = pc[OFF_PHI + 1] + f1;
        mc[OFF_PHI + 2] = pc[OFF_PHI + 2] + f2;
        for (int k = 0; k < 3; ++k) mc[OFF_BG + k] = E[k][0] * g0 + E[k][1] * g1 + E[k][2] * g2;
        for (int i = 0; i < NC; ++i) {
            T a0 = w[i][1] * f2 - w[i][2] * f1;
            T a1 = w[i][2] * f0 - w[i][0] * f2;
            T a2 = w[i][0] * f1 - w[i][1] * f0;
            for (int j = 0; j < NC; ++j) {
                const T fij = F[i][j];
                const T* x = pc + OFF_W + 3 * j;
                a0 += fij * x[0];
                a1 += fij * x[1];
                a2 += fij * x[2];
            }
            mc[OFF_W + 3 * i] = a0;
            mc[OFF_W + 3 * i + 1] = a1;
            mc[OFF_W + 3 * i + 2] = a2;
        }
        const T h0 = pc[OFF_BA], h1 = pc[OFF_BA + 1], h2 = pc[OFF_BA + 2];
        for (int k = 0; k < 3; ++k) mc[OFF_BA + k] = Eb[k][0] * h0 + Eb[k][1] * h1 + Eb[k][2] * h2;
    }

    // Upper block triangle of M Phi^T, row by row: row r of block rb needs the
    // column blocks cb >= rb only.
    for (int r = 0; r < NX; ++r) {
        const int rb = r / 3;
        const T* mr = m + r;                       // M(r, c) = mr[NX * c]
        T* pr = p + r;
        if (rb < BLK_BA) {
            const T g0 = mr[NX * OFF_BG], g1 = mr[NX * (OFF_BG + 1)], g2 = mr[NX * (OFF_BG + 2)];
            const T t0 = g0 * B[0][0] + g1 * B[0][1] + g2 * B[0][2];
            const T t1 = g0 * B[1][0] + g1 * B[1][1] + g2 * B[1][2];
            const T t2 = g0 * B[2][0] + g1 * B[2][1] + g2 * B[2][2];
            if (rb == BLK_PHI) {
                pr[NX * OFF_PHI] = mr[NX * OFF_PHI] + t0;
                pr[NX * (OFF_PHI + 1)] = mr[NX * (OFF_PHI + 1)] + t1;
                pr[NX * (OFF_PHI + 2)] = mr[NX * (OFF_PHI + 2)] + t2;
            }
            if (rb <= BLK_BG)
                for (int k = 0; k < 3; ++k) pr[NX * (OFF_BG + k)] = g0 * E[k][0] + g1 * E[k][1] + g2 * E[k][2];
            for (int j = rb > BLK_W ? rb - BLK_W : 0; j < NC; ++j) {
                T a0 = w[j][1] * t2 - w[j][2] * t1;
                T a1 = w[j][2] * t0 - w[j][0] * t2;
                T a2 = w[j][0] * t1 - w[j][1] * t0;
                for (int l = 0; l < NC; ++l) {
                    const T fjl = F[j][l];
                    const T* x = mr + NX * (OFF_W + 3 * l);
                    a0 += fjl * x[0];
                    a1 += fjl * x[NX];
                    a2 += fjl * x[2 * NX];
                }
                T* y = pr + NX * (OFF_W + 3 * j);
                y[0] = a0;
                y[NX] = a1;
                y[2 * NX] = a2;
            }
        }
        const T h0 = mr[NX * OFF_BA], h1 = mr[NX * (OFF_BA + 1)], h2 = mr[NX * (OFF_BA + 2)];
        for (int k = 0; k < 3; ++k) pr[NX * (OFF_BA + k)] = h0 * Eb[k][0] + h1 * Eb[k][1] + h2 * Eb[k][2];
    }

    // Qd on the same blocks (add_noise_()): a product with [w_j]x^T is a cross
    // product of every row with w_j.
    T Qpp[3][3], Qpb[3][3], Qbb[3][3], Qaa[3][3], Sre[3][3], Sim[3][3];
    d::rows3<T>(q.Qpp, Qpp);
    d::rows3<T>(q.Qpb, Qpb);
    d::rows3<T>(q.Qbb, Qbb);
    d::rows3<T>(q.Qaa, Qaa);
    if constexpr (kMotion) {
        d::rows3<T>(motion_.excitation().re, Sre);
        d::rows3<T>(motion_.excitation().im, Sim);
    } else {
        d::rows3<T>(Sigma_s_, Sre);
    }
    auto add3 = [p](int r, int c, const T (&M)[3][3]) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) p[r + i + NX * (c + j)] += M[i][j];
    };
    // out(r, .) = w x a(r, .), or w x a(., r) when `transposed`.
    auto cross_rows = [](const T (&a)[3][3], bool transposed, const T* wv, T (&out)[3][3]) {
        for (int r = 0; r < 3; ++r) {
            const T x = transposed ? a[0][r] : a[r][0];
            const T y = transposed ? a[1][r] : a[r][1];
            const T z = transposed ? a[2][r] : a[r][2];
            out[r][0] = wv[1] * z - wv[2] * y;
            out[r][1] = wv[2] * x - wv[0] * z;
            out[r][2] = wv[0] * y - wv[1] * x;
        }
    };
    add3(OFF_PHI, OFF_PHI, Qpp);
    add3(OFF_PHI, OFF_BG, Qpb);
    add3(OFF_BG, OFF_BG, Qbb);
    add3(OFF_BA, OFF_BA, Qaa);
    T GQp[NC][3][3], M3[3][3];                     // GQp[i] = [w_i]x Qpp
    for (int i = 0; i < NC; ++i) {
        cross_rows(Qpp, false, w[i], M3);          // Qpp [w_i]x^T = GQp[i]^T
        add3(OFF_PHI, OFF_W + 3 * i, M3);
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) GQp[i][r][c] = M3[c][r];
        cross_rows(Qpb, true, w[i], M3);           // Qpb^T [w_i]x^T
        add3(OFF_BG, OFF_W + 3 * i, M3);
    }
    for (int i = 0; i < NC; ++i)
        for (int j = i; j < NC; ++j) {
            cross_rows(GQp[i], false, w[j], M3);
            const T qr = Qc.data()[i + NC * j];
            if constexpr (kMotion) {
                const T qi = Qi->data()[i + NC * j];
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c) M3[r][c] += qr * Sre[r][c] + qi * Sim[r][c];
            } else {
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c) M3[r][c] += qr * Sre[r][c];
            }
            add3(OFF_W + 3 * i, OFF_W + 3 * j, M3);
        }

    for (int c = 1; c < NX; ++c)
        for (int r = 0; r < c; ++r) p[c + NX * r] = p[r + NX * c];
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
    // P <- J P J^T with J = blockdiag(Jr) + C e_phi^T (C_0 = 0).  In 3 x 3
    // blocks, with W_i = Jr P_{i,phi} + C_i P_{phi,phi} / 2,
    //     P'_ij = Jr P_ij Jr^T + W_i C_j^T + C_i W_j^T,
    // evaluated in place on the upper block triangle and mirrored.
    // AttitudeOnly is J = blockdiag(Jr, I, ..): the attitude block row alone.
    namespace d = tfg_carma_detail;
    const bool full = reset_mode_ == ResetMode::Full;
    const T* const av = a.data();
    const T nphi[3] = {-av[OFF_PHI], -av[OFF_PHI + 1], -av[OFF_PHI + 2]};
    T ja, jb, Jr[3][3], C[NBLK][3][3];
    d::left_jacobian_coefficients<T>(nphi[0] * nphi[0] + nphi[1] * nphi[1] + nphi[2] * nphi[2], ja, jb);
    d::so3_series<T>(nphi, ja, jb, Jr);            // J_l(-phi)
    T* const p = P_.data();
    // W in the measurement scratch, free once the gain is applied.
    T* const wm = PHt_scratch_.data();             // W_i(r, k) = wm[3 i + r + NX k]
    if (full) {
        const d::RightJacobianCross<T> cross(Vector3(a.template segment<3>(OFF_PHI)));
        for (int blk = 1; blk < NBLK; ++blk) cross(av + 3 * blk, C[blk]);
        T Pff[3][3];
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) Pff[r][c] = T(0.5) * p[r + NX * c];
        for (int i = 0; i < NBLK; ++i) {
            const T* const x = p + 3 * i;          // P block (i, phi)
            for (int r = 0; r < 3; ++r)
                for (int k = 0; k < 3; ++k) {
                    T v = Jr[r][0] * x[NX * k] + Jr[r][1] * x[1 + NX * k] + Jr[r][2] * x[2 + NX * k];
                    if (i != 0) v += C[i][r][0] * Pff[0][k] + C[i][r][1] * Pff[1][k] + C[i][r][2] * Pff[2][k];
                    wm[3 * i + r + NX * k] = v;
                }
        }
    }
    const int rows = full ? NBLK : 1;
    for (int i = 0; i < rows; ++i)
        for (int j = i; j < NBLK; ++j) {
            T* const x = p + 3 * i + NX * 3 * j;   // P block (i, j)
            T Y[3][3], X[3][3];
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    Y[r][c] = Jr[r][0] * x[NX * c] + Jr[r][1] * x[1 + NX * c] + Jr[r][2] * x[2 + NX * c];
            if (full || j == 0) {
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c)
                        X[r][c] = Y[r][0] * Jr[c][0] + Y[r][1] * Jr[c][1] + Y[r][2] * Jr[c][2];
            } else {
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c) X[r][c] = Y[r][c];
            }
            if (full && j != 0) {
                const T* const wi = wm + 3 * i;
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c)
                        X[r][c] += wi[r] * C[j][c][0] + wi[r + NX] * C[j][c][1] + wi[r + 2 * NX] * C[j][c][2];
            }
            if (full && i != 0) {
                const T* const wj = wm + 3 * j;
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c)
                        X[r][c] += C[i][r][0] * wj[c] + C[i][r][1] * wj[c + NX] + C[i][r][2] * wj[c + 2 * NX];
            }
            // The block and its mirror; a diagonal block keeps its upper triangle.
            T* const xt = p + 3 * j + NX * 3 * i;
            for (int r = 0; r < 3; ++r)
                for (int c = (i == j ? r : 0); c < 3; ++c) {
                    x[r + NX * c] = X[r][c];
                    xt[c + NX * r] = X[r][c];
                }
        }
}

TFGC3_TPL
void TFGC3_CLS::inject(const Tangent& a) {
    // Xhat <- Exp_G(-a) o Xhat (lie/TwoFrameGroup.h), on plain storage:
    //     R <- Exp(-phi) R,   w_i <- Exp(-phi) w_i - J_l(-phi) rho_i,
    //     b_j <- b_j - R^T J_l(phi) beta_j     (R before its update).
    namespace d = tfg_carma_detail;
    const T* const av = a.data();
    const T phi[3] = {av[OFF_PHI], av[OFF_PHI + 1], av[OFF_PHI + 2]};
    const T t2 = phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2];
    T ea, eb, ja, jb, Em[3][3], Jm[3][3], Jp[3][3];
    d::exp_coefficients<T>(t2, ea, eb);
    d::left_jacobian_coefficients<T>(t2, ja, jb);
    d::so3_series<T>(phi, -ea, eb, Em);            // Exp(-phi)
    d::so3_series<T>(phi, -ja, jb, Jm);            // J_l(-phi)
    d::so3_series<T>(phi, ja, jb, Jp);             // J_l(phi)
    T* const R = X_.R.data();                      // column-major
    T* const B = X_.B.data();
    for (int j = 0; j < Group::kBiasCols; ++j) {
        const T* const beta = av + Group::bias_offset(j);
        T v[3];
        for (int k = 0; k < 3; ++k) v[k] = Jp[k][0] * beta[0] + Jp[k][1] * beta[1] + Jp[k][2] * beta[2];
        for (int k = 0; k < 3; ++k) B[k + 3 * j] -= R[3 * k] * v[0] + R[1 + 3 * k] * v[1] + R[2 + 3 * k] * v[2];
    }
    T* const W = X_.X.data();
    for (int i = 0; i < NC; ++i) {
        const T* const rho = av + Group::world_offset(i);
        const T x = W[3 * i], y = W[3 * i + 1], z = W[3 * i + 2];
        for (int k = 0; k < 3; ++k)
            W[k + 3 * i] = Em[k][0] * x + Em[k][1] * y + Em[k][2] * z
                         - (Jm[k][0] * rho[0] + Jm[k][1] * rho[1] + Jm[k][2] * rho[2]);
    }
    for (int c = 0; c < 3; ++c) {
        const T x = R[3 * c], y = R[3 * c + 1], z = R[3 * c + 2];
        for (int k = 0; k < 3; ++k) R[k + 3 * c] = Em[k][0] * x + Em[k][1] * y + Em[k][2] * z;
    }
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
    {
        const T r0 = r(0), r1 = r(1), r2 = r(2);
        const T* const k = K.data();
        T* const av = a.data();
        for (int i = 0; i < NX; ++i) av[i] = k[i] * r0 + k[i + NX] * r1 + k[i + 2 * NX] * r2;
    }
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
    // Off-CoG IMU: remove the rotational term L at the sensor (see the header).
    const bool lever = lever_.enabled();
    const Vector3 f_cog = lever ? Vector3(acc_meas - lever_.acceleration()) : acc_meas;
    if constexpr (kMotion) motion_.observe(R * (f_cog - ba_term) + g_world);
    const Vector3 r = R * (f_cog - ba_term) - (Vector3(X_.X.col(IA)) - g_world);
    last_acc_diag_.r = r;
    const Matrix3 Rw = R * Racc_ * R.transpose();

    // H = [ [g]x | H_bg | 0 0 -I 0.. | -I ];  [g]x = [[0,-g,0],[g,0,0],[0,0,0]],
    // H_bg = -R J_L R^T with a lever arm, 0 without.
    const T g = gravity_magnitude_;
    MatrixNX3& PHt = PHt_scratch_;
    const T* const p = P_.data();
    T* const ht = PHt.data();
    for (int i = 0; i < NX; ++i) {
        ht[i] = -g * p[i + NX * (OFF_PHI + 1)] - p[i + NX * OFF_A] - p[i + NX * OFF_BA];
        ht[i + NX] = g * p[i + NX * OFF_PHI] - p[i + NX * (OFF_A + 1)] - p[i + NX * (OFF_BA + 1)];
        ht[i + 2 * NX] = -p[i + NX * (OFF_A + 2)] - p[i + NX * (OFF_BA + 2)];
    }
    Matrix3 H_bg = Matrix3::Zero();
    if (lever) {
        H_bg = -(R * lever_.bias_jacobian() * R.transpose());
        T Hb[3][3];
        tfg_carma_detail::rows3<T>(H_bg, Hb);
        for (int b = 0; b < 3; ++b)
            for (int i = 0; i < NX; ++i)
                ht[i + NX * b] += p[i + NX * OFF_BG] * Hb[b][0] + p[i + NX * (OFF_BG + 1)] * Hb[b][1] +
                                  p[i + NX * (OFF_BG + 2)] * Hb[b][2];
    }
    Matrix3 S = Rw;
    for (int b = 0; b < 3; ++b) {
        const T* const hb = ht + NX * b;
        S(0,b) += -g * hb[OFF_PHI + 1] - hb[OFF_A] - hb[OFF_BA];
        S(1,b) += g * hb[OFF_PHI] - hb[OFF_A + 1] - hb[OFF_BA + 1];
        S(2,b) += -hb[OFF_A + 2] - hb[OFF_BA + 2];
    }
    if (lever) {
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                S(a,b) += H_bg(a,0) * PHt(OFF_BG, b) + H_bg(a,1) * PHt(OFF_BG + 1, b) +
                          H_bg(a,2) * PHt(OFF_BG + 2, b);
    }
    for (int a = 0; a < 3; ++a)
        for (int b = a + 1; b < 3; ++b) {
            const T v = T(0.5) * (S(a,b) + S(b,a));
            S(a,b) = v;
            S(b,a) = v;
        }
    if (!use_ba) {
        // CARMA3's held bias: in S, not in the cross-covariance or the gain.
        for (int k = 0; k < 3; ++k)
            for (int i = 0; i < NX; ++i) ht[i + NX * k] += p[i + NX * (OFF_BA + k)];
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
    T Hp[3][3];
    tfg_carma_detail::rows3<T>(Hphi, Hp);
    const T* const p = P_.data();
    T* const ht = PHt.data();
    for (int b = 0; b < 3; ++b)
        for (int i = 0; i < NX; ++i)
            ht[i + NX * b] = p[i + NX * OFF_PHI] * Hp[b][0] + p[i + NX * (OFF_PHI + 1)] * Hp[b][1] +
                             p[i + NX * (OFF_PHI + 2)] * Hp[b][2];
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
