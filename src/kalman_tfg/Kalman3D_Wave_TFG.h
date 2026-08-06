#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    Right-invariant two-frame Lie-group error-state EKF for marine wave state.

    Same physics as Kalman3D_Wave_OU_III -- Ornstein-Uhlenbeck world
    acceleration, the v -> p -> S integral chain, random-walk biases -- but
    carried on the two-frame group in src/lie/TwoFrameGroup.h instead of a
    quaternion MEKF with additive translational states.

    CONVENTION CONTRACT (see docs/tfg-design.md; three parts, always together):

        1. rotation direction : R = R_bw, body to world
        2. error definition   : right-invariant, eta = Xhat o X^-1
        3. correction side    : left, Xhat+ = Exp_G(dxi) o Xhat-

    The public accessors keep OU-III's API. quaternion() returns body-to-world,
    exactly as Kalman3D_Wave_OU_III::quaternion() does (it conjugates its
    world-to-body qref to get there; this filter already stores that
    direction). So IW3dFusionAdapter, FilterSnapshot, the plots and the
    sketches see no change.

    STATE, in OU-III's ordering so the world block stays contiguous and the
    block structure lines up index-for-index:

        xi = [ phi | delta_bg | rho_v rho_p rho_S rho_a | delta_ba ]
               0     3          6     9     12    15      18

    ------------------------------------------------------------------------
    ERROR DYNAMICS (Level 1: world vectors on the group, additive bias errors)
    ------------------------------------------------------------------------

    Attitude. With Rdot = R [w_m - b_g - n_g]x and Rhatdot = Rhat [w_m - bhat_g]x,
    and delta_bg = bhat_g - b_g in the body frame,

        d(dR)/dt = Rhat ( [w_m - bhat_g]x - [w_m - b_g - n_g]x ) R^T
                 = Rhat [ -delta_bg + n_g ]x R^T

    and since Rhat [u]x Rhat^T = [Rhat u]x, to first order

        phidot = -Rhat delta_bg + Rhat n_g                                 (1)

    Note what is NOT there: no [w]x phi term. The attitude error is driftless
    in this parameterization, which is the structural gain over the MEKF -- it
    is why Phi_phi,phi = I here where OU-III needs an exact Rstep.

    World vectors. With rho_x = xhat - dR x and d(dR)/dt = [phidot]x dR,

        rho_xdot = xhatdot - [phidot]x x - dR xdot
                 = (chain term) + [x]x phidot                              (2)

    using -[u]x x = [x]x u. Substituting (1), every world vector picks up the
    gyro bias through -[x]x Rhat. The chain terms are the OU chain itself:

        rho_vdot = rho_a,  rho_pdot = rho_v,  rho_Sdot = rho_p,
        rho_adot = -lambda rho_a

    which is exactly what ou_detail::IntegratedOUChain<T,3> already
    discretizes. That is the reuse this design is built around: the world
    block of Phi and Qd is the OU-III code, unmodified, and
    ou_chain_identity-test asserts it bit for bit.

    ------------------------------------------------------------------------
    DISCRETIZATION
    ------------------------------------------------------------------------

    The continuous generator is block lower triangular in (phi, delta_bg,
    world), with A_phi,phi = 0 and A_bg,bg = 0. That makes several blocks
    exactly integrable rather than truncated:

        Phi_phi,phi   = I
        Phi_phi,bg    = -Rbar h                        exact
        Phi_bg,bg     = I
        Phi_world     = Phi_OU(h, tau)                 exact, IntegratedOUChain
        Phi_world,bg  = -Psi(h, tau) [x]x Rbar
        Phi_world,phi = 0

    where Psi(h,tau) = int_0^h Phi_OU(u) du, in closed form with Taylor
    branches (integral_transition_axis below).

    Rbar is the *average* body-to-world rotation over the step, not Rhat at
    its start. This matters and is easy to get wrong. The attitude propagates
    as Rhat(s) = Rhat Exp(w s) while (1) says phidot = -Rhat(s) delta_bg, so

        phi(h) = phi(0) - ( int_0^h Rhat Exp(w s) ds ) delta_bg
               = phi(0) - Rhat h J_l(w h) delta_bg

    using J_l(theta) = int_0^1 Exp(u theta) du. So Rbar = Rhat J_l(w h). Using
    Rhat alone costs a relative error of order |w| h -- 0.25% at 1 rad/s and a
    5 ms step, small but well above the tolerance a finite-difference check
    should be run at, and it grows with turn rate exactly when attitude and
    bias are hardest to separate.

    Doing Phi exactly rather than truncating at h^2 is what lets
    tfg_propagation-test run at a tight tolerance and at large h, instead of
    being loosened until it stops discriminating. Phi_world,bg keeps one
    documented approximation: it uses the step-averaged Rbar rather than
    correlating Rhat(s)'s variation with the Phi_OU weighting and with the
    drift of x(s) across the step. Both neglected effects are O((|w| h)^2) and
    O(h) on a term that is already a product of two small quantities.

    Qd is assembled blockwise. The world block carries the exact OU
    contribution from IntegratedOUChain; the attitude and bias blocks carry
    the exactly integrated bias random-walk terms (h^2/2 and h^3/3); the gyro
    white-noise leakage into the world vectors is first order in h. That last
    truncation is deliberate and matches the surrounding code's standard --
    Kalman3D_Wave_OU_III.h:1560 does the same, and Qd magnitudes are tuned
    quantities, unlike Phi which has to match the propagation it linearizes.

    ------------------------------------------------------------------------

    Scalar-templated: firmware runs float, tests run double.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>

#include "kalman_ou_common/KalmanOUCoreMath.h"
#include "lie/SO3Jacobians.h"
#include "lie/TwoFrameGroup.h"

namespace ocean_imu::kalman {

namespace tfg_detail {

using ocean_imu::kalman::ou_detail::safe_inv_tau;

/*
    Psi(h, tau) = int_0^h Phi_OU(u) du for one axis, where Phi_OU is the 4x4
    [v p S a] transition of IntegratedOUChain<T,3>.

    Closed forms, with x = h/tau:

        Psi_00 = h            Psi_10 = h^2/2        Psi_20 = h^3/6
        Psi_11 = h            Psi_21 = h^2/2        Psi_22 = h
        Psi_33 = -tau expm1(-x)
        Psi_03 = tau^2 ( expm1(-x) + x )
        Psi_13 = tau^3 ( x^2/2 - expm1(-x) - x )
        Psi_23 = tau^4 ( x^3/6 - x^2/2 + expm1(-x) + x )

    The last three cancel to O(x^2), O(x^3), O(x^4) respectively, so the
    direct forms lose all significance for small steps -- at h/tau = 1e-3,
    Psi_23's leading term is 1e-12 of the quantities being subtracted. Hence
    the Taylor branch, on the same 1e-2 threshold the rest of the codebase
    uses.
*/
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
        const T x2 = x * x, x3 = x2 * x, x4 = x3 * x, x5 = x4 * x;
        const T x6 = x5 * x, x7 = x6 * x, x8 = x7 * x;
        // Each series starts one power later than the one above it, so the
        // number of terms needed to reach a given relative accuracy is the
        // same in every case but the absolute order climbs. Psi_23's leading
        // term is x^4/24, and its first dropped term sets the error: stopping
        // at x^6 leaves x^3/210 relative, which is 4.6e-9 right at the 1e-2
        // seam -- measurable, so it is carried to x^8.
        // -tau expm1(-x) = tau (x - x^2/2 + x^3/6 - ...)
        Psi(3,3) = tau  * (x - x2/T(2) + x3/T(6) - x4/T(24) + x5/T(120) - x6/T(720));
        // expm1(-x) + x = x^2/2 - x^3/6 + x^4/24 - ...
        Psi(0,3) = tau2 * (x2/T(2) - x3/T(6) + x4/T(24) - x5/T(120) + x6/T(720)
                           - x7/T(5040));
        // x^2/2 - expm1(-x) - x = x^3/6 - x^4/24 + ...
        Psi(1,3) = tau3 * (x3/T(6) - x4/T(24) + x5/T(120) - x6/T(720)
                           + x7/T(5040) - x8/T(40320));
        // x^3/6 - x^2/2 + expm1(-x) + x = x^4/24 - x^5/120 + ...
        Psi(2,3) = tau4 * (x4/T(24) - x5/T(120) + x6/T(720) - x7/T(5040)
                           + x8/T(40320));
    } else {
        const T em1 = std::expm1(-x);
        const T x2 = x * x, x3 = x2 * x;
        Psi(3,3) = -tau  * em1;
        Psi(0,3) =  tau2 * (em1 + x);
        Psi(1,3) =  tau3 * (x2/T(2) - em1 - x);
        Psi(2,3) =  tau4 * (x3/T(6) - x2/T(2) + em1 + x);
    }
}

} // namespace tfg_detail

/*
    two_frame_bias selects the bias geometry:

      true  (Level 2) -- full two-frame group, beta_b = R (bhat - b). Default.
      false (Level 1) -- world vectors on SE_4(3), bias errors additive in the
                         body frame. Kept as a build flag: it isolates
                         group-retraction defects from bias-geometry defects,
                         and it is the ablation that answers whether the bias
                         geometry earns its keep.

    THE TWO LEVELS ARE THE SAME FILTER IN DIFFERENT COORDINATES, related by
    the exact linear map beta = Rhat delta_b:

        T     = blkdiag(I, Rhat, I, I, I, I, Rhat)
        xi_L2 = T xi_L1,   P_L2 = T P_L1 T^T,   H_L1 = H_L2 T

    Under pure propagation that relation is exact rather than approximate, and
    covariance_transport-test asserts it -- which cross-checks the whole
    Level-2 derivation against the already-verified Level-1 one. Note the
    conjugation is time-varying: Phi_L2 = T(h) Phi_L1 T(0)^-1, which is
    precisely why Phi_beta,beta is Exp(w_world h) here where Level 1 has I.

    Under measurement updates the two agree only to first order, because the
    Level-2 injection carries a J_l(-phi) that the Level-1 injection does not.
*/
template <typename T = float,
          bool with_gyro_bias = true,
          bool with_accel_bias = true,
          bool two_frame_bias = true>
class Kalman3D_Wave_TFG {
    static_assert(with_gyro_bias,
                  "the two-frame tangent layout pins the gyro bias at offset 3; "
                  "a no-gyro-bias variant would need its own layout");

  public:
    using Group   = ocean_imu::lie::TwoFrameGroup<T, 4, with_accel_bias ? 2 : 1>;
    using Scalar  = T;
    using Vector3 = Eigen::Matrix<T,3,1>;
    using Matrix3 = Eigen::Matrix<T,3,3>;

    static constexpr int NX = Group::NX;               // 21, or 18 without b_a

    static constexpr int OFF_PHI = Group::OFF_PHI;     // 0
    static constexpr int OFF_BG  = Group::OFF_BG;      // 3
    static constexpr int OFF_V   = Group::world_offset(0);   // 6
    static constexpr int OFF_P   = Group::world_offset(1);   // 9
    static constexpr int OFF_S   = Group::world_offset(2);   // 12
    static constexpr int OFF_AW  = Group::world_offset(3);   // 15
    static constexpr int OFF_BA  = with_accel_bias ? Group::OFF_BA : -1;  // 18

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

    // ---------------------------------------------------------------- state

    void initialize_identity() {
        X_ = Group::Identity();
        P_.setIdentity();
        P_ *= T(1e-4);
    }

    // Truth injection, for tests and for the paired simulator.
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

    // --------------------------------------------------------- accessors
    // Deliberately the same shapes and frames as Kalman3D_Wave_OU_III, so the
    // two filters are drop-in comparable in the simulator and the plots.

    // BODY -> WORLD, matching Kalman3D_Wave_OU_III::quaternion().
    [[nodiscard]] Eigen::Quaternion<T> quaternion() const {
        Eigen::Quaternion<T> q(X_.R);
        q.normalize();
        return q;
    }
    [[nodiscard]] Matrix3 R_bw() const { return X_.R; }
    [[nodiscard]] Matrix3 R_wb() const { return X_.R.transpose(); }

    [[nodiscard]] Vector3 get_velocity() const             { return X_.X.col(0); }
    [[nodiscard]] Vector3 get_position() const             { return X_.X.col(1); }
    [[nodiscard]] Vector3 get_integral_displacement() const{ return X_.X.col(2); }
    [[nodiscard]] Vector3 get_world_accel() const          { return X_.X.col(3); }
    [[nodiscard]] Vector3 gyroscope_bias() const           { return X_.B.col(0); }
    [[nodiscard]] Vector3 get_acc_bias() const {
        if constexpr (with_accel_bias) return X_.B.col(1);
        else return Vector3::Zero();
    }

    [[nodiscard]] const MatrixNX& covariance_full() const { return P_; }
    [[nodiscard]] MatrixNX& covariance_full() { return P_; }
    [[nodiscard]] const Group& state() const { return X_; }
    [[nodiscard]] Group& state() { return X_; }

    // ------------------------------------------------------------- tuning

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
        Q_gyro_ = Matrix3::Identity() * (d * d);
    }
    void set_Q_bgyro_rw(const Vector3& var_per_s) {
        Q_bg_ = var_per_s.cwiseAbs().asDiagonal();
    }
    void set_Q_bacc_rw(const Vector3& var_per_s) {
        Q_ba_ = var_per_s.cwiseAbs().asDiagonal();
    }

    void set_initial_linear_uncertainty(T std_v, T std_p, T std_S, T std_aw) {
        P_.template block<12,12>(OFF_V, OFF_V).setZero();
        const T d[4] = {std_v*std_v, std_p*std_p, std_S*std_S, std_aw*std_aw};
        for (int c = 0; c < 4; ++c) {
            for (int k = 0; k < 3; ++k) P_(OFF_V + 3*c + k, OFF_V + 3*c + k) = d[c];
        }
    }

    // --------------------------------------------------------- propagation

    /*
        Nominal (mean) propagation.

            Rhat+   = Rhat Exp((w_m - bhat_g) h)
            [v p S a_w]+ = Phi_OU applied per axis
            biases held constant

        R is body-to-world and the rate is body-frame, so the rotation
        increment multiplies on the RIGHT here. That is not in tension with
        the left-multiplicative *correction*: propagation moves the nominal
        state along the trajectory, while a correction acts on the error,
        which lives in the world frame. Mixing the two up is the single
        easiest way to get this filter subtly wrong.
    */
    void propagate_mean(const Vector3& gyro_meas, T Ts) {
        if (!(Ts > T(0)) || !std::isfinite(Ts)) return;

        const Vector3 omega = gyro_meas - X_.B.col(0);
        last_omega_ = omega;

        X_.R = X_.R * ocean_imu::lie::Exp<T>(Vector3(omega * Ts));
        reorthonormalize_();

        Eigen::Matrix<T,4,4> Phi_axis;
        ou_detail::IntegratedOUChain<T,3>::transition(tau_aw_, Ts, Phi_axis);
        // Row i of X is [v_i p_i S_i a_i], so applying the per-axis 4x4 to
        // every row at once is a right multiplication by its transpose.
        X_.X = (X_.X * Phi_axis.transpose()).eval();
    }

    /*
        Discrete transition of the right-invariant error, assembled blockwise.
        Exposed so tfg_propagation-test can check it against finite
        differences of the nonlinear propagation, and ou_chain_identity-test
        can check the world block against IntegratedOUChain directly.
    */
    void build_transition(T Ts, const Vector3& gyro_meas, MatrixNX& Phi) const {
        Phi.setIdentity();

        /*
            Step-averaged body-to-world rotation:

                (1/h) int_0^h Rhat Exp(w s) ds = Rhat J_l(w h) =: Rbar

            Using Rhat alone is wrong at order |w| h -- see the header note.

            Bmap is the map from a bias error onto phidot, integrated over the
            step. The two levels differ by exactly one transpose of Rhat,
            which is the whole content of the coordinate change beta = Rhat
            delta_b:

                Level 1:  phidot = -Rhat delta_b   ->  Bmap = Rbar
                Level 2:  phidot = -beta           ->  Bmap = Rbar Rhat^T
                                                            = J_l(w_world h)

            using the equivariance J_l(R u) = R J_l(u) R^T.
        */
        const Vector3 omega = gyro_meas - X_.B.col(0);
        const Matrix3 Rbar =
            X_.R * ocean_imu::lie::left_jacobian<T>(Vector3(omega * Ts));
        const Matrix3 Bmap = two_frame_bias ? Matrix3(Rbar * X_.R.transpose()) : Rbar;

        // Attitude <- bias. Exact: A_phi,phi = 0, and the bias block's own
        // dynamics are a pure rotation, handled below.
        Phi.template block<3,3>(OFF_PHI, OFF_BG) = -Bmap * Ts;

        /*
            Bias self-dynamics.

            Level 1: delta_b is a random walk in the body frame, so Phi = I.

            Level 2: beta = R delta_b inherits the frame's rotation --

                betadot = Rdot delta_b + R delta_bdot
                        = R [w]x delta_b - R w_b
                        = [R w]x beta - R w_b

            so Phi_beta,beta = Exp(w_world h), a rotation rather than the
            identity. This term has no Level-1 analogue and is the substantive
            new content of the two-frame geometry: the bias error is carried
            along by the vehicle's rotation instead of sitting still.
        */
        if constexpr (two_frame_bias) {
            const Matrix3 Rw_step = ocean_imu::lie::Exp<T>(Vector3(X_.R * omega * Ts));
            Phi.template block<3,3>(OFF_BG, OFF_BG) = Rw_step;
            if constexpr (with_accel_bias) {
                Phi.template block<3,3>(OFF_BA, OFF_BA) = Rw_step;
            }
        }

        // World block: the OU chain, exact, per axis.
        Eigen::Matrix<T,4,4> Phi_axis;
        ou_detail::IntegratedOUChain<T,3>::transition(tau_aw_, Ts, Phi_axis);
        Phi.template block<12,12>(OFF_V, OFF_V).setZero();
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (Phi_axis(i,j) == T(0)) continue;
                for (int k = 0; k < 3; ++k) {
                    Phi(OFF_V + 3*i + k, OFF_V + 3*j + k) = Phi_axis(i,j);
                }
            }
        }

        // World <- gyro bias:  Phi_world,bg = -Psi(h) [x]x Rbar, stacked over
        // the world columns x in {v, p, S, a_w}. This is where rho_xdot picks
        // up +[x]x phidot from (2).
        Eigen::Matrix<T,4,4> Psi;
        tfg_detail::integral_transition_axis<T>(tau_aw_, Ts, Psi);
        Matrix3 Gx[4];
        for (int c = 0; c < 4; ++c) {
            Gx[c] = -ocean_imu::lie::skew<T>(Vector3(X_.X.col(c))) * Bmap;
        }
        for (int i = 0; i < 4; ++i) {
            Matrix3 acc = Matrix3::Zero();
            for (int j = 0; j < 4; ++j) {
                if (Psi(i,j) != T(0)) acc.noalias() += Psi(i,j) * Gx[j];
            }
            Phi.template block<3,3>(OFF_V + 3*i, OFF_BG) = acc;
        }

        // Accelerometer bias is a pure random walk and does not enter the
        // propagation of anything else -- it only shows up in the
        // accelerometer measurement, which is Phase C.
    }

    /*
        Discrete process noise.

        World block: exact, from IntegratedOUChain, including the correlated
        vector-OU assembly when Sigma_aw is not diagonal.

        Attitude and gyro-bias blocks: exactly integrated, since
        Phi_phi,bg(s) = -Rhat (h - s) is exact --

            Q_phi,phi += W h + Rhat Q_bg Rhat^T h^3/3
            Q_phi,bg   = -Rhat Q_bg h^2/2
            Q_bg,bg    = Q_bg h                with W = Rhat Q_gyro Rhat^T

        Gyro white noise leaking into the world vectors is first order in h.
        See the header note on why that truncation is deliberate.
    */
    void build_process_noise(T Ts, const Vector3& gyro_meas, MatrixNX& Qd) const {
        Qd.setZero();
        if (!(Ts > T(0)) || !std::isfinite(Ts)) return;

        const Matrix3 Rhat = X_.R;
        /*
            At Level 2 the bias noise lands in the world frame of the END of
            the step, not its start.

            The bias random walk is body-frame, and beta carries it forward
            through the frame's rotation:

                Q_beta = int_0^h Exp(w_world (h-s)) Rhat(s) Q_b Rhat(s)^T
                                 Exp(w_world (h-s))^T ds

            and because Exp(w_world (h-s)) Rhat(s) = Rhat(h) identically, the
            integrand is the constant Rhat(h) Q_b Rhat(h)^T. Using Rhat(0)
            instead is wrong at order |w| h -- the same order, and for the same
            reason, as using Rhat rather than Rbar in Phi.
        */
        const Matrix3 R_end =
            two_frame_bias
                ? Matrix3(Rhat * ocean_imu::lie::Exp<T>(
                      Vector3((gyro_meas - X_.B.col(0)) * Ts)))
                : Rhat;
        const Matrix3 W  = Rhat * Q_gyro_ * Rhat.transpose();   // world-frame gyro noise
        const Matrix3 Wb = Rhat * Q_bg_   * Rhat.transpose();

        const T h2 = Ts * Ts;
        const T h3 = h2 * Ts;

        /*
            Attitude and gyro-bias.

            The bias random walk drives phi through Phi_phi,bg(s), which is
            exactly integrable, giving the h^3/3 and h^2/2 terms.

            The bias blocks themselves are expressed in whichever frame the
            error lives in: body at Level 1, world at Level 2. That is the
            T-conjugation again -- Q_L2 = T Q_L1 T^T -- and it is why the
            Level-2 forms carry Rhat on both sides where Level 1 carries it on
            one or neither.
        */
        // phi's own block is frame-independent (T's phi block is the
        // identity), so it uses Rhat at the start of the step in both levels.
        // The bias blocks land at R_end.
        const Matrix3 Q_bias_own   = two_frame_bias
            ? Matrix3(R_end * Q_bg_ * R_end.transpose()) : Q_bg_;
        const Matrix3 Q_bias_cross = two_frame_bias
            ? Matrix3(Rhat * Q_bg_ * R_end.transpose()) : Matrix3(Rhat * Q_bg_);

        Qd.template block<3,3>(OFF_PHI, OFF_PHI) = W * Ts + Wb * (h3 / T(3));
        Qd.template block<3,3>(OFF_PHI, OFF_BG)  = -Q_bias_cross * (h2 / T(2));
        Qd.template block<3,3>(OFF_BG, OFF_PHI)  =
            Qd.template block<3,3>(OFF_PHI, OFF_BG).transpose();
        Qd.template block<3,3>(OFF_BG, OFF_BG)   = Q_bias_own * Ts;

        // Gyro white noise into the world vectors, via rho_xdot = ... + [x]x phidot.
        Matrix3 Sx[4];
        for (int c = 0; c < 4; ++c) Sx[c] = ocean_imu::lie::skew<T>(Vector3(X_.X.col(c)));
        for (int i = 0; i < 4; ++i) {
            const Matrix3 SiW = Sx[i] * W;
            Qd.template block<3,3>(OFF_PHI, OFF_V + 3*i) += (W * Sx[i].transpose()) * Ts;
            Qd.template block<3,3>(OFF_V + 3*i, OFF_PHI) += (SiW) * Ts;
            for (int j = 0; j < 4; ++j) {
                Qd.template block<3,3>(OFF_V + 3*i, OFF_V + 3*j) +=
                    (SiW * Sx[j].transpose()) * Ts;
            }
        }

        // World block: the exact OU contribution, added on top.
        add_ou_process_noise_(Ts, Qd);

        if constexpr (with_accel_bias) {
            Qd.template block<3,3>(OFF_BA, OFF_BA) =
                (two_frame_bias ? Matrix3(R_end * Q_ba_ * R_end.transpose()) : Q_ba_) * Ts;
        }

        Qd = T(0.5) * (Qd + Qd.transpose()).eval();
        ou_detail::project_psd_ou_iii<T,NX>(Qd);
    }

    /*
        Right-invariant error of this state relative to a reference, in the
        same TANGENT coordinates that inject() consumes and that P is
        expressed in. Exact inverse of inject().

            phi        = Log(Rhat R^T)
            rho_group  = xhat - Exp(phi) x
            xi_rho     = J_l(phi)^-1 rho_group
            xi_db      = bhat - b               (additive, body frame)

        The J_l inverse is the part that is easy to drop, and dropping it is
        not harmless. There are two different objects here:

          - the group element eta = Xhat o X^-1, whose world components are
            rho_group = xhat - dR x;
          - its tangent coordinates xi = Log_G(eta), whose world components
            are J_l(phi)^-1 rho_group.

        The covariance lives on the tangent, and the retraction is Exp_G, so
        the injection multiplies by J_l. Returning rho_group here instead
        would make error_from and inject differ by exactly J_l -- agreeing to
        first order, so every small-perturbation check would still pass, while
        finite-difference Jacobians at usable step sizes and NEES at realistic
        error magnitudes would both be quietly wrong.
    */
    [[nodiscard]] Tangent error_from(const Kalman3D_Wave_TFG& reference) const {
        Tangent xi;
        xi.setZero();

        const Vector3 phi = ocean_imu::lie::Log<T>(
            Matrix3(X_.R * reference.X_.R.transpose()));
        const Matrix3 dR  = ocean_imu::lie::Exp<T>(phi);
        const Matrix3 Jli = ocean_imu::lie::inv_left_jacobian<T>(phi);
        xi.template segment<3>(OFF_PHI) = phi;

        for (int c = 0; c < 4; ++c) {
            xi.template segment<3>(OFF_V + 3*c) =
                Jli * (X_.X.col(c) - dR * reference.X_.X.col(c));
        }
        if constexpr (two_frame_bias) {
            /*
                The group element's bias part is eta_B = R_ref (Bhat - B_ref),
                and its tangent coordinates apply J_l(-phi)^-1, mirroring the
                J_l(-phi) that Exp_G puts there.

                R_ref, not Rhat: eta = Xhat o X^-1 evaluates the composition's
                R2^T B1 term with R2 = X^-1's rotation, i.e. the reference's.
                The two differ by Exp(phi), so using the estimate's rotation
                would again be right only at zero error.
            */
            const Matrix3 Jri = ocean_imu::lie::inv_left_jacobian<T>(Vector3(-phi));
            const Matrix3 JriR = Jri * reference.X_.R;
            xi.template segment<3>(OFF_BG) = JriR * (X_.B.col(0) - reference.X_.B.col(0));
            if constexpr (with_accel_bias) {
                xi.template segment<3>(OFF_BA) = JriR * (X_.B.col(1) - reference.X_.B.col(1));
            }
        } else {
            xi.template segment<3>(OFF_BG) = X_.B.col(0) - reference.X_.B.col(0);
            if constexpr (with_accel_bias) {
                xi.template segment<3>(OFF_BA) = X_.B.col(1) - reference.X_.B.col(1);
            }
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
        P_ = T(0.5) * (P_ + P_.transpose()).eval();
    }

    // ---------------------------------------------------------- correction

    /*
        State injection, Xhat+ = Exp_G(xi) o Xhat-.

        The world part is the group retraction at both levels. The bias line
        is what differs:

            Level 1:  Bhat+ = Bhat + delta_b                 additive, body frame
            Level 2:  Bhat+ = Bhat + Rhat^- ^T J_l(-phi) beta

        Level 2 is the two-frame group law, B_result = B2 + R2^T B1, with the
        PRE-update Rhat -- which is why R_pre is captured before the attitude
        is touched. Using the post-update rotation would be wrong by exactly
        Exp(phi), i.e. right only when the correction is zero.

        The split lives here rather than in TwoFrameGroup so that the group
        stays a pure Level-2 object with no filter-level switch in it.

        Note what is deliberately absent: there is no MEKF-style covariance
        reset (ou_detail::apply_left_error_reset). The invariant update already
        leaves P in the post-update tangent frame, and applying
        G = I + 0.5 [dtheta]x on top would double-correct it.
    */
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
            // E_B = J_l(-phi) beta, transported into the body frame by R_pre^T.
            const Matrix3 Jr = ocean_imu::lie::left_jacobian<T>(Vector3(-phi));
            const Matrix3 RtJr = R_pre.transpose() * Jr;
            X_.B.col(0) += RtJr * xi.template segment<3>(OFF_BG);
            if constexpr (with_accel_bias) {
                X_.B.col(1) += RtJr * xi.template segment<3>(OFF_BA);
            }
        } else {
            X_.B.col(0) += xi.template segment<3>(OFF_BG);
            if constexpr (with_accel_bias) X_.B.col(1) += xi.template segment<3>(OFF_BA);
        }
    }

    // ------------------------------------------------------- measurements
    /*
        SIGN CONVENTION, fixed here once and used everywhere.

            xi     is the error OF THE ESTIMATE:  eta = Xhat X^-1 = Exp_G(xi),
                   and P = E[xi xi^T].
            r      is the innovation, measured minus predicted.
            r      ~ H xi + noise                       (note: +, not -)
            K      = P H^T (H P H^T + R)^-1
            xihat  = K r                                 estimated error
            Xhat+  = Exp_G(-xihat) o Xhat                remove it

        The correction is injected NEGATED, because xi is the error rather
        than the correction. Half the sign bugs in filters like this come from
        mixing the two conventions -- deriving H against "r ~ -H xi" and then
        injecting +K r, which cancels out only if both mistakes are made
        together. Here H is exactly d r / d xi, which is what
        tfg_jacobians-test measures by finite differences, so there is nothing
        left to get backwards.
    */

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
    void set_Racc_matrix(const Matrix3& R_body) { Racc_ = T(0.5) * (R_body + R_body.transpose()); }
    void set_Rmag_std(const Vector3& std_body)  { Rmag_ = std_body.cwiseAbs2().asDiagonal(); }
    void set_Rmag(T variance) { Rmag_ = Matrix3::Identity() * std::abs(variance); }
    void set_RS_noise(T variance) { R_S_ = Matrix3::Identity() * std::abs(variance); }
    void set_RS_noise_vector(const Vector3& variances) { R_S_ = variances.cwiseAbs().asDiagonal(); }
    void set_RS_noise_matrix(const Matrix3& R_S) { R_S_ = T(0.5) * (R_S + R_S.transpose()); }

    // World-frame magnetic reference, NED. Stored as the fixed vector the
    // magnetometer Jacobian is built from -- not as a body-frame prediction.
    void set_magnetic_reference_world(const Vector3& B_w) {
        if (B_w.allFinite() && B_w.norm() > T(1e-9)) {
            B_w_ = B_w;
            have_mag_reference_ = true;
        }
    }
    [[nodiscard]] const Vector3& magnetic_reference_world() const { return B_w_; }
    [[nodiscard]] bool has_magnetic_reference() const { return have_mag_reference_; }

    void set_accel_bias_temp_coeff(const Vector3& k_a) { k_a_ = k_a; }
    void set_accel_bias_limit(T limit) { acc_bias_limit_ = std::abs(limit); }
    void set_acc_bias_updates_enabled(bool on) { acc_bias_updates_enabled_ = on; }

    // Chi-square gate on the 3-DoF NIS. Zero disables it, which is the
    // default and matches OU-III's behaviour of accepting every update.
    void set_meas_gate_nis(T threshold) { nis_gate_ = std::max(T(0), threshold); }

    [[nodiscard]] Vector3 gravity_world() const { return Vector3(T(0), T(0), gravity_magnitude_); }

    /*
        Accelerometer, as a world-frame invariant residual.

            f_m = R^T (a_w - g) + b_a + n_a
            r_a = Rhat (f_m - bhat_a) - (ahat_w - g)

        which is zero at the true state. Linearizing with Rhat = Exp(phi) R
        and ahat_w = rho_a + Exp(phi) a_w:

            r_a ~ [g]x phi - rho_a - Rhat delta_ba + Rhat n_a

        The Jacobian contains the FIXED gravity vector and nothing else --
        no estimated attitude, no estimated wave acceleration. That is the
        structural gain over OU-III, whose J_att = -[f_cog_b]x is rebuilt from
        the current specific-force estimate every step.

        WHAT IS APPROXIMATED, precisely. Differentiating the residual exactly
        gives

            d r_a / d phi = -[Rhat (f_m - bhat_a)]x + [ahat_w]x
                          = [g]x - [r_a]x

        because Rhat(f_m - bhat_a) = ahat_w - g + r_a by definition of r_a. So
        [g]x is the exact derivative only where the residual vanishes, and the
        neglected term is exactly -[r_a]x: first order in the innovation,
        vanishing as the filter converges.

        That term is dropped on purpose. Keeping it would put the measurement
        and the estimated attitude back into the Jacobian and forfeit the one
        property this formulation exists for. It is the standard invariant-EKF
        choice, and tfg_jacobians-test pins its exact size and structure --
        finite differences at a consistent state must match to 1e-8, and the
        discrepancy at an inconsistent state must equal -[r_a]x rather than
        merely being "small".

        Note the accelerometer-bias column: -Rhat, not -I. At Level 1
        delta_ba is an additive BODY-frame error, so it has to be rotated into
        the world frame the residual lives in. It collapses to -I only at
        Level 2, where beta_a = R(bhat_a - b_a) is already world-referred.
        That difference is exactly what the two-frame bias geometry buys, and
        tfg_jacobians-test pins both the value and the reason.
    */
    void accel_residual(const Vector3& acc_meas_body, T tempC,
                        Vector3& r, Eigen::Matrix<T,3,NX>& H, Matrix3& Rw) const {
        const Vector3 ba = accel_bias_at_(tempC);
        r = X_.R * (acc_meas_body - ba) - (X_.X.col(3) - gravity_world());

        H.setZero();
        H.template block<3,3>(0, OFF_PHI) = ocean_imu::lie::skew<T>(gravity_world());
        H.template block<3,3>(0, OFF_AW)  = -Matrix3::Identity();
        if constexpr (with_accel_bias) {
            if (acc_bias_updates_enabled_) {
                H.template block<3,3>(0, OFF_BA) =
                    two_frame_bias ? Matrix3(-Matrix3::Identity()) : Matrix3(-X_.R);
            }
        }
        // Body-frame sensor noise, expressed in the world frame the residual
        // lives in. For isotropic Racc this is a no-op, which is worth knowing
        // when reading the anisotropic case.
        Rw = X_.R * Racc_ * X_.R.transpose();
    }

    /*
        Magnetometer.

            m_m = R^T B_w + n_m
            r_m = Rhat m_m - B_w        ~  -[B_w]x phi + Rhat n_m

        Again the Jacobian is built from the fixed world reference rather than
        from the current predicted body-frame vector.

        Same approximation as the accelerometer, same structure: the exact
        derivative is -[Rhat m_m]x = -[B_w]x - [r_m]x, so the neglected term is
        -[r_m]x, first order in the innovation.
    */
    void mag_residual(const Vector3& mag_meas_body,
                      Vector3& r, Eigen::Matrix<T,3,NX>& H, Matrix3& Rw) const {
        r = X_.R * mag_meas_body - B_w_;
        H.setZero();
        H.template block<3,3>(0, OFF_PHI) = -ocean_imu::lie::skew<T>(B_w_);
        Rw = X_.R * Rmag_ * X_.R.transpose();
    }

    /*
        Integral pseudo-measurement, z_S = 0 = S + n_S.

            r_S = 0 - Shat = -Shat       ~  [Shat]x phi - rho_S + n_S

        The [Shat]x term is O(S) and vanishes at the regularization target, so
        it is routinely dropped. It is kept because it costs one cross product
        and makes the finite-difference check exact rather than approximate --
        a test that has to be loosened to accommodate a term you chose not to
        write stops discriminating against the terms you got wrong.

        R_S stays in world/NED coordinates. A NED-anisotropic R_S is a
        deliberate physical choice (heave regularized differently from surge
        and sway), and the right-invariant error is what preserves its meaning
        -- but it does break rotational symmetry, so the filter must not then
        be described as invariant under arbitrary world rotations.
    */
    void integral_residual(Vector3& r, Eigen::Matrix<T,3,NX>& H, Matrix3& Rw) const {
        r = -X_.X.col(2);
        H.setZero();
        H.template block<3,3>(0, OFF_PHI) = ocean_imu::lie::skew<T>(Vector3(X_.X.col(2)));
        H.template block<3,3>(0, OFF_S)   = -Matrix3::Identity();
        Rw = R_S_;
    }

    bool measurement_update_acc_only(const Vector3& acc_meas_body, T tempC = T(35)) {
        if (!acc_meas_body.allFinite()) { last_acc_diag_ = MeasDiag3{}; return false; }
        Vector3 r; Eigen::Matrix<T,3,NX> H; Matrix3 Rw;
        accel_residual(acc_meas_body, tempC, r, H, Rw);
        const bool ok = apply_update3_(r, H, Rw, last_acc_diag_);
        if (ok) project_acc_bias_();
        return ok;
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
        Vector3 r; Eigen::Matrix<T,3,NX> H; Matrix3 Rw;
        integral_residual(r, H, Rw);
        return apply_update3_(r, H, Rw, last_S_diag_);
    }

    // Exposed so the covariance-transport test can drive one update in
    // isolation and inspect what it did.
    bool apply_update3(const Vector3& r, const Eigen::Matrix<T,3,NX>& H,
                       const Matrix3& Rw, MeasDiag3& diag) {
        return apply_update3_(r, H, Rw, diag);
    }

  private:
    [[nodiscard]] Vector3 accel_bias_at_(T tempC) const {
        if constexpr (with_accel_bias) {
            return X_.B.col(1) + k_a_ * (tempC - kTempRefC);
        } else {
            (void)tempC;
            return Vector3::Zero();
        }
    }

    /*
        Rank-3 Joseph update in the invariant tangent frame.

        There is deliberately NO MEKF-style covariance reset afterwards.
        ou_detail::apply_left_error_reset's G = I + 0.5 [dtheta]x exists because
        OU-III's error state is zeroed after injection and its covariance has
        to be transported into the new linearization point. Here the update is
        already formulated in the post-update tangent frame, so applying that
        transport on top would double-correct the covariance. The
        no-double-reset assertion lives in covariance_transport-test.
    */
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

        const MatrixNX IKH = MatrixNX::Identity() - K * H;
        P_ = (IKH * P_ * IKH.transpose() + K * Rw * K.transpose()).eval();
        P_ = T(0.5) * (P_ + P_.transpose()).eval();

        // xi = K r is the estimated error; remove it.
        inject(Tangent(-(K * r)));

        diag.accepted = true;
        return true;
    }

    /*
        Confine b_a to a ball, as OU-III does. The rationale there is the ISS
        argument in doc/kalman_ou_iii/w3d-iss-stability.tex-part: tilt,
        horizontal acceleration and accelerometer bias are only weakly
        separable, so an unbounded bias state will absorb a persistent tilt
        error and walk away. A Lie-group formulation does not change that --
        it is a physical observability limit, not a coordinate artifact.
    */
    void project_acc_bias_() {
        if constexpr (with_accel_bias) {
            if (!(acc_bias_limit_ > T(0))) return;
            const T n = X_.B.col(1).norm();
            if (n > acc_bias_limit_) X_.B.col(1) *= (acc_bias_limit_ / n);
        }
    }

    void reorthonormalize_() {
        // Cheaper and better conditioned than a QR: round-trip through the
        // quaternion, which is what the surrounding filters store anyway.
        Eigen::Quaternion<T> q(X_.R);
        q.normalize();
        X_.R = q.toRotationMatrix();
    }

    // Exact OU process noise for [v p S a_w], reusing OU-III's assembly
    // including the correlated-axis Kronecker path.
    void add_ou_process_noise_(T Ts, MatrixNX& Qd) const {
        const bool correlated = !ou_detail::is_isotropic3<T>(Sigma_aw_) &&
                                !is_diagonal_(Sigma_aw_);

        if (correlated) {
            Eigen::Matrix<T,4,4> Q_unit;
            ou_detail::IntegratedOUChain<T,3>::process_covariance(tau_aw_, Ts, T(1), Q_unit);
            Q_unit = T(0.5) * (Q_unit + Q_unit.transpose()).eval();
            const Matrix3 Sig = T(0.5) * (Sigma_aw_ + Sigma_aw_.transpose());
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    Qd.template block<3,3>(OFF_V + 3*i, OFF_V + 3*j).noalias()
                        += Sig * Q_unit(i,j);
                }
            }
        } else {
            for (int axis = 0; axis < 3; ++axis) {
                Eigen::Matrix<T,4,4> Q_axis;
                ou_detail::IntegratedOUChain<T,3>::process_covariance(
                    tau_aw_, Ts, Sigma_aw_(axis,axis), Q_axis);
                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        Qd(OFF_V + 3*i + axis, OFF_V + 3*j + axis) += Q_axis(i,j);
                    }
                }
            }
        }
    }

    static bool is_diagonal_(const Matrix3& M) {
        Matrix3 off = M;
        off.diagonal().setZero();
        return off.cwiseAbs().maxCoeff() <= T(1e-12) * std::max(T(1), M.cwiseAbs().maxCoeff());
    }

    T gravity_magnitude_;

    Group   X_{};
    MatrixNX P_{MatrixNX::Identity()};

    T       tau_aw_{T(1.0)};
    Matrix3 Sigma_aw_{Matrix3::Identity()};
    Matrix3 Q_gyro_{Matrix3::Identity() * T(2.5e-5)};
    Matrix3 Q_bg_{Matrix3::Identity() * T(1e-10)};
    Matrix3 Q_ba_{Matrix3::Identity() * T(2.5e-7)};

    // Measurement noise. Racc_ and Rmag_ are BODY frame and get rotated into
    // the world frame per update; R_S_ is already world/NED.
    Matrix3 Racc_{Matrix3::Identity() * T(0.16)};
    Matrix3 Rmag_{Matrix3::Identity() * T(0.01)};
    Matrix3 R_S_{Matrix3::Identity() * T(1.5)};

    Vector3 B_w_{Vector3::UnitX()};
    bool    have_mag_reference_{false};

    Vector3 k_a_{Vector3::Zero()};
    T       acc_bias_limit_{T(0.5)};
    bool    acc_bias_updates_enabled_{true};
    T       nis_gate_{T(0)};

    static constexpr T kTempRefC = T(35);

    MeasDiag3 last_acc_diag_{};
    MeasDiag3 last_mag_diag_{};
    MeasDiag3 last_S_diag_{};

    Vector3 last_omega_{Vector3::Zero()};
};

} // namespace ocean_imu::kalman
