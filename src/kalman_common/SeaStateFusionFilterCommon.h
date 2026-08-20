#pragma once

#include <cmath>
#include <algorithm>

#include "tuner/SeaStateFusionTunerCommon.h"
#include "tuner/SeaStateAdaptationLimits.h"

namespace seastate::common {

using FreqInputLPF = seastate::tuner::common::FreqInputLPF;
using StillnessAdapter = seastate::tuner::common::StillnessAdapter;

struct StartupTiltObserver {
    Eigen::Vector3f s_body = Eigen::Vector3f(0.0f, 0.0f, -1.0f); // predicted specific-force dir at rest
    bool initialized = false;

    void reset() {
        s_body = Eigen::Vector3f(0.0f, 0.0f, -1.0f);
        initialized = false;
    }

    Eigen::Vector3f step(const Eigen::Vector3f& gyro_body_ned,
                         const Eigen::Vector3f& acc_body_ned,
                         float dt,
                         float acc_tau_sec,
                         float g_ref,
                         float norm_frac)
    {
        if (!(dt > 0.0f) || !std::isfinite(dt)) {
            return s_body;
        }

        // Initialize directly from first valid accel direction.
        if (!initialized) {
            if (acc_body_ned.allFinite()) {
                const float an = acc_body_ned.norm();
                if (std::isfinite(an) && an > 1e-6f) {
                    s_body = acc_body_ned / an;
                    initialized = true;
                }
            }
            return s_body;
        }

        // Propagate with gyro: s_dot = -omega x s
        Eigen::Vector3f s_pred = s_body;
        if (gyro_body_ned.allFinite()) {
            s_pred += dt * (-gyro_body_ned.cross(s_pred));
            const float sn = s_pred.norm();
            if (std::isfinite(sn) && sn > 1e-6f) {
                s_pred /= sn;
            } else {
                s_pred = s_body;
            }
        }

        // Correct slowly toward measured accel direction, with norm-based weighting.
        if (acc_body_ned.allFinite()) {
            const float an = acc_body_ned.norm();
            if (std::isfinite(an) && an > 1e-6f) {
                const Eigen::Vector3f u_a = acc_body_ned / an;

                const float rel_err =
                    std::fabs(an - g_ref) / std::max(g_ref, 1e-6f);

                float w = 1.0f - rel_err / std::max(norm_frac, 1e-3f);
                w = std::min(std::max(w, 0.0f), 1.0f);

                const float tau = std::max(acc_tau_sec, 1.0e-3f);
                const float alpha = w * (1.0f - std::exp(-dt / tau));

                Eigen::Vector3f s_upd = (1.0f - alpha) * s_pred + alpha * u_a;
                const float un = s_upd.norm();
                if (std::isfinite(un) && un > 1e-6f) {
                    s_body = s_upd / un;
                } else {
                    s_body = s_pred;
                }
                return s_body;
            }
        }

        s_body = s_pred;
        return s_body;
    }
};

inline float unitVecAlignSin(const Eigen::Vector3f& a,
                             const Eigen::Vector3f& b)
{
    if (!a.allFinite() || !b.allFinite()) return 1.0f;

    const float an = a.norm();
    const float bn = b.norm();
    if (!(an > 1e-6f) || !(bn > 1e-6f) ||
        !std::isfinite(an) || !std::isfinite(bn))
    {
        return 1.0f;
    }

    const Eigen::Vector3f ua = a / an;
    const Eigen::Vector3f ub = b / bn;

    const float s = ua.cross(ub).norm();
    if (!std::isfinite(s)) return 1.0f;

    return std::min(std::max(s, 0.0f), 1.0f);
}

// Companion to unitVecAlignSin.  A sine residual is symmetric about
// 90 deg, so it reads the same at an angle and at its supplement: it cannot
// tell the aligned branch from the antipodal one.  This returns the cosine,
// whose sign is exactly that branch, and fails closed at -1 so a degenerate
// input can never be read as aligned.
inline float unitVecAlignCos(const Eigen::Vector3f& a,
                             const Eigen::Vector3f& b)
{
    if (!a.allFinite() || !b.allFinite()) return -1.0f;

    const float an = a.norm();
    const float bn = b.norm();
    if (!(an > 1e-6f) || !(bn > 1e-6f) ||
        !std::isfinite(an) || !std::isfinite(bn))
    {
        return -1.0f;
    }

    const float c = (a / an).dot(b / bn);
    if (!std::isfinite(c)) return -1.0f;

    return std::min(std::max(c, -1.0f), 1.0f);
}

inline Eigen::Vector3f predictedGravityDirBody(const Eigen::Quaternionf& q_bw_in)
{
    if (!q_bw_in.coeffs().allFinite()) {
        return Eigen::Vector3f::Zero();
    }

    Eigen::Quaternionf q_bw = q_bw_in;
    q_bw.normalize();

    const Eigen::Vector3f world_down(0.0f, 0.0f, 1.0f);

    // Predicted specific-force direction at rest in body frame:
    // acc_body ~= -(world_down expressed in body)
    Eigen::Vector3f u_g_body = -(q_bw.conjugate() * world_down);

    const float n = u_g_body.norm();
    if (!(n > 1e-6f) || !u_g_body.allFinite()) {
        return Eigen::Vector3f::Zero();
    }

    return u_g_body / n;
}

inline float gravityAlignResidualSin(const Eigen::Quaternionf& q_bw_in,
                                     const Eigen::Vector3f& acc_body_ned)
{
    if (!acc_body_ned.allFinite()) return 1.0f;

    const float an = acc_body_ned.norm();
    if (!(an > 1e-6f) || !std::isfinite(an)) return 1.0f;

    const Eigen::Vector3f u_a = acc_body_ned / an;
    const Eigen::Vector3f u_g = predictedGravityDirBody(q_bw_in);

    const float gn = u_g.norm();
    if (!(gn > 1e-6f) || !u_g.allFinite()) return 1.0f;

    const Eigen::Vector3f r = u_a.cross(u_g);
    const float s = r.norm();
    if (!std::isfinite(s)) return 1.0f;

    return std::min(std::max(s, 0.0f), 1.0f);
}

// Branch companion to gravityAlignResidualSin: the cosine between the
// predicted and measured rest-specific-force directions.  Positive means the
// attitude is on the aligned branch; the sine residual alone is satisfied just
// as well by the attitude flipped through 180 deg.
inline float gravityAlignResidualCos(const Eigen::Quaternionf& q_bw_in,
                                     const Eigen::Vector3f& acc_body_ned)
{
    return unitVecAlignCos(predictedGravityDirBody(q_bw_in), acc_body_ned);
}

// The runtime branch certificate the startup gates are judged against:
// s_hat^T s_meas > 0.  Kept as a named predicate because it is a proof
// condition, not a tuning knob, and is therefore not exposed as a threshold.
inline bool gravityAlignedBranch(const Eigen::Quaternionf& q_bw_in,
                                 const Eigen::Vector3f& acc_body_ned)
{
    return gravityAlignResidualCos(q_bw_in, acc_body_ned) > 0.0f;
}

template<typename Vec3LPFType, typename OnReadyFn>
inline bool runStartupGravityInit(const Eigen::Vector3f& gyro_body_ned,
                                  const Eigen::Vector3f& acc_body_ned,
                                  float dt,
                                  float elapsed_sec,
                                  float gravity_std,
                                  float obs_acc_tau_sec,
                                  float gravity_slow_tau_sec,
                                  float gravity_align_max_sin,
                                  float gravity_hold_sec,
                                  float gravity_min_sec,
                                  float gravity_timeout_sec,
                                  float gravity_norm_frac,
                                  StartupTiltObserver& bootstrap_tilt_obs,
                                  Vec3LPFType& bootstrap_gravity_slow_lpf,
                                  float& bootstrap_gravity_good_sec,
                                  OnReadyFn&& on_ready)
{
    if (!acc_body_ned.allFinite()) return false;

    const Eigen::Vector3f s_obs =
        bootstrap_tilt_obs.step(
            gyro_body_ned,
            acc_body_ned,
            dt,
            obs_acc_tau_sec,
            gravity_std,
            gravity_norm_frac);

    const Eigen::Vector3f g_slow =
        bootstrap_gravity_slow_lpf.step(
            acc_body_ned, dt, gravity_slow_tau_sec);

    const float align_sin = unitVecAlignSin(s_obs, g_slow);

    // Branch certificate.  The sine residual above is blind to a 180 deg flip,
    // so on its own it would accept an observer that has settled on the
    // antipodal gravity direction and seed the filter upside down.  The
    // aligned branch is the sign of the dot product, and it is required both
    // for the quality path and for the timeout path below.
    const bool aligned_branch = (unitVecAlignCos(s_obs, g_slow) > 0.0f);

    const float g_slow_n = g_slow.norm();
    const float g_rel_err =
        (std::isfinite(g_slow_n) && gravity_std > 1e-6f)
            ? std::fabs(g_slow_n - gravity_std) / gravity_std
            : INFINITY;

    const bool gravity_good_now =
        std::isfinite(align_sin) &&
        (align_sin <= gravity_align_max_sin) &&
        aligned_branch &&
        std::isfinite(g_rel_err) &&
        (g_rel_err <= gravity_norm_frac);

    if (gravity_good_now) {
        bootstrap_gravity_good_sec += dt;
        if (bootstrap_gravity_good_sec > 10.0f) {
            bootstrap_gravity_good_sec = 10.0f;
        }
    } else {
        bootstrap_gravity_good_sec =
            std::max(0.0f, bootstrap_gravity_good_sec - 2.0f * dt);
    }

    const bool ready_by_quality =
        (elapsed_sec >= gravity_min_sec) &&
        (bootstrap_gravity_good_sec >= gravity_hold_sec);

    const bool ready_by_timeout =
        (elapsed_sec >= gravity_timeout_sec);

    if (!ready_by_quality && !ready_by_timeout) return false;

    // Both acceptance paths are held to the branch certificate on the sample
    // that actually seeds the filter.  The hold counter alone is not that
    // certificate: it can still stand above the hold time for one sample after
    // the branch flips.  A forced handoff is still a handoff -- the timeout may
    // bound how long startup takes, but it may not hand over an attitude that
    // is upside down with respect to the measurement it is derived from.  The
    // antipodal set is not attracting for the accel-corrected observer, so this
    // is a bounded extra wait rather than a stall.
    if (!aligned_branch) return false;

    Eigen::Vector3f g_init_dir = s_obs;

    if (std::isfinite(g_slow_n) && g_slow_n > 1e-6f) {
        const Eigen::Vector3f g_slow_u = g_slow / g_slow_n;
        Eigen::Vector3f blend = 0.90f * s_obs + 0.10f * g_slow_u;
        const float bn = blend.norm();
        if (std::isfinite(bn) && bn > 1e-6f) {
            g_init_dir = blend / bn;
        }
    }

    on_ready(gravity_std * g_init_dir);
    return true;
}

template<typename DetrenderConfig>
inline DetrenderConfig defaultDisplacementDetrenderConfig(float freq_guess_hz) {
    return seastate::tuner::common::defaultDisplacementDetrenderConfig<DetrenderConfig>(freq_guess_hz);
}

// Smoothing horizon for a drift-correction channel whose target is a power law
// in the OU operating point (r_S ~ sigma_aw tau^3, r_p0 ~ sigma_aw tau^2,
// r_v0 ~ sigma_aw tau).  Those targets are rebuilt every step from the raw
// tuner estimates, so the channel needs an exponential smoother that its
// inputs do not already provide.
//
// The horizon is `mult * tau`, i.e. a fixed number of wave periods rather than
// a fixed number of seconds, so it follows the sea state.  A single horizon has
// to serve two regimes that want opposite things: while the sea is stationary
// the target only jitters and a long horizon rejects that jitter, but while the
// sea state moves the same long horizon lags the target and the filter runs on
// a stale operating point.
//
// `slew_log` resolves that.  It is the size, in natural-log units of the
// target-to-applied ratio, of a discrepancy considered "real" rather than
// jitter.  The horizon is divided by 1 + (discrepancy / slew_log)^2, so jitter
// well inside the threshold leaves the long horizon intact while a sustained
// move shortens it quadratically.  slew_log <= 0 disables the term and leaves
// the plain proportional horizon.
inline float adaptiveSmoothingHorizonSec(float mult,
                                         float tau_sec,
                                         float target,
                                         float applied,
                                         float slew_log,
                                         float dt)
{
    const float safe_time_scale =
        seastate::tuner::limits::clampDynamicEmaTimeScaleSec(tau_sec);
    float horizon = mult * safe_time_scale;

    if (slew_log > 0.0f && target > 0.0f && applied > 0.0f &&
        std::isfinite(target) && std::isfinite(applied))
    {
        const float discrepancy = std::fabs(std::log(target / applied));
        if (std::isfinite(discrepancy)) {
            const float k = discrepancy / slew_log;
            horizon /= (1.0f + k * k);
        }
    }

    // The dynamic time base and the final EMA horizon have independent
    // guards: the latter also protects the optional discrepancy gate above
    // from collapsing a legitimate sea-scaled horizon into sample-to-sample
    // pass-through.
    return seastate::tuner::limits::clampDynamicEmaHorizonSec(horizon, dt);
}

template<typename MekfPtr, typename EnterColdFn, typename ApplyTuneFn>
inline void finalizeInitialization(MekfPtr& mekf, EnterColdFn&& enterCold, ApplyTuneFn&& applyTune) {
    enterCold();
    applyTune();
    // Initial tuning changes the stationary OU model before the linear block is
    // activated. Seed the matching posterior covariance once at initialization.
    mekf->reset_aw_covariance_to_stationary();
    mekf->set_exact_att_bias_Qd(true);
}

}  // namespace seastate::common
