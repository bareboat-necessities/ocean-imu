#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

// Startup-attitude and magnetic-acquisition primitives shared by the OU-II,
// OU-III and TFG orchestrators.
//
// Everything here is measurement-side plumbing: yaw bookkeeping on BODY->NED
// quaternions, the world-frame gravity gate that certifies a startup tilt, the
// MagAutoTuner / continuous hard-iron configuration, and the slew that moves
// an applied hard-iron offset toward the continuous estimate.  None of it reads
// or writes an estimator state; each wrapper decides *when* these run and what
// the results are written into, and that part stays in the wrapper because it
// differs (see the TFG north-frame acquisition, for example).

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>

#include "kalman_common/SeaStateFusionFilterCommon.h"
#include "tuner/ContinuousMagHardIronEstimator.h"
#include "tuner/MagAutoTuner.h"

namespace seastate::common {

inline float wrapPi(float a) {
    constexpr float PI_F = 3.14159265358979323846f;
    constexpr float TWO_PI_F = 2.0f * PI_F;
    if (!std::isfinite(a)) return NAN;
    while (a > PI_F) a -= TWO_PI_F;
    while (a <= -PI_F) a += TWO_PI_F;
    return a;
}

// Heading of a BODY->WORLD quaternion, atan2(R10, R00).
inline float yawFromBoatQuatRad(const Eigen::Quaternionf& q_bw_in) {
    if (!q_bw_in.coeffs().allFinite()) return NAN;

    Eigen::Quaternionf q_bw = q_bw_in;
    const float qn = q_bw.norm();
    if (!(qn > 1.0e-6f) || !std::isfinite(qn)) return NAN;
    q_bw.normalize();

    const Eigen::Matrix3f R = q_bw.toRotationMatrix();
    const float c = R(0, 0);
    const float s = R(1, 0);
    if (!std::isfinite(c) || !std::isfinite(s)) return NAN;
    return std::atan2(s, c);
}

// Tilt part of a BODY->WORLD quaternion, with the heading divided out.
// Invariant under q_bw -> Rz(psi) q_bw, which is what makes it usable as a
// magnetometer accumulation frame that cannot leak an arbitrary estimator yaw.
inline Eigen::Quaternionf yawRemovedBoatQuat(const Eigen::Quaternionf& q_bw_in) {
    if (!q_bw_in.coeffs().allFinite()) return Eigen::Quaternionf::Identity();

    Eigen::Quaternionf q_bw = q_bw_in;
    const float qn = q_bw.norm();
    if (!(qn > 1.0e-6f) || !std::isfinite(qn)) return Eigen::Quaternionf::Identity();
    q_bw.normalize();

    const float yaw = yawFromBoatQuatRad(q_bw);
    if (!std::isfinite(yaw)) return Eigen::Quaternionf::Identity();

    const Eigen::Quaternionf q_yaw_inv(Eigen::AngleAxisf(-yaw, Eigen::Vector3f::UnitZ()));
    Eigen::Quaternionf q_tilt = q_yaw_inv * q_bw;
    q_tilt.normalize();
    if (!q_tilt.coeffs().allFinite()) return Eigen::Quaternionf::Identity();
    return q_tilt;
}

// Rewrites heading only, keeping the tilt untouched.  The incoming heading is
// stripped first, so a drifted yaw cannot survive the write.
inline Eigen::Quaternionf boatQuatWithAbsoluteYaw(const Eigen::Quaternionf& q_bw_in,
                                                  float yaw_abs_rad) {
    if (!std::isfinite(yaw_abs_rad)) return q_bw_in;

    const Eigen::Quaternionf q_tilt = yawRemovedBoatQuat(q_bw_in);
    const Eigen::Quaternionf q_yaw(Eigen::AngleAxisf(yaw_abs_rad, Eigen::Vector3f::UnitZ()));
    Eigen::Quaternionf q_out = q_yaw * q_tilt;
    q_out.normalize();
    if (!q_out.coeffs().allFinite()) return q_bw_in;
    return q_out;
}

// First-order vector low-pass, seeded by its first finite input.
struct Vec3LPF {
    Eigen::Vector3f state = Eigen::Vector3f::Zero();
    bool initialized = false;

    void reset() {
        state.setZero();
        initialized = false;
    }

    Eigen::Vector3f step(const Eigen::Vector3f& x, float dt, float tau_sec) {
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

// Advance a sample clock and return the interval since its last sample, or
// the nominal interval on the first sample and after a non-advancing clock.
inline float advanceSampleClock(float& last_t, float t, float nominal_dt) {
    const float dt = (std::isfinite(last_t) && t > last_t) ? (t - last_t) : nominal_dt;
    last_t = t;
    return dt;
}

// World-frame gravity-agreement gate for a startup attitude.
//
// The residual is taken in the attitude's own world frame rather than in the
// body frame.  A body-frame average of the specific force is not gravity under
// way: the hull rolls and pitches through the window, so the orbital term the
// average is there to remove is smeared across it instead of cancelling.  What
// a body-frame gate then reports can be dominated by the sea state rather than
// the levelling error.  Rotating before averaging keeps the residual in a
// fixed frame.  See gravityAlignResidualSinWorld().
//
// The warmup exists because the average and the observer are seeded from the
// *same* accelerometer sample.  Until the average has moved off that seed, a
// small residual only says the two agree about the instant they both started
// from, which they do by construction even when the boat was mid-wave.
//
// The branch (sign of the world down component) is deliberately not held
// behind the warmup.  The sine residual is the same at an angle and at its
// supplement, so it accepts an attitude flipped through 180 deg as readily as
// the right one; the branch is the one part of the certificate an unaveraged
// sample can answer, and handoff timeouts are gated on it, so withholding it
// early would let a stalled startup sit unbranched rather than fail closed.
struct WorldGravityGate {
    Vec3LPF acc_world_lpf{};
    float   elapsed_sec    = 0.0f;
    float   good_sec       = 0.0f;
    bool    aligned_branch = false;

    void reset() {
        acc_world_lpf.reset();
        elapsed_sec = 0.0f;
        good_sec = 0.0f;
        aligned_branch = false;
    }

    // q_bw is the attitude being judged; only its tilt matters (see
    // accWorldFromBody).  The gyro only vetoes truly violent motion.
    void step(const Eigen::Quaternionf& q_bw,
              const Eigen::Vector3f& acc_body_ned,
              const Eigen::Vector3f& gyro_body_ned,
              float dt,
              float lpf_tau_sec,
              float warmup_sec,
              float max_sin,
              float extreme_gyro_dps)
    {
        acc_world_lpf.step(accWorldFromBody(q_bw, acc_body_ned), dt, lpf_tau_sec);
        const Eigen::Vector3f acc_world_lp = acc_world_lpf.state;
        elapsed_sec += dt;

        const bool average_warm = elapsed_sec >= warmup_sec;
        const float align_sin = average_warm
            ? gravityAlignResidualSinWorld(acc_world_lp)
            : 1.0f;
        aligned_branch = gravityAlignedBranchWorld(acc_world_lp);

        const float gyro_dps = gyro_body_ned.norm() * 57.295779513f;
        const bool extreme_motion =
            !std::isfinite(gyro_dps) || (gyro_dps > extreme_gyro_dps);

        const bool good_now =
            std::isfinite(align_sin) &&
            (align_sin <= max_sin) &&
            aligned_branch &&
            !extreme_motion;

        if (good_now) {
            good_sec += dt;
            if (good_sec > 10.0f) good_sec = 10.0f;
        } else {
            good_sec = std::max(0.0f, good_sec - 2.0f * dt);
        }
    }

    // Held agreement on the aligned branch.
    bool trusted(float hold_sec) const {
        return aligned_branch && (good_sec >= hold_sec);
    }
};

// MagAutoTuner configuration from a wrapper Config.  All three wrappers name
// these fields identically.
template <typename Config>
inline MagAutoTuner::Config magAutoTunerConfig(const Config& cfg) {
    MagAutoTuner::Config mag_cfg;
    mag_cfg.mag_norm_min             = cfg.mag_init_min_mag_norm;
    mag_cfg.min_samples              = cfg.mag_min_samples;
    mag_cfg.min_window_sec           = cfg.mag_min_window_sec;
    mag_cfg.max_window_sec           = cfg.mag_max_window_sec;
    mag_cfg.sample_dt_sec            = cfg.mag_sample_dt_sec;
    mag_cfg.gravity_ref              = cfg.gravity_magnitude;
    mag_cfg.enable_quality_weighting = cfg.mag_enable_quality_weighting;
    mag_cfg.estimate_hard_iron       = cfg.mag_estimate_hard_iron;
    mag_cfg.min_effective_weight     = cfg.mag_min_effective_weight;
    mag_cfg.acc_norm_rel_soft        = cfg.mag_acc_norm_rel_soft;
    mag_cfg.gyro_soft_dps            = cfg.mag_gyro_soft_dps;
    return mag_cfg;
}

// Second-stage acquisition: the same tuner, re-armed with the refinement
// window.
inline void beginMagRefinement(MagAutoTuner& tuner, float window_sec, int min_samples) {
    MagAutoTuner::Config refine_cfg = tuner.config();
    refine_cfg.min_window_sec = window_sec;
    refine_cfg.min_samples    = min_samples;
    tuner.setConfig(refine_cfg);
    tuner.reset();
}

// Continuous hard-iron estimation, and the reference that goes with it.
//
// The startup estimate is a single window, and a single window is where the
// offset is least identifiable: the excitation is whatever tilt the hull
// happened to take in fifteen seconds.  The continuous estimator never closes
// its accumulation, so the question stops being "can the offset be read out of
// these fifteen seconds" and becomes "keep watching, and correct when the data
// finally say something".
//
// Two things make that safe to leave running.  The estimator is exogenous --
// gravity-referenced tilt from the private Mahony observer and the raw
// magnetometer, never a filter state, so no loop is closed through the
// estimator.  And the applied offset and the magnetic reference move together,
// out of the same statistics, so the filter is never subtracting one offset
// while steering to a reference that belongs to another.
struct ContinuousHardIronTracker {
    ContinuousMagHardIronEstimator estimator{};

    Eigen::Vector3f startup_body_uT     = Eigen::Vector3f::Zero();
    Eigen::Vector3f applied_body_uT     = Eigen::Vector3f::Zero();
    Eigen::Vector3f anchor_bias_body_uT = Eigen::Vector3f::Zero();
    Eigen::Vector3f anchor_world_ref_uT = Eigen::Vector3f::Zero();
    bool  anchored      = false;
    float last_sample_t = NAN;
    float last_apply_t  = NAN;

    template <typename Config>
    void configure(const Config& cfg) {
        ContinuousMagHardIronEstimator::Config hi_cfg;
        hi_cfg.memory_sec           = cfg.mag_hi_memory_sec;
        hi_cfg.model_ridge          = cfg.mag_hi_model_ridge;
        hi_cfg.model_ridge_relative = cfg.mag_hi_model_ridge_relative;
        hi_cfg.min_information      = cfg.mag_hi_min_information;
        hi_cfg.min_effective_weight = cfg.mag_hi_min_effective_weight;
        hi_cfg.max_residual_rms_uT  = cfg.mag_hi_max_residual_rms_uT;
        hi_cfg.max_bias_fraction    = cfg.mag_hi_max_bias_fraction;
        hi_cfg.min_mag_norm_uT      = cfg.mag_init_min_mag_norm;
        estimator.setConfig(hi_cfg);
    }

    // Clears the applied offset and its anchors; the estimator's statistics
    // are left to configure().
    void resetApplied() {
        startup_body_uT.setZero();
        applied_body_uT.setZero();
        anchor_bias_body_uT.setZero();
        anchor_world_ref_uT.setZero();
        anchored = false;
        last_sample_t = NAN;
        last_apply_t  = NAN;
    }

    // Feed the exogenous accumulation.  Raw magnetometer -- not the corrected
    // stream -- because the estimator is fitting the offset itself and must
    // not be shown data with its own answer already subtracted.
    void accumulate(float t, float nominal_dt,
                    const Eigen::Quaternionf& q_tilt_bw,
                    const Eigen::Vector3f& mag_body_uT) {
        const float dt_mag = advanceSampleClock(last_sample_t, t, nominal_dt);
        estimator.update(dt_mag, q_tilt_bw, mag_body_uT);
    }

    // Move the applied offset toward the fit, and re-gauge the reference.
    //
    // The reference has to be rebuilt in MagAutoTuner's canonical form --
    // horizontal magnitude on +X, vertical below it -- and not merely shifted
    // by the same amount as the measurement.  A shift that tracks the offset
    // exactly is a no-op: subtracting b from every sample and subtracting the
    // matching mean(R) b from the reference leaves the innovation identical at
    // the attitude the filter already holds, so nothing moves and the standing
    // yaw error survives the correction that was meant to remove it.
    //
    // The standing error is a *gauge*: the startup acquisition put the world
    // frame's north along the average of the uncorrected field, which is
    // magnetic north rotated by whatever the offset contributes.  Leaving the
    // canonical reference in place while the offset comes out of the stream
    // asks the filter for the heading the corrected field implies, and the
    // magnetometer update walks the yaw there over its own time constant.  No
    // attitude state is written, so the correction remains a change of
    // measurement-model parameters.
    //
    // Only the horizontal magnitude and the vertical component move with the
    // offset, and only by the amount the offset changes them.  They are not
    // recomputed from the estimator's own window: that window is longer and
    // less selective than the one the startup acquisition gated, and simply
    // adopting its magnitude and dip costs a fifth of the roll accuracy on the
    // moderate seas while the offset correction itself costs none of it.
    //
    // Returns true with the new canonical reference when the applied offset
    // moved; applied_body_uT is updated only then.  The caller owns the
    // eligibility gates and writes the reference and the total offset.
    bool slewTowardEstimate(float t, float nominal_dt,
                            const Eigen::Vector3f& current_world_ref_uT,
                            float apply_fraction,
                            float slew_tau_sec,
                            float min_mag_norm,
                            Eigen::Vector3f& new_world_ref_uT)
    {
        const auto& est = estimator.estimate();
        if (!est.valid) return false;

        if (!anchored) {
            anchor_bias_body_uT = applied_body_uT;
            anchor_world_ref_uT = current_world_ref_uT;
            anchored = true;
        }

        const Eigen::Vector3f target = apply_fraction * est.bias_body_uT;
        if (!target.allFinite()) return false;

        const float dt_apply = advanceSampleClock(last_apply_t, t, nominal_dt);

        const float tau = slew_tau_sec;
        const float alpha = (std::isfinite(tau) && tau > 1.0e-3f)
                                ? (1.0f - std::exp(-dt_apply / tau))
                                : 1.0f;

        const Eigen::Vector3f applied =
            applied_body_uT + alpha * (target - applied_body_uT);
        if (!applied.allFinite()) return false;

        // Both evaluated against the statistics as they stand now, so the
        // difference is the offset's doing and nothing else.
        Eigen::Vector3f level_new;
        Eigen::Vector3f level_anchor;
        if (!estimator.levelReferenceForBias(applied, level_new) ||
            !estimator.levelReferenceForBias(anchor_bias_body_uT, level_anchor)) {
            return false;
        }

        const float h_new = level_new.head<2>().norm();
        const float h_anchor = level_anchor.head<2>().norm();
        if (!std::isfinite(h_new) || !std::isfinite(h_anchor)) return false;

        const float h = anchor_world_ref_uT.x() + (h_new - h_anchor);
        const float z = anchor_world_ref_uT.z() + (level_new.z() - level_anchor.z());
        if (!(h > min_mag_norm) || !std::isfinite(h) || !std::isfinite(z)) return false;

        const Eigen::Vector3f ref(h, 0.0f, z);
        if (!ref.allFinite()) return false;

        applied_body_uT = applied;
        new_world_ref_uT = ref;
        return true;
    }
};

}  // namespace seastate::common
