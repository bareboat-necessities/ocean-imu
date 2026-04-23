#pragma once

#include <cmath>
#include <algorithm>

namespace seastate::common {

struct FreqInputLPF {
    float state       = 0.0f;
    float fc_hz       = 1.0f;
    bool  initialized = false;

    void setCutoff(float fc) {
        if (std::isfinite(fc) && fc > 0.0f) {
            fc_hz = fc;
        }
    }

    float step(float x, float dt) {
        const float alpha = std::exp(-2.0f * static_cast<float>(M_PI) * fc_hz * dt);

        if (!initialized) {
            state       = x;
            initialized = true;
            return state;
        }

        state = (1.0f - alpha) * x + alpha * state;
        return state;
    }
};

struct StillnessAdapter {
    explicit StillnessAdapter(float gravity_std = 9.80665f,
                              float target_freq_hz_in = 0.2f,
                              float freq_guess = 0.3f)
        : gravity_std_(gravity_std),
          target_freq_hz(target_freq_hz_in),
          freq_state(freq_guess) {}

    float energy_ema      = 0.0f;
    float energy_alpha    = 0.05f;
    float energy_thresh   = 8e-4f;

    float still_time_sec  = 0.0f;
    float still_thresh_s  = 2.0f;

    float relax_tau_sec   = 1.0f;
    float target_freq_hz;

    bool  freq_init       = false;
    float freq_state;

    bool  last_is_still   = false;

    void setTargetFreqHz(float f) {
        if (std::isfinite(f) && f > 0.0f) {
            target_freq_hz = f;
        }
    }

    float step(float a_z_body_proxy_lp, float dt, float freq_in) {
        if (!(dt > 0.0f) || !std::isfinite(freq_in)) {
            return freq_in;
        }

        if (!freq_init || !std::isfinite(freq_state)) {
            freq_state = freq_in;
            freq_init  = true;
        }

        const float a_norm      = a_z_body_proxy_lp / gravity_std_;
        const float inst_energy = a_norm * a_norm;

        energy_ema = (1.0f - energy_alpha) * energy_ema + energy_alpha * inst_energy;
        const bool is_still = (energy_ema < energy_thresh);
        last_is_still = is_still;

        if (is_still) {
            still_time_sec += dt;
            if (still_time_sec > 60.0f) {
                still_time_sec = 60.0f;
            }

            if (still_time_sec > still_thresh_s) {
                const float relax_alpha = 1.0f - std::exp(-dt / relax_tau_sec);
                freq_state += relax_alpha * (target_freq_hz - freq_state);
            } else {
                freq_state = freq_in;
            }
        } else {
            still_time_sec = 0.0f;
            freq_state     = freq_in;
        }

        return freq_state;
    }

    bool  isStill()      const { return last_is_still; }
    float getStillTime() const { return still_time_sec; }
    float getEnergyEma() const { return energy_ema; }

private:
    float gravity_std_;
};



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

    const float g_slow_n = g_slow.norm();
    const float g_rel_err =
        (std::isfinite(g_slow_n) && gravity_std > 1e-6f)
            ? std::fabs(g_slow_n - gravity_std) / gravity_std
            : INFINITY;

    const bool gravity_good_now =
        std::isfinite(align_sin) &&
        (align_sin <= gravity_align_max_sin) &&
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
    DetrenderConfig dcfg;
    dcfg.init_wave_freq_hz = freq_guess_hz;
    dcfg.min_wave_freq_hz  = 0.02f;
    dcfg.max_wave_freq_hz  = 1.20f;
    dcfg.baseline_cutoff_fraction = 0.25f;
    dcfg.min_baseline_cutoff_hz   = 0.003f;
    dcfg.max_baseline_cutoff_hz   = 0.25f;
    dcfg.freq_smooth_tau_s = 12.0f;
    dcfg.slope_lpf_tau_s   = 0.20f;
    dcfg.slope_rms_tau_s   = 8.0f;
    dcfg.freq_learn_axis   = 2;
    dcfg.threshold_rms_fraction  = 0.15f;
    dcfg.min_slope_threshold_abs = 0.002f;
    dcfg.max_slope_threshold_abs = 1.0e9f;
    dcfg.startup_hold_s      = 2.0f;
    dcfg.freq_timeout_cycles = 3.0f;
    dcfg.enable_wave_cleanup     = true;
    dcfg.cleanup_cutoff_fraction = 1.0f;
    dcfg.min_cleanup_cutoff_hz   = 0.003f;
    dcfg.max_cleanup_cutoff_hz   = 0.50f;
    dcfg.cleanup_stages          = 2;
    dcfg.min_dt_s = 1.0e-4f;
    dcfg.max_dt_s = 0.25f;
    dcfg.output_abs_limit = 0.0f;
    return dcfg;
}

template<typename MekfPtr, typename EnterColdFn, typename ApplyTuneFn>
inline void finalizeInitialization(MekfPtr& mekf, EnterColdFn&& enterCold, ApplyTuneFn&& applyTune) {
    enterCold();
    applyTune();
    mekf->set_exact_att_bias_Qd(true);
}

}  // namespace seastate::common
