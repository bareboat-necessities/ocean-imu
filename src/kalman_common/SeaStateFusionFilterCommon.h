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
