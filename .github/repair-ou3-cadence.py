from pathlib import Path

p = Path('src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h')
s = p.read_text()

def once(old, new):
    global s
    n = s.count(old)
    if n != 1:
        raise SystemExit(f'expected one match, found {n}: {old[:80]!r}')
    s = s.replace(old, new, 1)

once(
    'constexpr float FREQ_SMOOTHER_DT = 1.0f / 200.0f;\n\nstruct TuneState {',
    '''constexpr float FREQ_SMOOTHER_DT = 1.0f / 200.0f;

// Self-similar S=0 pseudo-measurement cadence.  The deployed wrapper historically
// used T_S=15 ms while its initial applied OU time constant is tau=1.1 s, so
// this ratio preserves that exact operating point and scales cadence thereafter.
constexpr float PSEUDO_UPDATE_PERIOD_NOMINAL_S = 0.015f;
constexpr float PSEUDO_UPDATE_TAU_NOMINAL_S = 1.1f;
constexpr float PSEUDO_UPDATE_TAU_RATIO_DEFAULT =
    PSEUDO_UPDATE_PERIOD_NOMINAL_S / PSEUDO_UPDATE_TAU_NOMINAL_S;
// A pseudo update cannot occur more often than the nominal 200 Hz IMU schedule.
// The upper guard is inactive over the current tau <= 12 s operating envelope.
constexpr float PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT = FREQ_SMOOTHER_DT;
constexpr float PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT = 0.25f;

struct TuneState {''')

once(
    '    void setAwCovarianceSyncCongruent(bool flag) { congruent_aw_cov_sync_ = flag; }\n    bool awCovarianceSyncCongruent() const noexcept { return congruent_aw_cov_sync_; }\n\n    // Freeze the online tuner',
    '''    void setAwCovarianceSyncCongruent(bool flag) { congruent_aw_cov_sync_ = flag; }
    bool awCovarianceSyncCongruent() const noexcept { return congruent_aw_cov_sync_; }

    // Self-similar integral pseudo-measurement cadence T_S = c_T * tau_applied.
    // Enabled by default; disabling restores the historical fixed 15 ms cadence
    // for direct old-versus-new ablation.
    void setTauScaledPseudoUpdateCadence(bool flag) {
        tau_scaled_pseudo_cadence_ = flag;
        if (!mekf_) return;
        if (flag) apply_pseudo_update_cadence_();
        else mekf_->set_pseudo_update_period_s(pseudo_update_fixed_period_s_);
    }
    bool tauScaledPseudoUpdateCadence() const noexcept { return tau_scaled_pseudo_cadence_; }
    void setPseudoUpdateTauRatio(float ratio) {
        if (!(std::isfinite(ratio) && ratio > 0.0f)) return;
        pseudo_update_tau_ratio_ = ratio;
        if (tau_scaled_pseudo_cadence_) apply_pseudo_update_cadence_();
    }
    void setPseudoUpdatePeriodBounds(float min_s, float max_s) {
        if (!(std::isfinite(min_s) && std::isfinite(max_s) &&
              min_s > 0.0f && max_s >= min_s)) return;
        pseudo_update_period_min_s_ = min_s;
        pseudo_update_period_max_s_ = max_s;
        if (tau_scaled_pseudo_cadence_) apply_pseudo_update_cadence_();
    }
    float getPseudoUpdateTauRatio() const noexcept { return pseudo_update_tau_ratio_; }
    float getPseudoUpdatePeriodSec() const noexcept {
        return mekf_ ? mekf_->get_pseudo_update_period_s() : NAN;
    }

    // Freeze the online tuner''')

once(
    '    // sync_covariance is set only by discrete reconfiguration events. The\n',
    '''    void apply_pseudo_update_cadence_() {
        if (!mekf_ || !tau_scaled_pseudo_cadence_) return;
        const float tau = tune_.tau_applied;
        if (!(std::isfinite(tau) && tau > 0.0f)) return;
        const float requested = pseudo_update_tau_ratio_ * tau;
        const float period = std::min(
            std::max(requested, pseudo_update_period_min_s_),
            pseudo_update_period_max_s_);
        mekf_->set_pseudo_update_period_s(period);
    }

    // sync_covariance is set only by discrete reconfiguration events. The
''')

once(
    '        mekf_->set_aw_time_constant(tune_.tau_applied);\n\n        const float sigma_floor',
    '''        mekf_->set_aw_time_constant(tune_.tau_applied);
        // Commit the S=0 cadence with the same applied tau so T_S/tau remains
        // constant apart from explicit safety clamps.
        apply_pseudo_update_cadence_();

        const float sigma_floor''')

once(
    '    double last_aw_cov_sync_sec_ = 0.0;\n\n    bool  wave_band_tuning_',
    '''    double last_aw_cov_sync_sec_ = 0.0;

    bool  tau_scaled_pseudo_cadence_ = true;
    float pseudo_update_tau_ratio_ = PSEUDO_UPDATE_TAU_RATIO_DEFAULT;
    float pseudo_update_period_min_s_ = PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT;
    float pseudo_update_period_max_s_ = PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT;
    float pseudo_update_fixed_period_s_ = PSEUDO_UPDATE_PERIOD_NOMINAL_S;

    bool  wave_band_tuning_''')

p.write_text(s)
