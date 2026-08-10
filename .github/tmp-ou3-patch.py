from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text()
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{path}: expected one match, found {count}")
    p.write_text(text.replace(old, new, 1))


src = "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
test = "tests/kalman_ou_iii/tuner_schedule-test.cpp"

replace_once(
    src,
    '''    // Self-similar integral pseudo-measurement cadence T_S = c_T * tau_applied.
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
''',
    '''    // Self-similar integral pseudo-measurement cadence T_S = c_T * tau_applied.
    // Enabled by default; disabling restores the historical fixed 15 ms cadence
    // for direct old-versus-new ablation. Whenever cadence changes while Live,
    // reapply R_S so its per-update covariance stays information-rate matched.
    void setTauScaledPseudoUpdateCadence(bool flag) {
        tau_scaled_pseudo_cadence_ = flag;
        if (!mekf_) return;
        if (flag) apply_pseudo_update_cadence_();
        else mekf_->set_pseudo_update_period_s(pseudo_update_fixed_period_s_);
        if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
            apply_RS_tune_();
        }
    }
    bool tauScaledPseudoUpdateCadence() const noexcept { return tau_scaled_pseudo_cadence_; }
    void setPseudoUpdateTauRatio(float ratio) {
        if (!(std::isfinite(ratio) && ratio > 0.0f)) return;
        pseudo_update_tau_ratio_ = ratio;
        if (tau_scaled_pseudo_cadence_) {
            apply_pseudo_update_cadence_();
            if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
                apply_RS_tune_();
            }
        }
    }
    void setPseudoUpdatePeriodBounds(float min_s, float max_s) {
        if (!(std::isfinite(min_s) && std::isfinite(max_s) &&
              min_s > 0.0f && max_s >= min_s)) return;
        pseudo_update_period_min_s_ = min_s;
        pseudo_update_period_max_s_ = max_s;
        if (tau_scaled_pseudo_cadence_) {
            apply_pseudo_update_cadence_();
            if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
                apply_RS_tune_();
            }
        }
    }
''',
)

replace_once(
    src,
    '''    void apply_RS_tune_(float rs_scale = 1.0f) {
        if (!mekf_) return;
        const float s = (std::isfinite(rs_scale) && rs_scale > 0.0f)
                        ? std::min(rs_scale, 1.0f)
                        : 1.0f;
        const float RSb = std::min(std::max(tune_.RS_applied, min_R_S_), max_R_S_);
        const float rs_xy = RSb * s * R_S_xy_factor_;
        mekf_->set_RS_noise(Eigen::Vector3f(
            rs_xy,
            rs_xy,
            RSb * s
        ));
    }
''',
    '''    // set_RS_noise() accepts a standard deviation, so one S=0 update has
    // covariance r_S^2. With updates every T_S seconds, the continuous-equivalent
    // information rate is proportional to 1/(r_S^2 T_S). Preserve the historical
    // 15 ms information rate by normalizing the filter-input standard deviation:
    //     r_S,filter = r_S,base * sqrt(T_0/T_S).
    // The base tuner value remains clamped to [min_R_S_, max_R_S_]. Do not clamp
    // again after this normalization: the smallest-sea operating point is already
    // on the 0.4 base floor and must move below 0.4 when T_S > 15 ms to preserve
    // r_S^2 T_S. With unclamped T_S proportional to tau this turns the existing
    // base sigma_aw*tau^3 schedule into an effective sigma_aw*tau^(5/2) schedule
    // at the filter input.
    float pseudo_update_information_rate_scale_() const noexcept {
        if (!tau_scaled_pseudo_cadence_ || !mekf_) return 1.0f;
        const float period = mekf_->get_pseudo_update_period_s();
        if (!(std::isfinite(period) && period > 0.0f)) return 1.0f;
        return std::sqrt(pseudo_update_fixed_period_s_ / period);
    }

    void apply_RS_tune_(float rs_scale = 1.0f) {
        if (!mekf_) return;
        const float s = (std::isfinite(rs_scale) && rs_scale > 0.0f)
                        ? std::min(rs_scale, 1.0f)
                        : 1.0f;
        const float RSbase = std::min(std::max(tune_.RS_applied, min_R_S_), max_R_S_);
        const float RSb = RSbase * pseudo_update_information_rate_scale_();
        const float rs_xy = RSb * s * R_S_xy_factor_;
        mekf_->set_RS_noise(Eigen::Vector3f(
            rs_xy,
            rs_xy,
            RSb * s
        ));
    }
''',
)

replace_once(
    test,
    '''    const float staged_tau = f.tune_.tau_applied;
    const float staged_rs = std::min(std::max(f.tune_.RS_applied, f.min_R_S_), f.max_R_S_);
    const float staged_rs_var = staged_rs * staged_rs;
    const float staged_pseudo = std::min(
        std::max(f.pseudo_update_tau_ratio_ * staged_tau,
                 f.pseudo_update_period_min_s_),
        f.pseudo_update_period_max_s_);
''',
    '''    const float staged_tau = f.tune_.tau_applied;
    const float staged_rs_base = std::min(std::max(f.tune_.RS_applied, f.min_R_S_), f.max_R_S_);
    const float staged_pseudo = std::min(
        std::max(f.pseudo_update_tau_ratio_ * staged_tau,
                 f.pseudo_update_period_min_s_),
        f.pseudo_update_period_max_s_);
    const float staged_cadence_scale = std::sqrt(
        PSEUDO_UPDATE_PERIOD_NOMINAL_S / staged_pseudo);
    const float staged_rs = staged_rs_base * staged_cadence_scale;
    const float staged_rs_var = staged_rs * staged_rs;
''',
)

replace_once(
    test,
    '''    if (!near(f.getPseudoUpdatePeriodSec() / f.getTauApplied(),
              f.getPseudoUpdateTauRatio(), 1e-5f)) {
        std::cerr << "FAIL: T_S/tau is not invariant after the staged update\\n";
        return 1;
    }
''',
    '''    if (!near(f.getPseudoUpdatePeriodSec() / f.getTauApplied(),
              f.getPseudoUpdateTauRatio(), 1e-5f)) {
        std::cerr << "FAIL: T_S/tau is not invariant after the staged update\\n";
        return 1;
    }
    const float expected_info_product =
        staged_rs_base * staged_rs_base * PSEUDO_UPDATE_PERIOD_NOMINAL_S;
    if (!near(f.mekf_->R_S(2, 2) * f.getPseudoUpdatePeriodSec(),
              expected_info_product, 1e-5f)) {
        std::cerr << "FAIL: cadence compensation did not preserve R_S*T_S information rate\\n";
        return 1;
    }
''',
)

replace_once(
    test,
    '''    f.setTauScaledPseudoUpdateCadence(false);
    if (!near(f.getPseudoUpdatePeriodSec(), 0.015f)) {
        std::cerr << "FAIL: fixed-cadence ablation did not restore 15 ms\\n";
        return 1;
    }
    f.setTauScaledPseudoUpdateCadence(true);
    if (!near(f.getPseudoUpdatePeriodSec(),
              std::min(std::max(f.getPseudoUpdateTauRatio() * f.getTauApplied(),
                                f.pseudo_update_period_min_s_),
                       f.pseudo_update_period_max_s_))) {
        std::cerr << "FAIL: re-enabling tau-scaled cadence did not reapply T_S=c_T*tau\\n";
        return 1;
    }

    std::cout << "OU-III predictable tuner scheduling and tau-scaled pseudo cadence passed\\n";
''',
    '''    f.setTauScaledPseudoUpdateCadence(false);
    const float fixed_rs_var = staged_rs_base * staged_rs_base;
    if (!near(f.getPseudoUpdatePeriodSec(), PSEUDO_UPDATE_PERIOD_NOMINAL_S) ||
        !near(f.mekf_->R_S(2, 2), fixed_rs_var)) {
        std::cerr << "FAIL: fixed-cadence ablation did not restore 15 ms/base R_S\\n";
        return 1;
    }
    f.setTauScaledPseudoUpdateCadence(true);
    if (!near(f.getPseudoUpdatePeriodSec(), staged_pseudo) ||
        !near(f.mekf_->R_S(2, 2), staged_rs_var)) {
        std::cerr << "FAIL: re-enabling tau-scaled cadence did not restore compensated schedule\\n";
        return 1;
    }

    std::cout << "OU-III predictable tuner scheduling, tau-scaled cadence, and information-rate compensation passed\\n";
''',
)
