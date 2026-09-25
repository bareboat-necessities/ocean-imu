#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

// Adaptation mechanics shared by the OU-II and OU-III orchestrators.
//
// What is common is the *measurement* of the OU operating point and the
// *schedule* on which it reaches the estimator:
//
//   * the Cold -> TunerWarm -> TunerReady -> Live stage machine;
//   * the period-scaled sigma band ahead of SeaStateAutoTuner, and the
//     bench-noise floor referred through that band's own coefficients;
//   * the operating-point targets tau = c_tau T_z/2 and
//     sigma_aw = c_sigma sigma_a,B, with their safety clamps;
//   * the self-similar tau/sigma_aw EMA and the one-sample staging of the
//     smoothed candidate (exogeneity as a timing property, below);
//   * the tau-scaled pseudo-measurement cadence T_S = c_T tau;
//   * the periodic a_w covariance-maintenance clock.
//
// What is NOT here is anything a regularization law decides: which pseudo
// measurements exist, how their targets follow from (tau, sigma_aw, T_S), how
// they are smoothed, and how the OU prior and the pseudo covariances are
// written into the MEKF.  Those are the estimator-specific layers of
// SeaStateFusionFilter_OU_II (dual p/v PhysicalMSE) and
// SeaStateFusionFilter_OU_III (integral SpectralMSE), and they stay there.
//
// EXOGENEITY IS A TIMING PROPERTY, NOT JUST A SIGNAL CHOICE.  Feeding the
// tuner from the complementary observer keeps its *inputs* independent of the
// filter.  That is necessary and not sufficient: the schedule smoothed during
// step k is committed at the top of the next update(), before y_{k+1} reaches
// the MEKF, so the active schedule at step k+1 is measurable with respect to
// data through k.  stageOnlineTuneCommit_() sets the flag; each wrapper's
// apply_pending_online_tune_() consumes it.

#include <algorithm>
#include <cmath>

#include "kalman_common/SeaStateFusionDefaults.h"
#include "tuner/AdaptiveWaveBandPass.h"
#include "tuner/SeaStateAdaptationLimits.h"
#include "tuner/SeaStateAutoTuner.h"

namespace seastate::common {

class SeaStateAdaptationCommon {
public:
    enum class StartupStage {
        Cold,        // just booted or just had a big tilt reset
        TunerWarm,   // front end running, tuner collecting stats
        TunerReady,  // tuner trusted, MEKF still held by an external bootstrap
        Live         // MEKF owns the attitude; full adaptation & extras allowed
    };

    // Estimator-specific values of the shared mechanics.  These differ on
    // purpose between the families; each wrapper documents its own.
    struct EstimatorBounds {
        float max_tune_freq_hz;
        float min_tau_s;
        float max_tau_s;
        float max_sigma_a;
        float pseudo_update_period_max_s;
    };

    explicit SeaStateAdaptationCommon(const EstimatorBounds& b)
        : max_tune_freq_hz_(b.max_tune_freq_hz),
          min_tau_s_(b.min_tau_s),
          max_tau_s_(b.max_tau_s),
          max_sigma_a_(b.max_sigma_a),
          pseudo_update_period_max_s_(b.pseudo_update_period_max_s) {}

    StartupStage getStartupStage() const noexcept { return startup_stage_; }
    bool isAdaptiveLive() const noexcept { return startup_stage_ == StartupStage::Live; }

    // The operating point is trustworthy.  It is reached while the MEKF is
    // still held, and it is one of the conditions the caller waits on before
    // handing the attitude over.
    bool isTunerReady() const noexcept {
        return startup_stage_ == StartupStage::TunerReady ||
               startup_stage_ == StartupStage::Live;
    }

    void setAccNoiseFloorSigma(float s) {
        if (std::isfinite(s) && s > 0.0f) acc_noise_floor_sigma_ = s;
    }
    float getAccNoiseFloorSigma() const noexcept { return acc_noise_floor_sigma_; }

    // Dimensionless sigma-band shape.  Defaults correspond to the theorem's
    // measurement-only soft band [0.5,4] in units of f_tune.  Changing these
    // ratios preserves the similarity structure as long as they remain fixed
    // across sea states; absolute limits below are safety clamps only.
    void setSigmaWaveBandRatios(float low_ratio, float high_ratio) {
        sigma_wave_band_.setRatios(low_ratio, high_ratio);
    }
    void setSigmaWaveBandLimitsHz(float min_hz, float max_hz) {
        sigma_wave_band_.setLimitsHz(min_hz, max_hz);
    }
    float getSigmaWaveBandLowHz() const noexcept { return sigma_wave_band_.lowHz(); }
    float getSigmaWaveBandHighHz() const noexcept { return sigma_wave_band_.highHz(); }
    float getSigmaWaveBandLowRatio() const noexcept { return sigma_wave_band_.lowRatio(); }
    float getSigmaWaveBandHighRatio() const noexcept { return sigma_wave_band_.highRatio(); }
    float getSigmaBandNoiseStd() const noexcept { return band_noise_floor_sigma_(); }

    void enableClamp(bool flag = true) { enable_clamp_ = flag; }
    void enableTuner(bool flag = true) {
        enable_tuner_ = flag;
        if (!flag) online_tune_apply_pending_ = false;
    }
    bool tunerEnabled() const noexcept { return enable_tuner_; }

    void setTauBounds(float min_tau_s, float max_tau_s) {
        if (!std::isfinite(min_tau_s) || !std::isfinite(max_tau_s)) return;
        if (min_tau_s <= 0.0f || max_tau_s <= min_tau_s) return;
        min_tau_s_ = min_tau_s;
        max_tau_s_ = max_tau_s;
    }

    void setMaxSigmaA(float max_sigma_a) {
        if (!std::isfinite(max_sigma_a) || max_sigma_a <= 0.0f) return;
        max_sigma_a_ = max_sigma_a;
    }

    // Bounds on the wave-band tuning frequency.  Distinct from setFreqBounds(),
    // which bounds the acceleration-band tracker that drives wave direction.
    void setTuneFreqBounds(float min_hz, float max_hz) {
        if (!(std::isfinite(min_hz) && std::isfinite(max_hz))) return;
        if (!(min_hz > 0.0f && max_hz > min_hz)) return;
        min_tune_freq_hz_ = min_hz;
        max_tune_freq_hz_ = max_hz;
    }

    // Fixed-seconds compatibility mode.  Explicitly selecting it disables the
    // period-scaled common tau/sigma EMA.
    void setAdaptationTimeConstants(float tau_sec) {
        if (std::isfinite(tau_sec) && tau_sec > 0.0f) {
            adapt_tau_sec_ = tau_sec;
            adapt_tau_sea_periods_ = 0.0f;
        }
    }

    // Common tau/sigma EMA horizon as a fraction/multiple of measured
    // T_sea = T_z/2.  This changes only averaging memory; the operating-point
    // coefficients, the sigma variance K-period horizon, and the estimator's
    // own regularizer horizons are untouched.
    void setAdaptationSeaPeriods(float periods) {
        if (std::isfinite(periods) && periods > 0.0f) {
            adapt_tau_sea_periods_ = periods;
        }
    }
    float getAdaptationSeaPeriods() const noexcept { return adapt_tau_sea_periods_; }

    // Compatibility no-ops: frequency smoothing belongs to
    // WavePeriodEstimator's canonical log-period state, not this tuner.
    void setTunerFreqSmoothingSeaPeriods(float /*periods*/) {}
    void setTunerFreqSmoothingTimeConstant(float /*tau_sec*/) {}
    float getTunerFreqSmoothingSeaPeriods() const noexcept {
        return tuner_.getFrequencySmoothingSeaPeriods();
    }

    void setAdaptationUpdatePeriod(float every_sec) {
        if (std::isfinite(every_sec) && every_sec > 0.0f) adapt_every_secs_ = every_sec;
    }

    void setOnlineTuneWarmupSec(float warmup_sec) {
        if (std::isfinite(warmup_sec) && warmup_sec >= 0.0f) online_tune_warmup_sec_ = warmup_sec;
    }

    // sigma_a averaging horizon, in periods of the tuning frequency, and its
    // absolute clamps in seconds.  The requested horizon is K/f_tune.
    void setSigmaVarianceKPeriods(float k) { tuner_.setKPeriods(k); }
    float getSigmaVarianceKPeriods() const noexcept { return tuner_.getKPeriods(); }
    void setSigmaVarianceHorizonBounds(float min_s, float max_s) {
        tuner_.setVarianceHorizonBounds(min_s, max_s);
    }
    // Horizon currently in force [s] for the first- and second-moment EWMAs.
    float getSigmaVarianceHorizonSec() const noexcept {
        return tuner_.getVarianceHorizonSec();
    }

    // Variance measured after the period-scaled sigma band.
    inline float getAccelVariance() const noexcept { return tuner_.getAccelVariance(); }

    inline float getTauTarget()   const noexcept { return tau_target_; }
    inline float getSigmaTarget() const noexcept { return sigma_target_; }

    // Policy for the latent-acceleration marginal P_{a_w a_w}.
    //
    // Default (true): once per adaptation period the marginal is re-aligned
    // with the stationary OU covariance, keeping the cross-covariances the
    // filter has learned.  This is a deliberate bounded covariance inflation.
    // It stops the a_w marginal from settling far below the level the process
    // model considers stationary, which keeps the accelerometer gain
    // responsive when the sea state changes.  It is not free -- it discards
    // posterior information at the adaptation cadence -- so the alternative
    // is available and measured rather than assumed.
    //
    // With false, the marginal is aligned only at discrete reconfiguration
    // events (construction, the transition to Live, setFixedTuning()) and a
    // changed stationary scale reaches the filter solely through the discrete
    // OU process covariance.
    //
    // The policy is applied independently of the tuner so that fixed-tuning
    // modes run it too.  Otherwise an adaptive-versus-fixed comparison would
    // confound whether the parameters adapt with whether part of the
    // covariance is periodically re-aligned.
    void setPeriodicAwCovarianceSync(bool flag) {
        periodic_aw_cov_sync_ = flag;
        last_aw_cov_sync_sec_ = time_;
    }
    bool periodicAwCovarianceSync() const noexcept { return periodic_aw_cov_sync_; }

    bool tauScaledPseudoUpdateCadence() const noexcept { return tau_scaled_pseudo_cadence_; }
    float getPseudoUpdateTauRatio() const noexcept { return pseudo_update_tau_ratio_; }

protected:
    // Bench white-noise sigma referred to the current adaptive band.  The band
    // is time varying, so the gain has to come from the filter's own state
    // rather than from a closed-form transfer function.
    float band_noise_floor_sigma_() const noexcept {
        if (!sigma_wave_band_.isReady()) {
            return acc_noise_floor_sigma_;
        }
        const float gain = sigma_wave_band_.whiteNoiseVarianceGain();
        if (!(std::isfinite(gain) && gain >= 0.0f)) return acc_noise_floor_sigma_;
        return acc_noise_floor_sigma_ * std::sqrt(gain);
    }

    // Variance channel: the period-scaled band, then the tuner.
    //
    // The previous smoothed tuner frequency sets the current sample's band
    // corners.  This keeps the band motion smooth and one-sample predictable
    // while remaining measurement-only.  Until the frequency EMA is ready,
    // fall back to the current external frequency estimate.
    void stepVarianceChannel_(float dt, float a_vertical_measurement,
                              float freq_hz_for_tuner) {
        float f_for_sigma_band = tuner_.isFreqReady()
            ? tuner_.getFrequencyHz()
            : freq_hz_for_tuner;
        const float f_tune_floor = min_tune_freq_hz_;
        const float f_tune_ceil = max_tune_freq_hz_;
        if (!std::isfinite(f_for_sigma_band) || f_for_sigma_band < f_tune_floor) {
            f_for_sigma_band = f_tune_floor;
        }
        f_for_sigma_band = std::min(f_for_sigma_band, f_tune_ceil);

        const float a_for_variance =
            sigma_wave_band_.step(a_vertical_measurement, dt, f_for_sigma_band);

        tuner_.update(dt, a_for_variance, freq_hz_for_tuner);
    }

    // Startup stage logic.  Returns false while Cold: the tuner accumulates
    // but no target is formed until the warmup has elapsed.
    bool advanceStartupStage_(bool wave_period_usable) {
        switch (startup_stage_) {
            case StartupStage::Cold:
                if (startup_stage_t_ >= online_tune_warmup_sec_) {
                    startup_stage_   = StartupStage::TunerWarm;
                    startup_stage_t_ = 0.0f;
                }
                return false;

            case StartupStage::TunerWarm:
                if (!tuner_.isFreqReady()) return false;
                if (tuner_.isReady() && wave_period_usable) {
                    // The operating point and measured wave period are
                    // trusted, but the attitude is not this filter's to
                    // decide.  Park here and let the bootstrap call goLive()
                    // once it has tilt and north.
                    startup_stage_   = StartupStage::TunerReady;
                    startup_stage_t_ = 0.0f;
                }
                break;

            case StartupStage::TunerReady:
            case StartupStage::Live:
                break;
        }
        return true;
    }

    // Operating-point targets tau_target_ = tau_coeff T_z/2 and
    // sigma_target_ = sigma_coeff sigma_a,B from the tuner.
    //
    // The two coefficients and still_decay_sec -- the time constant with
    // which the wave variance is attenuated while the stillness detector
    // reports still water -- are estimator-specific fits and are passed in by
    // each wrapper.  Returns the measured sea time T_sea = T_z/2 the tau/sigma
    // EMA is scaled by.
    float measureOperatingPoint_(bool is_still, float still_time_sec,
                                 float still_decay_sec,
                                 float tau_coeff, float sigma_coeff) {
        // The tuning frequency is a wave-band quantity and is bounded by the
        // wave band, not by the tracker's bounds: a developed sea has
        // T_z = 8.6 s, i.e. 0.12 Hz, well under the 0.2 Hz the tracker is
        // bounded to, and setFreqBounds() moves the demodulator carrier rather
        // than the OU operating point.
        const float f_tune_floor = min_tune_freq_hz_;
        const float f_tune_ceil = max_tune_freq_hz_;
        float f_tune = tuner_.getFrequencyHz();
        if (!std::isfinite(f_tune) || f_tune < f_tune_floor) f_tune = f_tune_floor;
        if (f_tune > f_tune_ceil) f_tune = f_tune_ceil;

        const float band_noise_sigma = band_noise_floor_sigma_();
        const float var_noise = band_noise_sigma * band_noise_sigma;
        float var_total = var_noise;
        if (tuner_.isVarReady()) {
            var_total = std::max(0.0f, tuner_.getAccelVariance());
        }
        float var_wave = var_total - var_noise;
        if (var_wave < 0.0f) var_wave = 0.0f;

        if (is_still) {
            const float still_t = std::max(0.0f, still_time_sec);
            float atten = std::exp(-still_t / still_decay_sec);
            atten = std::min(std::max(atten, 0.0f), 1.0f);
            var_wave *= atten;
        }

        var_wave = std::max(var_wave, 1e-6f);
        float sigma_wave = std::sqrt(var_wave);
        float tau_raw = tau_coeff * 0.5f / f_tune;

        if (enable_clamp_) {
            tau_target_   = std::min(std::max(tau_raw, min_tau_s_), max_tau_s_);
            sigma_target_ = std::min(sigma_wave * sigma_coeff, max_sigma_a_);
        } else {
            tau_target_   = tau_raw;
            sigma_target_ = sigma_wave;
        }

        if (!tuner_.isVarReady()) {
            sigma_target_ = std::max(sigma_target_, std::max(0.05f, band_noise_sigma));
        }

        return 0.5f / f_tune;  // T_sea = T_z/2
    }

    // EMA factor of the common tau/sigma_aw channel: 0.40 T_sea by default, a
    // fixed number of seconds in the fixed-horizon ablation.
    float tauSigmaSmoothingAlpha_(float dt, float sea_time_sec) const {
        float adapt_sec = adapt_tau_sec_;
        if (adapt_tau_sea_periods_ > 0.0f &&
            std::isfinite(sea_time_sec) && sea_time_sec > 0.0f) {
            const float safe_sea_time =
                seastate::tuner::limits::clampDynamicEmaTimeScaleSec(sea_time_sec);
            adapt_sec = seastate::tuner::limits::clampDynamicEmaHorizonSec(
                adapt_tau_sea_periods_ * safe_sea_time, dt);
        }
        return 1.0f - std::exp(-dt / adapt_sec);
    }

    // Every valid physical sample contributes to every EMA.  The activation
    // cadence is deliberately kept separate from that sampling: it throttles
    // parameter commits, not physical measurements or EMA input.  y_k may
    // change the smoothed candidate, but the MEKF keeps the schedule that was
    // active before y_k arrived; the candidate is committed at the beginning
    // of the next sample.
    void stageOnlineTuneCommit_() {
        if (time_ - last_adapt_time_sec_ > adapt_every_secs_) {
            online_tune_apply_pending_ = true;
            last_adapt_time_sec_ = time_;
        }
    }

    // Pseudo-update period the cadence scheduler selects for a given tau:
    // T_S = c_T tau, clamped.  Used both to commit the cadence and by laws
    // that must refer their targets to the same operating point, including
    // the safety clamps.
    float pseudo_update_period_for_(float tau) const noexcept {
        if (!tau_scaled_pseudo_cadence_) return pseudo_update_fixed_period_s_;
        if (!(std::isfinite(tau) && tau > 0.0f)) return pseudo_update_fixed_period_s_;
        return std::min(std::max(pseudo_update_tau_ratio_ * tau,
                                 pseudo_update_period_min_s_),
                        pseudo_update_period_max_s_);
    }

    // Commit the pseudo-update cadence for the applied tau, so T_S/tau stays
    // constant apart from explicit safety clamps.
    template <typename Mekf>
    void apply_pseudo_update_cadence_(Mekf* mekf, float tau_applied) {
        if (!mekf || !tau_scaled_pseudo_cadence_) return;
        const float tau = tau_applied;
        if (!(std::isfinite(tau) && tau > 0.0f)) return;
        const float requested = pseudo_update_tau_ratio_ * tau;
        const float period = std::min(
            std::max(requested, pseudo_update_period_min_s_),
            pseudo_update_period_max_s_);
        mekf->set_pseudo_update_period_s(period);
    }

    // One pseudo update has covariance r^2; with updates every T_S seconds the
    // continuous-equivalent information rate is proportional to 1/(r^2 T_S).
    // A law whose base target does not contain the realized cadence preserves
    // the nominal 15 ms information rate by renormalizing the filter-input
    // standard deviation, r_filter = r_base sqrt(T_0/T_S).  Not re-clamped
    // afterwards: the smallest-sea operating point may already sit on a base
    // floor and must be allowed below it when T_S > T_0.
    float cadenceRenormalization_(float realized_period_s) const noexcept {
        if (!tau_scaled_pseudo_cadence_) return 1.0f;
        const float period = realized_period_s;
        if (!(std::isfinite(period) && period > 0.0f)) return 1.0f;
        return std::sqrt(pseudo_update_fixed_period_s_ / period);
    }

    // Covariance-maintenance clock; see setPeriodicAwCovarianceSync().
    bool periodicAwCovSyncDue_() const {
        if (!periodic_aw_cov_sync_) return false;
        if (startup_stage_ != StartupStage::Live) return false;
        return !(time_ - last_aw_cov_sync_sec_ <= adapt_every_secs_);
    }

    void resetAdaptation_() {
        sigma_wave_band_.reset();
        tuner_.reset();
        last_adapt_time_sec_ = time_;
        last_aw_cov_sync_sec_ = time_;
        online_tune_apply_pending_ = false;
    }

    StartupStage startup_stage_   = StartupStage::Cold;
    float        startup_stage_t_ = 0.0f;

    double time_                = 0.0;
    double last_adapt_time_sec_ = 0.0;

    bool enable_clamp_ = true;
    bool enable_tuner_ = true;
    bool online_tune_apply_pending_ = false;

    // Covariance-inflation policy; see setPeriodicAwCovarianceSync.
    bool   periodic_aw_cov_sync_ = true;
    double last_aw_cov_sync_sec_ = 0.0;

    bool  tau_scaled_pseudo_cadence_    = true;
    float pseudo_update_tau_ratio_      = defaults::PSEUDO_UPDATE_TAU_RATIO;
    float pseudo_update_period_min_s_   = defaults::PSEUDO_UPDATE_PERIOD_MIN_S;
    float pseudo_update_fixed_period_s_ = defaults::PSEUDO_UPDATE_PERIOD_NOMINAL_S;

    float min_tune_freq_hz_ = defaults::MIN_TUNE_FREQ_HZ;
    float max_tune_freq_hz_;
    float min_tau_s_;
    float max_tau_s_;
    float max_sigma_a_;
    float pseudo_update_period_max_s_;

    float adapt_tau_sec_          = defaults::ADAPT_TAU_SEC;
    float adapt_tau_sea_periods_  = defaults::ADAPT_TAU_SEA_PERIODS;
    float adapt_every_secs_       = defaults::ADAPT_EVERY_SECS;
    float online_tune_warmup_sec_ = 5.0f;

    SeaStateAutoTuner    tuner_;
    AdaptiveWaveBandPass sigma_wave_band_{
        defaults::SIGMA_BAND_LOW_RATIO,
        defaults::SIGMA_BAND_HIGH_RATIO,
        defaults::SIGMA_BAND_MIN_HZ,
        defaults::SIGMA_BAND_MAX_HZ};

    float tau_target_   = NAN;
    float sigma_target_ = NAN;

    float acc_noise_floor_sigma_ = defaults::ACC_NOISE_FLOOR_SIGMA;
};

}  // namespace seastate::common
