#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Level-1 (prescribed-attitude) building blocks shared by carma_bench and
  carma3_bench: the shipping adaptation mechanics on the vertical channel, the
  operating point they produce, the shipping residual-bias prior and a plain
  float Kalman filter with Joseph scalar updates.
*/

#include <cmath>

#include <Eigen/Dense>

#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

namespace carma_bench {

namespace ou3 = ocean_imu::ou3;

// The shipping adaptation mechanics on the vertical channel alone: the same
// WavePeriodEstimator, sigma band, SeaStateAutoTuner, operating-point targets,
// EMA and r_S schedule as SeaStateFusionFilter_OU_III, minus the Mahony
// levelling (the attitude is prescribed) and the stillness detector.
class VerticalTuner : public seastate::common::SeaStateAdaptationCommon {
public:
    VerticalTuner()
        : SeaStateAdaptationCommon(EstimatorBounds{
              ou3::MAX_TUNE_FREQ_HZ, ou3::MIN_TAU_S, ou3::MAX_TAU_S, ou3::MAX_SIGMA_A,
              ou3::PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT}) {
        setOnlineTuneWarmupSec(seastate::common::defaults::STARTUP_ONLINE_TUNE_WARMUP_SEC);
    }

    void step(float dt, float a_up) {
        time_ += dt;
        startup_stage_t_ += dt;
        if (pending_) {
            committed_tau_ = tau_applied_;
            committed_sigma_ = sigma_applied_;
            committed_RS_ = RS_applied_;
            pending_ = false;
        }
        wave_period_.update(dt, a_up);
        const float wave_hz = wave_period_.getFrequencyHz();
        const float f_tune = (wave_period_.hasUsablePeriod() && std::isfinite(wave_hz) && wave_hz > 0.0f)
                             ? wave_hz : seastate::common::defaults::TUNE_FREQ_PRIOR_HZ;
        stepVarianceChannel_(dt, a_up, f_tune);
        if (!advanceStartupStage_(wave_period_.hasUsablePeriod())) return;
        const float sea_time = measureOperatingPoint_(false, 0.0f, 1.0f, 1.0f, 0.9f);
        const float rs_t = std::min(std::max(rs_target_(tau_target_, sigma_target_), ou3::MIN_R_S), ou3::MAX_R_S);
        if (!std::isfinite(tau_applied_)) {
            tau_applied_ = tau_target_; sigma_applied_ = sigma_target_; RS_applied_ = rs_t;
        }
        const float alpha = tauSigmaSmoothingAlpha_(dt, sea_time);
        const float rs_sec = seastate::common::adaptiveSmoothingHorizonSec(
            ou3::ADAPT_RS_MULT, tau_target_, rs_t, RS_applied_, ou3::ADAPT_RS_SLEW_LOG, dt);
        tau_applied_ += alpha * (tau_target_ - tau_applied_);
        sigma_applied_ += alpha * (sigma_target_ - sigma_applied_);
        RS_applied_ += (1.0f - std::exp(-dt / rs_sec)) * (rs_t - RS_applied_);
        stageOnlineTuneCommit_();
        if (online_tune_apply_pending_) { pending_ = true; online_tune_apply_pending_ = false; }
        if (startup_stage_ == StartupStage::TunerReady && std::isfinite(committed_tau_)) ready_ = true;
    }

    bool ready() const { return ready_; }
    float tau() const { return committed_tau_; }
    float sigma() const { return committed_sigma_; }
    float rs() const { return committed_RS_; }
    float ts() const { return pseudo_update_period_for_(committed_tau_); }
    float sigma_floor() const { return std::max(0.05f, band_noise_floor_sigma_()); }
    float noise_floor() const { return band_noise_floor_sigma_(); }

    static float rs_target_(float tau, float sigma) {
        // SpectralMSE law of SeaStateFusionFilter_OU_III with its constants.
        const float TS = std::min(std::max(seastate::common::defaults::PSEUDO_UPDATE_TAU_RATIO * tau,
                                           seastate::common::defaults::PSEUDO_UPDATE_PERIOD_MIN_S),
                                  ou3::PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT);
        const float sigma_aB = std::max(sigma / 0.9f, 1e-6f);
        const float qpow = std::pow(2.0f * ou3::R_S_ACCEL_NOISE_DENSITY_DEFAULT, 1.0f / 14.0f);
        const float tau2 = tau * tau;
        return ou3::R_S_MSE_COEFF_DEFAULT * qpow * std::pow(sigma_aB * tau2 * tau2, 6.0f / 7.0f) / std::sqrt(TS);
    }

private:
    WavePeriodEstimator wave_period_;
    float tau_applied_ = NAN, sigma_applied_ = NAN, RS_applied_ = NAN;
    float committed_tau_ = NAN, committed_sigma_ = NAN, committed_RS_ = NAN;
    bool pending_ = false;
    bool ready_ = false;
};

struct OperatingPoint {
    bool valid = false;
    float tau = NAN, sigma_aw = NAN, rs = NAN, ts = NAN, sigma_floor = 0.05f, noise_floor = 0.0f;
};

// Residual accelerometer-bias OU, the shipping prior.
struct BiasPrior {
    float tau = 5000.0f, q = 2.5e-7f, sigma0 = 0.004f;
};

// Plain float Kalman filter with scalar updates in Joseph form.
template<int NS>
struct LinearKF {
    using Vec = Eigen::Matrix<float,NS,1>;
    using Mat = Eigen::Matrix<float,NS,NS>;
    Vec x = Vec::Zero();
    Mat P = Mat::Zero();
    float last_nis = NAN;

    void predict(const Mat& F, const Mat& Q) {
        x = (F * x).eval();
        P = (F * P * F.transpose() + Q).eval();
        P = 0.5f * (P + P.transpose()).eval();
    }
    // z = h^T x + n with h having entries at the given indices.
    template<int M>
    void update(const int (&idx)[M], const float (&coef)[M], float z, float R) {
        Vec Ph = Vec::Zero();
        float pred = 0.0f;
        for (int k = 0; k < M; ++k) {
            Ph += P.col(idx[k]) * coef[k];
            pred += coef[k] * x(idx[k]);
        }
        float S = R;
        for (int k = 0; k < M; ++k) S += coef[k] * Ph(idx[k]);
        if (!(S > 0.0f) || !std::isfinite(S)) return;
        const Vec K = Ph / S;
        const float r = z - pred;
        last_nis = r * r / S;
        x += K * r;
        // Joseph: (I - K h^T) P (I - K h^T)^T + K R K^T
        Mat IKH = Mat::Identity();
        for (int k = 0; k < M; ++k) IKH.col(idx[k]) -= K * coef[k];
        P = (IKH * P * IKH.transpose() + K * R * K.transpose()).eval();
        P = 0.5f * (P + P.transpose()).eval();
    }
};

}  // namespace carma_bench
