#pragma once
#include <algorithm>
#include <cmath>
#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

namespace ocean_imu::ou3 {

// Horizontal spectral moments in the world frame. For a single encounter
// frequency omega, the displacement primitive S and velocity satisfy
// S + v/omega² = 0. An irregular sea has a finite residual variance.
class HorizontalWaveBand {
public:
    void reset() { *this = {}; }

    void update(float dt, float period, const Eigen::Vector2f& integral,
                const Eigen::Vector2f& position, const Eigen::Vector2f& velocity) {
        if (!(dt > 0.0f && period > 0.0f) || !std::isfinite(dt) || !std::isfinite(period) ||
            !integral.allFinite() || !position.allFinite() || !velocity.allFinite()) {
            reset();
            return;
        }
        if (cycles_ == 0.0f) {
            mean_s_ = integral;
            mean_p_ = position;
            mean_v_ = velocity;
        }
        const float w = -std::expm1(-dt / (8.0f * period));
        const Eigen::Vector2f ds = integral - mean_s_;
        const Eigen::Vector2f dp = position - mean_p_;
        const Eigen::Vector2f dv = velocity - mean_v_;
        mean_s_ += w * ds;
        mean_p_ += w * dp;
        mean_v_ += w * dv;
        // Half the horizontal trace gives a per-axis variance and is
        // invariant to the choice of horizontal north.
        const auto moment = [w](float old, float product) {
            return (1.0f - w) * (old + (0.5f * w) * product);
        };
        var_s_ = moment(var_s_, ds.squaredNorm());
        var_p_ = moment(var_p_, dp.squaredNorm());
        var_v_ = moment(var_v_, dv.squaredNorm());
        cov_sv_ = moment(cov_sv_, ds.dot(dv));
        cycles_ = std::min(4.0f, cycles_ + dt / period);
        if (!(var_v_ > 1e-12f)) {
            beta_ = 0.0f;
            return;
        }
        constexpr float two_pi = 6.283185307179586f;
        const float inverse_omega = period / two_pi;
        const float nominal = inverse_omega * inverse_omega;
        beta_ = std::clamp(var_p_ / var_v_, 0.25f * nominal, 4.0f * nominal);
        residual_variance_ = std::max(0.0f, var_s_ + beta_ * (beta_ * var_v_ + 2.0f * cov_sv_));
    }

    float beta() const { return beta_; }
    float residualVariance() const { return residual_variance_; }
    float weight() const { return std::clamp(cycles_ - 3.0f, 0.0f, 1.0f); }

private:
    Eigen::Vector2f mean_s_ = Eigen::Vector2f::Zero();
    Eigen::Vector2f mean_p_ = Eigen::Vector2f::Zero();
    Eigen::Vector2f mean_v_ = Eigen::Vector2f::Zero();
    float var_s_ = 0.0f, var_p_ = 0.0f, var_v_ = 0.0f, cov_sv_ = 0.0f;
    float cycles_ = 0.0f, beta_ = 0.0f, residual_variance_ = 0.0f;
};

} // namespace ocean_imu::ou3
