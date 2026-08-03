#pragma once

/*
  WaveBandPeriodEstimator

  Estimates the zero-crossing period T_z from vertical acceleration without
  using any state from the navigation filter. Two cascaded leaky integrators
  approximate velocity and elevation over the wave band while keeping DC and
  very-low-frequency bias bounded. Long-horizon EWMA moments then give

      omega_z = sigma_v / sigma_eta,
      T_z     = 2*pi / omega_z.

  For frequencies well above the leak corner, this is the standard spectral
  moment relation omega_z = sqrt(m2/m0) for elevation.
*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

template<typename T = float>
class WaveBandPeriodEstimator {
public:
    explicit WaveBandPeriodEstimator(T leak_cutoff_hz = T(0.02),
                                     T moment_tau_sec = T(180),
                                     T warmup_sec = T(60))
    {
        setLeakCutoffHz(leak_cutoff_hz);
        setMomentTimeConstant(moment_tau_sec);
        setWarmupSec(warmup_sec);
        reset();
    }

    void reset() noexcept {
        velocity_ = T(0);
        elevation_ = T(0);
        velocity_mean_ = T(0);
        elevation_mean_ = T(0);
        velocity_second_ = T(0);
        elevation_second_ = T(0);
        elapsed_sec_ = T(0);
        accepted_samples_ = 0;
    }

    void setLeakCutoffHz(T cutoff_hz) noexcept {
        if (std::isfinite(cutoff_hz) && cutoff_hz > T(0)) {
            leak_cutoff_hz_ = cutoff_hz;
        }
    }

    void setMomentTimeConstant(T tau_sec) noexcept {
        if (std::isfinite(tau_sec) && tau_sec > T(0)) {
            moment_tau_sec_ = tau_sec;
        }
    }

    void setWarmupSec(T warmup_sec) noexcept {
        if (std::isfinite(warmup_sec) && warmup_sec >= T(0)) {
            warmup_sec_ = warmup_sec;
        }
    }

    void update(T acceleration, T dt) noexcept {
        if (!std::isfinite(acceleration) || !std::isfinite(dt) || !(dt > T(0))) {
            return;
        }

        const T lambda = T(2) * std::numbers::pi_v<T> * leak_cutoff_hz_;
        const T decay = std::exp(-lambda * dt);
        const T one_minus_decay = -std::expm1(-lambda * dt);
        const T integration_gain = one_minus_decay / lambda;
        const T old_velocity = velocity_;

        // Exact zero-order-hold update of
        //   v_dot   = a - lambda*v
        //   eta_dot = v - lambda*eta
        velocity_ = decay * old_velocity + integration_gain * acceleration;
        elevation_ = decay * (elevation_ + dt * old_velocity)
                   + acceleration * (one_minus_decay / (lambda * lambda)
                                   - dt * decay / lambda);

        const T alpha = -std::expm1(-dt / moment_tau_sec_);
        update_moment_(velocity_, alpha, velocity_mean_, velocity_second_);
        update_moment_(elevation_, alpha, elevation_mean_, elevation_second_);

        elapsed_sec_ += dt;
        ++accepted_samples_;
    }

    bool isReady() const noexcept {
        return elapsed_sec_ >= warmup_sec_
            && accepted_samples_ >= 2
            && velocityVariance() > variance_floor_
            && elevationVariance() > variance_floor_;
    }

    T getPeriodSec() const noexcept {
        const T omega = getAngularFrequencyRadPerSec();
        if (!(omega > T(0)) || !std::isfinite(omega)) {
            return std::numeric_limits<T>::quiet_NaN();
        }
        return T(2) * std::numbers::pi_v<T> / omega;
    }

    T getFrequencyHz() const noexcept {
        const T period = getPeriodSec();
        return (period > T(0) && std::isfinite(period))
            ? T(1) / period
            : std::numeric_limits<T>::quiet_NaN();
    }

    T getAngularFrequencyRadPerSec() const noexcept {
        const T var_v = velocityVariance();
        const T var_eta = elevationVariance();
        if (!(var_v > variance_floor_) || !(var_eta > variance_floor_)) {
            return std::numeric_limits<T>::quiet_NaN();
        }
        return std::sqrt(var_v / var_eta);
    }

    T velocityVariance() const noexcept {
        return std::max(T(0), velocity_second_ - velocity_mean_ * velocity_mean_);
    }

    T elevationVariance() const noexcept {
        return std::max(T(0), elevation_second_ - elevation_mean_ * elevation_mean_);
    }

    T velocityProxy() const noexcept { return velocity_; }
    T elevationProxy() const noexcept { return elevation_; }
    T elapsedSec() const noexcept { return elapsed_sec_; }
    T leakCutoffHz() const noexcept { return leak_cutoff_hz_; }
    T momentTimeConstantSec() const noexcept { return moment_tau_sec_; }

private:
    static void update_moment_(T sample, T alpha, T& mean, T& second) noexcept {
        mean += alpha * (sample - mean);
        second += alpha * (sample * sample - second);
    }

    T leak_cutoff_hz_ = T(0.02);
    T moment_tau_sec_ = T(180);
    T warmup_sec_ = T(60);
    T variance_floor_ = T(1e-12);

    T velocity_ = T(0);
    T elevation_ = T(0);
    T velocity_mean_ = T(0);
    T elevation_mean_ = T(0);
    T velocity_second_ = T(0);
    T elevation_second_ = T(0);
    T elapsed_sec_ = T(0);
    unsigned long accepted_samples_ = 0;
};
