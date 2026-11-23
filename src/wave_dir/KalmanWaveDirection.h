#ifndef KALMAN_WAVE_DIRECTION_H
#define KALMAN_WAVE_DIRECTION_H
#include "util/HotLoop.h"

/*
  Copyright 2025-2026, Mikhail Grushinskiy

  Phase-invariant estimator for the horizontal wave-propagation axis from IMU
  horizontal acceleration. The result is axial (modulo 180 degrees);
  propagation sense is resolved separately by WaveDirectionDetector.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

#ifdef KALMAN_WAVE_DIRECTION_TEST
#include <fstream>
#include <iostream>
#include <random>
#endif

class EIGEN_ALIGN_MAX KalmanWaveDirection {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    KalmanWaveDirection(float initialOmega, float directionReportAlpha = 0.0025f)
        : omega(initialOmega),
          direction_report_tau_sec(alphaToTau(directionReportAlpha, REFERENCE_DT_SEC))
    {
        reset();
    }

    void reset() {
        coefficient_estimate.setZero();
        p_cc = p_ss = 1.0f;
        p_cs = 0.0f;
        phase = 0.0f;
        invalidateCurrentAxis();

        lastStableAmplitude = 0.0f;
        lastStableConfidence = 0.0f;
        lastStableLinearity = 0.0f;
        lastStableVariance = 2.0f;
        lastStableDir = Eigen::Vector2f(1.0f, 0.0f);

        have_report_axis = false;
        report_cos_2theta = 1.0f;
        report_sin_2theta = 0.0f;
        raw_double_axis = Eigen::Vector2f::UnitX();
    }

    OCEAN_IMU_SAMPLE_IRAM void update(float ax, float ay, float currentOmega, float deltaT) {
        if (!(deltaT > 0.0f) || !std::isfinite(deltaT) ||
            !std::isfinite(ax) || !std::isfinite(ay)) {
            invalidateCurrentAxis();
            return;
        }

        if (std::isfinite(currentOmega) && currentOmega > 0.0f) {
            const float scale = std::max(std::fabs(omega), 1e-6f);
            if (!(omega > 0.0f) ||
                std::fabs(currentOmega - omega) > 0.01f * scale) {
                omega = currentOmega;
            }
        }
        if (!(omega > 0.0f) || !std::isfinite(omega)) {
            invalidateCurrentAxis();
            return;
        }

        updatePhase(deltaT);

        const float c = std::cos(phase);
        const float s = std::sin(phase);

        // z = C cos(phi) + S sin(phi), where C and S are two-dimensional
        // horizontal coefficient vectors. Estimating both quadratures removes
        // a carrier-phase singularity that arises with a cosine-only model.
        // The x and y measurements have independent, identical isotropic
        // noise. Starting from identity, their two quadrature covariance
        // blocks stay identical and every cross-axis entry stays zero.
        const float q = (deltaT / REFERENCE_DT_SEC) * process_noise;
        const float cc = (p_cc + q) + 1e-10f;
        const float ss = (p_ss + q) + 1e-10f;
        const float ch = cc * c + p_cs * s;
        const float sh = p_cs * c + ss * s;
        const float innovation_variance = c * ch + s * sh + measurement_noise;
        const float kc = ch / innovation_variance;
        const float ks = sh / innovation_variance;
        const Eigen::Vector2f residual = Eigen::Vector2f(ax, ay)
            - (c * cosineCoefficients() + s * sineCoefficients());
        coefficient_estimate.head<2>() += kc * residual;
        coefficient_estimate.tail<2>() += ks * residual;

        // Joseph correction of the common 2x2 block, on every sample.
        const float a = 1.0f - kc * c, b = -kc * s;
        const float d = -ks * c, e = 1.0f - ks * s;
        const float ac = a * cc + b * p_cs;
        const float as = a * p_cs + b * ss;
        const float dc = d * cc + e * p_cs;
        const float ds = d * p_cs + e * ss;
        p_cc = ac * a + as * b + kc * measurement_noise * kc;
        p_cs = ac * d + as * e + kc * measurement_noise * ks;
        p_ss = dc * d + ds * e + ks * measurement_noise * ks;

        updateStableAxis(deltaT);
    }

    bool isAxisReliable() const {
        return current_amplitude > AMP_THRESHOLD &&
               confidence > CONFIDENCE_THRESHOLD &&
               current_linearity > LINEARITY_THRESHOLD;
    }

    // Current usable axis. A zero vector explicitly means that the latest
    // motion is not sufficiently axial/confident for propagation-sense use.
    Eigen::Vector2f getAxis() const {
        if (isAxisReliable()) {
            return lastStableDir;
        }
        return Eigen::Vector2f::Zero();
    }

    // Last accepted representative retained for plots/diagnostics even while the
    // current validity gate is closed.
    Eigen::Vector2f getLastStableAxis() const { return lastStableDir; }

    // Current usable angle. NaN explicitly marks an invalid current axis so
    // callers that only consume the scalar angle still fail closed.
    float getAxisDegrees() const {
        return isAxisReliable()
            ? getLastStableAxisDegrees()
            : std::numeric_limits<float>::quiet_NaN();
    }

    float getAxisDegreesRaw() const {
        return isAxisReliable()
            ? wrapAxisDegrees(0.5f * std::atan2(raw_double_axis.y(), raw_double_axis.x())
                              * (180.0f / std::numbers::pi_v<float>))
            : std::numeric_limits<float>::quiet_NaN();
    }

    float getLastStableAxisDegrees() const {
        return wrapAxisDegrees(std::atan2(lastStableDir.y(), lastStableDir.x())
                               * (180.0f / std::numbers::pi_v<float>));
    }

    // Backward-compatible aliases. These return an axis modulo 180 degrees,
    // not a fully directed apparent propagation angle.
    Eigen::Vector2f getDirection() const { return getAxis(); }
    float getDirectionDegrees() const { return getAxisDegrees(); }
    float getDirectionDegreesRaw() const { return getAxisDegreesRaw(); }

    // Axial uncertainty at approximately 95% confidence (2 sigma).
    float getAxisUncertaintyDegrees() const {
        if (!isAxisReliable() || !(lastStableAmplitude > 1e-6f)) {
            return 90.0f;
        }

        // The cosine and sine covariances are isotropic in x/y, so
        // projecting their sum on a unit tangent gives p_cc + p_ss.
        const float angular_std_rad =
            std::sqrt(std::max(0.0f, lastStableVariance)) / lastStableAmplitude;
        const float uncertainty_deg =
            2.0f * angular_std_rad *
            (180.0f / std::numbers::pi_v<float>);
        return std::clamp(uncertainty_deg, 0.0f, 90.0f);
    }

    float getDirectionUncertaintyDegrees() const {
        return getAxisUncertaintyDegrees();
    }

    // Reports whether the retained axis
    // is valid for the current sample rather than holding stale confidence.
    float getLastStableConfidence() const { return confidence; }
    float getHistoricalStableConfidence() const {
        return lastStableConfidence;
    }

    float getAxisLinearity() const { return current_linearity; }
    float getLastStableLinearity() const { return lastStableLinearity; }

    Eigen::Vector2f getFilteredSignal() const {
        return cosineCoefficients() * std::cos(phase) +
               sineCoefficients() * std::sin(phase);
    }

    Eigen::Vector2f getAmplitudeVector() const {
        return getAxis() * current_amplitude;
    }

    float getAmplitude() const { return current_amplitude; }

    Eigen::Vector2f getCosineCoefficients() const {
        return cosineCoefficients();
    }

    Eigen::Vector2f getSineCoefficients() const {
        return sineCoefficients();
    }

    float getPhase() const { return phase; }
    float getConfidence() const { return confidence; }

    void setProcessNoise(float q) {
        if (!std::isfinite(q) || q < 0.0f) return;
        process_noise = q;
    }

    void setMeasurementNoise(float r) {
        if (!std::isfinite(r) || !(r > 0.0f)) return;
        measurement_noise = r;
    }

private:
    static constexpr float REFERENCE_DT_SEC = 1.0f / 200.0f;
    static constexpr float AMP_THRESHOLD = 0.08f;
    static constexpr float CONFIDENCE_THRESHOLD = 20.0f;
    static constexpr float LINEARITY_THRESHOLD = 0.25f;

    static float alphaToTau(float alpha, float reference_dt) {
        if (!(alpha > 0.0f) || alpha >= 1.0f || !std::isfinite(alpha)) {
            return 2.0f;
        }
        return -reference_dt / std::log(1.0f - alpha);
    }

    static float emaAlpha(float dt, float tau) {
        if (!(tau > 0.0f)) return 1.0f;
        return 1.0f - std::exp(-dt / tau);
    }

    static float wrapAxisDegrees(float deg) {
        deg = std::fmod(deg, 180.0f);
        if (deg < 0.0f) deg += 180.0f;
        return deg;
    }

    void invalidateCurrentAxis() {
        confidence = 0.0f;
        current_amplitude = 0.0f;
        current_linearity = 0.0f;
    }

    Eigen::Vector2f cosineCoefficients() const {
        return coefficient_estimate.segment<2>(0);
    }

    Eigen::Vector2f sineCoefficients() const {
        return coefficient_estimate.segment<2>(2);
    }

    void updatePhase(float deltaT) {
        phase = std::remainder(
            phase + omega * deltaT,
            2.0f * std::numbers::pi_v<float>);
    }

    void updateStableAxis(float deltaT) {
        // For a symmetric 2x2 moment matrix, (xx - yy, 2xy) points
        // along twice the principal-axis angle. Smooth it directly.
        const Eigen::Vector2f cosine = cosineCoefficients();
        const Eigen::Vector2f sine = sineCoefficients();
        const float xx = cosine.x() * cosine.x() + sine.x() * sine.x();
        const float yy = cosine.y() * cosine.y() + sine.y() * sine.y();
        const float xy = cosine.x() * cosine.y() + sine.x() * sine.y();
        const float u = xx - yy, v = 2.0f * xy;
        const float gap = std::hypot(u, v);
        const float trace = xx + yy;
        const float large = std::max(0.0f, 0.5f * (trace + gap));
        const float small = std::max(0.0f, 0.5f * (trace - gap));
        current_amplitude = std::sqrt(large);
        const float total = large + small;
        current_linearity = (total > 1e-12f)
            ? std::clamp((large - small) / total, 0.0f, 1.0f) : 0.0f;
        confidence = current_linearity / (2.0f * (p_cc + p_ss) + 1e-6f);
        if (!isAxisReliable()) return;

        raw_double_axis = Eigen::Vector2f(u, v);
        const float target_cos_2theta = u / gap;
        const float target_sin_2theta = v / gap;
        bool near_cancellation = false;
        if (!have_report_axis) {
            report_cos_2theta = target_cos_2theta;
            report_sin_2theta = target_sin_2theta;
            have_report_axis = true;
        } else {
            const float alpha = emaAlpha(deltaT, direction_report_tau_sec);
            report_cos_2theta +=
                alpha * (target_cos_2theta - report_cos_2theta);
            report_sin_2theta +=
                alpha * (target_sin_2theta - report_sin_2theta);
            const float norm =
                std::hypot(report_cos_2theta, report_sin_2theta);
            if (norm > 1e-12f) {
                report_cos_2theta /= norm;
                report_sin_2theta /= norm;
            } else {
                near_cancellation = true;
            }
        }

        // Recover the unit representative with the half-angle identity.
        // Pick the larger component for the division at either axial pole.
        Eigen::Vector2f smoothed_axis;
        if (near_cancellation) {
            const float theta = 0.5f * std::atan2(report_sin_2theta, report_cos_2theta);
            smoothed_axis = Eigen::Vector2f(std::cos(theta), std::sin(theta));
        } else if (report_cos_2theta >= 0.0f) {
            smoothed_axis.x() = std::sqrt(0.5f * (1.0f + report_cos_2theta));
            smoothed_axis.y() = report_sin_2theta / (2.0f * smoothed_axis.x());
        } else {
            smoothed_axis.y() = std::copysign(std::sqrt(0.5f * (1.0f - report_cos_2theta)),
                                              report_sin_2theta);
            smoothed_axis.x() = report_sin_2theta / (2.0f * smoothed_axis.y());
        }
        if (lastStableDir.dot(smoothed_axis) < 0.0f) {
            smoothed_axis = -smoothed_axis;
        }

        lastStableDir = smoothed_axis;
        lastStableAmplitude = current_amplitude;
        lastStableConfidence = confidence;
        lastStableLinearity = current_linearity;
        lastStableVariance = p_cc + p_ss;
    }

    // State ordering: [C_x, C_y, S_x, S_y].
    Eigen::Vector4f coefficient_estimate = Eigen::Vector4f::Zero();
    float p_cc = 1.0f, p_cs = 0.0f, p_ss = 1.0f;
    float process_noise = 1e-6f;
    float measurement_noise = 0.01f;

    float omega = 0.0f;
    float phase = 0.0f;
    float confidence = 0.0f;
    float current_amplitude = 0.0f;
    float current_linearity = 0.0f;

    Eigen::Vector2f lastStableDir = Eigen::Vector2f(1.0f, 0.0f);
    float lastStableAmplitude = 0.0f;
    float lastStableConfidence = 0.0f;
    float lastStableLinearity = 0.0f;
    float lastStableVariance = 2.0f;

    float direction_report_tau_sec = 2.0f;
    bool have_report_axis = false;
    float report_cos_2theta = 1.0f;
    float report_sin_2theta = 0.0f;
    Eigen::Vector2f raw_double_axis = Eigen::Vector2f::UnitX();
};

#ifdef KALMAN_WAVE_DIRECTION_TEST

inline void KalmanWaveDirection_test_signal(
    float t,
    float freq_hz,
    float& ax,
    float& ay,
    std::normal_distribution<float>& noise,
    std::default_random_engine& generator) {
  const float amp = 0.8f + 0.4f * std::sin(0.005f * t);
  const float omega =
      2.0f * std::numbers::pi_v<float> * freq_hz;
  const float carrier = amp * std::cos(omega * t);

  Eigen::Vector2f axis(1.0f, 1.5f);
  axis.normalize();
  ax = carrier * axis.x() + noise(generator);
  ay = carrier * axis.y() + noise(generator);
}

inline void KalmanWaveDirection_test_1() {
  const float delta_t = 0.02f;
  const float freq_hz = 0.5f;
  const float omega =
      2.0f * std::numbers::pi_v<float> * freq_hz;
  const int num_steps = 2000;

  std::default_random_engine generator;
  generator.seed(42u);
  std::normal_distribution<float> dist(0.0f, 0.08f);

  KalmanWaveDirection filter(omega);
  filter.setMeasurementNoise(0.01f);
  filter.setProcessNoise(1e-6f);

  std::ofstream out("wave_dir.csv");
  out << "t,ax,ay,filtered_ax,filtered_ay,freq_hz,amplitude,phase,"
         "confidence,deg,uncertaintyDeg\n";

  for (int i = 0; i < num_steps; ++i) {
    const float t = i * delta_t;
    float ax = 0.0f;
    float ay = 0.0f;
    KalmanWaveDirection_test_signal(
        t, freq_hz, ax, ay, dist, generator);
    filter.update(ax, ay, omega, delta_t);

    const Eigen::Vector2f filtered = filter.getFilteredSignal();
    out << t << "," << ax << "," << ay << ","
        << filtered.x() << "," << filtered.y() << ","
        << freq_hz << "," << filter.getAmplitude() << ","
        << filter.getPhase() << "," << filter.getConfidence() << ","
        << filter.getDirectionDegrees() << ","
        << filter.getDirectionUncertaintyDegrees() << "\n";
  }
}

#endif  // KALMAN_WAVE_DIRECTION_TEST

#endif  // KALMAN_WAVE_DIRECTION_H
