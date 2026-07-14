#ifndef KALMAN_WAVE_DIRECTION_H
#define KALMAN_WAVE_DIRECTION_H

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
        P = Eigen::Matrix4f::Identity();
        phase = 0.0f;
        invalidateCurrentAxis();

        lastStableAmplitude = 0.0f;
        lastStableConfidence = 0.0f;
        lastStableLinearity = 0.0f;
        lastStableCovariance = Eigen::Matrix4f::Identity();
        lastStableDir = Eigen::Vector2f(1.0f, 0.0f);

        have_report_axis = false;
        report_cos_2theta = 1.0f;
        report_sin_2theta = 0.0f;
        direction_deg_raw = 0.0f;
        direction_deg_smoothed = 0.0f;
    }

    void update(float ax, float ay, float currentOmega, float deltaT) {
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
        // the carrier-phase singularity of the previous cosine-only model.
        Eigen::Matrix<float, 2, 4> H;
        H << c, 0.0f, s, 0.0f,
             0.0f, c, 0.0f, s;

        const float q_scale = deltaT / REFERENCE_DT_SEC;
        Eigen::Matrix4f P_pred =
            P + std::max(0.0f, q_scale) * Q_reference_sample;
        P_pred = 0.5f * (P_pred + P_pred.transpose());
        P_pred += Eigen::Matrix4f::Identity() * 1e-10f;

        const Eigen::Vector2f z(ax, ay);
        const Eigen::Matrix2f innovation_covariance =
            H * P_pred * H.transpose() + R;
        const Eigen::Matrix<float, 4, 2> K =
            P_pred * H.transpose() *
            innovation_covariance.ldlt().solve(Eigen::Matrix2f::Identity());

        coefficient_estimate += K * (z - H * coefficient_estimate);

        const Eigen::Matrix4f I = Eigen::Matrix4f::Identity();
        const Eigen::Matrix4f KH = K * H;
        P = (I - KH) * P_pred * (I - KH).transpose() +
            K * R * K.transpose();
        P = 0.5f * (P + P.transpose());

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
        return isAxisReliable() ? lastStableDir : Eigen::Vector2f::Zero();
    }

    // Historical representative retained for plots/diagnostics even while the
    // current validity gate is closed.
    Eigen::Vector2f getLastStableAxis() const { return lastStableDir; }

    float getAxisDegrees() const { return direction_deg_smoothed; }
    float getAxisDegreesRaw() const { return direction_deg_raw; }

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

        const Eigen::Vector2f tangent(-lastStableDir.y(), lastStableDir.x());
        const Eigen::Matrix2f P_cos =
            lastStableCovariance.block<2, 2>(0, 0);
        const Eigen::Matrix2f P_sin =
            lastStableCovariance.block<2, 2>(2, 2);
        float tangent_variance =
            tangent.dot((P_cos + P_sin) * tangent);
        tangent_variance = std::max(0.0f, tangent_variance);

        const float angular_std_rad =
            std::sqrt(tangent_variance) / lastStableAmplitude;
        const float uncertainty_deg =
            2.0f * angular_std_rad *
            (180.0f / std::numbers::pi_v<float>);
        return std::clamp(uncertainty_deg, 0.0f, 90.0f);
    }

    float getDirectionUncertaintyDegrees() const {
        return getAxisUncertaintyDegrees();
    }

    // Kept under the original API name, but reports whether the retained axis
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
        Q_reference_sample = Eigen::Matrix4f::Identity() * q;
    }

    void setMeasurementNoise(float r) {
        if (!std::isfinite(r) || !(r > 0.0f)) return;
        R = Eigen::Matrix2f::Identity() * r;
    }

private:
    struct AxisMetrics {
        Eigen::Vector2f representative = Eigen::Vector2f(1.0f, 0.0f);
        float amplitude = 0.0f;
        float linearity = 0.0f;
    };

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

    AxisMetrics computeAxisMetrics() const {
        const Eigen::Vector2f c = cosineCoefficients();
        const Eigen::Vector2f s = sineCoefficients();
        const Eigen::Matrix2f moment =
            c * c.transpose() + s * s.transpose();

        const float trace = moment.trace();
        const float discriminant = std::hypot(
            moment(0, 0) - moment(1, 1),
            2.0f * moment(0, 1));
        const float lambda_max =
            std::max(0.0f, 0.5f * (trace + discriminant));
        const float lambda_min =
            std::max(0.0f, 0.5f * (trace - discriminant));

        const float theta = 0.5f * std::atan2(
            2.0f * moment(0, 1),
            moment(0, 0) - moment(1, 1));

        AxisMetrics metrics;
        metrics.representative =
            Eigen::Vector2f(std::cos(theta), std::sin(theta));
        metrics.amplitude = std::sqrt(lambda_max);
        const float total = lambda_max + lambda_min;
        metrics.linearity = (total > 1e-12f)
            ? std::clamp((lambda_max - lambda_min) / total, 0.0f, 1.0f)
            : 0.0f;
        return metrics;
    }

    void updatePhase(float deltaT) {
        phase = std::remainder(
            phase + omega * deltaT,
            2.0f * std::numbers::pi_v<float>);
    }

    void updateStableAxis(float deltaT) {
        const AxisMetrics metrics = computeAxisMetrics();
        current_amplitude = metrics.amplitude;
        current_linearity = metrics.linearity;

        const float covariance_confidence =
            1.0f / (P.trace() + 1e-6f);
        confidence = covariance_confidence * current_linearity;

        Eigen::Vector2f new_axis = metrics.representative;
        if (lastStableDir.dot(new_axis) < 0.0f) {
            new_axis = -new_axis;
        }
        direction_deg_raw = wrapAxisDegrees(
            std::atan2(new_axis.y(), new_axis.x()) *
            (180.0f / std::numbers::pi_v<float>));

        if (!isAxisReliable()) {
            return;
        }

        const float theta = std::atan2(new_axis.y(), new_axis.x());
        const float target_cos_2theta = std::cos(2.0f * theta);
        const float target_sin_2theta = std::sin(2.0f * theta);

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
            }
        }

        const float smoothed_theta = 0.5f * std::atan2(
            report_sin_2theta, report_cos_2theta);
        Eigen::Vector2f smoothed_axis(
            std::cos(smoothed_theta), std::sin(smoothed_theta));
        if (lastStableDir.dot(smoothed_axis) < 0.0f) {
            smoothed_axis = -smoothed_axis;
        }

        lastStableDir = smoothed_axis;
        direction_deg_smoothed = wrapAxisDegrees(
            std::atan2(lastStableDir.y(), lastStableDir.x()) *
            (180.0f / std::numbers::pi_v<float>));

        lastStableAmplitude = current_amplitude;
        lastStableConfidence = confidence;
        lastStableLinearity = current_linearity;
        lastStableCovariance = P;
    }

    // State ordering: [C_x, C_y, S_x, S_y].
    Eigen::Vector4f coefficient_estimate = Eigen::Vector4f::Zero();
    Eigen::Matrix4f P = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f Q_reference_sample =
        Eigen::Matrix4f::Identity() * 1e-6f;
    Eigen::Matrix2f R = Eigen::Matrix2f::Identity() * 0.01f;

    float omega = 0.0f;
    float phase = 0.0f;
    float confidence = 0.0f;
    float current_amplitude = 0.0f;
    float current_linearity = 0.0f;

    Eigen::Vector2f lastStableDir = Eigen::Vector2f(1.0f, 0.0f);
    float lastStableAmplitude = 0.0f;
    float lastStableConfidence = 0.0f;
    float lastStableLinearity = 0.0f;
    Eigen::Matrix4f lastStableCovariance =
        Eigen::Matrix4f::Identity();

    float direction_report_tau_sec = 2.0f;
    bool have_report_axis = false;
    float report_cos_2theta = 1.0f;
    float report_sin_2theta = 0.0f;
    float direction_deg_raw = 0.0f;
    float direction_deg_smoothed = 0.0f;
};

#ifdef KALMAN_WAVE_DIRECTION_TEST

void KalmanWaveDirection_test_signal(
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

void KalmanWaveDirection_test_1() {
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
