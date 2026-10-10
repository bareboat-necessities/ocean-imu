#pragma once

// Copyright 2026, Mikhail Grushinskiy
// Relative gyro rotations only: no attitude/reference/correction state.
// See docs/magnetic-timing.md for clocks, covariance and correlation limits.
#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>

namespace ocean_imu::magnetic {
using V3 = Eigen::Vector3f;
using M3 = Eigen::Matrix3f;
using Q = Eigen::Quaternionf;

inline Q increment(const V3& omega, float dt) {
    const V3 v = -omega * dt; // shipping qref is W->B; LEFT increment
    const float a2 = v.squaredNorm();
    const float a = std::sqrt(a2);
    const float s = a < 1e-4f ? 0.5f - a2 / 48.0f : std::sin(0.5f * a) / a;
    Q q(std::cos(0.5f * a), s * v.x(), s * v.y(), s * v.z());
    q.normalize();
    return q;
}

enum class Status : uint8_t { Ready, First, Duplicate, OutOfOrder, Invalid, Gap, NoHistory, TooOld, Future };

// micros() wrap is supported; comparisons require intervals below 2^31 us.
// 32 slots cover 160 ms at 200 Hz; queries are additionally capped at 120 ms.
template <unsigned Capacity = 32>
class RotationHistory {
    static_assert(Capacity > 1, "rotation history needs at least two slots");
public:
    static constexpr uint32_t MaxAgeUs = 120000;
    static constexpr uint32_t MaxGapUs = 20000;
    struct Segment {
        uint32_t begin = 0, end = 0;
        V3 omega = V3::Zero(); // calibrated physical-body rate minus residual bias
        float bias_variance = 0;
        Q dq = Q::Identity();
    };
    struct Rotation {
        Q q = Q::Identity(); // B_to <- B_from
        float seconds = 0;
        float start_speed = 0, end_speed = 0;
        float bias_variance = 0;
        unsigned visited = 0;
    };

    void reset() { count_ = next_ = 0; have_time_ = false; }
    Status push(uint32_t t, const V3& omega, float bias_variance = 0) {
        if (!omega.allFinite() || !std::isfinite(omega.squaredNorm()) ||
            !std::isfinite(bias_variance) || bias_variance < 0) {
            reset(); return Status::Invalid;
        }
        if (!have_time_) { last_ = t; have_time_ = true; return Status::First; }
        const int32_t delta = static_cast<int32_t>(t - last_);
        if (delta == 0) return Status::Duplicate;
        if (delta < 0) return Status::OutOfOrder;
        if (uint32_t(delta) > MaxGapUs) {
            count_ = next_ = 0; last_ = t; return Status::Gap;
        }
        const Q dq = increment(omega, float(delta) * 1e-6f);
        if (!dq.coeffs().allFinite()) { reset(); return Status::Invalid; }
        data_[next_] = Segment{last_, t, omega, bias_variance, dq};
        next_ = (next_ + 1) % Capacity;
        count_ = std::min(count_ + 1, Capacity);
        last_ = t;
        return Status::Ready;
    }

    Status between(uint32_t from, uint32_t to, Rotation& out) const {
        out = Rotation{};
        if (!have_time_) return Status::NoHistory;
        const int32_t age = static_cast<int32_t>(last_ - from);
        const int32_t end_age = static_cast<int32_t>(last_ - to);
        if (age < 0 || end_age < 0) return Status::Future;
        if (end_age > age) return Status::OutOfOrder;
        if (uint32_t(age) > MaxAgeUs) return Status::TooOld;
        if (!count_) return from == to && to == last_ ? Status::Ready : Status::NoHistory;
        const unsigned oldest = (next_ + Capacity - count_) % Capacity;
        if (static_cast<int32_t>(from - data_[oldest].begin) < 0) return Status::NoHistory;
        uint32_t covered = 0;
        for (unsigned j = 0; j < count_; ++j) {
            const auto& s = data_[(oldest + j) % Capacity];
            ++out.visited;
            if (from == to && static_cast<int32_t>(from-s.begin) >= 0 &&
                static_cast<int32_t>(s.end-from) >= 0) {
                out.start_speed = out.end_speed = s.omega.norm();
            }
            const int32_t start = std::max(int32_t(0), static_cast<int32_t>(from - s.begin));
            const int32_t stop = std::min(static_cast<int32_t>(s.end - s.begin),
                                          static_cast<int32_t>(to - s.begin));
            if (stop <= start) continue;
            const uint32_t us = uint32_t(stop - start);
            const Q dq = (start == 0 && stop == static_cast<int32_t>(s.end - s.begin))
                ? s.dq : increment(s.omega, float(us) * 1e-6f);
            out.q = dq * out.q;
            if (!covered) out.start_speed = s.omega.norm();
            out.end_speed = s.omega.norm();
            out.bias_variance = std::max(out.bias_variance, s.bias_variance);
            covered += us;
        }
        if (covered != uint32_t(to - from)) return Status::NoHistory;
        out.q.normalize();
        out.seconds = float(covered) * 1e-6f;
        return out.q.coeffs().allFinite() ? Status::Ready : Status::Invalid;
    }
private:
    std::array<Segment, Capacity> data_{};
    unsigned count_ = 0, next_ = 0;
    uint32_t last_ = 0;
    bool have_time_ = false;
};

struct Uncertainty {
    // Scalar upper bounds on eigenvalues, not per-sample gyro standard deviation.
    float gyro_density = 0;          // rad/sqrt(s)
    float bias_variance = 0;         // (rad/s)^2, held/common over the interval
    float from_time_sigma = 0, to_time_sigma = 0; // seconds; may be conservative bounds
    bool correlated = true;          // gyro-bias estimate can depend on magnetic samples
    bool valid() const {
        return std::isfinite(gyro_density) && gyro_density >= 0 &&
            std::isfinite(bias_variance) && bias_variance >= 0 &&
            std::isfinite(from_time_sigma) && from_time_sigma >= 0 &&
            std::isfinite(to_time_sigma) && to_time_sigma >= 0;
    }
};

inline float covarianceBound(const M3& p) {
    return p.cwiseAbs().rowwise().sum().maxCoeff(); // Gershgorin/spectral upper bound
}

template <class Rotation>
inline M3 rotationNoise(const V3& transported, const Rotation& r, const Uncertainty& u) {
    const float time_angle = r.start_speed * u.from_time_sigma + r.end_speed * u.to_time_sigma;
    const float gyro_sd = u.gyro_density * std::sqrt(r.seconds);
    const float bias_sd = std::sqrt(std::max(u.bias_variance, r.bias_variance)) * r.seconds;
    const float variance = u.correlated ? (gyro_sd + bias_sd + time_angle) * (gyro_sd + bias_sd + time_angle)
        : gyro_sd * gyro_sd + bias_sd * bias_sd + time_angle * time_angle;
    return variance * (transported.squaredNorm() * M3::Identity() - transported * transported.transpose());
}

inline bool covarianceValid(const M3& r) {
    return r.allFinite() && r.isApprox(r.transpose(), 1e-5f) && Eigen::LLT<M3>(r).info() == Eigen::Success;
}

template <class Rotation>
inline bool align(const V3& calibrated, const V3& body_hard_iron, const M3& sample_covariance,
                  const Rotation& r, const Uncertainty& u, V3& field, M3& covariance) {
    if (!calibrated.allFinite() || !body_hard_iron.allFinite() || !u.valid() ||
        !covarianceValid(sample_covariance) || !r.q.coeffs().allFinite() ||
        !std::isfinite(r.seconds) || r.seconds < 0 || r.seconds > .120001f ||
        !std::isfinite(r.start_speed) || r.start_speed < 0 ||
        !std::isfinite(r.end_speed) || r.end_speed < 0 ||
        !std::isfinite(r.bias_variance) || r.bias_variance < 0 ||
        std::abs(r.q.squaredNorm() - 1.0f) > 1e-3f) return false;
    field = r.q * (calibrated - body_hard_iron); // subtract body-fixed offset FIRST
    const M3 d = r.q.toRotationMatrix();
    covariance = d * sample_covariance * d.transpose() + rotationNoise(field, r, u);
    // Cov(x+y) <= 2 Cov(x)+2 Cov(y) for unknown cross-correlation.
    if (u.correlated && (r.seconds > 0 || u.from_time_sigma > 0 || u.to_time_sigma > 0)) covariance *= 2.0f;
    return field.allFinite() && covarianceValid(covariance);
}

struct Consistency {
    Status status = Status::First;
    float residual_uT = NAN, nis = NAN;
    float predicted_angle_rad = NAN, observed_angle_rad = NAN;
    float rotation_signal_to_noise = NAN;
    bool unexplained = false; // advisory ONLY, never estimator acceptance
};

// Attitude/IMU snapshots belong to the observation epoch. In particular the
// proxy snapshot must be exogenous: never rewind it with MEKF bias estimates.
struct Observation {
    RotationHistory<>::Rotation rotation{};
    Uncertainty uncertainty{};
    M3 covariance = M3::Identity();
    Q proxy_bw = Q::Identity(), attitude_bw = Q::Identity();
    V3 accel = V3::Zero(), gyro = V3::Zero();
    bool valid() const {
        return uncertainty.valid() && covarianceValid(covariance) &&
            proxy_bw.coeffs().allFinite() && attitude_bw.coeffs().allFinite() &&
            std::abs(proxy_bw.squaredNorm()-1.0f)<1e-3f &&
            std::abs(attitude_bw.squaredNorm()-1.0f)<1e-3f &&
            accel.allFinite() && gyro.allFinite() &&
            std::isfinite(rotation.seconds) && rotation.seconds>=0 && rotation.seconds<=.120001f &&
            std::isfinite(rotation.start_speed) && rotation.start_speed>=0 &&
            std::isfinite(rotation.end_speed) && rotation.end_speed>=0 &&
            std::isfinite(rotation.bias_variance) && rotation.bias_variance>=0 &&
            rotation.q.coeffs().allFinite() && std::abs(rotation.q.squaredNorm()-1.0f)<1e-3f;
    }
};

inline float angle(const V3& a, const V3& b) {
    const float n = a.norm() * b.norm();
    return n > 1e-12f ? std::acos(std::clamp(a.dot(b) / n, -1.0f, 1.0f)) : NAN;
}

class RotationConsistency {
public:
    void reset() { have_ = false; }
    template <unsigned N>
    Consistency observe(const RotationHistory<N>& history, uint32_t sequence, uint32_t time,
                        const V3& calibrated, const V3& body_hard_iron,
                        const M3& covariance, const Uncertainty& u) {
        Consistency c;
        if (!calibrated.allFinite() || calibrated.norm() < 1e-6f ||
            !body_hard_iron.allFinite() || !u.valid() || !covarianceValid(covariance)) {
            c.status = Status::Invalid; return c;
        }
        if (have_ && sequence == sequence_) { c.status = Status::Duplicate; return c; }
        if (have_ && (static_cast<int32_t>(sequence - sequence_) < 0 ||
                      static_cast<int32_t>(time - time_) <= 0)) {
            c.status = Status::OutOfOrder; return c;
        }
        if (have_) {
            typename RotationHistory<N>::Rotation r;
            c.status = history.between(time_, time, r);
            if (c.status == Status::Ready) {
                // Evaluate BOTH samples at the current applied offset. Do not
                // mistake a changed calibration parameter for field motion.
                const V3 previous = previous_ - body_hard_iron;
                const V3 current = calibrated - body_hard_iron;
                const V3 prediction = r.q * previous;
                const V3 residual = current - prediction;
                const M3 d = r.q.toRotationMatrix();
                M3 sigma = covariance + d * covariance_ * d.transpose() + rotationNoise(prediction, r, u);
                // Three possibly correlated terms: current mag, previous mag,
                // and gyro transport. This is a marginal upper bound; adjacent
                // residuals still share samples and must NOT be pooled as iid.
                if (u.correlated) sigma *= 3.0f;
                const Eigen::LLT<M3> llt(sigma);
                if (llt.info() == Eigen::Success) {
                    c.residual_uT = residual.norm();
                    c.nis = residual.dot(llt.solve(residual));
                    c.predicted_angle_rad = angle(previous, prediction);
                    c.observed_angle_rad = angle(previous, current);
                    c.rotation_signal_to_noise = (prediction - previous).norm() / std::sqrt(sigma.trace());
                    // chi-square(3), 0.999. A reference level, not a calibrated
                    // false-alarm probability for bounded/timestamp/model errors.
                    c.unexplained = c.nis > 16.266236f;
                } else c.status = Status::Invalid;
            }
        }
        previous_ = calibrated; covariance_ = covariance;
        time_ = time; sequence_ = sequence; have_ = true;
        return c;
    }
private:
    V3 previous_ = V3::Zero();
    M3 covariance_ = M3::Identity();
    uint32_t time_ = 0, sequence_ = 0;
    bool have_ = false;
};
} // namespace ocean_imu::magnetic
