#pragma once

// Copyright 2026, Mikhail Grushinskiy
// Relative gyro rotations only: no attitude/reference/correction state.
// See docs/magnetic-timing.md for clocks and history limits.
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
template <unsigned Capacity = 32, uint32_t QueryAgeUs = 120000>
class RotationHistory {
    static_assert(Capacity > 1, "rotation history needs at least two slots");
public:
    static constexpr uint32_t MaxAgeUs = QueryAgeUs;
    static constexpr uint32_t MaxGapUs = 20000;
    struct Segment {
        uint32_t begin = 0, end = 0;
        V3 omega = V3::Zero(); // calibrated physical-body rate minus residual bias
        float bias_variance = 0;
        float gyro_angle_variance = 0; // optional discrete angle-noise variance
        Q dq = Q::Identity();
    };
    struct Rotation {
        Q q = Q::Identity(); // B_to <- B_from
        float seconds = 0;
        float start_speed = 0, end_speed = 0;
        float bias_variance = 0;
        unsigned visited = 0;
        float gyro_angle_variance = 0;
    };

    void reset() { count_ = next_ = 0; have_time_ = false; }
    Status push(uint32_t t, const V3& omega, float bias_variance = 0, float gyro_angle_variance = 0) {
        if (!omega.allFinite() || !std::isfinite(omega.squaredNorm()) ||
            !std::isfinite(bias_variance) || bias_variance < 0 ||
            !std::isfinite(gyro_angle_variance) || gyro_angle_variance < 0) {
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
        data_[next_] = Segment{last_, t, omega, bias_variance, gyro_angle_variance, dq};
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
            if (s.gyro_angle_variance > 0) {
                const float fraction = us == s.end-s.begin ? 1.0f : float(us) / float(s.end-s.begin);
                out.gyro_angle_variance += fraction * fraction * s.gyro_angle_variance;
            }
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

} // namespace ocean_imu::magnetic
