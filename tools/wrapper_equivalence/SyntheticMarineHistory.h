#pragma once

// Deterministic synthetic marine IMU history for wrapper characterization.
//
// One continuous 200 Hz record of a moving hull: roll/pitch in the wave band,
// a slow heading sweep (so a hard-iron offset is identifiable), orbital
// acceleration from three wave components, gyro/accelerometer bias, a fixed
// hard-iron offset, white noise from a fixed-seed generator, and scheduled
// events that exercise the orchestration branches the characterization
// compares:
//
//   * an engine-vibration interval (accelerometer vibration guard and the
//     vibration-aware accelerometer covariance),
//   * a quiet-water interval (stillness detector and sigma attenuation),
//   * an optional capsize interval (tilt watchdog reset while Live).
//
// Everything is plain float/double arithmetic with no platform RNG, so the
// same history is produced by any build of the same compiler flags.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Geometry>
#endif

namespace characterization {

struct Lcg {
    uint64_t s;
    explicit Lcg(uint64_t seed) : s(seed) {}
    double uniform() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return (double)((s >> 11) + 1) / 9007199254740994.0;
    }
    double gauss() {
        const double u1 = uniform(), u2 = uniform();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    }
};

struct HistoryConfig {
    double duration_sec   = 420.0;
    double dt             = 0.005;
    double gravity        = 9.80665;
    bool   vibration      = true;   // 200..240 s engine interval
    bool   quiet_interval = true;   // 300..330 s still water
    bool   capsize        = false;  // rolled past the tilt watchdog for 3 s
    double capsize_start  = 360.0;
    uint64_t seed         = 0x5eed0fu;
};

struct Sample {
    double t = 0.0;
    Eigen::Vector3f gyro;
    Eigen::Vector3f acc;
    bool mag_due = false;
    Eigen::Vector3f mag;
};

class SyntheticMarineHistory {
public:
    explicit SyntheticMarineHistory(const HistoryConfig& c) : cfg_(c), rng_(c.seed) {
        q_prev_ = attitude_(0.0);
    }

    bool next(Sample& out) {
        if (k_ * cfg_.dt > cfg_.duration_sec) return false;
        const double t = (k_ + 1) * cfg_.dt;
        const Eigen::Quaterniond q = attitude_(t);
        // Body rate from consecutive attitudes: q_{k+1} = q_k * exp(w dt / 2).
        const Eigen::Quaterniond dq = q_prev_.conjugate() * q;
        Eigen::AngleAxisd aa(dq);
        Eigen::Vector3d w = aa.axis() * aa.angle() / cfg_.dt;
        if (!(aa.angle() > 0.0)) w.setZero();
        q_prev_ = q;

        const double amp = amplitude_(t);
        Eigen::Vector3d a_world = Eigen::Vector3d::Zero();
        const double comps[3][3] = {{0.9, 0.125, 0.3}, {0.5, 0.21, 1.4}, {0.25, 0.33, 2.2}};
        for (const auto& c : comps) {
            const double om = 6.283185307179586 * c[1];
            const double a = -amp * c[0] * om * om;
            a_world.z() += a * std::sin(om * t + c[2]);
            a_world.x() += 0.6 * a * std::cos(om * t + c[2]);
            a_world.y() += 0.35 * a * std::cos(om * t + c[2] + 0.4);
        }
        const Eigen::Vector3d g_world(0.0, 0.0, cfg_.gravity);
        Eigen::Vector3d f_body = q.conjugate() * (a_world - g_world);

        if (cfg_.vibration && t >= 200.0 && t < 240.0) {
            f_body.x() += 0.35 * std::sin(6.283185307179586 * 27.0 * t);
            f_body.z() += 0.60 * std::sin(6.283185307179586 * 31.0 * t + 0.3);
        }

        const Eigen::Vector3d gyro_bias(0.004, -0.003, 0.002);
        const Eigen::Vector3d acc_bias(0.03, -0.02, 0.04);
        out.t = t;
        out.gyro = (w + gyro_bias + 0.002 * noise3_()).cast<float>();
        out.acc = (f_body + acc_bias + 0.015 * noise3_()).cast<float>();

        out.mag_due = (k_ % 8) == 0;  // 25 Hz
        if (out.mag_due) {
            const Eigen::Vector3d B_world(20.0, 1.5, 44.0);
            const Eigen::Vector3d hard_iron(4.0, -2.5, 1.5);
            out.mag = (q.conjugate() * B_world + hard_iron + 0.3 * noise3_()).cast<float>();
        }
        ++k_;
        return true;
    }

private:
    Eigen::Vector3d noise3_() { return {rng_.gauss(), rng_.gauss(), rng_.gauss()}; }

    double amplitude_(double t) const {
        if (cfg_.quiet_interval && t >= 300.0 && t < 330.0) return 0.02;
        return 1.0 + 0.4 * std::sin(6.283185307179586 * t / 180.0);
    }

    Eigen::Quaterniond attitude_(double t) const {
        const double amp = amplitude_(t);
        double roll  = amp * 0.10 * std::sin(6.283185307179586 * 0.14 * t + 0.2)
                     + 0.03 * std::sin(6.283185307179586 * 0.31 * t);
        const double pitch = amp * 0.06 * std::sin(6.283185307179586 * 0.19 * t + 1.1);
        const double yaw   = 0.4 + 0.08 * std::sin(6.283185307179586 * t / 150.0);
        if (cfg_.capsize && t >= cfg_.capsize_start && t < cfg_.capsize_start + 3.0) {
            const double s = (t - cfg_.capsize_start) / 3.0;
            roll += 1.8 * std::sin(3.141592653589793 * s);
        }
        return Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ()) *
               Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY()) *
               Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX());
    }

    HistoryConfig cfg_;
    Lcg rng_;
    Eigen::Quaterniond q_prev_;
    uint64_t k_ = 0;
};

// Binary trace: every emitted value is written as its raw float bits, so two
// traces compare equal only when every value is bit-identical.
struct Trace {
    std::FILE* fp = nullptr;
    uint64_t hash = 1469598103934665603ULL;
    uint64_t values = 0;
    void put(float v) {
        uint32_t u;
        std::memcpy(&u, &v, 4);
        if (fp) std::fwrite(&u, 4, 1, fp);
        for (int i = 0; i < 4; ++i) { hash ^= (u >> (8 * i)) & 0xffu; hash *= 1099511628211ULL; }
        ++values;
    }
    void put(double v) { put(static_cast<float>(v)); }
    void put(bool v) { put(v ? 1.0f : 0.0f); }
    void put(int v) { put(static_cast<float>(v)); }
    template <typename Derived>
    void put(const Eigen::MatrixBase<Derived>& m) {
        for (int i = 0; i < m.size(); ++i) put(static_cast<float>(m(i)));
    }
    void put(const Eigen::Quaternionf& q) { put(q.coeffs()); }
};

}  // namespace characterization
