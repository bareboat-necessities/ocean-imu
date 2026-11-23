#pragma once
#include "util/HotLoop.h"
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include "tuner/VerticalAccelComplementary.h"

namespace seastate::common {

// Separate yaw-induced acceleration of steady forward travel from wave motion.
// Learns speed from IMU transport jerk: j_b = df_b/dt + omega_b x f_b.
// Gravity cancels in that derivative. In the level bow frame a steady speed U
// gives j_h = U * (-yaw_rate^2, yaw_acceleration). Identical low-pass filters
// on the regressor and observation suppress wave jerk without notching waves
// in the output. No position, attitude truth, speed or steering period is read.
class MarineTurnSeparator {
public:
    struct Config {
        bool enabled = true;
        // Sensor/body coordinates, after mounting calibration; not world north.
        Eigen::Vector3f forward_body = Eigen::Vector3f::UnitX();
    };
    void reset(Config config) {
        *this = MarineTurnSeparator();
        // Take Config by value so reset(cfg_) remains safe when resetting the
        // object invalidates the member that supplied the argument.
        cfg_ = std::move(config);
        if (!cfg_.forward_body.allFinite() || cfg_.forward_body.norm() < 0.5f)
            cfg_.enabled = false;
        else cfg_.forward_body.normalize();
    }
    OCEAN_IMU_SAMPLE_IRAM Eigen::Vector3f update(float dt, const Eigen::Vector3f& gyro,
                          const Eigen::Vector3f& force, float gravity) {
        if (!cfg_.enabled) return force;
        if (!(dt > 0 && dt <= 0.1f) || !std::isfinite(dt) ||
            !gyro.allFinite() || !force.allFinite()) {
            reset(cfg_);
            return force;
        }
        proxy_.update(dt, gyro, force, gravity);
        // Express the levelled bow axes in the body frame. Rotating the
        // observer around world down cancels from every dot product below.
        const Eigen::Vector3f down = proxy_.quaternion().conjugate() * Eigen::Vector3f::UnitZ();
        const float bow_down = down.dot(cfg_.forward_body);
        const float horizontal2 = 1.0f - bow_down*bow_down;
        if (horizontal2 < 0.25f) {
            reset(cfg_);
            return force;
        }
        const float inv_horizontal = 1.0f/std::sqrt(horizontal2);
        const Eigen::Vector3f forward = (cfg_.forward_body - bow_down*down)*inv_horizontal;
        const Eigen::Vector3f starboard = down.cross(cfg_.forward_body)*inv_horizontal;
        const float rate = (gyro.dot(down) - bow_down*gyro.dot(cfg_.forward_body))/horizontal2;
        elapsed_ += dt;
        if (!initialized_) {
            previous_force_ = force;
            previous_gyro_ = gyro;
            previous_rate_ = rate;
            initialized_ = true;
            return force;
        }
        const Eigen::Vector3f jerk = (force-previous_force_)/dt +
            (0.5f*(gyro+previous_gyro_)).cross(0.5f*(force+previous_force_));
        std::array<float,5> sample{{-rate*rate,
            (rate-previous_rate_)/dt, forward.dot(jerk), starboard.dot(jerk), rate}};
        previous_force_ = force;
        previous_gyro_ = gyro;
        previous_rate_ = rate;
        // At dt <= 0.1 s the prewarp argument is at most 0.011 radians.
        // Dropping terms after x^3/3 in tan(x) has relative error < 2e-9.
        const float x = (3.14159265358979323846f * 0.035f) * dt;
        const float k = x * (1.0f + x*x*(1.0f/3.0f));
        const float k2 = k*k;
        for (std::size_t section=0; section<2; ++section) {
            const float inv_q = section == 0 ? 1.847759065f : 0.765366865f;
            const float damping = k*inv_q;
            const float norm = 1.0f/(1.0f+damping+k2);
            const float b = k2*norm;
            const float d = 2.0f*(damping+2.0f*k2)*norm;
            // Transposed direct-form II, in the fixed coordinates z=s0,
            // v=s0+s1. Computing d=2+a1 and 4b=1+a1+a2 directly avoids
            // cancellation near the unit pole at 200 Hz. The transform is
            // independent of dt, so a jittered sample uses the same states.
            for (std::size_t i=0; i<sample.size(); ++i) {
                auto& state=lowpass_[section][i];
                const float input=sample[i];
                const float feed=b*input;
                const float y=state.z+feed;
                state.z += (state.v-d*y)+4.0f*feed;
                state.v += (4.0f*b)*(input-y);
                sample[i]=y;
            }
        }
        // Ignore startup derivative/filter transients. Centered regression
        // absorbs the constant jerk from residual gyro/accelerometer bias.
        if (elapsed_ > 60.0f) {
            if (!statistics_initialized_) { mean_=sample; statistics_initialized_=true; }
            const float a=decayGain(dt/120.0f);
            std::array<float,5> delta{};
            for (std::size_t i=0; i<sample.size(); ++i) {
                delta[i]=sample[i]-mean_[i]; mean_[i]+=a*delta[i];
            }
            const float retain=1.0f-a;
            xx_+=a*(retain*(delta[0]*delta[0]+delta[1]*delta[1])-xx_);
            xy_+=a*(retain*(delta[0]*delta[2]+delta[1]*delta[3])-xy_);
            yy_+=a*(retain*(delta[2]*delta[2]+delta[3]*delta[3])-yy_);
            rr_+=a*(retain*delta[4]*delta[4]-rr_);
            float target=0;
            // Constant heading gives no information about cruise speed. Do
            // not divide by a noise-only angular-acceleration regressor.
            if (elapsed_ > 120 && xx_ > 1e-9f && rr_ > 1e-7f) {
                const float speed=xy_/xx_;
                if (std::abs(speed) >= 0.2f && std::abs(speed) <= 15 &&
                    xy_*xy_ > 0.015f*xx_*yy_) target=speed;
            }
            speed_ += decayGain(dt/15.0f)*(target-speed_);
        }
        correction_ = starboard*(speed_*rate);
        return force-correction_;
    }
    // Apply in the estimator's current gravity plane. The nuisance speed is
    // learned independently, but a wave-only horizontal correction must have
    // zero vertical component in the navigation frame used by that estimator.
    Eigen::Vector3f correctionBodyMs2(const Eigen::Vector3f& down_body,
                                     const Eigen::Vector3f& gyro) const {
        if (!down_body.allFinite() || !gyro.allFinite() || down_body.norm() < 1e-3f)
            return Eigen::Vector3f::Zero();
        const Eigen::Vector3f down=down_body.normalized();
        const float c=down.dot(cfg_.forward_body);
        const float horizontal2=1-c*c;
        if (!(horizontal2 > 0.25f)) return Eigen::Vector3f::Zero();
        const float rate=(gyro.dot(down)-gyro.dot(cfg_.forward_body)*c)/horizontal2;
        return down.cross(cfg_.forward_body)*(speed_*rate/std::sqrt(horizontal2));
    }
    float speedMps() const { return speed_; }
    const Eigen::Vector3f& correctionBodyMs2() const { return correction_; }
private:
    // 1-exp(-x), with 0 < x <= 0.1/15. Relative truncation error < 1.3e-8.
    static float decayGain(float x) {
        return x*(1.0f+x*(-0.5f+x*(1.0f/6.0f)));
    }

    struct SectionState {
        float z=0, v=0;
    };
    Config cfg_{};
    VerticalAccelComplementary proxy_{0.2f,0.02f,0};
    bool initialized_=false, statistics_initialized_=false;
    float elapsed_=0, previous_rate_=0, speed_=0;
    Eigen::Vector3f previous_force_=Eigen::Vector3f::Zero();
    Eigen::Vector3f previous_gyro_=Eigen::Vector3f::Zero();
    Eigen::Vector3f correction_=Eigen::Vector3f::Zero();
    std::array<std::array<SectionState,5>,2> lowpass_{};
    std::array<float,5> mean_{};
    float xx_=0, xy_=0, yy_=0, rr_=0;
};
} // namespace seastate::common
