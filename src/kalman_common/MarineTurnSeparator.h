#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
#include <algorithm>
#include <array>
#include <cmath>
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
    void reset(const Config& config) {
        const Config saved = config;
        *this = MarineTurnSeparator();
        cfg_ = saved;
        if (!cfg_.forward_body.allFinite() || cfg_.forward_body.norm() < 0.5f)
            cfg_.enabled = false;
        else cfg_.forward_body.normalize();
    }
    Eigen::Vector3f update(float dt, const Eigen::Vector3f& gyro,
                          const Eigen::Vector3f& force, float gravity) {
        if (!cfg_.enabled) return force;
        if (!(dt > 0 && dt <= 0.1f) || !std::isfinite(dt) ||
            !gyro.allFinite() || !force.allFinite()) {
            reset(cfg_);
            return force;
        }
        proxy_.update(dt, gyro, force, gravity);
        const Eigen::Matrix3f R = proxy_.tiltQuaternion().toRotationMatrix();
        const Eigen::Vector3f bow = R * cfg_.forward_body;
        const float horizontal2 = bow.head<2>().squaredNorm();
        if (horizontal2 < 0.25f) {
            reset(cfg_);
            return force;
        }
        Eigen::Vector3f forward(bow.x(), bow.y(), 0);
        forward /= std::sqrt(horizontal2);
        const Eigen::Vector3f starboard(-forward.y(), forward.x(), 0);
        const Eigen::Vector3f omega = R * gyro;
        const float rate = omega.z() - bow.z() *
            (omega.x()*bow.x() + omega.y()*bow.y()) / horizontal2;
        elapsed_ += dt;
        if (!initialized_) {
            previous_force_ = force;
            previous_gyro_ = gyro;
            previous_rate_ = rate;
            initialized_ = true;
            return force;
        }
        const Eigen::Vector3f jerk = R * ((force-previous_force_)/dt +
            (0.5f*(gyro+previous_gyro_)).cross(0.5f*(force+previous_force_)));
        std::array<double,5> sample{{-double(rate)*rate,
            double(rate-previous_rate_)/dt, forward.dot(jerk), starboard.dot(jerk), rate}};
        previous_force_ = force;
        previous_gyro_ = gyro;
        previous_rate_ = rate;
        // Two cascaded Butterworth sections. Double precision is needed at
        // IMU rates because the learning bandwidth is only 0.035 Hz.
        const double k = std::tan(3.14159265358979323846 * 0.035 * dt);
        for (std::size_t section=0; section<2; ++section) {
            const double q = section == 0 ? 0.541196100146197 : 1.306562964876377;
            const double norm = 1/(1+k/q+k*k);
            const double b0=k*k*norm, a1=2*(k*k-1)*norm, a2=(1-k/q+k*k)*norm;
            for (std::size_t i=0; i<sample.size(); ++i) {
                auto& state=lowpass_[section][i];
                const double y=b0*sample[i]+state[0];
                state[0]=2*b0*sample[i]-a1*y+state[1];
                state[1]=b0*sample[i]-a2*y;
                sample[i]=y;
            }
        }
        // Ignore startup derivative/filter transients. Centered regression
        // absorbs the constant jerk from residual gyro/accelerometer bias.
        if (elapsed_ > 60.0f) {
            if (!statistics_initialized_) { mean_=sample; statistics_initialized_=true; }
            const double a=-std::expm1(-double(dt)/120.0);
            std::array<double,5> delta{};
            for (std::size_t i=0; i<sample.size(); ++i) {
                delta[i]=sample[i]-mean_[i]; mean_[i]+=a*delta[i];
            }
            xx_=(1-a)*(xx_+a*(delta[0]*delta[0]+delta[1]*delta[1]));
            xy_=(1-a)*(xy_+a*(delta[0]*delta[2]+delta[1]*delta[3]));
            yy_=(1-a)*(yy_+a*(delta[2]*delta[2]+delta[3]*delta[3]));
            rr_=(1-a)*(rr_+a*delta[4]*delta[4]);
            double target=0;
            // Constant heading gives no information about cruise speed. Do
            // not divide by a noise-only angular-acceleration regressor.
            if (elapsed_ > 120 && xx_ > 1e-9 && rr_ > 1e-7) {
                const double speed=xy_/xx_;
                if (std::abs(speed) >= 0.2 && std::abs(speed) <= 15 &&
                    xy_*xy_ > 0.015*xx_*yy_) target=speed;
            }
            speed_ += float(-std::expm1(-double(dt)/15.0)*(target-speed_));
        }
        correction_ = R.transpose()*(starboard*(speed_*rate));
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
    Config cfg_{};
    VerticalAccelComplementary proxy_{0.2f,0.02f,0};
    bool initialized_=false, statistics_initialized_=false;
    float elapsed_=0, previous_rate_=0, speed_=0;
    Eigen::Vector3f previous_force_=Eigen::Vector3f::Zero();
    Eigen::Vector3f previous_gyro_=Eigen::Vector3f::Zero();
    Eigen::Vector3f correction_=Eigen::Vector3f::Zero();
    std::array<std::array<std::array<double,2>,5>,2> lowpass_{};
    std::array<double,5> mean_{};
    double xx_=0, xy_=0, yy_=0, rr_=0;
};
} // namespace seastate::common
