#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif
#include <algorithm>
#include <array>
#include <cmath>
#include "kalman_common/GyroInterval.h"

namespace seastate::common {

// Startup-only variable projection. Transport both observations into the
// initial body frame. For a small constant body gyro bias b the transported
// vector y changes by [y]x A b, A = integral R(t) dt. Center magnetic
// observations to eliminate the unknown earth field. Polynomial moments
// integrate the tapered force and slope for any actual handover time.
class WaveStartup {
    using V = Eigen::Vector3f;
    using M = Eigen::Matrix3f;
    using Q = Eigen::Quaternionf;
    using ForceModel = Eigen::Matrix<double,3,4>;
    using Mode = ocean_imu::kalman::filter_detail::GyroSampleMode;

    static M derivative(const V& y, const M& a) {
        M j;
        for (int k = 0; k < 3; ++k) j.col(k) = y.cross(a.col(k));
        return j;
    }
    static Q rotation(const V& v) {
        const float angle = v.norm();
        if (angle > 1e-6f) return Q(Eigen::AngleAxisf(angle, v / angle));
        return Q(1.0f, v.x() * 0.5f, v.y() * 0.5f, v.z() * 0.5f).normalized();
    }

public:
    struct Seed {
        Q attitude = Q::Identity();
        V field = V::Zero();
        V gyro_bias = V::Zero();
    };
    WaveStartup() { reset(); }
    void reset() {
        time_ = 0; last_mag_ = -1; nmag_ = 0; have_imu_ = false;
        relative_.setIdentity(); integrated_rotation_.setZero();
        mean_field_.setZero(); mean_j_.setZero(); info_.setZero(); gradient_.setZero();
        scatter_ = 0; force_energy_ = 0;
        field_block_sum_.setZero(); field_previous_mean_.setZero();
        field_block_n_ = field_previous_n_ = 0;
        field_block_start_ = field_time_sum_ = field_previous_time_ = 0;
        field_block_norm_ = field_previous_norm_ = 0;
        field_discontinuity_ = false;
        field_diff_n_ = 0; field_diff_energy_ = 0;
        field_last_.setZero();
        for (auto& m : moment_) m.setZero();
    }

    void imu(float dt, const V& gyro, const V& accel, Mode mode) {
        if (!(dt > 0 && dt <= 0.1f) || !gyro.allFinite() || !accel.allFinite()) {
            reset();
            return;
        }
        // Start a fresh contiguous window if an application defers startup.
        if (time_ + dt > 35.001) reset();
        if (!have_imu_) {
            last_gyro_ = gyro; have_imu_ = true;
            previous_force_.col(0) = accel.cast<double>();
            previous_force_.rightCols<3>().setZero();
            return;
        }
        const V rate = mode == Mode::Instantaneous ? V(0.5f*(gyro + last_gyro_)) : gyro;
        last_gyro_ = gyro;
        const M old_r = relative_.toRotationMatrix();
        relative_ = (relative_ * rotation(dt * rate)).normalized();
        const M r = relative_.toRotationMatrix();
        integrated_rotation_ += (0.5f*dt)*(old_r+r);
        const V force = r*accel;
        ForceModel next;
        next.col(0) = force.cast<double>();
        next.rightCols<3>() = derivative(force, integrated_rotation_).cast<double>();

        // Midpoint quadrature of the linearized force. Normalize the time
        // coordinate for stable evaluation even at the 35 s boundary.
        const double u = (time_ + 0.5*dt)/35.0;
        const ForceModel increment = (0.5*dt)*(previous_force_+next);
        double power = 1;
        for (auto& m : moment_) { m += power*increment; power *= u; }
        force_energy_ += (0.5*dt)*(previous_force_.col(0).squaredNorm()+next.col(0).squaredNorm());
        previous_force_ = next;
        time_ += dt;
    }

    void mag(const V& body) {
        if (!have_imu_ || !body.allFinite() || body.squaredNorm() < 1e-8f || last_mag_ == time_) return;
        const V field = relative_ * body;
        checkFieldContinuity_(field);
        const M j = derivative(field, integrated_rotation_);
        const V delta = field-mean_field_;
        const M dj = j-mean_j_;
        ++nmag_;
        const float alpha = 1.0f/static_cast<float>(nmag_);
        mean_field_ += alpha*delta;
        mean_j_ += alpha*dj;
        info_.noalias() += (1.0f-alpha)*dj.transpose()*dj;
        gradient_.noalias() += (1.0f-alpha)*dj.transpose()*delta;
        scatter_ += (1.0f-alpha)*delta.squaredNorm();
        last_mag_ = time_;
    }

    bool solve(float gravity, bool scale_slope, Seed& out) const {
        if (field_discontinuity_ || time_ < 19.99 || nmag_ < 10 || time_-last_mag_ > 0.2 || !(gravity > 0)) return false;
        // For w(x)=30 x²(1-x)², E[x]=1/2 and Var[x]=1/28.
        const double k = 35/time_;
        const double c = 30*k*k/time_;
        const ForceModel f = c*(moment_[2]-2*k*moment_[3]+k*k*moment_[4]);
        const ForceModel xf = c*k*(moment_[3]-2*k*moment_[4]+k*k*moment_[5]);
        const ForceModel slope = (28/time_)*(xf-0.5*f);
        const V avg = (moment_[0].col(0)/time_).cast<float>();
        const double variance = std::max(0.0, force_energy_/time_ - (moment_[0].col(0)/time_).squaredNorm());
        const double scale = scale_slope ? std::clamp(std::sqrt(variance), 0.15, 1.0) : 1.0;
        const float sigma = static_cast<float>(0.005*scale*k*k*k);
        const M js = slope.rightCols<3>().cast<float>() / sigma;
        const V s = slope.col(0).cast<float>() / sigma;
        constexpr float mag_var = 0.8f*0.8f;
        M normal = info_/mag_var + js.transpose()*js;
        normal.diagonal().array() += 1e6f; // prior sigma 0.001 rad/s
        const V b = normal.ldlt().solve(-gradient_/mag_var-js.transpose()*s);
        const V drift = integrated_rotation_*b;
        if (!b.allFinite() || b.norm() > 0.004f || drift.norm() > 0.15f) return false;
        const float residual = scatter_+2*gradient_.dot(b)+b.dot(info_*b);
        if (residual > 12.0f*static_cast<float>(nmag_)) return false;
        const V corrected = f.col(0).cast<float>() + f.rightCols<3>().cast<float>()*b;
        const V uniform = avg + (moment_[0].rightCols<3>()/time_).cast<float>()*b;
        if (std::abs(corrected.norm()-gravity) > 0.5f || uniform.norm() < gravity*0.5f) return false;
        const V down = -corrected.normalized();
        if (down.dot(-uniform.normalized()) < 0.99938756f) return false; // cos(0.035)
        const V field = mean_field_ + mean_j_*b;
        const V horizontal = field-down*down.dot(field);
        if (!(horizontal.norm() > 5.0f)) return false;
        M frame;
        frame.row(0) = horizontal.normalized().transpose();
        frame.row(1) = down.cross(horizontal.normalized()).transpose();
        frame.row(2) = down.transpose();
        const Seed candidate{(Q(frame)*rotation(-drift)*relative_).normalized(), frame*field, b};
        if (!candidate.attitude.coeffs().allFinite() || !candidate.field.allFinite()) return false;
        out = candidate;
        return true;
    }

private:
    // Compare adjacent half-second means in the gyro-transported frame.
    // The existing 0.8 uT observation model sets the uncertainty of their
    // difference. Include rotation over the mean epochs at the accepted
    // four-sigma gyro-bias limit. A fit must not absorb a local field jump
    // as a slowly accumulating gyro bias. This remains a consistency gate,
    // not a certificate against every smoothly varying magnetic disturbance.
    void checkFieldContinuity_(const V& field) {
        if (field_diff_n_ > 0) field_diff_energy_ += (field-field_last_).squaredNorm();
        ++field_diff_n_; field_last_ = field;
        if (field_block_n_ == 0) field_block_start_ = time_;
        ++field_block_n_;
        field_block_sum_ += field;
        field_time_sum_ += time_;
        field_block_norm_ = std::max(field_block_norm_, field.norm());
        if (time_ - field_block_start_ < 0.5 || field_block_n_ < 2) return;
        const float count = static_cast<float>(field_block_n_);
        const V mean = field_block_sum_ / count;
        const double epoch = field_time_sum_ / count;
        if (field_previous_n_ > 0) {
            // Adjacent gyro-transported differences upper-estimate white
            // magnetic noise: genuine smooth motion/model residual adds energy.
            // Keep the existing prior until 25 samples, then allow a conservative
            // twofold sample standard deviation, capped by that prior.
            const float observed = field_diff_n_ >= 25
                ? std::min(0.8f, std::max(1e-4f, 2.0f*std::sqrt(field_diff_energy_ /
                    (6.0f*static_cast<float>(field_diff_n_-1))))) : 0.8f;
            const float sigma = observed * std::sqrt(1.0f/count + 1.0f/static_cast<float>(field_previous_n_));
            const float rotation = 0.004f * static_cast<float>(epoch-field_previous_time_) *
                std::max(field_previous_norm_, field_block_norm_);
            if ((mean-field_previous_mean_).norm() > 6.0f*sigma + rotation)
                field_discontinuity_ = true;
        }
        field_previous_mean_ = mean;
        field_previous_time_ = epoch;
        field_previous_n_ = field_block_n_;
        field_previous_norm_ = field_block_norm_;
        field_block_sum_.setZero(); field_time_sum_ = 0;
        field_block_n_ = 0; field_block_norm_ = 0;
    }
    V field_last_;
    unsigned field_diff_n_{};
    float field_diff_energy_{};
    V field_block_sum_, field_previous_mean_;
    double field_block_start_{}, field_time_sum_{}, field_previous_time_{};
    unsigned field_block_n_{}, field_previous_n_{};
    float field_block_norm_{}, field_previous_norm_{};
    bool field_discontinuity_{};
    double time_{}, last_mag_{}, force_energy_{};
    unsigned nmag_{};
    bool have_imu_{};
    Q relative_;
    M integrated_rotation_, mean_j_, info_;
    V last_gyro_, mean_field_, gradient_;
    float scatter_{};
    ForceModel previous_force_;
    std::array<ForceModel,6> moment_;
};
} // namespace seastate::common
