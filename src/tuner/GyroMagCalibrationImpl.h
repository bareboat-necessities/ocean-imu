#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Shared implementation: do not duplicate host and embedded equations.
#if !defined(OCEAN_IMU_GYRO_MAG_CALIBRATION_DECLARED)
#include "GyroMagCalibration.h"
#endif
#if defined(ARDUINO) && !defined(EIGEN_NON_ARDUINO)
#define OCEAN_IMU_MAG_INLINE EIGEN_DONT_INLINE
#else
#define OCEAN_IMU_MAG_INLINE inline EIGEN_DONT_INLINE
#endif

OCEAN_IMU_MAG_INLINE void GyroMagCalibration::predict(float dt, const Vector3& gyro) {
        if (!std::isfinite(dt) || dt <= 0.0f || dt > 0.1f || !gyro.allFinite()) {
            reset(); return;
        }
        const auto step = interval_.prepare(gyro, dt);
        if (!initialized_ || !step.predict) return;
        const Vector3 rate = step.rate - x_.segment<3>(6);
        const Eigen::Quaternionf dq = ocean_imu::magnetic::increment(rate, dt);
        const Matrix3 R = dq.toRotationMatrix();
        // Left-error derivative of exp(-omega dt) with respect to gyro bias.
        const Vector3 v = -rate * dt;
        const Matrix3 K = skew_(v);
        const float a2 = v.squaredNorm();
        float c1 = 0.5f - a2 / 24.0f, c2 = 1.0f / 6.0f - a2 / 120.0f;
        if (a2 > 1.0e-6f) {
            const float a = std::sqrt(a2);
            c1 = (1.0f - std::cos(a)) / a2;
            c2 = (a - std::sin(a)) / (a2 * a);
        }
        bias_rotation_ = R * bias_rotation_ + dt * (Matrix3::Identity() + c1*K + c2*K*K);
        rotation_ = dq * rotation_; rotation_.normalize();
        elapsed_ += dt;
        // Endpoint samples are correlated after interval averaging. Retain
        // the raw sample exposure variance instead of halving it.
        angle_variance_ += sq_(cfg_.gyro_white_std_radps) *
            std::max(dt, cfg_.gyro_sample_period_sec) * dt;
        if (elapsed_ > 1.0f) reset();
    }

OCEAN_IMU_MAG_INLINE bool GyroMagCalibration::update(const Vector3& mag) {
        if (!mag.allFinite() || mag.norm() < 1.0e-3f) return false;
        if (!initialized_) { initialize_(mag); return false; }
        if (!(elapsed_ > 0.0f)) return false; // duplicate mag without IMU time
        const float dt = elapsed_;
        const Matrix3 R = rotation_.toRotationMatrix();
        x_.head<3>() = R * x_.head<3>();
        const Matrix3 J = -skew_(x_.head<3>()) * bias_rotation_;
        const Eigen::Matrix<float,3,18> rows = R*P_.topRows<3>() + J*P_.middleRows<3>(6);
        const Matrix3 field_cov = rows.leftCols<3>()*R.transpose() + rows.middleCols<3>(6)*J.transpose();
        P_.topRows<3>() = rows;
        P_.leftCols<3>() = rows.transpose();
        P_.block<3,3>(0,0) = field_cov;
        const Matrix3 Sfield = skew_(x_.head<3>());
        P_.block<3,3>(0,0) += angle_variance_ * Sfield * Sfield.transpose();
        P_.block<3,3>(3,3).diagonal().array() += sq_(cfg_.mag_bias_rw_uT) * dt;
        P_.block<3,3>(6,6).diagonal().array() += sq_(cfg_.gyro_bias_rw_radps) * dt;
        clearRotation_();
        const Measurement H = measurement_();
        const Vector3 residual = mag - (matrix_() * x_.head<3>() + biasBodyUT());
        const float variance = sq_(cfg_.mag_white_std_uT);
        const Matrix3 S = H * P_ * H.transpose() + variance * Matrix3::Identity();
        const Eigen::LDLT<Matrix3> solve(S);
        if (solve.info() != Eigen::Success || !solve.isPositive()) { reset(); return false; }
        const float nis = residual.dot(solve.solve(residual));
        if (!std::isfinite(nis) || nis > 25.0f) return false;
        const Eigen::Matrix<float,18,3> gain = P_ * H.transpose() * solve.solve(Matrix3::Identity());
        const State next = x_ + gain * residual;
        const Matrix3 A = matrixFor_(next);
        if (!next.allFinite() || (A - Matrix3::Identity()).norm() > 0.25f ||
            A.determinant() < 0.5f || next.segment<3>(3).norm() > 0.35f*mag.norm() ||
            next.segment<3>(6).norm() > 0.1f) return false;
        // Expanded Joseph update; only the three measurement rows are dense.
        // Avoid several 18x18 stack temporaries on the embedded target.
        const Eigen::Matrix<float,3,18> HP = H * P_;
        Covariance Pnext = P_;
        Pnext.noalias() -= gain * HP;
        Pnext.noalias() -= HP.transpose() * gain.transpose();
        Pnext.noalias() += gain * S * gain.transpose();
        if (!Pnext.allFinite()) { reset(); return false; }
        x_ = next; P_ = (0.5f * (Pnext + Pnext.transpose())).eval();
        const float alpha = accepted_ ? -std::expm1(-dt/30.0f) : 1.0f;
        residual_sq_ += alpha * (residual.squaredNorm() - residual_sq_);
        ++accepted_;
        return true;
    }

OCEAN_IMU_MAG_INLINE void GyroMagCalibration::initialize_(const Vector3& mag) {
        x_.setZero(); x_.head<3>() = mag; P_.setZero();
        P_.block<3,3>(3,3).diagonal().setConstant(sq_(cfg_.mag_bias_std_uT));
        P_.block<3,3>(6,6).diagonal().setConstant(sq_(cfg_.gyro_bias_std_radps));
        for(int i=0;i<3;++i) P_(9+4*i,9+4*i) = sq_(cfg_.scale_std);
        for(int i=0;i<3;++i) for(int j=i+1;j<3;++j) {
            const int a=9+3*j+i,b=9+3*i+j;
            const float vs=sq_(cfg_.cross_axis_std),vr=sq_(cfg_.alignment_std_rad);
            P_(a,a)=P_(b,b)=vs+vr; P_(a,b)=P_(b,a)=vs-vr;
        }
        // The first raw reading seeds the field, so its calibration and
        // offset uncertainty is correlated with that field. Do not pretend
        // the initial field is an independently known magnetic reference.
        const auto Hc = measurement_().block<3,15>(0,3).eval();
        P_.block<3,15>(0,3) = -Hc * P_.block<15,15>(3,3);
        P_.block<15,3>(3,0) = P_.block<3,15>(0,3).transpose();
        P_.block<3,3>(0,0) = sq_(cfg_.mag_white_std_uT)*Matrix3::Identity() +
            Hc*P_.block<15,15>(3,3)*Hc.transpose();
        initialized_ = true; accepted_ = 0; residual_sq_ = 0.0f; clearRotation_();
    }

OCEAN_IMU_MAG_INLINE bool GyroMagCalibration::hasSignificantDistortion() const {
        if (!valid()) return false;
        const Vector3 field = x_.head<3>();
        const Matrix3 distortion = matrix_() - Matrix3::Identity();
        const Vector3 shift = distortion * field + biasBodyUT();
        // Jacobian of (A-I)f+b at the last posterior, including field, offset
        // and gain cross-covariances. No incoming measurement or truth enters.
        auto jacobian = measurement_();
        jacobian.leftCols<3>() = distortion;
        const Matrix3 covariance = jacobian * P_ * jacobian.transpose();
        const Eigen::LDLT<Matrix3> solve(covariance);
        const bool valid = solve.info() == Eigen::Success && (solve.vectorD().array() > 0.f).all();
        const float score = valid ? shift.dot(solve.solve(shift)) : 0.f;
        return std::isfinite(score) && score > 16.27f;
}

#undef OCEAN_IMU_MAG_INLINE
