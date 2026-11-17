#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Residual magnetic calibration from raw magnetometer and relative gyro motion.
// No gravity, attitude estimate, geographical field or heading is observed.
#include "util/MagneticRotation.h"
#include "kalman_common/GyroInterval.h"
#include <cmath>

class GyroMagCalibration {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using Vector3 = Eigen::Vector3f;
    using Matrix3 = Eigen::Matrix3f;
    using State = Eigen::Matrix<float, 18, 1>;
    using Covariance = Eigen::Matrix<float, 18, 18>;
    using Measurement = Eigen::Matrix<float, 3, 18>;
    struct Config {
        float mag_white_std_uT = 0.8f;
        float mag_bias_std_uT = 1.15470054f;
        float mag_bias_rw_uT = 0.01f;
        float gyro_bias_std_radps = 0.00050383316f;
        float gyro_bias_rw_radps = 0.00001f;
        float gyro_white_std_radps = 0.00157f;
        float gyro_sample_period_sec = 0.005f;
        float scale_std = 0.008660254f;
        float cross_axis_std = 0.005773503f;
        float alignment_std_rad = 0.010076663f;
    };
    void reset() {
        x_.setZero(); P_.setZero(); initialized_ = false;
        accepted_ = 0; residual_sq_ = 0.0f;
        interval_.reset(); clearRotation_();
    }
    void setGyroMode(ocean_imu::kalman::filter_detail::GyroSampleMode mode) {
        interval_.set_mode(mode); interval_.reset();
    }
    void predict(float dt, const Vector3& gyro) {
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
    bool update(const Vector3& mag) {
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
    Vector3 biasBodyUT() const { return x_.segment<3>(3); }
    Matrix3 correction() const { return matrix_().inverse(); }
    bool valid() const { return initialized_ && accepted_ >= 500; }
    float residualRmsUT() const { return std::sqrt(std::max(0.0f,residual_sq_)); }
    const Covariance& covariance() const { return P_; }
    const State& state() const { return x_; }
private:
    static float sq_(float x) { return x*x; }
    static Matrix3 skew_(const Vector3& v) {
        Matrix3 s; s << 0,-v.z(),v.y(),v.z(),0,-v.x(),-v.y(),v.x(),0; return s;
    }
    static Matrix3 matrixFor_(const State& x) {
        Matrix3 A = Matrix3::Identity();
        for (int j=0;j<3;++j) A.col(j) += x.segment<3>(9+3*j);
        return A;
    }
    Matrix3 matrix_() const { return matrixFor_(x_); }
    Measurement measurement_() const {
        Measurement H = Measurement::Zero();
        H.block<3,3>(0,0) = matrix_(); H.block<3,3>(0,3).setIdentity();
        for (int j=0;j<3;++j) H.block<3,3>(0,9+3*j) = x_[j] * Matrix3::Identity();
        return H;
    }
    void initialize_(const Vector3& mag) {
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
    void clearRotation_() { rotation_.setIdentity(); bias_rotation_.setZero(); elapsed_=angle_variance_=0.0f; }
    Config cfg_{};
    State x_ = State::Zero();
    Covariance P_ = Covariance::Zero();
    Eigen::Quaternionf rotation_ = Eigen::Quaternionf::Identity();
    Matrix3 bias_rotation_ = Matrix3::Zero();
    float elapsed_ = 0.0f, angle_variance_ = 0.0f, residual_sq_ = 0.0f;
    unsigned accepted_ = 0;
    bool initialized_ = false;
    ocean_imu::kalman::filter_detail::GyroInterval<Vector3> interval_;
};
