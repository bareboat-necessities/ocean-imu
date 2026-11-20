#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Residual magnetic calibration from raw magnetometer and relative gyro motion.
// No gravity, attitude estimate, geographical field or heading is observed.
#define OCEAN_IMU_GYRO_MAG_CALIBRATION_DECLARED 1
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
    void predict(float dt, const Vector3& gyro);
    bool update(const Vector3& mag);
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
    void initialize_(const Vector3& mag);
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

// Host tests keep header-only integration; Arduino builds the same bodies
// once in GyroMagCalibration.cpp to bound Xtensa literal-pool expansion.
#if !defined(ARDUINO) || defined(EIGEN_NON_ARDUINO)
#include "GyroMagCalibrationImpl.h"
#endif
