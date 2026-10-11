#pragma once

// Per-observation covariance validation/installation. The caller retains its
// original measurement equation, coordinate transform and correction method.
namespace ocean_imu {
template <class Matrix>
bool validMeasurementCovariance(const Matrix& covariance) {
    using T = typename Matrix::Scalar;
    return covariance.allFinite() &&
        covariance.isApprox(covariance.transpose(), T(1e-5)) &&
        Eigen::LLT<Matrix>(covariance).info() == Eigen::Success;
}

template <class Matrix, class Correction>
bool withMeasurementCovariance(Matrix& configured, const Matrix& observation, Correction correct) {
    if (!validMeasurementCovariance(observation)) return false;
    const Matrix saved = configured;
    configured = observation;
    correct();
    configured = saved;
    return true;
}
} // namespace ocean_imu
