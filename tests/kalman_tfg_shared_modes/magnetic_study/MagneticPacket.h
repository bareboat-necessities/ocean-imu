#pragma once
#include <Eigen/Dense>

namespace ocean_imu::research {
using SensorCovariance = Eigen::Matrix<float,18,18>;

// Forecast the noise of a corrected packet using the prior to that packet.
// Keep the offset/affine cross covariance and allow for unapplied slew.
inline Eigen::Matrix3f correctedMagneticCovariance(const Eigen::Matrix3f& applied,
        const Eigen::Vector3f& corrected_field, const Eigen::Matrix3f& noise,
        const Eigen::Vector3f& lag, const SensorCovariance& posterior) {
    Eigen::Matrix<float,3,18> j = Eigen::Matrix<float,3,18>::Zero();
    j.middleCols<3>(3) = -applied;
    for (int k = 0; k < 3; ++k) j.middleCols<3>(9+3*k) = -corrected_field[k]*applied;
    return applied*noise*applied.transpose() + j*posterior*j.transpose() + lag*lag.transpose();
}
inline bool quarantinePacket(const Eigen::Vector3f& received, const Eigen::Vector3f& predicted) {
    return (received-predicted).norm() > 6.f;
}
} // namespace ocean_imu::research
