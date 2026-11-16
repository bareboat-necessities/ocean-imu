#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif
#include <cmath>

namespace ocean_imu::kalman {
// Steady horizontal travel along the projected vessel bow. The wave state
// excludes this acceleration, but the physical accelerometer model includes
// it, together with its attitude and gyro-bias sensitivities.
template<class T> struct MarineTurnModel {
    using V = Eigen::Matrix<T,3,1>;
    using M = Eigen::Matrix<T,3,3>;
    struct Evaluation {
        V force = V::Zero();
        M down_jacobian = M::Zero();
        M rate_jacobian = M::Zero();
    };
    T speed = T(0);
    V bow = V::UnitX();
    V endpoint_gyro = V::Zero();
    void set(T speed_mps, const V& forward_body) {
        speed = T(0);
        if (!std::isfinite(speed_mps) || !forward_body.allFinite() ||
            forward_body.norm() < T(0.5)) return;
        bow = forward_body.normalized();
        speed = speed_mps;
    }
    static M skew(const V& v) {
        M result;
        result << T(0),-v.z(),v.y(),v.z(),T(0),-v.x(),-v.y(),v.x(),T(0);
        return result;
    }
    // down and bow are unit vectors in the same body frame; omega is the
    // endpoint angular rate after subtracting the estimated gyro bias.
    Evaluation evaluate(const V& down, const V& omega, const V& forward) const {
        Evaluation result;
        if (speed == T(0)) return result;
        const T c = down.dot(forward), h2 = T(1)-c*c;
        if (!(h2 > T(0.25))) return result;
        const T wb = omega.dot(forward);
        const T n = omega.dot(down)-c*wb;
        const T scale = speed/(h2*std::sqrt(h2));
        const V side = down.cross(forward);
        result.force = scale*n*side;
        result.rate_jacobian = scale*side*(down-c*forward).transpose();
        result.down_jacobian = scale*(side*(omega-wb*forward).transpose()
            + (T(3)*n*c/h2)*side*forward.transpose() - n*skew(forward));
        return result;
    }
};
} // namespace ocean_imu::kalman
