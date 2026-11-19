#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>

// Research-only gate. Six nuisance columns per independent 20 s block:
// local starting field (uT) and a constant gyro bias (millirad/s).
// Retain only the five calibration coordinates shared across blocks.
struct SelectiveMag {
    using M = Eigen::MatrixXd;
    using V = Eigen::Vector3d;
    M normal = M::Zero(11,11), information = M::Zero(5,5);
    Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity(), bias = Eigen::Matrix3d::Zero();
    V start = V::Zero();
    int packets = 0;
    double hard = 0, scale = 0;

    static Eigen::Matrix3d skew(const V& v) {
        Eigen::Matrix3d k;
        k << 0,-v.z(),v.y(),v.z(),0,-v.x(),-v.y(),v.x(),0;
        return k;
    }
    static M profile(const M& n, int nuisance);
    void observe(const Eigen::Matrix3f& transport, const Eigen::Matrix3f& bias_transport,
                 const Eigen::Vector3f& predicted);
};
