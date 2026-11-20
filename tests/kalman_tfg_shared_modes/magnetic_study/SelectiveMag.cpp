#include "SelectiveMag.h"

SelectiveMag::M SelectiveMag::profile(const M& n, int nuisance) {
    const int p = static_cast<int>(n.rows()) - nuisance;
    const M cross = n.topRightCorner(nuisance,p);
    const M result = n.bottomRightCorner(p,p) -
        cross.transpose() * n.topLeftCorner(nuisance,nuisance).completeOrthogonalDecomposition().solve(cross);
    return (result + result.transpose()) * 0.5;
}

void SelectiveMag::observe(const Eigen::Matrix3f& transport, const Eigen::Matrix3f& bias_transport,
                           const Eigen::Vector3f& predicted) {
    if (packets == 0) {
        start = predicted.cast<double>();
        rotation.setIdentity(); bias.setZero(); normal.setZero();
    } else {
        const Eigen::Matrix3d step = transport.cast<double>();
        rotation = step * rotation;
        bias = step * bias + bias_transport.cast<double>();
    }
    // One row block per 0.4 s; ignore absolute heading and magnitude truth.
    // Two orthonormal diagonal, trace-free gain directions per 1% change.
    if (packets % 10 == 0) {
        const V field = rotation * start;
        M j = M::Zero(3,11);
        j.leftCols(3) = rotation;
        j.middleCols(3,3) = -0.001 * skew(field) * bias;
        j.middleCols(6,3).setIdentity();
        j.block<3,1>(0,9) = 0.01 / std::sqrt(2.) * V(field.x(),-field.y(),0.);
        j.block<3,1>(0,10) = 0.01 / std::sqrt(6.) * V(field.x(),field.y(),-2.*field.z());
        normal.noalias() += j.transpose() * j / (0.8 * 0.8);
    }
    if (++packets != 500) return;
    packets = 0;
    information += profile(normal,6);
    scale = profile(information,3).selfadjointView<Eigen::Lower>().eigenvalues().minCoeff();
    M reordered(5,5);
    reordered << information.bottomRightCorner(2,2), information.bottomLeftCorner(2,3),
                 information.topRightCorner(3,2), information.topLeftCorner(3,3);
    hard = profile(reordered,2).selfadjointView<Eigen::Lower>().eigenvalues().minCoeff();
}
