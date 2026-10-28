#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

// Covariance hygiene shared by the wave MEKFs: symmetrisation, the PSD
// projection and the PSD regulariser of the process-noise blocks.
// KalmanOUCoreMath.h re-exports these names into ou_detail; its
// project_psd_ou_ii and project_psd_ou_iii are project_psd_regularize_large
// and project_psd_symmetrize_large.

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>

namespace ocean_imu::kalman::filter_detail {

// S = (S + S^T)/2. The .eval() matters: without it Eigen writes S(i,j) before
// reading S(j,i) for the mirrored entry, and the result is not symmetric.
template<typename T, int N>
inline void symmetrize(Eigen::Matrix<T,N,N>& S) {
    S = (T(0.5) * (S + S.transpose())).eval();
}

// Returns a symmetric matrix whose eigenvalues are non-negative up to the
// rounding tolerance psd_roundoff_tol(); eigenvalues that have to be repaired
// are lifted to eps. Non-finite entries are replaced (eps on the diagonal,
// zero elsewhere) before the check, as before. N <= 4 always goes through the
// eigensolver and lifts every non-positive eigenvalue. Larger matrices take
// the cheap path first: Eigen's pivoted LDLT reports Success for any finite
// symmetric matrix without an exact zero pivot, indefinite or not, so
// info() == Success proves nothing on its own; by Sylvester's law of inertia
// the signs of its D carry the PSD test. When the LDLT stops at an exact zero
// pivot (singular PSD matrices do this) its D is no longer a factorisation
// and the eigenvalues decide. Only an indefinite matrix is modified. N > 6
// without regularize_large (the TFG full state) is symmetrised only.
template<typename T, int N>
inline T psd_roundoff_tol(const Eigen::Matrix<T,N,N>& S) {
    // Backward error of a diagonally pivoted LDLT: O(N eps ||S||_max).
    return T(4 * N) * std::numeric_limits<T>::epsilon() * S.cwiseAbs().maxCoeff();
}

// Lifts every eigenvalue that is not positive to eps. With keep_if_psd, a
// matrix whose smallest eigenvalue is within -tol of zero is left unchanged.
template<typename T, int N>
inline void clamp_eigenvalues_to_floor(Eigen::Matrix<T,N,N>& S, T eps,
                                       bool keep_if_psd = false, T tol = T(0)) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,N,N>> es(S);
    if (es.info() != Eigen::Success) {
        S.diagonal().array() += eps;
        return;
    }
    Eigen::Matrix<T,N,1> lam = es.eigenvalues();
    if (keep_if_psd && lam.minCoeff() >= -tol) return;
    for (int i = 0; i < N; ++i) {
        if (!(lam(i) > T(0))) lam(i) = eps;
    }
    S = es.eigenvectors() * lam.asDiagonal() * es.eigenvectors().transpose();
}

template<typename T, int N, bool regularize_large>
inline void project_psd_impl(Eigen::Matrix<T,N,N>& S, T eps) {
    symmetrize<T,N>(S);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (!std::isfinite(S(i,j))) S(i,j) = (i == j) ? eps : T(0);
        }
    }

    if constexpr (N <= 4) {
        clamp_eigenvalues_to_floor<T,N>(S, eps);
    } else if constexpr (regularize_large || N <= 6) {
        const T tol = psd_roundoff_tol<T,N>(S);
        Eigen::LDLT<Eigen::Matrix<T,N,N>> ldlt(S);
        // Congruence preserves inertia, not eigenvalue magnitudes. A small
        // negative D pivot can correspond to an eigenvalue below -tol, so
        // apply the roundoff tolerance to S's eigenvalues, not to its pivots.
        if (ldlt.info() != Eigen::Success || !(ldlt.vectorD().minCoeff() >= T(0))) {
            clamp_eigenvalues_to_floor<T,N>(S, eps, true, tol);
        }
    }
    symmetrize<T,N>(S);
}

// Every size is repaired.
template<typename T, int N>
inline void project_psd_regularize_large(Eigen::Matrix<T,N,N>& S, T eps = T(1e-12)) {
    project_psd_impl<T,N,true>(S, eps);
}

// N > 6 is symmetrised only.
template<typename T, int N>
inline void project_psd_symmetrize_large(Eigen::Matrix<T,N,N>& S, T eps = T(1e-12)) {
    project_psd_impl<T,N,false>(S, eps);
}

template<typename T, int N>
inline void regularize_psd_if_needed(Eigen::Matrix<T,N,N>& S) {
    symmetrize<T,N>(S);
    T scale = std::max(T(1), S.cwiseAbs().maxCoeff());
    const T tol = T(64) * std::numeric_limits<T>::epsilon() * scale;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (!std::isfinite(S(i,j))) S(i,j) = (i == j) ? tol : T(0);
        }
    }

    Eigen::LDLT<Eigen::Matrix<T,N,N>> ldlt(S);
    if (ldlt.info() == Eigen::Success && ldlt.vectorD().minCoeff() >= T(0)) return;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,N,N>> es(S);
    if (es.info() != Eigen::Success) {
        S.diagonal().array() += tol;
        symmetrize<T,N>(S);
        return;
    }

    Eigen::Matrix<T,N,1> eigenvalues = es.eigenvalues();
    if (eigenvalues.minCoeff() >= -tol) return;
    for (int i = 0; i < N; ++i) eigenvalues(i) = std::max(T(0), eigenvalues(i));
    S = es.eigenvectors() * eigenvalues.asDiagonal() * es.eigenvectors().transpose();
    symmetrize<T,N>(S);
}

} // namespace ocean_imu::kalman::filter_detail
