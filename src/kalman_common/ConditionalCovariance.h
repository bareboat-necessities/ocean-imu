#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// Shrink only the horizontal/vertical Schur cross block conditional on r.
// P_rr and P_rm stay fixed. The conditional covariance becomes the convex
// combination alpha C + (1-alpha) blockdiag(C_hh,C_zz), preserving PSD.
// Indices allow axis-major additive and column-major invariant coordinates.
namespace ocean_imu::kalman::filter_detail {
template<typename T, int NR, int NHk, int NBk, class Matrix>
bool shrink_conditional_cross(Matrix& P, const int (&ridx)[NR],
                              const int (&hidx)[NHk], const int (&zidx)[NBk], T alpha) {
        if (!(alpha >= T(0)) || !(alpha < T(1))) return false;
        Eigen::Matrix<T,NR,NR> Prr;
        Eigen::Matrix<T,NR,NBk> Prz;
        for (int i = 0; i < NR; ++i) {
            for (int j = 0; j < NR; ++j) Prr(i,j) = P(ridx[i], ridx[j]);
            for (int j = 0; j < NBk; ++j) Prz(i,j) = P(ridx[i], zidx[j]);
        }
        // Correlation scaling keeps the factorization well conditioned.
        Eigen::Matrix<T,NR,1> d;
        for (int i = 0; i < NR; ++i) {
            if (!(Prr(i,i) > T(0))) return false;
            d(i) = T(1) / std::sqrt(Prr(i,i));
        }
        const Eigen::LDLT<Eigen::Matrix<T,NR,NR>> ldlt(d.asDiagonal() * Prr * d.asDiagonal());
        if (ldlt.info() != Eigen::Success || !(ldlt.vectorD().minCoeff() > T(0))) return false;
        const Eigen::Matrix<T,NR,NBk> X = d.asDiagonal() * ldlt.solve(d.asDiagonal() * Prz);   // P_rr^-1 P_rz
        for (int i = 0; i < NHk; ++i)
            for (int j = 0; j < NBk; ++j) {
                T e = T(0);
                for (int k = 0; k < NR; ++k) e += P(hidx[i], ridx[k]) * X(k,j);
                const int r = hidx[i], c = zidx[j];
                const T v = e + alpha * (T(0.5) * (P(r,c) + P(c,r)) - e);
                P(r,c) = v;
                P(c,r) = v;
            }
        return true;
}
}  // namespace ocean_imu::kalman::filter_detail
