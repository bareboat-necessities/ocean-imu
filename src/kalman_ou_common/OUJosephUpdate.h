#pragma once

// Full three-channel Joseph update for a possibly constrained gain. Keep the
// original nine-term K S K^T accumulation order. Each K(i,a)*S(a,b) factor
// is shared by all the output columns in row i.
namespace ocean_imu::kalman::ou_detail {

template<class Covariance, class Gain, class Innovation>
inline void joseph_update3(Covariance& P, const Gain& K, const Innovation& S,
                           const Gain& PCt) {
    using T = typename Covariance::Scalar;
    constexpr int n = Covariance::RowsAtCompileTime;
    for (int i = 0; i < n; ++i) {
        T factors[3][3];
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                factors[a][b] = K(i,a) * S(a,b);
        for (int j = i; j < n; ++j) {
            T kcp = 0, kcp_transposed = 0, ksk = 0;
            for (int a = 0; a < 3; ++a) {
                kcp += K(i,a) * PCt(j,a);
                kcp_transposed += K(j,a) * PCt(i,a);
                for (int b = 0; b < 3; ++b)
                    ksk += factors[a][b] * K(j,b);
            }
            const T updated = P(i,j) + (-(kcp + kcp_transposed) + ksk);
            P(i,j) = updated;
            P(j,i) = updated;
        }
    }
}

}  // namespace ocean_imu::kalman::ou_detail
