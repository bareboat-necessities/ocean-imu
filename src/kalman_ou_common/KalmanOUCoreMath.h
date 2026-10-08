#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace ocean_imu::kalman::ou_detail {

// The estimated state is residual sensor bias, not vessel angular rate.
// This deliberately generous fixed ball is separate from the commissioned
// physical residual-bias bound (0.02 rad/s). See docs/ou-gyro-bias-projection.md.
inline constexpr double gyro_bias_radius_rad_s = 0.5;

template<typename T, class Derived>
inline void project_gyro_bias(Eigen::MatrixBase<Derived>& bias) {
    const T radius = T(gyro_bias_radius_rad_s);
    if (!bias.allFinite()) { bias.setZero(); return; }
    const T largest = bias.cwiseAbs().maxCoeff();
    if (largest == T(0)) return;

    // The ordinary path is bit-for-bit unchanged. Bound the rounding of the
    // three squared components before accepting a point close to the sphere.
    if (largest <= radius) {
        if (bias.squaredNorm() < radius*radius*(T(1)-T(8)*std::numeric_limits<T>::epsilon())) return;
        // Exact axial boundary points must remain unchanged.
        if (largest == radius && (bias.array() != T(0)).count() == 1) return;
        const double infinity = std::numeric_limits<double>::infinity();
        double upper = 0;
        for (int i = 0; i < 3; ++i) {
            const double x = static_cast<double>(bias(i));
            if (x != 0) {
                const double square = std::nextafter(x*x, infinity);
                upper = std::nextafter(upper+square, infinity);
            }
        }
        if (upper <= double(radius)*double(radius)) return;
    }

    // Scale first: even finite near-FLT_MAX/DBL_MAX vectors must retain their
    // direction rather than overflow their norm. The tiny inward rounding
    // allowance keeps the stored vector inside the Euclidean sphere.
    const Eigen::Matrix<T,3,1> unit = bias / largest;
    const T inward = radius*(T(1)-T(8)*std::numeric_limits<T>::epsilon());
    bias = unit * (inward / unit.norm());
    // Mean only. In particular, no covariance clipping, scaling or reset.
}

template<typename T>
inline T safe_inv_tau(T tau) {
    return T(1) / ((std::abs(tau) >= T(1e-8)) ? tau : std::copysign(T(1e-8), tau));
}

template<typename T>
struct OUPrims {
    T alpha;
    T em1;
};

template<typename T>
inline OUPrims<T> make_prims(T h, T tau) {
    const T x = h * safe_inv_tau(tau);
    return {std::exp(-x), std::expm1(-x)};
}

template<typename T>
inline Eigen::Matrix<T,3,3> skew(const Eigen::Matrix<T,3,1>& vec) {
    Eigen::Matrix<T,3,3> M;
    M << T(0), -vec(2), vec(1),
         vec(2), T(0), -vec(0),
        -vec(1), vec(0), T(0);
    return M;
}

template<typename T>
inline Eigen::Quaternion<T> quat_from_delta_theta(const Eigen::Matrix<T,3,1>& dtheta) {
    const T theta = dtheta.norm();
    const T half_theta = T(0.5) * theta;

    T w, k;
    if (theta < T(1e-2)) {
        const T t2 = theta * theta;
        const T t4 = t2 * t2;
        w = T(1);
        w = std::fma(-t2, T(1)/T(8), w);
        w = std::fma( t4, T(1)/T(384), w);
        k = T(0.5);
        k = std::fma(-t2, T(1)/T(48), k);
        k = std::fma( t4, T(1)/T(3840), k);
    } else {
        w = std::cos(half_theta);
        k = std::sin(half_theta) / theta;
    }

    const Eigen::Matrix<T,3,1> v = k * dtheta;
    Eigen::Quaternion<T> q(w, v.x(), v.y(), v.z());
    q.normalize();
    return q;
}

// Dimensionless integrated-OU coefficients.  With x = h/tau and a = exp(-x),
// every transition and process-noise entry of the integrated OU chain is
// tau^m g(x) with g(x) = x^p G(x) entire, p up to 7.  The closed forms of g
// subtract O(1) terms to leave O(x^p), so their rounding error grows like
// x^-p ulp and is large in float well inside the deployed range of x.  For
// |x| < 1, G is evaluated from its Taylor series instead: each table holds the
// exact rational coefficients rounded to double, long enough that truncation
// is below 2^-60 relative at |x| = 1.  At |x| >= 1 the closed forms lose at
// most a few hundred ulp, decreasing with x.
template<typename T, std::size_t N>
inline T ou_series(const double (&c)[N], T x) {
    T acc = T(c[N-1]);
    for (std::size_t k = N-1; k-- > 0;) acc = T(c[k]) + x*acc;
    return acc;
}

// vv: x^3 * sum_k c_k x^k
inline constexpr double ou_series_vv[23] = {
    0.3333333333333333, -0.25, 0.11666666666666667, -0.041666666666666664,
    0.012301587301587301, -0.003125, 0.0006999559082892416, -0.00014054232804232804,
    2.5603254769921436e-05, -4.2713844797178134e-06, 6.574572546794769e-07,
    -9.394540644540644e-08, 1.2527583625467223e-08, -1.5660435427300506e-09,
    1.8424603970627936e-10, -2.0472094573725967e-11, 2.1549737648205952e-12,
    -2.1549819854558418e-13, 2.052367710260443e-14, -1.865790606867772e-15,
    1.622427388301662e-16, -1.3520231459322325e-17, 1.0816186456847918e-18};
// vp: x^4 * sum_k c_k x^k
inline constexpr double ou_series_vp[23] = {
    0.125, -0.08333333333333333, 0.034722222222222224, -0.011111111111111112,
    0.002951388888888889, -0.000677910052910053, 0.0001380621693121693,
    -2.535273368606702e-05, 4.248420047031158e-06, -6.555301694190583e-07,
    9.379628675263596e-08, -1.2516877596242676e-08, 1.5653266211301926e-09,
    -1.8420105639020982e-10, 2.0469439308541306e-11, -2.154825793386156e-12,
    2.1549038894209988e-13, -2.052328564378316e-14, 1.8657719236058477e-15,
    -1.6224188783272867e-16, 1.3520194389358189e-17, -1.0816170984167235e-18,
    8.320137725273146e-20};
// va: x^2 * sum_k c_k x^k
inline constexpr double ou_series_va[24] = {
    0.5, -0.5, 0.2916666666666667, -0.125, 0.043055555555555555, -0.0125,
    0.0031498015873015874, -0.0007027116402116402, 0.0001408179012345679,
    -2.5628306878306877e-05, 4.2734721554166e-06, -6.576178451178452e-07,
    9.395687719100417e-08, -1.2528348341840405e-08, 1.5660913375033746e-09,
    -1.842488511635337e-10, 2.0472250765795654e-11, -2.1549819854558417e-12,
    2.154986095773465e-13, -2.0523696675545493e-14, 1.8657914965469113e-15,
    -1.622427775118679e-16, 1.3520233071059896e-17, -1.0816187101542946e-18};
// pp: x^5 * sum_k c_k x^k
inline constexpr double ou_series_pp[22] = {
    0.05, -0.027777777777777776, 0.00992063492063492, -0.002777777777777778,
    0.0006558641975308642, -0.00013558201058201059, 2.5102212602212603e-05,
    -4.225455614344503e-06, 6.536030841586397e-07, -9.364716705986547e-08,
    1.2506171567018128e-08, -1.5646096995303344e-09, 1.841560730741403e-10,
    -2.046678404335665e-11, 2.1546778219517165e-12, -2.1548257933861558e-13,
    2.0522894184961894e-14, -1.865753240343924e-15, 1.6224103683529112e-16,
    -1.3520157319394055e-17, 1.081615551148655e-18, -8.320131526282488e-20};
// pa: x^3 * sum_k c_k x^k
inline constexpr double ou_series_pa[24] = {
    0.16666666666666666, -0.16666666666666666, 0.09166666666666666,
    -0.03611111111111111, 0.01130952380952381, -0.002976190476190476,
    0.0006806657848324515, -0.00013833774250440918, 2.537778579445246e-05,
    -4.250507722729945e-06, 6.556907598574265e-07, -9.38077574982337e-08,
    1.2517642312615858e-08, -1.5653744159035164e-09, 1.8420386784746418e-10,
    -2.0469595500610993e-11, 2.1548340140214024e-12, -2.154907999738622e-13,
    2.0523305216724225e-14, -1.865772813284987e-15, 1.6224192651443036e-16,
    -1.352019600109576e-17, 1.0816171628862263e-18, -8.320137973232772e-20};
// aa: x^1 * sum_k c_k x^k
inline constexpr double ou_series_aa[25] = {
    1.0, -1.0, 0.6666666666666666, -0.3333333333333333, 0.13333333333333333,
    -0.044444444444444446, 0.012698412698412698, -0.0031746031746031746,
    0.0007054673721340388, -0.00014109347442680775, 2.565335898669232e-05,
    -4.275559831115387e-06, 6.577784355562133e-07, -9.39683479366019e-08,
    1.2529113058213587e-08, -1.5661391322766984e-09, 1.8425166262078804e-10,
    -2.0472406957865337e-11, 2.1549902060910882e-12, -2.1549902060910883e-13,
    2.0523716248486557e-14, -1.8657923862260506e-15, 1.622428161935696e-16,
    -1.3520234682797467e-17, 1.0816187746237974e-18};
// vS: x^5 * sum_k c_k x^k
inline constexpr double ou_series_vS[22] = {
    0.03333333333333333, -0.020833333333333332, 0.008134920634920636,
    -0.0024305555555555556, 0.0006007495590828924, -0.00012814153439153438,
    2.4225388808722142e-05, -4.133597883597884e-06, 6.44931200486756e-07,
    -9.290156859601304e-08, 1.2447288406283115e-08, -1.5603081699311859e-09,
    1.8386368151968837e-10, -2.044819718706403e-11, 2.1535680361934223e-12,
    -2.1542010251074124e-13, 2.0519566784981117e-14, -1.8655850909866063e-15,
    1.622329523596345e-16, -1.3519786619752702e-17, 1.0815993048339384e-18,
    -8.32006333738525e-20};
// pS: x^6 * sum_k c_k x^k
inline constexpr double ou_series_pS[21] = {
    0.013888888888888888, -0.006944444444444444, 0.0022569444444444442,
    -0.0005787037037037037, 0.00012566137566137565, -2.3974867724867726e-05,
    4.110633450911229e-06, -6.430041152263375e-07, 9.275244890324255e-08,
    -1.2436582377058568e-08, 1.5595912483313277e-09, -1.8381869820361883e-10,
    2.044554192187937e-11, -2.153420064758983e-12, 2.1541229290725694e-13,
    -2.051917532615985e-14, 1.865566407724682e-15, -1.6223210136219697e-16,
    1.3519749549788567e-17, -1.0815977575658702e-18, 8.320057138394592e-20};
// SS: x^7 * sum_k c_k x^k
inline constexpr double ou_series_SS[21] = {
    0.003968253968253968, -0.001736111111111111, 0.0005015432098765432,
    -0.00011574074074074075, 2.284752284752285e-05, -3.995811287477954e-06,
    6.324051462940352e-07, -9.185773074661964e-08, 1.2366993187099007e-08,
    -1.554572797132321e-09, 1.8348132333309738e-10, -2.0424299800402093e-11,
    2.1521623075662494e-12, -2.153420064758983e-13, 2.0515456467357805e-14,
    -1.8653795751054408e-15, 1.622231658891028e-16, -1.351934178018308e-17,
    1.0815799639830853e-18, -8.319982750506694e-20, 6.163005287699697e-21};
// Sa: x^4 * sum_k c_k x^k
inline constexpr double ou_series_Sa[23] = {
    0.041666666666666664, -0.041666666666666664, 0.022222222222222223,
    -0.008333333333333333, 0.002455357142857143, -0.000603505291005291,
    0.00012841710758377425, -2.4250440917107583e-05, 4.135685559296671e-06,
    -6.450917909251242e-07, 9.291303934161077e-08, -1.2448053122656297e-08,
    1.5603559647045097e-09, -1.8386649297694272e-10, 2.0448353379133717e-11,
    -2.153576256828669e-12, 2.1542051354250356e-13, -2.051958635792218e-14,
    1.8655859806657457e-15, -1.622329910413362e-16, 1.3519788231490273e-17,
    -1.0815993693034415e-18, 8.320063585344876e-20};
// phi_pa: x^2 * sum_k c_k x^k
inline constexpr double ou_series_phi_pa[19] = {
    0.5, -0.16666666666666666, 0.041666666666666664, -0.008333333333333333,
    0.001388888888888889, -0.0001984126984126984, 2.48015873015873e-05,
    -2.7557319223985893e-06, 2.755731922398589e-07, -2.505210838544172e-08,
    2.08767569878681e-09, -1.6059043836821613e-10, 1.1470745597729725e-11,
    -7.647163731819816e-13, 4.779477332387385e-14, -2.8114572543455206e-15,
    1.5619206968586225e-16, -8.22063524662433e-18, 4.110317623312165e-19};
// phi_Sa: x^3 * sum_k c_k x^k
inline constexpr double ou_series_phi_Sa[18] = {
    0.16666666666666666, -0.041666666666666664, 0.008333333333333333,
    -0.001388888888888889, 0.0001984126984126984, -2.48015873015873e-05,
    2.7557319223985893e-06, -2.755731922398589e-07, 2.505210838544172e-08,
    -2.08767569878681e-09, 1.6059043836821613e-10, -1.1470745597729725e-11,
    7.647163731819816e-13, -4.779477332387385e-14, 2.8114572543455206e-15,
    -1.5619206968586225e-16, 8.22063524662433e-18, -4.110317623312165e-19};
// phi_Ia (TFG integral of phi_Sa): x^4 * sum_k c_k x^k
inline constexpr double ou_series_phi_Ia[17] = {
    0.041666666666666664, -0.008333333333333333, 0.001388888888888889,
    -0.0001984126984126984, 2.48015873015873e-05, -2.7557319223985893e-06,
    2.755731922398589e-07, -2.505210838544172e-08, 2.08767569878681e-09,
    -1.6059043836821613e-10, 1.1470745597729725e-11, -7.647163731819816e-13,
    4.779477332387385e-14, -2.8114572543455206e-15, 1.5619206968586225e-16,
    -8.22063524662433e-18, 4.110317623312165e-19};

template<typename T>
struct OUDiscreteCoeffs {
    T phi_pa;
    T phi_Sa;
};

template<typename T>
inline OUDiscreteCoeffs<T> safe_phi_A_coeffs(T h, T tau) {
    OUDiscreteCoeffs<T> c{};
    const T x = h * safe_inv_tau(tau);
    const T tau2 = tau * tau;
    const T tau3 = tau2 * tau;

    if (std::abs(x) < T(1)) {
        const T x2 = x*x;
        c.phi_pa = tau2 * x2 * ou_series<T>(ou_series_phi_pa, x);
        c.phi_Sa = tau3 * (x2*x) * ou_series<T>(ou_series_phi_Sa, x);
    } else {
        const T em1 = std::expm1(-x);
        c.phi_pa = tau2 * (x + em1);
        c.phi_Sa = tau3 * (T(0.5)*x*x - x - em1);
    }
    return c;
}

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

template<typename T, int N>
inline void project_psd_ou_ii(Eigen::Matrix<T,N,N>& S, T eps = T(1e-12)) {
    project_psd_impl<T,N,true>(S, eps);
}

template<typename T, int N>
inline void project_psd_ou_iii(Eigen::Matrix<T,N,N>& S, T eps = T(1e-12)) {
    project_psd_impl<T,N,false>(S, eps);
}

// Scalar coefficients of the constant-rate SO(3) integrals, written as
// functions of the dimensionless angle x = |omega|*t so that, with W = [omega]x,
//
//   R(t)          = I - t g0(x) W + t^2 g1(x) W^2
//   B(t)          = int_0^t R = I t - t^2 g1(x) W + t^3 g2(x) W^2
//   int_0^t B(s)  = I t^2/2 - t^3 g2(x) W + t^4 g3(x) W^2
//
// g0 = sin(x)/x, g1 = (1 - cos x)/x^2, g2 = (x - sin x)/x^3 and
// g3 = (x^2/2 - 1 + cos x)/x^4. The closed forms of g1 (as written), g2 and
// g3 cancel: their rounding error grows like 1/x^2 ulp, so float keeps no
// correct digit below x ~ 3e-4 whatever |omega| is. For |x| < 1 all four are
// evaluated from their alternating even Taylor series through x^18; the
// first omitted term, at most 1/21! ~ 2e-20 at |x| = 1, is far below half a
// double ulp of every coefficient. At |x| >= 1 the closed-form
// amplification is bounded by ~12 ulp (g3) for float and double alike.
template<typename T>
struct SO3IntegralCoeffs {
    T g0, g1, g2, g3;
};

template<typename T>
inline SO3IntegralCoeffs<T> so3_integral_coeffs(T x) {
    SO3IntegralCoeffs<T> c{};
    if (std::abs(x) < T(1)) {
        const T y = x * x;
        // sum_{k=0..9} (-1)^k y^k / (2k+n)!, Horner from the highest order.
        auto alternating = [y](const T* inv_fact) {
            T acc = inv_fact[9];
            for (int k = 8; k >= 0; --k) acc = inv_fact[k] - y * acc;
            return acc;
        };
        static constexpr T n1[10] = {
            T(1.0/1.0), T(1.0/6.0), T(1.0/120.0), T(1.0/5040.0),
            T(1.0/362880.0), T(1.0/39916800.0), T(1.0/6227020800.0), T(1.0/1307674368000.0),
            T(1.0/355687428096000.0), T(1.0/1.21645100408832e+17)};
        static constexpr T n2[10] = {
            T(1.0/2.0), T(1.0/24.0), T(1.0/720.0), T(1.0/40320.0),
            T(1.0/3628800.0), T(1.0/479001600.0), T(1.0/87178291200.0), T(1.0/20922789888000.0),
            T(1.0/6402373705728000.0), T(1.0/2.43290200817664e+18)};
        static constexpr T n3[10] = {
            T(1.0/6.0), T(1.0/120.0), T(1.0/5040.0), T(1.0/362880.0),
            T(1.0/39916800.0), T(1.0/6227020800.0), T(1.0/1307674368000.0), T(1.0/355687428096000.0),
            T(1.0/1.21645100408832e+17), T(1.0/5.109094217170944e+19)};
        static constexpr T n4[10] = {
            T(1.0/24.0), T(1.0/720.0), T(1.0/40320.0), T(1.0/3628800.0),
            T(1.0/479001600.0), T(1.0/87178291200.0), T(1.0/20922789888000.0), T(1.0/6402373705728000.0),
            T(1.0/2.43290200817664e+18), T(1.0/1.1240007277776077e+21)};
        c.g0 = alternating(n1);
        c.g1 = alternating(n2);
        c.g2 = alternating(n3);
        c.g3 = alternating(n4);
        return c;
    }
    const T s = std::sin(x);
    const T half_s = std::sin(T(0.5) * x);
    const T x2 = x * x;
    c.g0 = s / x;
    c.g1 = T(2) * half_s * half_s / x2;
    c.g2 = (x - s) / (x2 * x);
    c.g3 = (T(0.5) * x2 - T(1) + std::cos(x)) / (x2 * x2);
    return c;
}

template<typename T>
inline void rot_and_B_from_wt(const Eigen::Matrix<T,3,1>& w, T t,
                              Eigen::Matrix<T,3,3>& R, Eigen::Matrix<T,3,3>& B) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    const Matrix3 W = skew(w);
    const Matrix3 W2 = W * W;
    const SO3IntegralCoeffs<T> c = so3_integral_coeffs(w.norm() * t);
    const T t2 = t * t;
    R = Matrix3::Identity() - (t*c.g0)*W + (t2*c.g1)*W2;
    B = Matrix3::Identity()*t - (t2*c.g1)*W + (t2*t*c.g2)*W2;
}

template<typename T>
inline void integral_B_ds(const Eigen::Matrix<T,3,1>& w, T step,
                          Eigen::Matrix<T,3,3>& IB) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    const Matrix3 W = skew(w);
    const SO3IntegralCoeffs<T> c = so3_integral_coeffs(w.norm() * step);
    const T T2 = step*step, T3 = T2*step;
    IB = Matrix3::Identity()*(T(0.5)*T2)
       - (T3*c.g2)*W
       + (T3*step*c.g3)*(W*W);
}

template<typename T>
inline Eigen::Matrix<T,3,3> simpson_R_Q_RT(const Eigen::Matrix<T,3,1>& w, T step,
                                            const Eigen::Matrix<T,3,3>& Q) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    Matrix3 R0, Btmp, Rm, R1;
    rot_and_B_from_wt(w, T(0), R0, Btmp);
    rot_and_B_from_wt(w, T(0.5)*step, Rm, Btmp);
    rot_and_B_from_wt(w, step, R1, Btmp);
    return (step/T(6)) * (R0*Q*R0.transpose() + T(4)*Rm*Q*Rm.transpose() + R1*Q*R1.transpose());
}

template<typename T>
inline Eigen::Matrix<T,3,3> simpson_B_Q_BT(const Eigen::Matrix<T,3,1>& w, T step,
                                            const Eigen::Matrix<T,3,3>& Q) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    Matrix3 Rtmp, B0, Bm, B1;
    rot_and_B_from_wt(w, T(0), Rtmp, B0);
    rot_and_B_from_wt(w, T(0.5)*step, Rtmp, Bm);
    rot_and_B_from_wt(w, step, Rtmp, B1);
    return (step/T(6)) * (B0*Q*B0.transpose() + T(4)*Bm*Q*Bm.transpose() + B1*Q*B1.transpose());
}

template<typename T>
inline bool is_isotropic3(const Eigen::Matrix<T,3,3>& S, T tol = T(1e-9)) {
    const T a = S(0,0), b = S(1,1), c = S(2,2);
    Eigen::Matrix<T,3,3> off_matrix = S;
    off_matrix.diagonal().setZero();
    const T off = off_matrix.cwiseAbs().sum();
    const T mean = (a+b+c)/T(3);
    return (std::abs(a-mean)+std::abs(b-mean)+std::abs(c-mean)+off)
        <= tol*(T(1)+std::abs(mean));
}

template<typename T, int NX>
inline void apply_left_error_reset(Eigen::Matrix<T,NX,NX>& covariance,
                                   const Eigen::Matrix<T,3,1>& dtheta) {
    if (!dtheta.allFinite() || !(dtheta.squaredNorm() > T(0))) return;
    const Eigen::Matrix<T,3,3> G = Eigen::Matrix<T,3,3>::Identity() + T(0.5)*skew(dtheta);

    Eigen::Matrix<T,3,3> Paa_old;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) Paa_old(i,j) = covariance(i,j);
    }

    Eigen::Matrix<T,3,3> GP;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            T sum = T(0);
            for (int k = 0; k < 3; ++k) sum += G(i,k)*Paa_old(k,j);
            GP(i,j) = sum;
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            T sum = T(0);
            for (int k = 0; k < 3; ++k) sum += GP(i,k)*G(j,k);
            covariance(i,j) = sum;
        }
    }

    for (int col = 3; col < NX; ++col) {
        const T old0 = covariance(0,col);
        const T old1 = covariance(1,col);
        const T old2 = covariance(2,col);
        const T new0 = G(0,0)*old0 + G(0,1)*old1 + G(0,2)*old2;
        const T new1 = G(1,0)*old0 + G(1,1)*old1 + G(1,2)*old2;
        const T new2 = G(2,0)*old0 + G(2,1)*old1 + G(2,2)*old2;
        covariance(0,col)=new0; covariance(1,col)=new1; covariance(2,col)=new2;
        covariance(col,0)=new0; covariance(col,1)=new1; covariance(col,2)=new2;
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = i+1; j < 3; ++j) {
            const T v = T(0.5)*(covariance(i,j)+covariance(j,i));
            covariance(i,j)=v; covariance(j,i)=v;
        }
    }
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

// Retarget a periodic scheduler without discarding elapsed service credit.
// If the new deadline has already passed, park immediately below it so the
// next valid sample services the owed update through periodic_update_due().
// Otherwise preserve elapsed bit-for-bit.
template<typename T>
inline T retarget_period_elapsed_progress_preserving(T elapsed, T period) {
    if (!(period > T(0)) || !std::isfinite(period)) return T(0);
    if (!(elapsed >= T(0)) || !std::isfinite(elapsed)) return T(0);
    if (elapsed < period) return elapsed;
    return std::nextafter(period, T(0));
}

template<typename T>
inline bool periodic_update_due(T dt, T period, T& elapsed) {
    if (!(dt > T(0)) || !std::isfinite(dt) || !(period > T(0)) || !std::isfinite(period)) return false;
    const T total = elapsed + dt;
    const T tol = T(16) * std::numeric_limits<T>::epsilon() * std::max(T(1), period);
    if (total + tol < period) {
        elapsed = total;
        return false;
    }
    elapsed = (total >= period) ? std::fmod(total, period) : T(0);
    if (!(elapsed >= T(0)) || !std::isfinite(elapsed) || elapsed >= period) elapsed = T(0);
    return true;
}

template<typename T, int Integrations>
struct IntegratedOUChain;

template<typename T>
struct IntegratedOUChain<T,2> {
    using MatrixAxis = Eigen::Matrix<T,3,3>;

    static void transition(T tau, T h, MatrixAxis& Phi) {
        const auto P = make_prims<T>(h, tau);
        const T phi_va = -tau * P.em1;
        const auto coeffs = safe_phi_A_coeffs<T>(h, tau);
        Phi.setZero();
        Phi(0,0)=T(1); Phi(0,2)=phi_va;
        Phi(1,0)=h;    Phi(1,1)=T(1); Phi(1,2)=coeffs.phi_pa;
        Phi(2,2)=std::min(P.alpha, T(1));
    }

    static void process_covariance(T tau, T h, T sigma2, MatrixAxis& Qd) {
        const T tau_eff = std::max(tau, T(1e-7));
        const T inv = T(1) / tau_eff;
        const T x = h * inv;
        if (std::abs(x) < T(1)) {
            // Qd = 2 sigma2 tau^(m-1) x^p G(x); see ou_series.
            const T s2=T(2)*sigma2;
            const T x2=x*x, x3=x2*x, x4=x3*x, x5=x4*x;
            const T t2=tau_eff*tau_eff, t3=t2*tau_eff, t4=t3*tau_eff;
            Qd(0,0)=s2*t2*x3*ou_series<T>(ou_series_vv,x);
            Qd(0,1)=s2*t3*x4*ou_series<T>(ou_series_vp,x);
            Qd(0,2)=s2*tau_eff*x2*ou_series<T>(ou_series_va,x);
            Qd(1,1)=s2*t4*x5*ou_series<T>(ou_series_pp,x);
            Qd(1,2)=s2*t2*x3*ou_series<T>(ou_series_pa,x);
            Qd(2,2)=s2*x*ou_series<T>(ou_series_aa,x);
            Qd(1,0)=Qd(0,1); Qd(2,0)=Qd(0,2); Qd(2,1)=Qd(1,2);
            regularize_psd_if_needed<T,3>(Qd);
            return;
        }

        const T a=std::exp(-x), a2=a*a, qc=T(2)*sigma2*inv;
        const T t2=tau_eff*tau_eff, t3=t2*tau_eff, t4=t3*tau_eff, t5=t4*tau_eff;
        const T x2=x*x, x3=x2*x;
        const T K00=t3*(-a2+T(4)*a+T(2)*x-T(3))/T(2);
        const T K01=t4*(a2+T(2)*a*(x-T(1))+x2-T(2)*x+T(1))/T(2);
        const T K02=t2*(a2-T(2)*a+T(1))/T(2);
        const T K11=t5*(-a2/T(2)-T(2)*a*x+x3/T(3)-x2+x+T(1)/T(2));
        const T K12=t3*(-a2-T(2)*a*x+T(1))/T(2);
        const T K22=tau_eff*(T(1)-a2)/T(2);
        Qd << qc*K00,qc*K01,qc*K02,
              qc*K01,qc*K11,qc*K12,
              qc*K02,qc*K12,qc*K22;
        regularize_psd_if_needed<T,3>(Qd);
    }
};

template<typename T>
struct IntegratedOUChain<T,3> {
    using MatrixAxis = Eigen::Matrix<T,4,4>;

    static void transition(T tau, T h, MatrixAxis& Phi) {
        const auto P = make_prims<T>(h, tau);
        const T phi_va = -tau*P.em1;
        const auto coeffs = safe_phi_A_coeffs<T>(h, tau);
        Phi.setZero();
        Phi(0,0)=T(1); Phi(0,3)=phi_va;
        Phi(1,0)=h; Phi(1,1)=T(1); Phi(1,3)=coeffs.phi_pa;
        Phi(2,0)=T(0.5)*h*h; Phi(2,1)=h; Phi(2,2)=T(1); Phi(2,3)=coeffs.phi_Sa;
        Phi(3,3)=std::min(P.alpha,T(1));
    }

    static void process_covariance(T tau, T h, T sigma2, MatrixAxis& Qd) {
        const T tau_eff = std::max(tau, T(1e-7));
        const T inv = T(1) / tau_eff;
        const T x = h * inv;

        Eigen::Matrix<T,3,3> marginal;
        IntegratedOUChain<T,2>::process_covariance(tau_eff, h, sigma2, marginal);
        Qd.setZero();
        const int idx[3] = {0,1,3};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) Qd(idx[i],idx[j]) = marginal(i,j);
        }

        T qvS, qpS, qSS, qSa;
        if (std::abs(x) < T(1)) {
            const T s2=T(2)*sigma2;
            const T x2=x*x, x4=x2*x2;
            const T t3=tau_eff*tau_eff*tau_eff, t4=t3*tau_eff, t5=t4*tau_eff, t6=t5*tau_eff;
            qvS=s2*t4*(x4*x)*ou_series<T>(ou_series_vS,x);
            qpS=s2*t5*(x4*x2)*ou_series<T>(ou_series_pS,x);
            qSS=s2*t6*(x4*x2*x)*ou_series<T>(ou_series_SS,x);
            qSa=s2*t3*x4*ou_series<T>(ou_series_Sa,x);
        } else {
            const T a=std::exp(-x), a2=a*a, qc=T(2)*sigma2*inv;
            const T t4=tau_eff*tau_eff*tau_eff*tau_eff;
            const T t5=t4*tau_eff, t6=t5*tau_eff, t7=t6*tau_eff;
            const T x2=x*x, x3=x2*x, x4=x3*x, x5=x4*x;
            const T K02=t5*(-T(3)*a2+T(3)*a*(x2+T(4))+x3-T(3)*x2+T(6)*x-T(9))/T(6);
            const T K12=t6*(a2/T(2)+a*(-x2+T(2)*x-T(2))/T(2)+x4/T(8)-x3/T(2)+x2-x+T(1)/T(2));
            const T K22=t7*(-a2/T(2)+a*x2+T(2)*a+x5/T(20)-x4/T(4)+T(2)*x3/T(3)-x2+x-T(3)/T(2));
            const T K23=t4*(a2-a*(x2+T(2))+T(1))/T(2);
            qvS=qc*K02;
            qpS=qc*K12;
            qSS=qc*K22;
            qSa=qc*K23;
        }

        Qd(0,2)=qvS; Qd(2,0)=qvS;
        Qd(1,2)=qpS; Qd(2,1)=qpS;
        Qd(2,2)=qSS;
        Qd(2,3)=qSa; Qd(3,2)=qSa;
        regularize_psd_if_needed<T,4>(Qd);
    }
};

} // namespace ocean_imu::kalman::ou_detail
