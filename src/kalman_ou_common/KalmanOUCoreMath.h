#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

// Integrated Ornstein-Uhlenbeck chain math of the OU filters (OU-II, OU-III,
// TFG): the OU primitives and their series, the discrete transition and the
// process noise of the integrated chain.

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

#include "kalman_common/CovarianceProjection.h"
#include "kalman_common/MekfAttitudeMath.h"
#include "kalman_common/PeriodicSchedule.h"

namespace ocean_imu::kalman::ou_detail {

// The attitude, covariance and schedule helpers the OU filters share with the
// CARMA filters, under the names the OU filters, their tests and the OU-III
// proof tools use.
using filter_detail::gyro_bias_radius_rad_s;
using filter_detail::project_gyro_bias;
using filter_detail::skew;
using filter_detail::quat_from_delta_theta;
using filter_detail::SO3IntegralCoeffs;
using filter_detail::so3_integral_coeffs;
using filter_detail::rot_and_B_from_wt;
using filter_detail::integral_B_ds;
using filter_detail::simpson_R_Q_RT;
using filter_detail::simpson_B_Q_BT;
using filter_detail::is_isotropic3;
using filter_detail::apply_left_error_reset;
using filter_detail::symmetrize;
using filter_detail::psd_roundoff_tol;
using filter_detail::clamp_eigenvalues_to_floor;
using filter_detail::project_psd_impl;
using filter_detail::regularize_psd_if_needed;
using filter_detail::retarget_period_elapsed_progress_preserving;
using filter_detail::periodic_update_due;

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

template<typename T, int N>
inline void project_psd_ou_ii(Eigen::Matrix<T,N,N>& S, T eps = T(1e-12)) {
    project_psd_impl<T,N,true>(S, eps);
}

template<typename T, int N>
inline void project_psd_ou_iii(Eigen::Matrix<T,N,N>& S, T eps = T(1e-12)) {
    project_psd_impl<T,N,false>(S, eps);
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
