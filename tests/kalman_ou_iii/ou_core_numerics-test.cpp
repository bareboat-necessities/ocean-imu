// Copyright 2026, Mikhail Grushinskiy
#define EIGEN_NON_ARDUINO

// Numerical contracts of the shared OU-II/OU-III core math:
//
//  * rot_and_B_from_wt() and integral_B_ds() must be accurate to a few ulp of
//    each entry's own magnitude for every dimensionless angle |omega|*dt,
//    in float (the deployed type) and double. The reference is the Van Loan
//    block exponential of [-W I 0; 0 0 I; 0 0 0] t in long double, which is
//    independent of the closed forms and their series.
//  * project_psd_ou_ii/iii() must never return an indefinite matrix, while
//    leaving already-PSD matrices (including singular ones and roundoff-level
//    negatives) unchanged apart from symmetrisation.
//
// Both groups fail against the previous implementation: the float integral
// returned IB(0,0) = 0 for omega = (0,0,0.01), dt = 0.005, and the six-state
// projection accepted diag(-1,1,1,1,1,1) because its LDLT "succeeded".

#include "kalman_ou_common/KalmanOUCoreMath.h"

#include <unsupported/Eigen/MatrixFunctions>

#include <array>
#include <cmath>
#include <iostream>
#include <limits>

namespace detail = ocean_imu::kalman::ou_detail;

namespace {

int failures = 0;

bool check(bool condition, const char* what, double got = 0, double limit = 0) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << what << " (value " << got << ", limit " << limit << ")\n";
    }
    return condition;
}

using LD = long double;
using Matrix3L = Eigen::Matrix<LD,3,3>;
using Vector3L = Eigen::Matrix<LD,3,1>;

struct Reference {
    Matrix3L R, B, IB;
    // Entry magnitudes of the three terms of each closed form; used as the
    // per-entry scale so a tiny entry obtained by cancellation is still
    // checked to relative precision.
    Matrix3L R_scale, B_scale, IB_scale;
};

Reference van_loan_reference(const Vector3L& w, LD t) {
    Eigen::Matrix<LD,9,9> A = Eigen::Matrix<LD,9,9>::Zero();
    A.block<3,3>(0,0) = -detail::skew<LD>(w);
    A.block<3,3>(0,3) = Matrix3L::Identity();
    A.block<3,3>(3,6) = Matrix3L::Identity();
    const Eigen::Matrix<LD,9,9> E = (A * t).exp();
    Reference ref;
    ref.R = E.block<3,3>(0,0);
    ref.B = E.block<3,3>(0,3);
    ref.IB = E.block<3,3>(0,6);

    const Matrix3L W = detail::skew<LD>(w).cwiseAbs();
    const Matrix3L W2 = (detail::skew<LD>(w) * detail::skew<LD>(w)).cwiseAbs();
    const LD t2 = t*t;
    // Leading-order magnitudes of the W and W^2 coefficients (1, 1/2, 1/6,
    // 1/24 times powers of t), which bound the true ones for all angles.
    ref.R_scale = Matrix3L::Identity() + t*W + t2*LD(0.5)*W2;
    ref.B_scale = Matrix3L::Identity()*t + t2*LD(0.5)*W + t2*t*W2/LD(6);
    ref.IB_scale = Matrix3L::Identity()*(t2/LD(2)) + t2*t*W/LD(6) + t2*t2*W2/LD(24);
    return ref;
}

template<typename T>
LD rel_error(const Eigen::Matrix<T,3,3>& got, const Matrix3L& ref, const Matrix3L& scale) {
    LD worst = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const LD err = std::abs(LD(got(i,j)) - ref(i,j));
            // Entries whose scale is zero must come out exactly zero.
            worst = std::max(worst, scale(i,j) > 0 ? err / scale(i,j) : (err > 0 ? LD(1) : LD(0)));
        }
    }
    return worst;
}

template<typename T>
bool check_integrals(const Eigen::Matrix<T,3,1>& w, T dt, const char* label) {
    using Matrix3 = Eigen::Matrix<T,3,3>;
    Matrix3 R, B, IB;
    detail::rot_and_B_from_wt(w, dt, R, B);
    detail::integral_B_ds(w, dt, IB);
    if (!R.allFinite() || !B.allFinite() || !IB.allFinite()) {
        return check(false, label);
    }
    const Reference ref = van_loan_reference(w.template cast<LD>(), LD(dt));
    // A handful of ulp of the working type; the reference error is ~1e-18.
    const LD limit = LD(32) * LD(std::numeric_limits<T>::epsilon());
    const LD eR = rel_error(R, ref.R, ref.R_scale);
    const LD eB = rel_error(B, ref.B, ref.B_scale);
    const LD eI = rel_error(IB, ref.IB, ref.IB_scale);
    bool ok = true;
    if (!check(eR <= limit, label, double(eR), double(limit))) { std::cerr << "  in R\n"; ok = false; }
    if (!check(eB <= limit, label, double(eB), double(limit))) { std::cerr << "  in B\n"; ok = false; }
    if (!check(eI <= limit, label, double(eI), double(limit))) { std::cerr << "  in IB\n"; ok = false; }
    return ok;
}

// Independent long-double evaluation of g0..g3 (series with 20 terms for
// |x| < 2, closed form above) for the branch-boundary continuity checks.
std::array<LD,4> reference_coeffs(LD x) {
    std::array<LD,4> g{};
    if (std::abs(x) < LD(2)) {
        for (int n = 1; n <= 4; ++n) {
            LD sum = 0, term = 1;
            for (int m = 2; m <= n; ++m) term /= LD(m);   // 1/n!
            LD y = 1;
            for (int k = 0; k < 20; ++k) {
                sum += ((k % 2) ? -term : term) * y;
                y *= x*x;
                term /= LD((2*k+n+1) * (2*k+n+2));
            }
            g[n-1] = sum;
        }
        return g;
    }
    const LD x2 = x*x;
    g[0] = std::sin(x)/x;
    g[1] = (LD(1)-std::cos(x))/x2;
    g[2] = (x-std::sin(x))/(x2*x);
    g[3] = (x2/LD(2)-LD(1)+std::cos(x))/(x2*x2);
    return g;
}

template<typename T>
void check_coefficients_at(T x, const char* label) {
    const auto c = detail::so3_integral_coeffs(x);
    const auto ref = reference_coeffs(LD(x));
    const T got[4] = {c.g0, c.g1, c.g2, c.g3};
    const LD limit = LD(16) * LD(std::numeric_limits<T>::epsilon());
    for (int i = 0; i < 4; ++i) {
        const LD rel = std::abs(LD(got[i]) - ref[i]) / std::abs(ref[i]);
        check(rel <= limit, label, double(rel), double(limit));
    }
}

template<typename T>
void test_so3_integrals(const char* type_name) {
    using Vector3 = Eigen::Matrix<T,3,1>;
    std::cerr << "SO(3) integrals, " << type_name << '\n';
    const Vector3 axis = Vector3(T(0.3), T(-0.5), T(0.81)).normalized();

    // Exactly zero angular rate: IB = I dt^2/2, B = I dt, R = I exactly.
    check_integrals<T>(Vector3::Zero(), T(0.005), "zero angular rate");
    {
        Eigen::Matrix<T,3,3> R, B, IB;
        detail::rot_and_B_from_wt<T>(Vector3::Zero(), T(0.005), R, B);
        detail::integral_B_ds<T>(Vector3::Zero(), T(0.005), IB);
        check(R == Eigen::Matrix<T,3,3>::Identity(), "zero rate R is not exactly I");
        check(IB(0,0) == T(0.5)*T(0.005)*T(0.005), "zero rate IB(0,0) is not dt^2/2");
    }

    // The reproduced float defect: x = 5e-5.
    {
        const Vector3 w(T(0), T(0), T(0.01));
        const T dt = T(0.005);
        Eigen::Matrix<T,3,3> IB;
        detail::integral_B_ds(w, dt, IB);
        const LD expected = LD(0.5)*LD(dt)*LD(dt) - LD(1e-4)*std::pow(LD(dt),4)/LD(24);
        const LD rel = std::abs(LD(IB(0,0)) - expected) / expected;
        check(rel <= LD(4)*LD(std::numeric_limits<T>::epsilon()),
              "omega=(0,0,0.01), dt=0.005: IB(0,0) != 1.25e-5", double(IB(0,0)), double(expected));
        check_integrals<T>(w, dt, "omega=(0,0,0.01), dt=0.005");
    }

    // Very small through large dimensionless angles, deployed dt and others.
    const std::array<T,3> dts = {T(0.005), T(0.02), T(0.1)};
    const std::array<T,14> rates = {T(1e-9), T(1e-6), T(1e-4), T(1e-3), T(0.01), T(0.05),
                                    T(0.2), T(0.5), T(1.0), T(3.0), T(10.0), T(31.4), T(60.0), T(150.0)};
    for (T dt : dts) {
        for (T rate : rates) check_integrals<T>((rate*axis).eval(), dt, "angle sweep");
    }

    // Immediately below, at and above the series/closed-form switch |x| = 1,
    // for both the full matrices and the scalar coefficients.
    const T one = T(1);
    const std::array<T,5> xs = {std::nextafter(one, T(0)) , one, std::nextafter(one, T(2)),
                                T(0.999), T(1.001)};
    const T dt = T(0.1);
    for (T x : xs) {
        check_integrals<T>((x/dt*axis).eval(), dt, "series/closed-form boundary");
        check_coefficients_at<T>(x, "coefficient accuracy at the boundary");
        check_coefficients_at<T>(-x, "coefficient accuracy at the negative boundary");
    }
    // Continuity: the jump across the switch is only rounding.
    {
        const auto lo = detail::so3_integral_coeffs(std::nextafter(one, T(0)));
        const auto hi = detail::so3_integral_coeffs(one);
        const T lov[4] = {lo.g0, lo.g1, lo.g2, lo.g3};
        const T hiv[4] = {hi.g0, hi.g1, hi.g2, hi.g3};
        for (int i = 0; i < 4; ++i) {
            const T jump = std::abs(lov[i]-hiv[i]) / std::abs(hiv[i]);
            check(jump <= T(32)*std::numeric_limits<T>::epsilon(), "coefficient jump at the switch",
                  double(jump), double(T(32)*std::numeric_limits<T>::epsilon()));
        }
    }
    // Small-angle coefficients themselves.
    for (T x : {T(0), T(1e-8), T(5e-5), T(1e-3), T(0.1), T(0.5), T(2.0), T(6.0)}) {
        if (x != T(0)) check_coefficients_at<T>(x, "coefficient accuracy");
    }
}

template<typename T, int N>
T min_eigenvalue(const Eigen::Matrix<T,N,N>& S) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<LD,N,N>> es(S.template cast<LD>());
    return T(es.eigenvalues().minCoeff());
}

template<typename T, int N>
bool returned_psd(const Eigen::Matrix<T,N,N>& S, const char* label) {
    if (!check(S.allFinite(), label)) return false;
    if (!check(S == S.transpose(), label)) return false;
    const T scale = S.cwiseAbs().maxCoeff();
    const T tol = T(4*N) * std::numeric_limits<T>::epsilon() * scale;
    const T lam = min_eigenvalue<T,N>(S);
    return check(lam >= -tol, label, double(lam), double(-tol));
}

template<typename T, int N, bool ou_iii>
void project(Eigen::Matrix<T,N,N>& S, T eps) {
    if constexpr (ou_iii) detail::project_psd_ou_iii<T,N>(S, eps);
    else detail::project_psd_ou_ii<T,N>(S, eps);
}

template<typename T, int N, bool ou_iii>
void test_psd_projection(const char* label) {
    using Matrix = Eigen::Matrix<T,N,N>;
    std::cerr << "PSD projection, " << label << '\n';
    const T eps = T(1e-12);

    // Valid SPD: unchanged bit for bit.
    Matrix A = Matrix::Random();
    Matrix spd = A * A.transpose() + T(0.1) * Matrix::Identity();
    detail::symmetrize<T,N>(spd);
    {
        Matrix S = spd;
        project<T,N,ou_iii>(S, eps);
        if constexpr (N > 4) check(S == spd, "SPD matrix changed");
        returned_psd<T,N>(S, "SPD result");
    }

    // PSD with zero eigenvalues (rank 2).
    {
        Eigen::Matrix<T,N,2> V = Eigen::Matrix<T,N,2>::Random();
        Matrix S = V * V.transpose();
        detail::symmetrize<T,N>(S);
        const Matrix before = S;
        project<T,N,ou_iii>(S, eps);
        if constexpr (N > 4) check(S == before, "singular PSD matrix changed");
        returned_psd<T,N>(S, "singular PSD result");
    }

    // Roundoff-level negative eigenvalue: harmless noise is left alone.
    {
        Eigen::SelfAdjointEigenSolver<Matrix> es(spd);
        Eigen::Matrix<T,N,1> lam = es.eigenvalues();
        lam(0) = -T(0.5) * std::numeric_limits<T>::epsilon() * lam.cwiseAbs().maxCoeff();
        Matrix S = es.eigenvectors() * lam.asDiagonal() * es.eigenvectors().transpose();
        detail::symmetrize<T,N>(S);
        project<T,N,ou_iii>(S, eps);
        returned_psd<T,N>(S, "roundoff-negative result");
    }

    // The counterexample diag(-1,1,...,1) and a hidden indefinite matrix with
    // a positive diagonal.
    {
        Matrix S = Matrix::Identity();
        S(0,0) = T(-1);
        project<T,N,ou_iii>(S, eps);
        returned_psd<T,N>(S, "diag(-1,1,...) result");
        check(S(0,0) >= T(0) && S(0,0) <= T(2)*eps + T(8)*std::numeric_limits<T>::epsilon(),
              "diag(-1,...) not floored at eps", double(S(0,0)), double(eps));
        for (int i = 1; i < N; ++i) {
            check(std::abs(S(i,i) - T(1)) <= T(8)*std::numeric_limits<T>::epsilon(),
                  "diag(-1,...) positive eigenvalue altered");
        }
    }
    {
        Matrix S = Matrix::Identity();
        S(0,1) = S(1,0) = T(2);  // eigenvalue -1 with positive diagonal
        project<T,N,ou_iii>(S, eps);
        returned_psd<T,N>(S, "positive-diagonal indefinite result");
    }

    // LDLT preserves inertia, not eigenvalue magnitudes. The correlated
    // negative block has one eigenvalue -(N-1)*amount*tol but a first pivot
    // of only -amount*tol; comparing D with -tol incorrectly accepts it.
    if constexpr (N > 4) {
        for (T scale : {T(1e-4), T(1), T(1e4)}) {
            const T tol = T(4*N)*std::numeric_limits<T>::epsilon()*scale;
            for (T amount : {T(0.25)/T(N-1), T(0.75)}) {
                Matrix S = Matrix::Zero();
                S(0,0) = scale;
                S.template bottomRightCorner<N-1,N-1>().setConstant(-amount*tol);
                const Matrix before = S;
                Eigen::LDLT<Matrix> ldlt(S);
                check(ldlt.info() == Eigen::Success && ldlt.vectorD().minCoeff() >= -tol,
                      "correlated negative block reproduces pivot-tolerance shortcut");
                project<T,N,ou_iii>(S, eps);
                returned_psd<T,N>(S, "correlated negative block result");
                if (amount < T(0.5)) check(S == before, "harmless correlated roundoff changed");
            }
        }
    }

    // Asymmetric input is symmetrised.
    {
        Matrix S = spd;
        S(0,N-1) += T(0.01);
        project<T,N,ou_iii>(S, eps);
        returned_psd<T,N>(S, "asymmetric result");
        check(std::abs(S(0,N-1) - (spd(0,N-1) + T(0.005))) <= T(8)*std::numeric_limits<T>::epsilon()*spd.cwiseAbs().maxCoeff(),
              "asymmetric input not symmetrised to the mean");
    }

    // Non-finite entries: existing contract replaces them (eps on the
    // diagonal, zero off it) and the result is finite and PSD.
    {
        Matrix S = spd;
        S(1,1) = std::numeric_limits<T>::quiet_NaN();
        S(0,2) = std::numeric_limits<T>::infinity();
        project<T,N,ou_iii>(S, eps);
        returned_psd<T,N>(S, "non-finite input result");
    }
}

// The same invariant applies to the shared OU-chain regularizer. Exercise
// the three- and four-state sizes used by the shipping OU process blocks.
template<typename T, int N>
void test_psd_regularizer() {
    using Matrix = Eigen::Matrix<T,N,N>;
    for (T scale : {T(1e-4), T(1), T(1e4)}) {
        const T tol = T(64)*std::numeric_limits<T>::epsilon()*std::max(T(1),scale);
        for (T amount : {T(0.25)/T(N-1), T(0.75)}) {
            Matrix S = Matrix::Zero();
            S(0,0) = scale;
            S.template bottomRightCorner<N-1,N-1>().setConstant(-amount*tol);
            const Matrix before = S;
            detail::regularize_psd_if_needed<T,N>(S);
            check(S.allFinite() && S == S.transpose(), "regularizer finite and symmetric");
            check(min_eigenvalue<T,N>(S) >= -tol, "regularizer correlated negative block");
            if (amount < T(0.5)) check(S == before, "regularizer harmless roundoff changed");
        }
    }
}

} // namespace

int main() {
    test_so3_integrals<float>("float");
    test_so3_integrals<double>("double");

    test_psd_projection<float,6,false>("OU-II float 6x6");
    test_psd_projection<float,6,true>("OU-III float 6x6");
    test_psd_projection<double,6,false>("OU-II double 6x6");
    test_psd_projection<double,6,true>("OU-III double 6x6");
    test_psd_projection<float,3,false>("OU-II float 3x3");
    test_psd_projection<double,3,true>("OU-III double 3x3");

    test_psd_regularizer<float,3>();
    test_psd_regularizer<double,3>();
    test_psd_regularizer<float,4>();
    test_psd_regularizer<double,4>();

    if (failures) {
        std::cerr << failures << " OU core numerics check(s) failed\n";
        return 1;
    }
    std::cout << "OU core numerics: all checks passed\n";
    return 0;
}
