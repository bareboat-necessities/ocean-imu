/*
  Copyright (c) 2026 Mikhail Grushinskiy

  Independent checks of the experimental CARMA vertical prior:

    * the realization reproduces g^2 |B|^2/|P|^2 from the polynomial factors;
    * relative degree three (C G = C A G = 0, C A^2 G = g) and the h, v, a rows;
    * every pole strictly in the open left half plane;
    * F and Q against a double-precision Van Loan matrix exponential, for
      uniform and non-uniform steps, and their semigroup composition law;
    * Q positive semidefinite;
    * the doubling stationary covariance against a Kronecker Lyapunov solve and
      against quadrature of the analytic h, v and a spectra;
    * discrete covariance recursion converging to the continuous stationary
      covariance, in double and in float;
    * the kappa scaling law used to set the driving gain at run time;
    * a long exact-discretization simulation staying finite and bounded.
*/

#define EIGEN_NON_ARDUINO

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>

#include "kalman_carma/CarmaVerticalModel.h"

using namespace ocean_imu::kalman::carma;

namespace {

int g_failures = 0;
int g_checks = 0;

void check(bool ok, const char* what, double value = NAN, double limit = NAN) {
    ++g_checks;
    if (!ok) {
        ++g_failures;
        std::printf("FAIL: %s (value=%.6g limit=%.6g)\n", what, value, limit);
    }
}

constexpr double TWO_PI = 6.283185307179586476925286766559;

template<int N>
using MatD = Eigen::Matrix<double,N,N>;

template<int N>
MatD<N> unit_W() {
    MatD<N> W = MatD<N>::Zero();
    W(2,2) = 1.0;
    return W;
}

// Van Loan: exp([[-A, W],[0, A^T]] h) = [[*, E12],[0, E22]], F = E22^T, Q = F E12.
//
// Evaluated in the diagonally scaled coordinates x = D z, D = diag(h^d_i),
// with d_i the number of couplings between state i and the noise input.  An
// entry two integrations from the noise is O(h^5) in Q; unscaled, the Pade
// rounding of the 2N x 2N exponential is relative to its largest entry and
// would swamp it.  Scaled, every entry of the reference carries its own
// relative precision.  Used for h <= 0.05 s, where exp(-A h) stays O(1).
template<int N>
void van_loan(const MatD<N>& A, const MatD<N>& W, double h, MatD<N>& F, MatD<N>& Q) {
    int dist[N];
    for (int i = 0; i < N; ++i) dist[i] = (W(i,i) != 0.0) ? 0 : 1000;
    for (int pass = 0; pass < N; ++pass)
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                if (A(i,j) != 0.0 && dist[j] + 1 < dist[i]) dist[i] = dist[j] + 1;
    Eigen::Matrix<double,N,1> d;
    for (int i = 0; i < N; ++i) d(i) = std::pow(h, dist[i]);
    const MatD<N> Az = d.cwiseInverse().asDiagonal() * A * d.asDiagonal();
    const MatD<N> Wz = d.cwiseInverse().asDiagonal() * W * d.cwiseInverse().asDiagonal();
    Eigen::Matrix<double,2*N,2*N> M = Eigen::Matrix<double,2*N,2*N>::Zero();
    M.template topLeftCorner<N,N>() = -Az * h;
    M.template topRightCorner<N,N>() = Wz * h;
    M.template bottomRightCorner<N,N>() = Az.transpose() * h;
    const Eigen::Matrix<double,2*N,2*N> E = M.exp();
    const MatD<N> Fz = E.template bottomRightCorner<N,N>().transpose();
    MatD<N> Qz = Fz * E.template topRightCorner<N,N>();
    Qz = 0.5 * (Qz + Qz.transpose()).eval();
    F = d.asDiagonal() * Fz * d.cwiseInverse().asDiagonal();
    Q = d.asDiagonal() * Qz * d.asDiagonal();
}

// Reference for any step: scaled Van Loan up to 0.05 s, beyond it the
// stationary identity Q = P - F P F^T with F = exp(A h) (cancellation in that
// identity is at most F ~ I, i.e. small h, where Van Loan is used instead).
template<int N>
void reference_discretization(const MatD<N>& A, const MatD<N>& W, const MatD<N>& Pinf,
                              double h, MatD<N>& F, MatD<N>& Q) {
    if (h <= 0.05) {
        van_loan<N>(A, W, h, F, Q);
        return;
    }
    F = (A * h).exp();
    Q = Pinf - F * Pinf * F.transpose();
    Q = 0.5 * (Q + Q.transpose()).eval();
}

// Kronecker solve of A P + P A^T + W = 0.
template<int N>
MatD<N> lyapunov_kron(const MatD<N>& A, const MatD<N>& W) {
    Eigen::Matrix<double,N*N,N*N> K = Eigen::Matrix<double,N*N,N*N>::Zero();
    const MatD<N> I = MatD<N>::Identity();
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k)
                for (int q = 0; q < N; ++q)
                    // vec column-major: (I (x) A + A (x) I) vec(P)
                    K(i + N * j, k + N * q) = I(j,q) * A(i,k) + A(j,q) * I(i,k);
    Eigen::Matrix<double,N*N,1> w;
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) w(i + N * j) = -W(i,j);
    const Eigen::Matrix<double,N*N,1> p = K.fullPivLu().solve(w);
    MatD<N> P;
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) P(i,j) = p(i + N * j);
    return 0.5 * (P + P.transpose());
}

// Error measured against each entry's natural scale sqrt(Q_ii Q_jj).
template<int N, typename M1, typename M2>
double normalized_error(const M1& X, const M2& Ref) {
    double e = 0.0;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            const double scale = std::sqrt(std::abs(Ref(i,i) * Ref(j,j))) + 1e-300;
            e = std::max(e, std::abs(static_cast<double>(X(i,j)) - Ref(i,j)) / scale);
        }
    return e;
}

template<int N>
std::vector<std::complex<double>> shape_poles(const CarmaShape& s, double kappa) {
    std::vector<std::complex<double>> p;
    auto quad = [&](double f, double z) {
        const double w = TWO_PI * kappa * f;
        const std::complex<double> disc = std::sqrt(std::complex<double>(z * z - 1.0, 0.0));
        p.push_back(w * (-z + disc));
        p.push_back(w * (-z - disc));
    };
    quad(s.f1_hz, s.zeta1);
    if (N >= 4) quad(s.f2_hz, s.zeta2);
    if (N != 4) p.push_back(-TWO_PI * kappa * s.f_lambda_hz);
    return p;
}

template<int N>
void test_structure(const char* name) {
    const CarmaShape sh = CarmaDefaultShape<N>::value;
    for (double kappa : {0.45, 1.0, 2.3}) {
        MatD<N> A;
        carma_dynamics<N>(sh, kappa, A);

        // h' = v and v' = a exactly; noise enters a' only.
        MatD<N> W = unit_W<N>();
        const Eigen::Matrix<double,N,1> G = Eigen::Matrix<double,N,1>::Unit(2);
        const Eigen::Matrix<double,1,N> C = Eigen::Matrix<double,1,N>::Unit(0);
        check(std::abs((C * G)(0)) == 0.0, "C G = 0");
        check(std::abs((C * A * G)(0)) == 0.0, "C A G = 0");
        check(std::abs((C * A * A * G)(0) - 1.0) < 1e-15, "C A^2 G = g");
        check(((C * A) - Eigen::Matrix<double,1,N>::Unit(1)).cwiseAbs().maxCoeff() == 0.0,
              "C A = velocity row");
        check(((C * A * A) - Eigen::Matrix<double,1,N>::Unit(2)).cwiseAbs().maxCoeff() == 0.0,
              "C A^2 = acceleration row");
        (void)W;

        // Poles: strictly stable and equal to the factor roots.
        Eigen::EigenSolver<MatD<N>> es(A);
        std::vector<std::complex<double>> ev(es.eigenvalues().data(),
                                             es.eigenvalues().data() + N);
        std::vector<std::complex<double>> ref = shape_poles<N>(sh, kappa);
        double max_re = -1e300;
        for (auto& e : ev) max_re = std::max(max_re, e.real());
        check(max_re < 0.0, "poles in open left half plane", max_re, 0.0);
        auto by = [](const std::complex<double>& a, const std::complex<double>& b) {
            return a.real() != b.real() ? a.real() < b.real() : a.imag() < b.imag();
        };
        std::sort(ev.begin(), ev.end(), by);
        std::sort(ref.begin(), ref.end(), by);
        double perr = 0.0;
        for (int i = 0; i < N; ++i) perr = std::max(perr, std::abs(ev[i] - ref[i]) / std::abs(ref[i]));
        check(perr < 1e-8, "eigenvalues equal factor poles", perr, 1e-8);

        // Transfer function against the polynomial PSD.
        double terr = 0.0;
        for (int k = 0; k <= 80; ++k) {
            const double f = 0.005 * std::pow(10.0, 3.6 * k / 80.0);  // 5 mHz .. 20 Hz
            const std::complex<double> s(0.0, TWO_PI * f);
            Eigen::Matrix<std::complex<double>,N,N> M =
                s * Eigen::Matrix<std::complex<double>,N,N>::Identity() - A.template cast<std::complex<double>>();
            const Eigen::Matrix<std::complex<double>,N,1> x =
                M.fullPivLu().solve(G.template cast<std::complex<double>>());
            const double ss = std::norm(x(0));
            const double ref_psd = carma_unit_psd<N>(sh, kappa, f);
            terr = std::max(terr, std::abs(ss - ref_psd) / ref_psd);
        }
        check(terr < 1e-9, "realization PSD equals g^2|B|^2/|P|^2", terr, 1e-9);
    }
    std::printf("[%s] structure, poles and transfer function checked\n", name);
}

// Integral of S over the real line, Var = 2 int_0^inf S(2 pi f) df, on a log grid.
template<int N>
double variance_quadrature(const CarmaShape& sh, int moment) {
    // Composite Simpson in u = log f over [1e-7, 1e7] Hz plus [0, 1e-7].
    const double lo = std::log(1e-7), hi = std::log(1e7);
    const int n = 400000;
    const double du = (hi - lo) / n;
    double s = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double u = lo + i * du;
        const double f = std::exp(u);
        const double w = TWO_PI * f;
        const double val = carma_unit_psd<N>(sh, 1.0, f) * std::pow(w, 2 * moment) * f;
        const double c = (i == 0 || i == n) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
        s += c * val;
    }
    s *= du / 3.0;
    s += carma_unit_psd<N>(sh, 1.0, 0.0) * (moment == 0 ? 1e-7 : 0.0);
    return 2.0 * s;
}

template<int N>
void test_stationary(const char* name) {
    const CarmaShape sh = CarmaDefaultShape<N>::value;
    MatD<N> A;
    carma_dynamics<N>(sh, 1.0, A);
    const MatD<N> W = unit_W<N>();
    MatD<N> P;
    check(stationary_covariance_doubling<double,N>(A, W, P), "doubling converged");
    const MatD<N> Pk = lyapunov_kron<N>(A, W);
    const double e_k = normalized_error<N>(P, Pk);
    check(e_k < 1e-9, "doubling equals Kronecker Lyapunov solve", e_k, 1e-9);
    const double resid = (A * P + P * A.transpose() + W).cwiseAbs().maxCoeff();
    check(resid < 1e-10, "Lyapunov residual", resid, 1e-10);
    Eigen::SelfAdjointEigenSolver<MatD<N>> es(P);
    check(es.eigenvalues().minCoeff() > 0.0, "stationary covariance positive definite",
          es.eigenvalues().minCoeff(), 0.0);
    for (int m = 0; m < 3; ++m) {
        const double q = variance_quadrature<N>(sh, m);
        const double rel = std::abs(q - P(m,m)) / P(m,m);
        check(rel < 1e-5, m == 0 ? "Var h = int S_h" : (m == 1 ? "Var v = int w^2 S_h"
                                                              : "Var a = int w^4 S_h"),
              rel, 1e-5);
    }
    // Finite acceleration variance, no displacement/acceleration cross term
    // beyond the stationary identities E[h v] = 0, E[v a] = 0, E[h a] = -Var v.
    check(std::abs(P(0,1)) < 1e-12 * std::sqrt(P(0,0) * P(1,1)), "E[h v] = 0", P(0,1));
    check(std::abs(P(1,2)) < 1e-12 * std::sqrt(P(1,1) * P(2,2)), "E[v a] = 0", P(1,2));
    check(std::abs(P(0,2) + P(1,1)) < 1e-9 * P(1,1), "E[h a] = -Var v", P(0,2) + P(1,1));
    std::printf("[%s] Var h=%.6e Var v=%.6e Var a=%.6e (unit gain); f_z=%.5f Hz\n",
                name, P(0,0), P(1,1), P(2,2),
                std::sqrt(P(1,1) / P(0,0)) / TWO_PI);
}

template<int N>
void test_discretization(const char* name) {
    const CarmaShape sh = CarmaDefaultShape<N>::value;
    for (double kappa : {0.6, 1.0, 1.6}) {
        MatD<N> A;
        carma_dynamics<N>(sh, kappa, A);
        const MatD<N> W = unit_W<N>() * 0.37;
        MatD<N> Pinf;
        stationary_covariance_doubling<double,N>(A, W, Pinf);
        for (double h : {1e-4, 2.5e-3, 5e-3, 0.0137, 0.05, 0.4, 3.0, 40.0}) {
            MatD<N> F, Q, Fr, Qr;
            check(discretize_lti<double,N>(A, W, h, F, Q), "discretize ok");
            reference_discretization<N>(A, W, Pinf, h, Fr, Qr);
            const double ef = (F - Fr).cwiseAbs().maxCoeff() / std::max(1.0, Fr.cwiseAbs().maxCoeff());
            const double eq = normalized_error<N>(Q, Qr);
            check(ef < 1e-11, "F equals Van Loan exp(A h)", ef, 1e-11);
            check(eq < 1e-9, "Q equals Van Loan integral", eq, 1e-9);

            // PSD of Q.
            Eigen::SelfAdjointEigenSolver<MatD<N>> es(Q);
            const double lmin = es.eigenvalues().minCoeff();
            const double lmax = es.eigenvalues().maxCoeff();
            check(lmin >= -1e-14 * lmax, "Q positive semidefinite (double)", lmin / lmax, -1e-14);

            // Float implementation against the double reference.  F is
            // compared in stationary-whitened coordinates S^-1 F S, S =
            // diag(sqrt(Pinf_ii)), where its entries are dimensionless and
            // O(1); Q entrywise against its own diagonal scale.
            Eigen::Matrix<float,N,N> Ff, Qf;
            discretize_lti<float,N>(A.template cast<float>(), W.template cast<float>(),
                                    static_cast<float>(h), Ff, Qf);
            const Eigen::Matrix<double,N,1> sd = Pinf.diagonal().cwiseSqrt();
            const MatD<N> dFw = sd.cwiseInverse().asDiagonal()
                                * (Ff.template cast<double>() - Fr) * sd.asDiagonal();
            const double eff = dFw.cwiseAbs().maxCoeff();
            const double eqf = normalized_error<N>(Qf.template cast<double>(), Qr);
            // Float rounding of A itself already moves the whitened F by up to
            // ~1e-5 for CARMA(5,2), whose jerk row mixes units over three
            // decades; these bounds characterize that, they are not slack.
            const double tol_f = (h <= 0.05) ? 1e-5 : 1e-4;
            const double tol_q = (h <= 0.05) ? 1e-4 : 2e-4;
            check(eff < tol_f, "float F (whitened)", eff, tol_f);
            check(eqf < tol_q, "float Q", eqf, tol_q);
            // Positive semidefinite to float rounding after whitening by its
            // own diagonal.
            const Eigen::Matrix<double,N,1> qd = Qf.template cast<double>().diagonal().cwiseSqrt();
            const MatD<N> Qfw = qd.cwiseInverse().asDiagonal() * Qf.template cast<double>()
                                * qd.cwiseInverse().asDiagonal();
            Eigen::SelfAdjointEigenSolver<MatD<N>> esf(Qfw);
            check(esf.eigenvalues().minCoeff() > -1e-6, "float Q positive semidefinite (correlation form)",
                  esf.eigenvalues().minCoeff(), -1e-6);
            if (std::getenv("CARMA_TEST_VERBOSE")) {
                std::printf("  kappa=%.2f h=%-7g double eF=%.1e eQ=%.1e | float eFw=%.1e eQ=%.1e corr-min-eig=%.1e\n",
                            kappa, h, ef, eq, eff, eqf, esf.eigenvalues().minCoeff());
            }
        }

        // Non-uniform steps compose exactly (semigroup).
        const double h1 = 0.0031, h2 = 0.0071;
        MatD<N> F1, Q1, F2, Q2, F12, Q12;
        discretize_lti<double,N>(A, W, h1, F1, Q1);
        discretize_lti<double,N>(A, W, h2, F2, Q2);
        discretize_lti<double,N>(A, W, h1 + h2, F12, Q12);
        const double esf = (F12 - F2 * F1).cwiseAbs().maxCoeff();
        const double esq = normalized_error<N>(MatD<N>(F2 * Q1 * F2.transpose() + Q2), Q12);
        check(esf < 1e-13, "F(h1+h2) = F(h2) F(h1)", esf, 1e-13);
        check(esq < 1e-11, "Q(h1+h2) = F(h2) Q(h1) F(h2)^T + Q(h2)", esq, 1e-11);
    }
    std::printf("[%s] F, Q against Van Loan; PSD; semigroup; float precision checked\n", name);
}

template<int N>
void test_recursion_converges(const char* name) {
    CarmaVerticalModel<double,N> md;
    CarmaVerticalModel<float,N> mf;
    md.configure(1.3, 0.8);
    mf.configure(1.3f, 0.8f);
    MatD<N> Pinf;
    md.stationary_covariance(Pinf);

    const double dt = 0.005;
    MatD<N> F, Q;
    md.discretize(dt, F, Q);
    Eigen::Matrix<float,N,N> Ff, Qf;
    mf.discretize(static_cast<float>(dt), Ff, Qf);
    MatD<N> P = MatD<N>::Zero();
    Eigen::Matrix<float,N,N> Pf = Eigen::Matrix<float,N,N>::Zero();
    for (int k = 0; k < 120000; ++k) {  // 600 s
        P = (F * P * F.transpose() + Q).eval();
        Pf = (Ff * Pf * Ff.transpose() + Qf).eval();
    }
    const double e_d = normalized_error<N>(P, Pinf);
    const double e_f = normalized_error<N>(Pf.template cast<double>(), Pinf);
    check(e_d < 1e-9, "discrete recursion -> continuous stationary covariance (double)", e_d, 1e-9);
    check(e_f < 2e-3, "discrete recursion -> continuous stationary covariance (float)", e_f, 2e-3);

    // Float stationary covariance from the doubling routine.
    Eigen::Matrix<float,N,N> Pfd;
    check(mf.stationary_covariance(Pfd), "float doubling converged");
    const double e_fd = normalized_error<N>(Pfd.template cast<double>(), Pinf);
    check(e_fd < 1e-4, "float doubling stationary covariance", e_fd, 1e-4);
    std::printf("[%s] recursion error double=%.2e float=%.2e, float doubling=%.2e\n",
                name, e_d, e_f, e_fd);
}

template<int N>
void test_scaling(const char* name) {
    CarmaVerticalModel<double,N> m;
    const double fz1 = m.base_zero_crossing_hz();
    for (double kappa : {0.5, 1.0, 1.7, 3.1}) {
        for (double sa : {0.02, 0.6, 3.0}) {
            check(m.configure(kappa, sa), "configure");
            MatD<N> P;
            m.stationary_covariance(P);
            const double ea = std::abs(std::sqrt(P(2,2)) - sa) / sa;
            check(ea < 1e-9, "configured stationary accel std", ea, 1e-9);
            const double fz = std::sqrt(P(1,1) / P(0,0)) / TWO_PI;
            const double ez = std::abs(fz - kappa * fz1) / (kappa * fz1);
            check(ez < 1e-9, "zero-crossing frequency scales with kappa", ez, 1e-9);
            const double eh = std::abs(m.displacement_std() - std::sqrt(P(0,0))) / std::sqrt(P(0,0));
            check(eh < 1e-9, "displacement std scaling law", eh, 1e-9);
        }
    }
    check(!m.configure(-1.0, 1.0) && !m.configure(1.0, 0.0) && !m.configure(NAN, 1.0),
          "configure rejects invalid scale or amplitude");
    std::printf("[%s] base displacement zero-crossing %.5f Hz; scaling law checked\n", name, fz1);
}

template<int N>
void test_long_simulation(const char* name) {
    CarmaVerticalModel<double,N> m;
    m.configure(1.0, 1.0);
    MatD<N> P, F, Q;
    m.stationary_covariance(P);
    const double dt = 0.005;
    m.discretize(dt, F, Q);
    Eigen::LDLT<MatD<N>> ldlt(Q);
    const MatD<N> L = ldlt.matrixL();
    Eigen::Matrix<double,N,1> D = ldlt.vectorD().cwiseMax(0.0).cwiseSqrt();
    std::mt19937_64 rng(20261025u);
    std::normal_distribution<double> n01;
    Eigen::LLT<MatD<N>> llt0(P);
    Eigen::Matrix<double,N,1> z, x;
    for (int i = 0; i < N; ++i) z(i) = n01(rng);
    x = llt0.matrixL() * z;
    double max_h = 0.0;
    Eigen::Matrix<double,3,1> sumsq = Eigen::Matrix<double,3,1>::Zero();
    const int steps = 2000000;  // 10000 s
    bool finite = true;
    for (int k = 0; k < steps; ++k) {
        for (int i = 0; i < N; ++i) z(i) = D(i) * n01(rng);
        Eigen::Matrix<double,N,1> w = ldlt.transpositionsP().transpose() * (L * z);
        x = (F * x + w).eval();
        if (!x.allFinite()) { finite = false; break; }
        max_h = std::max(max_h, std::abs(x(0)));
        for (int i = 0; i < 3; ++i) sumsq(i) += x(i) * x(i);
    }
    check(finite, "long simulation finite");
    check(max_h < 8.0 * std::sqrt(P(0,0)), "long simulation bounded (8 sigma)", max_h / std::sqrt(P(0,0)), 8.0);
    for (int i = 0; i < 3; ++i) {
        const double ratio = sumsq(i) / steps / P(i,i);
        check(std::abs(ratio - 1.0) < 0.15, "sample variance near stationary variance", ratio, 1.0);
    }
    std::printf("[%s] 10000 s simulation: max|h|=%.2f sigma, var ratios h=%.3f v=%.3f a=%.3f\n",
                name, max_h / std::sqrt(P(0,0)), sumsq(0) / steps / P(0,0),
                sumsq(1) / steps / P(1,1), sumsq(2) / steps / P(2,2));
}

template<int N>
void run_all(const char* name) {
    test_structure<N>(name);
    test_stationary<N>(name);
    test_discretization<N>(name);
    test_recursion_converges<N>(name);
    test_scaling<N>(name);
    test_long_simulation<N>(name);
}

}  // namespace

int main() {
    run_all<3>("CARMA(3,0)");
    run_all<4>("CARMA(4,1)");
    run_all<5>("CARMA(5,2)");

    // A zero step is the identity with no process noise; invalid steps are rejected.
    {
        MatD<5> A, F, Q;
        carma_dynamics<5>(CARMA52_JONSWAP_SHAPE, 1.0, A);
        const MatD<5> W = unit_W<5>();
        check(discretize_lti<double,5>(A, W, 0.0, F, Q) && F.isIdentity() && Q.isZero(),
              "zero step is identity");
        check(!discretize_lti<double,5>(A, W, -1.0, F, Q) && F.isIdentity() && Q.isZero(),
              "negative step rejected");
        check(!discretize_lti<double,5>(A, W, NAN, F, Q), "NaN step rejected");
    }

    std::printf("carma_model-test: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
