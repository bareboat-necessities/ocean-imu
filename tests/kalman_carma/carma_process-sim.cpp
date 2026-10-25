/*
  Copyright (c) 2026 Mikhail Grushinskiy

  carma_process-sim -- long exact-discretization realizations of the
  experimental CARMA vertical priors, for the independent spectral validation
  in tools/carma_process_validation.py.

  For each order the library model (CarmaVerticalModel<float, N>, kappa = 1,
  unit stationary acceleration std) supplies F and Q; the state is propagated
  in double with Gaussian increments w ~ N(0, Q).  Writes little-endian
  float64 triples (h, v, a):

    <out>/carma<N>_coarse.f64   dt = 0.02 s, 2^21 samples (11.6 h)
    <out>/carma<N>_fine.f64     dt = 0.001 s, 400000 samples (400 s)

  Usage: carma_process-sim OUTDIR [SEED]
*/

#define EIGEN_NON_ARDUINO

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "kalman_carma/CarmaVerticalModel.h"

using namespace ocean_imu::kalman::carma;

namespace {

template<int N>
bool simulate(const std::string& path, double dt, long n, unsigned seed) {
    CarmaVerticalModel<float, N> m;
    m.configure(1.0f, 1.0f);
    Eigen::Matrix<float, N, N> Ff, Qf, Pf;
    if (!m.discretize(static_cast<float>(dt), Ff, Qf) || !m.stationary_covariance(Pf)) return false;
    const Eigen::Matrix<double, N, N> F = Ff.template cast<double>();
    const Eigen::Matrix<double, N, N> Q = Qf.template cast<double>();
    Eigen::LDLT<Eigen::Matrix<double, N, N>> ldlt(Q);
    const Eigen::Matrix<double, N, N> L = ldlt.matrixL();
    const Eigen::Matrix<double, N, 1> D = ldlt.vectorD().cwiseMax(0.0).cwiseSqrt();
    Eigen::LLT<Eigen::Matrix<double, N, N>> llt(Pf.template cast<double>());
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> n01;
    Eigen::Matrix<double, N, 1> z, x;
    for (int i = 0; i < N; ++i) z(i) = n01(rng);
    x = llt.matrixL() * z;
    std::vector<double> buf;
    buf.reserve(static_cast<size_t>(3 * n));
    for (long k = 0; k < n; ++k) {
        buf.push_back(x(0)); buf.push_back(x(1)); buf.push_back(x(2));
        for (int i = 0; i < N; ++i) z(i) = D(i) * n01(rng);
        x = (F * x + ldlt.transpositionsP().transpose() * (L * z)).eval();
        if (!x.allFinite()) return false;
    }
    std::ofstream out(path, std::ios::binary);
    out.write(reinterpret_cast<const char*>(buf.data()),
              static_cast<std::streamsize>(buf.size() * sizeof(double)));
    return static_cast<bool>(out);
}

template<int N>
bool both(const std::string& dir, unsigned seed) {
    const std::string stem = dir + "/carma" + std::to_string(N);
    return simulate<N>(stem + "_coarse.f64", 0.02, 1L << 21, seed) &&
           simulate<N>(stem + "_fine.f64", 0.001, 400000, seed + 1u);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "Usage: %s OUTDIR [SEED]\n", argv[0]);
        return 2;
    }
    const std::string dir = argv[1];
    const unsigned seed = argc > 2 ? static_cast<unsigned>(std::strtoul(argv[2], nullptr, 10)) : 20261025u;
    const bool ok = both<3>(dir, seed) && both<4>(dir, seed + 10u) && both<5>(dir, seed + 20u);
    std::printf("carma_process-sim: %s\n", ok ? "ok" : "FAILED");
    return ok ? 0 : 1;
}
