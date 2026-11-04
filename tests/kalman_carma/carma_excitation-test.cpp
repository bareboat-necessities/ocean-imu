/*
  Copyright (c) 2026 Mikhail Grushinskiy

  carma_excitation-test -- the wave-excitation evidence of the CARMA research
  wrappers (src/kalman_carma/CarmaWaveEvidence.h,
  docs/carma-wave-excitation.md):

    * the sensor-noise model of the evidence statistic against the shipping
      measurement chain (AdaptiveWaveBandPass -> SeaStateAutoTuner, and
      WavePeriodEstimator) driven by simulated white noise plus a bias random
      walk: mean band variance, its effective sample size, and the
      elevation-proxy variance;
    * the false-alarm rate of the evidence weight on pure noise;
    * the CARMA process covariance as the driving noise goes to zero: Q and
      the stationary covariance scale exactly with sigma^2, stay symmetric
      and positive semidefinite, W = 0 gives Q = 0 with an unchanged
      transition, and the numerical floor keeps Q out of float subnormals at
      constant and jittered steps;
    * the Legacy mode's floor and the Evidence mode on a still IMU;
    * repeated wave / no-wave transitions through both wrappers (quaternion
      MEKF and two-frame Lie group): finite, symmetric, positive
      semidefinite covariance at every check, the excitation leaving and
      re-entering, the heave following the waves, no reset.
*/

#define EIGEN_NON_ARDUINO

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "kalman_carma/CarmaWaveEvidence.h"
#include "kalman_carma/CarmaVerticalModel.h"
#include "tuner/AdaptiveWaveBandPass.h"
#include "tuner/SeaStateAutoTuner.h"
#include "tuner/WavePeriodEstimator.h"
#include "CarmaBenchTruth.h"
#include "carma3/FullAdapters.h"
#include "carma3/FullAdaptersTFG.h"

namespace kc = ocean_imu::kalman::carma;
using ocean_imu::carma::CarmaEvidenceConfig;
using ocean_imu::carma::CarmaExcitationMode;
using ocean_imu::carma::CarmaSensorNoise;
using ocean_imu::carma::CarmaWaveEvidence;

namespace {

int g_fail = 0;

void check(const char* what, bool ok) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok) ++g_fail;
}

// ---------------------------------------------------------------------------
// Sensor-noise model against the measurement chain
// ---------------------------------------------------------------------------

struct NoiseRun {
    double var_mean = 0, var_relstd = 0, var0 = 0, n_eff = 0, elev_var = 0, elev0 = 0;
    double p_frac = 0, exc_mean = 0;
};

NoiseRun run_noise(float sw, float qrw, float fref, double T, unsigned seed) {
    const float dt = 0.005f;
    std::mt19937 rng(seed);
    std::normal_distribution<float> n01(0.0f, 1.0f);
    WavePeriodEstimator wpe;
    AdaptiveWaveBandPass band(0.5f, 4.0f, 0.01f, 6.0f);
    SeaStateAutoTuner tuner;
    CarmaWaveEvidence ev;
    const CarmaSensorNoise model{sw, qrw, 0.5f};
    const CarmaEvidenceConfig cfg{};
    float b = 0.0f;
    std::vector<double> V;
    double se = 0, sv0 = 0, sn = 0, spf = 0, sexc = 0;
    long n = 0;
    const long steps = static_cast<long>(T / dt);
    for (long k = 0; k < steps; ++k) {
        b += qrw * std::sqrt(dt) * n01(rng);
        const float x = sw * n01(rng) + b;
        wpe.update(dt, x);
        tuner.update(dt, band.step(x, dt, fref), fref);
        CarmaWaveEvidence::Inputs in;
        in.dt = dt;
        in.accel = x;
        in.band_low_hz = band.lowHz();
        in.band_high_hz = band.highHz();
        in.band_white_gain = band.whiteNoiseVarianceGain();
        in.var_horizon_s = tuner.getVarianceHorizonSec();
        ev.update(in, model, cfg);
        const float es = wpe.getElevationStd();
        if (static_cast<double>(k) * dt > 600.0 && k % 20 == 0) {
            V.push_back(tuner.getAccelVariance());
            se += static_cast<double>(es) * es;
            sv0 += ev.noiseVariance();
            sn += ev.effectiveDof();
            spf += ev.waveProbability() >= 0.5f ? 1.0 : 0.0;
            sexc += ev.excitationVariance(true);
            ++n;
        }
    }
    NoiseRun r;
    double m = 0, m2 = 0;
    for (double v : V) { m += v; m2 += v * v; }
    m /= static_cast<double>(V.size());
    m2 = m2 / static_cast<double>(V.size()) - m * m;
    r.var_mean = m;
    r.var_relstd = std::sqrt(std::max(0.0, m2)) / m;
    r.var0 = sv0 / static_cast<double>(n);
    r.n_eff = sn / static_cast<double>(n);
    r.elev_var = se / static_cast<double>(n);
    // WavePeriodEstimator's default leak corner, 0.02 Hz.
    r.elev0 = CarmaWaveEvidence::elevationNoise(model, dt, 0.02f);
    r.p_frac = spf / static_cast<double>(n);
    r.exc_mean = sexc / static_cast<double>(n);
    return r;
}

void test_noise_model() {
    const float sw = 0.0148f;
    for (float qrw : {0.0f, 5e-4f}) {
        for (float fref : {0.05f, 0.3f, 1.0f}) {
            const NoiseRun r = run_noise(sw, qrw, fref, 6000.0, 17u + static_cast<unsigned>(fref * 100));
            const double nemp = 2.0 / (r.var_relstd * r.var_relstd);
            std::printf("[noise] q=%.0e f=%.2f: E[V] %.3e model %.3e (%.1f%%), N_eff %.0f model %.0f, "
                        "elev %.3e model %.3e, p>=0.5 %.2f%%, mean p E+ / E[V|0] %.3f\n",
                        qrw, fref, r.var_mean, r.var0, 100.0 * (r.var_mean / r.var0 - 1.0), nemp, r.n_eff,
                        r.elev_var, r.elev0, 100.0 * r.p_frac, r.exc_mean / r.var0);
            char buf[160];
            std::snprintf(buf, sizeof buf, "band variance of noise matches the model within 8%% (q=%.0e, f=%.2f)", qrw, fref);
            check(buf, std::abs(r.var_mean / r.var0 - 1.0) < 0.08);
            std::snprintf(buf, sizeof buf, "effective sample size matches the model within a factor 1.5 (q=%.0e, f=%.2f)", qrw, fref);
            check(buf, nemp / r.n_eff < 1.5 && r.n_eff / nemp < 1.5);
            std::snprintf(buf, sizeof buf, "evidence weight >= 0.5 on <= 2%% of pure-noise samples (q=%.0e, f=%.2f)", qrw, fref);
            check(buf, r.p_frac <= 0.02);
            if (fref == 0.3f) {
                std::snprintf(buf, sizeof buf, "elevation-proxy noise variance matches the model within 15%% (q=%.0e)", qrw);
                check(buf, std::abs(r.elev_var / r.elev0 - 1.0) < 0.15);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Process covariance as the driving noise goes to zero
// ---------------------------------------------------------------------------

template<int N>
double min_eig_rel(const Eigen::Matrix<float,N,N>& Q) {
    const Eigen::Matrix<double,N,N> Qd = Q.template cast<double>();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,N,N>> es(Qd);
    const double mx = es.eigenvalues().cwiseAbs().maxCoeff();
    return mx > 0 ? es.eigenvalues().minCoeff() / mx : 0.0;
}

template<int N>
bool no_subnormal(const Eigen::Matrix<float,N,N>& Q) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            const float v = std::abs(Q(i,j));
            if (v != 0.0f && v < FLT_MIN) return false;
        }
    return true;
}

template<int N>
void test_q_scaling(const char* label) {
    kc::CarmaVerticalModel<float,N> m1, ms;
    const double fz = m1.base_zero_crossing_hz();
    const float dt = 0.005f;
    double worst_rel = 0, worst_eig = 0, worst_stat = 0;
    bool sym = true, finite = true, subn = true;
    for (float f_tune : {0.03f, 0.05f, 0.3f, 1.2f}) {
        const float kappa = static_cast<float>(f_tune / fz);
        m1.configure(kappa, 1.0f);
        Eigen::Matrix<float,N,N> F1, Q1, P1;
        m1.discretize(dt, F1, Q1);
        m1.stationary_covariance(P1);
        for (float sigma : {1e-2f, 1e-3f, 1e-5f}) {
            ms.configure(kappa, sigma);
            Eigen::Matrix<float,N,N> Fs, Qs, Ps;
            ms.discretize(dt, Fs, Qs);
            ms.stationary_covariance(Ps);
            const float s2 = sigma * sigma;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    sym = sym && Qs(i,j) == Qs(j,i);
                    finite = finite && std::isfinite(Qs(i,j)) && std::isfinite(Ps(i,j));
                    const double ref = static_cast<double>(Q1(i,j)) * s2;
                    const double scale = std::sqrt(static_cast<double>(Q1(i,i)) * Q1(j,j)) * s2;
                    if (scale > 0) worst_rel = std::max(worst_rel, std::abs(Qs(i,j) - ref) / scale);
                    const double pref = static_cast<double>(P1(i,j)) * s2;
                    const double ps = std::sqrt(static_cast<double>(P1(i,i)) * P1(j,j)) * s2;
                    if (ps > 0) worst_stat = std::max(worst_stat, std::abs(Ps(i,j) - pref) / ps);
                }
            worst_eig = std::min(worst_eig, min_eig_rel<N>(Qs));
            subn = subn && no_subnormal<N>(Qs);
            if (!(Fs - F1).isZero(0.0f)) finite = false;   // Phi does not depend on sigma
        }
    }
    std::printf("[Q->0] %s: max |Q(s)/s^2 - Q(1)| %.2e, max |P(s)/s^2 - P(1)| %.2e (correlation units), "
                "min eig / max %.2e\n", label, worst_rel, worst_stat, worst_eig);
    char buf[160];
    std::snprintf(buf, sizeof buf, "%s: Q and P_inf scale with sigma^2 (float, 1e-2 .. 1e-5 m/s^2)", label);
    check(buf, worst_rel < 1e-4 && worst_stat < 1e-3);
    std::snprintf(buf, sizeof buf, "%s: Q symmetric, finite, PSD, Phi independent of sigma", label);
    check(buf, sym && finite && worst_eig > -1e-5);
    std::snprintf(buf, sizeof buf, "%s: no subnormal Q entry down to the 1e-5 m/s^2 numerical floor", label);
    check(buf, subn);

    // W = 0 exactly: a deterministic, exponentially stable prior.
    m1.configure(static_cast<float>(0.05 / fz), 1.0f);
    const Eigen::Matrix<float,N,N> A = m1.A(), Z = Eigen::Matrix<float,N,N>::Zero();
    Eigen::Matrix<float,N,N> F0, Q0, Fw, Qw, P0;
    const bool ok0 = kc::discretize_lti<float,N>(A, Z, dt, F0, Q0);
    kc::discretize_lti<float,N>(A, m1.W(), dt, Fw, Qw);
    const bool st0 = kc::stationary_covariance_doubling<float,N>(A, Z, P0);
    std::snprintf(buf, sizeof buf, "%s: W = 0 gives Q = 0, the same Phi, and P_inf = 0", label);
    check(buf, ok0 && st0 && Q0.isZero(0.0f) && P0.isZero(0.0f) && (F0 - Fw).cwiseAbs().maxCoeff() < 1e-6f);
    // Covariance under W = 0 decays from any PSD start: P <- Phi P Phi^T.
    Eigen::Matrix<float,N,N> P = Qw * 1e4f + Eigen::Matrix<float,N,N>::Identity() * 1e-3f;
    double rho = 0;
    {
        Eigen::EigenSolver<Eigen::Matrix<double,N,N>> es(F0.template cast<double>());
        rho = es.eigenvalues().cwiseAbs().maxCoeff();
    }
    const float tr0 = P.trace();
    for (int k = 0; k < 200 * 600; ++k) P = (F0 * P * F0.transpose()).eval();
    std::printf("[Q->0] %s: spectral radius of Phi %.9f, trace P after 600 s of W = 0: %.3e of start\n",
                label, rho, P.trace() / tr0);
    std::snprintf(buf, sizeof buf, "%s: rho(Phi) < 1 and P decays under W = 0", label);
    check(buf, rho < 1.0 && P.allFinite() && P.trace() < 1e-3f * tr0);
}

// Jittered sample clock at the numerical floor.
void test_step_discretization_floor() {
    kc::CarmaVerticalModel<float,4> m;
    m.configure(static_cast<float>(0.03 / m.base_zero_crossing_hz()), 1e-5f);
    kc::StepDiscretization<float,4> sd;
    std::mt19937 rng(5);
    std::uniform_real_distribution<float> u(-0.25f, 0.25f);
    bool ok = true;
    double worst = 0;
    for (int k = 0; k < 20000; ++k) {
        const float h = 0.005f * (1.0f + u(rng));
        if (!sd.advance(m.generation(), h)) sd.refresh(m.A(), m.W(), m.generation(), h);
        const auto& Q = sd.Q();
        ok = ok && Q.allFinite() && (Q - Q.transpose()).isZero(0.0f) && no_subnormal<4>(Q);
        worst = std::min(worst, min_eig_rel<4>(Q));
    }
    std::printf("[Q->0] jittered +-25 %% clock at sigma = 1e-5 m/s^2, f = 0.03 Hz: min eig / max %.2e\n", worst);
    check("jittered steps at the numerical floor: Q finite, symmetric, PSD, no subnormals", ok && worst > -1e-4);
}

// ---------------------------------------------------------------------------
// Wrappers
// ---------------------------------------------------------------------------

// One synthetic stream (bench truth and sensor models) through an adapter.
struct StreamCfg {
    std::string scenario;
    double T;
    unsigned seed;
};

template<class Adapter>
struct WrapperRun {
    std::vector<double> t, z_est, z_true, sig, pw;
    bool cov_ok = true;
    bool always_live = true;
    double max_jump = 0;   // largest one-sample change of the heave error
    double min_eig = 0;
};

template<class Adapter>
WrapperRun<Adapter> run_wrapper(Adapter& f, const carma_bench::Scenario& sc, unsigned seed, int check_every) {
    using namespace carma_bench;
    using namespace carma3_bench;
    const double dt = 1.0 / 200.0;
    ImuNoiseModel acc_noise = make_imu_noise_model(ACC_SIGMA, 8e-3f * 9.80665f, 0.0005f, w3d_expand_seed(seed, 11u),
                                                   w3d_expand_seed(seed, 12u));
    ImuNoiseModel gyr_noise = make_imu_noise_model(GYR_SIGMA, 0.05f * 3.14159265f / 180.0f, 0.00001f,
                                                   w3d_expand_seed(seed, 13u), w3d_expand_seed(seed, 14u));
    MagNoiseModel mag_noise = make_mag_noise_model(MAG_SIGMA_UT, 2.0f, 0.01f, 0.015f, 0.010f, 1.0f,
                                                   w3d_expand_seed(seed, 15u), w3d_expand_seed(seed, 16u));
    const Eigen::Vector3f mag_world = MagSim_WMM::mag_world_aero();
    WaveTruth truth(sc, w3d_expand_seed(seed, 21u), dt);
    WrapperRun<Adapter> r;
    float mag_phase = 0.0f;
    bool was_live = false;
    double z_prev = NAN;
    const long steps = std::lround(sc.duration_s / dt);
    for (long k = 0; k < steps; ++k) {
        const TruthSample ts = truth.next();
        Sensors s;
        s.dt = static_cast<float>(dt);
        s.acc_ned = apply_imu_noise(ts.f_b.cast<float>(), acc_noise, s.dt);
        s.gyr_ned = apply_imu_noise(ts.omega_b.cast<float>(), gyr_noise, s.dt);
        mag_phase += s.dt;
        if (mag_phase >= 0.04f) {
            mag_phase -= 0.04f;
            s.mag_tick = true;
            s.mag_ned = apply_mag_noise(ts.R_nb.transpose().cast<float>() * mag_world, mag_noise, 0.04f);
        }
        f.step(s);
        const Estimate3 e = f.estimate();
        if (was_live && !e.live) r.always_live = false;
        if (!e.live) continue;
        was_live = true;
        const double z = e.p.z();
        const double err = z - ts.p_wave.z();
        if (std::isfinite(z_prev)) r.max_jump = std::max(r.max_jump, std::abs(err - z_prev));
        z_prev = err;
        if (k % check_every == 0) {
            const auto& P = f.fusion().raw().mekf().covariance_full();
            const auto Pd = P.template cast<double>().eval();
            r.cov_ok = r.cov_ok && Pd.allFinite() && (Pd - Pd.transpose()).cwiseAbs().maxCoeff() <= 1e-12 * Pd.cwiseAbs().maxCoeff();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Pd, Eigen::EigenvaluesOnly);
            r.min_eig = std::min(r.min_eig, es.eigenvalues().minCoeff() / es.eigenvalues().cwiseAbs().maxCoeff());
            double d[8];
            f.excitation_diag(d);
            r.t.push_back(ts.t);
            r.z_est.push_back(z);
            r.z_true.push_back(ts.p_wave.z());
            r.sig.push_back(d[0]);
            r.pw.push_back(d[1]);
        }
    }
    return r;
}

template<class R>
void segment(const R& r, double t0, double t1, double& rms_err, double& rms_true, double& sig_max, double& sig_min) {
    double se = 0, st = 0;
    long n = 0;
    sig_max = 0;
    sig_min = INFINITY;
    for (size_t i = 0; i < r.t.size(); ++i) {
        if (r.t[i] < t0 || r.t[i] > t1) continue;
        se += (r.z_est[i] - r.z_true[i]) * (r.z_est[i] - r.z_true[i]);
        st += r.z_true[i] * r.z_true[i];
        sig_max = std::max(sig_max, r.sig[i]);
        sig_min = std::min(sig_min, r.sig[i]);
        ++n;
    }
    rms_err = n ? std::sqrt(se / static_cast<double>(n)) : NAN;
    rms_true = n ? std::sqrt(st / static_cast<double>(n)) : NAN;
}

template<class Adapter>
void test_transitions(Adapter& f, const char* label) {
    // Calm 0-300 s, waves 300-600 s, calm 600-900 s, waves 900-1200 s,
    // calm 1200-1500 s: four 10 s envelope transitions, one stream.
    carma_bench::Scenario sc;
    carma_bench::base_scenario("jonswap_fp0.29_hs0.5_g3.3", sc);
    sc.hmodel.physical = true;
    sc.duration_s = 1500.0;
    sc.env.level0 = 0.0;
    sc.env.steps = {{300.0, 10.0, 1.0}, {600.0, 10.0, 0.0}, {900.0, 10.0, 1.0}, {1200.0, 10.0, 0.0}};
    f.set_options("C");
    const auto r = run_wrapper(f, sc, 3u, 200);
    char buf[200];
    std::snprintf(buf, sizeof buf, "%s: covariance finite, symmetric and PSD at every 1 s check (min eig / max %.1e)",
                  label, r.min_eig);
    check(buf, r.cov_ok && r.min_eig > -1e-6);
    // A covariance collapse or a state reset would step the error by the
    // wave amplitude (~0.1 m); onset corrections stay at millimetres.
    std::snprintf(buf, sizeof buf, "%s: Live throughout, heave error never steps by more than 5 mm in one sample (max %.2f mm)",
                  label, 1e3 * r.max_jump);
    check(buf, r.always_live && r.max_jump < 5e-3);
    // Calm windows end before the next transition starts (t_c - w/2).
    const double calm[][2] = {{200, 294}, {750, 894}, {1350, 1500}};
    const double wave[][2] = {{420, 600}, {1020, 1200}};
    for (const auto& c : calm) {
        double e, tr, smax, smin;
        segment(r, c[0], c[1], e, tr, smax, smin);
        // The excitation decays with the tuner's EW memory after waves stop
        // and is back at the numerical floor by the end of the window.
        std::printf("[transitions] %s calm %4.0f-%4.0f s: heave RMS %.2f mm, CARMA sigma <= %.2e, at the end %.2e m/s^2\n",
                    label, c[0], c[1], 1e3 * e, smax, smin);
        std::snprintf(buf, sizeof buf, "%s: calm %.0f-%.0f s heave RMS < 2 mm, sigma < 0.01 and at the 1e-5 floor by the end",
                      label, c[0], c[1]);
        check(buf, e < 2e-3 && smax < 0.01 && smin <= 1.0001e-5);
    }
    for (const auto& w : wave) {
        double e, tr, smax, smin;
        segment(r, w[0], w[1], e, tr, smax, smin);
        std::printf("[transitions] %s waves %4.0f-%4.0f s: heave RMS error %.1f mm of %.1f mm, CARMA sigma >= %.3f m/s^2\n",
                    label, w[0], w[1], 1e3 * e, 1e3 * tr, smin);
        std::snprintf(buf, sizeof buf, "%s: waves %.0f-%.0f s excitation re-entered (sigma > 0.05) and heave error < 25 %% of truth RMS",
                      label, w[0], w[1]);
        check(buf, smin > 0.05 && e < 0.25 * tr);
    }
}

template<class Adapter>
void test_still(Adapter& legacy, Adapter& evidence, const char* label) {
    carma_bench::Scenario sc;
    carma_bench::base_scenario("stationary", sc);
    sc.duration_s = 400.0;
    legacy.set_options("A");
    evidence.set_options("C");
    const auto a = run_wrapper(legacy, sc, 2u, 200);
    const auto c = run_wrapper(evidence, sc, 2u, 200);
    double ea, ta, sa_max, sa_min, ec, tc, sc_max, sc_min;
    segment(a, 200, 400, ea, ta, sa_max, sa_min);
    segment(c, 200, 400, ec, tc, sc_max, sc_min);
    double da[8], dc[8];
    legacy.excitation_diag(da);
    evidence.excitation_diag(dc);
    std::printf("[still] %s: Legacy sigma %.4f f %.3f Hz heave RMS %.1f mm; Evidence sigma %.1e f %.3f Hz heave RMS %.2f mm\n",
                label, sa_min, da[6], 1e3 * ea, sc_max, dc[6], 1e3 * ec);
    char buf[200];
    std::snprintf(buf, sizeof buf, "%s: Legacy still IMU is driven at the 0.05 m/s^2 floor at the 0.05 Hz tuner bound", label);
    check(buf, std::abs(sa_min - 0.05) < 1e-6 && std::abs(sa_max - 0.05) < 1e-6 && std::abs(da[6] - 0.05) < 1e-3);
    std::snprintf(buf, sizeof buf, "%s: Evidence still IMU: sigma at the numerical floor, frequency held at the 0.2 Hz prior", label);
    check(buf, sc_max < 1e-4 && std::abs(dc[6] - 0.2) < 0.01);
    std::snprintf(buf, sizeof buf, "%s: Evidence still-IMU heave RMS below a tenth of Legacy's", label);
    check(buf, ec < 0.1 * ea);
}

// Keel-frame horizontal factors: the inner filter's world-frame covariance is
// the keel-frame one rotated by the filter's own heading (yaw-equivariant),
// and the default stays exactly isotropic.
template <class Adapter>
double keel_frame_case(Adapter& f, const char* label, double heading_deg, const char* opts) {
    carma_bench::Scenario sc;
    carma_bench::base_scenario("jonswap_fp0.29_hs0.5_g3.3", sc);
    sc.duration_s = 300.0;
    sc.att.heading_deg = heading_deg;
    f.set_options(opts);
    (void)run_wrapper(f, sc, 4u, 300);
    const auto& m = f.fusion().raw().mekf();
    const Eigen::Matrix2f S = m.horizontal_sigma();
    const Eigen::Matrix3f R = m.quaternion_boat().toRotationMatrix();
    const double psi = std::atan2(static_cast<double>(R(1, 0)), static_cast<double>(R(0, 0)));
    const double major = 0.5 * std::atan2(2.0 * S(0, 1), static_cast<double>(S(0, 0) - S(1, 1)));
    double dpsi = std::remainder(major - psi, 3.14159265358979323846);
    const double tr = S(0, 0) + S(1, 1), det = S(0, 0) * S(1, 1) - S(0, 1) * S(0, 1);
    const double disc = std::sqrt(std::max(0.0, tr * tr / 4.0 - det));
    const double ratio = (tr / 2.0 - disc) / (tr / 2.0 + disc);
    char buf[220];
    if (std::string(opts).find("sxy") == std::string::npos) {
        std::snprintf(buf, sizeof buf, "%s heading %.0f deg: default factors give an exactly isotropic prior", label, heading_deg);
        check(buf, S(0, 1) == 0.0f && S(0, 0) == S(1, 1));
        return psi;
    }
    std::printf("[keel] %s heading %.0f: filter heading %.2f deg, major axis - filter heading %.4f deg, "
                "eigenvalue ratio %.5f\n", label, heading_deg, psi * 180.0 / 3.14159265358979323846,
                dpsi * 180.0 / 3.14159265358979323846, ratio);
    std::snprintf(buf, sizeof buf, "%s heading %.0f deg: surge axis of the prior on the filter heading, sway/surge variance 0.25",
                  label, heading_deg);
    check(buf, std::abs(dpsi) < 2e-3 && std::abs(ratio - 0.25) < 1e-3);
    return psi;
}

void test_keel_frame_factors() {
    // The filters' world frame is magnetic NED, so their absolute heading
    // carries the simulated declination; turning the hull by 90 degrees must
    // turn the prior with it.
    double psi_q[2], psi_t[2];
    for (int i = 0; i < 2; ++i) {
        const double hdg = i ? 90.0 : 0.0;
        {
            carma3_bench::FullCARMA3<4, false, 4, false> f("A41@C");
            psi_q[i] = keel_frame_case(f, "CARMA3 (4,1)", hdg, "C/sxy1_0.5");
        }
        {
            using M4 = typename carma3_bench::FullTFGCARMA3<4>::M;
            carma3_bench::FullTFGCARMA3<4> f("T41@C", M4::ResetMode::Full, M4::CovarianceWriteMode::Physical);
            psi_t[i] = keel_frame_case(f, "TFG-CARMA3 (4,1)", hdg, "C/sxy1_0.5");
        }
    }
    const double d2r = 3.14159265358979323846 / 180.0;
    const double tq = std::remainder(psi_q[1] - psi_q[0] - 90.0 * d2r, 2.0 * 3.14159265358979323846) / d2r;
    const double tt = std::remainder(psi_t[1] - psi_t[0] - 90.0 * d2r, 2.0 * 3.14159265358979323846) / d2r;
    std::printf("[keel] hull turned 90 deg: prior turned 90 %+.2f deg (CARMA3), 90 %+.2f deg (TFG)\n", tq, tt);
    check("turning the hull by 90 deg turns the keel-frame prior by 90 deg (within 3 deg), both inner filters",
          std::abs(tq) < 3.0 && std::abs(tt) < 3.0);
    carma3_bench::FullCARMA3<4, false, 4, false> f("A41@C");
    keel_frame_case(f, "CARMA3 (4,1)", 90.0, "C");
}

}  // namespace

int main() {
    test_keel_frame_factors();
    test_noise_model();
    test_q_scaling<3>("CARMA(3,0)");
    test_q_scaling<4>("CARMA(4,1)");
    test_q_scaling<5>("CARMA(5,2)");
    test_step_discretization_floor();
    {
        carma3_bench::FullCARMA3<4, false, 4, false> a("A41"), c("A41@C");
        test_still(a, c, "CARMA3 (4,1)");
    }
    {
        carma3_bench::FullCARMA3<4, false, 4, false> f("A41@C");
        test_transitions(f, "CARMA3 (4,1)");
    }
    {
        carma3_bench::FullCARMA3<5, false, 5, false> f("A52@C");
        test_transitions(f, "CARMA3 (5,2)");
    }
    {
        using M4 = typename carma3_bench::FullTFGCARMA3<4>::M;
        carma3_bench::FullTFGCARMA3<4> f("T41@C", M4::ResetMode::Full, M4::CovarianceWriteMode::Physical);
        test_transitions(f, "TFG-CARMA3 (4,1)");
    }
    std::printf("%s: %d failure(s)\n", g_fail ? "FAILED" : "OK", g_fail);
    return g_fail ? 1 : 0;
}
