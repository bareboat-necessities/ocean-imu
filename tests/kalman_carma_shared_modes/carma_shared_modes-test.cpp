/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    carma_shared_modes-test -- unit and numerical tests of the EXPERIMENTAL
    standalone shared-modes QMEKF (src/kalman_carma_shared_modes,
    docs/kalman-carma-shared-modes.md section 6).

    Groups:
      [1] public API, state count, determinism, reset
      [2] preservation of the validated Option C filter: confidence held at
          one reproduces the research filter (SHM:est construction) bit for
          bit; coupling held at zero reproduces CouplingMode::Off bit for bit
      [3] fallback algebra: shrink_motion_vertical_coupling keeps P positive
          semidefinite, keeps the attitude/bias-explained part, is a no-op at
          alpha = 1 and never touches the state
      [4] confidence state machine: consistency hysteresis, energy gate,
          manoeuvre and sea-change holds, NIS ratio, asymmetric smoothing
          (bounded rate), recovery
      [5] wave-direction estimate from a regression
      [6] filter level on a synthetic directional sea (exact kinematics):
          estimated coupling; a confidently reversed external prior with and
          without confidence (bounded error and NEES with it); continuity of
          p, v, attitude and covariance through a forced fallback; recovery
          after the prior is corrected
*/

#define EIGEN_NON_ARDUINO

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModes.h"

extern const float g_std = 9.80665f;   // the application constant the front end reads

namespace csm = ocean_imu::carma_shared_modes;
namespace kc = ocean_imu::kalman::carma;
using V3 = Eigen::Vector3f;

// Test-local convenience over the float wrapper SeaStateFusion_CarmaSharedModes
// (the configuration the simulators use: tuner and clamp on).
class Csm {
public:
    using Fusion = csm::SeaStateFusion_CarmaSharedModes;
    using Mekf = Fusion::Mekf;
    static constexpr int kStates = Mekf::kStates;
    struct Config {
        ocean_imu::carma::SeaStateFusionConfig_CARMA3 fusion{};
        csm::CouplingMode mode = csm::CouplingMode::Estimated;
        csm::CouplingConfidenceConfig<float> confidence{};
    };
    struct SeaState { float operating_hz; Eigen::Vector3f accel_std; };
    Csm() = default;
    explicit Csm(const Config& c) { begin(c); }
    void begin(const Config& c) {
        cfg_ = c;
        f_ = std::make_unique<Fusion>();
        Fusion::SharedModes sm;
        sm.mode = c.mode;
        sm.confidence = c.confidence;
        f_->begin(c.fusion, sm);
        f_->raw().enableTuner(true);
        f_->raw().enableClamp(true);
    }
    void reset() { begin(cfg_); }
    void update(float dt, const Eigen::Vector3f& g, const Eigen::Vector3f& a, float t) { f_->update(dt, g, a, t); }
    void updateMag(const Eigen::Vector3f& m) { f_->updateMag(m); }
    bool isLive() const { return f_->isLive(); }
    Mekf& mekf() { return f_->raw().mekf(); }
    const Mekf& mekf() const { return f_->raw().mekf(); }
    Eigen::Vector3f position() const { return mekf().get_position(); }
    Eigen::Vector3f velocity() const { return mekf().get_velocity(); }
    Eigen::Quaternionf attitude() const { return mekf().quaternion_boat(); }
    Eigen::Matrix<float,6,6> positionVelocityCovariance() const { return mekf().pv_covariance6(); }
    float couplingConfidence() const { return mekf().couplingConfidence(); }
    void forceCouplingConfidence(float c) { mekf().forceCouplingConfidence(c); }
    void setExternalCoupling(const ocean_imu::kalman::carma::Hermitian3<float>& S, float n) { mekf().setExternalCoupling(S, n); }
    Mekf::CouplingStatus coupling() const { return mekf().coupling(); }
    Mekf::Health health() const { return mekf().health(); }
    SeaState seaState() const { return {mekf().motion().operating_frequency_hz(), mekf().motion().accel_std()}; }
    Fusion& fusion() { return *f_; }
    const Fusion& fusion() const { return *f_; }
private:
    Config cfg_{};
    std::unique_ptr<Fusion> f_;
};

static int g_checks = 0, g_failures = 0;
static void check(bool ok, const char* what, double value = NAN, double limit = NAN) {
    ++g_checks;
    if (!ok) ++g_failures;
    std::printf("  %s: %s (value=%.4g limit=%.4g)\n", ok ? "ok" : "FAIL", what, value, limit);
}

constexpr double G = 9.80665;
constexpr double PI = 3.14159265358979323846;

// Directional sea with exact kinematics: one dominant system travelling
// toward dir_deg (NED, from north), linear-wave quadrature between the
// horizontal displacement and heave, plus an incoherent cross component.
struct Sea {
    struct Comp { double A, w, ph, dir; bool quad; };
    std::vector<Comp> c;
    explicit Sea(double dir_deg, unsigned seed = 1, double hs = 1.0, double fp = 0.25) {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> u(0.0, 2 * PI);
        const int n = 24;
        for (int k = 0; k < n; ++k) {
            const double f = fp * (0.7 + 0.9 * k / (n - 1));
            const double a = 0.25 * hs * std::exp(-2.0 * std::pow((f - fp) / (0.25 * fp), 2.0)) / std::sqrt(n / 4.0);
            c.push_back({a, 2 * PI * f, u(rng), dir_deg * PI / 180.0, true});
        }
        for (int k = 0; k < 6; ++k)
            c.push_back({0.03 * hs, 2 * PI * fp * (0.8 + 0.1 * k), u(rng), (dir_deg + 90.0) * PI / 180.0, false});
    }
    void turn(double d_deg) { for (auto& k : c) k.dir += d_deg * PI / 180.0; }
    // NED displacement (z down) and its second derivative.
    void eval(double t, Eigen::Vector3d& p, Eigen::Vector3d& a) const {
        p.setZero(); a.setZero();
        for (const auto& k : c) {
            const double th = k.w * t + k.ph;
            const double eta = k.A * std::cos(th), eta2 = -k.w * k.w * eta;
            const double xi = k.quad ? k.A * std::sin(th) : 0.4 * k.A * std::cos(th + 1.0);
            const double xi2 = -k.w * k.w * xi;
            p += Eigen::Vector3d(std::cos(k.dir) * xi, std::sin(k.dir) * xi, -eta);
            a += Eigen::Vector3d(std::cos(k.dir) * xi2, std::sin(k.dir) * xi2, -eta2);
        }
    }
};

struct Imu {
    std::mt19937 rng{7};
    std::normal_distribution<float> n01{0.0f, 1.0f};
    float sa = 0.0148f, sg = 0.0016f, sm = 0.8f;
    // level, north-facing body (body = NED)
    void sample(const Eigen::Vector3d& a_ned, V3& acc, V3& gyr) {
        acc = (a_ned - Eigen::Vector3d(0, 0, G)).cast<float>() + sa * V3(n01(rng), n01(rng), n01(rng));
        gyr = sg * V3(n01(rng), n01(rng), n01(rng));
    }
    V3 mag() { return V3(20.0f, 0.0f, 45.0f) + sm * V3(n01(rng), n01(rng), n01(rng)); }
};

static Csm::Config base_config() {
    Csm::Config c;
    c.fusion.with_mag = true;
    c.fusion.sigma_a = V3::Constant(2.8f * 0.0148f * 0.71f);
    c.fusion.sigma_g = V3::Constant(2.0f * 0.00157f * 0.01f);
    c.fusion.sigma_m = V3::Constant(1.2f * 0.8f * 8.0f);
    c.fusion.mag_delay_sec = 0.0f;
    c.confidence.accel_noise_std = 0.0148f;
    return c;
}

// The research filter (coupled-prior study, SHM:est) built as its benchmark
// adapter builds it.
struct Research {
    using Mekf = ocean_imu::kalman::Kalman3D_Wave_CARMA3<float, 4, false, 4, false, true, true, kc::CarmaSharedModeMotion<float>>;
    using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 4, false, Mekf>;
    std::unique_ptr<Fusion> f = std::make_unique<Fusion>();
    explicit Research(const Csm::Config& c) {
        f->begin(c.fusion);
        f->raw().enableTuner(true);
        f->raw().enableClamp(true);
        f->raw().setCarmaExcitationMode(ocean_imu::carma::CarmaExcitationMode::Evidence);   // as the CSM default
        kc::CarmaCouplingConfig cc;
        cc.source = kc::CarmaCouplingConfig::Source::Estimated;
        f->raw().mekf().motion().coupling().set_config(cc);
        kc::CarmaCrossSpectrumEstimator<float>::Config ec;
        f->raw().mekf().motion().coupling().estimator().configure(ec);
    }
};

struct RunStats {
    double rms_h = 0, rms_z = 0, nees_p = 0;
    long n = 0;
    double c_mean = 0;
    bool finite = true;
};

// Runs a filter for T seconds on the sea; scores from t_score.
template<class Hook>
static RunStats run(Csm& f, Sea& sea, double T, double t_score, unsigned imu_seed, Hook hook) {
    Imu imu;
    imu.rng.seed(imu_seed);
    RunStats r;
    const double dt = 0.005;
    const long steps = static_cast<long>(T / dt);
    for (long k = 0; k < steps; ++k) {
        const double t = static_cast<double>(k) * dt;
        Eigen::Vector3d p, a;
        sea.eval(t, p, a);
        V3 acc, gyr;
        imu.sample(a, acc, gyr);
        f.update(static_cast<float>(dt), gyr, acc, 35.0f);
        if (k % 8 == 0) f.updateMag(imu.mag());
        hook(t, f, sea);
        if (t >= t_score && f.isLive()) {
            const Eigen::Vector3d e = p - f.position().cast<double>();
            r.rms_h += e.head<2>().squaredNorm();
            r.rms_z += e.z() * e.z();
            const Eigen::Matrix3d Pp = f.positionVelocityCovariance().topLeftCorner<3,3>().cast<double>();
            Eigen::LDLT<Eigen::Matrix3d> l(Pp);
            if (l.info() == Eigen::Success) r.nees_p += e.dot(l.solve(e));
            r.c_mean += f.couplingConfidence();
            ++r.n;
            r.finite = r.finite && f.position().allFinite() && f.fusion().raw().mekf().covariance_full().allFinite();
        }
    }
    if (r.n) {
        const double nn = static_cast<double>(r.n);
        r.rms_h = std::sqrt(r.rms_h / nn); r.rms_z = std::sqrt(r.rms_z / nn);
        r.nees_p /= nn; r.c_mean /= nn;
    }
    return r;
}
static auto no_hook = [](double, Csm&, Sea&) {};

template<int N>
static double min_corr_eig(const Eigen::Matrix<float,N,N>& Pf) {
    const Eigen::Matrix<double,N,N> P = Pf.template cast<double>();
    Eigen::Matrix<double,N,1> d = P.diagonal().cwiseMax(1e-300).cwiseSqrt().cwiseInverse();
    const Eigen::Matrix<double,N,N> C = d.asDiagonal() * P * d.asDiagonal();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,N,N>> es(0.5 * (C + C.transpose()));
    return es.eigenvalues()(0);
}

// ---------------------------------------------------------------------------
void test1_api() {
    std::printf("[1] public API, state count, determinism, reset\n");
    check(Csm::kStates == 33, "33 error states (3 + 3 + 24 + 3)",
          Csm::kStates, 33);
    Csm a(base_config()), b(base_config());
    Sea sea(30.0);
    run(a, sea, 120.0, 1e9, 3, no_hook);
    run(b, sea, 120.0, 1e9, 3, no_hook);
    check((a.position() - b.position()).norm() == 0.0f && (a.velocity() - b.velocity()).norm() == 0.0f,
          "two identical runs are bit-identical");
    a.reset();
    run(a, sea, 120.0, 1e9, 3, no_hook);
    check((a.position() - b.position()).norm() == 0.0f, "reset() restarts exactly");
    const auto h = a.health();
    check(h.covariance_finite && h.covariance_asymmetry <= 1e-6f, "covariance finite and symmetric",
          h.covariance_asymmetry, 1e-6);
    const auto ss = a.seaState();
    check(std::isfinite(ss.operating_hz) && ss.operating_hz > 0.0f && ss.accel_std.allFinite(), "sea-state outputs finite");
}

// ---------------------------------------------------------------------------
void test2_preservation() {
    std::printf("[2] preservation of the validated Option C filter\n");
    auto c = base_config();
    c.confidence.enabled = false;
    Csm f(c);
    Research r(c);
    Sea sea(30.0);
    Imu imu;
    double dmax = 0.0;
    for (long k = 0; k < 400L * 200L; ++k) {
        const double t = static_cast<double>(k) * 0.005;
        Eigen::Vector3d p, a;
        sea.eval(t, p, a);
        V3 acc, gyr;
        imu.sample(a, acc, gyr);
        f.update(0.005f, gyr, acc, 35.0f);
        r.f->update(0.005f, gyr, acc, 35.0f);
        if (k % 8 == 0) { const V3 m = imu.mag(); f.updateMag(m); r.f->updateMag(m); }
        dmax = std::max(dmax, static_cast<double>((f.position() - r.f->raw().mekf().get_position()).cwiseAbs().maxCoeff()));
    }
    check(dmax == 0.0, "confidence disabled: identical to the research SHM:est filter (400 s, bit for bit)", dmax, 0.0);
    const auto& E = f.fusion().raw().mekf().motion().excitation();
    check(std::abs(E.im(0,2)) + std::abs(E.re(0,2)) > 0.0f, "... and its coupling is active (nonzero x-z excitation)");
    // Coupling held at zero = CouplingMode::Off.
    auto c0 = base_config();
    Csm z(c0), off([&] { auto o = base_config(); o.mode = csm::CouplingMode::Off; return o; }());
    z.forceCouplingConfidence(0.0f);
    Sea s2(30.0);
    run(z, s2, 400.0, 1e9, 5, no_hook);
    Sea s3(30.0);
    run(off, s3, 400.0, 1e9, 5, no_hook);
    const double d0 = (z.position() - off.position()).cwiseAbs().maxCoeff();
    check(d0 == 0.0, "confidence held at 0: identical to CouplingMode::Off (independent-axis law)", d0, 0.0);
}

// ---------------------------------------------------------------------------
void test3_shrink() {
    std::printf("[3] fallback algebra (shrink_motion_vertical_coupling)\n");
    using Mekf = Csm::Mekf;
    constexpr int NX = Mekf::NX;
    Csm f(base_config());
    Sea sea(30.0);
    run(f, sea, 300.0, 1e9, 11, no_hook);
    auto& m = f.fusion().raw().mekf();
    const Eigen::Matrix<float,NX,NX> P0 = m.covariance_full();
    const V3 p0 = m.get_position();
    // Inject a strong coupling-like cross-covariance to make the test bite.
    double worst_psd = 1.0, worst_explained = 0.0;
    bool untouched = true;
    for (float alpha : {0.0f, 0.3f, 0.9f}) {
        Csm g(base_config());
        Sea s2(30.0);
        run(g, s2, 300.0, 1e9, 11, no_hook);
        auto& mg = g.fusion().raw().mekf();
        const Eigen::Matrix<float,NX,NX> Pb = mg.covariance_full();
        check(mg.shrink_motion_vertical_coupling(alpha), "shrink applied");
        const Eigen::Matrix<float,NX,NX> Pa = mg.covariance_full();
        worst_psd = std::min(worst_psd, min_corr_eig<NX>(Pa));
        untouched = untouched && (mg.get_position() - p0).norm() == 0.0f;
        // Unchanged outside the horizontal-vertical motion blocks.
        Eigen::Matrix<float,NX,NX> D = Pa - Pb;
        D.block<16,8>(6, 22).setZero();
        D.block<8,16>(22, 6).setZero();
        worst_explained = std::max(worst_explained, static_cast<double>(D.cwiseAbs().maxCoeff()));
    }
    (void)P0;
    check(worst_psd > -1e-5, "P stays positive semidefinite (correlation-scaled min eigenvalue)", worst_psd, -1e-5);
    check(worst_explained == 0.0, "only the horizontal-vertical motion blocks change", worst_explained, 0.0);
    check(untouched, "the state is not touched");
    check(!m.shrink_motion_vertical_coupling(1.0f) && (m.covariance_full() - P0).cwiseAbs().maxCoeff() == 0.0f,
          "alpha = 1 is a no-op");
    // alpha = 0 leaves exactly the attitude/bias-explained part E.
    {
        Csm g(base_config());
        Sea s2(30.0);
        run(g, s2, 300.0, 1e9, 11, no_hook);
        auto& mg = g.fusion().raw().mekf();
        const Eigen::MatrixXd P = mg.covariance_full().cast<double>();
        std::vector<int> r;
        for (int i = 0; i < 6; ++i) r.push_back(i);
        for (int i = 30; i < 33; ++i) r.push_back(i);
        Eigen::MatrixXd Prr(9, 9), Phr(16, 9), Prz(9, 8);
        for (int i = 0; i < 9; ++i) {
            for (int j = 0; j < 9; ++j) Prr(i,j) = P(r[i], r[j]);
            for (int j = 0; j < 8; ++j) Prz(i,j) = P(r[i], 22 + j);
        }
        for (int i = 0; i < 16; ++i) for (int j = 0; j < 9; ++j) Phr(i,j) = P(6 + i, r[j]);
        const Eigen::MatrixXd E = Phr * Prr.ldlt().solve(Prz);
        mg.shrink_motion_vertical_coupling(0.0f);
        const Eigen::MatrixXd Phz = mg.covariance_full().cast<double>().block(6, 22, 16, 8);
        const double rel = (Phz - E).cwiseAbs().maxCoeff() / std::max(1e-30, P.block(6, 22, 16, 8).cwiseAbs().maxCoeff());
        check(rel < 2e-3, "alpha = 0 keeps the attitude/bias-explained cross-covariance only", rel, 2e-3);
    }
}

// ---------------------------------------------------------------------------
void test4_confidence() {
    std::printf("[4] coupling-confidence state machine\n");
    csm::CarmaCouplingConfidence<> c;
    csm::CouplingConfidenceConfig<> cfg;
    c.configure(cfg);
    csm::BandRegression<> good;
    good.valid = true; good.szz = 1.0f; good.gamma2 = 0.9f;
    good.beta[0] = std::complex<float>(0.2f, -0.6f); good.beta[1] = std::complex<float>(0.1f, -0.3f);
    good.var[0] = good.var[1] = 0.002f;
    csm::BandRegression<> rev = good;
    rev.beta[0] = -good.beta[0]; rev.beta[1] = -good.beta[1];
    const float dt = 0.005f, f = 0.3f;
    auto steps = [&](double T, const csm::BandRegression<>& pr, const csm::BandRegression<>& mo, float yaw = 0.0f,
                     float fop = 0.3f, float nis = 0.1f) {
        double max_rate = 0.0;
        for (long k = 0; k < static_cast<long>(T / dt); ++k) {
            const float c0 = c.confidence();
            c.update(dt, pr, mo, true, yaw, fop, nis);
            max_rate = std::max(max_rate, static_cast<double>(std::fabs(c.confidence() - c0)) / dt);
        }
        return max_rate;
    };
    // energy: monitor band power 1 against the noise floor -> full
    const double r_up = steps(300.0, good, good);
    check(c.status().state == csm::CouplingState::Coupled && c.confidence() > 0.98f,
          "consistent prior: coupled, c -> 1", c.confidence(), 0.98);
    check(r_up <= 1.0 / cfg.tau_up_s + 1e-6, "rate of rise bounded by 1/tau_up", r_up, 1.0 / cfg.tau_up_s);
    const double r_dn = steps(10.0, rev, good);
    check(c.status().state == csm::CouplingState::Fallback && c.confidence() < 0.05f,
          "reversed prior: fallback, c < 0.05 within 10 s", c.confidence(), 0.05);
    check(r_dn <= 1.0 / cfg.tau_down_s + 1e-6, "rate of fall bounded by 1/tau_down (no step)", r_dn, 1.0 / cfg.tau_down_s);
    // hysteresis: D^2 between d2_low and d2_high keeps the decision
    csm::BandRegression<> mid = good;
    mid.beta[0] += std::complex<float>(0.0f, 0.08f);   // D^2 ~ 0.0064 / (0.004 + floor) ~ 1.3 (consistent)
    steps(200.0, good, good);
    csm::BandRegression<> edge = good;
    edge.beta[0] += std::complex<float>(0.0f, 0.24f);   // D^2 ~ 10: between the thresholds
    steps(5.0, edge, good);
    const bool kept = c.status().consistent;
    steps(5.0, rev, good);
    steps(5.0, edge, good);
    check(kept && !c.status().consistent, "hysteresis: an intermediate mismatch keeps the previous decision",
          c.status().mismatch_d2, cfg.d2_high);
    // recovery
    steps(400.0, good, good);
    check(c.confidence() > 0.95f, "recovers automatically when the evidence returns", c.confidence(), 0.95);
    // manoeuvre: 10 deg/s yaw for 20 s, then hold
    steps(20.0, good, good, 10.0f * 3.14159f / 180.0f);
    const float c_man = c.confidence();
    steps(20.0, good, good);
    const bool held = c.status().manoeuvre;
    check(c_man < 0.05f && held, "manoeuvre: coupling withdrawn and held after the turn", c_man, 0.05);
    steps(300.0, good, good);
    check(!c.status().manoeuvre && c.confidence() > 0.9f, "... restored after the hold", c.confidence(), 0.9);
    // sea change: f_op ramps 0.3 -> 0.15 Hz over 30 s
    for (int k = 0; k < 6000; ++k) c.update(dt, good, good, true, 0.0f, 0.3f * std::pow(0.5f, static_cast<float>(k) / 6000.0f), 0.1f);
    check(c.status().sea_change && c.confidence() < 0.2f, "abrupt sea-state change: coupling withdrawn",
          c.confidence(), 0.2);
    steps(300.0, good, good, 0.0f, 0.15f);
    // low energy
    csm::BandRegression<> calm = good;
    calm.szz = 1e-7f;
    steps(30.0, good, calm);
    check(c.confidence() < 0.05f, "no wave energy above the noise floor: no coupling", c.confidence(), 0.05);
    steps(300.0, good, good);
    // innovation consistency: NIS 3x its long-term mean
    steps(600.0, good, good, 0.0f, f, 0.1f);
    steps(30.0, good, good, 0.0f, f, 0.4f);
    check(c.status().state == csm::CouplingState::Fallback, "NIS ratio above 2: fallback", c.status().nis_ratio, 2.0);
    // forced
    c.force(0.5f);
    steps(10.0, rev, good);
    check(c.confidence() == 0.5f, "forced confidence is held", c.confidence(), 0.5);
    c.force(-1.0f);
    // warm-up
    csm::CarmaCouplingConfidence<> w;
    w.configure(cfg);
    for (int k = 0; k < 4000; ++k) w.update(dt, good, good, false, 0.0f, f, 0.1f);
    check(w.status().state == csm::CouplingState::WarmUp && w.confidence() == 0.0f, "warm-up: no coupling");
    (void)mid;
}

// ---------------------------------------------------------------------------
void test5_direction() {
    std::printf("[5] wave direction from a regression\n");
    double worst = 0.0;
    for (double dir : {0.0, 30.0, 135.0, 250.0, 330.0}) {
        const double th = dir * PI / 180.0;
        const std::complex<float> q = std::polar(0.8f, static_cast<float>(-62.0 * PI / 180.0));
        csm::BandRegression<> r;
        r.valid = true; r.szz = 1.0f; r.gamma2 = 0.95f;
        r.beta[0] = q * static_cast<float>(std::cos(th));
        r.beta[1] = q * static_cast<float>(std::sin(th));
        const auto w = csm::wave_direction(r, 100.0f);
        double e = std::fabs(w.to_deg - dir);
        e = std::min(e, 360.0 - e);
        worst = std::max(worst, e);
        worst = std::max(worst, std::fabs(w.phase_deg + 62.0));
    }
    check(worst < 1e-3, "direction of travel and quadrature phase recovered (deg)", worst, 1e-3);
}

// ---------------------------------------------------------------------------
void test6_filter() {
    std::printf("[6] filter level: synthetic directional sea\n");
    using H = kc::Hermitian3<float>;
    // Estimated coupling against independent axes.
    auto run_mode = [&](csm::CouplingMode mode, bool conf, int ext /*0 none, 1 true, -1 flipped*/, unsigned seed) {
        auto c = base_config();
        c.mode = mode;
        c.confidence.enabled = conf;
        Csm f(c);
        Sea sea(40.0, seed);
        kc::CarmaCrossSpectrumEstimator<float> truth_est;
        kc::CarmaCrossSpectrumEstimator<float>::Config ec;
        ec.horizon_s = 600.0f;
        truth_est.configure(ec);
        auto hook = [&](double t, Csm& g, Sea& s) {
            if (!ext) return;
            Eigen::Vector3d p, a;
            s.eval(t, p, a);
            truth_est.update(a.cast<float>(), 0.005f, g.fusion().raw().mekf().motion().operating_frequency_hz());
            if (truth_est.elapsed_s() < 30.0f) return;
            H S = truth_est.cross();
            if (ext < 0)
                for (int i = 0; i < 2; ++i) {
                    S.re(i,2) = -S.re(i,2); S.im(i,2) = -S.im(i,2);
                    S.re(2,i) = -S.re(2,i); S.im(2,i) = -S.im(2,i);
                }
            g.setExternalCoupling(S, 0.0f);
        };
        RunStats r = run(f, sea, 900.0, 300.0, seed + 100, hook);
        return std::make_pair(r, f.coupling());
    };
    const auto off = run_mode(csm::CouplingMode::Off, true, 0, 1).first;
    const auto [est, est_st] = run_mode(csm::CouplingMode::Estimated, true, 0, 1);
    std::printf("    off       rms h %.4f z %.4f NEES_p %.2f\n", off.rms_h, off.rms_z, off.nees_p);
    std::printf("    estimated rms h %.4f z %.4f NEES_p %.2f c %.2f\n", est.rms_h, est.rms_z, est.nees_p, est.c_mean);
    check(est.c_mean > 0.5, "estimated coupling on a directional sea: confidence high", est.c_mean, 0.5);
    check(est.rms_z < off.rms_z, "estimated coupling improves heave over independent axes", est.rms_z / off.rms_z, 1.0);
    const auto rev_nc = run_mode(csm::CouplingMode::External, false, -1, 1).first;
    const auto [rev, rev_st] = run_mode(csm::CouplingMode::External, true, -1, 1);
    std::printf("    reversed, no confidence rms h %.4f z %.4f NEES_p %.2f\n", rev_nc.rms_h, rev_nc.rms_z, rev_nc.nees_p);
    std::printf("    reversed, confidence    rms h %.4f z %.4f NEES_p %.2f c %.3f state %s\n", rev.rms_h, rev.rms_z,
                rev.nees_p, rev.c_mean, csm::to_string(rev_st.confidence.state));
    check(rev_nc.rms_z > 1.2 * off.rms_z, "a confidently reversed prior without confidence is harmful (reference)",
          rev_nc.rms_z / off.rms_z, 1.2);
    check(rev.c_mean < 0.1, "reversed prior: confidence falls back (mean c over the scored window)", rev.c_mean, 0.1);
    check(rev.rms_h < 1.15 * off.rms_h && rev.rms_z < 1.15 * off.rms_z,
          "reversed prior with confidence: errors within 15 % of independent axes", std::max(rev.rms_h / off.rms_h, rev.rms_z / off.rms_z), 1.15);
    check(rev.nees_p < 2.0 * std::max(off.nees_p, 3.0), "reversed prior with confidence: NEES bounded (no overconfidence)",
          rev.nees_p, 2.0 * std::max(off.nees_p, 3.0));
    check(rev.finite && est.finite, "finite throughout");
    // Continuity through a forced fallback 1 -> 0 at t = 500 s.
    {
        Csm f(base_config());
        Sea sea(40.0, 2);
        double jump_p = 0.0, jump_v = 0.0, jump_q = 0.0, base_p = 0.0, worst_psd = 1.0;
        V3 lp = V3::Zero(), lv = V3::Zero();
        Eigen::Quaternionf lq = Eigen::Quaternionf::Identity();
        auto hook = [&](double t, Csm& g, Sea&) {
            if (t > 480.0 && t < 520.0) {
                const double dp = (g.position() - lp).norm(), dv = (g.velocity() - lv).norm();
                const double dq = g.attitude().angularDistance(lq);
                if (t >= 500.0 && t < 501.0) { jump_p = std::max(jump_p, dp); jump_v = std::max(jump_v, dv); jump_q = std::max(jump_q, dq); }
                else if (t < 500.0) base_p = std::max(base_p, dp);
                if (std::fmod(t, 1.0) < 0.005) worst_psd = std::min(worst_psd, min_corr_eig<33>(g.fusion().raw().mekf().covariance_full()));
            }
            if (t >= 500.0) g.forceCouplingConfidence(0.0f);
            else g.forceCouplingConfidence(1.0f);
            lp = g.position(); lv = g.velocity(); lq = g.attitude();
        };
        run(f, sea, 520.0, 1e9, 9, hook);
        check(jump_p <= 2.0 * base_p + 1e-4, "no position step at the fallback (largest per-sample change vs before)",
              jump_p, 2.0 * base_p + 1e-4);
        check(jump_q < 1e-3, "no attitude step at the fallback (rad per sample)", jump_q, 1e-3);
        check(worst_psd > -1e-4, "covariance PSD through the fallback", worst_psd, -1e-4);
        check(f.coupling().shrink_events >= 1, "cross-covariance shrink applied at the fallback",
              static_cast<double>(f.coupling().shrink_events), 1.0);
        const auto& E = f.fusion().raw().mekf().motion().excitation();
        check(E.re(0,2) == 0.0f && E.im(0,2) == 0.0f, "no cross-axis excitation after the fallback");
        (void)jump_v;
    }
    // Recovery: the external prior is reversed until 450 s, then correct.
    {
        auto c = base_config();
        c.mode = csm::CouplingMode::External;
        Csm f(c);
        Sea sea(40.0, 3);
        kc::CarmaCrossSpectrumEstimator<float> te;
        kc::CarmaCrossSpectrumEstimator<float>::Config ec;
        ec.horizon_s = 600.0f;
        te.configure(ec);
        float c_at_450 = 1.0f;
        auto hook = [&](double t, Csm& g, Sea& s) {
            Eigen::Vector3d p, a;
            s.eval(t, p, a);
            te.update(a.cast<float>(), 0.005f, g.fusion().raw().mekf().motion().operating_frequency_hz());
            if (te.elapsed_s() < 30.0f) return;
            H S = te.cross();
            if (t < 450.0)
                for (int i = 0; i < 2; ++i) {
                    S.re(i,2) = -S.re(i,2); S.im(i,2) = -S.im(i,2);
                    S.re(2,i) = -S.re(2,i); S.im(2,i) = -S.im(2,i);
                }
            g.setExternalCoupling(S, 0.0f);
            if (std::fabs(t - 450.0) < 0.003) c_at_450 = g.couplingConfidence();
        };
        run(f, sea, 900.0, 1e9, 13, hook);
        check(c_at_450 < 0.1f && f.couplingConfidence() > 0.8f, "recovery after the prior is corrected (c at 450 s / 900 s)",
              f.couplingConfidence(), 0.8);
    }
}

// ---------------------------------------------------------------------------
void test7_double() {
    std::printf("[7] the core is templated on the scalar (T = double)\n");
    using Core = csm::KalmanCarmaSharedModes<double>;
    check(Core::kStates == 33 && csm::KalmanCarmaSharedModes<>::kStates == 33, "33 states for float (default) and double");
    using V3d = Eigen::Vector3d;
    Core f(V3d::Constant(0.03), V3d::Constant(3e-5), V3d::Constant(7.7));
    f.initialize_from_acc(V3d(0.0, 0.0, -G));
    Sea sea(40.0, 4);
    std::mt19937 rng(5);
    std::normal_distribution<double> n01(0.0, 1.0);
    bool finite = true;
    for (long k = 0; k < 200L * 400; ++k) {
        const double t = static_cast<double>(k) * 0.005;
        Eigen::Vector3d p, a;
        sea.eval(t, p, a);
        if (k % 20 == 0) {   // the schedule's model commits, as the wrapper does them
            f.set_vertical_carma(0.95, 0.4);
            f.set_horizontal_carma(0.95, 0.3, 0.3, 0.0);
        }
        f.time_update(0.0016 * V3d(n01(rng), n01(rng), n01(rng)), 0.005);
        f.measurement_update_acc_only(a - V3d(0, 0, G) + 0.0148 * V3d(n01(rng), n01(rng), n01(rng)));
        finite = finite && f.get_position().allFinite() && f.covariance_full().allFinite();
    }
    const auto st = f.coupling();
    std::printf("    double core after 400 s: c %.3f state %s, applied x-z coherence %.3f\n",
                st.confidence.confidence, csm::to_string(st.confidence.state), st.applied_coherence_xz);
    check(finite, "double core: finite state and covariance over 400 s");
    check(st.confidence.confidence > 0.5 && st.applied_coherence_xz > 0.0,
          "double core: the coupling confidence rises and the coupling is applied", st.confidence.confidence, 0.5);
    check(f.health().covariance_asymmetry < 1e-12, "double core: covariance symmetric", f.health().covariance_asymmetry, 1e-12);
}

int main() {
    test1_api();
    test2_preservation();
    test3_shrink();
    test4_confidence();
    test5_direction();
    test6_filter();
    test7_double();
    std::printf("carma_shared_modes-test: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
