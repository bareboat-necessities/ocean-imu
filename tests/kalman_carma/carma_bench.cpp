/*
  Copyright (c) 2026 Mikhail Grushinskiy

  carma_bench -- paired CARMA-versus-OU-III heave benchmark.

  One process runs one (scenario, seed): it generates an independent wave
  truth (CarmaBenchTruth.h, random-phase spectral synthesis with explicit
  RAOs; no CARMA model is involved), corrupts it once with the repository's
  shared IMU and magnetometer noise models, and feeds the identical sample
  stream to every estimator under test.  Each estimator therefore sees the
  same sea realization, sensor noise, biases, sample timing, magnetometer
  cadence, duration and scoring windows.

  --level vertical   Reduced vertical subproblem under the prescribed
                     (true) attitude: y = a_z + b_z + n, with the true
                     body-z bias history projected to the vertical.  Filters:
                       OU3v     the OU-III vertical chain [v,p,S,a_w] + b
                                (shared IntegratedOUChain code, SpectralMSE
                                r_S, tau-scaled S cadence)
                       C30v C41v C52v   CARMA(3,0)/(4,1)/(5,2) + b
                       C52Sv    CARMA(5,2) + S_z with the pole-placement law
                       C52Fv    CARMA(5,2) with the frequency scale fixed at
                                the 0.29 Hz design (amplitude adapts)
                       C52Nv    mechanism ablation: CARMA(5,2) amplitude
                                floored at the band noise floor only
                       C52Mv    oracle-shape diagnostic: CARMA(5,2) with the
                                shape given by --shape52 (e.g. refitted to the
                                true vessel spectrum), kappa = 1
                     --params tuner (default): one shared measurement-only
                     tuner (WavePeriodEstimator + sigma band + auto-tuner,
                     shipping code) drives every filter.  --params oracle:
                     tau, sigma from the realized truth record.
  --level full       Complete INS: shipping SeaStateFusion_OU_III
                     (unchanged) and SeaStateFusion_CARMA variants, behind
                     the same proxy startup and magnetic acquisition.
                       OU3 C30 C41 C52 C52S C52F C52N
  --level truth      No estimator: prints the realized motion levels of the
                     scenario (heave, vertical acceleration, roll, pitch,
                     body rate, incident steepness H_s / L_p).

  Scenario suffixes: _T<s> duration, _norao unit heave RAO, _calbias
  post-calibration accelerometer turn-on bias (+-0.01 m/s^2 instead of the
  shared harness's +-8 mg), _jit +-25 % uniform sample-interval jitter.

  Output: one CSV row per estimator (header with --header).  --timeseries
  writes a decimated trace for plotting.
*/

#define EIGEN_NON_ARDUINO

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "util/W3dSimCommon.h"
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#include "kalman_carma/SeaStateFusionFilter_CARMA.h"
#include "CarmaBenchTruth.h"
#include "CarmaBenchLevel1.h"

using namespace carma_bench;
using Eigen::Vector3f;
namespace ou3 = ocean_imu::ou3;
namespace carma = ocean_imu::carma;

namespace {

// ---------------------------------------------------------------------------
// Scenarios
// ---------------------------------------------------------------------------

// Scenario name without the modifier suffixes handled by the caller.
std::string scenario_base(const std::string& name) {
    std::string base = name;
    for (const char* suffix : {"_calbias", "_jit", "_norao"}) {
        const auto pos = base.find(suffix);
        if (pos != std::string::npos) base.erase(pos, std::strlen(suffix));
    }
    const auto tpos = base.find("_T");
    if (tpos != std::string::npos) base.erase(tpos);
    return base;
}

// n comma-separated numbers (strtod; false on a malformed list).
bool parse_numbers(const std::string& s, double* out, int n) {
    const char* p = s.c_str();
    for (int i = 0; i < n; ++i) {
        char* end = nullptr;
        out[i] = std::strtod(p, &end);
        if (end == p) return false;
        p = end;
        if (i + 1 < n) {
            if (*p != ',') return false;
            ++p;
        }
    }
    return true;
}

Scenario make_scenario(const std::string& full_name) {
    Scenario sc;
    sc.name = full_name;
    const std::string name = scenario_base(full_name);
    if (!base_scenario(name, sc)) throw std::invalid_argument("unknown scenario " + name);
    if (full_name.find("_norao") != std::string::npos) sc.heave_rao = false;
    {
        const auto pos = full_name.find("_T");
        sc.duration_s = (pos == std::string::npos) ? 1200.0 : std::strtod(full_name.c_str() + pos + 2, nullptr);
    }
    return sc;
}

// ---------------------------------------------------------------------------
// Common estimator interface
// ---------------------------------------------------------------------------

struct Sensors {
    float dt;
    Vector3f gyr_ned, acc_ned;
    bool mag_tick;
    Vector3f mag_ned;
    double y_vertical;   // level 1: vertical specific force + g (NED down)
};

struct Estimate {
    bool live = false;
    double h_up = 0, v_up = 0, a_up = 0;
    double var_p = NAN, var_v = NAN, cov_pv = NAN;   // NED [p_z, v_z]
    double nis = NAN;
    int nis_dof = 3;
};

class IEstimator {
public:
    virtual ~IEstimator() = default;
    virtual std::string name() const = 0;
    virtual void step(const Sensors& s) = 0;
    virtual Estimate estimate() const = 0;
    virtual size_t memory_bytes() const = 0;
    virtual int states() const = 0;
};

// Shipping simulator configuration (tests/kalman_ou_iii/kalman_ou_iii-sim.cpp
// FusionAdapter_OU_III), without its dataset-specific direction-RAO and
// low-wave noise weighting, which name the pinned 28 ft hull: the firmware
// sketch sets neither.
constexpr float SIGMA_A_RESCALE = 0.71f;
constexpr float SIGMA_G_RESCALE = 0.01f;
constexpr float SIGMA_M_RESCALE = 8.0f;

const float ACC_SIGMA = 1.51e-3f * 9.80665f;
const float GYR_SIGMA = 0.00157f;
const float MAG_SIGMA_UT = 0.80f;

template<class Cfg>
void fill_config(Cfg& cfg) {
    cfg.with_mag = true;
    cfg.sigma_a = Vector3f::Constant(2.8f * ACC_SIGMA * SIGMA_A_RESCALE);
    cfg.sigma_g = Vector3f::Constant(2.0f * GYR_SIGMA * SIGMA_G_RESCALE);
    cfg.sigma_m = Vector3f::Constant(1.2f * MAG_SIGMA_UT * SIGMA_M_RESCALE);
    cfg.mag_delay_sec = MAG_DELAY_SEC;
}

class FullOU3 final : public IEstimator {
public:
    FullOU3() {
        SeaStateFusionConfig_OU_III cfg;
        fill_config(cfg);
        fusion_.begin(cfg);
        auto& f = fusion_.raw();
        f.setPeriodicAwCovarianceSync(true);
        f.setAwCovarianceSyncCongruent(false);
        f.enableTuner(true);
        f.enableClamp(true);
    }
    std::string name() const override { return "OU3"; }
    void step(const Sensors& s) override {
        fusion_.update(s.dt, s.gyr_ned, s.acc_ned, 35.0f);
        if (s.mag_tick) fusion_.updateMag(s.mag_ned);
    }
    Estimate estimate() const override {
        Estimate e;
        const auto& m = fusion_.raw().mekf();
        e.live = fusion_.isLive();
        e.h_up = -m.get_position().z();
        e.v_up = -m.get_velocity().z();
        e.a_up = -m.get_world_accel().z();
        const auto& P = m.covariance_full();
        constexpr int IV = 8, IP = 11;   // v_z, p_z in [dtheta, b_g, v, p, S, a_w, b_a]
        e.var_p = P(IP, IP); e.var_v = P(IV, IV); e.cov_pv = P(IP, IV);
        e.nis = m.lastAccDiag().nis;
        return e;
    }
    size_t memory_bytes() const override {
        return sizeof(fusion_) + sizeof(fusion_.raw().mekf());
    }
    int states() const override { return 21; }

private:
    SeaStateFusion_OU_III<TrackerType::KALMANF> fusion_;
};

template<int N, bool VI>
class FullCARMA final : public IEstimator {
public:
    FullCARMA(std::string label, bool fixed_shape, float sigma_floor = -1.0f) : label_(std::move(label)) {
        carma::SeaStateFusionConfig_CARMA cfg;
        fill_config(cfg);
        fusion_.begin(cfg);
        auto& f = fusion_.raw();
        f.setPeriodicAwCovarianceSync(true);
        f.enableTuner(true);
        f.enableClamp(true);
        if (fixed_shape) f.setFixedVerticalCarma(1.0f, 0.0f);
        if (sigma_floor >= 0.0f) f.setCarmaSigmaFloor(sigma_floor);
    }
    std::string name() const override { return label_; }
    void step(const Sensors& s) override {
        fusion_.update(s.dt, s.gyr_ned, s.acc_ned, 35.0f);
        if (s.mag_tick) fusion_.updateMag(s.mag_ned);
    }
    Estimate estimate() const override {
        Estimate e;
        const auto& m = fusion_.raw().mekf();
        e.live = fusion_.isLive();
        e.h_up = -m.get_position().z();
        e.v_up = -m.get_velocity().z();
        e.a_up = -m.get_world_accel().z();
        const auto Ppv = m.vertical_pv_covariance();
        e.var_p = Ppv(0,0); e.var_v = Ppv(1,1); e.cov_pv = Ppv(0,1);
        e.nis = m.lastAccDiag().nis;
        return e;
    }
    size_t memory_bytes() const override {
        return sizeof(fusion_) + sizeof(fusion_.raw().mekf());
    }
    int states() const override { return Fusion::Impl::Mekf::NX; }

private:
    using Fusion = carma::SeaStateFusion_CARMA<TrackerType::KALMANF, N, VI>;
    std::string label_;
    Fusion fusion_;
};

// ---------------------------------------------------------------------------
// Level 1: shared measurement-only operating point and reduced filters
// ---------------------------------------------------------------------------

class IVertical : public IEstimator {
public:
    virtual void set_operating_point(const OperatingPoint& op) = 0;
    int nis_dof() const { return 1; }
};

// OU-III vertical chain [v, p, S, a_w] + b.
class VerticalOU3 final : public IVertical {
public:
    VerticalOU3(float R_acc, BiasPrior bp) : R_(R_acc), bp_(bp) {}
    std::string name() const override { return "OU3v"; }
    void set_operating_point(const OperatingPoint& op) override {
        if (!started_) {
            started_ = true;
            kf_.P.setZero();
            kf_.P(0,0) = 1.0f; kf_.P(1,1) = 400.0f; kf_.P(2,2) = 2500.0f;
            kf_.P(3,3) = op.sigma_aw * op.sigma_aw; kf_.P(4,4) = bp_.sigma0 * bp_.sigma0;
        }
        op_ = op;
    }
    void step(const Sensors& s) override {
        if (!started_) return;
        const float dt = s.dt;
        const float tau = op_.tau;
        const float sig = std::max(op_.sigma_floor, op_.sigma_aw);
        Eigen::Matrix<float,4,4> Phi, Qd;
        ocean_imu::kalman::ou_detail::IntegratedOUChain<float,3>::transition(tau, dt, Phi);
        ocean_imu::kalman::ou_detail::IntegratedOUChain<float,3>::process_covariance(tau, dt, sig * sig, Qd);
        Eigen::Matrix<float,5,5> F = Eigen::Matrix<float,5,5>::Zero(), Q = F;
        F.topLeftCorner<4,4>() = Phi;
        Q.topLeftCorner<4,4>() = Qd;
        const float phi_b = std::exp(-dt / bp_.tau);
        F(4,4) = phi_b;
        Q(4,4) = bp_.q * (-0.5f * bp_.tau * std::expm1(-2.0f * dt / bp_.tau));
        kf_.predict(F, Q);
        if (ocean_imu::kalman::ou_detail::periodic_update_due(dt, op_.ts, s_elapsed_)) {
            const int idx[1] = {2};
            const float c[1] = {1.0f};
            const float nis = kf_.last_nis;
            kf_.update(idx, c, 0.0f, op_.rs * op_.rs);
            kf_.last_nis = nis;
        }
        const int idx[2] = {3, 4};
        const float c[2] = {1.0f, 1.0f};
        kf_.update(idx, c, static_cast<float>(s.y_vertical), R_);
    }
    Estimate estimate() const override {
        Estimate e;
        e.live = started_;
        e.h_up = -kf_.x(1); e.v_up = -kf_.x(0); e.a_up = -kf_.x(3);
        e.var_p = kf_.P(1,1); e.var_v = kf_.P(0,0); e.cov_pv = kf_.P(1,0);
        e.nis = kf_.last_nis; e.nis_dof = 1;
        return e;
    }
    size_t memory_bytes() const override { return sizeof(*this); }
    int states() const override { return 5; }

private:
    float R_;
    BiasPrior bp_;
    OperatingPoint op_;
    bool started_ = false;
    float s_elapsed_ = 0.0f;
    LinearKF<5> kf_;
};

// CARMA(N) [+ S_z] + b.
template<int N, bool VI>
class VerticalCARMA final : public IVertical {
    static constexpr int NZ = N + (VI ? 1 : 0);
    static constexpr int NS = NZ + 1;
public:
    VerticalCARMA(std::string label, float R_acc, BiasPrior bp, bool fixed_shape, bool noise_floor_only = false)
        : label_(std::move(label)), R_(R_acc), bp_(bp), fixed_shape_(fixed_shape),
          noise_floor_only_(noise_floor_only) {}
    // Oracle-shape diagnostic: replace the default shape; with fixed_shape the
    // frequency scale stays at kappa = 1, i.e. the shape is used as given.
    void set_shape(const ocean_imu::kalman::carma::CarmaShape& shape) { model_.set_shape(shape); }
    std::string name() const override { return label_; }
    void set_operating_point(const OperatingPoint& op) override {
        op_ = op;
        const double fz = model_.base_zero_crossing_hz();
        const float kappa = fixed_shape_ ? 1.0f : static_cast<float>((0.5 / op.tau) / fz);
        const float floor = noise_floor_only_ ? std::max(op.noise_floor, 1e-4f) : op.sigma_floor;
        const float sa = std::max(floor, op.sigma_aw / 0.9f);
        model_.configure(kappa, sa);
        if (!started_) {
            started_ = true;
            Eigen::Matrix<float,N,N> P0;
            model_.stationary_covariance(P0);
            kf_.P.setZero();
            kf_.P.template topLeftCorner<N,N>() = P0;
            if constexpr (VI) kf_.P(N,N) = 2500.0f;
            kf_.P(NS - 1, NS - 1) = bp_.sigma0 * bp_.sigma0;
        }
        if constexpr (VI) {
            const double c0 = static_cast<double>(model_.noise_gain2()) *
                ocean_imu::kalman::carma::carma_unit_psd<N>(model_.shape(), model_.kappa(), 0.0);
            rs_z_ = static_cast<float>(std::sqrt(c0)) * op.tau / (ou3::R_S_POLE_KAPPA_DEFAULT * std::sqrt(op.ts));
        }
    }
    void step(const Sensors& s) override {
        if (!started_) return;
        const float dt = s.dt;
        Eigen::Matrix<float,NZ,NZ> Fz, Qz;
        if constexpr (VI) {
            Eigen::Matrix<float,NZ,NZ> A = Eigen::Matrix<float,NZ,NZ>::Zero(), W = A;
            A.template topLeftCorner<N,N>() = model_.A();
            A(N, 0) = 1.0f;
            W.template topLeftCorner<N,N>() = model_.W();
            ocean_imu::kalman::carma::discretize_lti<float,NZ>(A, W, dt, Fz, Qz);
        } else {
            model_.discretize(dt, Fz, Qz);
        }
        Eigen::Matrix<float,NS,NS> F = Eigen::Matrix<float,NS,NS>::Zero(), Q = F;
        F.template topLeftCorner<NZ,NZ>() = Fz;
        Q.template topLeftCorner<NZ,NZ>() = Qz;
        F(NS - 1, NS - 1) = std::exp(-dt / bp_.tau);
        Q(NS - 1, NS - 1) = bp_.q * (-0.5f * bp_.tau * std::expm1(-2.0f * dt / bp_.tau));
        kf_.predict(F, Q);
        if constexpr (VI) {
            if (ocean_imu::kalman::ou_detail::periodic_update_due(dt, op_.ts, s_elapsed_) &&
                std::isfinite(rs_z_) && rs_z_ > 0.0f) {
                const int idx[1] = {N};
                const float c[1] = {1.0f};
                const float nis = kf_.last_nis;
                kf_.update(idx, c, 0.0f, rs_z_ * rs_z_);
                kf_.last_nis = nis;
            }
        }
        const int idx[2] = {2, NS - 1};
        const float c[2] = {1.0f, 1.0f};
        kf_.update(idx, c, static_cast<float>(s.y_vertical), R_);
    }
    Estimate estimate() const override {
        Estimate e;
        e.live = started_;
        // NED states; h = p_z.
        e.h_up = -kf_.x(0); e.v_up = -kf_.x(1); e.a_up = -kf_.x(2);
        e.var_p = kf_.P(0,0); e.var_v = kf_.P(1,1); e.cov_pv = kf_.P(0,1);
        e.nis = kf_.last_nis; e.nis_dof = 1;
        return e;
    }
    size_t memory_bytes() const override { return sizeof(*this); }
    int states() const override { return NS; }

private:
    std::string label_;
    float R_;
    BiasPrior bp_;
    bool fixed_shape_;
    bool noise_floor_only_;
    OperatingPoint op_;
    bool started_ = false;
    float s_elapsed_ = 0.0f;
    float rs_z_ = NAN;
    ocean_imu::kalman::carma::CarmaVerticalModel<float,N> model_;
    LinearKF<NS> kf_;
};

// ---------------------------------------------------------------------------
// Metrics
// ---------------------------------------------------------------------------

struct Series {
    std::vector<float> t, eh, ev, ea, nees2, nees1, nis, h_true, h_est;
    std::vector<unsigned char> live;
    double runtime_ns = 0.0;
    long runtime_count = 0;
    bool failed = false;
    std::string fail_reason;
};

struct WindowStats {
    double rms_h = NAN, mean_h = NAN, peak_h = NAN, rms_v = NAN, rms_a = NAN;
    double lf_rms_h = NAN;      // RMS of the 60 s moving average of the heave error
    double trend_m_per_hr = NAN;
    double nees2 = NAN, nees1 = NAN, nis = NAN, nis_in95 = NAN;
    double truth_std_h = NAN;
};

WindowStats window_stats(const Series& s, const std::vector<double>& h_true_all, double t0, double t1,
                         double dt, int nis_dof) {
    WindowStats w;
    double se = 0, sm = 0, peak = 0, sv = 0, sa = 0, sn2 = 0, sn1 = 0, snis = 0, nin = 0, sht = 0, sht2 = 0;
    long n = 0, nn = 0, nnis = 0;
    double st = 0, stt = 0, ste = 0;
    const size_t N = s.t.size();
    for (size_t i = 0; i < N; ++i) {
        if (s.t[i] < t0 || s.t[i] > t1) continue;
        const double e = s.eh[i];
        se += e * e; sm += e; peak = std::max(peak, std::abs(e));
        sv += double(s.ev[i]) * s.ev[i]; sa += double(s.ea[i]) * s.ea[i];
        st += s.t[i]; stt += double(s.t[i]) * s.t[i]; ste += s.t[i] * e;
        sht += h_true_all[i]; sht2 += h_true_all[i] * h_true_all[i];
        ++n;
        if (std::isfinite(s.nees2[i])) { sn2 += s.nees2[i]; sn1 += s.nees1[i]; ++nn; }
        if (std::isfinite(s.nis[i])) {
            snis += s.nis[i]; ++nnis;
            const double lo = nis_dof == 1 ? 0.000982 : 0.2158, hi = nis_dof == 1 ? 5.0239 : 9.3484;
            if (s.nis[i] >= lo && s.nis[i] <= hi) nin += 1.0;
        }
    }
    if (n == 0) return w;
    const double dn = static_cast<double>(n);
    w.rms_h = std::sqrt(se / dn); w.mean_h = sm / dn; w.peak_h = peak;
    w.rms_v = std::sqrt(sv / dn); w.rms_a = std::sqrt(sa / dn);
    const double denom = dn * stt - st * st;
    if (denom > 0) w.trend_m_per_hr = (dn * ste - st * sm) / denom * 3600.0;
    if (nn > 0) { w.nees2 = sn2 / static_cast<double>(nn); w.nees1 = sn1 / static_cast<double>(nn); }
    if (nnis > 0) { w.nis = snis / static_cast<double>(nnis); w.nis_in95 = nin / static_cast<double>(nnis); }
    w.truth_std_h = std::sqrt(std::max(0.0, sht2 / dn - (sht / dn) * (sht / dn)));
    // 60 s moving average of the error.
    const int L = static_cast<int>(std::lround(60.0 / dt));
    double acc = 0, slf = 0;
    long nlf = 0;
    std::vector<size_t> idx;
    for (size_t i = 0; i < N; ++i) if (s.t[i] >= t0 && s.t[i] <= t1) idx.push_back(i);
    for (size_t j = 0; j < idx.size(); ++j) {
        acc += s.eh[idx[j]];
        if (j >= static_cast<size_t>(L)) acc -= s.eh[idx[j - L]];
        if (j + 1 >= static_cast<size_t>(L)) {
            const double m = acc / L;
            slf += m * m;
            ++nlf;
        }
    }
    if (nlf > 0) w.lf_rms_h = std::sqrt(slf / static_cast<double>(nlf));
    return w;
}

double convergence_time(const Series& s, double t_live, double rms_ref, double dt) {
    // First time after which the 30 s rolling RMS heave error stays within
    // 1.5 x the scored-window RMS for the rest of the record.  The window
    // starts at Live: the error is NaN before it, and a NaN in the running
    // sum would make every later window NaN, which std::max(0.0, NaN) used to
    // turn into a rolling RMS of 0 and a convergence time of 0 for every
    // tuner-driven run.
    if (!std::isfinite(t_live) || !(rms_ref > 0)) return NAN;
    const int L = static_cast<int>(std::lround(30.0 / dt));
    const size_t N = s.t.size();
    std::vector<double> roll(N, NAN);
    double acc = 0;
    size_t n = 0;
    for (size_t i = 0; i < N; ++i) {
        if (!std::isfinite(s.eh[i])) { acc = 0; n = 0; continue; }
        acc += double(s.eh[i]) * s.eh[i];
        ++n;
        if (n > static_cast<size_t>(L)) acc -= double(s.eh[i - L]) * s.eh[i - L];
        if (n >= static_cast<size_t>(L)) roll[i] = std::sqrt(std::max(0.0, acc / L));
    }
    double t_conv = NAN;
    for (size_t i = N; i-- > 0;) {
        if (s.t[i] < t_live) break;
        if (!(roll[i] <= 1.5 * rms_ref)) break;
        t_conv = s.t[i] - 30.0;
    }
    return std::isfinite(t_conv) ? std::max(0.0, t_conv - t_live) : NAN;
}

// ---------------------------------------------------------------------------

struct Args {
    std::string scenario = "jonswap_fp0.29_hs0.8_g3.3";
    unsigned seed = 1;
    std::string level = "vertical";
    std::string params = "tuner";
    std::string filters;
    std::string shape52;
    std::string timeseries;
    int decimate = 10;
    bool header = false;
    double score_start = 300.0;
};

Args parse(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        const std::string k = argv[i];
        auto val = [&]() -> std::string {
            if (++i >= argc) throw std::invalid_argument("missing value for " + k);
            return argv[i];
        };
        if (k == "--scenario") a.scenario = val();
        else if (k == "--seed") a.seed = static_cast<unsigned>(std::stoul(val()));
        else if (k == "--level") a.level = val();
        else if (k == "--params") a.params = val();
        else if (k == "--filters") a.filters = val();
        else if (k == "--shape52") a.shape52 = val();
        else if (k == "--timeseries") a.timeseries = val();
        else if (k == "--decimate") a.decimate = std::stoi(val());
        else if (k == "--score-start") a.score_start = std::stod(val());
        else if (k == "--header") a.header = true;
        else throw std::invalid_argument("unknown argument " + k);
    }
    return a;
}

std::vector<std::string> split(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) if (!item.empty()) out.push_back(item);
    return out;
}

const char* HEADER =
    "level,params,scenario,seed,filter,states,memory_bytes,runtime_ns,hs_incident,hs_heave,"
    "t_live,conv_s,failed,"
    "rms_h,rms_h_pct_hs,mean_h,peak_h,rms_v,rms_a,lf_rms_h,trend_m_per_hr,nees2,nees1,nis,nis_in95,"
    "f60_rms_h,f60_mean_h,f60_peak_h,f60_rms_v,f60_rms_a,f60_nees2,f60_nis";

// Motion levels of one realization, for checking that a scenario is
// physically plausible.  Steepness uses the deep-water peak wavelength of the
// steepest spectral part.
int truth_stats(const Scenario& sc, const Args& args, double dt, long steps) {
    WaveTruth truth(sc, w3d_expand_seed(args.seed, 21u), dt);
    double sh = 0, sa = 0, sr = 0, sp = 0, amax = 0, rmax = 0, pmax = 0, wmax = 0, fmax = 0;
    for (long k = 0; k < steps; ++k) {
        const TruthSample ts = truth.next();
        const double r = (ts.roll - sc.att.heel_deg * PI / 180.0) * 180.0 / PI;
        const double q = (ts.pitch - sc.att.trim_deg * PI / 180.0) * 180.0 / PI;
        sh += ts.h_up * ts.h_up; sa += ts.a_up * ts.a_up; sr += r * r; sp += q * q;
        amax = std::max(amax, std::abs(ts.a_up));
        rmax = std::max(rmax, std::abs(r)); pmax = std::max(pmax, std::abs(q));
        wmax = std::max(wmax, ts.omega_b.norm());
        fmax = std::max(fmax, ts.f_b.norm());
    }
    double steep = 0.0;
    for (const auto* parts : {&sc.parts, &sc.parts_b})
        for (const auto& pt : *parts) {
            const double fp = pt.type == SpectrumType::Jonswap ? pt.fp_hz : pt.f_lo_hz;
            steep = std::max(steep, pt.hs_m / (G / (TWO_PI * fp * fp)));
        }
    const double n = static_cast<double>(steps);
    std::printf("scenario=%s seed=%u hs_heave_m=%.3f a_rms=%.3f a_max=%.3f roll_rms_deg=%.2f roll_max_deg=%.2f "
                "pitch_rms_deg=%.2f pitch_max_deg=%.2f rate_max_rad_s=%.3f f_max_mps2=%.2f steepness=%.4f\n",
                args.scenario.c_str(), args.seed, 4.0 * std::sqrt(sh / n), std::sqrt(sa / n), amax,
                std::sqrt(sr / n), rmax, std::sqrt(sp / n), pmax, wmax, fmax, steep);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    Args args;
    try {
        args = parse(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    if (args.header) {
        std::cout << HEADER << "\n";
        if (argc == 2) return 0;
    }

    const Scenario sc = make_scenario(args.scenario);
    const double dt = 1.0 / 200.0;
    const long steps = static_cast<long>(std::lround(sc.duration_s / dt));
    const bool full = args.level == "full";
    if (args.level == "truth") return truth_stats(sc, args, dt, steps);
    if (!full && args.level != "vertical") {
        std::cerr << "ERROR: --level must be vertical, full or truth\n";
        return 2;
    }

    // Sensor models: the shared harness parameters, independent streams.
    const unsigned base = args.seed;
    const bool calbias = args.scenario.find("_calbias") != std::string::npos;
    const bool jitter = args.scenario.find("_jit") != std::string::npos;
    ImuNoiseModel acc_noise = make_imu_noise_model(ACC_SIGMA, calbias ? 0.01f : 8e-3f * 9.80665f, 0.0005f,
                                                   w3d_expand_seed(base, 11u), w3d_expand_seed(base, 12u));
    ImuNoiseModel gyr_noise = make_imu_noise_model(GYR_SIGMA, 0.05f * 3.14159265f / 180.0f, 0.00001f,
                                                   w3d_expand_seed(base, 13u), w3d_expand_seed(base, 14u));
    MagNoiseModel mag_noise = make_mag_noise_model(MAG_SIGMA_UT, 2.0f, 0.01f, 0.015f, 0.010f, 1.0f,
                                                   w3d_expand_seed(base, 15u), w3d_expand_seed(base, 16u));
    const Vector3f mag_world_ned = MagSim_WMM::mag_world_aero();
    const float mag_dt = 1.0f / 25.0f;
    float mag_phase = 0.0f;

    WaveTruth truth(sc, w3d_expand_seed(base, 21u), dt);
    const double hs_inc = sc.stationary ? 0.0
        : (sc.fade.enabled ? WaveTruth::incident_hs(sc.parts_b) : WaveTruth::incident_hs(sc.parts));

    // Estimators.
    std::vector<std::unique_ptr<IEstimator>> est;
    std::vector<IVertical*> vert;
    const float R_acc = std::pow(2.8f * ACC_SIGMA * SIGMA_A_RESCALE, 2.0f);
    const BiasPrior bp{};
    std::vector<std::string> names = split(args.filters);
    if (names.empty()) {
        names = full ? std::vector<std::string>{"OU3", "C30", "C41", "C52", "C52S", "C52F", "C52N"}
                     : std::vector<std::string>{"OU3v", "C30v", "C41v", "C52v", "C52Sv", "C52Fv", "C52Nv"};
    }
    for (const auto& n : names) {
        if (full) {
            if (n == "OU3") est.push_back(std::make_unique<FullOU3>());
            else if (n == "C30") est.push_back(std::make_unique<FullCARMA<3,false>>(n, false));
            else if (n == "C41") est.push_back(std::make_unique<FullCARMA<4,false>>(n, false));
            else if (n == "C52") est.push_back(std::make_unique<FullCARMA<5,false>>(n, false));
            else if (n == "C52S") est.push_back(std::make_unique<FullCARMA<5,true>>(n, false));
            else if (n == "C52F") est.push_back(std::make_unique<FullCARMA<5,false>>(n, true));
            else if (n == "C52N") est.push_back(std::make_unique<FullCARMA<5,false>>(n, false, 0.0f));
            else { std::cerr << "ERROR: unknown full filter " << n << "\n"; return 2; }
        } else {
            std::unique_ptr<IVertical> v;
            if (n == "OU3v") v = std::make_unique<VerticalOU3>(R_acc, bp);
            else if (n == "C30v") v = std::make_unique<VerticalCARMA<3,false>>(n, R_acc, bp, false);
            else if (n == "C41v") v = std::make_unique<VerticalCARMA<4,false>>(n, R_acc, bp, false);
            else if (n == "C52v") v = std::make_unique<VerticalCARMA<5,false>>(n, R_acc, bp, false);
            else if (n == "C52Sv") v = std::make_unique<VerticalCARMA<5,true>>(n, R_acc, bp, false);
            else if (n == "C52Fv") v = std::make_unique<VerticalCARMA<5,false>>(n, R_acc, bp, true);
            else if (n == "C52Nv") v = std::make_unique<VerticalCARMA<5,false>>(n, R_acc, bp, false, true);
            else if (n == "C52Mv") {
                double q[7];
                if (!parse_numbers(args.shape52, q, 7)) {
                    std::cerr << "ERROR: C52Mv needs --shape52 f1,z1,f2,z2,fl,fz,zz\n";
                    return 2;
                }
                auto m = std::make_unique<VerticalCARMA<5,false>>(n, R_acc, bp, true);
                m->set_shape(ocean_imu::kalman::carma::CarmaShape{q[0], q[1], q[2], q[3], q[4], q[5], q[6]});
                v = std::move(m);
            }
            else { std::cerr << "ERROR: unknown vertical filter " << n << "\n"; return 2; }
            vert.push_back(v.get());
            est.push_back(std::move(v));
        }
    }

    // Oracle operating point from the realized truth (vertical level only):
    // T_z = 2 pi std(h)/std(v), sigma_a,B = std(a).  A pre-pass over an
    // identical generator (same seed) leaves the main pass untouched.
    OperatingPoint oracle;
    if (!full && args.params == "oracle") {
        WaveTruth pre(sc, w3d_expand_seed(base, 21u), dt);
        double sh = 0, sv = 0, sa = 0;
        for (long k = 0; k < steps; ++k) {
            const TruthSample ts = pre.next();
            sh += ts.h_up * ts.h_up; sv += ts.v_up * ts.v_up; sa += ts.a_up * ts.a_up;
        }
        const double tz = TWO_PI * std::sqrt(sh / std::max(sv, 1e-30));
        oracle.valid = true;
        oracle.tau = static_cast<float>(std::min(std::max(0.5 * tz, double(ou3::MIN_TAU_S)), double(ou3::MAX_TAU_S)));
        oracle.sigma_aw = static_cast<float>(0.9 * std::sqrt(sa / static_cast<double>(steps)));
        oracle.rs = std::min(std::max(VerticalTuner::rs_target_(oracle.tau, oracle.sigma_aw), ou3::MIN_R_S), ou3::MAX_R_S);
        oracle.ts = std::min(std::max(seastate::common::defaults::PSEUDO_UPDATE_TAU_RATIO * oracle.tau,
                                      seastate::common::defaults::PSEUDO_UPDATE_PERIOD_MIN_S),
                             ou3::PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT);
        oracle.sigma_floor = 0.05f;
        oracle.noise_floor = 0.0f;
    } else if (!full && args.params != "tuner") {
        std::cerr << "ERROR: --params must be tuner or oracle\n";
        return 2;
    }

    VerticalTuner tuner;
    std::vector<Series> series(est.size());
    for (auto& s : series) {
        s.t.reserve(steps); s.eh.reserve(steps); s.ev.reserve(steps); s.ea.reserve(steps);
        s.nees2.reserve(steps); s.nees1.reserve(steps); s.nis.reserve(steps); s.live.reserve(steps);
    }
    std::vector<double> h_true_all;
    h_true_all.reserve(steps);
    std::vector<double> t_live(est.size(), NAN);

    std::ofstream ts_out;
    if (!args.timeseries.empty()) {
        ts_out.open(args.timeseries);
        ts_out << "t,eta_up,h_true,v_true,a_true";
        for (const auto& e : est) ts_out << ",h_" << e->name() << ",v_" << e->name() << ",a_" << e->name();
        ts_out << "\n";
    }

    std::mt19937 jitter_rng(w3d_expand_seed(base, 31u));
    std::uniform_real_distribution<double> jitter_u(-0.25, 0.25);
    double t_now = 0.0;
    for (long k = 0; k < steps; ++k) {
        // Non-uniform timing: every consumer receives the same jittered step.
        const double dt_k = (jitter && k > 0) ? dt * (1.0 + jitter_u(jitter_rng)) : dt;
        if (k > 0) t_now += dt_k;
        const TruthSample ts = jitter ? truth.at(t_now) : truth.next();
        // Sensor corruption, once, shared by every estimator.
        const Vector3f f_true = ts.f_b.cast<float>();
        const Vector3f w_true = ts.omega_b.cast<float>();
        Vector3f acc = apply_imu_noise(f_true, acc_noise, static_cast<float>(dt_k));
        acc += ts.extra_acc_bias_b.cast<float>();
        const Vector3f gyr = apply_imu_noise(w_true, gyr_noise, static_cast<float>(dt_k));
        Sensors s;
        s.dt = static_cast<float>(dt_k);
        s.acc_ned = acc;
        s.gyr_ned = gyr;
        mag_phase += s.dt;
        s.mag_tick = false;
        if (mag_phase >= mag_dt) {
            while (mag_phase >= mag_dt) mag_phase -= mag_dt;
            s.mag_tick = true;
            const Vector3f m_true = (ts.R_nb.transpose().cast<float>() * mag_world_ned);
            s.mag_ned = apply_mag_noise(m_true, mag_noise, mag_dt);
        }
        // Vertical subproblem: the prescribed attitude levels the measured
        // specific force exactly, so the vertical channel carries the true
        // vertical acceleration, the levelled sensor bias and noise.
        const Eigen::Vector3d f_lvl = ts.R_nb * acc.cast<double>();
        s.y_vertical = f_lvl.z() + G;

        if (!full) {
            tuner.step(s.dt, static_cast<float>(-s.y_vertical));
            OperatingPoint op;
            if (oracle.valid) op = oracle;
            else if (tuner.ready()) {
                op.valid = true;
                op.tau = tuner.tau(); op.sigma_aw = tuner.sigma(); op.rs = tuner.rs(); op.ts = tuner.ts();
                op.sigma_floor = tuner.sigma_floor();
                op.noise_floor = tuner.noise_floor();
            }
            if (op.valid) for (auto* v : vert) v->set_operating_point(op);
        }

        h_true_all.push_back(ts.h_up);
        for (size_t i = 0; i < est.size(); ++i) {
            const auto t0 = std::chrono::steady_clock::now();
            est[i]->step(s);
            const auto t1 = std::chrono::steady_clock::now();
            const Estimate e = est[i]->estimate();
            Series& se = series[i];
            if (e.live) {
                se.runtime_ns += std::chrono::duration<double, std::nano>(t1 - t0).count();
                se.runtime_count++;
                if (!std::isfinite(t_live[i])) t_live[i] = ts.t;
            }
            const double eh = e.live ? (e.h_up - ts.h_up) : NAN;
            const double ev = e.live ? (e.v_up - ts.v_up) : NAN;
            const double ea = e.live ? (e.a_up - ts.a_up) : NAN;
            double n2 = NAN, n1 = NAN;
            if (e.live && std::isfinite(e.var_p)) {
                // NED errors for [p_z, v_z].
                const double ep = -eh, evv = -ev;
                const double det = e.var_p * e.var_v - e.cov_pv * e.cov_pv;
                if (det > 0) n2 = (e.var_v * ep * ep - 2 * e.cov_pv * ep * evv + e.var_p * evv * evv) / det;
                if (e.var_p > 0) n1 = ep * ep / e.var_p;
            }
            if (!se.failed && e.live &&
                (!std::isfinite(eh) || !std::isfinite(e.var_p) || e.var_p < 0 ||
                 std::abs(eh) > 10.0 * std::max(1.0, hs_inc))) {
                se.failed = true;
                se.fail_reason = !std::isfinite(eh) ? "nonfinite" : (e.var_p < 0 ? "negative_variance" : "diverged");
            }
            se.t.push_back(static_cast<float>(ts.t));
            se.eh.push_back(static_cast<float>(eh));
            se.ev.push_back(static_cast<float>(ev));
            se.ea.push_back(static_cast<float>(ea));
            se.nees2.push_back(static_cast<float>(n2));
            se.nees1.push_back(static_cast<float>(n1));
            se.nis.push_back(static_cast<float>(e.live ? e.nis : NAN));
            se.live.push_back(e.live ? 1 : 0);
        }
        if (ts_out.is_open() && k % args.decimate == 0) {
            char buf[128];
            std::snprintf(buf, sizeof buf, "%.3f,%.5f,%.5f,%.5f,%.5f", ts.t, ts.eta_up, ts.h_up, ts.v_up, ts.a_up);
            ts_out << buf;
            for (size_t i = 0; i < est.size(); ++i) {
                const Estimate e = est[i]->estimate();
                std::snprintf(buf, sizeof buf, ",%.5f,%.5f,%.5f",
                              e.live ? e.h_up : NAN, e.live ? e.v_up : NAN, e.live ? e.a_up : NAN);
                ts_out << buf;
            }
            ts_out << "\n";
        }
    }

    // Vessel significant heave over the scored window.
    double sh = 0, sh2 = 0;
    long nh = 0;
    for (long k = 0; k < steps; ++k) {
        const double t = static_cast<double>(k) * dt;
        if (t < args.score_start) continue;
        sh += h_true_all[k]; sh2 += h_true_all[k] * h_true_all[k]; ++nh;
    }
    const double dnh = static_cast<double>(nh);
    const double hs_heave = nh ? 4.0 * std::sqrt(std::max(0.0, sh2 / dnh - (sh / dnh) * (sh / dnh))) : NAN;

    for (size_t i = 0; i < est.size(); ++i) {
        const Series& se = series[i];
        const int dof = full ? 3 : 1;
        const WindowStats w = window_stats(se, h_true_all, args.score_start, sc.duration_s, dt, dof);
        const WindowStats f = window_stats(se, h_true_all, sc.duration_s - 60.0, sc.duration_s, dt, dof);
        const double conv = convergence_time(se, t_live[i], w.rms_h, dt);
        const double hs_ref = hs_inc > 0 ? hs_inc : NAN;
        const double runtime = se.runtime_count ? se.runtime_ns / static_cast<double>(se.runtime_count) : NAN;
        std::printf("%s,%s,%s,%u,%s,%d,%zu,%.1f,%.4f,%.4f,%.2f,%.2f,%d,"
                    "%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,"
                    "%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n",
                    args.level.c_str(), full ? "tuner" : args.params.c_str(), sc.name.c_str(), args.seed,
                    est[i]->name().c_str(), est[i]->states(), est[i]->memory_bytes(), runtime, hs_inc, hs_heave,
                    t_live[i], conv, se.failed ? 1 : 0,
                    w.rms_h, 100.0 * w.rms_h / hs_ref, w.mean_h, w.peak_h, w.rms_v, w.rms_a, w.lf_rms_h,
                    w.trend_m_per_hr, w.nees2, w.nees1, w.nis, w.nis_in95,
                    f.rms_h, f.mean_h, f.peak_h, f.rms_v, f.rms_a, f.nees2, f.nis);
        if (se.failed) std::cerr << "FAILURE " << est[i]->name() << ": " << se.fail_reason << "\n";
    }
    return 0;
}
