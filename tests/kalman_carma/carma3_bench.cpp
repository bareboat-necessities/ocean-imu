/*
  Copyright (c) 2026 Mikhail Grushinskiy

  carma3_bench -- paired three-axis CARMA benchmark (EXPERIMENTAL).

  One process runs one (scenario, seed): an independent random-phase
  directional sea through vessel RAOs (CarmaBenchTruth.h, horizontal model
  "physical" by default; no CARMA model in the truth), corrupted once by the
  shared IMU and magnetometer noise models, fed identically to every filter.
  Surge, sway and heave are scored against the WAVE coordinate (oscillatory
  motion about the local equilibrium); a low-frequency translation of that
  equilibrium, when the scenario has one, reaches the IMU but is not part of
  the scored truth.

  --level full     complete INS behind the shared proxy startup and magnetic
                   acquisition (carma3/FullAdapters.h):
                     OU3   shipping OU-III (unchanged)
                     Z41 Z52      PR #704 hybrids: OU-III x/y, CARMA(4,1)/(5,2) z
                     A41 A52 A30  CARMA on all three axes
                     M45          CARMA(4,1) x/y, CARMA(5,2) z
                     M54          CARMA(5,2) x/y, CARMA(4,1) z
                     A52S         CARMA(5,2) all axes + S_x, S_y (pole law)
                     A41S         CARMA(4,1) all axes + S_x, S_y with OU-III's
                                  r_S carried over (negative control)
                     H41          CARMA(4,1) all axes, horizontal shape fitted
                                  to the benchmark's own sway RAO (hull-dependent)
                     A41C:<r>     A41 with horizontal noise correlation r
                     Q41 Q52      the three-axis class with the OU-III
                                  horizontal chain (equivalence with Z41, Z52)
  --level l1       prescribed attitude, per-axis reduced filters on
                   y = R_nb f + g e_z = a + R_nb b + n, one shared
                   measurement-only operating point (the vertical channel's
                   shipping tuner, as the INS); same filter names.
  --level truth    realized motion levels of the scenario.

  Scenario: <base>[+modifier...], base as carma_bench (CarmaBenchTruth.h
  base_scenario).  Modifiers: T<s> duration, rel<deg> wave propagation
  direction relative to the bow (0 following, 90 beam, 180 head; default 30),
  U<m/s> steady forward speed (encounter-frequency Doppler), anchor, drift,
  maneuver (low-frequency translation), maggap (no magnetometer for 300 s from
  T/3), jit (+-25 % sample-interval jitter), calbias (+-0.01 m/s^2 turn-on
  bias), norao (unit heave RAO), hlegacy (legacy 0.8 x heave RAO horizontal
  model).

  The complete INS has no declination input, so its horizontal axes are
  magnetic; estimates are rotated into the truth's geographic frame before
  scoring (as W3dSimulationRunner does), and the yaw error is taken against
  the magnetic heading.
*/

#define EIGEN_NON_ARDUINO

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "util/W3dSimCommon.h"
#include "kalman_carma/CarmaVerticalModel.h"
#include "CarmaBenchTruth.h"
#include "CarmaBenchLevel1.h"
#include "carma3/Bench3.h"
#include "carma3/FullAdapters.h"

using namespace carma_bench;
using namespace carma3_bench;
using Eigen::Vector3d;
using Eigen::Vector3f;
namespace kc = ocean_imu::kalman::carma;

namespace {

// ---------------------------------------------------------------------------
// Scenarios
// ---------------------------------------------------------------------------

struct ScenarioX {
    Scenario sc;
    bool jitter = false, calbias = false, maggap = false;
};

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, sep)) if (!item.empty()) out.push_back(item);
    return out;
}

ScenarioX make_scenario(const std::string& full) {
    ScenarioX x;
    Scenario& sc = x.sc;
    sc.name = full;
    const auto tok = split(full, '+');
    if (tok.empty() || !base_scenario(tok[0], sc)) throw std::invalid_argument("unknown scenario " + full);
    sc.hmodel.physical = true;
    double rel_deg = NAN;
    for (size_t i = 1; i < tok.size(); ++i) {
        const std::string& m = tok[i];
        if (m[0] == 'T' && std::isdigit(static_cast<unsigned char>(m[1]))) sc.duration_s = std::atof(m.c_str() + 1);
        else if (m.rfind("rel", 0) == 0) rel_deg = std::atof(m.c_str() + 3);
        else if (m[0] == 'U' && (std::isdigit(static_cast<unsigned char>(m[1])) || m[1] == '-')) sc.speed_mps = std::atof(m.c_str() + 1);
        else if (m == "anchor") sc.lf.kind = Translation::Kind::Anchor;
        else if (m == "drift") sc.lf.kind = Translation::Kind::Drift;
        else if (m == "maneuver") sc.lf.kind = Translation::Kind::Maneuver;
        else if (m == "maggap") x.maggap = true;
        else if (m == "jit") x.jitter = true;
        else if (m == "calbias") x.calbias = true;
        else if (m == "norao") sc.heave_rao = false;
        else if (m == "hlegacy") sc.hmodel.physical = false;
        else throw std::invalid_argument("unknown scenario modifier " + m);
    }
    // An anchored hull lies head to the waves unless told otherwise.
    if (sc.lf.kind == Translation::Kind::Anchor && !std::isfinite(rel_deg)) rel_deg = 180.0;
    if (std::isfinite(rel_deg) && !sc.parts.empty()) sc.att.heading_deg = sc.parts[0].mean_dir_deg - rel_deg;
    return x;
}

// ---------------------------------------------------------------------------
// Level 1: prescribed attitude, per-axis reduced filters
// ---------------------------------------------------------------------------

class IAxis {
public:
    virtual ~IAxis() = default;
    virtual void set_operating_point(const OperatingPoint& op) = 0;
    virtual void step(float dt, float y) = 0;
    virtual bool started() const = 0;
    virtual double p() const = 0;
    virtual double v() const = 0;
    virtual double a() const = 0;
    virtual Eigen::Matrix2d Ppv() const = 0;
    virtual double nis() const = 0;
    virtual int states() const = 0;
};

// OU-III chain [v, p, S, a_w] + b with the SpectralMSE r_S times the axis factor.
class AxisOU3 final : public IAxis {
public:
    AxisOU3(float R, BiasPrior bp, float rs_factor) : R_(R), bp_(bp), rs_factor_(rs_factor) {}
    void set_operating_point(const OperatingPoint& op) override {
        if (!started_) {
            started_ = true;
            kf_.P.setZero();
            kf_.P(0,0) = 1.0f; kf_.P(1,1) = 400.0f; kf_.P(2,2) = 2500.0f;
            kf_.P(3,3) = op.sigma_aw * op.sigma_aw; kf_.P(4,4) = bp_.sigma0 * bp_.sigma0;
        }
        op_ = op;
    }
    void step(float dt, float y) override {
        if (!started_) return;
        const float sig = std::max(op_.sigma_floor, op_.sigma_aw);
        Eigen::Matrix<float,4,4> Phi, Qd;
        ocean_imu::kalman::ou_detail::IntegratedOUChain<float,3>::transition(op_.tau, dt, Phi);
        ocean_imu::kalman::ou_detail::IntegratedOUChain<float,3>::process_covariance(op_.tau, dt, sig * sig, Qd);
        Eigen::Matrix<float,5,5> F = Eigen::Matrix<float,5,5>::Zero(), Q = F;
        F.topLeftCorner<4,4>() = Phi;
        Q.topLeftCorner<4,4>() = Qd;
        F(4,4) = std::exp(-dt / bp_.tau);
        Q(4,4) = bp_.q * (-0.5f * bp_.tau * std::expm1(-2.0f * dt / bp_.tau));
        kf_.predict(F, Q);
        if (ocean_imu::kalman::ou_detail::periodic_update_due(dt, op_.ts, s_elapsed_)) {
            const int idx[1] = {2};
            const float c[1] = {1.0f};
            const float nis = kf_.last_nis;
            const float rs = rs_factor_ * op_.rs;
            kf_.update(idx, c, 0.0f, rs * rs);
            kf_.last_nis = nis;
        }
        const int idx[2] = {3, 4};
        const float c[2] = {1.0f, 1.0f};
        kf_.update(idx, c, y, R_);
    }
    bool started() const override { return started_; }
    double p() const override { return kf_.x(1); }
    double v() const override { return kf_.x(0); }
    double a() const override { return kf_.x(3); }
    Eigen::Matrix2d Ppv() const override {
        Eigen::Matrix2d M; M << kf_.P(1,1), kf_.P(1,0), kf_.P(0,1), kf_.P(0,0); return M;
    }
    double nis() const override { return kf_.last_nis; }
    int states() const override { return 5; }

private:
    float R_;
    BiasPrior bp_;
    float rs_factor_;
    OperatingPoint op_;
    bool started_ = false;
    float s_elapsed_ = 0.0f;
    LinearKF<5> kf_;
};

// CARMA(N) [+ S] + b.  Scale and amplitude by the INS mapping; S by the
// pole-placement law when the shape has a DC level, else OU-III's r_S times
// the axis factor (the negative control).
template<int N, bool VI>
class AxisCARMA final : public IAxis {
    static constexpr int NZ = N + (VI ? 1 : 0);
    static constexpr int NS = NZ + 1;
public:
    AxisCARMA(float R, BiasPrior bp, float rs_factor, const kc::CarmaShape& shape)
        : R_(R), bp_(bp), rs_factor_(rs_factor), model_(shape) {}
    void set_operating_point(const OperatingPoint& op) override {
        op_ = op;
        const double fz = model_.base_zero_crossing_hz();
        const float kappa = static_cast<float>((0.5 / op.tau) / fz);
        const float sa = std::max(op.sigma_floor, op.sigma_aw / 0.9f);
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
                              kc::carma_unit_psd<N>(model_.shape(), model_.kappa(), 0.0);
            rs_ = (c0 > 1e-14) ? static_cast<float>(std::sqrt(c0)) * op.tau /
                                     (ocean_imu::ou3::R_S_POLE_KAPPA_DEFAULT * std::sqrt(op.ts))
                               : rs_factor_ * op.rs;
        }
    }
    void step(float dt, float y) override {
        if (!started_) return;
        Eigen::Matrix<float,NZ,NZ> Fz, Qz;
        if constexpr (VI) {
            Eigen::Matrix<float,NZ,NZ> A = Eigen::Matrix<float,NZ,NZ>::Zero(), W = A;
            A.template topLeftCorner<N,N>() = model_.A();
            A(N, 0) = 1.0f;
            W.template topLeftCorner<N,N>() = model_.W();
            kc::discretize_lti<float,NZ>(A, W, dt, Fz, Qz);
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
                std::isfinite(rs_) && rs_ > 0.0f) {
                const int idx[1] = {N};
                const float c[1] = {1.0f};
                const float nis = kf_.last_nis;
                kf_.update(idx, c, 0.0f, rs_ * rs_);
                kf_.last_nis = nis;
            }
        }
        const int idx[2] = {2, NS - 1};
        const float c[2] = {1.0f, 1.0f};
        kf_.update(idx, c, y, R_);
    }
    bool started() const override { return started_; }
    double p() const override { return kf_.x(0); }
    double v() const override { return kf_.x(1); }
    double a() const override { return kf_.x(2); }
    Eigen::Matrix2d Ppv() const override {
        Eigen::Matrix2d M; M << kf_.P(0,0), kf_.P(0,1), kf_.P(1,0), kf_.P(1,1); return M;
    }
    double nis() const override { return kf_.last_nis; }
    int states() const override { return NS; }

private:
    float R_;
    BiasPrior bp_;
    float rs_factor_;
    OperatingPoint op_;
    bool started_ = false;
    float s_elapsed_ = 0.0f;
    float rs_ = NAN;
    kc::CarmaVerticalModel<float,N> model_;
    LinearKF<NS> kf_;
};

class Level1Set final : public IEstimator3 {
public:
    Level1Set(std::string name, std::unique_ptr<IAxis> x, std::unique_ptr<IAxis> y, std::unique_ptr<IAxis> z)
        : name_(std::move(name)) { ax_[0] = std::move(x); ax_[1] = std::move(y); ax_[2] = std::move(z); }
    std::string name() const override { return name_; }
    void set_operating_point(const OperatingPoint& op) { for (auto& a : ax_) a->set_operating_point(op); }
    void step(const Sensors& s) override {
        for (int k = 0; k < 3; ++k) ax_[k]->step(s.dt, static_cast<float>(s.y_level(k)));
    }
    Estimate3 estimate() const override {
        Estimate3 e;
        e.live = ax_[0]->started() && ax_[1]->started() && ax_[2]->started();
        e.Ppv.setZero();
        for (int k = 0; k < 3; ++k) {
            e.p(k) = ax_[k]->p(); e.v(k) = ax_[k]->v(); e.a(k) = ax_[k]->a();
            const Eigen::Matrix2d M = ax_[k]->Ppv();
            e.Ppv(k, k) = M(0,0); e.Ppv(k, 3 + k) = M(0,1); e.Ppv(3 + k, k) = M(1,0); e.Ppv(3 + k, 3 + k) = M(1,1);
        }
        e.nis = ax_[2]->nis();
        e.nis_dof = 1;
        return e;
    }
    size_t memory_bytes() const override { return sizeof(*this); }
    int states() const override { return ax_[0]->states() + ax_[1]->states() + ax_[2]->states(); }

private:
    std::string name_;
    std::unique_ptr<IAxis> ax_[3];
};

std::unique_ptr<Level1Set> make_level1(const std::string& n, float R, BiasPrior bp) {
    const float fx = 0.72f, fy = 0.50f;
    const kc::CarmaShape c30 = kc::CARMA30_JONSWAP_SHAPE, c41 = kc::CARMA41_JONSWAP_SHAPE,
                         c52 = kc::CARMA52_JONSWAP_SHAPE, h41 = CARMA41_SWAY_SHAPE;
    auto ou = [&](float f) { return std::make_unique<AxisOU3>(R, bp, f); };
    auto c = [&](int order, bool vi, float f, const kc::CarmaShape& sh) -> std::unique_ptr<IAxis> {
        if (order == 3) return std::make_unique<AxisCARMA<3,false>>(R, bp, f, sh);
        if (order == 4) return vi ? std::unique_ptr<IAxis>(std::make_unique<AxisCARMA<4,true>>(R, bp, f, sh))
                                  : std::unique_ptr<IAxis>(std::make_unique<AxisCARMA<4,false>>(R, bp, f, sh));
        return vi ? std::unique_ptr<IAxis>(std::make_unique<AxisCARMA<5,true>>(R, bp, f, sh))
                  : std::unique_ptr<IAxis>(std::make_unique<AxisCARMA<5,false>>(R, bp, f, sh));
    };
    if (n == "OU3") return std::make_unique<Level1Set>(n, ou(fx), ou(fy), ou(1.0f));
    if (n == "Z41") return std::make_unique<Level1Set>(n, ou(fx), ou(fy), c(4, false, 1, c41));
    if (n == "Z52") return std::make_unique<Level1Set>(n, ou(fx), ou(fy), c(5, false, 1, c52));
    if (n == "A30") return std::make_unique<Level1Set>(n, c(3, false, fx, c30), c(3, false, fy, c30), c(3, false, 1, c30));
    if (n == "A41") return std::make_unique<Level1Set>(n, c(4, false, fx, c41), c(4, false, fy, c41), c(4, false, 1, c41));
    if (n == "A52") return std::make_unique<Level1Set>(n, c(5, false, fx, c52), c(5, false, fy, c52), c(5, false, 1, c52));
    if (n == "A52S") return std::make_unique<Level1Set>(n, c(5, true, fx, c52), c(5, true, fy, c52), c(5, false, 1, c52));
    if (n == "A41S") return std::make_unique<Level1Set>(n, c(4, true, fx, c41), c(4, true, fy, c41), c(4, false, 1, c41));
    if (n == "M45") return std::make_unique<Level1Set>(n, c(4, false, fx, c41), c(4, false, fy, c41), c(5, false, 1, c52));
    if (n == "M54") return std::make_unique<Level1Set>(n, c(5, false, fx, c52), c(5, false, fy, c52), c(4, false, 1, c41));
    if (n == "H41") return std::make_unique<Level1Set>(n, c(4, false, fx, h41), c(4, false, fy, h41), c(4, false, 1, c41));
    return nullptr;
}

// ---------------------------------------------------------------------------
// Metrics
// ---------------------------------------------------------------------------

struct Series {
    std::vector<float> t, e[3], ev[3], pp[3], nees[3], nees3, nees6, nis, droll, dpitch, dyaw;
    double runtime_ns = 0.0;
    long runtime_count = 0;
    bool failed = false;
    std::string fail_reason;
    double t_live = NAN;
};

struct AxisStats {
    double rms = NAN, mean = NAN, peak = NAN, lf = NAN, trend = NAN, rmsv = NAN, nees = NAN, predstd = NAN;
};

AxisStats axis_stats(const Series& s, int k, double t0, double t1, double dt) {
    AxisStats w;
    double se = 0, sm = 0, peak = 0, sv = 0, sn = 0, spp = 0, st = 0, stt = 0, ste = 0;
    long n = 0, nn = 0;
    std::vector<double> win;
    for (size_t i = 0; i < s.t.size(); ++i) {
        if (s.t[i] < t0 || s.t[i] > t1 || !std::isfinite(s.e[k][i])) continue;
        const double e = s.e[k][i];
        se += e * e; sm += e; peak = std::max(peak, std::abs(e));
        sv += double(s.ev[k][i]) * s.ev[k][i];
        st += s.t[i]; stt += double(s.t[i]) * s.t[i]; ste += s.t[i] * e;
        if (std::isfinite(s.nees[k][i])) { sn += s.nees[k][i]; spp += s.pp[k][i]; ++nn; }
        win.push_back(e);
        ++n;
    }
    if (n == 0) return w;
    const double dn = static_cast<double>(n);
    w.rms = std::sqrt(se / dn); w.mean = sm / dn; w.peak = peak; w.rmsv = std::sqrt(sv / dn);
    const double den = dn * stt - st * st;
    if (den > 0) w.trend = (dn * ste - st * sm) / den * 3600.0;
    if (nn > 0) { w.nees = sn / static_cast<double>(nn); w.predstd = std::sqrt(spp / static_cast<double>(nn)); }
    const size_t L = static_cast<size_t>(std::lround(60.0 / dt));
    double acc = 0, slf = 0;
    long nlf = 0;
    for (size_t j = 0; j < win.size(); ++j) {
        acc += win[j];
        if (j >= L) acc -= win[j - L];
        if (j + 1 >= L) { const double m = acc / static_cast<double>(L); slf += m * m; ++nlf; }
    }
    if (nlf > 0) w.lf = std::sqrt(slf / static_cast<double>(nlf));
    return w;
}

// First time after Live from which every 30 s rolling RMS error stays within
// 1.5 x the scored RMS (NaN-safe: the rolling window starts at Live).
double convergence_time(const Series& s, int k, double rms_ref, double dt) {
    if (!std::isfinite(s.t_live) || !(rms_ref > 0)) return NAN;
    const size_t L = static_cast<size_t>(std::lround(30.0 / dt));
    const size_t N = s.t.size();
    std::vector<double> roll(N, NAN);
    double acc = 0;
    size_t n = 0;
    for (size_t i = 0; i < N; ++i) {
        if (!std::isfinite(s.e[k][i])) { acc = 0; n = 0; continue; }
        acc += double(s.e[k][i]) * s.e[k][i];
        ++n;
        if (n > L) acc -= double(s.e[k][i - L]) * s.e[k][i - L];
        if (n >= L) roll[i] = std::sqrt(std::max(0.0, acc / static_cast<double>(L)));
    }
    double t_conv = NAN;
    for (size_t i = N; i-- > 0;) {
        if (s.t[i] < s.t_live) break;
        if (!(roll[i] <= 1.5 * rms_ref)) break;
        t_conv = s.t[i] - 30.0;
    }
    return std::isfinite(t_conv) ? std::max(0.0, t_conv - s.t_live) : NAN;
}

double window_mean(const std::vector<float>& x, const std::vector<float>& t, double t0, double t1, bool square) {
    double s = 0; long n = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        if (t[i] < t0 || t[i] > t1 || !std::isfinite(x[i])) continue;
        s += square ? double(x[i]) * x[i] : x[i]; ++n;
    }
    return n ? s / static_cast<double>(n) : NAN;
}

struct Args {
    std::string scenario = "jonswap_fp0.29_hs0.8_g3.3";
    unsigned seed = 1;
    std::string level = "full";
    std::string filters;
    std::string timeseries;
    int decimate = 10;
    bool header = false;
    double score_start = 300.0;
    double startup_window = 120.0;
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
        else if (k == "--filters") a.filters = val();
        else if (k == "--timeseries") a.timeseries = val();
        else if (k == "--decimate") a.decimate = std::stoi(val());
        else if (k == "--score-start") a.score_start = std::stod(val());
        else if (k == "--header") a.header = true;
        else throw std::invalid_argument("unknown argument " + k);
    }
    return a;
}

void print_header() {
    std::printf("level,scenario,seed,filter,states,memory_bytes,runtime_ns,hs_incident,hs_heave,"
                "std_x,std_y,std_z,lf_std_x,lf_std_y,t_live,failed");
    for (const char* ax : {"x", "y", "z"})
        std::printf(",rms_%s,pct_%s,mean_%s,peak_%s,lf_%s,trend_%s,rmsv_%s,nees_%s,predstd_%s,conv_%s,"
                    "su_rms_%s,su_peak_%s,f60_%s", ax, ax, ax, ax, ax, ax, ax, ax, ax, ax, ax, ax, ax);
    std::printf(",rms_h,nees_p3,nees_pv6,nis,nis_in95,roll_rms_deg,pitch_rms_deg,yaw_rms_deg,yaw_mean_deg\n");
}

int truth_stats(const ScenarioX& sx, const Args& args, double dt, long steps) {
    const Scenario& sc = sx.sc;
    WaveTruth truth(sc, w3d_expand_seed(args.seed, 21u), dt);
    double s[3] = {}, sa[3] = {}, slf[2] = {}, slfa[2] = {}, sr = 0, sp = 0, sy = 0, wmax = 0;
    for (long k = 0; k < steps; ++k) {
        const TruthSample ts = truth.next();
        for (int i = 0; i < 3; ++i) { s[i] += ts.p_wave(i) * ts.p_wave(i); sa[i] += ts.a_wave(i) * ts.a_wave(i); }
        for (int i = 0; i < 2; ++i) {
            const double d = ts.p_ned(i) - ts.p_wave(i), da = ts.a_ned(i) - ts.a_wave(i);
            slf[i] += d * d; slfa[i] += da * da;
        }
        const double r = (ts.roll - sc.att.heel_deg * PI / 180.0) * 180.0 / PI;
        const double q = (ts.pitch - sc.att.trim_deg * PI / 180.0) * 180.0 / PI;
        sr += r * r; sp += q * q;
        sy += std::pow((ts.yaw - sc.att.heading_deg * PI / 180.0) * 180.0 / PI, 2);
        wmax = std::max(wmax, ts.omega_b.norm());
    }
    const double n = static_cast<double>(steps);
    std::printf("scenario=%s seed=%u std_wave_m=%.3f,%.3f,%.3f accel_wave_rms=%.3f,%.3f,%.3f "
                "lf_disp_rms_m=%.2f,%.2f lf_accel_rms=%.4f,%.4f roll_rms_deg=%.2f pitch_rms_deg=%.2f "
                "yaw_dev_rms_deg=%.2f rate_max=%.3f heading_deg=%.1f\n",
                args.scenario.c_str(), args.seed, std::sqrt(s[0] / n), std::sqrt(s[1] / n), std::sqrt(s[2] / n),
                std::sqrt(sa[0] / n), std::sqrt(sa[1] / n), std::sqrt(sa[2] / n),
                std::sqrt(slf[0] / n), std::sqrt(slf[1] / n), std::sqrt(slfa[0] / n), std::sqrt(slfa[1] / n),
                std::sqrt(sr / n), std::sqrt(sp / n), std::sqrt(sy / n), wmax, sc.att.heading_deg);
    return 0;
}

double wrap_pi(double a) {
    while (a > PI) a -= TWO_PI;
    while (a < -PI) a += TWO_PI;
    return a;
}

void euler_zyx(const Eigen::Quaterniond& q, double& roll, double& pitch, double& yaw) {
    const Eigen::Matrix3d R = q.toRotationMatrix();
    pitch = std::asin(std::max(-1.0, std::min(1.0, -R(2,0))));
    roll = std::atan2(R(2,1), R(2,2));
    yaw = std::atan2(R(1,0), R(0,0));
}

}  // namespace

int main(int argc, char** argv) {
    Args args;
    ScenarioX sx;
    try {
        args = parse(argc, argv);
        if (args.header) {
            print_header();
            if (argc == 2) return 0;
        }
        sx = make_scenario(args.scenario);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
    const Scenario& sc = sx.sc;
    const double dt = 1.0 / 200.0;
    const long steps = static_cast<long>(std::lround(sc.duration_s / dt));
    if (args.level == "truth") return truth_stats(sx, args, dt, steps);
    const bool full = args.level == "full";
    if (!full && args.level != "l1") {
        std::cerr << "ERROR: --level must be full, l1 or truth\n";
        return 2;
    }

    const unsigned base = args.seed;
    ImuNoiseModel acc_noise = make_imu_noise_model(ACC_SIGMA, sx.calbias ? 0.01f : 8e-3f * 9.80665f, 0.0005f,
                                                   w3d_expand_seed(base, 11u), w3d_expand_seed(base, 12u));
    ImuNoiseModel gyr_noise = make_imu_noise_model(GYR_SIGMA, 0.05f * 3.14159265f / 180.0f, 0.00001f,
                                                   w3d_expand_seed(base, 13u), w3d_expand_seed(base, 14u));
    MagNoiseModel mag_noise = make_mag_noise_model(MAG_SIGMA_UT, 2.0f, 0.01f, 0.015f, 0.010f, 1.0f,
                                                   w3d_expand_seed(base, 15u), w3d_expand_seed(base, 16u));
    const Vector3f mag_world_ned = MagSim_WMM::mag_world_aero();
    const double dec = MagSim_WMM::default_declination_deg * PI / 180.0;
    // Magnetic -> geographic NED: the filter's x axis points at azimuth dec.
    Eigen::Matrix3d R_mg = Eigen::Matrix3d::Identity();
    R_mg(0,0) = std::cos(dec); R_mg(0,1) = -std::sin(dec);
    R_mg(1,0) = std::sin(dec); R_mg(1,1) = std::cos(dec);
    const float mag_dt = 1.0f / 25.0f;
    float mag_phase = 0.0f;
    const double gap0 = sc.duration_s / 3.0, gap1 = gap0 + 300.0;

    WaveTruth truth(sc, w3d_expand_seed(base, 21u), dt);
    const double hs_inc = sc.stationary ? 0.0
        : (sc.fade.enabled ? WaveTruth::incident_hs(sc.parts_b) : WaveTruth::incident_hs(sc.parts));

    std::vector<std::unique_ptr<IEstimator3>> est;
    std::vector<Level1Set*> l1;
    const float R_acc = std::pow(2.8f * ACC_SIGMA * SIGMA_A_RESCALE, 2.0f);
    std::vector<std::string> names = split(args.filters, ',');
    if (names.empty()) names = {"OU3", "Z41", "Z52", "A41", "A52", "M45", "A52S", "A41S", "H41"};
    for (const auto& n : names) {
        if (full) {
            auto e = make_full(n);
            if (!e) { std::cerr << "ERROR: unknown full filter " << n << "\n"; return 2; }
            est.push_back(std::move(e));
        } else {
            auto e = make_level1(n, R_acc, BiasPrior{});
            if (!e) { std::cerr << "ERROR: unknown level-1 filter " << n << "\n"; return 2; }
            l1.push_back(e.get());
            est.push_back(std::move(e));
        }
    }

    VerticalTuner tuner;
    std::vector<Series> series(est.size());
    for (auto& s : series) {
        s.t.reserve(steps);
        for (int k = 0; k < 3; ++k) { s.e[k].reserve(steps); s.ev[k].reserve(steps); s.pp[k].reserve(steps); s.nees[k].reserve(steps); }
    }
    std::vector<Vector3d> pw_all;
    pw_all.reserve(steps);
    double lf2[2] = {0, 0};
    long nlf = 0;

    std::ofstream ts_out;
    if (!args.timeseries.empty()) {
        ts_out.open(args.timeseries);
        ts_out << "t,x_true,y_true,z_true,x_total,y_total";
        for (const auto& e : est) ts_out << ",x_" << e->name() << ",y_" << e->name() << ",z_" << e->name();
        for (const auto& e : est) ts_out << ",fz_" << e->name();
        ts_out << "\n";
    }

    std::mt19937 jitter_rng(w3d_expand_seed(base, 31u));
    std::uniform_real_distribution<double> jitter_u(-0.25, 0.25);
    double t_now = 0.0;
    for (long k = 0; k < steps; ++k) {
        const double dt_k = (sx.jitter && k > 0) ? dt * (1.0 + jitter_u(jitter_rng)) : dt;
        if (k > 0) t_now += dt_k;
        const TruthSample ts = sx.jitter ? truth.at(t_now) : truth.next();
        Vector3f acc = apply_imu_noise(ts.f_b.cast<float>(), acc_noise, static_cast<float>(dt_k));
        acc += ts.extra_acc_bias_b.cast<float>();
        const Vector3f gyr = apply_imu_noise(ts.omega_b.cast<float>(), gyr_noise, static_cast<float>(dt_k));
        Sensors s;
        s.dt = static_cast<float>(dt_k);
        s.acc_ned = acc;
        s.gyr_ned = gyr;
        mag_phase += s.dt;
        if (mag_phase >= mag_dt) {
            while (mag_phase >= mag_dt) mag_phase -= mag_dt;
            const Vector3f m_true = ts.R_nb.transpose().cast<float>() * mag_world_ned;
            const Vector3f m = apply_mag_noise(m_true, mag_noise, mag_dt);
            s.mag_tick = !(sx.maggap && ts.t >= gap0 && ts.t < gap1);
            s.mag_ned = m;
        }
        s.y_level = ts.R_nb * acc.cast<double>();
        s.y_level.z() += G;

        if (!full) {
            tuner.step(s.dt, static_cast<float>(-s.y_level.z()));
            if (tuner.ready()) {
                OperatingPoint op;
                op.valid = true;
                op.tau = tuner.tau(); op.sigma_aw = tuner.sigma(); op.rs = tuner.rs(); op.ts = tuner.ts();
                op.sigma_floor = tuner.sigma_floor();
                op.noise_floor = tuner.noise_floor();
                for (auto* e : l1) e->set_operating_point(op);
            }
        }

        pw_all.push_back(ts.p_wave);
        if (ts.t >= args.score_start) {
            for (int i = 0; i < 2; ++i) lf2[i] += std::pow(ts.p_ned(i) - ts.p_wave(i), 2);
            ++nlf;
        }
        for (size_t i = 0; i < est.size(); ++i) {
            const auto t0 = std::chrono::steady_clock::now();
            est[i]->step(s);
            const auto t1 = std::chrono::steady_clock::now();
            Estimate3 e = est[i]->estimate();
            Series& se = series[i];
            if (e.live) {
                se.runtime_ns += std::chrono::duration<double, std::nano>(t1 - t0).count();
                se.runtime_count++;
                if (!std::isfinite(se.t_live)) se.t_live = ts.t;
            }
            se.t.push_back(static_cast<float>(ts.t));
            double rr = NAN, pr = NAN, yr = NAN;
            if (e.live) {
                Vector3d p = e.p, v = e.v;
                Eigen::Matrix<double,6,6> P = e.Ppv;
                if (full) {
                    // Into the geographic frame.
                    Eigen::Matrix<double,6,6> T6 = Eigen::Matrix<double,6,6>::Zero();
                    T6.topLeftCorner<3,3>() = R_mg; T6.bottomRightCorner<3,3>() = R_mg;
                    p = R_mg * p; v = R_mg * v;
                    P = (T6 * P * T6.transpose()).eval();
                }
                const Vector3d ep = p - ts.p_wave, ev = v - ts.v_wave;
                for (int a = 0; a < 3; ++a) {
                    se.e[a].push_back(static_cast<float>(ep(a)));
                    se.ev[a].push_back(static_cast<float>(ev(a)));
                    se.pp[a].push_back(static_cast<float>(P(a, a)));
                    Eigen::Matrix2d M;
                    M << P(a, a), P(a, 3 + a), P(3 + a, a), P(3 + a, 3 + a);
                    const Eigen::Vector2d x(ep(a), ev(a));
                    const double det = M.determinant();
                    se.nees[a].push_back(det > 0 ? static_cast<float>(x.dot(M.inverse() * x)) : NAN);
                }
                const Eigen::Matrix3d P3 = P.topLeftCorner<3,3>();
                Eigen::LDLT<Eigen::Matrix3d> l3(P3);
                se.nees3.push_back(l3.info() == Eigen::Success && l3.isPositive() ? static_cast<float>(ep.dot(l3.solve(ep))) : NAN);
                Eigen::Matrix<double,6,1> x6; x6 << ep, ev;
                Eigen::LDLT<Eigen::Matrix<double,6,6>> l6(P);
                se.nees6.push_back(l6.info() == Eigen::Success && l6.isPositive() ? static_cast<float>(x6.dot(l6.solve(x6))) : NAN);
                se.nis.push_back(static_cast<float>(e.nis));
                if (e.has_att) {
                    double r1, p1, y1;
                    euler_zyx(e.q_bw, r1, p1, y1);
                    rr = wrap_pi(r1 - ts.roll); pr = wrap_pi(p1 - ts.pitch);
                    yr = wrap_pi(y1 + dec - ts.yaw);
                }
                if (!se.failed && (!ep.allFinite() || !P.allFinite() || P(2,2) < 0 || P(0,0) < 0 || P(1,1) < 0 ||
                                   ep.cwiseAbs().maxCoeff() > 10.0 * std::max(1.0, hs_inc))) {
                    se.failed = true;
                    se.fail_reason = !ep.allFinite() ? "nonfinite" : (ep.cwiseAbs().maxCoeff() > 10.0 * std::max(1.0, hs_inc)
                                                                          ? "diverged" : "negative_variance");
                }
            } else {
                for (int a = 0; a < 3; ++a) {
                    se.e[a].push_back(NAN); se.ev[a].push_back(NAN); se.pp[a].push_back(NAN); se.nees[a].push_back(NAN);
                }
                se.nees3.push_back(NAN); se.nees6.push_back(NAN); se.nis.push_back(NAN);
            }
            se.droll.push_back(static_cast<float>(rr));
            se.dpitch.push_back(static_cast<float>(pr));
            se.dyaw.push_back(static_cast<float>(yr));
        }
        if (ts_out.is_open() && k % args.decimate == 0) {
            char buf[160];
            std::snprintf(buf, sizeof buf, "%.3f,%.5f,%.5f,%.5f,%.4f,%.4f", ts.t, ts.p_wave.x(), ts.p_wave.y(),
                          ts.p_wave.z(), ts.p_ned.x(), ts.p_ned.y());
            ts_out << buf;
            for (size_t i = 0; i < est.size(); ++i) {
                const Series& se = series[i];
                const size_t j = se.t.size() - 1;
                std::snprintf(buf, sizeof buf, ",%.5f,%.5f,%.5f",
                              se.e[0][j] + ts.p_wave.x(), se.e[1][j] + ts.p_wave.y(), se.e[2][j] + ts.p_wave.z());
                ts_out << buf;
            }
            for (const auto& e : est) {
                std::snprintf(buf, sizeof buf, ",%.4f", e->op_fz_hz());
                ts_out << buf;
            }
            ts_out << "\n";
        }
    }

    // Truth statistics over the scored window.
    double m[3] = {}, m2[3] = {};
    long nw = 0;
    for (long k = 0; k < steps; ++k) {
        if (static_cast<double>(k) * dt < args.score_start) continue;
        for (int a = 0; a < 3; ++a) { m[a] += pw_all[k](a); m2[a] += pw_all[k](a) * pw_all[k](a); }
        ++nw;
    }
    double std_w[3];
    for (int a = 0; a < 3; ++a) {
        const double mu = nw ? m[a] / static_cast<double>(nw) : 0.0;
        std_w[a] = nw ? std::sqrt(std::max(0.0, m2[a] / static_cast<double>(nw) - mu * mu)) : NAN;
    }
    const double hs_heave = 4.0 * std_w[2];
    const double lfx = nlf ? std::sqrt(lf2[0] / static_cast<double>(nlf)) : 0.0;
    const double lfy = nlf ? std::sqrt(lf2[1] / static_cast<double>(nlf)) : 0.0;

    for (size_t i = 0; i < est.size(); ++i) {
        const Series& se = series[i];
        const double runtime = se.runtime_count ? se.runtime_ns / static_cast<double>(se.runtime_count) : NAN;
        std::printf("%s,%s,%u,%s,%d,%zu,%.1f,%.4f,%.4f,%.5f,%.5f,%.5f,%.3f,%.3f,%.2f,%d",
                    args.level.c_str(), sc.name.c_str(), args.seed, est[i]->name().c_str(), est[i]->states(),
                    est[i]->memory_bytes(), runtime, hs_inc, hs_heave, std_w[0], std_w[1], std_w[2], lfx, lfy,
                    se.t_live, se.failed ? 1 : 0);
        const double hs_ref = hs_inc > 0 ? hs_inc : NAN;
        double rms_h2 = 0;
        for (int a = 0; a < 3; ++a) {
            const AxisStats w = axis_stats(se, a, args.score_start, sc.duration_s, dt);
            const AxisStats f60 = axis_stats(se, a, sc.duration_s - 60.0, sc.duration_s, dt);
            const AxisStats su = axis_stats(se, a, se.t_live, se.t_live + args.startup_window, dt);
            const double conv = convergence_time(se, a, w.rms, dt);
            if (a < 2) rms_h2 += w.rms * w.rms;
            std::printf(",%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.2f,%.6g,%.6g,%.6g",
                        w.rms, 100.0 * w.rms / hs_ref, w.mean, w.peak, w.lf, w.trend, w.rmsv, w.nees, w.predstd,
                        conv, su.rms, su.peak, f60.rms);
        }
        double nis_in = 0; long nn = 0;
        const int dof = full ? 3 : 1;
        const double lo = dof == 1 ? 0.000982 : 0.2158, hi = dof == 1 ? 5.0239 : 9.3484;
        for (size_t j = 0; j < se.t.size(); ++j) {
            if (se.t[j] < args.score_start || !std::isfinite(se.nis[j])) continue;
            nis_in += (se.nis[j] >= lo && se.nis[j] <= hi) ? 1.0 : 0.0; ++nn;
        }
        const double t0 = args.score_start, t1 = sc.duration_s;
        std::printf(",%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n",
                    std::sqrt(rms_h2), window_mean(se.nees3, se.t, t0, t1, false),
                    window_mean(se.nees6, se.t, t0, t1, false), window_mean(se.nis, se.t, t0, t1, false),
                    nn ? nis_in / static_cast<double>(nn) : NAN,
                    std::sqrt(window_mean(se.droll, se.t, t0, t1, true)) * 180.0 / PI,
                    std::sqrt(window_mean(se.dpitch, se.t, t0, t1, true)) * 180.0 / PI,
                    std::sqrt(window_mean(se.dyaw, se.t, t0, t1, true)) * 180.0 / PI,
                    window_mean(se.dyaw, se.t, t0, t1, false) * 180.0 / PI);
        if (se.failed) std::cerr << "FAILURE " << est[i]->name() << ": " << se.fail_reason << "\n";
    }
    return 0;
}
