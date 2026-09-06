// Canonical non-promoting complete-SEA3 same-driver startup bridge.
//
// The mathematical source is tools/stability/ou3_sea3_fixed_history_source_core.py:
// one continuum Hilbert-ball coefficient field, one admissible SEA3 partition,
// one admissible continuum RAO member, no replay, no finite harmonic source,
// and no phase reseed.  Simpson nodes below evaluate the continuum integral;
// they are not source modes.
#define EIGEN_NON_ARDUINO

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

const float g_std = 9.80665f;

namespace {
using Filter = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double DT = 0.005;
constexpr double PREHISTORY_S = 60.0;
constexpr double TP_S = 6.0;
constexpr double H_M = 1.5;
constexpr double GAMMA = 3.3;
constexpr double FC_HZ = 0.5;
constexpr double DRIVER_CENTER_HZ = 1.0 / TP_S;
constexpr double DRIVER_SIGMA_HZ = 0.005;
constexpr double DRIVER_BETA = 0.5;
constexpr double SIGMA_LO = 0.07;
constexpr double SIGMA_HI = 0.09;
constexpr double PM_EXPONENT = 1.25;

struct ContinuumPacket {
    std::vector<double> omega;
    std::vector<double> coeff;

    static double shape(double f) {
        if (!(f > 0.0)) return 0.0;
        const double fp = 1.0 / TP_S;
        const double x = f / fp;
        const double sigma = x <= 1.0 ? SIGMA_LO : SIGMA_HI;
        const double peak = std::exp(-((x - 1.0) * (x - 1.0)) /
                                     (2.0 * sigma * sigma));
        return std::pow(f, -5.0) * std::exp(-PM_EXPONENT * std::pow(x, -4.0)) *
               std::pow(GAMMA, peak);
    }
    static double driverRaw(double f) {
        const double z = (f - DRIVER_CENTER_HZ) / DRIVER_SIGMA_HZ;
        return std::exp(-0.5 * z * z);
    }
    static double rao(double f) {
        return f <= FC_HZ ? 1.0 : std::pow(FC_HZ / f, 2.0);
    }

    explicit ContinuumPacket(int panels) {
        if (panels <= 0 || panels % 2) std::abort();
        const double lo = (1.0 / TP_S) / 64.0;
        const double hi = (1.0 / TP_S) * 256.0;
        const double a = std::log(lo), b = std::log(hi);
        const double h = (b - a) / static_cast<double>(panels);
        std::vector<double> f, measure;
        f.reserve(static_cast<std::size_t>(panels + 1));
        measure.reserve(static_cast<std::size_t>(panels + 1));
        double raw_spectrum = 0.0, raw_driver_norm2 = 0.0;
        for (int i = 0; i <= panels; ++i) {
            const double fi = std::exp(a + static_cast<double>(i) * h);
            const double sw = (i == 0 || i == panels) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
            const double m = (h / 3.0) * sw * fi;
            f.push_back(fi); measure.push_back(m);
            raw_spectrum += m * shape(fi);
            raw_driver_norm2 += m * driverRaw(fi) * driverRaw(fi);
        }
        const double spectrum_scale = (H_M * H_M / 16.0) / raw_spectrum;
        const double driver_c = 1.0 / std::sqrt(raw_driver_norm2);
        omega.resize(f.size()); coeff.resize(f.size());
        for (std::size_t i = 0; i < f.size(); ++i) {
            omega[i] = 2.0 * PI * f[i];
            const double acc_transfer = -(omega[i] * omega[i]) * rao(f[i]);
            coeff[i] = measure[i] * std::sqrt(std::max(0.0, spectrum_scale * shape(f[i]))) *
                       (DRIVER_BETA * driver_c * driverRaw(f[i])) * acc_transfer;
        }
    }
    double acceleration(double t) const {
        double s = 0.0;
        for (std::size_t i = 0; i < coeff.size(); ++i) s += coeff[i] * std::cos(omega[i] * t);
        return s;
    }
};

bool require(bool cond, const char* message) {
    if (!cond) std::cerr << "FAIL: " << message << '\n';
    return cond;
}
}

int main() {
    bool ok = true;
    ContinuumPacket coarse(1024), fine(2048);

    double max_delta = 0.0, max_word_abs = 0.0;
    for (int k = 0; k <= 600; k += 10) {
        const double t = PREHISTORY_S + static_cast<double>(k) * DT;
        const double a = coarse.acceleration(t), b = fine.acceleration(t);
        max_delta = std::max(max_delta, std::fabs(a - b));
        max_word_abs = std::max(max_word_abs, std::fabs(b));
    }
    const double quadrature_rel = max_delta / std::max(max_word_abs, 1e-30);
    ok &= require(quadrature_rel < 2e-4, "continuum quadrature did not converge");

    Filter filter(false);
    filter.setWithMag(false);
    filter.setOnlineTuneWarmupSec(10.0f);
    filter.initialize(Eigen::Vector3f::Constant(0.2f),
                      Eigen::Vector3f::Constant(0.00157f),
                      Eigen::Vector3f::Constant(0.3f));
    const Eigen::Vector3f gyro = Eigen::Vector3f::Zero();

    double first_ready_s = -1.0;
    double max_prehistory_abs = 0.0;
    const int pre_steps = static_cast<int>(PREHISTORY_S / DT);
    for (int k = 0; k < pre_steps; ++k) {
        const double t = static_cast<double>(k) * DT;
        const float az = static_cast<float>(fine.acceleration(t));
        max_prehistory_abs = std::max(max_prehistory_abs, std::fabs(static_cast<double>(az)));
        filter.updateFrontEnd(static_cast<float>(DT), gyro,
                              Eigen::Vector3f(0.0f, 0.0f, -g_std + az));
        if (first_ready_s < 0.0 && filter.isTunerReady())
            first_ready_s = (static_cast<double>(k) + 1.0) * DT;
    }

    ok &= require(first_ready_s > 0.0 && first_ready_s < PREHISTORY_S,
                  "same-history source did not reach TunerReady before handoff");
    ok &= require(filter.wavePeriodUsable(), "WPE is not usable at same-history handoff");
    ok &= require(filter.startupProxyInitialized(), "startup proxy is not initialized");
    ok &= require(filter.accelVibrationGuardEngagement() == 0.0f,
                  "vibration guard must remain dormant-transparent");
    ok &= require(max_prehistory_abs <= 4.0,
                  "same-history source exceeded Normal-Live acceleration cap");

    const double period = filter.getWavePeriodSec();
    const double tau = filter.getTauApplied();
    const double sigma = filter.getSigmaApplied();
    const double rs_entry = filter.getRSApplied();
    ok &= require(std::isfinite(period) && period > 0.0, "invalid measured period");
    ok &= require(std::isfinite(tau) && tau > 0.0, "invalid applied tau");
    ok &= require(std::isfinite(sigma) && sigma > 0.0, "invalid applied sigma");
    ok &= require(std::isfinite(rs_entry) && rs_entry > 0.0, "invalid applied R_S");

    filter.goLive(filter.startupProxyQuat(), 0.035f, 1.5708f);
    ok &= require(filter.isAdaptiveLive(), "shipping goLive did not enter Live");
    const auto P0 = filter.mekf().covariance_full();
    Eigen::SelfAdjointEigenSolver<decltype(P0)> es(P0);
    ok &= require(P0.allFinite() && es.info() == Eigen::Success && es.eigenvalues().minCoeff() > 0.0f,
                  "shipping A21 Live covariance seed is not SPD");

    double min_rs = rs_entry, max_rs = rs_entry;
    std::size_t rs_change_count = 0;
    double previous_rs = rs_entry;
    for (int k = 0; k < 601; ++k) {
        const double t = PREHISTORY_S + static_cast<double>(k) * DT;
        const float az = static_cast<float>(fine.acceleration(t));
        filter.updateTime(static_cast<float>(DT), gyro,
                          Eigen::Vector3f(0.0f, 0.0f, -g_std + az));
        const double rs = filter.getRSApplied();
        min_rs = std::min(min_rs, rs); max_rs = std::max(max_rs, rs);
        if (std::fabs(rs - previous_rs) > 1e-9 * std::max(1.0, std::fabs(previous_rs))) ++rs_change_count;
        previous_rs = rs;
    }
    ok &= require(std::isfinite(min_rs) && min_rs > 0.0 && std::isfinite(max_rs),
                  "applied R_S became invalid in the 601-sample word");

    std::cout << std::setprecision(17)
              << "SEA3_SAME_DRIVER_STARTUP"
              << " first_ready_s=" << first_ready_s
              << " period_s=" << period
              << " tau=" << tau
              << " sigma=" << sigma
              << " RS_entry=" << rs_entry
              << " RS_word_min=" << min_rs
              << " RS_word_max=" << max_rs
              << " RS_change_count=" << rs_change_count
              << " P0_min_eig=" << es.eigenvalues().minCoeff()
              << " quadrature_rel=" << quadrature_rel
              << " max_source_accel=" << std::max(max_prehistory_abs, max_word_abs)
              << '\n';
    return ok ? 0 : 1;
}
