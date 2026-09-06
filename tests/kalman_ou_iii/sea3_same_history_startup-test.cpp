// Executes one legal complete-SEA3 continuum-driver member through the real
// shipping startup path.  This is a feasibility witness, not a source model:
// the mathematical source is the continuum Hilbert-ball member documented in
// tools/stability/ou3_sea3_fixed_history_source_core.py.  The quadrature here
// evaluates that continuum integral; its nodes are never source modes.
#define EIGEN_NON_ARDUINO

#include <algorithm>
#include <cmath>
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
constexpr double H_M = 1.5;
constexpr double TP_S = 6.0;
constexpr double GAMMA = 3.3;
constexpr double FC_HZ = 0.5;
constexpr double DRIVER_CENTER_HZ = 1.0 / TP_S;
constexpr double DRIVER_SIGMA_HZ = 0.005;
constexpr double DRIVER_BETA = 0.5;
constexpr double SIGMA_LO = 0.07;
constexpr double SIGMA_HI = 0.09;
constexpr double PM_EXPONENT = 1.25;

struct ContinuumPacket {
    std::vector<double> f;
    std::vector<double> coeff;

    static double shape(double hz) {
        if (!(hz > 0.0)) return 0.0;
        const double fp = 1.0 / TP_S;
        const double x = hz / fp;
        const double sigma = x <= 1.0 ? SIGMA_LO : SIGMA_HI;
        const double peak = std::exp(-((x - 1.0) * (x - 1.0)) /
                                     (2.0 * sigma * sigma));
        return std::pow(hz, -5.0) *
               std::exp(-PM_EXPONENT * std::pow(x, -4.0)) *
               std::pow(GAMMA, peak);
    }

    static double driverRaw(double hz) {
        const double z = (hz - DRIVER_CENTER_HZ) / DRIVER_SIGMA_HZ;
        return std::exp(-0.5 * z * z);
    }

    static double rao(double hz) {
        if (hz <= FC_HZ) return 1.0;
        return std::pow(FC_HZ / hz, 2.0);
    }

    explicit ContinuumPacket(int panels) {
        if (panels <= 0 || (panels % 2) != 0) std::abort();
        const double lo = (1.0 / TP_S) / 64.0;
        const double hi = (1.0 / TP_S) * 256.0;
        const double a = std::log(lo);
        const double b = std::log(hi);
        const double h = (b - a) / static_cast<double>(panels);

        std::vector<double> base_weight;
        f.reserve(static_cast<std::size_t>(panels + 1));
        base_weight.reserve(static_cast<std::size_t>(panels + 1));
        double raw_spectrum = 0.0;
        double raw_driver_norm2 = 0.0;
        for (int i = 0; i <= panels; ++i) {
            const double logf = a + static_cast<double>(i) * h;
            const double hz = std::exp(logf);
            const double sw = (i == 0 || i == panels) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
            const double measure = (h / 3.0) * sw * hz;
            f.push_back(hz);
            base_weight.push_back(measure);
            raw_spectrum += measure * shape(hz);
            const double dr = driverRaw(hz);
            raw_driver_norm2 += measure * dr * dr;
        }
        const double spectrum_scale = (H_M * H_M / 16.0) / raw_spectrum;
        const double driver_c = 1.0 / std::sqrt(raw_driver_norm2);

        coeff.resize(f.size());
        for (std::size_t i = 0; i < f.size(); ++i) {
            const double omega = 2.0 * PI * f[i];
            const double acc_transfer = -(omega * omega) * rao(f[i]);
            const double spectral = std::sqrt(std::max(0.0, spectrum_scale * shape(f[i])));
            const double driver = DRIVER_BETA * driver_c * driverRaw(f[i]);
            coeff[i] = base_weight[i] * spectral * acc_transfer * driver;
        }
    }

    double acceleration(double t) const {
        double total = 0.0;
        for (std::size_t i = 0; i < f.size(); ++i) {
            total += coeff[i] * std::cos(2.0 * PI * f[i] * t);
        }
        return total;
    }
};

bool check(bool condition, const char* message) {
    if (!condition) std::cerr << "FAIL: " << message << '\n';
    return condition;
}

}  // namespace

int main() {
    bool ok = true;

    // The source is the exact continuum member; these two quadratures only
    // verify that its numerical evaluation is stable at the handoff/window.
    ContinuumPacket coarse(1024);
    ContinuumPacket fine(2048);
    double max_eval_delta = 0.0;
    double max_eval_mag = 0.0;
    for (int k = 0; k <= 600; k += 20) {
        const double t = PREHISTORY_S + static_cast<double>(k) * DT;
        const double yc = coarse.acceleration(t);
        const double yf = fine.acceleration(t);
        max_eval_delta = std::max(max_eval_delta, std::fabs(yc - yf));
        max_eval_mag = std::max(max_eval_mag, std::fabs(yf));
    }
    ok &= check(max_eval_delta / std::max(max_eval_mag, 1e-30) < 2e-4,
                "continuum source quadrature must converge at the Live word");

    Filter f(false);
    f.setWithMag(false);
    f.setOnlineTuneWarmupSec(10.0f);
    f.initialize(Eigen::Vector3f::Constant(0.2f),
                 Eigen::Vector3f::Constant(0.00157f),
                 Eigen::Vector3f::Constant(0.3f));

    const Eigen::Vector3f gyro = Eigen::Vector3f::Zero();
    double first_ready_s = -1.0;
    double max_abs_accel = 0.0;
    const int pre_steps = static_cast<int>(PREHISTORY_S / DT);
    for (int k = 0; k < pre_steps; ++k) {
        const double t = static_cast<double>(k) * DT;
        const float az = static_cast<float>(fine.acceleration(t));
        max_abs_accel = std::max(max_abs_accel, std::fabs(static_cast<double>(az)));
        const Eigen::Vector3f acc(0.0f, 0.0f, -g_std + az);
        f.updateFrontEnd(static_cast<float>(DT), gyro, acc);
        if (first_ready_s < 0.0 && f.isTunerReady()) {
            first_ready_s = (static_cast<double>(k) + 1.0) * DT;
        }
    }

    ok &= check(first_ready_s > 0.0,
                "same-history continuum SEA3 prehistory must reach TunerReady");
    ok &= check(first_ready_s < PREHISTORY_S,
                "TunerReady must be reached before the declared Live handoff");
    ok &= check(f.wavePeriodUsable(),
                "shipping WavePeriodEstimator must be usable at Live entry");
    ok &= check(f.startupProxyInitialized(),
                "shipping startup proxy must be initialized at Live entry");
    ok &= check(max_abs_accel <= 4.0,
                "same-history source must remain inside Normal-Live acceleration cap");
    ok &= check(f.accelVibrationGuardEngagement() == 0.0f,
                "proof branch requires dormant-transparent vibration guard");

    const float tau_live = f.getTauApplied();
    const float sigma_live = f.getSigmaApplied();
    const float rs_live = f.getRSApplied();
    const float wave_period_live = f.getWavePeriodSec();
    ok &= check(std::isfinite(tau_live) && tau_live > 0.0f,
                "same-history startup must produce positive applied tau");
    ok &= check(std::isfinite(sigma_live) && sigma_live > 0.0f,
                "same-history startup must produce positive applied sigma");
    ok &= check(std::isfinite(rs_live) && rs_live > 0.0f,
                "same-history startup must produce positive applied R_S");
    ok &= check(std::isfinite(wave_period_live) && wave_period_live > 0.0f,
                "same-history startup must produce a measured wave period");

    // The explicit handoff is the shipping operation.  No MEKF propagation has
    // occurred before this point.  Use the ungauged-yaw seed for this H-mode
    // feasibility member; magnetic/vector events remain source events in the
    // subsequent complete word rather than an invented startup gauge.
    f.goLive(f.startupProxyQuat(), 0.035f, 1.5708f);
    ok &= check(f.isAdaptiveLive(), "shipping goLive must enter Live");

    const auto P0A = f.mekf().covariance_full();
    ok &= check(P0A.allFinite(), "shipping Live covariance seed must be finite");
    Eigen::SelfAdjointEigenSolver<decltype(P0A)> es(P0A);
    ok &= check(es.info() == Eigen::Success && es.eigenvalues().minCoeff() > 0.0f,
                "shipping A21 Live covariance seed must be SPD");

    // Execute the same continuum phase immediately after handoff.  This is not
    // yet the nonlinear P4 ratio; it verifies that the real shipping front end
    // carries its committed schedule and actual R_S into the word with no
    // source reseed.
    double min_rs = rs_live;
    double max_rs = rs_live;
    for (int k = 0; k < 601; ++k) {
        const double source_t = PREHISTORY_S + static_cast<double>(k) * DT;
        const float az = static_cast<float>(fine.acceleration(source_t));
        const Eigen::Vector3f acc(0.0f, 0.0f, -g_std + az);
        f.updateTime(static_cast<float>(DT), gyro, acc);
        const double rs = f.getRSApplied();
        min_rs = std::min(min_rs, rs);
        max_rs = std::max(max_rs, rs);
    }
    ok &= check(std::isfinite(min_rs) && min_rs > 0.0 && std::isfinite(max_rs),
                "every shipping sample must retain a finite positive applied R_S");

    std::cout << "SEA3_SAME_HISTORY_STARTUP first_ready_s=" << first_ready_s
              << " wave_period_s=" << wave_period_live
              << " tau=" << tau_live
              << " sigma=" << sigma_live
              << " RS_entry=" << rs_live
              << " RS_word_min=" << min_rs
              << " RS_word_max=" << max_rs
              << " max_abs_accel=" << max_abs_accel
              << " quadrature_rel="
              << max_eval_delta / std::max(max_eval_mag, 1e-30)
              << '\n';

    if (!ok) return 1;
    std::cout << "complete SEA3 same-history startup/live seed check passed\n";
    return 0;
}
