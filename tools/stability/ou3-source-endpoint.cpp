// External-generator CI probe: one physical source and one unmodified shipping observer, from power-on.
// This records an ACTUAL (generally forced) finite endpoint, not a frozen-gain
// shadow or a replacement homogeneous observer. No window/direction search.
#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "util/W3dSimCommon.h"
#include "kalman_ou_common/KalmanOUCoreMath.h"

// Read-only host access, as in ou3-operation-ledger-sim.cpp. No writes to
// private estimator state, no reseeding P or intervention at window entrance.
#define private public
#include "PiersonMoskowitzStokes3D_Waves.h"
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

extern const float g_std = 9.80665f;

namespace {
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;
using Sea = PMStokesN3dWaves<128, 3>;
using V3 = Eigen::Vector3d;
using M3 = Eigen::Matrix3d;
using V12 = Eigen::Matrix<double, 12, 1>;
constexpr float dt = 1.0f / 200.0f;

template<class Derived>
void array(std::ostream& out, const Eigen::MatrixBase<Derived>& x)
{
    out << '[';
    for (Eigen::Index i = 0; i < x.rows(); ++i) {
        for (Eigen::Index j = 0; j < x.cols(); ++j) {
            if (i || j) out << ',';
            const double value = static_cast<double>(x(i,j));
            if (!std::isfinite(value)) throw std::runtime_error("nonfinite capture");
            out << value;
        }
    }
    out << ']';
}

struct Truth {
    M3 R;
    V12 linear;
    V3 acc;
    V3 gyro;
};

Truth source(const Sea& sea, double t)
{
    const M3 basis = basis_zu_to_ned().cast<double>();
    const auto state = sea.getLagrangianState(t);
    const auto imu = sea.getIMUReadings(0.0, 0.0, t, 0.0, double(dt));
    V3 integral = V3::Zero();
    V3 drift = V3::Zero();
    for (int i = 0; i < 128; ++i) {
        const double a = sea.A1_(i), k = sea.k_(i), omega = sea.omega_(i);
        const V3 direction(sea.dir_x_(i), sea.dir_y_(i), 0.0);
        drift += omega * k * a * a * direction;
        const std::array<double, 3> coeff = {a, k*a*a/2, 3*k*k*a*a*a/8};
        for (int n = 1; n <= 3; ++n) {
            const double theta = double(n) * (-omega*t + sea.phi_(i));
            integral += coeff[std::size_t(n-1)] / (double(n)*omega)
                * (std::sin(theta)*direction + std::cos(theta)*V3::UnitZ());
        }
    }
    Truth truth;
    truth.R = basis * sea.rotationMatrixAt(0.0, 0.0, t) * basis.transpose();
    // The source's constant Stokes drift is not the derivative of its
    // centered displacement. Use the same centered orbital v/p/S graph.
    truth.linear.segment<3>(0) = basis * (state.velocity - drift);
    truth.linear.segment<3>(3) = basis * state.displacement;
    truth.linear.segment<3>(6) = basis * integral;
    truth.linear.segment<3>(9) = basis * state.acceleration;
    truth.acc = basis * imu.accel_body;
    truth.gyro = basis * imu.gyro_body;
    return truth;
}

void checkpoint(std::ostream& out, const Fusion& fusion, const Truth& truth,
                const char* word, const char* event, unsigned index, double t,
                bool s_due, const Eigen::Matrix3f& applied_rs,
                const V12& prediction_forcing)
{
    const auto& raw = fusion.raw();
    const auto& m = raw.mekf();
    M3 unheel;
    for (int j = 0; j < 3; ++j)
        unheel.col(j) = m.deheel_vector_(Eigen::Vector3f::Unit(j)).cast<double>();
    const M3 Rtrue = unheel * truth.R;
    out << "{\"word\":\"" << word << "\",\"event\":\"" << event
        << "\",\"index\":" << index << ",\"source_time\":" << t
        << ",\"live\":" << fusion.isLive()
        << ",\"active\":" << m.acc_bias_updates_enabled()
        << ",\"mag_lock\":" << fusion.hasMagNorthLock()
        << ",\"mag_refined\":" << fusion.hasRefinedMagReference()
        << ",\"acc_accepted\":" << m.lastAccDiag().accepted
        << ",\"s_due\":" << s_due
        << ",\"tau_aw\":" << raw.getTauApplied()
        << ",\"tau_b\":" << m.tau_bacc_
        << ",\"S_period\":" << m.pseudo_update_period_s_
        << ",\"S_elapsed\":" << m.pseudo_update_elapsed_s_
        << ",\"wind_heel\":" << m.wind_heel_rad_
        << ",\"R_true\":";
    array(out, Rtrue);
    out << ",\"R_hat\":"; array(out, m.R_wb());
    out << ",\"x_hat\":"; array(out, m.xext);
    out << ",\"linear_true\":"; array(out, truth.linear);
    out << ",\"P\":"; array(out, m.covariance_full());
    out << ",\"R_S\":"; array(out, applied_rs);
    out << ",\"prediction_forcing\":"; array(out, prediction_forcing);
    out << ",\"true_bias\":[0,0,0],\"physical_acceleration\":";
    array(out, truth.linear.tail<3>());
    out << ",\"physical_gyro\":"; array(out, truth.gyro);
    out << ",\"mag_reference\":"; array(out, m.v2ref);
    out << ",\"acc_innovation\":"; array(out, m.lastAccDiag().r);
    out << "}\n";
}
} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc != 2) throw std::runtime_error("usage: source-endpoint OUTPUT_PREFIX");
        const std::string prefix = argv[1];
        std::ofstream root(prefix + ".root.json");
        std::ofstream inputs(prefix + ".inputs.csv");
        std::ofstream points(prefix + ".prefixes.jsonl");
        if (!root || !inputs || !points) throw std::runtime_error("cannot open capture");
        root << std::setprecision(17);
        inputs << std::setprecision(17);
        points << std::setprecision(17) << std::boolalpha;
        auto spread = std::make_shared<Cosine2sRandomizedDistribution>(
            -30.0f * std::numbers::pi / 180.0, 10.0, 42u);
        Sea sea(1.5f, 5.7f, spread, 0.02, 0.8, g_std, 42u);
        root << "{\"branch\":\"STOKES_WAVE_FOLLOWING\",\"seed\":42,"
             << "\"gravity\":" << g_std << ",\"dt\":" << dt
             << ",\"bias_root\":[0,0,0],\"bias_driver\":\"ZERO\","
             << "\"temperature\":35,\"mag_world\":[20,0,40],\"atoms\":[";
        for (int i = 0; i < 128; ++i) {
            if (i) root << ',';
            root << '[' << sea.omega_(i) << ',' << sea.k_(i) << ',' << sea.A1_(i)
                 << ',' << sea.phi_(i) << ',' << sea.dir_x_(i) << ',' << sea.dir_y_(i) << ']';
        }
        root << "]}\n";

        Fusion fusion;
        Fusion::Config cfg;
        cfg.sigma_a = Eigen::Vector3f::Constant(2.8f * 1.51e-3f * g_std * 0.71f);
        cfg.sigma_g = Eigen::Vector3f::Constant(2.0f * 0.00157f * 0.05f);
        cfg.sigma_m = Eigen::Vector3f::Constant(1.2f * 0.80f * 2.0f);
        fusion.begin(cfg);
        fusion.raw().setPeriodicAwCovarianceSync(true);
        fusion.raw().setAwCovarianceSyncCongruent(false);
        fusion.raw().enableTuner(true);
        fusion.raw().enableClamp(true);

        // Fixed original observer-clock entrances; no rescan if this newly
        // captured source does not meet the requested mode/domain there.
        constexpr std::array<float, 2> starts = {32.99784469604492f, 1143.96044921875f};
        constexpr std::array<const char*, 2> names = {"H18", "A21"};
        std::array<unsigned, 2> counts = {0, 0};
        std::array<bool, 2> started = {false, false};
        double source_time = 0.0;
        float observer_time = 0.0f, mag_elapsed = 0.0f;
        Truth previous = source(sea, -double(dt));
        inputs << "index,source_time,observer_time,ax,ay,az,gx,gy,gz,mag_tick,mx,my,mz\n";
        for (unsigned index = 0; index < 250000u; ++index) {
            const Truth truth = source(sea, source_time);
            const Eigen::Vector3f acc = truth.acc.cast<float>();
            const Eigen::Vector3f gyro = truth.gyro.cast<float>();
            const Eigen::Vector3f mag = (truth.R * V3(20,0,40)).cast<float>();
            auto& raw = fusion.raw();
            auto& m = raw.mekf();
            const bool live_before = fusion.isLive();
            // Exact same stage/commit boundary used by the single observer.
            // The subsequent update sees pending=false; no second schedule.
            if (live_before) raw.apply_pending_online_tune_();
            const Eigen::Matrix3f rs = m.R_S;
            float elapsed = m.pseudo_update_elapsed_s_;
            const bool due = live_before && ocean_imu::kalman::ou_detail::periodic_update_due(
                dt, m.pseudo_update_period_s_, elapsed);
            for (std::size_t w = 0; w < 2; ++w) {
                if (!started[w] && observer_time >= starts[w] - 2.0e-6f) {
                    if (std::abs(observer_time - starts[w]) > 2.0e-4f)
                        throw std::runtime_error("fixed entrance clock missed");
                    started[w] = true;
                    checkpoint(points, fusion, previous, names[w], "root", index,
                               source_time-double(dt), false, rs, V12::Zero());
                }
            }
            fusion.update(dt, gyro, acc, 35.0f);
            // Handoff may replace the core; do not retain a pre-update core
            // reference across the literal wrapper update.
            const V12 forcing = truth.linear
                - raw.mekf().F_LL_scratch_.cast<double>() * previous.linear;
            for (std::size_t w = 0; w < 2; ++w)
                if (started[w] && counts[w] < 600u)
                    checkpoint(points, fusion, truth, names[w], "imu", index,
                               source_time, due, rs, forcing);

            mag_elapsed += dt;
            const bool mag_tick = mag_elapsed >= 1.0f/25.0f;
            if (mag_tick) {
                while (mag_elapsed >= 1.0f/25.0f) mag_elapsed -= 1.0f/25.0f;
                fusion.updateMag(mag);
                for (std::size_t w = 0; w < 2; ++w)
                    if (started[w] && counts[w] < 600u)
                        checkpoint(points, fusion, truth, names[w], "mag", index,
                                   source_time, false, rs, V12::Zero());
            }
            inputs << index << ',' << source_time << ',' << observer_time;
            for (int j = 0; j < 3; ++j) inputs << ',' << acc(j);
            for (int j = 0; j < 3; ++j) inputs << ',' << gyro(j);
            inputs << ',' << mag_tick;
            for (int j = 0; j < 3; ++j) inputs << ',' << mag(j);
            inputs << '\n';
            for (std::size_t w = 0; w < 2; ++w)
                if (started[w] && counts[w] < 600u) ++counts[w];
            if (counts[1] == 600u) break;
            previous = truth;
            source_time += double(dt);
            observer_time += dt;
        }
        if (counts[0] != 600u || counts[1] != 600u)
            throw std::runtime_error("incomplete fixed word");
        if (!root || !inputs || !points) throw std::runtime_error("capture write failed");
        std::cout << "SOURCE_ENDPOINT_CAPTURE H18=600 A21=600 single_shipping_observer=1\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
