// Coupled analytic physical source from construction; passive finite storage probe.
#define EIGEN_NON_ARDUINO 1
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include "util/W3dSimCommon.h"
#include "kalman_ou_common/KalmanOUCoreMath.h"
#define private public
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private
const float g_std = 9.80665f;
using Wrapper = SeaStateFusion_OU_III<TrackerType::KALMANF>;

// One scalar harmonic in a fixed direction, not independently selected samples.
struct Physical {
    Eigen::Vector3d p, v, a, primitive, beta;
};
Physical physical(double t, double omega, int family) {
    const Eigen::Vector3d direction(0, .6, .8);
    Physical x;
    x.p = .25*std::cos(omega*t)*direction;
    x.v = -.25*omega*std::sin(omega*t)*direction;
    x.a = -.25*omega*omega*std::cos(omega*t)*direction;
    x.primitive = .25/omega*std::sin(omega*t)*direction;
    const double pi = std::acos(-1.0);
    double beta;
    if (family == 0) beta = .02 + .005*std::sin(2*pi*t/600);
    else if (family == 1) beta = .08*std::exp(-t/1200)+.015*std::sin(2*pi*t/600);
    else beta = .05+.01*std::sin(2*pi*t/600);
    x.beta = beta*Eigen::Vector3d(0, 1, 0);
    return x;
}
int main(int argc, char** argv) {
    try {
        if (argc != 5) return 1;
        const std::string mode(argv[1]);
        const int family = std::stoi(argv[2]);
        const double omega = std::stod(argv[3]);
        if ((mode != "H" && mode != "A" && mode != "HA") || family < 0 || family > 2 ||
            (omega != .5 && omega != 1)) return 1;
        std::ofstream out(argv[4]);
        if (!out) return 1;
        out << std::setprecision(17);
        Wrapper f;
        Wrapper::Config config;
        f.begin(config);
        f.raw().setAccBiasHold(mode != "A");
        const double h = double(.005f);
        const Eigen::Vector3f mag(32, 0, 0);
        int tick = 0, live_tick = 0;
        Eigen::Vector3d origin = Eigen::Vector3d::Zero();
        double sensor_rounding_max = 0;
        auto advance = [&]() {
            ++tick;
            const auto truth = physical(tick*h, omega, family);
            Eigen::Vector3d acc = truth.a + truth.beta;
            acc.z() -= double(g_std);
            const Eigen::Vector3f rounded = acc.cast<float>();
            sensor_rounding_max = std::max(sensor_rounding_max, (rounded.cast<double>()-acc).norm());
            const bool was_live = f.isLive();
            f.update(.005f, Eigen::Vector3f::Zero(), rounded);
            if (!was_live && f.isLive()) {
                if (live_tick) throw std::runtime_error("second Live origin");
                live_tick = tick;
                origin = truth.primitive;
            }
            if (tick % 8 == 0) f.updateMag(mag);
        };
        while (tick < 70000) {
            advance();
            if (f.isLive() && f.mag_ref_set_ && !f.raw().accel_bias_locked_) break;
        }
        if (tick == 70000) throw std::runtime_error("no diagnostic root before experiment cap (not universal timeout)");
        const int startup = tick;
        if (!live_tick) throw std::runtime_error("missing physical Live origin");
        for (int step = 0; step <= 1800; ++step) {
            if (step) {
                if (mode == "HA" && step == 601) f.raw().setAccBiasHold(false);
                advance();
            }
            const auto& c = f.raw().mekf();
            if (!f.isLive() || f.raw().accelVibrationGuardEngagement() != 0 ||
                c.use_imu_lever_arm_ || !c.Pext.allFinite() || !c.xext.allFinite() ||
                !c.qref.coeffs().allFinite() || c.qref.w() == 0 ||
                c.acc_bias_updates_enabled() != (mode == "A" || (mode == "HA" && step >= 601)))
                throw std::runtime_error("left declared diagnostic branch");
            const auto truth = physical(tick*h, omega, family);
            Eigen::Matrix<double,21,1> target = Eigen::Matrix<double,21,1>::Zero();
            target.segment<3>(6) = truth.v;
            target.segment<3>(9) = truth.p;
            target.segment<3>(12) = truth.primitive-origin;
            target.segment<3>(15) = truth.a;
            target.segment<3>(18) = truth.beta;
            out << step << ' ' << c.acc_bias_updates_enabled();
            for (int i=0; i<3; ++i) out << ' ' << -2.0*double(c.qref.vec()(i))/double(c.qref.w());
            for (int i=3; i<21; ++i) out << ' ' << target(i)-double(c.xext(i));
            for (int i=0; i<3; ++i) out << ' ' << truth.beta(i);
            for (int i=0; i<21; ++i)
                for (int j=0; j<21; ++j) out << ' ' << double(c.Pext(i,j));
            // Physical-only candidate supply coordinates: p,v,S,a,beta.
            for (int i=6; i<21; ++i) out << ' ' << target(i);
            out << '\n';
        }
        std::cout << "{\"startup_samples\":" << startup << ",\"live_sample\":" << live_tick
                  << ",\"dt_s\":" << std::setprecision(17) << h
                  << ",\"max_acc_serialization_residual\":" << sensor_rounding_max << "}\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
}
