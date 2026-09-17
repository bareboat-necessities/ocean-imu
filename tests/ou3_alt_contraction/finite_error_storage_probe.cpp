// Actual finite errors generated through a sensor prefix, never state installs.
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
int main(int argc, char** argv) {
    try {
        if (argc != 5) return 1;
        const std::string mode(argv[1]);
        const int axis = std::stoi(argv[2]);
        const float amplitude = std::stof(argv[3]);
        if ((mode != "H" && mode != "A" && mode != "HA") || axis < 0 || axis > 2 ||
            !std::isfinite(amplitude) || amplitude <= 0) return 1;
        std::ofstream out(argv[4]);
        if (!out) return 1;
        out << std::setprecision(17);
        Wrapper f;
        Wrapper::Config config;
        f.begin(config);
        f.raw().setAccBiasHold(mode != "A");
        const Eigen::Vector3f quiet(0, 0, -g_std), mag(32, 0, 0);
        int tick = 0;
        auto advance = [&](const Eigen::Vector3f& gyro) {
            ++tick;
            f.update(.005f, gyro, quiet);
            if (tick % 8 == 0) f.updateMag(mag);
        };
        while (tick < 70000) {
            advance(Eigen::Vector3f::Zero());
            if (f.isLive() && f.mag_ref_set_ && !f.raw().accel_bias_locked_) break;
        }
        if (tick == 70000) throw std::runtime_error("no startup root");
        const int startup = tick;
        const auto field = f.raw().mekf().v2ref.eval();
        const bool refined = f.mag_refine_done_;
        Eigen::Vector3f pulse = Eigen::Vector3f::Zero();
        pulse(axis) = amplitude;
        // Predeclared half-second gyro disturbance; physical source stays quiet.
        for (int i = 0; i < 100; ++i) advance(pulse);
        for (int step = 0; step <= 1800; ++step) {
            if (step) {
                if (mode == "HA" && step == 601) f.raw().setAccBiasHold(false);
                advance(Eigen::Vector3f::Zero());
            }
            const auto& c = f.raw().mekf();
            if (!f.isLive() || f.raw().accelVibrationGuardEngagement() != 0 ||
                c.use_imu_lever_arm_ || !c.Pext.allFinite() || !c.xext.allFinite() ||
                !c.qref.coeffs().allFinite() || c.qref.w() == 0 ||
                f.mag_refine_done_ != refined || !(c.v2ref-field).isZero(0))
                throw std::runtime_error("left declared finite-error diagnostic branch");
            // Physical q_WB=identity. c=2*vec(q_true*conj(q_hat))/scalar.
            // All other true states and true BIAS0 beta are zero.
            out << step << ' ' << c.acc_bias_updates_enabled();
            for (int i = 0; i < 3; ++i) out << ' ' << -2.0*double(c.qref.vec()(i))/double(c.qref.w());
            for (int i = 3; i < 21; ++i) out << ' ' << -double(c.xext(i));
            for (int i = 0; i < 21; ++i)
                for (int j = 0; j < 21; ++j) out << ' ' << double(c.Pext(i,j));
            out << '\n';
        }
        std::cout << "{\"startup_samples\":" << startup << ",\"pulse_samples\":100,\"amplitude\":"
                  << std::setprecision(17) << double(amplitude) << "}\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
}
