// Passive zero-residual tangent diagnostic on an analytic stationary source.
// The shipping wrapper, tuner, covariance, scheduler and magnetic gates run.
// This does not perturb or install an estimator state or certify finite errors.
#define EIGEN_NON_ARDUINO 1
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
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

namespace ou3_alt_probe {
template<class F> void record(int, const F&);
template<class F, class V, class N, class S, class K, class Y>
void innovation(int, const F&, const V&, const N&, const S&, const K&, const Y&);
template<class F> struct Scope {
    using C = std::remove_reference_t<F>;
    int k;
    const C& f;
    Scope(int k, const C& f): k(k), f(f) { record(k, f); }
    ~Scope() { record(k + 1, f); }
};
}
#define private public
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std = 9.80665f;
using Wrapper = SeaStateFusion_OU_III<TrackerType::KALMANF>;
using TangentMatrix = Eigen::Matrix<double, 21, 21>;
using Rows = Eigen::Matrix<double, 3, 21>;
static bool observing = false;
static int current_step = 0;
static std::ofstream stream;

template<class M> void emit(const M& m) {
    for (int i = 0; i < m.rows(); ++i)
        for (int j = 0; j < m.cols(); ++j) stream << ' ' << double(m(i, j));
}

template<class F> void require_zero(const F& f) {
    if (!f.Pext.allFinite() || !f.xext.isZero(0) ||
        !f.qref.vec().isZero(0) || f.qref.w() != 1 ||
        f.use_imu_lever_arm_ || f.wind_heel_rad_ != 0)
        throw std::runtime_error("source left the exact zero-residual tangent branch");
}

template<class F> void covariance(int kind, const F& f) {
    require_zero(f);
    stream << current_step << ' ' << kind << ' ' << f.acc_bias_updates_enabled();
    emit(f.Pext);
    stream << '\n';
}

namespace ou3_alt_probe {
template<class F> void record(int kind, const F& f) {
    if (!observing) return;
    require_zero(f);
    if (kind == 2) {
        // These are the factors computed by this very prediction, before
        // pending covariance inflation, symmetrization and S service.
        TangentMatrix A = TangentMatrix::Identity();
        A.block<6, 6>(0, 0) = f.F_AA_scratch_.template cast<double>();
        A.block<12, 12>(6, 6) = f.F_LL_scratch_.template cast<double>();
        const float phi = f.acc_bias_updates_enabled()
            ? std::exp(-f.last_dt_ / std::max(1e-3f, f.tau_bacc_)) : 1.f;
        A.block<3, 3>(18, 18) *= double(phi);
        stream << current_step << " 2";
        emit(A);
        stream << '\n';
    }
    // Covariance-only floors and held/released bias transitions are retained
    // even though their mean tangent is identity at this zero state.
    if (kind == 2 || kind == 3 || kind == 4 || kind == 11 || kind == 21 ||
        kind == 31 || kind == 61) covariance(1000 + kind, f);
}

template<class F, class V, class N, class S, class K, class Y>
void innovation(int kind, const F& f, const V& r, const N&, const S&,
                const K& gain, const Y&) {
    if (!observing) return;
    require_zero(f);
    if (!r.isZero(0)) throw std::runtime_error("nonzero innovation: tangent requires reset terms");
    Rows H = Rows::Zero();
    Eigen::Matrix3d R;
    if (kind == 12) {
        const auto force = (f.R_wb() * Eigen::Vector3f(0, 0, -f.gravity_magnitude_)).eval();
        H.block<3, 3>(0, 0) = -f.skew_symmetric_matrix(force).template cast<double>();
        H.block<3, 3>(0, 15) = f.R_wb().template cast<double>();
        // Physical held bias still enters the residual, while shipping K has
        // zero held-bias rows. Never discard the held covariance block.
        H.block<3, 3>(0, 18).setIdentity();
        R = f.Racc.template cast<double>();
    } else if (kind == 22) {
        const auto field = (f.R_wb() * f.v2ref).eval();
        H.block<3, 3>(0, 0) = -f.skew_symmetric_matrix(field).template cast<double>();
        R = f.Rmag.template cast<double>();
    } else if (kind == 32) {
        H.block<3, 3>(0, 12).setIdentity();
        R = f.R_S.template cast<double>();
    } else throw std::runtime_error("unknown accepted measurement");
    TangentMatrix A = TangentMatrix::Identity() - gain.template cast<double>() * H;
    stream << current_step << ' ' << kind;
    emit(A);
    emit(H);
    emit(R);
    stream << '\n';
}
}

int main(int argc, char** argv) {
    try {
        if (argc != 5) return 1;
        const std::string mode = argv[1];
        const int samples = std::stoi(argv[2]);
        if ((mode != "H" && mode != "A" && mode != "HA") || samples < 600) return 1;
        stream.open(argv[3]);
        if (!stream) return 1;
        stream << std::setprecision(17);
        std::ofstream states(argv[4], std::ios::binary);
        if (!states) return 1;
        Wrapper f;
        Wrapper::Config config;
        f.begin(config);
        // Public external hold is an admitted control; no state is installed.
        f.raw().setAccBiasHold(mode != "A");
        const Eigen::Vector3f quiet(0, 0, -g_std), mag(32, 0, 0);
        int startup = 0;
        for (; startup < 70000; ) {
            ++startup;
            f.update(.005f, Eigen::Vector3f::Zero(), quiet);
            if (startup % 8 == 0) f.updateMag(mag);
            // Begin both branches at the same actual unlocked/gauged root.
            if (f.isLive() && f.mag_ref_set_ && !f.raw().accel_bias_locked_) break;
        }
        if (startup == 70000) throw std::runtime_error("startup did not reach the diagnostic root");
        if (f.raw().mekf().acc_bias_updates_enabled() != (mode == "A"))
            throw std::runtime_error("requested initial hold was not preserved");
        const bool initial_refine_done = f.mag_refine_done_;
        const Eigen::Vector3f initial_field = f.raw().mekf().v2ref;
        covariance(0, f.raw().mekf());
        observing = true;
        for (current_step = 1; current_step <= samples; ++current_step) {
            if (mode == "HA" && current_step == 601) f.raw().setAccBiasHold(false);
            f.update(.005f, Eigen::Vector3f::Zero(), quiet);
            if ((startup + current_step) % 8 == 0) f.updateMag(mag);
            // A later magnetic regauge has its own tangent. Even a zero-angle
            // nominal write cannot be silently represented by identity.
            if (f.mag_refine_done_ != initial_refine_done ||
                !(f.raw().mekf().v2ref - initial_field).isZero(0))
                throw std::runtime_error("magnetic regauge requires a separate event tangent");
            if (f.raw().accelVibrationGuardEngagement() != 0)
                throw std::runtime_error("vibration guard left dormant branch");
            covariance(100, f.raw().mekf());
            const auto& core = f.raw().mekf();
            states.write(reinterpret_cast<const char*>(core.xext.data()), 21 * sizeof(float));
            states.write(reinterpret_cast<const char*>(core.Pext.data()), 441 * sizeof(float));
            states.write(reinterpret_cast<const char*>(core.qref.coeffs().data()), 4 * sizeof(float));
        }
        std::cout << "{\"mode\":\"" << mode << "\",\"startup_samples\":" << startup
                  << ",\"samples\":" << samples << ",\"dt_s\":" << std::setprecision(17)
                  << double(.005f) << ",\"source\":\"stationary_zero_BIAS0_horizontal_field\"}\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 2;
    }
}
