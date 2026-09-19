// Literal magnetically informed H18 service-superword export.
//
// This harness runs the shipping SeaStateFusionFilter_OU_III through its own
// startup, hands over at the deployed goLive() call, keeps the accelerometer
// bias held with the shipping external hold, and then records one service
// superword rooted at an actually reached state. It reseeds nothing: the root
// is whatever the shipping execution arrived at.
//
// What it exports is measurement, not proof:
//
//   * the shipping covariance at every prefix of the superword,
//   * central differences of the complete closed-loop finite-error map taken
//     through the shipping code itself, so the map is the deployed one rather
//     than a re-derivation of it,
//   * the actual innovation covariance and magnetic sensitivity of every
//     correction the shipping update really applied,
//   * the reference finite-error trajectory carrying the physical bias
//     predecessors and the wave-model mismatch.
//
// Every structural claim is left to the reader of the export. Nothing here
// evaluates the target inequality; that is done separately, at high precision,
// and it is not promoted to a certificate either.
#define EIGEN_NON_ARDUINO

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Geometry>

// The export reads the shipping covariance, gain and magnetic reference that
// the deployed update actually used. Nothing is written back.
#define private public
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#undef private

const float g_std = 9.80665f;

namespace {

using Filter = SeaStateFusionFilter_OU_III<TrackerType::KALMANF>;
using V3d = Eigen::Vector3d;
using V3f = Eigen::Vector3f;
using Qd = Eigen::Quaterniond;

constexpr int NE = 21;                    // shipping extended error dimension
constexpr double DT = 1.0 / 200.0;        // deployed IMU sample period
constexpr int MAG_STRIDE = 20;            // 10 Hz magnetometer service
constexpr double GRAVITY = 9.80665;

// One deterministic marine history. Displacement is a finite sum of zero-mean
// sinusoids, so its primitive is uniformly bounded and permanent displacement
// DC is excluded by construction. Angular rate, residual accelerometer bias
// and residual gyro bias are bounded and rate bounded.
struct Harmonic {
    double amp, w, phase;
};

constexpr Harmonic DISP[3][2] = {
    {{0.60, 0.80, 0.00}, {0.30, 1.30, 0.70}},
    {{0.50, 0.90, 0.30}, {0.20, 1.60, 1.20}},
    {{0.90, 0.70, 1.00}, {0.40, 1.10, 0.20}},
};
constexpr Harmonic RATE[3] = {{0.10, 0.90, 0.00}, {0.08, 0.70, 0.40}, {0.05, 0.30, 1.10}};
constexpr Harmonic ACC_BIAS[3] = {{0.050, 0.015, 0.00}, {0.040, 0.012, 1.57}, {0.030, 0.010, 0.50}};
constexpr Harmonic GYR_BIAS[3] = {{0.0020, 0.0040, 0.00}, {0.0015, 0.0030, 1.57}, {0.0010, 0.0035, 0.90}};
constexpr Harmonic MAG_RESIDUAL[3] = {{0.50, 2.10, 0.00}, {0.40, 1.70, 0.60}, {0.30, 1.30, 1.40}};

const V3d MAG_WORLD(20.0, 2.0, 43.0);   // micro tesla, NED

V3d harmonic_value(const Harmonic (&h)[3], double t) {
    return V3d(h[0].amp * std::sin(h[0].w * t + h[0].phase),
               h[1].amp * std::sin(h[1].w * t + h[1].phase),
               h[2].amp * std::sin(h[2].w * t + h[2].phase));
}

V3d displacement(double t) {
    V3d out = V3d::Zero();
    for (int axis = 0; axis < 3; ++axis)
        for (const Harmonic& h : DISP[axis])
            out(axis) += h.amp * std::sin(h.w * t + h.phase);
    return out;
}

V3d velocity(double t) {
    V3d out = V3d::Zero();
    for (int axis = 0; axis < 3; ++axis)
        for (const Harmonic& h : DISP[axis])
            out(axis) += h.amp * h.w * std::cos(h.w * t + h.phase);
    return out;
}

V3d acceleration(double t) {
    V3d out = V3d::Zero();
    for (int axis = 0; axis < 3; ++axis)
        for (const Harmonic& h : DISP[axis])
            out(axis) -= h.amp * h.w * h.w * std::sin(h.w * t + h.phase);
    return out;
}

// Primitive of the displacement, so that the integral state has one fixed
// physical origin chosen once at certified-tail entry.
V3d displacement_primitive(double t) {
    V3d out = V3d::Zero();
    for (int axis = 0; axis < 3; ++axis)
        for (const Harmonic& h : DISP[axis])
            out(axis) += (h.amp / h.w) * (std::cos(h.phase) - std::cos(h.w * t + h.phase));
    return out;
}

struct Sample {
    Qd q_bw;              // body to world, integrated from the same rate the gyro reports
    V3d rate;             // body angular rate
    V3d acc_body;         // specific force with the physical residual bias
    V3d mag_body;         // magnetometer with a bounded hard-iron residual
    V3d disp, vel, acc;   // world wave coordinate, its rate and its acceleration
    V3d primitive;        // displacement primitive, before the tail origin shift
    V3d acc_bias, gyr_bias;
};

std::vector<Sample> build_history(int samples) {
    std::vector<Sample> history(static_cast<size_t>(samples));
    Qd q_bw = Qd::Identity();
    for (int k = 0; k < samples; ++k) {
        const double t = static_cast<double>(k) * DT;
        Sample& s = history[static_cast<size_t>(k)];
        s.q_bw = q_bw;
        s.rate = harmonic_value(RATE, t);
        s.disp = displacement(t);
        s.vel = velocity(t);
        s.acc = acceleration(t);
        s.primitive = displacement_primitive(t);
        s.acc_bias = harmonic_value(ACC_BIAS, t);
        s.gyr_bias = harmonic_value(GYR_BIAS, t);

        const Eigen::Matrix3d R_wb = q_bw.conjugate().toRotationMatrix();
        const V3d gravity_world(0.0, 0.0, GRAVITY);
        s.acc_body = R_wb * (s.acc - gravity_world) + s.acc_bias;
        s.mag_body = R_wb * MAG_WORLD + harmonic_value(MAG_RESIDUAL, t);

        // The reported rate is held over the step, so the physical attitude and
        // the gyroscope stream describe the same motion exactly.
        const double angle = s.rate.norm() * DT;
        Qd step = Qd::Identity();
        if (angle > 1e-14) step = Qd(Eigen::AngleAxisd(angle, s.rate.normalized()));
        q_bw = q_bw * step;
        q_bw.normalize();
    }
    return history;
}

V3f to_float(const V3d& v) {
    return V3f(static_cast<float>(v.x()), static_cast<float>(v.y()), static_cast<float>(v.z()));
}

using ErrorVec = std::array<double, NE>;

// Finite estimation error in the shipping error coordinates:
// attitude (left-multiplicative on the shipping WORLD->BODY' reference),
// gyro bias, velocity, position, integral state, wave acceleration and
// accelerometer bias. Physical truth minus the deployed estimate throughout.
ErrorVec finite_error(const Filter& f, const Sample& s, const V3d& primitive_origin) {
    ErrorVec e{};
    const Eigen::Quaternionf qf = f.mekf().quaternion_boat();
    const Qd q_hat(static_cast<double>(qf.w()), static_cast<double>(qf.x()),
                   static_cast<double>(qf.y()), static_cast<double>(qf.z()));
    Qd delta = s.q_bw.conjugate() * q_hat;   // reference-frame error, world to body'
    if (delta.w() < 0.0) delta.coeffs() *= -1.0;
    const double vec_norm = delta.vec().norm();
    const double angle = 2.0 * std::atan2(vec_norm, delta.w());
    V3d theta = V3d::Zero();
    if (vec_norm > 1e-300) theta = (angle / vec_norm) * delta.vec();

    const auto& x = f.mekf().xext;
    const V3d truth_tail[6] = {
        s.gyr_bias, s.vel, s.disp, s.primitive - primitive_origin, s.acc, s.acc_bias,
    };
    for (int i = 0; i < 3; ++i) e[static_cast<size_t>(i)] = theta(i);
    for (int block = 0; block < 6; ++block) {
        for (int i = 0; i < 3; ++i) {
            const int idx = 3 + 3 * block + i;
            e[static_cast<size_t>(idx)] =
                truth_tail[block](i) - static_cast<double>(x(idx));
        }
    }
    return e;
}

struct MagEvent {
    int prefix = 0;                  // superword prefix index the correction lands on
    double time = 0.0;
    std::array<double, 3> sensitivity_axis{};   // predicted field in body coordinates
    std::array<double, 9> innovation_covariance{};
};

struct RunRecord {
    std::vector<ErrorVec> prefix_error;                 // one per superword prefix
    std::vector<std::array<double, NE * (NE + 1) / 2>> prefix_covariance;
    std::vector<MagEvent> applied_mag;
    std::array<double, 9> root_rotation{};              // world to body' at the root
    double root_time = 0.0;
    double live_time = 0.0;
    int live_prefix_samples = 0;
    int mag_attempts = 0;
    bool acc_bias_held = true;
    bool attitude_injection_valid = true;
};

// One complete shipping execution. The perturbation, when present, displaces
// the deployed estimate at the superword root and nothing else: the physical
// history, the measurement stream and every later shipping event are shared.
struct Perturbation {
    int coordinate = -1;
    double size = 0.0;
};

void perturb_estimate(Filter& f, const Perturbation& p) {
    if (p.coordinate < 0) return;
    auto& m = f.mekf();
    if (p.coordinate < 3) {
        V3f dtheta = V3f::Zero();
        dtheta(p.coordinate) = static_cast<float>(p.size);
        const Eigen::Quaternionf corr =
            ocean_imu::kalman::ou_detail::quat_from_delta_theta(dtheta);
        m.qref = corr * m.qref;
        m.qref.normalize();
    } else {
        m.xext(p.coordinate) += static_cast<float>(p.size);
    }
}

RunRecord run_execution(const std::vector<Sample>& history, int settle_samples,
                        int superword_samples, const Perturbation& perturbation,
                        bool export_events) {
    RunRecord rec;
    Filter f(true);
    f.setWithMag(true);
    f.setOnlineTuneWarmupSec(5.0f);
    f.initialize(V3f::Constant(0.0148f), V3f::Constant(0.00157f), V3f::Constant(0.25f));
    f.setAccBiasHold(true);   // deployed external hold: the superword stays H18/held
    f.setMagDelaySec(0.0f);
    f.mekf().set_mag_world_ref(to_float(MAG_WORLD));
    f.initialize_from_acc(to_float(history[0].acc_body));

    int live_index = -1;
    int root_index = -1;
    V3d primitive_origin = V3d::Zero();
    bool rooted = false;

    const int samples = static_cast<int>(history.size());
    for (int k = 0; k < samples; ++k) {
        const Sample& s = history[static_cast<size_t>(k)];
        const V3f gyro = to_float(s.rate + s.gyr_bias);
        const V3f acc = to_float(s.acc_body);

        if (f.isAdaptiveLive()) {
            f.updateTime(static_cast<float>(DT), gyro, acc);
        } else {
            f.updateFrontEnd(static_cast<float>(DT), gyro, acc);
            if (f.isTunerReady()) {
                f.goLive(f.startupProxyQuat(), 0.035f, 1.5708f);
                live_index = k;
                rec.live_time = static_cast<double>(k) * DT;
                // One fixed physical origin for the integral state, taken once
                // at certified-tail entry and never restarted afterwards.
                primitive_origin = s.primitive;
            }
        }

        const bool mag_due = f.isAdaptiveLive() && (k % MAG_STRIDE) == 0;
        if (mag_due) f.updateMag(to_float(s.mag_body));

        if (live_index >= 0 && k == live_index + settle_samples) {
            root_index = k;
            rooted = true;
            perturb_estimate(f, perturbation);
            rec.root_time = static_cast<double>(k) * DT;
            rec.live_prefix_samples = k - live_index;
            const Eigen::Matrix3f R = f.mekf().R_wb();
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    rec.root_rotation[static_cast<size_t>(3 * i + j)] =
                        static_cast<double>(R(i, j));
        }

        if (rooted) {
            const int prefix = k - root_index;
            // The error is read against the physical state the step advanced to.
            // A startup that reached Live unusually late can run the history out
            // before the superword closes; stop here so the caller reports an
            // incomplete superword rather than reading past the end.
            if (k + 1 >= samples) break;
            const Sample& next = history[static_cast<size_t>(k + 1)];
            rec.prefix_error.push_back(finite_error(f, next, primitive_origin));
            if (export_events) {
                std::array<double, NE * (NE + 1) / 2> upper{};
                size_t at = 0;
                for (int i = 0; i < NE; ++i)
                    for (int j = i; j < NE; ++j)
                        upper[at++] = static_cast<double>(f.mekf().Pext(i, j));
                rec.prefix_covariance.push_back(upper);
                // A correction applied at the root precedes the displaced estimate,
                // so it belongs to the inherited prefix rather than to the superword.
                if (prefix > 0 && mag_due && f.mekf().lastMagDiag().accepted) {
                    MagEvent ev;
                    ev.prefix = prefix;
                    ev.time = static_cast<double>(k) * DT;
                    const V3f axis = f.mekf().R_wb() * f.mekf().v2ref;
                    for (int i = 0; i < 3; ++i)
                        ev.sensitivity_axis[static_cast<size_t>(i)] =
                            static_cast<double>(axis(i));
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j)
                            ev.innovation_covariance[static_cast<size_t>(3 * i + j)] =
                                static_cast<double>(f.mekf().lastMagDiag().S(i, j));
                    rec.applied_mag.push_back(ev);
                }
            }
            if (prefix >= superword_samples) break;
        }
    }

    rec.mag_attempts = f.mag_updates_applied_;
    rec.acc_bias_held = !f.mekf().acc_bias_updates_enabled();
    rec.attitude_injection_valid = f.mekf().xext.allFinite();
    return rec;
}

// Per-block central-difference sizes. They are small against the reached
// finite error in each block and large against single-precision resolution.
double coordinate_step(int coordinate) {
    if (coordinate < 3) return 2.0e-3;        // rad
    if (coordinate < 6) return 2.0e-5;        // rad/s
    if (coordinate < 9) return 2.0e-3;        // m/s
    if (coordinate < 12) return 2.0e-3;       // m
    if (coordinate < 15) return 2.0e-3;       // m s
    if (coordinate < 18) return 2.0e-3;       // m/s^2
    return 2.0e-4;                            // m/s^2
}

void print_vector(std::FILE* out, const double* values, size_t count) {
    std::fputc('[', out);
    for (size_t i = 0; i < count; ++i) {
        if (i) std::fputc(',', out);
        std::fprintf(out, "%.17g", values[i]);
    }
    std::fputc(']', out);
}

}  // namespace

int main(int argc, char** argv) {
    std::string output;
    int settle_seconds = 30;
    double superword_seconds = 1.0;
    double conditioning_scale = 4.0;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&](double fallback) {
            return i + 1 < argc ? std::atof(argv[++i]) : fallback;
        };
        if (arg == "--output" && i + 1 < argc) output = argv[++i];
        else if (arg == "--settle-seconds") settle_seconds = static_cast<int>(next(30.0));
        else if (arg == "--superword-seconds") superword_seconds = next(1.0);
        else if (arg == "--conditioning-scale") conditioning_scale = next(4.0);
        else {
            std::fprintf(stderr, "unknown argument: %s\n", arg.c_str());
            return 2;
        }
    }

    const int settle_samples = static_cast<int>(std::lround(static_cast<double>(settle_seconds) / DT));
    const int superword_samples = static_cast<int>(std::lround(superword_seconds / DT));
    const int history_samples = 200 * 600 + settle_samples + superword_samples + 4;
    const std::vector<Sample> history = build_history(history_samples);

    const Perturbation none;
    const RunRecord reference = run_execution(history, settle_samples, superword_samples, none, true);
    if (reference.prefix_error.size() != static_cast<size_t>(superword_samples) + 1) {
        std::fprintf(stderr, "the shipping execution did not reach a complete superword\n");
        return 1;
    }
    if (reference.applied_mag.empty()) {
        std::fprintf(stderr, "no magnetic correction was actually applied in the superword\n");
        return 1;
    }

    // Central differences of the complete closed-loop map, taken through the
    // shipping code. Column i is the response to a displaced estimate in
    // coordinate i; the root response is measured rather than assumed.
    const size_t prefixes = reference.prefix_error.size();
    std::vector<std::vector<double>> difference(prefixes,
                                                std::vector<double>(NE * NE, 0.0));
    std::vector<double> endpoint_coarse(NE * NE, 0.0);
    std::vector<double> root_coarse(NE * NE, 0.0);
    for (int coordinate = 0; coordinate < NE; ++coordinate) {
        const double step = coordinate_step(coordinate);
        for (const double scale : {1.0, conditioning_scale}) {
            const Perturbation plus{coordinate, step * scale};
            const Perturbation minus{coordinate, -step * scale};
            const RunRecord up = run_execution(history, settle_samples, superword_samples, plus, false);
            const RunRecord down = run_execution(history, settle_samples, superword_samples, minus, false);
            if (up.prefix_error.size() != prefixes || down.prefix_error.size() != prefixes) {
                std::fprintf(stderr, "a displaced shipping execution did not reach the superword end\n");
                return 1;
            }
            if (scale == 1.0) {
                for (size_t prefix = 0; prefix < prefixes; ++prefix)
                    for (int row = 0; row < NE; ++row)
                        difference[prefix][static_cast<size_t>(row * NE + coordinate)] =
                            0.5 * (up.prefix_error[prefix][static_cast<size_t>(row)] -
                                   down.prefix_error[prefix][static_cast<size_t>(row)]);
            } else {
                for (int row = 0; row < NE; ++row) {
                    endpoint_coarse[static_cast<size_t>(row * NE + coordinate)] =
                        0.5 * (up.prefix_error[prefixes - 1][static_cast<size_t>(row)] -
                               down.prefix_error[prefixes - 1][static_cast<size_t>(row)]) / scale;
                    root_coarse[static_cast<size_t>(row * NE + coordinate)] =
                        0.5 * (up.prefix_error[0][static_cast<size_t>(row)] -
                               down.prefix_error[0][static_cast<size_t>(row)]) / scale;
                }
            }
        }
    }

    std::FILE* out = output.empty() ? stdout : std::fopen(output.c_str(), "w");
    if (!out) {
        std::fprintf(stderr, "cannot write %s\n", output.c_str());
        return 1;
    }

    std::fprintf(out, "{\n");
    std::fprintf(out, "  \"qualification\": \"OU3_H18_SERVICE_SUPERWORD_EXPORT_V1\",\n");
    std::fprintf(out, "  \"role\": \"literal shipping measurement; no inequality is evaluated or asserted here\",\n");
    std::fprintf(out, "  \"error_dimension\": %d,\n", NE);
    std::fprintf(out, "  \"sample_period_s\": %.17g,\n", DT);
    std::fprintf(out, "  \"magnetometer_stride_samples\": %d,\n", MAG_STRIDE);
    std::fprintf(out, "  \"live_time_s\": %.17g,\n", reference.live_time);
    std::fprintf(out, "  \"root_time_s\": %.17g,\n", reference.root_time);
    std::fprintf(out, "  \"root_after_live_samples\": %d,\n", reference.live_prefix_samples);
    std::fprintf(out, "  \"superword_samples\": %d,\n", superword_samples);
    std::fprintf(out, "  \"magnetic_attempts_total\": %d,\n", reference.mag_attempts);
    std::fprintf(out, "  \"acc_bias_held_through_superword\": %s,\n",
                 reference.acc_bias_held ? "true" : "false");
    std::fprintf(out, "  \"attitude_injection_finite\": %s,\n",
                 reference.attitude_injection_valid ? "true" : "false");
    std::fprintf(out, "  \"inherited_state\": true,\n");
    std::fprintf(out, "  \"conditioning_scale\": %.17g,\n", conditioning_scale);

    std::fprintf(out, "  \"root_rotation_world_to_body\": ");
    print_vector(out, reference.root_rotation.data(), reference.root_rotation.size());
    std::fprintf(out, ",\n");

    std::fprintf(out, "  \"reference_error\": [\n");
    for (size_t prefix = 0; prefix < prefixes; ++prefix) {
        std::fprintf(out, "    ");
        print_vector(out, reference.prefix_error[prefix].data(), NE);
        std::fprintf(out, prefix + 1 < prefixes ? ",\n" : "\n");
    }
    std::fprintf(out, "  ],\n");

    std::fprintf(out, "  \"covariance_upper\": [\n");
    for (size_t prefix = 0; prefix < prefixes; ++prefix) {
        std::fprintf(out, "    ");
        print_vector(out, reference.prefix_covariance[prefix].data(),
                     reference.prefix_covariance[prefix].size());
        std::fprintf(out, prefix + 1 < prefixes ? ",\n" : "\n");
    }
    std::fprintf(out, "  ],\n");

    std::fprintf(out, "  \"error_difference\": [\n");
    for (size_t prefix = 0; prefix < prefixes; ++prefix) {
        std::fprintf(out, "    ");
        print_vector(out, difference[prefix].data(), difference[prefix].size());
        std::fprintf(out, prefix + 1 < prefixes ? ",\n" : "\n");
    }
    std::fprintf(out, "  ],\n");

    std::fprintf(out, "  \"root_difference_coarse\": ");
    print_vector(out, root_coarse.data(), root_coarse.size());
    std::fprintf(out, ",\n");
    std::fprintf(out, "  \"endpoint_difference_coarse\": ");
    print_vector(out, endpoint_coarse.data(), endpoint_coarse.size());
    std::fprintf(out, ",\n");

    std::fprintf(out, "  \"applied_magnetic_corrections\": [\n");
    for (size_t i = 0; i < reference.applied_mag.size(); ++i) {
        const MagEvent& ev = reference.applied_mag[i];
        std::fprintf(out, "    {\"prefix\": %d, \"time_s\": %.17g, \"sensitivity_axis\": ",
                     ev.prefix, ev.time);
        print_vector(out, ev.sensitivity_axis.data(), ev.sensitivity_axis.size());
        std::fprintf(out, ", \"innovation_covariance\": ");
        print_vector(out, ev.innovation_covariance.data(), ev.innovation_covariance.size());
        std::fprintf(out, "}%s", i + 1 < reference.applied_mag.size() ? ",\n" : "\n");
    }
    std::fprintf(out, "  ]\n}\n");

    if (out != stdout) std::fclose(out);
    return 0;
}
