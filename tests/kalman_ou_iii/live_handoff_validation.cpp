// Truth-based validation of the Phase-3 Live-entrance certificate.
//
// The certificate constructs a bound on each block of the handoff error from a
// declared physical envelope and a handful of startup observables.  In
// simulation the true error is available, so the one property that matters can
// be checked directly rather than argued:
//
//     certified == true  =>  every true handoff error is inside the bound it
//                            was certified against.
//
// A single violation is a false certification and fails this test.  The
// converse -- that a rejected handoff really was outside some bound -- is not
// required and is not checked: the certificate is a sufficient condition, and
// it is allowed to be conservative.  How conservative is reported, because a
// certificate with no false certifications and no acceptances would be safe
// and useless, and only the two numbers together say which one this is.
//
// The seas here are synthesised at the eight committed reference (Hs, Tp)
// pairs rather than replayed from the committed records.  That is deliberate:
// what this test needs is truth for v, p, a_w, the body rate and the attitude
// at one specific sample, which a generator has exactly and a recorded file
// does not.  The record replay in kalman_ou_iii-sim remains the scored
// evidence for filter accuracy; this is evidence for the certificate.
#define EIGEN_NON_ARDUINO

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

namespace {

using Wrapper = SeaStateFusion_OU_III<TrackerType::KALMANF>;
using Vec3 = Eigen::Vector3f;
using Quat = Eigen::Quaternionf;

constexpr float kPi = 3.14159265358979f;
constexpr float kG  = 9.80665f;
constexpr float kDt = 0.005f;              // 200 Hz IMU
constexpr float kMagDt = 0.04f;            // 25 Hz magnetometer
constexpr float kRunSec = 200.0f;

// The eight committed reference seas, by significant height and peak period.
struct Sea {
    const char* name;
    float Hs;    // m
    float Tp;    // s
};

const std::array<Sea, 8> kSeas{{
    {"J0.27", 0.27f, 2.4f},
    {"J1.50", 1.50f, 5.0f},
    {"J4.00", 4.00f, 7.4f},
    {"J8.50", 8.50f, 8.7f},
    {"P0.27", 0.27f, 2.5f},
    {"P1.50", 1.50f, 5.2f},
    {"P4.00", 4.00f, 7.2f},
    {"P8.50", 8.50f, 8.5f},
}};

// A wave orbit and the platform attitude riding it, both analytic, so every
// truth channel is exact at every sample rather than differenced.
//
// Elevation is a small sum of components around the peak period; the orbital
// horizontal motion is taken in the deep-water limit, where the horizontal and
// vertical displacements of a component have equal amplitude and are in
// quadrature.  Roll and pitch follow the surface slope, which is what a hull
// short against the wavelength does.
class TruthSea {
public:
    TruthSea(const Sea& sea, unsigned seed, float phase_offset)
        : sea_(sea)
    {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<float> uni(0.0f, 2.0f * kPi);
        const float wp = 2.0f * kPi / sea.Tp;
        // Amplitudes summing to the right significant height: Hs = 4 sigma.
        float acc = 0.0f;
        for (int i = 0; i < kNComp; ++i) {
            const float r = 0.7f + 0.15f * float(i);
            w_[size_t(i)] = wp * r;
            const float shape = std::exp(-6.0f * (r - 1.0f) * (r - 1.0f));
            a_[size_t(i)] = shape;
            acc += shape * shape;
            phi_[size_t(i)] = uni(rng) + phase_offset;
        }
        const float sigma_target = sea.Hs / 4.0f;
        const float scale = sigma_target / std::sqrt(0.5f * acc);
        for (int i = 0; i < kNComp; ++i) a_[size_t(i)] *= scale;
        dir_ = uni(rng);
    }

    void setSeaChange(float t_change, float hs_scale) {
        t_change_ = t_change;
        hs_scale_ = hs_scale;
    }

    // World NED position, velocity and acceleration of the platform.
    void kinematics(float t, Vec3& p, Vec3& v, Vec3& a) const {
        const float s = amplitudeScale(t);
        float eta = 0, deta = 0, d2eta = 0;      // up-positive elevation
        float hx = 0, dhx = 0, d2hx = 0;         // along-crest horizontal
        for (int i = 0; i < kNComp; ++i) {
            const float w = w_[size_t(i)];
            const float A = a_[size_t(i)] * s;
            const float ph = w * t + phi_[size_t(i)];
            eta   += A * std::sin(ph);
            deta  += A * w * std::cos(ph);
            d2eta -= A * w * w * std::sin(ph);
            hx    -= A * std::cos(ph);
            dhx   += A * w * std::sin(ph);
            d2hx  += A * w * w * std::cos(ph);
        }
        const float cd = std::cos(dir_), sd = std::sin(dir_);
        p = Vec3(hx * cd, hx * sd, -eta);
        v = Vec3(dhx * cd, dhx * sd, -deta);
        a = Vec3(d2hx * cd, d2hx * sd, -d2eta);
    }

    // Truth attitude and body rate.  Euler angles are analytic, so the body
    // rate follows from the Euler-rate transformation exactly.
    void attitude(float t, Quat& q_bw, Vec3& omega_body) const {
        float r = roll0_, pch = pitch0_, yaw = yaw0_;
        float dr = 0, dp = 0, dy = 0;
        const float s = amplitudeScale(t);
        for (int i = 0; i < kNComp; ++i) {
            const float w = w_[size_t(i)];
            const float slope = a_[size_t(i)] * s * w * w / kG;  // surface slope
            const float ph = w * t + phi_[size_t(i)];
            r   += slope * std::cos(ph) * std::sin(dir_);
            pch += slope * std::cos(ph) * std::cos(dir_);
            dr  -= slope * w * std::sin(ph) * std::sin(dir_);
            dp  -= slope * w * std::sin(ph) * std::cos(dir_);
            yaw += 0.05f * slope * std::sin(ph);
            dy  += 0.05f * slope * w * std::cos(ph);
        }
        q_bw = Quat(Eigen::AngleAxisf(yaw, Vec3::UnitZ())) *
               Quat(Eigen::AngleAxisf(pch, Vec3::UnitY())) *
               Quat(Eigen::AngleAxisf(r,   Vec3::UnitX()));
        q_bw.normalize();

        // Body rates from ZYX Euler rates.
        const float sr = std::sin(r), cr = std::cos(r);
        const float sp = std::sin(pch), cp = std::cos(pch);
        omega_body = Vec3(dr - dy * sp,
                          dp * cr + dy * sr * cp,
                          -dp * sr + dy * cr * cp);
    }

    void setInitialAttitude(float roll, float pitch, float yaw) {
        roll0_ = roll; pitch0_ = pitch; yaw0_ = yaw;
    }

private:
    float amplitudeScale(float t) const {
        if (!(t_change_ > 0.0f) || t < t_change_) return 1.0f;
        const float ramp = std::min(1.0f, (t - t_change_) / 20.0f);
        return 1.0f + (hs_scale_ - 1.0f) * ramp;
    }

    static constexpr int kNComp = 7;
    Sea  sea_;
    std::array<float, size_t(kNComp)> w_{}, a_{}, phi_{};
    float dir_ = 0.0f;
    float roll0_ = 0.0f, pitch0_ = 0.0f, yaw0_ = 0.0f;
    float t_change_ = -1.0f;
    float hs_scale_ = 1.0f;
};

struct Scenario {
    std::string label;
    size_t sea_index = 0;
    unsigned seed = 1;
    float phase = 0.0f;
    Vec3 gyro_bias = Vec3::Zero();     // rad/s, true
    Vec3 accel_bias = Vec3::Zero();    // m/s^2, true
    float mag_ref_tilt_rad = 0.0f;     // world-reference direction error
    float init_roll = 0.0f, init_pitch = 0.0f, init_yaw = 0.0f;
    float mag_delay_sec = 0.0f;
    float sea_change_at = -1.0f;
    float sea_change_scale = 1.0f;
    bool  with_mag = true;
};

struct Outcome {
    std::string label;
    bool handed_off = false;
    float handoff_time = std::numeric_limits<float>::quiet_NaN();
    bool certified = false;
    std::string reason;

    // true errors at the handoff sample
    float e_theta = 0, e_bg = 0, e_S = 0, e_aw = 0, e_ba = 0, e_v = 0, e_p = 0;
    float e_tilt = 0, e_yaw = 0;   // the two halves of the attitude error
    float b_tilt = 0, b_yaw = 0;   // and the two assumptions they answer to
    // the bounds those were certified against
    float b_theta = 0, b_bg = 0, b_S = 0, b_aw = 0, b_ba = 0, b_v = 0, b_p = 0;

    float margin = std::numeric_limits<float>::quiet_NaN();
    float sensitive_metric = std::numeric_limits<float>::quiet_NaN();
    float kinematic_metric = std::numeric_limits<float>::quiet_NaN();
    // The largest small-gain slope that would still have certified this
    // handoff.  The budget (1 - c r) r peaks at 1/(4c), so the inequality
    // closes exactly when c_eff < 1 / (4 * basin_lhs).  Reporting it turns a
    // rejection into a number: it says how much of the gap is the handoff and
    // how much is the interval.
    float required_c_eff = std::numeric_limits<float>::quiet_NaN();

    bool violated = false;
    bool bounds_checked = false;
    float worst_ratio = std::numeric_limits<float>::infinity();
    std::string violated_block;
};

// Tilt and heading parts of the error between two BODY->NED attitudes.  The
// certificate bounds them from two different assumptions, so a validation that
// only reported their sum could not say which one was wrong.
void tiltYawError(const Quat& est, const Quat& truth, float& tilt, float& yaw) {
    const Eigen::Matrix3f Re = est.toRotationMatrix();
    const Eigen::Matrix3f Rt = truth.toRotationMatrix();
    // Third row of R(body->NED) is world down expressed in body coordinates.
    const Vec3 de = Re.transpose() * Vec3::UnitZ();
    const Vec3 dt = Rt.transpose() * Vec3::UnitZ();
    float c = de.dot(dt);
    c = std::min(1.0f, std::max(-1.0f, c));
    tilt = std::acos(c);
    const float ye = std::atan2(Re(1,0), Re(0,0));
    const float yt = std::atan2(Rt(1,0), Rt(0,0));
    float d = ye - yt;
    while (d > kPi) d -= 2.0f * kPi;
    while (d <= -kPi) d += 2.0f * kPi;
    yaw = std::fabs(d);
}

float angleBetween(const Quat& a, const Quat& b) {
    Quat d = a.conjugate() * b;
    d.normalize();
    float w = std::fabs(d.w());
    w = std::min(1.0f, std::max(-1.0f, w));
    return 2.0f * std::acos(w);
}

Outcome runScenario(const Scenario& sc) {
    Outcome out;
    out.label = sc.label;

    TruthSea sea(kSeas[sc.sea_index], sc.seed, sc.phase);
    sea.setInitialAttitude(sc.init_roll, sc.init_pitch, sc.init_yaw);
    if (sc.sea_change_at > 0.0f) sea.setSeaChange(sc.sea_change_at, sc.sea_change_scale);

    // World magnetic field, and the field the platform actually sees.  The
    // reference the filter learns comes out of the measurements, so a tilt
    // applied to the true field is what a real reference error looks like.
    // Zero declination on purpose.  The filter gauges heading to the learned
    // magnetic reference, while the truth yaw here is the generator's; a field
    // with an east component would show up as a constant heading error that is
    // the declination and not the gauge.  Declination is a chart correction,
    // not something this test is about.
    const Vec3 mag_world(20.5f, 0.0f, 44.0f);
    const Vec3 mag_world_seen =
        Quat(Eigen::AngleAxisf(sc.mag_ref_tilt_rad, Vec3::UnitZ())) * mag_world;

    std::mt19937 rng(sc.seed * 7919u + 13u);
    std::normal_distribution<float> n_acc(0.0f, 0.0294f);
    std::normal_distribution<float> n_gyr(0.0f, 0.00157f);
    std::normal_distribution<float> n_mag(0.0f, 0.36f);

    Wrapper filter;
    Wrapper::Config cfg;
    cfg.with_mag = sc.with_mag;
    cfg.mag_delay_sec = sc.mag_delay_sec;
    // The declared envelope the certificate is conditional on.
    //
    // It is a declaration, not a fit, and a deployment has to declare one that
    // covers the conditions it will actually see.  This test draws its biases,
    // reference errors and sea states from a family it chose, so it declares
    // an envelope that covers that family with headroom -- which is exactly
    // what a deployment does, and the reason the numbers below are wider than
    // the library defaults rather than equal to them.  Declaring an envelope
    // the scenarios then leave is not a stricter test; it is applying the
    // theorem outside its own hypotheses.
    cfg.live_envelope.gravity_direction_rad   = 0.12f;
    cfg.live_envelope.mag_reference_rel_error = 0.25f;
    cfg.live_envelope.gyro_bias_axial_rad_s   = 3.0e-4f;
    cfg.live_envelope.gyro_bias_perp_rad_s    = 3.0e-4f;
    cfg.live_envelope.accel_bias_ms2          = 0.12f;
    cfg.live_envelope.accel_meas_ms2          = 0.10f;
    cfg.live_envelope.translational_velocity_ms = 4.0f;
    cfg.live_envelope.translational_position_m  = 8.0f;

    cfg.sigma_a = Vec3::Constant(0.0294f);
    cfg.sigma_g = Vec3::Constant(0.00157f);
    cfg.sigma_m = Vec3::Constant(0.36f);
    if (std::getenv("NO_BG_SEED")) cfg.seed_gyro_bias_at_handoff = false;
    if (std::getenv("NO_AW_SEED")) cfg.seed_world_accel_at_handoff = false;
    if (std::getenv("NO_ENV_COV")) cfg.seed_covariance_from_envelope = false;
    filter.begin(cfg);

    float t = 0.0f;
    float t_next_mag = 0.0f;
    bool was_live = false;

    const int steps = int(kRunSec / kDt);
    for (int k = 0; k < steps; ++k) {
        Quat q_bw;
        Vec3 omega;
        sea.attitude(t, q_bw, omega);
        Vec3 p_t, v_t, a_t;
        sea.kinematics(t, p_t, v_t, a_t);

        const Vec3 g_ned(0.0f, 0.0f, kG);
        const Vec3 f_body = q_bw.conjugate() * (a_t - g_ned);

        Vec3 gyro_meas = omega + sc.gyro_bias;
        Vec3 acc_meas  = f_body + sc.accel_bias;
        for (int i = 0; i < 3; ++i) {
            gyro_meas[i] += n_gyr(rng);
            acc_meas[i]  += n_acc(rng);
        }

        filter.update(kDt, gyro_meas, acc_meas, 35.0f);

        if (sc.with_mag && t >= t_next_mag) {
            Vec3 mag_meas = q_bw.conjugate() * mag_world_seen;
            for (int i = 0; i < 3; ++i) mag_meas[i] += n_mag(rng);
            filter.updateMag(mag_meas);
            t_next_mag += kMagDt;
        }

        if (!was_live && filter.isLive()) {
            was_live = true;
            out.handed_off = true;
            out.handoff_time = t;

            const auto& c = filter.liveEntranceCertificate();
            out.certified = filter.isLiveCertified();
            out.reason = filter.liveCertificationReason();
            out.margin = c.margin;
            out.sensitive_metric = c.sensitive_metric_norm;
            out.kinematic_metric = c.kinematic_metric_norm;
            if (std::isfinite(c.basin_lhs) && c.basin_lhs > 0.0f) {
                out.required_c_eff = 1.0f / (4.0f * c.basin_lhs);
            }

            const auto& m = filter.raw().mekf();
            out.e_theta = angleBetween(m.quaternion_boat(), q_bw);
            tiltYawError(m.quaternion_boat(), q_bw, out.e_tilt, out.e_yaw);
            out.e_bg = (m.gyroscope_bias() - sc.gyro_bias).norm();
            // S is exactly zero on both sides: the filter resets its
            // integration epoch here, and the truth integral measured from the
            // same epoch is zero at the epoch by definition.
            out.e_S = m.get_integral_displacement().norm();
            out.e_aw = (m.get_world_accel() - a_t).norm();
            out.e_ba = (m.get_acc_bias_at_temperature(35.0f) - sc.accel_bias).norm();
            out.e_v = (m.get_velocity() - v_t).norm();
            out.e_p = (m.get_position() - p_t).norm();

            out.b_theta = c.attitude.value;
            // The attitude bound is a sum of two separate assumptions, and the
            // sum can hold while one of them does not.  Check them apart, or
            // the derivation is not being tested at all.
            out.b_tilt = std::asin(cfg.mag_gravity_align_max_sin) +
                         c.eta_gravity_rad;
            out.b_yaw = c.eta_heading_rad;
            out.b_bg = c.gyro_bias.value;
            out.b_S = c.integral_S.value;
            out.b_aw = c.world_accel.value;
            out.b_ba = c.accel_bias.value;
            out.b_v = c.velocity.value;
            out.b_p = c.position.value;

            {
                struct Check { const char* n; float e, b; };
                const Check checks[] = {
                    {"theta", out.e_theta, out.b_theta},
                    {"tilt",  out.e_tilt,  out.b_tilt},
                    {"yaw",   out.e_yaw,   out.b_yaw},
                    {"b_g",   out.e_bg,    out.b_bg},
                    {"S",     out.e_S,     out.b_S},
                    {"a_w",   out.e_aw,    out.b_aw},
                    {"b_a",   out.e_ba,    out.b_ba},
                    {"v",     out.e_v,     out.b_v},
                    {"p",     out.e_p,     out.b_p},
                };
                // Checked on every handoff, not only certified ones.  If the
                // certificate never certifies, checking only certified
                // handoffs is vacuous, and the property that actually has to
                // hold -- that a constructed bound really bounds the truth --
                // is testable regardless of whether the basin inequality
                // closed.  A block with no finite bound (no magnetic gauge,
                // say) has nothing to check.
                for (const Check& ch : checks) {
                    if (!std::isfinite(ch.b)) continue;
                    // A tolerance of exactly zero on the S block: it is an
                    // identity, not a bound, and a nonzero value there means
                    // the epoch reset did not happen.
                    const float tol = (ch.b == 0.0f) ? 1e-6f : 0.0f;
                    if (!(ch.e <= ch.b + tol)) {
                        out.violated = true;
                        out.violated_block = ch.n;
                        break;
                    }
                    if (ch.e > 0.0f) {
                        out.worst_ratio = std::min(out.worst_ratio, ch.b / ch.e);
                    }
                    out.bounds_checked = true;
                }
            }
        }

        t += kDt;
    }

    return out;
}

std::vector<Scenario> buildScenarios() {
    std::vector<Scenario> v;
    std::mt19937 rng(20260821u);
    std::uniform_real_distribution<float> ph(0.0f, 2.0f * kPi);

    // Bias draws are uniform inside a ball whose radius is strictly below the
    // declared envelope, rather than Gaussian.  A Gaussian draw leaves any
    // finite envelope with positive probability, so a test that checked
    // envelope-conditional bounds against Gaussian draws would fail
    // occasionally for a reason that is not a defect -- it would be the
    // theorem correctly declining to cover a case outside its hypotheses.
    // Bounded draws keep the CI signal about the bound and not about the tail.
    std::normal_distribution<float> dir(0.0f, 1.0f);
    std::uniform_real_distribution<float> unit(0.0f, 1.0f);
    auto ball = [&](float radius) {
        Vec3 d(dir(rng), dir(rng), dir(rng));
        const float n = d.norm();
        if (!(n > 1e-9f)) return Vec3::Zero().eval();
        return (d / n * (radius * std::cbrt(unit(rng)))).eval();
    };
    constexpr float kGyroBiasBall  = 4.0e-4f;   // envelope declares 4.24e-4
    constexpr float kAccelBiasBall = 0.10f;     // envelope declares 0.12

    // Base sweep: eight seas x three seeds x three wave phases, with gyro- and
    // accelerometer-bias draws inside the declared envelope.
    for (size_t s = 0; s < kSeas.size(); ++s) {
        for (unsigned seed = 1; seed <= 3; ++seed) {
            for (int phase = 0; phase < 3; ++phase) {
                Scenario sc;
                char buf[96];
                std::snprintf(buf, sizeof buf, "%s/seed%u/ph%d",
                              kSeas[s].name, seed, phase);
                sc.label = buf;
                sc.sea_index = s;
                sc.seed = seed * 31u + unsigned(phase);
                sc.phase = ph(rng);
                sc.gyro_bias = ball(kGyroBiasBall);
                sc.accel_bias = ball(kAccelBiasBall);
                sc.mag_ref_tilt_rad = 0.02f * (float(phase) - 1.0f);
                v.push_back(sc);
            }
        }
    }

    // Large initial roll/pitch/yaw: the proxy has to capture from well away
    // from level before any of this means anything.
    const float big[][3] = {
        {0.7f, 0.0f, 1.0f}, {0.0f, 0.7f, -2.0f}, {-0.6f, 0.5f, 3.0f},
        {1.2f, 0.0f, 0.5f}, {0.0f, -1.1f, -1.5f},
    };
    for (size_t i = 0; i < sizeof(big) / sizeof(big[0]); ++i) {
        Scenario sc;
        char buf[96];
        std::snprintf(buf, sizeof buf, "large-attitude-%zu", i);
        sc.label = buf;
        sc.sea_index = i % kSeas.size();
        sc.seed = 900u + unsigned(i);
        sc.phase = ph(rng);
        sc.init_roll = big[i][0];
        sc.init_pitch = big[i][1];
        sc.init_yaw = big[i][2];
        sc.gyro_bias = ball(kGyroBiasBall);
        sc.accel_bias = ball(kAccelBiasBall);
        v.push_back(sc);
    }

    // Near, but not on, the antipodal set the accel-corrected proxy is not
    // attracted to.  Inverted by 170 deg rather than 180.
    for (int i = 0; i < 3; ++i) {
        Scenario sc;
        char buf[96];
        std::snprintf(buf, sizeof buf, "near-antipodal-%d", i);
        sc.label = buf;
        sc.sea_index = size_t(i);
        sc.seed = 1200u + unsigned(i);
        sc.phase = ph(rng);
        sc.init_roll = 2.967f;                 // 170 deg
        sc.init_pitch = 0.05f * float(i);
        sc.init_yaw = 0.4f * float(i);
        v.push_back(sc);
    }

    // Delayed magnetometer availability.
    for (int i = 0; i < 3; ++i) {
        Scenario sc;
        char buf[96];
        std::snprintf(buf, sizeof buf, "mag-delay-%d", i);
        sc.label = buf;
        sc.sea_index = size_t(i + 1);
        sc.seed = 1500u + unsigned(i);
        sc.phase = ph(rng);
        sc.mag_delay_sec = 20.0f + 20.0f * float(i);
        v.push_back(sc);
    }

    // Sea state changing during startup.
    for (int i = 0; i < 3; ++i) {
        Scenario sc;
        char buf[96];
        std::snprintf(buf, sizeof buf, "sea-change-%d", i);
        sc.label = buf;
        sc.sea_index = size_t(i);
        sc.seed = 1800u + unsigned(i);
        sc.phase = ph(rng);
        sc.sea_change_at = 25.0f;
        sc.sea_change_scale = 3.0f + float(i);
        v.push_back(sc);
    }

    // No magnetometer at all: the handoff has to fall through to the timeout,
    // and the certificate has to refuse it for want of a heading bound.
    for (int i = 0; i < 2; ++i) {
        Scenario sc;
        char buf[96];
        std::snprintf(buf, sizeof buf, "timeout-nomag-%d", i);
        sc.label = buf;
        sc.sea_index = size_t(i * 3);
        sc.seed = 2100u + unsigned(i);
        sc.phase = ph(rng);
        sc.with_mag = false;
        v.push_back(sc);
    }

    // Magnetic reference error at the edge of the declared envelope.
    for (int i = 0; i < 3; ++i) {
        Scenario sc;
        char buf[96];
        std::snprintf(buf, sizeof buf, "mag-ref-error-%d", i);
        sc.label = buf;
        sc.sea_index = size_t(i + 4);
        sc.seed = 2400u + unsigned(i);
        sc.phase = ph(rng);
        sc.mag_ref_tilt_rad = 0.02f + 0.015f * float(i);
        v.push_back(sc);
    }

    return v;
}

}  // namespace

int main(int argc, char** argv) {
    bool verbose = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-v") verbose = true;
    }

    const std::vector<Scenario> scenarios = buildScenarios();
    std::vector<Outcome> outcomes;
    outcomes.reserve(scenarios.size());
    for (const Scenario& sc : scenarios) {
        outcomes.push_back(runScenario(sc));
    }

    int n_total = 0, n_handed = 0, n_certified = 0, n_false = 0, n_uncertified = 0;
    int n_violations = 0, n_bounds_checked = 0;
    float min_margin = std::numeric_limits<float>::infinity();
    float max_handoff_certified = 0.0f;
    float max_handoff_uncertified = 0.0f;
    float worst_ratio_theta = std::numeric_limits<float>::infinity();
    float worst_ratio_any = std::numeric_limits<float>::infinity();
    float best_required_c_eff = 0.0f;
    float worst_required_c_eff = std::numeric_limits<float>::infinity();
    float worst_tilt = 0.0f, worst_yaw = 0.0f, worst_S = 0.0f;

    std::cout << std::setprecision(6);
    for (const Outcome& o : outcomes) {
        ++n_total;
        if (o.handed_off) ++n_handed;
        if (o.certified) {
            ++n_certified;
            min_margin = std::min(min_margin, o.margin);
            max_handoff_certified = std::max(max_handoff_certified, o.handoff_time);
            if (o.e_theta > 0.0f) {
                worst_ratio_theta = std::min(worst_ratio_theta, o.b_theta / o.e_theta);
            }
        } else if (o.handed_off) {
            ++n_uncertified;
            max_handoff_uncertified = std::max(max_handoff_uncertified, o.handoff_time);
        }
        if (std::isfinite(o.required_c_eff)) {
            best_required_c_eff = std::max(best_required_c_eff, o.required_c_eff);
            worst_required_c_eff = std::min(worst_required_c_eff, o.required_c_eff);
        }
        worst_tilt = std::max(worst_tilt, o.e_tilt);
        worst_yaw = std::max(worst_yaw, o.e_yaw);
        worst_S = std::max(worst_S, o.e_S);
        if (o.bounds_checked) {
            ++n_bounds_checked;
            worst_ratio_any = std::min(worst_ratio_any, o.worst_ratio);
        }
        if (o.violated) {
            ++n_violations;
            if (o.certified) ++n_false;
            std::cerr << "BOUND VIOLATION " << o.label
                      << " block=" << o.violated_block << '\n';
        }
        if (verbose || o.violated) {
            std::cout << "  " << o.label
                      << " live=" << (o.handed_off ? 1 : 0)
                      << " t=" << o.handoff_time
                      << " cert=" << (o.certified ? 1 : 0)
                      << " reason=" << o.reason
                      << " margin=" << o.margin
                      << " |zxi|=" << o.sensitive_metric
                      << " |zl|=" << o.kinematic_metric
                      << " theta " << o.e_theta << "/" << o.b_theta
                      << " (tilt " << o.e_tilt << " yaw " << o.e_yaw << ")"
                      << " bg " << o.e_bg << "/" << o.b_bg
                      << " aw " << o.e_aw << "/" << o.b_aw
                      << " ba " << o.e_ba << "/" << o.b_ba
                      << " v " << o.e_v << "/" << o.b_v
                      << " p " << o.e_p << "/" << o.b_p
                      << '\n';
        }
    }

    std::cout << "HANDOFF_VALIDATION scenarios=" << n_total
              << " handed_off=" << n_handed
              << " certified=" << n_certified
              << " uncertified=" << n_uncertified
              << " false_certifications=" << n_false
              << " bounds_checked=" << n_bounds_checked
              << " bound_violations=" << n_violations
              << " acceptance_rate="
              << (n_handed ? double(n_certified) / double(n_handed) : 0.0)
              << " min_margin=" << min_margin
              << " worst_conservatism_theta=" << worst_ratio_theta
              << " worst_conservatism_any=" << worst_ratio_any
              << " max_handoff_sec_certified=" << max_handoff_certified
              << " max_handoff_sec_uncertified=" << max_handoff_uncertified
              << " required_c_eff_best=" << best_required_c_eff
              << " required_c_eff_worst=" << worst_required_c_eff
              << " worst_true_tilt_rad=" << worst_tilt
              << " worst_true_yaw_rad=" << worst_yaw
              << " worst_true_S=" << worst_S
              << '\n';

    // The safety property.  Nothing else in this program is allowed to fail
    // the build: an acceptance rate of zero is a reportable result, a false
    // certification is a broken theorem.
    if (n_false != 0) {
        std::cerr << "FAIL: " << n_false
                  << " certified handoff(s) were outside the bound they were "
                     "certified against\n";
        return 1;
    }

    // The unconditional form of the same property.  A constructed bound has to
    // bound the truth whether or not the basin inequality closed; this is what
    // keeps the safety statement from being vacuous while the acceptance rate
    // is zero.
    if (n_violations != 0) {
        std::cerr << "FAIL: " << n_violations
                  << " handoff(s) left a constructed bound while inside the "
                     "declared envelope\n";
        return 1;
    }
    if (n_bounds_checked == 0) {
        std::cerr << "FAIL: no handoff had its constructed bounds checked\n";
        return 1;
    }

    // Every scenario must still reach Live.  The certificate is not allowed to
    // hold the filter back, and a bootstrap that never hands over would make
    // the numbers above vacuous.
    if (n_handed != n_total) {
        std::cerr << "FAIL: " << (n_total - n_handed)
                  << " scenario(s) never reached Live; the certificate must not"
                     " withhold the handoff\n";
        return 1;
    }

    // The integration-epoch reset is an identity, not a bound, and it applies
    // to every handoff rather than only to certified ones.
    if (!(worst_S <= 1e-6f)) {
        std::cerr << "FAIL: the integral state was not exactly zero at handoff"
                     " (worst " << worst_S << " m s); the epoch reset did not"
                     " happen\n";
        return 1;
    }

    // A handoff with no magnetometer has no heading bound, so it must not be
    // certified.  This is the one direction of the converse that is a
    // correctness property rather than conservatism.
    for (const Outcome& o : outcomes) {
        if (o.label.rfind("timeout-nomag", 0) == 0 && o.certified) {
            std::cerr << "FAIL: " << o.label
                      << " was certified without a magnetic heading gauge\n";
            return 1;
        }
    }

    return 0;
}
