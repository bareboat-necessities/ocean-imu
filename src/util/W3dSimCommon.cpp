/*
  Copyright 2025-2026, Mikhail Grushinskiy
*/

extern const float g_std = 9.80665f;

#ifdef ARDUINO
#else

#include "util/W3dSimCommon.h"

#include <cerrno>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <sstream>
#include <stdexcept>

static bool g_w3d_any_gate_failed = false;

ImuNoiseModel make_imu_noise_model(float sigma_white,
                                   float bias_half_range,
                                   float sigma_bias_rw,
                                   unsigned seed)
{
    ImuNoiseModel m;
    m.rng = std::mt19937(seed);
    m.w = std::normal_distribution<float>(0.0f, sigma_white);
    m.n01 = std::normal_distribution<float>(0.0f, 1.0f);
    m.bias0.setZero();
    m.bias_rw.setZero();
    m.sigma_bias_rw = sigma_bias_rw;

    std::uniform_real_distribution<float> ub(-bias_half_range, bias_half_range);
    m.bias0 = Vector3f(ub(m.rng), ub(m.rng), ub(m.rng));
    return m;
}

ImuNoiseModel make_imu_noise_model(float sigma_white,
                                   float bias_half_range,
                                   float sigma_bias_rw,
                                   unsigned noise_seed,
                                   unsigned initialization_seed)
{
    // Preserve the original draw order for the historical single-seed path.
    if (noise_seed == initialization_seed) {
        return make_imu_noise_model(
            sigma_white, bias_half_range, sigma_bias_rw, noise_seed);
    }

    ImuNoiseModel m;
    m.rng = std::mt19937(noise_seed);
    m.w = std::normal_distribution<float>(0.0f, sigma_white);
    m.n01 = std::normal_distribution<float>(0.0f, 1.0f);
    m.bias_rw.setZero();
    m.sigma_bias_rw = sigma_bias_rw;

    std::mt19937 initialization_rng(initialization_seed);
    std::uniform_real_distribution<float> ub(-bias_half_range, bias_half_range);
    m.bias0 = Vector3f(
        ub(initialization_rng), ub(initialization_rng), ub(initialization_rng));
    return m;
}

Vector3f apply_imu_noise(const Vector3f& truth, ImuNoiseModel& m, float dt)
{
    if (m.sigma_bias_rw > 0.0f) {
        const float s = m.sigma_bias_rw * std::sqrt(dt);
        m.bias_rw += Vector3f(s * m.n01(m.rng), s * m.n01(m.rng), s * m.n01(m.rng));
    }
    Vector3f white(m.w(m.rng), m.w(m.rng), m.w(m.rng));
    return truth + (m.bias0 + m.bias_rw) + white;
}

MagNoiseModel make_mag_noise_model(float sigma_white_uT,
                                   float bias_residual_range_uT,
                                   float sigma_bias_rw_uT_sqrt_s,
                                   float scale_err_max,
                                   float cross_axis_max,
                                   float misalign_deg_max,
                                   unsigned seed)
{
    MagNoiseModel m;
    m.rng = std::mt19937(seed);
    m.w_uT = std::normal_distribution<float>(0.0f, sigma_white_uT);
    m.n01 = std::normal_distribution<float>(0.0f, 1.0f);

    std::uniform_real_distribution<float> ub(-bias_residual_range_uT, bias_residual_range_uT);
    m.bias0_uT = Vector3f(ub(m.rng), ub(m.rng), ub(m.rng));
    m.bias_rw_uT.setZero();
    m.sigma_bias_rw_uT_sqrt_s = sigma_bias_rw_uT_sqrt_s;

    std::uniform_real_distribution<float> us(1.0f - scale_err_max, 1.0f + scale_err_max);
    std::uniform_real_distribution<float> uc(-cross_axis_max, cross_axis_max);

    Eigen::Matrix3f A = Eigen::Matrix3f::Identity();
    A(0, 0) = us(m.rng);
    A(1, 1) = us(m.rng);
    A(2, 2) = us(m.rng);

    float a01 = uc(m.rng), a02 = uc(m.rng), a12 = uc(m.rng);
    A(0, 1) = A(1, 0) = a01;
    A(0, 2) = A(2, 0) = a02;
    A(1, 2) = A(2, 1) = a12;

    auto deg2rad = [](float d) { return d * float(std::numbers::pi_v<float> / 180.0); };
    std::uniform_real_distribution<float> ua(-misalign_deg_max, misalign_deg_max);
    float rx = deg2rad(ua(m.rng));
    float ry = deg2rad(ua(m.rng));
    float rz = deg2rad(ua(m.rng));

    auto Rx = [&](float a) {
        Eigen::Matrix3f R;
        float c = std::cos(a), s = std::sin(a);
        R << 1, 0, 0,
             0, c, -s,
             0, s, c;
        return R;
    };
    auto Ry = [&](float a) {
        Eigen::Matrix3f R;
        float c = std::cos(a), s = std::sin(a);
        R << c, 0, s,
             0, 1, 0,
            -s, 0, c;
        return R;
    };
    auto Rz = [&](float a) {
        Eigen::Matrix3f R;
        float c = std::cos(a), s = std::sin(a);
        R << c, -s, 0,
             s, c, 0,
             0, 0, 1;
        return R;
    };

    Eigen::Matrix3f R = Rz(rz) * Ry(ry) * Rx(rx);
    m.Mis = R * A;
    return m;
}

MagNoiseModel make_mag_noise_model(float sigma_white_uT,
                                   float bias_residual_range_uT,
                                   float sigma_bias_rw_uT_sqrt_s,
                                   float scale_err_max,
                                   float cross_axis_max,
                                   float misalign_deg_max,
                                   unsigned noise_seed,
                                   unsigned initialization_seed)
{
    // Preserve the original draw order for the historical single-seed path.
    if (noise_seed == initialization_seed) {
        return make_mag_noise_model(
            sigma_white_uT,
            bias_residual_range_uT,
            sigma_bias_rw_uT_sqrt_s,
            scale_err_max,
            cross_axis_max,
            misalign_deg_max,
            noise_seed);
    }

    MagNoiseModel m;
    m.rng = std::mt19937(noise_seed);
    m.w_uT = std::normal_distribution<float>(0.0f, sigma_white_uT);
    m.n01 = std::normal_distribution<float>(0.0f, 1.0f);
    m.bias_rw_uT.setZero();
    m.sigma_bias_rw_uT_sqrt_s = sigma_bias_rw_uT_sqrt_s;

    std::mt19937 initialization_rng(initialization_seed);
    std::uniform_real_distribution<float> ub(
        -bias_residual_range_uT, bias_residual_range_uT);
    m.bias0_uT = Vector3f(
        ub(initialization_rng), ub(initialization_rng), ub(initialization_rng));

    std::uniform_real_distribution<float> us(
        1.0f - scale_err_max, 1.0f + scale_err_max);
    std::uniform_real_distribution<float> uc(-cross_axis_max, cross_axis_max);

    Eigen::Matrix3f A = Eigen::Matrix3f::Identity();
    A(0, 0) = us(initialization_rng);
    A(1, 1) = us(initialization_rng);
    A(2, 2) = us(initialization_rng);

    const float a01 = uc(initialization_rng);
    const float a02 = uc(initialization_rng);
    const float a12 = uc(initialization_rng);
    A(0, 1) = A(1, 0) = a01;
    A(0, 2) = A(2, 0) = a02;
    A(1, 2) = A(2, 1) = a12;

    auto deg2rad = [](float d) {
        return d * float(std::numbers::pi_v<float> / 180.0);
    };
    std::uniform_real_distribution<float> ua(
        -misalign_deg_max, misalign_deg_max);
    const float rx = deg2rad(ua(initialization_rng));
    const float ry = deg2rad(ua(initialization_rng));
    const float rz = deg2rad(ua(initialization_rng));

    auto Rx = [](float a) {
        Eigen::Matrix3f R;
        const float c = std::cos(a), s = std::sin(a);
        R << 1, 0, 0,
             0, c, -s,
             0, s, c;
        return R;
    };
    auto Ry = [](float a) {
        Eigen::Matrix3f R;
        const float c = std::cos(a), s = std::sin(a);
        R << c, 0, s,
             0, 1, 0,
            -s, 0, c;
        return R;
    };
    auto Rz = [](float a) {
        Eigen::Matrix3f R;
        const float c = std::cos(a), s = std::sin(a);
        R << c, -s, 0,
             s, c, 0,
             0, 0, 1;
        return R;
    };

    m.Mis = Rz(rz) * Ry(ry) * Rx(rx) * A;
    return m;
}

namespace {

constexpr float kTwoPi = 6.28318530717958647692f;

// Crank-order weights for a four-stroke.  The firing order carries the
// dominant gas-torque line; the non-firing half and whole orders come from
// reciprocating unbalance and cylinder-to-cylinder variability, which an
// inline engine with an odd cylinder count never fully cancels.
struct EngineOrderWeight {
    float order;
    float weight;
};

constexpr EngineOrderWeight kFiringHarmonicWeights[] = {
    {1.0f, 1.00f},
    {2.0f, 0.55f},
    {3.0f, 0.30f},
    {4.0f, 0.18f},
};

constexpr EngineOrderWeight kNonFiringOrderWeights[] = {
    {0.5f, 0.25f},
    {1.0f, 0.50f},
    {2.0f, 0.20f},
};

// Driveline weights, on shaft-rate orders.  The blade-rate line is the
// propeller's pressure pulse on the hull above the aperture; shaft order one
// is residual shaft unbalance and coupling misalignment.
constexpr float kDrivelineShaftWeight = 0.35f;
constexpr float kDrivelineBladeWeight = 0.45f;
constexpr float kDrivelineBladeSecondWeight = 0.18f;

// RMS split at the reference speed, relative to the engine harmonic RMS.
constexpr float kDrivelineRmsFraction = 0.45f;
constexpr float kBroadbandRmsFraction = 0.55f;

// Body-frame ZU anisotropy: x is athwartships, y is fore-and-aft, z is
// vertical.  The vertical path through the mounts and the hull bottom is the
// stiffest and carries the most; the fore-and-aft axis, aligned with the
// shaft, carries the least.
constexpr float kAxisGainX = 0.75f;
constexpr float kAxisGainY = 0.55f;
constexpr float kAxisGainZ = 1.00f;

// High-pass corner of the broadband structural floor.
constexpr float kBroadbandHighPassHz = 5.0f;

// Reference sensor bandwidth at which level_mps2 is stated, so the default
// configuration records exactly the configured broadband level and a narrower
// or wider anti-alias filter records proportionally less or more power.
constexpr float kReferenceBandwidthHz = 80.0f;

float engine_mount_transmissibility(float freq_hz, float mount_hz, float zeta)
{
    if (!(mount_hz > 0.0f)) return 1.0f;
    const float r = freq_hz / mount_hz;
    const float num = 1.0f + (2.0f * zeta * r) * (2.0f * zeta * r);
    const float one_minus = 1.0f - r * r;
    const float den = one_minus * one_minus + (2.0f * zeta * r) * (2.0f * zeta * r);
    if (!(den > 0.0f)) return 1.0f;
    return std::sqrt(num / den);
}

// Two-pole magnitude response of the sensor's anti-alias path.
float engine_antialias_gain(float freq_hz, float bandwidth_hz)
{
    if (!(bandwidth_hz > 0.0f)) return 1.0f;
    const float r = freq_hz / bandwidth_hz;
    const float r2 = r * r;
    return 1.0f / std::sqrt(1.0f + r2 * r2);
}

struct EngineHarmonicSpec {
    float freq_ratio = 0.0f;   // multiplier on crank frequency
    float weight = 0.0f;
    bool driveline = false;
};

std::vector<EngineHarmonicSpec> engine_harmonic_specs(const EngineVibrationConfig& cfg,
                                                      float firing_order)
{
    std::vector<EngineHarmonicSpec> specs;

    auto push_engine = [&specs](float order, float weight) {
        for (auto& spec : specs) {
            if (!spec.driveline && std::fabs(spec.freq_ratio - order) < 1e-4f) {
                spec.weight = std::max(spec.weight, weight);
                return;
            }
        }
        specs.push_back(EngineHarmonicSpec{order, weight, false});
    };

    for (const auto& entry : kFiringHarmonicWeights) {
        push_engine(entry.order * firing_order, entry.weight);
    }
    for (const auto& entry : kNonFiringOrderWeights) {
        push_engine(entry.order, entry.weight);
    }

    const float gear = (cfg.gear_ratio > 0.0f) ? cfg.gear_ratio : 1.0f;
    const float shaft_ratio = 1.0f / gear;
    const float blades = static_cast<float>(std::max(2, cfg.blades));
    specs.push_back(EngineHarmonicSpec{shaft_ratio, kDrivelineShaftWeight, true});
    specs.push_back(EngineHarmonicSpec{shaft_ratio * blades, kDrivelineBladeWeight, true});
    specs.push_back(
        EngineHarmonicSpec{shaft_ratio * 2.0f * blades, kDrivelineBladeSecondWeight, true});

    return specs;
}

// Unit-gain amplitude of one engine harmonic: inertial f^2 excitation through
// the flexible mount, referred to the firing line at the reference speed.
float engine_unit_amplitude(const EngineVibrationConfig& cfg, float freq_hz, float ref_firing_hz)
{
    if (!(ref_firing_hz > 0.0f)) return 0.0f;
    const float ratio = freq_hz / ref_firing_hz;
    return ratio * ratio *
        engine_mount_transmissibility(freq_hz, cfg.mount_hz, cfg.mount_zeta);
}

float engine_harmonic_rms(const EngineVibrationConfig& cfg,
                          const std::vector<EngineHarmonicSpec>& specs,
                          float crank_hz,
                          float ref_firing_hz,
                          bool driveline,
                          float driveline_speed_ratio)
{
    double sum = 0.0;
    for (const auto& spec : specs) {
        if (spec.driveline != driveline) continue;
        const float freq_hz = spec.freq_ratio * crank_hz;
        const float amp = driveline
            ? (spec.weight * driveline_speed_ratio * driveline_speed_ratio)
            : (spec.weight * engine_unit_amplitude(cfg, freq_hz, ref_firing_hz));
        sum += 0.5 * static_cast<double>(amp) * static_cast<double>(amp);
    }
    return static_cast<float>(std::sqrt(sum));
}

} // namespace

EngineVibrationModel make_engine_vibration_model(const EngineVibrationConfig& cfg, float dt)
{
    EngineVibrationModel m;
    m.cfg = cfg;
    m.rng = std::mt19937(cfg.seed);
    m.n01 = std::normal_distribution<float>(0.0f, 1.0f);

    if (!(cfg.rpm > 0.0f) || !(dt > 0.0f)) {
        return m;
    }

    const float cylinders = static_cast<float>(std::max(1, cfg.cylinders));
    m.firing_order = 0.5f * cylinders;
    m.crank_hz = cfg.rpm / 60.0f;
    const float gear = (cfg.gear_ratio > 0.0f) ? cfg.gear_ratio : 1.0f;
    m.shaft_hz = m.crank_hz / gear;

    const float reference_rpm = (cfg.reference_rpm > 0.0f) ? cfg.reference_rpm : cfg.rpm;
    const float reference_crank_hz = reference_rpm / 60.0f;
    const float reference_firing_hz = m.firing_order * reference_crank_hz;

    const auto specs = engine_harmonic_specs(cfg, m.firing_order);

    // One overall gain, calibrated so the hull broadband RMS is level_mps2 at
    // the reference speed.  Everything else follows from the physics above.
    const float engine_rms_ref =
        engine_harmonic_rms(cfg, specs, reference_crank_hz, reference_firing_hz, false, 1.0f);
    if (!(engine_rms_ref > 0.0f)) {
        return m;
    }
    const float shape = std::sqrt(1.0f +
                                  kDrivelineRmsFraction * kDrivelineRmsFraction +
                                  kBroadbandRmsFraction * kBroadbandRmsFraction);
    const float gain = cfg.level_mps2 / (shape * engine_rms_ref);

    const float speed_ratio = m.crank_hz / reference_crank_hz;
    const float driveline_rms_ref =
        engine_harmonic_rms(cfg, specs, reference_crank_hz, reference_firing_hz, true, 1.0f);
    const float driveline_gain = (driveline_rms_ref > 0.0f)
        ? (kDrivelineRmsFraction * engine_rms_ref / driveline_rms_ref)
        : 0.0f;

    const Vector3f axis_gain_raw(kAxisGainX, kAxisGainY, kAxisGainZ);
    const Vector3f axis_gain = axis_gain_raw / axis_gain_raw.norm();

    std::uniform_real_distribution<float> uphase(0.0f, kTwoPi);

    Vector3f hull_ms = Vector3f::Zero();
    Vector3f recorded_ms = Vector3f::Zero();

    for (const auto& spec : specs) {
        EngineVibrationLine line;
        line.freq_ratio = spec.freq_ratio;
        line.freq_hz = spec.freq_ratio * m.crank_hz;

        const float hull_amp = spec.driveline
            ? (gain * driveline_gain * spec.weight * speed_ratio * speed_ratio)
            : (gain * spec.weight *
               engine_unit_amplitude(cfg, line.freq_hz, reference_firing_hz));
        if (!(hull_amp > 0.0f)) continue;

        // The mechanical vibration reaches the sensing element in full; the
        // anti-alias path is what limits the part that survives to the sample.
        const float recorded_amp =
            hull_amp * engine_antialias_gain(line.freq_hz, cfg.sensor_bandwidth_hz);

        line.amp = axis_gain * recorded_amp;
        line.phase = uphase(m.rng);
        line.phase_offset = Vector3f(uphase(m.rng), uphase(m.rng), uphase(m.rng));
        // Anchor the per-axis offsets on the line phase so the axes stay
        // physically correlated rather than independent.
        line.phase_offset -= Vector3f::Constant(line.phase_offset.x());
        line.mod = 0.0f;
        m.lines.push_back(line);

        const Vector3f hull_axis = axis_gain * hull_amp;
        hull_ms += 0.5f * hull_axis.cwiseProduct(hull_axis);
        recorded_ms += 0.5f * line.amp.cwiseProduct(line.amp);
    }

    // Broadband structural floor.  Its recorded power scales with the sensor's
    // noise bandwidth, so a narrower anti-alias filter folds in less of it.
    const float engine_rms_now =
        engine_harmonic_rms(cfg, specs, m.crank_hz, reference_firing_hz, false, 1.0f);
    const float broadband_hull = gain * kBroadbandRmsFraction * engine_rms_now;
    const float bandwidth_ratio = (cfg.sensor_bandwidth_hz > 0.0f)
        ? std::sqrt(cfg.sensor_bandwidth_hz / kReferenceBandwidthHz)
        : 1.0f;
    m.broadband_sigma = axis_gain * (broadband_hull * bandwidth_ratio);
    hull_ms += (axis_gain * broadband_hull).cwiseAbs2();
    recorded_ms += m.broadband_sigma.cwiseAbs2();

    m.hp_alpha = 1.0f - std::exp(-kTwoPi * kBroadbandHighPassHz * dt);

    m.hull_rms = hull_ms.cwiseSqrt();
    m.recorded_rms = recorded_ms.cwiseSqrt();

    // Vibration rectification is a mechanical effect in the proof mass, so it
    // responds to the full hull vibration rather than to the filtered signal.
    // The specification is in mg of offset per g^2 of vibration.
    for (int i = 0; i < 3; ++i) {
        const float a_g = m.hull_rms[i] / g_std;
        m.vre_offset[i] = cfg.vre_mg_per_g2 * a_g * a_g * 1e-3f * g_std;
    }

    // Governor: 0.4 percent RMS speed deviation with a 2 s correlation time,
    // plus a slow periodic hunt.  This is what widens each order from a line
    // into a narrow band, and it is why an alias that lands near DC wanders
    // through the wave band instead of sitting at exactly zero.
    const float ou_tau_s = 2.0f;
    m.rpm_ou_alpha = 1.0f - std::exp(-dt / ou_tau_s);
    m.rpm_ou_sigma = 0.004f * std::sqrt(2.0f * m.rpm_ou_alpha);
    m.hunt_phase = uphase(m.rng);

    const float mod_tau_s = 1.5f;
    m.mod_alpha = 1.0f - std::exp(-dt / mod_tau_s);
    m.mod_sigma = std::sqrt(2.0f * m.mod_alpha);

    return m;
}

void apply_engine_vibration(Vector3f& acc_body_zu,
                            Vector3f& gyr_body_zu,
                            EngineVibrationModel& m,
                            float dt)
{
    if (m.lines.empty() && !(m.broadband_sigma.norm() > 0.0f)) return;

    // Optional shutdown part-way through the record.  Everything downstream
    // then sees a quiet sensor again, which is what exercises the release
    // side of an estimator's vibration handling rather than only the onset.
    m.time_sec += dt;
    if (m.cfg.stop_sec > 0.0f && m.time_sec >= m.cfg.stop_sec) return;

    // Instantaneous crank frequency: nominal, plus governor wander.
    m.rpm_ou += -m.rpm_ou_alpha * m.rpm_ou + m.rpm_ou_sigma * m.n01(m.rng);
    m.hunt_phase += kTwoPi * m.hunt_hz * dt;
    if (m.hunt_phase > kTwoPi) m.hunt_phase -= kTwoPi;
    const float speed_factor =
        1.0f + m.rpm_ou + m.hunt_amplitude * std::sin(m.hunt_phase);
    const float crank_hz_now = m.crank_hz * speed_factor;

    Vector3f vib = Vector3f::Zero();

    for (auto& line : m.lines) {
        // Advancing the phase at the sample rate is exactly what a sampled
        // sensor does, so a line above Nyquist folds on its own.  The wrap has
        // to be a real modulus rather than one conditional subtraction: a line
        // above Nyquist advances by more than 2 pi per sample, so subtracting
        // at most one turn would let the phase run away and lose float
        // precision over a long record.
        line.phase = std::fmod(
            line.phase + kTwoPi * line.freq_ratio * crank_hz_now * dt, kTwoPi);
        if (line.phase < 0.0f) line.phase += kTwoPi;

        line.mod += -m.mod_alpha * line.mod + m.mod_sigma * m.n01(m.rng);
        const float envelope = std::max(0.0f, 1.0f + m.mod_depth * line.mod);

        for (int i = 0; i < 3; ++i) {
            vib[i] += line.amp[i] * envelope *
                std::sin(line.phase + line.phase_offset[i]);
        }
    }

    if (m.broadband_sigma.norm() > 0.0f) {
        const Vector3f white(m.broadband_sigma.x() * m.n01(m.rng),
                             m.broadband_sigma.y() * m.n01(m.rng),
                             m.broadband_sigma.z() * m.n01(m.rng));
        // One-pole high pass, so the engine adds little directly in the wave
        // band and what does arrive there comes from the folded lines.
        m.hp_out = (1.0f - m.hp_alpha) * (m.hp_out + white - m.hp_prev_in);
        m.hp_prev_in = white;
        vib += m.hp_out;
    }

    acc_body_zu += vib + m.vre_offset;

    // The gyroscope sees hull angular vibration about the mounts, modeled by
    // an effective lever arm (a line of linear amplitude A at frequency f
    // corresponds to an angular rate A / (r 2 pi f)), plus the part of the
    // linear vibration its own g-sensitivity converts into rate.
    Vector3f rate = m.cfg.gyro_g_sensitivity * vib;
    if (m.cfg.gyro_lever_m > 0.0f) {
        for (const auto& line : m.lines) {
            const float omega = kTwoPi * line.freq_hz;
            if (!(omega > 0.0f)) continue;
            const float scale = 1.0f / (m.cfg.gyro_lever_m * omega);
            // Angular vibration is about the axes orthogonal to the linear
            // motion that drives it: vertical shake rocks the boat in roll and
            // pitch, athwartships shake rolls it.
            rate.x() += line.amp.z() * scale *
                std::sin(line.phase + line.phase_offset.z());
            rate.y() += line.amp.x() * scale *
                std::sin(line.phase + line.phase_offset.x());
            rate.z() += line.amp.y() * scale *
                std::sin(line.phase + line.phase_offset.y());
        }
    }
    gyr_body_zu += rate;
}

namespace {

bool w3d_engine_float_from_env(const char* name, float& out)
{
    const char* text = std::getenv(name);
    if (!text || *text == '\0') return false;

    errno = 0;
    char* end = nullptr;
    const double value = std::strtod(text, &end);
    if (errno == ERANGE || end == text || *end != '\0' || !std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be a finite number");
    }
    out = static_cast<float>(value);
    return true;
}

bool w3d_engine_int_from_env(const char* name, int& out)
{
    float value = 0.0f;
    if (!w3d_engine_float_from_env(name, value)) return false;
    if (value < 1.0f || value != std::floor(value)) {
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    }
    out = static_cast<int>(value);
    return true;
}

} // namespace

std::optional<EngineVibrationConfig> w3d_engine_vibration_from_env()
{
    EngineVibrationConfig cfg;
    if (!w3d_engine_float_from_env("W3D_ENGINE_RPM", cfg.rpm)) {
        return std::nullopt;
    }
    if (!(cfg.rpm > 0.0f)) {
        return std::nullopt;
    }

    w3d_engine_int_from_env("W3D_ENGINE_CYLINDERS", cfg.cylinders);
    w3d_engine_int_from_env("W3D_ENGINE_BLADES", cfg.blades);
    w3d_engine_float_from_env("W3D_ENGINE_GEAR_RATIO", cfg.gear_ratio);
    w3d_engine_float_from_env("W3D_ENGINE_LEVEL_MPS2", cfg.level_mps2);
    w3d_engine_float_from_env("W3D_ENGINE_REFERENCE_RPM", cfg.reference_rpm);
    w3d_engine_float_from_env("W3D_ENGINE_MOUNT_HZ", cfg.mount_hz);
    w3d_engine_float_from_env("W3D_ENGINE_MOUNT_ZETA", cfg.mount_zeta);
    w3d_engine_float_from_env("W3D_ENGINE_BANDWIDTH_HZ", cfg.sensor_bandwidth_hz);
    w3d_engine_float_from_env("W3D_ENGINE_GYRO_LEVER_M", cfg.gyro_lever_m);
    w3d_engine_float_from_env("W3D_ENGINE_GYRO_G_SENS", cfg.gyro_g_sensitivity);
    w3d_engine_float_from_env("W3D_ENGINE_VRE_MG_PER_G2", cfg.vre_mg_per_g2);
    w3d_engine_float_from_env("W3D_ENGINE_STOP_SEC", cfg.stop_sec);

    float seed = static_cast<float>(cfg.seed);
    if (w3d_engine_float_from_env("W3D_ENGINE_SEED", seed)) {
        if (seed < 0.0f || seed != std::floor(seed)) {
            throw std::invalid_argument("W3D_ENGINE_SEED must be an unsigned integer");
        }
        cfg.seed = static_cast<unsigned>(seed);
    }

    if (!(cfg.level_mps2 >= 0.0f)) {
        throw std::invalid_argument("W3D_ENGINE_LEVEL_MPS2 must be non-negative");
    }
    return cfg;
}

void w3d_install_engine_vibration_from_env(SimulationNoiseModels& noise_models, float dt)
{
    const auto cfg = w3d_engine_vibration_from_env();
    if (!cfg) return;

    auto model = std::make_shared<EngineVibrationModel>(
        make_engine_vibration_model(*cfg, dt));

    const float nyquist_hz = (dt > 0.0f) ? (0.5f / dt) : 0.0f;
    auto alias_hz = [nyquist_hz](float freq_hz) {
        if (!(nyquist_hz > 0.0f)) return freq_hz;
        const float fs = 2.0f * nyquist_hz;
        float folded = std::fmod(freq_hz, fs);
        if (folded < 0.0f) folded += fs;
        if (folded > nyquist_hz) folded = fs - folded;
        return folded;
    };

    std::cout << "ENGINE_VIBRATION rpm=" << cfg->rpm
              << " cylinders=" << cfg->cylinders
              << " firing_hz=" << (model->firing_order * model->crank_hz)
              << " shaft_hz=" << model->shaft_hz
              << " blade_hz=" << (model->shaft_hz * static_cast<float>(cfg->blades))
              << " level_mps2=" << cfg->level_mps2
              << " bandwidth_hz=" << cfg->sensor_bandwidth_hz
              << " hull_rms_mps2=" << model->hull_rms.norm()
              << " recorded_rms_mps2=" << model->recorded_rms.norm()
              << " vre_offset_mps2=" << model->vre_offset.norm()
              << " lines=" << model->lines.size()
              << " stop_sec=" << cfg->stop_sec
              << "\n";
    for (const auto& line : model->lines) {
        std::cout << "ENGINE_LINE order=" << (line.freq_ratio)
                  << " freq_hz=" << line.freq_hz
                  << " alias_hz=" << alias_hz(line.freq_hz)
                  << " amp_mps2=" << line.amp.norm()
                  << "\n";
    }

    noise_models.extra_imu_noise_models.push_back(
        [model](Vector3f& acc_body_zu, Vector3f& gyr_body_zu, float step) {
            apply_engine_vibration(acc_body_zu, gyr_body_zu, *model, step);
        });
}

// ---------------------------------------------------------------------------
// IMU installation lever arm
// ---------------------------------------------------------------------------

Vector3f W3dRateDerivative::update(const Vector3f& omega)
{
    Vector3f alpha = Vector3f::Zero();
    if (dt_ > 0.0f) {
        if (seen_ >= 2) {
            // Second-order backward difference.  Causal, so the stage can run
            // inside the streaming record loop without buffering the future,
            // and its O(dt^2) error is negligible next to a rotation whose
            // energy sits below 1 Hz at a 200 Hz sample rate.
            alpha = (3.0f * omega - 4.0f * prev_ + prev2_) / (2.0f * dt_);
        } else if (seen_ == 1) {
            alpha = (omega - prev_) / dt_;
        }
    }
    prev2_ = prev_;
    prev_ = omega;
    if (seen_ < 2) ++seen_;
    return alpha;
}

W3dLeverArm::W3dLeverArm(const W3dLeverArmConfig& cfg, float dt)
    : cfg_(cfg),
      offset_(cfg.offset_body_zu),
      offset_norm_(cfg.offset_body_zu.norm()),
      truth_derivative_(dt),
      model_derivative_(dt)
{
    if (!(dt > 0.0f) || !std::isfinite(dt)) {
        throw std::invalid_argument("lever-arm stage needs a positive dt");
    }
    if (!offset_.allFinite()) {
        throw std::invalid_argument("W3D_IMU_LEVER_ARM_M must be finite");
    }
    // Two cascaded one-pole sections, i.e. a critically damped second-order
    // low-pass with corner cfg.derivative_cutoff_hz.  A cutoff at or above
    // Nyquist degenerates to a pass-through, which is the honest behaviour
    // for "do not band-limit".
    const float cutoff = cfg_.derivative_cutoff_hz;
    if (cutoff > 0.0f && std::isfinite(cutoff)) {
        const float tau = 1.0f / (2.0f * float(std::numbers::pi_v<float>) * cutoff);
        lp_gain_ = dt / (tau + dt);
        if (lp_gain_ > 1.0f) lp_gain_ = 1.0f;
    } else {
        lp_gain_ = 1.0f;
    }
}

Vector3f W3dLeverArm::install(const Vector3f& acc_cg_body_zu,
                              const Vector3f& gyr_truth_body_zu)
{
    const Vector3f alpha = truth_derivative_.update(gyr_truth_body_zu);
    last_installed_ = installs()
        ? w3d_lever_acceleration(gyr_truth_body_zu, alpha, offset_)
        : Vector3f::Zero();
    installed_sumsq_ += double(last_installed_.squaredNorm());
    ++samples_;
    return acc_cg_body_zu + last_installed_;
}

Vector3f W3dLeverArm::compensate(const Vector3f& acc_meas_body_zu,
                                 const Vector3f& gyr_meas_body_zu)
{
    Vector3f modelled = Vector3f::Zero();
    if (cfg_.model == W3dLeverArmConfig::Model::Exact) {
        // The oracle reuses the term the installation stage applied, which is
        // the exact rigid-body acceleration at r.  It answers how much of the
        // penalty is deterministically recoverable, nothing more.
        modelled = last_installed_;
    } else if (cfg_.model == W3dLeverArmConfig::Model::MeasuredGyro) {
        // The deployable model sees only the corrupted rate.  Band-limit it
        // to the rigid-body band before differencing, then evaluate the same
        // rigid-body expression on what is left.
        if (!lp_primed_) {
            lp_stage1_ = gyr_meas_body_zu;
            lp_stage2_ = gyr_meas_body_zu;
            lp_primed_ = true;
        } else {
            lp_stage1_ += lp_gain_ * (gyr_meas_body_zu - lp_stage1_);
            lp_stage2_ += lp_gain_ * (lp_stage1_ - lp_stage2_);
        }
        const Vector3f alpha = model_derivative_.update(lp_stage2_);
        modelled = w3d_lever_acceleration(lp_stage2_, alpha, offset_);
    }

    // With no model this is the whole installed term, which is exactly what
    // the unmodeled arm is meant to report.
    residual_sumsq_ += double((last_installed_ - modelled).squaredNorm());
    return acc_meas_body_zu - modelled;
}

float W3dLeverArm::installed_rms_mps2() const
{
    if (samples_ == 0) return 0.0f;
    return float(std::sqrt(installed_sumsq_ / double(samples_)));
}

float W3dLeverArm::residual_rms_mps2() const
{
    if (samples_ == 0) return 0.0f;
    return float(std::sqrt(residual_sumsq_ / double(samples_)));
}

namespace {

W3dLeverArmConfig::Model w3d_lever_arm_model_from_text(const std::string& text)
{
    if (text == "none" || text == "off" || text == "unmodeled") {
        return W3dLeverArmConfig::Model::None;
    }
    if (text == "exact" || text == "oracle") {
        return W3dLeverArmConfig::Model::Exact;
    }
    if (text == "gyro" || text == "measured") {
        return W3dLeverArmConfig::Model::MeasuredGyro;
    }
    throw std::invalid_argument(
        "W3D_IMU_LEVER_ARM_MODEL must be none, exact, or gyro");
}

Vector3f w3d_vector_from_text(const char* name, const std::string& text)
{
    Vector3f out = Vector3f::Zero();
    std::stringstream stream(text);
    std::string field;
    int count = 0;
    while (std::getline(stream, field, ',')) {
        if (count >= 3) {
            throw std::invalid_argument(std::string(name) + " must be \"x,y,z\"");
        }
        errno = 0;
        char* end = nullptr;
        const double value = std::strtod(field.c_str(), &end);
        if (errno == ERANGE || end == field.c_str() || *end != '\0' ||
            !std::isfinite(value))
        {
            throw std::invalid_argument(
                std::string(name) + " must be three finite numbers");
        }
        out[count] = static_cast<float>(value);
        ++count;
    }
    if (count != 3) {
        throw std::invalid_argument(std::string(name) + " must be \"x,y,z\"");
    }
    return out;
}

const char* w3d_lever_arm_model_name(W3dLeverArmConfig::Model model)
{
    switch (model) {
        case W3dLeverArmConfig::Model::Exact: return "exact";
        case W3dLeverArmConfig::Model::MeasuredGyro: return "gyro";
        case W3dLeverArmConfig::Model::None: break;
    }
    return "none";
}

} // namespace

std::optional<W3dLeverArmConfig> w3d_lever_arm_config_from_env()
{
    const char* offset = std::getenv("W3D_IMU_LEVER_ARM_M");
    if (!offset || *offset == '\0') return std::nullopt;

    W3dLeverArmConfig cfg;
    cfg.offset_body_zu = w3d_vector_from_text("W3D_IMU_LEVER_ARM_M", offset);

    if (const char* model = std::getenv("W3D_IMU_LEVER_ARM_MODEL")) {
        if (*model != '\0') {
            cfg.model = w3d_lever_arm_model_from_text(model);
        }
    }
    // Shares the engine model's numeric env reader; the name is historical,
    // the parsing and the error text are not engine specific.
    (void)w3d_engine_float_from_env("W3D_IMU_LEVER_ARM_CUTOFF_HZ",
                                    cfg.derivative_cutoff_hz);
    if (!(cfg.derivative_cutoff_hz > 0.0f)) {
        throw std::invalid_argument(
            "W3D_IMU_LEVER_ARM_CUTOFF_HZ must be positive");
    }
    return cfg;
}

std::shared_ptr<W3dLeverArm> w3d_lever_arm_from_env(float dt)
{
    const auto cfg = w3d_lever_arm_config_from_env();
    if (!cfg) return nullptr;

    auto stage = std::make_shared<W3dLeverArm>(*cfg, dt);
    std::cout << "IMU_LEVER_ARM offset_m=[" << cfg->offset_body_zu.transpose()
              << "] norm_m=" << stage->offset_norm_m()
              << " model=" << w3d_lever_arm_model_name(cfg->model)
              << " derivative_cutoff_hz=" << cfg->derivative_cutoff_hz
              << "\n";
    return stage;
}

void w3d_report_lever_arm(const std::shared_ptr<W3dLeverArm>& lever_arm)
{
    if (!lever_arm) return;
    std::cout << "IMU_LEVER_ARM_RESULT norm_m=" << lever_arm->offset_norm_m()
              << " model=" << w3d_lever_arm_model_name(lever_arm->config().model)
              << " samples=" << lever_arm->samples()
              << " installed_rms_mps2=" << lever_arm->installed_rms_mps2()
              << " residual_rms_mps2=" << lever_arm->residual_rms_mps2()
              << "\n";
}

unsigned w3d_expand_seed(unsigned base_seed, unsigned stream_id)
{
    // SplitMix32 finalizer. Stream IDs make the expanded streams independent
    // while keeping the mapping stable across platforms and runs.
    std::uint32_t z = static_cast<std::uint32_t>(base_seed) +
        0x9e3779b9u * (static_cast<std::uint32_t>(stream_id) + 1u);
    z = (z ^ (z >> 16u)) * 0x85ebca6bu;
    z = (z ^ (z >> 13u)) * 0xc2b2ae35u;
    z ^= z >> 16u;
    return static_cast<unsigned>(z);
}

namespace {

std::optional<unsigned> w3d_seed_from_env(const char* name)
{
    const char* text = std::getenv(name);
    if (!text) return std::nullopt;
    if (*text == '\0' || *text == '-') {
        throw std::invalid_argument(std::string(name) + " must be an unsigned integer");
    }

    errno = 0;
    char* end = nullptr;
    const unsigned long value = std::strtoul(text, &end, 10);
    if (errno == ERANGE || end == text || *end != '\0' ||
        value > std::numeric_limits<unsigned>::max())
    {
        throw std::invalid_argument(std::string(name) + " must be an unsigned integer");
    }
    return static_cast<unsigned>(value);
}

} // namespace

W3dRandomSeeds w3d_random_seeds_from_env()
{
    W3dRandomSeeds seeds;
    const auto combined = w3d_seed_from_env("W3D_SEED");
    const auto imu = w3d_seed_from_env("W3D_IMU_SEED");
    const auto initialization = w3d_seed_from_env("W3D_INIT_SEED");

    if (!combined && !imu && !initialization) {
        return seeds;
    }

    if (combined) {
        seeds.accel_noise = w3d_expand_seed(*combined, 0u);
        seeds.gyro_noise = w3d_expand_seed(*combined, 1u);
        seeds.mag_noise = w3d_expand_seed(*combined, 2u);
        seeds.accel_initialization = w3d_expand_seed(*combined, 3u);
        seeds.gyro_initialization = w3d_expand_seed(*combined, 4u);
        seeds.mag_initialization = w3d_expand_seed(*combined, 5u);
    }
    if (imu) {
        seeds.accel_noise = w3d_expand_seed(*imu, 0u);
        seeds.gyro_noise = w3d_expand_seed(*imu, 1u);
        seeds.mag_noise = w3d_expand_seed(*imu, 2u);
    }
    if (initialization) {
        seeds.accel_initialization = w3d_expand_seed(*initialization, 3u);
        seeds.gyro_initialization = w3d_expand_seed(*initialization, 4u);
        seeds.mag_initialization = w3d_expand_seed(*initialization, 5u);
    }
    return seeds;
}

Vector3f apply_mag_noise(const Vector3f& ideal_mag_uT_body, MagNoiseModel& m, float dt_mag)
{
    if (m.sigma_bias_rw_uT_sqrt_s > 0.0f) {
        const float s = m.sigma_bias_rw_uT_sqrt_s * std::sqrt(dt_mag);
        m.bias_rw_uT += Vector3f(s * m.n01(m.rng), s * m.n01(m.rng), s * m.n01(m.rng));
    }
    Vector3f white(m.w_uT(m.rng), m.w_uT(m.rng), m.w_uT(m.rng));
    return (m.Mis * ideal_mag_uT_body) + (m.bias0_uT + m.bias_rw_uT) + white;
}

static void write_tvg_nlo_csv_header(std::ofstream& ofs)
{
    ofs << ",tvg_k1"
        << ",tvg_k2"
        << ",tvg_kI"
        << ",tvg_vartheta"
        << ",tvg_theta"
        << ",tvg_p0z_hat"
        << ",tvg_wave_freq_hz"
        << ",tvg_wave_freq_confidence"

        << ",tvg_xi_n_x"
        << ",tvg_xi_n_y"
        << ",tvg_xi_n_z"
        << ",tvg_xi_norm"

        << ",tvg_fhat_n_x"
        << ",tvg_fhat_n_y"
        << ",tvg_fhat_n_z"
        << ",tvg_fhat_norm"

        << ",tvg_sigma_b_x"
        << ",tvg_sigma_b_y"
        << ",tvg_sigma_b_z"
        << ",tvg_sigma_norm"

        << ",tvg_gyro_bias_b_x"
        << ",tvg_gyro_bias_b_y"
        << ",tvg_gyro_bias_b_z"
        << ",tvg_gyro_bias_norm";
}

static void write_tvg_nlo_csv_row(std::ofstream& ofs, const TvgNloFilterSnapshot& snap)
{
    const auto& d = snap.tvg;

    ofs << "," << d.k1
        << "," << d.k2
        << "," << d.kI
        << "," << d.vartheta
        << "," << d.theta
        << "," << d.p0z_hat
        << "," << d.wave_freq_hz
        << "," << d.wave_freq_confidence

        << "," << d.xi_n.x()
        << "," << d.xi_n.y()
        << "," << d.xi_n.z()
        << "," << d.xi_norm

        << "," << d.fhat_n.x()
        << "," << d.fhat_n.y()
        << "," << d.fhat_n.z()
        << "," << d.fhat_norm

        << "," << d.sigma_b.x()
        << "," << d.sigma_b.y()
        << "," << d.sigma_b.z()
        << "," << d.sigma_norm

        << "," << d.gyro_bias_b.x()
        << "," << d.gyro_bias_b.y()
        << "," << d.gyro_bias_b.z()
        << "," << d.gyro_bias_norm;
}

W3dSimulationRunner::W3dSimulationRunner(W3dSimulationOptions options,
                                         SimulationNoiseModels noise_models,
                                         IW3dFusionAdapter& fusion_adapter)
    : options_(std::move(options)),
      noise_models_(std::move(noise_models)),
      fusion_adapter_(fusion_adapter)
{
}

std::string W3dSimulationRunner::make_output_name(const std::string& filename) const
{
    std::string outname = filename;
    auto pos_prefix = outname.find("wave_data_");
    if (pos_prefix != std::string::npos) {
        outname.replace(pos_prefix, std::string("wave_data_").size(), "w3d_");
    } else {
        outname = "w3d_" + outname;
    }

    auto pos_ext = outname.rfind(".csv");
    const std::string& suffix = options_.with_mag ? options_.output_suffix_with_mag
                                                  : options_.output_suffix_no_mag;
    if (pos_ext != std::string::npos) {
        outname.insert(pos_ext, suffix);
    } else {
        outname += suffix + std::string(".csv");
    }
    return outname;
}

TvgNloSimulationRunner::TvgNloSimulationRunner(W3dSimulationOptions options,
                                               SimulationNoiseModels noise_models,
                                               Adapter& fusion_adapter)
    : options_(std::move(options)),
      noise_models_(std::move(noise_models)),
      fusion_adapter_(fusion_adapter)
{
}

std::string TvgNloSimulationRunner::make_output_name(const std::string& filename) const
{
    std::string outname = filename;
    auto pos_prefix = outname.find("wave_data_");
    if (pos_prefix != std::string::npos) {
        outname.replace(pos_prefix, std::string("wave_data_").size(), "w3d_");
    } else {
        outname = "w3d_" + outname;
    }

    auto pos_ext = outname.rfind(".csv");
    const std::string& suffix = options_.with_mag ? options_.output_suffix_with_mag
                                                  : options_.output_suffix_no_mag;
    if (pos_ext != std::string::npos) {
        outname.insert(pos_ext, suffix);
    } else {
        outname += suffix + std::string(".csv");
    }
    return outname;
}

std::optional<W3dSimulationRunResult> W3dSimulationRunner::run(const std::string& filename)
{
    auto parsed = WaveFileNaming::parse_to_params(filename);
    if (!parsed) return std::nullopt;

    auto [kind, type, wp] = *parsed;
    if (kind != FileKind::Data) return std::nullopt;
    if (!(type == WaveType::JONSWAP || type == WaveType::PMSTOKES)) return std::nullopt;

    W3dSimulationRunResult result;
    result.input_name = filename;
    result.output_name = make_output_name(filename);
    result.wave_type = type;
    result.wave_params = wp;
    result.with_mag = options_.with_mag;

    std::cout << "Processing " << filename << " (type="
              << EnumTraits<WaveType>::to_string(type)
              << ")\n";

    std::ofstream ofs;
    if (options_.write_timeseries) {
        ofs.open(result.output_name);
        ofs << "time,roll_ref,pitch_ref,yaw_ref,"
        << "disp_ref_x,disp_ref_y,disp_ref_z,"
        << "vel_ref_x,vel_ref_y,vel_ref_z,"
        << "acc_ref_x,acc_ref_y,acc_ref_z,"
        << "roll_est,pitch_est,yaw_est,"
        << "disp_est_x,disp_est_y,disp_est_z,"
        << "vel_est_x,vel_est_y,vel_est_z,"
        << "acc_est_x,acc_est_y,acc_est_z,"
        << "acc_bias_x,acc_bias_y,acc_bias_z,"
        << "gyro_bias_x,gyro_bias_y,gyro_bias_z,"
        << "acc_bias_est_x,acc_bias_est_y,acc_bias_est_z,"
        << "gyro_bias_est_x,gyro_bias_est_y,gyro_bias_est_z,"
        << "mag_bias_x,mag_bias_y,mag_bias_z,"
        << "mag_bias_est_x,mag_bias_est_y,mag_bias_est_z,"
        << "mag_bias_err_x,mag_bias_err_y,mag_bias_err_z,"
        << "tau_applied,sigma_a_applied,R_p0_applied,"
        << "freq_tracker_hz,Tp_tuner_s,accel_var_tuner,"
        << "disp_scale_m,vel_scale_mps,"
        << "dir_phase,"
        << "dir_axis_deg,dir_apparent_to_deg,dir_apparent_from_deg,"
        << "dir_sense_coherence,dir_uncert_deg,dir_conf,dir_amp,"
        << "dir_sign,dir_sign_num,"
        << "dir_vec_x,dir_vec_y,"
            << "dfilt_ax,dfilt_ay\n";
    }

    const Vector3f mag_world_zu = MagSim_WMM::mag_world_nautical();

    const float mag_dt = (options_.mag_odr_hz > 0.0f)
        ? (1.0f / options_.mag_odr_hz)
        : std::numeric_limits<float>::infinity();

    float mag_phase_s = 0.0f;
    Vector3f mag_body_ned_hold = Vector3f::Zero();

    auto quat_from_csv = [](const IMU_Sample& imu) -> Quaternionf {
        Quaternionf q(imu.q_wb_zu_w, imu.q_wb_zu_x, imu.q_wb_zu_y, imu.q_wb_zu_z);
        const bool finite =
            std::isfinite(q.w()) && std::isfinite(q.x()) &&
            std::isfinite(q.y()) && std::isfinite(q.z());
        if (!finite || q.norm() < 1e-6f) {
            return Quaternionf::Identity();
        }
        q.normalize();
        return q;
    };

    auto reference_euler_from_csv = [&](const IMU_Sample& imu,
                                        const Quaternionf& q_wb_zu,
                                        float& roll_deg,
                                        float& pitch_deg,
                                        float& yaw_deg) {
        float r_w = 0.0f, p_w = 0.0f, y_w = 0.0f;
        quat_wb_zu_to_euler_nautical(q_wb_zu, r_w, p_w, y_w);
        (void)imu;
        roll_deg = r_w;
        pitch_deg = p_w;
        yaw_deg = y_w;
    };

    WaveDataCSVReader reader(filename);
    reader.for_each_record([&](const Wave_Data_Sample& rec) {
        Vector3f acc_b(rec.imu.acc_bx, rec.imu.acc_by, rec.imu.acc_bz);
        Vector3f gyr_b(rec.imu.gyro_x, rec.imu.gyro_y, rec.imu.gyro_z);

        // Installation physics, on the record's own truth: an IMU at r sees
        // the CG specific force plus the rigid-body rotational terms.
        if (options_.lever_arm) {
            acc_b = options_.lever_arm->install(acc_b, gyr_b);
        }

        if (options_.add_noise) {
            if (noise_models_.accel_noise) {
                acc_b = apply_imu_noise(acc_b, *noise_models_.accel_noise, options_.dt);
            }
            if (noise_models_.gyro_noise) {
                gyr_b = apply_imu_noise(gyr_b, *noise_models_.gyro_noise, options_.dt);
            }
            for (auto& model : noise_models_.extra_imu_noise_models) {
                model(acc_b, gyr_b, options_.dt);
            }
        }

        // Filter-side lever-arm model, on the corrupted signals the estimator
        // actually receives.  Nothing downstream of this point knows the IMU
        // is off the CG.
        if (options_.lever_arm) {
            acc_b = options_.lever_arm->compensate(acc_b, gyr_b);
        }

        const Vector3f acc_meas_ned = zu_to_ned(acc_b);
        const Vector3f gyr_meas_ned = zu_to_ned(gyr_b);

        const Quaternionf q_ref_wb_zu = quat_from_csv(rec.imu);
        const Matrix3f C_wb_zu = q_ref_wb_zu.toRotationMatrix();

        float r_ref_out = 0.0f;
        float p_ref_out = 0.0f;
        float y_ref_out = 0.0f;

        reference_euler_from_csv(rec.imu, q_ref_wb_zu, r_ref_out, p_ref_out, y_ref_out);

        // Vessel heading in the record's own convention, before the magnetic
        // declination offset below moves the yaw reference into the frame the
        // filter actually learns.  Travel-sense scoring needs the physical
        // heading so the boat frame can be removed from the directed estimate.
        const float heading_ref_deg = y_ref_out;

        if (options_.with_mag) {
            y_ref_out = wrapDeg(y_ref_out + MagSim_WMM::default_declination_deg);
        }

        fusion_adapter_.update(options_.dt, gyr_meas_ned, acc_meas_ned, options_.temperature_c);

        if (options_.with_mag) {
            mag_phase_s += options_.dt;
            bool mag_tick = false;
            if (mag_phase_s >= mag_dt) {
                while (mag_phase_s >= mag_dt) mag_phase_s -= mag_dt;
                mag_tick = true;
            }
            if (mag_tick) {
                Vector3f mag_b_enu = C_wb_zu * mag_world_zu;

                if (options_.add_noise && noise_models_.mag_noise) {
                    mag_b_enu = apply_mag_noise(mag_b_enu, *noise_models_.mag_noise, mag_dt);
                }
                if (options_.add_noise) {
                    for (auto& model : noise_models_.extra_mag_noise_models) {
                        model(mag_b_enu, mag_dt);
                    }
                }
                mag_body_ned_hold = zu_to_ned(mag_b_enu);
                fusion_adapter_.updateMag(mag_body_ned_hold);
            }
        }

        FilterSnapshot snap = fusion_adapter_.snapshot();

        // The filter has no declination input, so with the magnetometer live
        // its world horizontal axes are magnetic, while the record's disp/vel/
        // acc columns are geographic.  Scoring x and y across that offset
        // charges every horizontal metric a fixed rotation the vertical axis
        // cannot have.  The yaw reference above is already moved into the
        // filter's frame; do the mirror move here and bring the translational
        // estimate back into the record's frame, before anything reads it.
        if (options_.with_mag) {
            snap.disp_est_zu = MagSim_WMM::mag_frame_to_geographic_nautical(snap.disp_est_zu);
            snap.vel_est_zu  = MagSim_WMM::mag_frame_to_geographic_nautical(snap.vel_est_zu);
            snap.acc_est_zu  = MagSim_WMM::mag_frame_to_geographic_nautical(snap.acc_est_zu);
        }

        Vector3f disp_ref(rec.wave.disp_x, rec.wave.disp_y, rec.wave.disp_z);
        Vector3f vel_ref(rec.wave.vel_x, rec.wave.vel_y, rec.wave.vel_z);
        Vector3f acc_ref(rec.wave.acc_x, rec.wave.acc_y, rec.wave.acc_z);

        Vector3f disp_err = snap.disp_est_zu - disp_ref;
        result.errs_x.push_back(disp_err.x());
        result.errs_y.push_back(disp_err.y());
        result.errs_z.push_back(disp_err.z());
        result.ref_x.push_back(disp_ref.x());
        result.ref_y.push_back(disp_ref.y());
        result.ref_z.push_back(disp_ref.z());
        result.errs_roll.push_back(diffDeg(snap.euler_nautical_deg.x(), r_ref_out));
        result.errs_pitch.push_back(diffDeg(snap.euler_nautical_deg.y(), p_ref_out));
        result.errs_yaw.push_back(diffDeg(snap.euler_nautical_deg.z(), y_ref_out));

        const Vector3f acc_bias_true_zu = (options_.add_noise && noise_models_.accel_noise)
            ? (noise_models_.accel_noise->bias0 + noise_models_.accel_noise->bias_rw).eval()
            : Vector3f::Zero().eval();
        const Vector3f gyro_bias_true_zu = (options_.add_noise && noise_models_.gyro_noise)
            ? (noise_models_.gyro_noise->bias0 + noise_models_.gyro_noise->bias_rw).eval()
            : Vector3f::Zero().eval();
        const Vector3f acc_bias_true_ned = zu_to_ned(acc_bias_true_zu);
        const Vector3f gyro_bias_true_ned = zu_to_ned(gyro_bias_true_zu);

        const Vector3f acc_bias_err = snap.acc_bias_est_ned - acc_bias_true_ned;
        const Vector3f gyro_bias_err = snap.gyro_bias_est_ned - gyro_bias_true_ned;

        const Vector3f mag_bias_true_zu = (options_.add_noise && options_.with_mag && noise_models_.mag_noise)
            ? (noise_models_.mag_noise->bias0_uT + noise_models_.mag_noise->bias_rw_uT).eval()
            : Vector3f::Zero().eval();
        const Vector3f mag_bias_true_ned = zu_to_ned(mag_bias_true_zu);
        const Vector3f mag_bias_err = snap.mag_bias_est_ned_uT - mag_bias_true_ned;

        result.accb_err_x.push_back(acc_bias_err.x());
        result.accb_err_y.push_back(acc_bias_err.y());
        result.accb_err_z.push_back(acc_bias_err.z());
        result.gyrb_err_x.push_back(gyro_bias_err.x());
        result.gyrb_err_y.push_back(gyro_bias_err.y());
        result.gyrb_err_z.push_back(gyro_bias_err.z());
        result.magb_err_x.push_back(mag_bias_err.x());
        result.magb_err_y.push_back(mag_bias_err.y());
        result.magb_err_z.push_back(mag_bias_err.z());

        result.accb_true_x.push_back(acc_bias_true_ned.x());
        result.accb_true_y.push_back(acc_bias_true_ned.y());
        result.accb_true_z.push_back(acc_bias_true_ned.z());
        result.gyrb_true_x.push_back(gyro_bias_true_ned.x());
        result.gyrb_true_y.push_back(gyro_bias_true_ned.y());
        result.gyrb_true_z.push_back(gyro_bias_true_ned.z());
        result.magb_true_x.push_back(mag_bias_true_ned.x());
        result.magb_true_y.push_back(mag_bias_true_ned.y());
        result.magb_true_z.push_back(mag_bias_true_ned.z());

        result.freq_hist.push_back(snap.freq_hz);
        result.dir_phase_hist.push_back(snap.direction.phase);
        // The axis arrives in the boat frame; adding the vessel heading puts it
        // in the generator frame the record azimuth lives in.  Every shipped
        // record has heading 0, which is why this was invisible until a
        // heading-rotated record was scored.
        result.dir_deg_hist.push_back(wrapAxialDeg90(
            snap.direction.direction_deg_generator_signed + heading_ref_deg));
        result.dir_unc_hist.push_back(snap.direction.uncertainty_deg);
        result.dir_conf_hist.push_back(snap.direction.confidence);
        result.dir_amp_hist.push_back(snap.direction.amplitude);
        result.dir_sign_num_hist.push_back(snap.direction.sign_num);
        result.dir_travel_deg_hist.push_back(
            travelDegGeneratorFromVec(snap.direction.travel_vec_boat, heading_ref_deg));

        if (options_.write_timeseries) {
            ofs << rec.time << ","
                << r_ref_out << "," << p_ref_out << "," << y_ref_out << ","
                << disp_ref.x() << "," << disp_ref.y() << "," << disp_ref.z() << ","
                << vel_ref.x() << "," << vel_ref.y() << "," << vel_ref.z() << ","
                << acc_ref.x() << "," << acc_ref.y() << "," << acc_ref.z() << ","
                << snap.euler_nautical_deg.x() << "," << snap.euler_nautical_deg.y() << "," << snap.euler_nautical_deg.z() << ","
                << snap.disp_est_zu.x() << "," << snap.disp_est_zu.y() << "," << snap.disp_est_zu.z() << ","
                << snap.vel_est_zu.x() << "," << snap.vel_est_zu.y() << "," << snap.vel_est_zu.z() << ","
                << snap.acc_est_zu.x() << "," << snap.acc_est_zu.y() << "," << snap.acc_est_zu.z() << ","
                << acc_bias_true_ned.x() << "," << acc_bias_true_ned.y() << "," << acc_bias_true_ned.z() << ","
                << gyro_bias_true_ned.x() << "," << gyro_bias_true_ned.y() << "," << gyro_bias_true_ned.z() << ","
                << snap.acc_bias_est_ned.x() << "," << snap.acc_bias_est_ned.y() << "," << snap.acc_bias_est_ned.z() << ","
                << snap.gyro_bias_est_ned.x() << "," << snap.gyro_bias_est_ned.y() << "," << snap.gyro_bias_est_ned.z() << ","
                << mag_bias_true_ned.x() << "," << mag_bias_true_ned.y() << "," << mag_bias_true_ned.z() << ","
                << snap.mag_bias_est_ned_uT.x() << "," << snap.mag_bias_est_ned_uT.y() << "," << snap.mag_bias_est_ned_uT.z() << ","
                << mag_bias_err.x() << "," << mag_bias_err.y() << "," << mag_bias_err.z() << ","
                << snap.tau_applied << ","
                << snap.sigma_applied << ","
                << snap.tuning_applied << ","
                << snap.freq_hz << ","
                << snap.period_sec << ","
                << snap.accel_variance << ","
                << snap.displacement_scale_m << ","
                << snap.velocity_scale_mps << ","
                << snap.direction.phase << "," << snap.direction.direction_deg << ","
                << snap.direction.apparent_to_deg << "," << snap.direction.apparent_from_deg << ","
                << snap.direction.sense_coherence << "," << snap.direction.uncertainty_deg << ","
                << snap.direction.confidence << "," << snap.direction.amplitude << ","
                << (snap.direction.sign == FORWARD ? "POSITIVE_AXIS" : snap.direction.sign == BACKWARD ? "NEGATIVE_AXIS" : "UNCERTAIN") << ","
                << snap.direction.sign_num << ","
                << snap.direction.direction_vec.x() << "," << snap.direction.direction_vec.y() << ","
                << snap.direction.filtered_signal.x() << "," << snap.direction.filtered_signal.y() << "\n";
        }

        result.final_tau_target = snap.tau_target;
        result.final_sigma_target = snap.sigma_target;
        result.final_tuning_target = snap.tuning_target;
        result.final_tau_applied = snap.tau_applied;
        result.final_sigma_applied = snap.sigma_applied;
        result.final_tuning_applied = snap.tuning_applied;
        result.final_wave_period_sec = snap.wave_period_sec;
        result.final_freq_hz = snap.freq_hz;
        result.final_period_sec = snap.period_sec;
        result.final_accel_variance = snap.accel_variance;
    });

    if (options_.write_timeseries) {
        ofs.close();
        std::cout << "Wrote " << result.output_name << "\n";
    }
    return result;
}

std::optional<TvgNloSimulationRunResult> TvgNloSimulationRunner::run(const std::string& filename)
{
    auto parsed = WaveFileNaming::parse_to_params(filename);
    if (!parsed) return std::nullopt;

    auto [kind, type, wp] = *parsed;
    if (kind != FileKind::Data) return std::nullopt;
    if (!(type == WaveType::JONSWAP || type == WaveType::PMSTOKES)) return std::nullopt;

    TvgNloSimulationRunResult result;
    result.input_name = filename;
    result.output_name = make_output_name(filename);
    result.wave_type = type;
    result.wave_params = wp;
    result.with_mag = options_.with_mag;

    std::cout << "Processing " << filename << " (type="
              << EnumTraits<WaveType>::to_string(type)
              << ")\n";

    std::ofstream ofs(result.output_name);
    ofs << "time,roll_ref,pitch_ref,yaw_ref,"
        << "disp_ref_x,disp_ref_y,disp_ref_z,"
        << "vel_ref_x,vel_ref_y,vel_ref_z,"
        << "acc_ref_x,acc_ref_y,acc_ref_z,"
        << "roll_est,pitch_est,yaw_est,"
        << "disp_est_x,disp_est_y,disp_est_z,"
        << "vel_est_x,vel_est_y,vel_est_z,"
        << "acc_est_x,acc_est_y,acc_est_z,"
        << "acc_bias_x,acc_bias_y,acc_bias_z,"
        << "gyro_bias_x,gyro_bias_y,gyro_bias_z,"
        << "acc_bias_est_x,acc_bias_est_y,acc_bias_est_z,"
        << "gyro_bias_est_x,gyro_bias_est_y,gyro_bias_est_z,"
        << "mag_bias_x,mag_bias_y,mag_bias_z,"
        << "mag_bias_est_x,mag_bias_est_y,mag_bias_est_z,"
        << "mag_bias_err_x,mag_bias_err_y,mag_bias_err_z,"
        << "dir_phase,"
        << "dir_axis_deg,dir_apparent_to_deg,dir_apparent_from_deg,"
        << "dir_sense_coherence,dir_uncert_deg,dir_conf,dir_amp,"
        << "dir_sign,dir_sign_num,"
        << "dir_vec_x,dir_vec_y,"
        << "dfilt_ax,dfilt_ay";

    write_tvg_nlo_csv_header(ofs);
    ofs << "\n";

    const Vector3f mag_world_zu = MagSim_WMM::mag_world_nautical();

    const float mag_dt = (options_.mag_odr_hz > 0.0f)
        ? (1.0f / options_.mag_odr_hz)
        : std::numeric_limits<float>::infinity();

    float mag_phase_s = 0.0f;
    Vector3f mag_body_ned_hold = Vector3f::Zero();

    auto quat_from_csv = [](const IMU_Sample& imu) -> Quaternionf {
        Quaternionf q(imu.q_wb_zu_w, imu.q_wb_zu_x, imu.q_wb_zu_y, imu.q_wb_zu_z);
        const bool finite =
            std::isfinite(q.w()) && std::isfinite(q.x()) &&
            std::isfinite(q.y()) && std::isfinite(q.z());
        if (!finite || q.norm() < 1e-6f) {
            return Quaternionf::Identity();
        }
        q.normalize();
        return q;
    };

    auto reference_euler_from_csv = [&](const IMU_Sample& imu,
                                        const Quaternionf& q_wb_zu,
                                        float& roll_deg,
                                        float& pitch_deg,
                                        float& yaw_deg) {
        float r_w = 0.0f, p_w = 0.0f, y_w = 0.0f;
        quat_wb_zu_to_euler_nautical(q_wb_zu, r_w, p_w, y_w);
        (void)imu;
        roll_deg = r_w;
        pitch_deg = p_w;
        yaw_deg = y_w;
    };

    WaveDataCSVReader reader(filename);
    reader.for_each_record([&](const Wave_Data_Sample& rec) {
        Vector3f acc_b(rec.imu.acc_bx, rec.imu.acc_by, rec.imu.acc_bz);
        Vector3f gyr_b(rec.imu.gyro_x, rec.imu.gyro_y, rec.imu.gyro_z);

        // Installation physics, on the record's own truth: an IMU at r sees
        // the CG specific force plus the rigid-body rotational terms.
        if (options_.lever_arm) {
            acc_b = options_.lever_arm->install(acc_b, gyr_b);
        }

        if (options_.add_noise) {
            if (noise_models_.accel_noise) {
                acc_b = apply_imu_noise(acc_b, *noise_models_.accel_noise, options_.dt);
            }
            if (noise_models_.gyro_noise) {
                gyr_b = apply_imu_noise(gyr_b, *noise_models_.gyro_noise, options_.dt);
            }
            for (auto& model : noise_models_.extra_imu_noise_models) {
                model(acc_b, gyr_b, options_.dt);
            }
        }

        // Filter-side lever-arm model, on the corrupted signals the estimator
        // actually receives.  Nothing downstream of this point knows the IMU
        // is off the CG.
        if (options_.lever_arm) {
            acc_b = options_.lever_arm->compensate(acc_b, gyr_b);
        }

        const Vector3f acc_meas_ned = zu_to_ned(acc_b);
        const Vector3f gyr_meas_ned = zu_to_ned(gyr_b);

        const Quaternionf q_ref_wb_zu = quat_from_csv(rec.imu);
        const Matrix3f C_wb_zu = q_ref_wb_zu.toRotationMatrix();

        float r_ref_out = 0.0f;
        float p_ref_out = 0.0f;
        float y_ref_out = 0.0f;

        reference_euler_from_csv(rec.imu, q_ref_wb_zu, r_ref_out, p_ref_out, y_ref_out);

        // Vessel heading in the record's own convention, before the magnetic
        // declination offset below moves the yaw reference into the frame the
        // filter actually learns.  Travel-sense scoring needs the physical
        // heading so the boat frame can be removed from the directed estimate.
        const float heading_ref_deg = y_ref_out;

        if (options_.with_mag) {
            y_ref_out = wrapDeg(y_ref_out + MagSim_WMM::default_declination_deg);
        }

        fusion_adapter_.update(options_.dt, gyr_meas_ned, acc_meas_ned, options_.temperature_c);

        if (options_.with_mag) {
            mag_phase_s += options_.dt;
            bool mag_tick = false;
            if (mag_phase_s >= mag_dt) {
                while (mag_phase_s >= mag_dt) mag_phase_s -= mag_dt;
                mag_tick = true;
            }
            if (mag_tick) {
                Vector3f mag_b_enu = C_wb_zu * mag_world_zu;

                if (options_.add_noise && noise_models_.mag_noise) {
                    mag_b_enu = apply_mag_noise(mag_b_enu, *noise_models_.mag_noise, mag_dt);
                }
                if (options_.add_noise) {
                    for (auto& model : noise_models_.extra_mag_noise_models) {
                        model(mag_b_enu, mag_dt);
                    }
                }
                mag_body_ned_hold = zu_to_ned(mag_b_enu);
                fusion_adapter_.updateMag(mag_body_ned_hold);
            }
        }

        TvgNloFilterSnapshot snap = fusion_adapter_.snapshot();

        // Same magnetic-versus-geographic horizontal offset as the W3D runner
        // above; see the comment there.
        if (options_.with_mag) {
            snap.disp_est_zu = MagSim_WMM::mag_frame_to_geographic_nautical(snap.disp_est_zu);
            snap.vel_est_zu  = MagSim_WMM::mag_frame_to_geographic_nautical(snap.vel_est_zu);
            snap.acc_est_zu  = MagSim_WMM::mag_frame_to_geographic_nautical(snap.acc_est_zu);
        }

        result.snapshots.push_back(snap);
        result.final_snapshot = snap;

        Vector3f disp_ref(rec.wave.disp_x, rec.wave.disp_y, rec.wave.disp_z);
        Vector3f vel_ref(rec.wave.vel_x, rec.wave.vel_y, rec.wave.vel_z);
        Vector3f acc_ref(rec.wave.acc_x, rec.wave.acc_y, rec.wave.acc_z);

        Vector3f disp_err = snap.disp_est_zu - disp_ref;
        result.errs_x.push_back(disp_err.x());
        result.errs_y.push_back(disp_err.y());
        result.errs_z.push_back(disp_err.z());
        result.ref_x.push_back(disp_ref.x());
        result.ref_y.push_back(disp_ref.y());
        result.ref_z.push_back(disp_ref.z());
        result.errs_roll.push_back(diffDeg(snap.euler_nautical_deg.x(), r_ref_out));
        result.errs_pitch.push_back(diffDeg(snap.euler_nautical_deg.y(), p_ref_out));
        result.errs_yaw.push_back(diffDeg(snap.euler_nautical_deg.z(), y_ref_out));

        const Vector3f acc_bias_true_zu = (options_.add_noise && noise_models_.accel_noise)
            ? (noise_models_.accel_noise->bias0 + noise_models_.accel_noise->bias_rw).eval()
            : Vector3f::Zero().eval();
        const Vector3f gyro_bias_true_zu = (options_.add_noise && noise_models_.gyro_noise)
            ? (noise_models_.gyro_noise->bias0 + noise_models_.gyro_noise->bias_rw).eval()
            : Vector3f::Zero().eval();
        const Vector3f acc_bias_true_ned = zu_to_ned(acc_bias_true_zu);
        const Vector3f gyro_bias_true_ned = zu_to_ned(gyro_bias_true_zu);

        const Vector3f acc_bias_err = snap.acc_bias_est_ned - acc_bias_true_ned;
        const Vector3f gyro_bias_err = snap.gyro_bias_est_ned - gyro_bias_true_ned;

        const Vector3f mag_bias_true_zu = (options_.add_noise && options_.with_mag && noise_models_.mag_noise)
            ? (noise_models_.mag_noise->bias0_uT + noise_models_.mag_noise->bias_rw_uT).eval()
            : Vector3f::Zero().eval();
        const Vector3f mag_bias_true_ned = zu_to_ned(mag_bias_true_zu);
        const Vector3f mag_bias_err = snap.mag_bias_est_ned_uT - mag_bias_true_ned;

        result.accb_err_x.push_back(acc_bias_err.x());
        result.accb_err_y.push_back(acc_bias_err.y());
        result.accb_err_z.push_back(acc_bias_err.z());
        result.gyrb_err_x.push_back(gyro_bias_err.x());
        result.gyrb_err_y.push_back(gyro_bias_err.y());
        result.gyrb_err_z.push_back(gyro_bias_err.z());
        result.magb_err_x.push_back(mag_bias_err.x());
        result.magb_err_y.push_back(mag_bias_err.y());
        result.magb_err_z.push_back(mag_bias_err.z());

        result.accb_true_x.push_back(acc_bias_true_ned.x());
        result.accb_true_y.push_back(acc_bias_true_ned.y());
        result.accb_true_z.push_back(acc_bias_true_ned.z());
        result.gyrb_true_x.push_back(gyro_bias_true_ned.x());
        result.gyrb_true_y.push_back(gyro_bias_true_ned.y());
        result.gyrb_true_z.push_back(gyro_bias_true_ned.z());
        result.magb_true_x.push_back(mag_bias_true_ned.x());
        result.magb_true_y.push_back(mag_bias_true_ned.y());
        result.magb_true_z.push_back(mag_bias_true_ned.z());

        result.freq_hist.push_back(NAN);
        result.dir_phase_hist.push_back(snap.direction.phase);
        // The axis arrives in the boat frame; adding the vessel heading puts it
        // in the generator frame the record azimuth lives in.  Every shipped
        // record has heading 0, which is why this was invisible until a
        // heading-rotated record was scored.
        result.dir_deg_hist.push_back(wrapAxialDeg90(
            snap.direction.direction_deg_generator_signed + heading_ref_deg));
        result.dir_unc_hist.push_back(snap.direction.uncertainty_deg);
        result.dir_conf_hist.push_back(snap.direction.confidence);
        result.dir_amp_hist.push_back(snap.direction.amplitude);
        result.dir_sign_num_hist.push_back(snap.direction.sign_num);
        result.dir_travel_deg_hist.push_back(
            travelDegGeneratorFromVec(snap.direction.travel_vec_boat, heading_ref_deg));

        ofs << rec.time << ","
            << r_ref_out << "," << p_ref_out << "," << y_ref_out << ","
            << disp_ref.x() << "," << disp_ref.y() << "," << disp_ref.z() << ","
            << vel_ref.x() << "," << vel_ref.y() << "," << vel_ref.z() << ","
            << acc_ref.x() << "," << acc_ref.y() << "," << acc_ref.z() << ","
            << snap.euler_nautical_deg.x() << "," << snap.euler_nautical_deg.y() << "," << snap.euler_nautical_deg.z() << ","
            << snap.disp_est_zu.x() << "," << snap.disp_est_zu.y() << "," << snap.disp_est_zu.z() << ","
            << snap.vel_est_zu.x() << "," << snap.vel_est_zu.y() << "," << snap.vel_est_zu.z() << ","
            << snap.acc_est_zu.x() << "," << snap.acc_est_zu.y() << "," << snap.acc_est_zu.z() << ","
            << acc_bias_true_ned.x() << "," << acc_bias_true_ned.y() << "," << acc_bias_true_ned.z() << ","
            << gyro_bias_true_ned.x() << "," << gyro_bias_true_ned.y() << "," << gyro_bias_true_ned.z() << ","
            << snap.acc_bias_est_ned.x() << "," << snap.acc_bias_est_ned.y() << "," << snap.acc_bias_est_ned.z() << ","
            << snap.gyro_bias_est_ned.x() << "," << snap.gyro_bias_est_ned.y() << "," << snap.gyro_bias_est_ned.z() << ","
            << mag_bias_true_ned.x() << "," << mag_bias_true_ned.y() << "," << mag_bias_true_ned.z() << ","
            << snap.mag_bias_est_ned_uT.x() << "," << snap.mag_bias_est_ned_uT.y() << "," << snap.mag_bias_est_ned_uT.z() << ","
            << mag_bias_err.x() << "," << mag_bias_err.y() << "," << mag_bias_err.z() << ","
            << snap.direction.phase << "," << snap.direction.direction_deg << ","
            << snap.direction.apparent_to_deg << "," << snap.direction.apparent_from_deg << ","
            << snap.direction.sense_coherence << "," << snap.direction.uncertainty_deg << ","
            << snap.direction.confidence << "," << snap.direction.amplitude << ","
            << (snap.direction.sign == FORWARD ? "POSITIVE_AXIS" : snap.direction.sign == BACKWARD ? "NEGATIVE_AXIS" : "UNCERTAIN") << ","
            << snap.direction.sign_num << ","
            << snap.direction.direction_vec.x() << "," << snap.direction.direction_vec.y() << ","
            << snap.direction.filtered_signal.x() << "," << snap.direction.filtered_signal.y();

        write_tvg_nlo_csv_row(ofs, snap);
        ofs << "\n";

        result.final_tau_target = NAN;
        result.final_sigma_target = NAN;
        result.final_tuning_target = NAN;
        result.final_tau_applied = NAN;
        result.final_sigma_applied = NAN;
        result.final_tuning_applied = NAN;
        result.final_wave_period_sec = NAN;
        result.final_freq_hz = NAN;
        result.final_period_sec = NAN;
        result.final_accel_variance = NAN;
    });

    ofs.close();
    std::cout << "Wrote " << result.output_name << "\n";
    return result;
}

std::vector<std::string> collect_wave_data_files(const std::filesystem::path& directory)
{
    std::vector<std::string> files;
    for (auto& entry : std::filesystem::directory_iterator(directory)) {
        if (!entry.is_regular_file()) continue;
        std::string fname = entry.path().string();
        if (fname.find("wave_data_") == std::string::npos) continue;
        if (auto kind = WaveFileNaming::parse_kind_only(fname); kind && *kind == FileKind::Data) {
            files.push_back(std::move(fname));
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

namespace {

// Absolute-time segments requested through W3D_VALIDATION_SEGMENTS, given as
// "name:t0:t1,name:t0:t1,..." with t0/t1 in seconds from the start of the
// record.  The transition case needs these because a single trailing window
// mixes the pure start sea, the crossfade, and the pure endpoint sea, and one
// normalized number over that mixture is not comparable with a stationary
// score.
struct NamedSegment {
    std::string name;
    float start_sec = 0.0f;
    float end_sec = 0.0f;
};

std::vector<NamedSegment> parse_validation_segments()
{
    std::vector<NamedSegment> segments;
    const char* raw = std::getenv("W3D_VALIDATION_SEGMENTS");
    if (!raw || *raw == '\0') return segments;

    std::stringstream spec(raw);
    std::string item;
    while (std::getline(spec, item, ',')) {
        if (item.empty()) continue;
        const size_t first = item.find(':');
        const size_t second = (first == std::string::npos)
            ? std::string::npos
            : item.find(':', first + 1);
        if (second == std::string::npos) {
            throw std::runtime_error(
                "W3D_VALIDATION_SEGMENTS entries must be name:t0:t1");
        }
        NamedSegment segment;
        segment.name = item.substr(0, first);
        segment.start_sec = std::stof(item.substr(first + 1, second - first - 1));
        segment.end_sec = std::stof(item.substr(second + 1));
        if (!(segment.end_sec > segment.start_sec) || segment.start_sec < 0.0f) {
            throw std::runtime_error(
                "W3D_VALIDATION_SEGMENTS requires 0 <= t0 < t1");
        }
        segments.push_back(std::move(segment));
    }
    return segments;
}

// Axial (mod 180 deg) difference, reduced to [-90, 90].
float axial_difference_deg(float estimated_deg, float truth_deg)
{
    return wrapAxialDeg90(estimated_deg - truth_deg);
}

void emit_window_metrics(const W3dSimulationRunResult& result,
                         float dt,
                         size_t start,
                         size_t stop,
                         const char* family,
                         const char* record_tag,
                         const std::string& segment_name)
{
    if (stop <= start || stop > result.errs_z.size()) return;
    const size_t count = stop - start;

    RMSReport x, y, z, roll, pitch, yaw;
    RMSReport acc_x, acc_y, acc_z, gyro_x, gyro_y, gyro_z;
    RMSReport ref_z_rms;
    // Mean error alongside the RMS.  An RMS cannot tell a zero-mean
    // fluctuation from a constant offset, and several deployed error sources
    // -- attitude rectification under vibration among them -- show up almost
    // entirely as an offset.  The angles use a circular mean so a yaw error
    // near the +/-180 wrap does not average to zero.
    double mean_x = 0.0, mean_y = 0.0, mean_z = 0.0;
    double roll_c = 0.0, roll_s = 0.0;
    double pitch_c = 0.0, pitch_s = 0.0;
    double yaw_c = 0.0, yaw_s = 0.0;
    float ref_max_3d = 0.0f;
    for (size_t i = start; i < stop; ++i) {
        ref_z_rms.add(result.ref_z[i]);
        x.add(result.errs_x[i]);
        y.add(result.errs_y[i]);
        z.add(result.errs_z[i]);
        roll.add(result.errs_roll[i]);
        pitch.add(result.errs_pitch[i]);
        yaw.add(result.errs_yaw[i]);
        mean_x += static_cast<double>(result.errs_x[i]);
        mean_y += static_cast<double>(result.errs_y[i]);
        mean_z += static_cast<double>(result.errs_z[i]);
        const double roll_rad = static_cast<double>(deg_to_rad(result.errs_roll[i]));
        const double pitch_rad = static_cast<double>(deg_to_rad(result.errs_pitch[i]));
        const double yaw_rad = static_cast<double>(deg_to_rad(result.errs_yaw[i]));
        roll_c += std::cos(roll_rad);
        roll_s += std::sin(roll_rad);
        pitch_c += std::cos(pitch_rad);
        pitch_s += std::sin(pitch_rad);
        yaw_c += std::cos(yaw_rad);
        yaw_s += std::sin(yaw_rad);
        acc_x.add(result.accb_err_x[i]);
        acc_y.add(result.accb_err_y[i]);
        acc_z.add(result.accb_err_z[i]);
        gyro_x.add(result.gyrb_err_x[i]);
        gyro_y.add(result.gyrb_err_y[i]);
        gyro_z.add(result.gyrb_err_z[i]);
        const float rx = result.ref_x[i];
        const float ry = result.ref_y[i];
        const float rz = result.ref_z[i];
        ref_max_3d = std::max(ref_max_3d, std::sqrt(rx * rx + ry * ry + rz * rz));
    }

    const float x_rms = x.rms();
    const float y_rms = y.rms();
    const float z_rms = z.rms();
    const double inv_count = 1.0 / static_cast<double>(count);
    const float x_mean = static_cast<float>(mean_x * inv_count);
    const float y_mean = static_cast<float>(mean_y * inv_count);
    const float z_mean = static_cast<float>(mean_z * inv_count);
    auto circular_mean_deg = [](double cos_sum, double sin_sum) -> float {
        if (std::fabs(cos_sum) < 1e-12 && std::fabs(sin_sum) < 1e-12) return NAN;
        return rad_to_deg(std::atan2(sin_sum, cos_sum));
    };
    const float roll_mean = circular_mean_deg(roll_c, roll_s);
    const float pitch_mean = circular_mean_deg(pitch_c, pitch_s);
    const float yaw_mean = circular_mean_deg(yaw_c, yaw_s);
    const float disp_3d = std::sqrt(
        x_rms * x_rms + y_rms * y_rms + z_rms * z_rms);
    const float acc_3d = std::sqrt(
        acc_x.rms() * acc_x.rms() +
        acc_y.rms() * acc_y.rms() +
        acc_z.rms() * acc_z.rms());
    const float gyro_3d = std::sqrt(
        gyro_x.rms() * gyro_x.rms() +
        gyro_y.rms() * gyro_y.rms() +
        gyro_z.rms() * gyro_z.rms());
    const float hs = result.wave_params.height;
    const float pct_hs = (hs > 0.0f) ? 100.0f / hs : NAN;
    const float disp_3d_pct_refmax = (ref_max_3d > 1e-12f)
        ? 100.0f * disp_3d / ref_max_3d
        : NAN;

    // Reference-RMS normalization.  Unlike %H_s it uses the motion actually
    // present in the scored interval, so it stays meaningful when the sea
    // state changes inside the window and is comparable between segments of
    // the transition record.
    const float ref_z = ref_z_rms.rms();
    const float disp_z_pct_refrms = (ref_z > 1e-9f) ? 100.0f * z_rms / ref_z : NAN;

    // Direction accuracy over the same window, against the generator azimuth
    // recovered from the record.  dir_deg_hist is already expressed in the
    // generator convention, and the axis is defined mod 180 deg.
    const float truth_deg = result.wave_params.direction;
    float dir_axis_mean_deg = NAN;
    float dir_axis_error_deg = NAN;
    float dir_axis_rmse_deg = NAN;
    float dir_axis_circ_std_deg = NAN;
    float sense_forward_pct = NAN;
    float sense_reverse_pct = NAN;
    float sense_uncertain_pct = NAN;
    float sense_dominant_pct = NAN;
    if (stop <= result.dir_deg_hist.size() && stop <= result.dir_sign_num_hist.size()) {
        std::vector<float> axis;
        axis.reserve(count);
        RMSReport axis_error;
        for (size_t i = start; i < stop; ++i) {
            const float d = result.dir_deg_hist[i];
            if (!std::isfinite(d)) continue;
            axis.push_back(d);
            axis_error.add(axial_difference_deg(d, truth_deg));
        }
        if (!axis.empty()) {
            const auto stats = circular_stats_180(axis);
            dir_axis_mean_deg = stats.mean_deg;
            dir_axis_circ_std_deg = stats.std_deg;
            dir_axis_error_deg = axial_difference_deg(stats.mean_deg, truth_deg);
            dir_axis_rmse_deg = axis_error.rms();
        }

        size_t forward = 0, reverse = 0, uncertain = 0;
        for (size_t i = start; i < stop; ++i) {
            const int sign = result.dir_sign_num_hist[i];
            if (sign > 0) ++forward;
            else if (sign < 0) ++reverse;
            else ++uncertain;
        }
        const float scale = 100.0f / static_cast<float>(count);
        sense_forward_pct = scale * static_cast<float>(forward);
        sense_reverse_pct = scale * static_cast<float>(reverse);
        sense_uncertain_pct = scale * static_cast<float>(uncertain);

        // FORWARD/BACKWARD are relative to the axis representative the axis
        // estimator happens to return, and that representative is not tied to
        // the generator azimuth: records of opposite nominal heading put the
        // same physical travel sense in opposite classes.  The share of
        // samples in whichever class dominates therefore measures how
        // consistently the classifier commits to one sense, which is what the
        // raw counts above can support; it is not a correctness rate against
        // the generator, and must not be reported as one.
        sense_dominant_pct = std::max(sense_forward_pct, sense_reverse_pct);
    }

    // Travel-sense correctness. Unlike the class shares above this scores the
    // directed propagation vector, with the vessel heading removed, against the
    // physical propagation direction of the record, so it is a genuine error
    // rate and is invariant to which end of the axis the estimator returns.
    float travel_error_deg = NAN;
    float travel_rmse_deg = NAN;
    float travel_correct_pct = NAN;
    float travel_wrong_pct = NAN;
    float travel_unresolved_pct = NAN;
    if (stop <= result.dir_travel_deg_hist.size()) {
        const float travel_truth_deg =
            travelTruthDegFromGeneratorAzimuth(truth_deg);
        double sum_sin = 0.0;
        double sum_cos = 0.0;
        RMSReport travel_error;
        size_t resolved = 0;
        size_t correct = 0;
        for (size_t i = start; i < stop; ++i) {
            const float estimate = result.dir_travel_deg_hist[i];
            if (!std::isfinite(estimate)) continue;
            ++resolved;
            const float error = wrapDeg(estimate - travel_truth_deg);
            sum_sin += std::sin(static_cast<double>(deg_to_rad(error)));
            sum_cos += std::cos(static_cast<double>(deg_to_rad(error)));
            travel_error.add(error);
            if (std::abs(error) < 90.0f) ++correct;
        }
        const float scale = 100.0f / static_cast<float>(count);
        travel_correct_pct = scale * static_cast<float>(correct);
        travel_wrong_pct = scale * static_cast<float>(resolved - correct);
        travel_unresolved_pct = scale * static_cast<float>(count - resolved);
        if (resolved > 0) {
            travel_error_deg = rad_to_deg(std::atan2(sum_sin, sum_cos));
            travel_rmse_deg = travel_error.rms();
        }
    }

    const char* mode = std::getenv("W3D_TUNING_MODE");
    if (!mode || *mode == '\0') mode = "adaptive";
    const char* aw_cov_sync = std::getenv("W3D_AW_COV_SYNC");
    if (!aw_cov_sync || *aw_cov_sync == '\0') aw_cov_sync = "periodic";

    const auto old_precision = std::cout.precision();
    std::cout << std::setprecision(9)
              << record_tag
              << " family=" << family
              << " tuning_mode=" << mode
              << " aw_cov_sync=" << aw_cov_sync
              << " input=" << std::filesystem::path(result.input_name).filename().string()
              << " segment=" << (segment_name.empty() ? "window" : segment_name)
              << " start_s=" << (static_cast<float>(start) * dt)
              << " window_s=" << (static_cast<float>(count) * dt)
              << " samples=" << count
              << " disp_x_rms_m=" << x_rms
              << " disp_y_rms_m=" << y_rms
              << " disp_z_rms_m=" << z_rms
              << " disp_3d_rms_m=" << disp_3d
              << " disp_x_pct_hs=" << x_rms * pct_hs
              << " disp_y_pct_hs=" << y_rms * pct_hs
              << " disp_z_pct_hs=" << z_rms * pct_hs
              << " disp_3d_pct_refmax=" << disp_3d_pct_refmax
              << " disp_z_ref_rms_m=" << ref_z
              << " disp_z_pct_refrms=" << disp_z_pct_refrms
              << " roll_rms_deg=" << roll.rms()
              << " pitch_rms_deg=" << pitch.rms()
              << " yaw_rms_deg=" << yaw.rms()
              << " dir_travel_error_deg=" << travel_error_deg
              << " dir_travel_rmse_deg=" << travel_rmse_deg
              << " dir_travel_correct_pct=" << travel_correct_pct
              << " dir_travel_wrong_pct=" << travel_wrong_pct
              << " dir_travel_unresolved_pct=" << travel_unresolved_pct
              << " dir_axis_mean_deg=" << dir_axis_mean_deg
              << " dir_axis_truth_deg=" << truth_deg
              << " dir_axis_error_deg=" << dir_axis_error_deg
              << " dir_axis_rmse_deg=" << dir_axis_rmse_deg
              << " dir_axis_circ_std_deg=" << dir_axis_circ_std_deg
              << " dir_sense_forward_pct=" << sense_forward_pct
              << " dir_sense_reverse_pct=" << sense_reverse_pct
              << " dir_sense_uncertain_pct=" << sense_uncertain_pct
              << " dir_sense_dominant_pct=" << sense_dominant_pct
              << " accel_bias_3d_rms_mps2=" << acc_3d
              << " gyro_bias_3d_rms_radps=" << gyro_3d
              << " tau_applied_s=" << result.final_tau_applied
              << " sigma_applied_mps2=" << result.final_sigma_applied
              << " tuning_applied=" << result.final_tuning_applied
              << " wave_period_s=" << result.final_wave_period_sec
              << " frequency_hz=" << result.final_freq_hz
              << " period_s=" << result.final_period_sec
              << " accel_variance_m2ps4=" << result.final_accel_variance
              << "\n";

    // Mean error on its own line.  An RMS cannot separate a zero-mean
    // fluctuation from a constant offset, and some deployed error sources --
    // attitude rectification under vibration among them -- are almost purely
    // an offset.  This is deliberately *not* folded into VALIDATION_METRICS:
    // that line is parsed into the committed validation and robustness
    // evidence, and adding columns to it would rewrite files the evidence
    // fingerprint covers without any of their numbers having changed.
    std::cout << record_tag << "_MEANS"
              << " family=" << family
              << " input=" << std::filesystem::path(result.input_name).filename().string()
              << " segment=" << (segment_name.empty() ? "window" : segment_name)
              << " samples=" << count
              << " disp_x_mean_m=" << x_mean
              << " disp_y_mean_m=" << y_mean
              << " disp_z_mean_m=" << z_mean
              << " roll_mean_deg=" << roll_mean
              << " pitch_mean_deg=" << pitch_mean
              << " yaw_mean_deg=" << yaw_mean
              << "\n";
    std::cout.precision(old_precision);
}

}  // namespace

void print_validation_metrics(const W3dSimulationRunResult& result,
                              float dt,
                              float window_seconds,
                              const char* family)
{
    if (!(dt > 0.0f) || !(window_seconds > 0.0f) || result.errs_z.empty()) return;

    const size_t total = result.errs_z.size();
    const size_t requested = static_cast<size_t>(window_seconds / dt);
    const size_t count = std::min(total, std::max<size_t>(requested, 1u));
    emit_window_metrics(
        result, dt, total - count, total, family, "VALIDATION_METRICS", "");

    for (const auto& segment : parse_validation_segments()) {
        const size_t i0 = std::min(
            total, static_cast<size_t>(segment.start_sec / dt));
        const size_t i1 = std::min(
            total, static_cast<size_t>(segment.end_sec / dt));
        emit_window_metrics(
            result, dt, i0, i1, family, "VALIDATION_SEGMENT", segment.name);
    }
}

void print_summary_and_fail_if_needed(const W3dSimulationRunResult& result,
                                      float dt,
                                      const W3dFailureLimits& limits,
                                      const W3dSummaryLabels& labels)
{
    // The scored trailing window for the OU-II/OU-III simulators.  It matches
    // the window the statistical validation runner scores, so the executable
    // gates and the ensemble study describe the same stretch of a 20-minute
    // replay instead of two windows an order of magnitude apart.
    constexpr float RMS_WINDOW_SEC = 900.0f;
    constexpr int RMS_WINDOW_SEC_LABEL = static_cast<int>(RMS_WINDOW_SEC);
    const int N_last = static_cast<int>(RMS_WINDOW_SEC / dt);
    // A partial window is not the scored window, so a record shorter than it is
    // left ungated rather than scored against sentinels fitted to the full one.
    // Said out loud, because a caller that scrapes QUALITY_GATE otherwise reads
    // the silence as a failure.
    if (result.errs_z.size() <= static_cast<size_t>(N_last)) {
        std::cout << "QUALITY_GATE: SKIPPED REASON=record_shorter_than_"
                  << RMS_WINDOW_SEC_LABEL << "s_window RECORD="
                  << result.output_name << "\n";
        return;
    }

    const size_t start = result.errs_z.size() - N_last;
    RMSReport rms_x, rms_y, rms_z, rms_roll, rms_pitch, rms_yaw;
    RMSReport rms_accb_x, rms_accb_y, rms_accb_z;
    RMSReport rms_gyrb_x, rms_gyrb_y, rms_gyrb_z;
    RMSReport rms_magb_x, rms_magb_y, rms_magb_z;

    float acc_true_max_x = 0.f, acc_true_max_y = 0.f, acc_true_max_z = 0.f, acc_true_max_3d = 0.f;
    float gyr_true_max_x = 0.f, gyr_true_max_y = 0.f, gyr_true_max_z = 0.f, gyr_true_max_3d = 0.f;
    float mag_true_max_x = 0.f, mag_true_max_y = 0.f, mag_true_max_z = 0.f, mag_true_max_3d = 0.f;
    float disp_true_max_3d = 0.f;

    for (size_t i = start; i < result.errs_z.size(); ++i) {
        rms_x.add(result.errs_x[i]);
        rms_y.add(result.errs_y[i]);
        rms_z.add(result.errs_z[i]);
        rms_roll.add(result.errs_roll[i]);
        rms_pitch.add(result.errs_pitch[i]);
        rms_yaw.add(result.errs_yaw[i]);
        rms_accb_x.add(result.accb_err_x[i]);
        rms_accb_y.add(result.accb_err_y[i]);
        rms_accb_z.add(result.accb_err_z[i]);
        rms_gyrb_x.add(result.gyrb_err_x[i]);
        rms_gyrb_y.add(result.gyrb_err_y[i]);
        rms_gyrb_z.add(result.gyrb_err_z[i]);
        rms_magb_x.add(result.magb_err_x[i]);
        rms_magb_y.add(result.magb_err_y[i]);
        rms_magb_z.add(result.magb_err_z[i]);

        const float dx = result.ref_x[i];
        const float dy = result.ref_y[i];
        const float dz = result.ref_z[i];
        disp_true_max_3d = std::max(disp_true_max_3d, std::sqrt(dx * dx + dy * dy + dz * dz));

        const float ax = result.accb_true_x[i], ay = result.accb_true_y[i], az = result.accb_true_z[i];
        acc_true_max_x = std::max(acc_true_max_x, std::abs(ax));
        acc_true_max_y = std::max(acc_true_max_y, std::abs(ay));
        acc_true_max_z = std::max(acc_true_max_z, std::abs(az));
        acc_true_max_3d = std::max(acc_true_max_3d, std::sqrt(ax * ax + ay * ay + az * az));

        const float gx = result.gyrb_true_x[i], gy = result.gyrb_true_y[i], gz = result.gyrb_true_z[i];
        gyr_true_max_x = std::max(gyr_true_max_x, std::abs(gx));
        gyr_true_max_y = std::max(gyr_true_max_y, std::abs(gy));
        gyr_true_max_z = std::max(gyr_true_max_z, std::abs(gz));
        gyr_true_max_3d = std::max(gyr_true_max_3d, std::sqrt(gx * gx + gy * gy + gz * gz));

        const float mx = result.magb_true_x[i], my = result.magb_true_y[i], mz = result.magb_true_z[i];
        mag_true_max_x = std::max(mag_true_max_x, std::abs(mx));
        mag_true_max_y = std::max(mag_true_max_y, std::abs(my));
        mag_true_max_z = std::max(mag_true_max_z, std::abs(mz));
        mag_true_max_3d = std::max(mag_true_max_3d, std::sqrt(mx * mx + my * my + mz * mz));
    }

    const float x_rms = rms_x.rms(), y_rms = rms_y.rms(), z_rms = rms_z.rms();
    const float x_pct = 100.f * x_rms / result.wave_params.height;
    const float y_pct = 100.f * y_rms / result.wave_params.height;
    const float z_pct = 100.f * z_rms / result.wave_params.height;
    const float rms_3d_err = std::sqrt(x_rms * x_rms + y_rms * y_rms + z_rms * z_rms);
    const float pct_3d = (disp_true_max_3d > 1e-12f && std::isfinite(rms_3d_err))
        ? 100.f * rms_3d_err / disp_true_max_3d
        : NAN;

    std::cout << "=== Last " << RMS_WINDOW_SEC_LABEL << " s RMS summary for "
              << result.output_name << " ===\n";
    std::cout << "XYZ RMS (m): X=" << x_rms << " Y=" << y_rms << " Z=" << z_rms << "\n";
    std::cout << "XYZ RMS (%Hs): X=" << x_pct << "% Y=" << y_pct << "% Z=" << z_pct
              << "% (Hs=" << result.wave_params.height << ")\n";
    std::cout << "3D RMS (m): " << rms_3d_err
              << " (3D % of max |disp_ref|_3D = " << pct_3d
              << "%, max |disp_ref|_3D = " << disp_true_max_3d << " m)\n";
    std::cout << "Angles RMS (deg): Roll=" << rms_roll.rms()
              << " Pitch=" << rms_pitch.rms()
              << " Yaw=" << rms_yaw.rms() << "\n";

    auto vec_rms = [](float rx, float ry, float rz) { return std::sqrt(rx * rx + ry * ry + rz * rz); };
    const float accb_rx = rms_accb_x.rms(), accb_ry = rms_accb_y.rms(), accb_rz = rms_accb_z.rms();
    const float gyrb_rx = rms_gyrb_x.rms(), gyrb_ry = rms_gyrb_y.rms(), gyrb_rz = rms_gyrb_z.rms();
    const float magb_rx = rms_magb_x.rms(), magb_ry = rms_magb_y.rms(), magb_rz = rms_magb_z.rms();
    const float accb_r3 = vec_rms(accb_rx, accb_ry, accb_rz);
    const float gyrb_r3 = vec_rms(gyrb_rx, gyrb_ry, gyrb_rz);
    const float magb_r3 = vec_rms(magb_rx, magb_ry, magb_rz);

    std::cout << "Bias error RMS (acc, m/s^2): X=" << accb_rx << " Y=" << accb_ry << " Z=" << accb_rz
              << " |3D|=" << accb_r3 << "\n";
    std::cout << "Bias error RMS (gyro, rad/s): X=" << gyrb_rx << " Y=" << gyrb_ry << " Z=" << gyrb_rz
              << " |3D|=" << gyrb_r3 << "\n";

    const float rad2deg = 180.0f / float(std::numbers::pi_v<float>);
    std::cout << "Bias error RMS (gyro, deg/s): X=" << (gyrb_rx * rad2deg)
              << " Y=" << (gyrb_ry * rad2deg)
              << " Z=" << (gyrb_rz * rad2deg)
              << " |3D|=" << (gyrb_r3 * rad2deg) << "\n";
    std::cout << "Bias error RMS (mag, uT): X=" << magb_rx << " Y=" << magb_ry << " Z=" << magb_rz
              << " |3D|=" << magb_r3 << "\n";

    auto pct_of_max = [](float rms, float maxv) -> float {
        return (maxv > 1e-12f && std::isfinite(rms)) ? (100.f * rms / maxv) : NAN;
    };

    std::cout << "Max TRUE bias in window (acc, m/s^2): X=" << acc_true_max_x << " Y=" << acc_true_max_y
              << " Z=" << acc_true_max_z << " |3D|=" << acc_true_max_3d << "\n";
    std::cout << "Max TRUE bias in window (gyro, rad/s): X=" << gyr_true_max_x << " Y=" << gyr_true_max_y
              << " Z=" << gyr_true_max_z << " |3D|=" << gyr_true_max_3d << "\n";

    const float accb_r3_pct = pct_of_max(accb_r3, acc_true_max_3d);
    const float gyrb_r3_pct = pct_of_max(gyrb_r3, gyr_true_max_3d);
    std::cout << "Bias error RMS (% of max TRUE bias) (acc): X=" << pct_of_max(accb_rx, acc_true_max_x)
              << "% Y=" << pct_of_max(accb_ry, acc_true_max_y)
              << "% Z=" << pct_of_max(accb_rz, acc_true_max_z)
              << "% |3D|=" << accb_r3_pct << "%\n";
    std::cout << "Bias error RMS (% of max TRUE bias) (gyro): X=" << pct_of_max(gyrb_rx, gyr_true_max_x)
              << "% Y=" << pct_of_max(gyrb_ry, gyr_true_max_y)
              << "% Z=" << pct_of_max(gyrb_rz, gyr_true_max_z)
              << "% |3D|=" << gyrb_r3_pct << "%\n";

    std::cout << "tau_target=" << result.final_tau_target
              << ", sigma_target=" << result.final_sigma_target
              << ", " << labels.target << "=" << result.final_tuning_target << "\n";
    std::cout << "tau_applied=" << result.final_tau_applied
              << ", sigma_applied=" << result.final_sigma_applied
              << ", " << labels.applied << "=" << result.final_tuning_applied << "\n";
    std::cout << "f_hz=" << result.final_freq_hz
              << ", Tp_tuner=" << result.final_period_sec
              << ", accel_var=" << result.final_accel_variance << "\n";

    if (start < result.dir_deg_hist.size()) {
        const size_t i0 = start;
        const size_t i1 = result.errs_z.size();
        std::vector<float> vf(result.freq_hist.begin() + i0, result.freq_hist.begin() + i1);
        std::vector<float> vd(result.dir_deg_hist.begin() + i0, result.dir_deg_hist.begin() + i1);
        std::vector<float> vu(result.dir_unc_hist.begin() + i0, result.dir_unc_hist.begin() + i1);
        std::vector<float> vc(result.dir_conf_hist.begin() + i0, result.dir_conf_hist.begin() + i1);
        vd.erase(std::remove_if(vd.begin(), vd.end(), [](float a){ return !std::isfinite(a); }), vd.end());
        auto cs = circular_stats_180(vd);

        int nToward = 0, nAway = 0, nUnc = 0;
        size_t good = 0;
        constexpr float CONF_THRESH = 20.0f;
        constexpr float AMP_THRESH = 0.08f;
        for (size_t k = i0; k < i1; ++k) {
            const int s = result.dir_sign_num_hist[k];
            if (s > 0) ++nToward;
            else if (s < 0) ++nAway;
            else ++nUnc;
            if (result.dir_conf_hist[k] > CONF_THRESH && result.dir_amp_hist[k] > AMP_THRESH) ++good;
        }
        const int nWin = int(i1 - i0);
        auto pct = [&](int n){ return (nWin > 0) ? (100.0 * double(n) / double(nWin)) : 0.0; };

        std::cout << "=== Direction Report (last " << RMS_WINDOW_SEC_LABEL
                  << " s only) for " << result.output_name << " ===\n";
        std::cout << "window_s: " << (float(i1 - i0) * dt) << " samples: " << (i1 - i0) << "\n";
        std::cout << "freq_hz: mean=" << mean_vec(vf)
                  << " median=" << median_vec(vf)
                  << " p05=" << percentile_vec(vf, 0.05)
                  << " p95=" << percentile_vec(vf, 0.95) << "\n";
        std::cout << "dir_deg_gen ([-90,90], 0=+Y CW): mean_circ=" << cs.mean_deg
                  << " circ_std≈" << cs.std_deg << " deg\n";
        std::cout << "uncert_deg: mean=" << mean_vec(vu)
                  << " median=" << median_vec(vu)
                  << " p95=" << percentile_vec(vu, 0.95) << "\n";
        std::cout << "confidence: mean=" << mean_vec(vc)
                  << " >" << CONF_THRESH << " count=" << good
                  << " (" << (100.0 * double(good) / double(i1 - i0)) << "%)\n";
        std::cout << "sense: +AXIS=" << nToward << " (" << pct(nToward) << "%)"
                  << " -AXIS=" << nAway << " (" << pct(nAway) << "%)"
                  << " UNCERTAIN=" << nUnc << " (" << pct(nUnc) << "%)\n";
        std::cout << "=============================================\n\n";
    }

    const float limit_z = (result.wave_type == WaveType::JONSWAP)
        ? limits.err_limit_percent_z_jonswap
        : limits.err_limit_percent_z_pmstokes;
    const float limit_3d = (result.wave_type == WaveType::JONSWAP)
        ? limits.err_limit_percent_3d_jonswap
        : limits.err_limit_percent_3d_pmstokes;
    const bool collect_all_gates = (std::getenv("W3D_COLLECT_ALL_GATES") != nullptr);
    bool failed = false;
    std::string fail_reason = "ok";
    auto fail_if = [&](const char* label, float pct, float limit) {
        if (pct > limit) {
            failed = true;
            fail_reason = std::string(label) + "_limit_exceeded";
            std::cerr << "ERROR: " << label << " RMS above limit (" << pct << "% > " << limit << "%). Failing.\n";
            if (!collect_all_gates) std::exit(EXIT_FAILURE);
        }
    };

    fail_if("Z", z_pct, limit_z);
    fail_if("3D", pct_3d, limit_3d);

    // The attitude bars are fitted to the deployed protocol, which runs the
    // magnetometer, so an IMU-only run is not scored against them -- the same
    // rule yaw has always followed here, and for the same reason.  Roll and
    // pitch are observable without a magnetometer, but not equally well:
    // worst-case pitch over the eight records rises 18 percent for OU-III and
    // 97 percent for OU-II under --nomag, while worst-case roll rises 13
    // percent for OU-II and *falls* 19 percent for OU-III.  A bar fitted with
    // the magnetometer is measuring a different filter without it.
    auto fail_attitude_if = [&](const char* label, float deg, float limit) {
        if (!result.with_mag || !(limit > 0.0f)) return;
        if (deg > limit) {
            failed = true;
            fail_reason = std::string(label) + "_limit_exceeded";
            std::cerr << "ERROR: " << label << " RMS above limit (" << deg
                      << " deg > " << limit << " deg). Failing.\n";
            if (!collect_all_gates) std::exit(EXIT_FAILURE);
        }
    };

    if (result.with_mag && rms_yaw.rms() > limits.err_limit_yaw_deg) {
        failed = true;
        fail_reason = "yaw_limit_exceeded";
        std::cerr << "ERROR: Yaw RMS above limit (" << rms_yaw.rms() << " deg > "
                  << limits.err_limit_yaw_deg << " deg). Failing.\n";
        if (!collect_all_gates) std::exit(EXIT_FAILURE);
    }
    fail_attitude_if("Roll", rms_roll.rms(), limits.err_limit_roll_deg);
    fail_attitude_if("Pitch", rms_pitch.rms(), limits.err_limit_pitch_deg);

    const float accb_z_pct = pct_of_max(accb_rz, acc_true_max_z);
    if (std::isfinite(accb_z_pct) && accb_z_pct > limits.acc_z_bias_percent) {
        failed = true;
        fail_reason = "acc_z_bias_percent_exceeded";
        std::cerr << "ERROR: accel Z bias error RMS above limit ("
                  << accb_z_pct << "% > " << limits.acc_z_bias_percent
                  << "% of max TRUE Z bias). Failing.\n";
        if (!collect_all_gates) std::exit(EXIT_FAILURE);
    }
    if (std::isfinite(accb_r3_pct) && accb_r3_pct > limits.bias_3d_percent) {
        failed = true;
        fail_reason = "acc_bias_3d_percent_exceeded";
        std::cerr << "ERROR: 3D accel bias error RMS above limit ("
                  << accb_r3_pct << "% > " << limits.bias_3d_percent
                  << "% of max TRUE bias). Failing.\n";
        if (!collect_all_gates) std::exit(EXIT_FAILURE);
    }
    const float gyro_bias_limit = (limits.gyro_bias_3d_percent > 0.0f)
        ? limits.gyro_bias_3d_percent
        : limits.bias_3d_percent;
    if (std::isfinite(gyrb_r3_pct) && gyrb_r3_pct > gyro_bias_limit) {
        failed = true;
        fail_reason = "gyro_bias_3d_percent_exceeded";
        std::cerr << "ERROR: 3D gyro bias error RMS above limit ("
                  << gyrb_r3_pct << "% > " << gyro_bias_limit
                  << "% of max TRUE bias). Failing.\n";
        if (!collect_all_gates) std::exit(EXIT_FAILURE);
    }
    if (failed) g_w3d_any_gate_failed = true;
    std::cout << "QUALITY_GATE: PASS=" << (failed ? 0 : 1)
              << " REASON=" << fail_reason << "\n";
}

bool w3d_any_quality_gate_failed()
{
    return g_w3d_any_gate_failed;
}

#endif
