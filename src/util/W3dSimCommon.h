#pragma once

/*
  Copyright 2025-2026, Mikhail Grushinskiy
*/

#ifdef ARDUINO
#else

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <numbers>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "util/WaveFilesSupport.h"
#include "ahrs/FrameConversions.h"
#include "wave_dir/WaveDirectionDetector.h"

using Eigen::Vector2f;
using Eigen::Vector3f;
using Eigen::Matrix3f;

inline float wrapDeg(float a) {
    a = std::fmod(a + 180.0f, 360.0f);
    if (a < 0) a += 360.0f;
    return a - 180.0f;
}

inline float diffDeg(float est_deg, float ref_deg) {
    return wrapDeg(est_deg - ref_deg);
}

static inline float deg_to_rad(float d) { return d * (std::numbers::pi_v<float> / 180.0f); }
static inline float rad_to_deg(float r) { return r * (180.0f / std::numbers::pi_v<float>); }
static inline float rad_to_deg(double r) { return static_cast<float>(r * (180.0 / std::numbers::pi)); }

inline float wrapAxialDeg90(float a) {
    a = std::fmod(a + 180.0f, 360.0f);
    if (a < 0) a += 360.0f;
    a -= 180.0f;
    if (a > 90.0f) a -= 180.0f;
    if (a < -90.0f) a += 180.0f;
    return a;
}

inline float dirDegGeneratorSignedFromVec(const Vector2f& v) {
    float deg = rad_to_deg(std::atan2(v.x(), v.y()));
    return wrapAxialDeg90(deg);
}

// Directed counterpart of dirDegGeneratorSignedFromVec, in [0, 360).
//
// heading_deg is the vessel heading in the record's own convention, which is 0
// for every shipped record; adding it removes the boat frame so the result is
// directly comparable with the record azimuth.  A zero vector means the travel
// sense was not resolved for this sample and yields NaN.
inline float travelDegGeneratorFromVec(const Vector2f& v, float heading_deg) {
    if (!(v.squaredNorm() > 1e-12f) || !std::isfinite(heading_deg)) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    float deg = rad_to_deg(std::atan2(v.x(), v.y())) + heading_deg;
    deg = std::fmod(deg, 360.0f);
    if (deg < 0.0f) deg += 360.0f;
    return deg;
}

// The generator azimuth in a record name is the direction the waves come from:
// the propagation-to vector recovered from the truth channels lies at
// azimuth + 180 in every shipped record.  Travel-sense scoring compares against
// this, not against the raw azimuth.
inline float travelTruthDegFromGeneratorAzimuth(float azimuth_deg) {
    float deg = std::fmod(azimuth_deg + 180.0f, 360.0f);
    if (deg < 0.0f) deg += 360.0f;
    return deg;
}

inline float p0_s_from_sigma_tau(float sigma_a, float tau) {
    if (!std::isfinite(sigma_a) || !std::isfinite(tau)) return NAN;
    return sigma_a * tau * tau;
}

template<typename T>
inline T mean_vec(const std::vector<T>& v) {
    if (v.empty()) return T(NAN);
    T s = 0;
    for (const auto& x : v) s += x;
    return s / T(v.size());
}

template<typename T>
inline T median_vec(std::vector<T> v) {
    if (v.empty()) return T(NAN);
    size_t n = v.size();
    std::nth_element(v.begin(), v.begin() + n / 2, v.end());
    if (n % 2) return v[n / 2];
    auto lo = *std::max_element(v.begin(), v.begin() + n / 2);
    auto hi = v[n / 2];
    return (lo + hi) / T(2);
}

template<typename T>
inline T percentile_vec(std::vector<T> v, double p01) {
    if (v.empty()) return T(NAN);
    if (p01 <= 0) return *std::min_element(v.begin(), v.end());
    if (p01 >= 1) return *std::max_element(v.begin(), v.end());
    std::sort(v.begin(), v.end());
    double idx = p01 * static_cast<double>(v.size() - 1);
    size_t i = size_t(std::floor(idx));
    double frac = idx - double(i);
    if (i + 1 >= v.size()) return v[i];
    return T(v[i] * (1.0 - frac) + v[i + 1] * frac);
}

struct CircStats {
    float mean_deg = NAN;
    float std_deg = NAN;
};

inline CircStats circular_stats_180(const std::vector<float>& degs) {
    CircStats cs;
    if (degs.empty()) return cs;

    double C = 0, S = 0;
    for (float d : degs) {
        const double a2 = 2.0 * deg_to_rad(d);
        C += std::cos(a2);
        S += std::sin(a2);
    }
    C /= double(degs.size());
    S /= double(degs.size());

    const double R = std::sqrt(C * C + S * S);
    const double a2_mean = std::atan2(S, C);

    float md = float(rad_to_deg(0.5 * a2_mean));
    md = wrapAxialDeg90(md);
    cs.mean_deg = md;
    cs.std_deg = (R > 1e-12)
        ? float(rad_to_deg(0.5 * std::sqrt(std::max(0.0, -2.0 * std::log(R)))))
        : 90.0f;
    cs.std_deg = std::min(cs.std_deg, 90.0f);
    return cs;
}

class RMSReport {
public:
    void add(float value) { sum_sq_ += value * value; count_++; }
    float rms() const { return count_ ? std::sqrt(sum_sq_ / float(count_)) : NAN; }

private:
    float sum_sq_ = 0.0f;
    size_t count_ = 0;
};

extern const float g_std;

struct ImuNoiseModel {
    std::mt19937 rng;
    std::normal_distribution<float> w;
    std::normal_distribution<float> n01;
    Vector3f bias0;
    Vector3f bias_rw;
    float sigma_bias_rw = 0.0f;
};

ImuNoiseModel make_imu_noise_model(float sigma_white,
                                   float bias_half_range,
                                   float sigma_bias_rw,
                                   unsigned seed);

ImuNoiseModel make_imu_noise_model(float sigma_white,
                                   float bias_half_range,
                                   float sigma_bias_rw,
                                   unsigned noise_seed,
                                   unsigned initialization_seed);

Vector3f apply_imu_noise(const Vector3f& truth, ImuNoiseModel& m, float dt);

struct MagNoiseModel {
    std::mt19937 rng;
    std::normal_distribution<float> w_uT;
    std::normal_distribution<float> n01;

    Vector3f bias0_uT;
    Vector3f bias_rw_uT;
    float sigma_bias_rw_uT_sqrt_s = 0.0f;

    Eigen::Matrix3f Mis;
};

MagNoiseModel make_mag_noise_model(float sigma_white_uT,
                                   float bias_residual_range_uT,
                                   float sigma_bias_rw_uT_sqrt_s,
                                   float scale_err_max,
                                   float cross_axis_max,
                                   float misalign_deg_max,
                                   unsigned seed);

MagNoiseModel make_mag_noise_model(float sigma_white_uT,
                                   float bias_residual_range_uT,
                                   float sigma_bias_rw_uT_sqrt_s,
                                   float scale_err_max,
                                   float cross_axis_max,
                                   float misalign_deg_max,
                                   unsigned noise_seed,
                                   unsigned initialization_seed);

Vector3f apply_mag_noise(const Vector3f& ideal_mag_uT_body, MagNoiseModel& m, float dt_mag);

template <class MekfT>
inline Vector3f get_mag_bias_est_uT(const MekfT& mekf)
{
    if constexpr (requires { mekf.get_mag_bias_uT(); }) {
        return mekf.get_mag_bias_uT();
    } else if constexpr (requires { mekf.get_mag_bias(); }) {
        return mekf.get_mag_bias();
    } else if constexpr (requires { mekf.magnetometer_bias_uT(); }) {
        return mekf.magnetometer_bias_uT();
    } else if constexpr (requires { mekf.magnetometer_bias(); }) {
        return mekf.magnetometer_bias();
    } else {
        return Vector3f::Zero();
    }
}

struct DirectionTelemetry {
    float phase = NAN;
    float direction_deg = NAN;              // propagation axis, modulo 180 deg
    float apparent_to_deg = NAN;             // directed encounter propagation
    float apparent_from_deg = NAN;           // opposite marine "waves from" angle
    float sense_coherence = NAN;
    float direction_deg_generator_signed = NAN;
    float uncertainty_deg = NAN;
    float confidence = NAN;
    float amplitude = NAN;
    Vector2f direction_vec = Vector2f::Zero();
    // Directed propagation unit vector in the boat frame (forward, starboard).
    // Zero means the travel sense is not resolved for this sample.  Unlike the
    // FORWARD/BACKWARD class this is a physical quantity: it is invariant to
    // the 180-degree representative the axis estimator happens to return.
    Vector2f travel_vec_boat = Vector2f::Zero();
    Vector2f filtered_signal = Vector2f::Zero();
    WaveDirection sign = UNCERTAIN;
    int sign_num = 0;
};

struct FilterSnapshot {
    Vector3f disp_est_zu = Vector3f::Zero();
    Vector3f vel_est_zu = Vector3f::Zero();
    Vector3f acc_est_zu = Vector3f::Zero();
    Vector3f euler_nautical_deg = Vector3f::Zero();

    Vector3f acc_bias_est_ned = Vector3f::Zero();
    Vector3f gyro_bias_est_ned = Vector3f::Zero();
    Vector3f mag_bias_est_ned_uT = Vector3f::Zero();

    float tau_target = NAN;
    float sigma_target = NAN;
    float tuning_target = NAN;
    float tau_applied = NAN;
    float sigma_applied = NAN;
    float tuning_applied = NAN;
    float freq_hz = NAN;
    // Zero-crossing wave period from the accelerometer-only estimator, NaN when
    // the filter family does not provide one.
    float wave_period_sec = NAN;
    float period_sec = NAN;
    float accel_variance = NAN;
    float displacement_scale_m = NAN;
    float velocity_scale_mps = NAN;

    DirectionTelemetry direction;
};

/*
  Additive typed snapshot for TimeVarGain NLO.
  Existing FilterSnapshot is unchanged.
*/
struct TvgNloTelemetry {
    float k1 = NAN;
    float k2 = NAN;
    float kI = NAN;
    float vartheta = NAN;
    float theta = NAN;
    float p0z_hat = NAN;

    float wave_freq_hz = NAN;
    float wave_freq_confidence = NAN;

    Vector3f xi_n = Vector3f::Constant(NAN);
    Vector3f fhat_n = Vector3f::Constant(NAN);
    Vector3f sigma_b = Vector3f::Constant(NAN);
    Vector3f gyro_bias_b = Vector3f::Constant(NAN);

    float xi_norm = NAN;
    float fhat_norm = NAN;
    float sigma_norm = NAN;
    float gyro_bias_norm = NAN;
};

struct TvgNloFilterSnapshot {
    Vector3f disp_est_zu = Vector3f::Zero();
    Vector3f vel_est_zu = Vector3f::Zero();
    Vector3f acc_est_zu = Vector3f::Zero();
    Vector3f euler_nautical_deg = Vector3f::Zero();

    Vector3f acc_bias_est_ned = Vector3f::Zero();
    Vector3f gyro_bias_est_ned = Vector3f::Zero();
    Vector3f mag_bias_est_ned_uT = Vector3f::Zero();

    DirectionTelemetry direction;

    TvgNloTelemetry tvg;
};

class IW3dFusionAdapter {
public:
    virtual ~IW3dFusionAdapter() = default;
    virtual void updateMag(const Vector3f& mag_body_ned) = 0;
    virtual void update(float dt,
                        const Vector3f& gyr_meas_ned,
                        const Vector3f& acc_meas_ned,
                        float temperature_c) = 0;
    virtual FilterSnapshot snapshot() const = 0;
};

template <typename SnapshotT>
class IW3dFusionAdapterTyped {
public:
    using Snapshot = SnapshotT;

    virtual ~IW3dFusionAdapterTyped() = default;

    virtual void updateMag(const Vector3f& mag_body_ned) = 0;

    virtual void update(float dt,
                        const Vector3f& gyr_meas_ned,
                        const Vector3f& acc_meas_ned,
                        float temperature_c) = 0;

    virtual SnapshotT snapshot() const = 0;
};

using ImuNoiseInjector = std::function<void(Vector3f& acc_body_zu,
                                            Vector3f& gyr_body_zu,
                                            float dt)>;
using MagNoiseInjector = std::function<void(Vector3f& mag_body_enu, float dt_mag)>;

struct SimulationNoiseModels {
    std::optional<ImuNoiseModel> accel_noise;
    std::optional<ImuNoiseModel> gyro_noise;
    std::optional<MagNoiseModel> mag_noise;
    std::vector<ImuNoiseInjector> extra_imu_noise_models;
    std::vector<MagNoiseInjector> extra_mag_noise_models;
};

// ---------------------------------------------------------------------------
// Inboard-diesel engine vibration
// ---------------------------------------------------------------------------
//
// Models what a hull-mounted IMU records on a mid-size recreational cruising
// sailboat (35-45 ft) motoring under its inboard auxiliary diesel.  The
// archetype is a naturally aspirated three-cylinder four-stroke of the
// Yanmar 3YM30 / Volvo D1-30 / Beta 25 class on flexible mounts, driving a
// three-blade fixed propeller through a 2.6:1 reduction gear.
//
// The model is a sensor-path model.  It leaves the vessel's rigid-body wave
// response untouched and adds only what the engine and driveline put into the
// accelerometer and gyroscope channels:
//
//   * discrete crank orders (half order, first order, the firing order
//     n_cyl/2 and its harmonics) plus driveline shaft-rate and
//     propeller-blade-rate lines,
//   * an elevated broadband structural floor above a few Hz,
//   * governor hunting and combustion variability, which turn each line into
//     a narrow band rather than a pure tone,
//   * the flexible-mount transmissibility, which amplifies near the mount
//     resonance and rolls off above it,
//   * the sensor's finite anti-alias bandwidth, after which every remaining
//     component is folded into [0, 1/(2 dt)] by the sample-rate phase
//     accumulation -- this is the mechanism that puts engine energy into the
//     wave band,
//   * accelerometer vibration rectification and gyroscope g-sensitivity.
//
// Amplitudes are physical rather than fitted: for a line of order k at crank
// frequency f_c the excitation is inertial and grows as f^2, and the flexible
// mount transmits it with
//
//   T(f) = sqrt( (1 + (2 zeta r)^2) / ((1 - r^2)^2 + (2 zeta r)^2) ),  r = f/f_n
//
// so a single overall gain, calibrated once so the hull broadband RMS equals
// level_mps2 at reference_rpm, fixes the level at every other speed.  The
// resulting speed dependence is the familiar one: a diesel auxiliary is rough
// at idle, where the low orders sit near the mount resonance, and does not
// simply fall as rpm^2.
struct EngineVibrationConfig {
    // Engine speed in rpm.  Zero or negative disables the model entirely.
    float rpm = 0.0f;
    // Four-stroke cylinder count; the firing order is cylinders/2.
    int cylinders = 3;
    // Reduction-gear ratio, engine rev per shaft rev.
    float gear_ratio = 2.6f;
    // Propeller blade count, for the blade-rate line.
    int blades = 3;
    // Hull broadband vibration RMS over all three axes at reference_rpm,
    // before the sensor's anti-alias filter.  0.60 m/s^2 is about 0.061 g and
    // corresponds to roughly 3 mm/s RMS velocity near 30 Hz, i.e. the ISO 6954
    // comfort range for small-craft accommodation adjacent to the engine bay.
    float level_mps2 = 0.60f;
    float reference_rpm = 2400.0f;
    // Flexible engine mount: natural frequency and damping ratio.  Marine
    // mounts for this engine class are chosen around 8-12 Hz.
    float mount_hz = 10.0f;
    float mount_zeta = 0.12f;
    // Sensor bandwidth ahead of the sample rate, as a two-pole rolloff.  A
    // consumer MEMS IMU delivering 200 Hz typically leaves 50-100 Hz of
    // usable bandwidth, so the high crank orders are attenuated but not
    // removed before they fold.
    float sensor_bandwidth_hz = 80.0f;
    // Effective lever arm converting hull linear vibration into the angular
    // vibration the gyroscope sees.
    float gyro_lever_m = 1.5f;
    // Gyroscope linear-acceleration sensitivity.  1.78e-4 (rad/s)/(m/s^2) is
    // 0.1 deg/s/g, the typical figure for a BMI270-class part.
    float gyro_g_sensitivity = 1.78e-4f;
    // Accelerometer vibration rectification, in mg per g^2 of vibration.
    float vre_mg_per_g2 = 1.0f;
    // Seconds into the record at which the engine is shut down; the model
    // contributes nothing from then on.  Zero or negative means it runs for
    // the whole record, which is the default.  This exists so a study can
    // watch the estimator's vibration guard release rather than only engage.
    float stop_sec = 0.0f;
    unsigned seed = 20260828u;
};

struct EngineVibrationLine {
    // True, unaliased line frequency in Hz.
    float freq_hz = 0.0f;
    // Multiplier taking the instantaneous crank frequency to this line's
    // frequency, so governor wander moves every line coherently.
    float freq_ratio = 0.0f;
    // Recorded per-axis amplitude in m/s^2, after the anti-alias rolloff.
    Vector3f amp = Vector3f::Zero();
    // Fixed per-axis phase offset from the structural transfer path.
    Vector3f phase_offset = Vector3f::Zero();
    // Running phase, advanced at the instantaneous frequency.
    float phase = 0.0f;
    // Slow multiplicative amplitude modulation state.
    float mod = 0.0f;
};

struct EngineVibrationModel {
    EngineVibrationConfig cfg;
    std::mt19937 rng;
    std::normal_distribution<float> n01;

    float crank_hz = 0.0f;
    float shaft_hz = 0.0f;
    float firing_order = 1.5f;

    std::vector<EngineVibrationLine> lines;

    // Governor: an Ornstein-Uhlenbeck speed deviation plus a periodic hunt.
    float rpm_ou = 0.0f;
    float rpm_ou_alpha = 0.0f;
    float rpm_ou_sigma = 0.0f;
    float hunt_phase = 0.0f;
    float hunt_hz = 0.35f;
    float hunt_amplitude = 0.002f;

    // Amplitude modulation from cycle-to-cycle combustion variability.
    float mod_alpha = 0.0f;
    float mod_sigma = 0.0f;
    float mod_depth = 0.12f;

    // Broadband structural floor: white noise through a one-pole high pass,
    // so the engine contributes little directly below a few Hz.
    Vector3f broadband_sigma = Vector3f::Zero();
    Vector3f hp_prev_in = Vector3f::Zero();
    Vector3f hp_out = Vector3f::Zero();
    float hp_alpha = 0.0f;

    // Elapsed record time, for the optional shutdown.
    float time_sec = 0.0f;

    // Constant accelerometer offset from vibration rectification, and the
    // hull vibration RMS per axis it was computed from.
    Vector3f vre_offset = Vector3f::Zero();
    Vector3f hull_rms = Vector3f::Zero();
    Vector3f recorded_rms = Vector3f::Zero();
};

EngineVibrationModel make_engine_vibration_model(const EngineVibrationConfig& cfg, float dt);

void apply_engine_vibration(Vector3f& acc_body_zu,
                            Vector3f& gyr_body_zu,
                            EngineVibrationModel& model,
                            float dt);

// Reads W3D_ENGINE_RPM and the optional W3D_ENGINE_* overrides.  Returns
// nullopt when no engine is configured, which is the default and reproduces
// the historical noise realization bit for bit.
std::optional<EngineVibrationConfig> w3d_engine_vibration_from_env();

// Builds the model from the environment, prints an ENGINE_VIBRATION banner
// describing it, and appends it to the simulator's extra IMU noise models.
// A no-op when no engine is configured.
void w3d_install_engine_vibration_from_env(SimulationNoiseModels& noise_models, float dt);

// Default values preserve the historical deterministic validation realization.
// When W3D_SEED, W3D_IMU_SEED, or W3D_INIT_SEED is supplied, the simulator
// expands the corresponding base seed into independent sensor streams.
struct W3dRandomSeeds {
    unsigned accel_noise = 1234u;
    unsigned gyro_noise = 5678u;
    unsigned mag_noise = 9012u;
    unsigned accel_initialization = 1234u;
    unsigned gyro_initialization = 5678u;
    unsigned mag_initialization = 9012u;
};

unsigned w3d_expand_seed(unsigned base_seed, unsigned stream_id);
W3dRandomSeeds w3d_random_seeds_from_env();

// ---------------------------------------------------------------------------
// IMU installation lever arm and its filter-side model
// ---------------------------------------------------------------------------
//
// The versioned wave records carry specific force at the vessel centre of
// gravity (CG).  A real IMU is bolted somewhere else on the hull.  For a
// body-fixed offset r the rigid-body relation between the two locations is
//
//   a_IMU = a_CG + omega_dot x r + omega x (omega x r),
//
// a tangential term and a centripetal term.  Both are deterministic and both
// are correlated with attitude, so a few decimetres of installation offset
// inject a structured error into the one channel the OU filters use for
// attitude and wave acceleration at the same time.
//
// The two halves of the experiment sit on opposite sides of the sensor and
// are therefore modelled as separate stages:
//
//   * install()    runs on the record's noiseless truth, before sensor
//                  corruption.  It is the physics: what a sensor at r really
//                  measures.
//   * compensate() runs after sensor corruption, immediately before fusion.
//                  It is the filter's own lever-arm model, and it is the
//                  stage a firmware implementation would own.
//
// Keeping them apart is what makes the comparison meaningful.  The oracle
// model removes exactly what the installation stage added and bounds the
// recoverable error; the deployable model has to rebuild the same term from
// the noisy, biased rate the filter actually receives.
struct W3dLeverArmConfig {
    enum class Model {
        None,          // no filter-side model: the installation error stands
        Exact,         // oracle: the record's own angular kinematics
        MeasuredGyro,  // deployable: band-limited derivative of the measured rate
    };

    // Offset of the IMU from the CG, in the body z-up frame, metres.
    Vector3f offset_body_zu = Vector3f::Zero();
    Model model = Model::None;
    // Corner frequency of the two-pole low-pass the MeasuredGyro model runs
    // ahead of its derivative.  This is the model's one real design choice
    // and it is a two-sided one.  Too narrow and the low-pass phase lag
    // misaligns a term that is already the right size, so the "correction"
    // adds error; too wide and differentiating the raw 200 Hz rate hands the
    // filter more white noise than lever-arm signal.  15 Hz sits in the flat
    // basin between the two over the whole Hs range measured here.
    float derivative_cutoff_hz = 15.0f;
};

// Rigid-body acceleration of a point at r relative to the body origin.
inline Vector3f w3d_lever_acceleration(const Vector3f& omega_radps,
                                       const Vector3f& alpha_radps2,
                                       const Vector3f& r_body_m)
{
    return alpha_radps2.cross(r_body_m) +
        omega_radps.cross(omega_radps.cross(r_body_m));
}

// Causal second-order backward difference of a body rate, with a first-order
// start-up so the first two samples of a record stay finite and bounded.
class W3dRateDerivative {
public:
    explicit W3dRateDerivative(float dt) : dt_(dt) {}

    Vector3f update(const Vector3f& omega);

private:
    float dt_;
    Vector3f prev_ = Vector3f::Zero();
    Vector3f prev2_ = Vector3f::Zero();
    int seen_ = 0;
};

class W3dLeverArm {
public:
    W3dLeverArm(const W3dLeverArmConfig& cfg, float dt);

    const W3dLeverArmConfig& config() const { return cfg_; }
    bool installs() const { return offset_norm_ > 0.0f; }
    bool compensates() const { return cfg_.model != W3dLeverArmConfig::Model::None; }
    float offset_norm_m() const { return offset_norm_; }

    // Physics stage: move the record's CG specific force to the sensor
    // location using the record's own body rate.  Call once per sample,
    // before any sensor corruption.
    Vector3f install(const Vector3f& acc_cg_body_zu,
                     const Vector3f& gyr_truth_body_zu);

    // Filter stage: remove the modelled term from the measured specific
    // force.  Call exactly once per install(), after sensor corruption and
    // immediately before fusion; the residual accounting below assumes the
    // pairing.  With Model::None it removes nothing, which is what makes the
    // unmodeled arm report the whole installed term as its residual.
    Vector3f compensate(const Vector3f& acc_meas_body_zu,
                        const Vector3f& gyr_meas_body_zu);

    // RMS magnitude of the term the installation added, and of the residual
    // the filter is left with after its own model.  Both are accumulated over
    // the whole record so a study can chart the mechanism and not only the
    // score.  Equal values mean the model removed nothing; a residual near
    // zero means it removed everything.
    float installed_rms_mps2() const;
    float residual_rms_mps2() const;
    std::size_t samples() const { return samples_; }

private:
    W3dLeverArmConfig cfg_;
    Vector3f offset_ = Vector3f::Zero();
    float offset_norm_ = 0.0f;
    float lp_gain_ = 1.0f;

    W3dRateDerivative truth_derivative_;
    W3dRateDerivative model_derivative_;
    Vector3f lp_stage1_ = Vector3f::Zero();
    Vector3f lp_stage2_ = Vector3f::Zero();
    bool lp_primed_ = false;

    Vector3f last_installed_ = Vector3f::Zero();
    double installed_sumsq_ = 0.0;
    double residual_sumsq_ = 0.0;
    std::size_t samples_ = 0;
};

// Reads W3D_IMU_LEVER_ARM_M ("x,y,z" in metres, body z-up) and the optional
// W3D_IMU_LEVER_ARM_MODEL (none|exact|gyro) and W3D_IMU_LEVER_ARM_CUTOFF_HZ.
// Returns nullopt when no lever arm is configured, which is the default and
// reproduces the historical realization bit for bit.
std::optional<W3dLeverArmConfig> w3d_lever_arm_config_from_env();

// Builds the stage from the environment and prints an IMU_LEVER_ARM banner
// describing it.  Returns nullptr when no lever arm is configured.
std::shared_ptr<W3dLeverArm> w3d_lever_arm_from_env(float dt);

// Prints the post-record IMU_LEVER_ARM_RESULT diagnostic line.  A no-op for a
// null stage, so callers need no branch of their own.
void w3d_report_lever_arm(const std::shared_ptr<W3dLeverArm>& lever_arm);

struct W3dSimulationOptions {
    float dt = 0.005f;
    bool with_mag = true;
    bool add_noise = true;
    float mag_odr_hz = 25.0f;
    float temperature_c = 35.0f;
    bool write_timeseries = true;
    std::string output_suffix_with_mag = "_fusion";
    std::string output_suffix_no_mag = "_fusion_nomag";
    // Optional IMU installation lever arm.  Null is the default and leaves
    // both the truth path and the fusion input exactly as they were.
    std::shared_ptr<W3dLeverArm> lever_arm;
};

struct W3dSimulationRunResult {
    std::string input_name;
    std::string output_name;
    WaveType wave_type = WaveType::JONSWAP;
    WaveParameters wave_params{};
    bool with_mag = true;

    std::vector<float> errs_x, errs_y, errs_z, errs_roll, errs_pitch, errs_yaw;
    std::vector<float> ref_x, ref_y, ref_z;
    std::vector<float> accb_err_x, accb_err_y, accb_err_z;
    std::vector<float> gyrb_err_x, gyrb_err_y, gyrb_err_z;
    std::vector<float> magb_err_x, magb_err_y, magb_err_z;
    std::vector<float> accb_true_x, accb_true_y, accb_true_z;
    std::vector<float> gyrb_true_x, gyrb_true_y, gyrb_true_z;
    std::vector<float> magb_true_x, magb_true_y, magb_true_z;
    std::vector<float> freq_hist;
    std::vector<float> dir_deg_hist, dir_unc_hist, dir_conf_hist, dir_amp_hist, dir_phase_hist;
    // Directed propagation angle in the generator convention, with the vessel
    // heading removed, so it is directly comparable with the record azimuth.
    // NaN where the travel sense is unresolved.
    std::vector<float> dir_travel_deg_hist;
    std::vector<int> dir_sign_num_hist;

    float final_tau_target = NAN;
    float final_sigma_target = NAN;
    float final_tuning_target = NAN;
    float final_tau_applied = NAN;
    float final_sigma_applied = NAN;
    float final_tuning_applied = NAN;
    float final_freq_hz = NAN;
    float final_wave_period_sec = NAN;
    float final_period_sec = NAN;
    float final_accel_variance = NAN;
};

template <typename SnapshotT>
struct W3dSimulationRunResultTyped : public W3dSimulationRunResult {
    using Snapshot = SnapshotT;

    std::vector<SnapshotT> snapshots;
    SnapshotT final_snapshot{};
};

using TvgNloSimulationRunResult = W3dSimulationRunResultTyped<TvgNloFilterSnapshot>;

struct W3dFailureLimits {
    float err_limit_percent_z_jonswap = 0.0f;
    float err_limit_percent_z_pmstokes = 0.0f;
    float err_limit_yaw_deg = 0.0f;
    // Roll and pitch, like yaw, in degrees RMS over the scored window.  Zero
    // means the channel is not gated for this family, which is what the
    // families that have never fitted a bar for it get: a sentinel of zero
    // would fail every run, and silence is the correct behaviour for a bar
    // nobody has measured.
    float err_limit_roll_deg = 0.0f;
    float err_limit_pitch_deg = 0.0f;
    float err_limit_percent_3d_jonswap = 0.0f;
    float err_limit_percent_3d_pmstokes = 0.0f;
    float acc_z_bias_percent = 0.0f;
    // Fitted to the accelerometer, which is the larger of the two 3D bias
    // errors by a factor of four on every family that has measured both.  A
    // single bar over both channels is therefore the accelerometer's bar, and
    // the gyro rides four times below it, catching nothing.
    float bias_3d_percent = 0.0f;
    // The gyro's own 3D bias bar.  Zero keeps the shared behaviour above, so a
    // family that has not fitted one is gated exactly as before.
    float gyro_bias_3d_percent = 0.0f;
};

struct W3dSummaryLabels {
    const char* target = "RS_target";
    const char* applied = "RS_applied";
};

class W3dSimulationRunner {
public:
    W3dSimulationRunner(W3dSimulationOptions options,
                        SimulationNoiseModels noise_models,
                        IW3dFusionAdapter& fusion_adapter);

    std::optional<W3dSimulationRunResult> run(const std::string& filename);

private:
    std::string make_output_name(const std::string& filename) const;

    W3dSimulationOptions options_;
    SimulationNoiseModels noise_models_;
    IW3dFusionAdapter& fusion_adapter_;
};

class TvgNloSimulationRunner {
public:
    using Adapter = IW3dFusionAdapterTyped<TvgNloFilterSnapshot>;

    TvgNloSimulationRunner(W3dSimulationOptions options,
                           SimulationNoiseModels noise_models,
                           Adapter& fusion_adapter);

    std::optional<TvgNloSimulationRunResult> run(const std::string& filename);

private:
    std::string make_output_name(const std::string& filename) const;

    W3dSimulationOptions options_;
    SimulationNoiseModels noise_models_;
    Adapter& fusion_adapter_;
};

void print_summary_and_fail_if_needed(const W3dSimulationRunResult& result,
                                      float dt,
                                      const W3dFailureLimits& limits,
                                      const W3dSummaryLabels& labels = {});

void print_validation_metrics(const W3dSimulationRunResult& result,
                              float dt,
                              float window_seconds,
                              const char* family);

std::vector<std::string> collect_wave_data_files(const std::filesystem::path& directory);
bool w3d_any_quality_gate_failed();

template <typename AdapterT>
inline std::optional<W3dSimulationRunResult> process_wave_file_for_tracker(const std::string& filename,
                                                                           float dt,
                                                                           bool with_mag,
                                                                           bool add_noise,
                                                                           float mag_odr_hz,
                                                                           std::string output_suffix_with_mag = "_fusion",
                                                                           std::string output_suffix_no_mag = "_fusion_nomag",
                                                                           W3dRandomSeeds seeds = {},
                                                                           bool write_timeseries = true)
{
    const float acc_sigma = 1.51e-3f * g_std;
    const float gyr_sigma = 0.00157f;
    const float acc_bias_range = 8e-3f * g_std;
    const float gyr_bias_range = 0.05f * float(std::numbers::pi_v<float> / 180.0f);
    const float acc_bias_rw = 0.0005f;
    const float gyr_bias_rw = 0.00001f;
    const float mag_sigma_uT = (mag_odr_hz <= 20.0f) ? 0.40f : 0.80f;

    SimulationNoiseModels noise_models;
    noise_models.accel_noise = make_imu_noise_model(
        acc_sigma, acc_bias_range, acc_bias_rw,
        seeds.accel_noise, seeds.accel_initialization);
    noise_models.gyro_noise = make_imu_noise_model(
        gyr_sigma, gyr_bias_range, gyr_bias_rw,
        seeds.gyro_noise, seeds.gyro_initialization);
    noise_models.mag_noise = make_mag_noise_model(
        mag_sigma_uT, 2.0f, 0.01f, 0.015f, 0.010f, 1.0f,
        seeds.mag_noise, seeds.mag_initialization);

    // Optional inboard-diesel vibration, on top of the sensor noise models
    // and identical for every filter family.  Absent unless W3D_ENGINE_RPM is
    // set, so the historical realization is untouched by default.
    w3d_install_engine_vibration_from_env(noise_models, dt);

    const Vector3f sigma_a_init(2.8f * acc_sigma, 2.8f * acc_sigma, 2.8f * acc_sigma);
    const Vector3f sigma_g(2.0f * gyr_sigma, 2.0f * gyr_sigma, 2.0f * gyr_sigma);
    const float sigma_m_uT = 1.2f * mag_sigma_uT;
    const Vector3f sigma_m(sigma_m_uT, sigma_m_uT, sigma_m_uT);

    AdapterT adapter(with_mag, sigma_a_init, sigma_g, sigma_m);
    W3dSimulationOptions options;
    options.dt = dt;
    options.with_mag = with_mag;
    options.add_noise = add_noise;
    options.mag_odr_hz = mag_odr_hz;
    options.temperature_c = 35.0f;
    options.write_timeseries = write_timeseries;
    options.output_suffix_with_mag = std::move(output_suffix_with_mag);
    options.output_suffix_no_mag = std::move(output_suffix_no_mag);

    // Optional IMU installation lever arm and its filter-side model.  Absent
    // unless W3D_IMU_LEVER_ARM_M is set, so the historical realization is
    // untouched by default.  The runner copies the options but shares this
    // stage, so the diagnostics it accumulates are readable here once the
    // record is done.
    options.lever_arm = w3d_lever_arm_from_env(dt);

    W3dSimulationRunner runner(options, std::move(noise_models), adapter);
    auto result = runner.run(filename);
    w3d_report_lever_arm(options.lever_arm);
    return result;
}

#endif
