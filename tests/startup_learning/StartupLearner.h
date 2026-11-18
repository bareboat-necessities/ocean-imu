#pragma once
// Host-only startup experiment. No shipping filter includes this file.
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <string>
#include <vector>

namespace startup_learning {
using V = Eigen::Vector3d;
using M = Eigen::Matrix3d;
using Q = Eigen::Quaterniond;
constexpr double gravity = 9.80665;

// All vectors use body NED axes. Rates are instantaneous endpoint samples.
// A magnetic observation belongs to this exact acquisition epoch, not its
// later delivery time. Calibrated sensor units: m/s^2, rad/s, microtesla.
struct Sample {
    double time = 0;
    V accel = V::Zero(), gyro = V::Zero(), mag = V::Zero();
    bool magnetic = false;
};

struct Handover {
    double time = 0;
    Q attitude = Q::Identity(); // body -> magnetic NED, at time
    V field = V::Zero();       // magnetic NED: [horizontal, 0, down], uT
    V down_body = V::UnitZ();  // unit gravity direction, fixed magnitude g0
    V gyro_bias = V::Zero(), accel_bias = V::Zero();
    V mag_bias = V::Zero(); // residual hard iron; joint_cal ablation only
    V angular_rate = V::Zero(); // instantaneous body rate, gyro - gyro_bias
    V last_raw_gyro = V::Zero(); // prime next interval; do not replay this one
    M gyro_bias_cov = M::Identity();
    double magnetic_rms = 0, gravity_norm_error = 0, horizontal_field = 0;
    double gravity_sensitivity_rad = 0;
    bool valid = false, ready = false, accel_bias_learned = false;
};

struct Config {
    V gyro_bias_prior = V::Zero(), accel_bias_prior = V::Zero();
    double gyro_bias_sigma = 0.001; // residual calibration uncertainty, rad/s
    double accel_bias_sigma = 0.08; // retained prior, NOT estimated in waves
    double mag_sigma = 0.8;
    double mag_bias_sigma = 1.2; // microtesla, for joint_cal only
    // Model uncertainty of a tapered force slope over the full window. This
    // admits residual wave/turn acceleration; it is not per-sample white noise.
    double force_slope_sigma = 0.005; // m/s^3 at a 35 s horizon
    double min_seconds = 20, max_seconds = 35;
    bool scale_force_slope = false; // experimental joint fit, fixed at the start of each solve
};

class Learner {
public:
    explicit Learner(Config config = {});
    void add(const Sample& sample);
    // uniform/tapered: gyro propagation using calibration prior only.
    // magfit: gyro/field fit followed by tapered gravity-direction estimate.
    // joint: also constrain the long-window transported force slope, to
    // inform gyro drift about the otherwise weak magnetic-field axis.
    Handover estimate(const std::string& method) const;
    const std::vector<Sample>& samples() const { return samples_; }
private:
    Config config_;
    std::vector<Sample> samples_;
};

Q rotation(const V& angle);
V rotation_log(const Q& q);
} // namespace startup_learning
