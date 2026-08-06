/*
  Trajectory transition matrices for the OU-III performance error.

  The pointwise certificate (generate_jacobians.cpp) freezes the operating
  point, settles the covariance there, and linearizes. That measures a system
  that never moves, and it badly understates this estimator: the frozen mode is
  hours, while the transition matrix over one wave period contracts by a factor
  of two to five.

  The switched certificate improved on that by multiplying frozen-point
  Jacobians along an attitude orbit, but each factor still carried a covariance
  settled at its own operating point and held for a dwell interval. This
  program removes that approximation. It flies the filter along a continuous
  wave trajectory, takes the one-step Jacobian at every step against the
  filter's actual running covariance, and accumulates the product over
  windows of one wave period.

  Output is one transition matrix Phi(k+N, k) per window, which is what a
  linear-time-varying stability argument needs. verify_trajectory_certificate.py
  consumes it.

  Build:
    g++ -O2 -std=c++20 -DEIGEN_NON_ARDUINO -I src -isystem /usr/include/eigen3 \
        tools/iss_certificate/generate_trajectory_jacobians.cpp -o gen_traj
*/
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#define private public
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#undef private

namespace {
using Filter = Kalman3D_Wave_OU_III<double, true, true>;
using Vec3 = Eigen::Vector3d;

constexpr int kNx = 21;
constexpr int kNe = 18;   // certified performance block
constexpr int kBg = 3;
constexpr int kV = 6;
constexpr int kP = 9;
constexpr int kS = 12;
constexpr int kAw = 15;
constexpr int kBa = 18;

struct Tuning {
    double dt = 1.0 / 240.0;
    double tau = 4.0;
    double sigma_aw = 2.8;
    double sigma_s = 6.0;

    /*
      Exogenous schedule for the adaptive parameters.

      The tuner in SeaStateFusionFilter_OU_III is driven by
      a_body_z_up_proxy = -(acc.z() + g) and by a frequency tracker running on
      the same signal, i.e. by the raw accelerometer only. Nothing in the
      estimator state feeds it, apart from a tilt gate that re-locks attitude
      above 70 degrees, which is outside the certified envelope (worst
      combined tilt there is 47 degrees, and 54 degrees on this trajectory).

      So on the envelope theta is a causal functional of the measurements and
      not of the state: there is no estimator-to-tuner feedback, and any
      bounded schedule in Theta is admissible. Rather than replay one tuner
      trace, sweep the whole box, out of phase per axis and faster than the
      real tuner slews, which is the conservative choice.
    */
    double sweep_period = 0.0;   // 0 disables the sweep

    void apply_sweep(double t) {
        if (!(sweep_period > 0.0)) return;
        const double w = 2.0 * M_PI / sweep_period;
        auto lerp_log = [](double lo, double hi, double u) {
            return std::exp(std::log(lo) + (std::log(hi) - std::log(lo)) * u);
        };
        tau      = lerp_log(0.55, 4.0,  0.5 * (1.0 + std::sin(w * t)));
        sigma_aw = lerp_log(0.12, 2.8,  0.5 * (1.0 + std::sin(w * t + 2.1)));
        sigma_s  = lerp_log(0.35, 6.0,  0.5 * (1.0 + std::sin(w * t + 4.2)));
    }
};

/*
  A wave trajectory, as a small sum of components so the motion is
  quasi-periodic rather than a single tone. Amplitudes are chosen to sweep the
  same attitude envelope the pointwise certificate samples.
*/
struct WaveTrajectory {
    double period = 10.0;

    double roll(double t) const {
        return 0.70 * std::sin(w(1.0) * t)
             + 0.12 * std::sin(w(1.7) * t + 0.9);
    }
    double pitch(double t) const {
        return 0.45 * std::sin(w(1.0) * t + 1.3)
             + 0.08 * std::sin(w(2.3) * t + 0.2);
    }
    double yaw(double t) const {
        return 0.20 * std::sin(w(0.5) * t + 0.4);
    }
    // Vertical wave acceleration the filter is meant to track.
    double heave_acc(double t) const {
        const double a = w(1.0);
        return -1.8 * a * a * std::sin(a * t);
    }

private:
    double w(double harmonic) const { return harmonic * 2.0 * M_PI / period; }
};

Eigen::Quaterniond body_to_world(double roll, double pitch, double yaw) {
    return Eigen::AngleAxisd(yaw, Vec3::UnitZ())
         * Eigen::AngleAxisd(pitch, Vec3::UnitY())
         * Eigen::AngleAxisd(roll, Vec3::UnitX());
}

// Body rate from the analytic attitude path, by central difference.
Vec3 body_rate(const WaveTrajectory& w, double t, double h) {
    const Eigen::Quaterniond q0 = body_to_world(w.roll(t - h), w.pitch(t - h), w.yaw(t - h));
    const Eigen::Quaterniond q1 = body_to_world(w.roll(t + h), w.pitch(t + h), w.yaw(t + h));
    Eigen::Quaterniond dq = q0.conjugate() * q1;
    if (dq.w() < 0.0) dq.coeffs() *= -1.0;
    return dq.vec() * (1.0 / h);
}

Filter make_filter(const Tuning& g, const WaveTrajectory& w) {
    const Vec3 sigma_acc = Vec3::Constant(0.0148);
    const Vec3 sigma_gyro = Vec3::Constant(0.00157);
    const Vec3 sigma_mag = Vec3::Constant(0.25);
    Filter f(sigma_acc, sigma_gyro, sigma_mag);
    f.set_aw_time_constant(g.tau);
    f.set_aw_stationary_std(Vec3::Constant(g.sigma_aw));
    f.set_RS_noise(Vec3::Constant(g.sigma_s));
    f.set_Racc_std(sigma_acc);
    f.set_Rmag(sigma_mag);
    f.set_initial_linear_uncertainty(0.25, 0.25, 0.5);
    f.set_initial_acc_bias_std(0.02);
    f.set_Q_bacc_rw(Vec3::Constant(5.0e-4));
    f.set_mag_world_ref(Filter::ned_field_from_decl_incl(0.0, 1.05, 48.0));
    const Eigen::Quaterniond qbw = body_to_world(w.roll(0.0), w.pitch(0.0), w.yaw(0.0));
    f.qref = qbw.conjugate();
    f.qref.normalize();
    f.xext.setZero();
    f.Pext = 0.5 * (f.Pext + f.Pext.transpose());
    return f;
}

Vec3 predicted_acc(const Filter& f) {
    const Vec3 gravity_world(0.0, 0.0, 9.80665);
    return f.qref * (f.xext.segment<3>(kAw) - gravity_world) + f.xext.segment<3>(kBa);
}

Vec3 predicted_mag(const Filter& f) {
    return f.qref * f.v2ref;
}

void retune(Filter& f, const Tuning& g) {
    f.set_aw_time_constant(g.tau);
    f.set_aw_stationary_std(Vec3::Constant(g.sigma_aw));
    f.set_RS_noise(Vec3::Constant(g.sigma_s));
}

void apply_cycle(Filter& f, const Vec3& omega, double dt,
                 const Vec3& acc, const Vec3& mag) {
    f.time_update(omega, dt);
    f.measurement_update_acc_only(acc, 35.0);
    f.measurement_update_mag_only(mag);
    f.applyIntegralZeroPseudoMeas();
}

void inject(Filter& f, int j, double h) {
    if (j < 3) {
        Vec3 d = Vec3::Zero();
        d(j) = h;
        Eigen::Quaterniond dq(1.0, 0.5 * d.x(), 0.5 * d.y(), 0.5 * d.z());
        dq.normalize();
        f.qref = dq * f.qref;
        f.qref.normalize();
    } else {
        f.xext(j) += h;
    }
}

Eigen::Matrix<double, kNx, 1> error_between(const Filter& pert, const Filter& nom) {
    Eigen::Matrix<double, kNx, 1> e;
    e.setZero();
    Eigen::Quaterniond dq = pert.qref * nom.qref.conjugate();
    if (dq.w() < 0.0) dq.coeffs() *= -1.0;
    e.segment<3>(0) = 2.0 * dq.vec();
    e.segment<3>(kBg) = pert.xext.segment<3>(kBg) - nom.xext.segment<3>(kBg);
    e.segment<3>(kV) = pert.xext.segment<3>(kV) - nom.xext.segment<3>(kV);
    e.segment<3>(kP) = pert.xext.segment<3>(kP) - nom.xext.segment<3>(kP);
    e.segment<3>(kS) = pert.xext.segment<3>(kS) - nom.xext.segment<3>(kS);
    e.segment<3>(kAw) = pert.xext.segment<3>(kAw) - nom.xext.segment<3>(kAw);
    e.segment<3>(kBa) = pert.xext.segment<3>(kBa) - nom.xext.segment<3>(kBa);
    return e;
}

/*
  One-step Jacobian at the filter's current state and covariance.

  The perturbed runs are copies of the same filter, so they share its running
  covariance; only the state is displaced. That is the whole point of doing
  this along a trajectory rather than at a settled operating point.
*/
Eigen::Matrix<double, kNx, kNx> step_jacobian(const Filter& base, const Vec3& omega,
                                              double dt, const Vec3& acc, const Vec3& mag) {
    Filter nominal = base;
    apply_cycle(nominal, omega, dt, acc, mag);

    Eigen::Matrix<double, kNx, kNx> F;
    const double h_att = 2.0e-6;
    const double h_add = 1.0e-6;
    for (int j = 0; j < kNx; ++j) {
        const double h = j < 3 ? h_att : h_add;
        Filter plus = base;
        Filter minus = base;
        inject(plus, j, h);
        inject(minus, j, -h);
        apply_cycle(plus, omega, dt, acc, mag);
        apply_cycle(minus, omega, dt, acc, mag);
        F.col(j) = (error_between(plus, nominal) - error_between(minus, nominal)) / (2.0 * h);
    }
    return F;
}
}  // namespace

int main(int argc, char** argv) {
    /*
      usage: gen_traj OUT [period_s] [instrumented_windows]
                          [warmup_s] [tau] [sigma_aw] [sigma_s]

      The covariance takes a few hundred seconds to settle onto the trajectory,
      and rho(Phi) drifts the whole time it is settling. Running the warmup
      without building Jacobians costs one filter cycle per step instead of
      43, so the settled window is reached about seven times faster and a
      sweep over tuning cells becomes affordable.
    */
    const std::string path = argc > 1 ? argv[1] : "iss_trajectory.txt";
    const double period = argc > 2 ? std::atof(argv[2]) : 10.0;
    const int windows = argc > 3 ? std::atoi(argv[3]) : 5;
    const double warmup_s = argc > 4 ? std::atof(argv[4]) : 340.0;

    Tuning tuning;
    if (argc > 5) tuning.tau = std::atof(argv[5]);
    if (argc > 6) tuning.sigma_aw = std::atof(argv[6]);
    if (argc > 7) tuning.sigma_s = std::atof(argv[7]);
    if (argc > 8) tuning.sweep_period = std::atof(argv[8]);

    WaveTrajectory wave;
    wave.period = period;

    const double dt = tuning.dt;
    const int steps_per_window = static_cast<int>(std::lround(period / dt));
    const int warmup = static_cast<int>(std::lround(warmup_s / dt));

    Filter f = make_filter(tuning, wave);

    double t = 0.0;
    auto advance_nominal = [&](Filter& g, double time) {
        tuning.apply_sweep(time);
        retune(g, tuning);
        const Vec3 omega = body_rate(wave, time, 1e-4);
        Filter truth = g;
        truth.time_update(omega, dt);
        const Vec3 acc = predicted_acc(truth);
        const Vec3 mag = predicted_mag(truth);
        apply_cycle(g, omega, dt, acc, mag);
    };

    for (int i = 0; i < warmup; ++i, t += dt) advance_nominal(f, t);

    std::ofstream os(path);
    if (!os) {
        std::cerr << "cannot open " << path << "\n";
        return 2;
    }
    os << std::setprecision(17);
    os << "OCEAN_IMU_ISS_TRAJECTORY_V1 " << kNe << " " << windows << " "
       << period << " " << dt << " " << steps_per_window << "\n";
    std::cout << "tau=" << tuning.tau << " sigma_aw=" << tuning.sigma_aw
              << " sigma_s=" << tuning.sigma_s
              << "  warmup=" << warmup_s << "s  period=" << period << "s"
              << "  theta_sweep=" << tuning.sweep_period << "s\n";

    for (int w = 0; w < windows; ++w) {
        Eigen::Matrix<double, kNe, kNe> Phi = Eigen::Matrix<double, kNe, kNe>::Identity();
        double peak = 1.0;

        for (int i = 0; i < steps_per_window; ++i, t += dt) {
            tuning.apply_sweep(t);
            retune(f, tuning);
            const Vec3 omega = body_rate(wave, t, 1e-4);
            Filter truth = f;
            truth.time_update(omega, dt);
            const Vec3 acc = predicted_acc(truth);
            const Vec3 mag = predicted_mag(truth);

            const auto F = step_jacobian(f, omega, dt, acc, mag);
            Phi = F.topLeftCorner<kNe, kNe>() * Phi;
            peak = std::max(peak, Phi.operatorNorm());

            apply_cycle(f, omega, dt, acc, mag);
        }

        os << "WINDOW " << w << " " << t << " " << peak << "\n";
        for (int r = 0; r < kNe; ++r) {
            for (int c = 0; c < kNe; ++c) os << Phi(r, c) << (c + 1 == kNe ? '\n' : ' ');
        }
        std::cout << "window " << w << "  t=" << t
                  << "  rho=" << std::abs(Phi.eigenvalues().array().abs().maxCoeff())
                  << "  peak||Phi||=" << peak << "\n";
    }

    std::cout << "wrote " << windows << " trajectory transition matrices to " << path << "\n";
    return 0;
}
