// ImuAccuracy (src/imu_calibrate/ImuAccuracy.h) on simulated still poses and
// rotations with known sensor errors: each error must show up, at its size,
// in the quantity meant to reveal it.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "imu_calibrate/ImuAccuracy.h"

using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::Vector3f;

static int failures = 0;
static void require(bool ok, const char* what) {
  std::printf("  %-72s %s\n", what, ok ? "ok" : "FAIL");
  if (!ok) ++failures;
}

constexpr double PI = 3.14159265358979323846, D2R = PI / 180.0, G = 9.8025605, G_STD = 9.80665;
constexpr double DT = 0.005;

struct Sensor {
  double accel_scale = 1.0, gyro_scale = 1.0;
  Vector3d gyro_bias = Vector3d::Zero();
  Matrix3d mag_map = Matrix3d::Identity();
  double accel_noise = 0.0015 * G_STD, gyro_noise = 0.07 * D2R, mag_noise = 0.15;
  double hand_tremor = 0;  // m/s^2, added to every accelerometer sample
};

// Field in NED (uT) for F, inclination, declination.
static Vector3d fieldNed(double f, double inc_deg, double dec_deg) {
  const double i = inc_deg * D2R, d = dec_deg * D2R;
  return f * Vector3d(std::cos(i) * std::cos(d), std::cos(i) * std::sin(d), std::sin(i));
}

struct Segment { Vector3d axis; double angle_deg; double seconds; };

struct Run {
  std::vector<imu_cal::StillPose> poses;
  std::vector<imu_cal::RotationCheck> checks;
  imu_cal::AccuracySummary summary;
};

// Still for `still_s`, then each motion followed by `still_s` still.
static Run simulate(const Sensor& s, const std::vector<Segment>& motions, double still_s, const Vector3d& m_ned,
                    unsigned seed) {
  std::mt19937 rng(seed);
  std::normal_distribution<double> n01(0.0, 1.0);
  imu_cal::ImuAccuracy acc;
  imu_cal::AccuracyReference ref;
  ref.gravity_mps2 = G;
  ref.have_field = true;
  ref.f_uT = m_ned.norm();
  ref.h_uT = std::hypot(m_ned.x(), m_ned.y());
  ref.z_uT = m_ned.z();
  ref.inclination_deg = std::atan2(ref.z_uT, ref.h_uT) / D2R;
  acc.setReference(ref);
  Run run;
  Quaterniond q = Quaterniond::Identity();  // body to NED
  long k = 0;
  auto step = [&](const Vector3d& w_body) {
    q = (q * Quaterniond(Eigen::AngleAxisd(w_body.norm() * DT, w_body.norm() > 0 ? w_body.normalized()
                                                                                    : Vector3d::UnitX())))
            .normalized();
    const Matrix3d Rt = q.toRotationMatrix().transpose();
    const Vector3d f = Rt * Vector3d(0, 0, -G);
    Vector3d a = s.accel_scale * f, w = s.gyro_scale * w_body + s.gyro_bias;
    for (int i = 0; i < 3; ++i) {
      a[i] += s.accel_noise * n01(rng) + s.hand_tremor * n01(rng);
      w[i] += s.gyro_noise * n01(rng);
    }
    acc.addFrame(a.cast<float>(), w.cast<float>(), a.cast<float>(), w.cast<float>(), DT, 30.0f);
    if (k++ % 8 == 0) {
      Vector3d m = s.mag_map * (Rt * m_ned);
      for (int i = 0; i < 3; ++i) m[i] += s.mag_noise * n01(rng);
      acc.addMag(m.cast<float>(), m.cast<float>());
    }
    imu_cal::StillPose p;
    if (acc.takePose(p)) run.poses.push_back(p);
    imu_cal::RotationCheck c;
    if (acc.takeRotationCheck(c)) run.checks.push_back(c);
  };
  auto still = [&](double sec) { for (int i = 0; i < int(sec / DT); ++i) step(Vector3d::Zero()); };
  still(still_s);
  for (const Segment& m : motions) {
    const int n = int(m.seconds / DT);
    const Vector3d w = m.axis.normalized() * (m.angle_deg * D2R / m.seconds);
    for (int i = 0; i < n; ++i) step(w);
    still(still_s);
  }
  // A pose ends at the first block with motion: finish the last one.
  for (int i = 0; i < 60; ++i) step(Vector3d(1.0, 0, 0));
  run.summary = acc.summary();
  return run;
}

static const std::vector<Segment> TOUR = {
    {Vector3d::UnitX(), 60, 2.0},   {Vector3d::UnitY(), -45, 1.5}, {Vector3d::UnitZ(), 90, 3.0},
    {Vector3d::UnitX(), -100, 3.0}, {Vector3d::UnitY(), 80, 2.5},  {Vector3d(1, 1, 0), 120, 3.0},
    {Vector3d::UnitX(), 180, 4.0},
};

int main() {
  const Vector3d m_ned = fieldNed(50.886, 65.72, -12.45);

  std::printf("calibrated sensor\n");
  {
    Sensor s;
    const Run r = simulate(s, TOUR, 6.0, m_ned, 1);
    std::printf("  %zu poses, %u distinct, %zu rotation checks, |m| std %.3f %% over %d faces\n", r.poses.size(),
                r.summary.distinct_poses, r.checks.size(), r.summary.mag_norm_std_pct, r.summary.mag_faces);
    require(r.poses.size() == TOUR.size() + 1, "every still pose found, none in motion");
    bool long_enough = true;
    for (const auto& p : r.poses) long_enough = long_enough && p.duration_s > 4.5 && p.duration_s < 6.0;
    require(long_enough, "each pose keeps 4.5-6 s of its 6 s (settling and motion edges dropped)");
    require(std::fabs(r.summary.accel_error_mg_min) < 0.3 && std::fabs(r.summary.accel_error_mg_max) < 0.3,
            "|a| matches gravity within 0.3 mg in every pose");
    require(r.summary.accel_noise_mg_max > 1.2 && r.summary.accel_noise_mg_max < 1.9, "accelerometer noise ~1.5 mg");
    require(r.summary.gyro_bias_dps_max < 0.01, "no gyro bias reported");
    require(r.summary.gyro_noise_dps_max > 0.055 && r.summary.gyro_noise_dps_max < 0.09, "gyro noise ~0.07 deg/s");
    require(std::fabs(r.summary.f_error_pct_min) < 0.05 && std::fabs(r.summary.f_error_pct_max) < 0.05,
            "F matches the expected field within 0.05 %");
    require(r.summary.inclination_error_deg_max < 0.1, "dip matches the expected 65.72 deg within 0.1 deg");
    require(r.summary.z_error_uT_max < 0.1, "vertical component matches within 0.1 uT");
    require(r.checks.size() == TOUR.size() && r.summary.rotation_accel_error_deg_max < 0.15,
            "gyro rotation predicts the gravity direction within 0.15 deg");
    require(r.summary.rotation_mag_error_deg_max < 0.3, "and the field direction within 0.3 deg");
    require(r.summary.scale_fits >= 5 && std::fabs(r.summary.gyro_scale_pct_min) < 0.15 &&
                std::fabs(r.summary.gyro_scale_pct_max) < 0.15,
            "fitted gyro scale error within 0.15 %, the 180 deg flip included");
    require(r.summary.mag_norm_std_pct < 0.5 && r.summary.mag_faces >= 4, "|m| constant in motion over 4+ faces");
    require(r.summary.distinct_poses >= 5, "distinct gravity directions counted");
  }

  std::printf("sensor with known errors\n");
  {
    Sensor s;
    s.accel_scale = 1.004;                          // +4 mg at 1 g
    s.gyro_scale = 1.02;                            // 2 % scale error
    s.gyro_bias = Vector3d(0.3, -0.2, 0.1) * D2R;   // |b| = 0.374 deg/s
    s.mag_map = Matrix3d::Identity() * 1.03;        // field reads 3 % high
    const Run r = simulate(s, TOUR, 6.0, m_ned, 2);
    std::printf("  |a| error %.2f..%.2f mg, bias %.3f deg/s, F error %.2f..%.2f %%, gyro scale %.2f..%.2f %% (%u fits)\n",
                r.summary.accel_error_mg_min, r.summary.accel_error_mg_max, r.summary.gyro_bias_dps_max,
                r.summary.f_error_pct_min, r.summary.f_error_pct_max, r.summary.gyro_scale_pct_min,
                r.summary.gyro_scale_pct_max, r.summary.scale_fits);
    const double expect_mg = 1000.0 * (1.004 * G - G) / G_STD;
    require(std::fabs(r.summary.accel_error_mg_min - expect_mg) < 0.3 &&
                std::fabs(r.summary.accel_error_mg_max - expect_mg) < 0.3,
            "accelerometer scale error shows as |a| error (+4 mg)");
    require(std::fabs(r.summary.gyro_bias_dps_max - 0.374) < 0.01, "gyro bias shows as the still mean rate");
    require(std::fabs(r.summary.f_error_pct_min - 3.0) < 0.06 && std::fabs(r.summary.f_error_pct_max - 3.0) < 0.06,
            "magnetometer scale error shows as F error (+3 %)");
    require(r.summary.inclination_error_deg_max < 0.1, "an isotropic field scale leaves the dip unchanged");
    require(r.summary.scale_fits >= 5 && std::fabs(r.summary.gyro_scale_pct_min - 2.0) < 0.15 &&
                std::fabs(r.summary.gyro_scale_pct_max - 2.0) < 0.15,
            "gyro scale error fitted as +2 % on every rotation (bias removed between poses)");
    double worst_mag = 0;
    for (const auto& c : r.checks)
      if (std::isfinite(c.gyro_scale_mag_pct)) worst_mag = std::max(worst_mag, std::fabs(c.gyro_scale_mag_pct - 2.0));
    require(r.summary.mag_scale_fits == TOUR.size() && worst_mag < 0.3,
            "and from the field direction on every rotation, the turn about the vertical too");
  }

  std::printf("soft-iron residual\n");
  {
    Sensor s;
    s.mag_map = Eigen::DiagonalMatrix<double, 3>(1.06, 0.95, 1.0);
    const Run r = simulate(s, TOUR, 6.0, m_ned, 3);
    std::printf("  |m| std %.2f %%, dip %.2f..%.2f deg\n", r.summary.mag_norm_std_pct, r.summary.inclination_deg_min,
                r.summary.inclination_deg_max);
    require(r.summary.mag_norm_std_pct > 1.5, "|m| varies with orientation");
    require(r.summary.inclination_deg_max - r.summary.inclination_deg_min > 1.0, "dip varies between poses");
  }

  std::printf("noisy gyro\n");
  {
    Sensor s;
    s.gyro_noise = 0.3 * D2R;  // above the worst 0.25 s block of an AtomS3R at rest (0.26 deg/s, GYR_CONF 0xA9)
    const Run r = simulate(s, TOUR, 6.0, m_ned, 5);
    std::printf("  %zu poses, gyro noise %.3f deg/s\n", r.poses.size(), r.summary.gyro_noise_dps_max);
    require(r.poses.size() == TOUR.size() + 1, "every still pose found");
    require(r.summary.gyro_noise_dps_max > 0.27 && r.summary.gyro_noise_dps_max < 0.33,
            "and the gyro noise reported at its size (above the 0.2 deg/s limit)");
  }

  std::printf("hand-held\n");
  {
    Sensor s;
    s.hand_tremor = 0.2;  // ~20 mg rms tremor
    const Run r = simulate(s, TOUR, 6.0, m_ned, 4);
    require(r.poses.empty(), "a trembling board yields no still pose");
  }

  std::printf("imu_accuracy: %s\n", failures ? "FAIL" : "PASS");
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
