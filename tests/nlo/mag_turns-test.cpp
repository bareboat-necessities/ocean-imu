// Device configuration of the NLO (no GNSS, 3D magnetometer, as the
// atomS3R_ins_nlo sketch runs it) on a bench history the wave simulator does
// not exercise: a gyro residual bias of 0.3 deg/s per axis, an accelerometer
// bias, and the body yawing. Checks both startups the adapter sees:
//
//   still boot : the device lies still while the bootstrap runs, so the
//                stillness-gated gyro- and accelerometer-bias seeds apply;
//                then it is turned 90 deg every 40 s, and afterwards spun
//                continuously at 5 deg/s. Heave is also checked from the
//                moment the observer starts reporting, where an unseeded
//                vertical accelerometer bias swings it by metres.
//   moving boot: it is already yawing at 3 deg/s at boot, so no seed is taken
//                and the bias has to be learned by the observer while it
//                turns. That takes minutes, so only boundedness is checked.
//   handled boot: switched on while rocked by hand or by waves, with linear
//                acceleration, so the first accelerometer sample is not
//                gravity and the bootstrap's gyro-rate gate keeps its tilt
//                correction off. Heading must hold from the handoff on.
//   still boot, AtomS3R gyro: the still boot with the gyro noise measured on
//                the device at rest, 0.2 deg/s white plus a 10 Hz wobble of
//                0.3 deg/s. Its per-sample spread is above the still seed's
//                gate, which used to reject the seed and leave the residual
//                to be integrated into tilt.
//   slow roll boot: a 0.5 deg roll at an 8 s period, which over the
//                bootstrap is close to a constant rate. It must not be seeded
//                as a gyro bias.
//
// Before the horizontal theta floor and the bias seeds, the still-boot history
// reached 11 deg of tilt, 24 deg of heading error and 23 m of heave, with a
// startup heave swing of several metres, and the moving-boot history 30 deg,
// 180 deg and 200 m, within ten minutes. While the still seed's spread was
// taken over single samples, the AtomS3R gyro noise rejected it: the still
// boot then swung heave by 8.8 m, and the slow roll was seeded as a bias
// 0.24 deg/s wrong.
#include "nlo/TimeVarGainNLO_Adapter.h"
#include <cmath>
#include <cstdio>
#include <random>
extern const float g_std = 9.80665f;  // generic default of the shared headers

using Fusion = TimeVarGainNloAdapter<false, NloMagType::Magnetometer, float, TrackerType::PLL>;
using Vec3 = Eigen::Vector3f;
using Mat3 = Eigen::Matrix3f;

static float wrap180(float d) { return std::remainder(d, 360.0f); }

struct Limits {
  float tilt_deg, hdg_deg, heave_m, startup_heave_m;
};

static bool run(const char* name, bool still_boot, Limits lim, bool atoms3r_gyro = false) {
  Fusion f;
  auto& cfg = f.config();
  cfg.filter.expected_mag_norm = 0.0f;
  f.reset();

  std::mt19937 rng(7);
  std::normal_distribution<float> n01(0.0f, 1.0f);
  auto noise = [&](float s) { return Vec3(s * n01(rng), s * n01(rng), s * n01(rng)); };

  const float dt = 0.005f;
  const float deg = float(M_PI / 180.0);
  const Vec3 f_n(0.0f, 0.0f, -g_std);
  const Vec3 mag_n(20.0f, 0.0f, 45.0f);  // uT, 66 deg inclination
  const Vec3 acc_bias(0.05f, -0.05f, 0.08f);
  const Vec3 gyro_bias(0.005f, -0.005f, 0.005f);

  Mat3 R_nb = Eigen::AngleAxisf(1.0f, Vec3::UnitZ()).toRotationMatrix();
  float max_tilt = 0.0f, max_hdg_err = 0.0f, max_heave = 0.0f, max_startup_heave = 0.0f;

  for (int k = 0; k < 200 * 600; ++k) {
    const float t = k * dt;
    float yaw_rate = 0.0f;
    if (still_boot) {
      if (t > 30.0f && t < 300.0f && std::fmod(t, 40.0f) < 2.0f) yaw_rate = 45.0f * deg;
      if (t >= 300.0f) yaw_rate = 5.0f * deg;
    } else {
      yaw_rate = 3.0f * deg;
    }
    R_nb = R_nb * Eigen::AngleAxisf(yaw_rate * dt, Vec3::UnitZ()).toRotationMatrix();

    f.setMagBody(R_nb.transpose() * mag_n + noise(0.4f), true);
    Vec3 gyro_noise = noise(0.0016f);
    if (atoms3r_gyro) {
      const float w = 2.0f * float(M_PI) * t;
      gyro_noise = noise(0.2f * deg) + 0.3f * deg * Vec3(std::sin(10.0f * w), std::sin(10.3f * w + 1.0f),
                                                          std::sin(9.7f * w + 2.0f));
    }
    f.update(dt, Vec3(0.0f, 0.0f, yaw_rate) + gyro_bias + gyro_noise,
             R_nb.transpose() * f_n + acc_bias + noise(0.015f));

    const auto s = f.snapshot();
    if (t < 150.0f) {  // startup: attitude gains still ramping
      max_startup_heave = std::max(max_startup_heave, std::fabs(s.disp_zu.z()));
      continue;
    }
    const float tilt = std::hypot(s.euler_rad.x(), s.euler_rad.y()) / deg;
    const float hdg_err = wrap180((s.euler_rad.z() - std::atan2(R_nb(1, 0), R_nb(0, 0))) / deg);
    max_tilt = std::max(max_tilt, tilt);
    max_hdg_err = std::max(max_hdg_err, std::fabs(hdg_err));
    max_heave = std::max(max_heave, std::fabs(s.disp_zu.z()));
  }

  const bool ok = std::isfinite(max_tilt) && max_tilt < lim.tilt_deg &&
                  max_hdg_err < lim.hdg_deg && max_heave < lim.heave_m &&
                  max_startup_heave < lim.startup_heave_m;
  std::printf("NLO mag turns %-11s max|tilt|=%.2f deg max|hdg err|=%.2f deg max|heave|=%.3f m "
              "startup max|heave|=%.3f m %s\n",
              name, max_tilt, max_hdg_err, max_heave, max_startup_heave, ok ? "OK" : "FAIL");
  return ok;
}

// Before the bootstrap's motion-tolerant handoff attitude, this history
// handed over with about 20 deg of tilt error; the startup heading gain then
// tumbled the attitude and the heading spun through 360 deg.
static bool runHandled(float phase) {
  Fusion f;
  f.config().filter.expected_mag_norm = 0.0f;
  f.reset();

  std::mt19937 rng(11);
  std::normal_distribution<float> n01(0.0f, 1.0f);
  auto noise = [&](float s) { return Vec3(s * n01(rng), s * n01(rng), s * n01(rng)); };

  const float dt = 0.005f;
  const float deg = float(M_PI / 180.0);
  const Vec3 f_n(0.0f, 0.0f, -g_std);
  const Vec3 mag_n(20.0f, 0.0f, 45.0f);
  const float yaw = 1.0f;
  auto attitude = [&](float t) {
    const float roll = (25.0f * std::sin(2 * float(M_PI) * t / 1.5f + phase) +
                        8.0f * std::sin(2 * float(M_PI) * t / 0.7f + 0.5f)) * deg;
    const float pitch = (10.0f * std::sin(2 * float(M_PI) * t / 1.2f + 0.7f + phase) +
                         5.0f * std::sin(2 * float(M_PI) * t / 0.5f + 1.1f)) * deg;
    return Mat3(Eigen::AngleAxisf(yaw, Vec3::UnitZ()) * Eigen::AngleAxisf(pitch, Vec3::UnitY()) *
                Eigen::AngleAxisf(roll, Vec3::UnitX()));
  };
  auto accel = [&](float t) {  // horizontal handling acceleration, about 2 m/s^2
    const float w1 = 2 * float(M_PI) / 1.2f, w2 = 2 * float(M_PI) / 0.6f;
    const float a = -0.08f * w1 * w1 * std::sin(w1 * t + 0.3f + phase) -
                    0.03f * w2 * w2 * std::sin(w2 * t + 2.0f);
    const float az = 0.10f * 7.9f * 7.9f * std::sin(4.19f * t + phase) / 4.0f;
    return Vec3(a * std::cos(yaw), a * std::sin(yaw), az);
  };

  float max_hdg_err = 0.0f;
  bool saw_handoff = false;
  for (int k = 0; k < 200 * 120; ++k) {
    const float t = k * dt;
    const Mat3 R_nb = attitude(t);
    const Eigen::AngleAxisf step(R_nb.transpose() * attitude(t + 1e-4f));
    const Vec3 w = step.axis() * (step.angle() / 1e-4f);

    f.setMagBody(R_nb.transpose() * mag_n + noise(0.4f), true);
    f.update(dt, w + noise(0.0016f), R_nb.transpose() * (accel(t) + f_n) + noise(0.02f));

    const auto s = f.snapshot();
    if (!s.initialized) continue;
    saw_handoff = true;
    const Vec3 x_hat = Eigen::AngleAxisf(s.euler_rad.z(), Vec3::UnitZ()) *
                       Eigen::AngleAxisf(s.euler_rad.y(), Vec3::UnitY()) *
                       Eigen::AngleAxisf(s.euler_rad.x(), Vec3::UnitX()) * Vec3::UnitX();
    const Vec3 x_true = R_nb * Vec3::UnitX();
    const float hdg_err =
        wrap180((std::atan2(x_hat.y(), x_hat.x()) - std::atan2(x_true.y(), x_true.x())) / deg);
    max_hdg_err = std::max(max_hdg_err, std::fabs(hdg_err));
  }

  const bool ok = saw_handoff && std::isfinite(max_hdg_err) && max_hdg_err < 25.0f;
  std::printf("NLO mag turns handled boot phase=%.1f max|hdg err| after handoff=%.2f deg %s\n",
              phase, max_hdg_err, ok ? "OK" : "FAIL");
  return ok;
}

// A seed taken from this roll is wrong by about 0.26 deg/s, more than half
// the residual it is meant to remove.
static bool runSlowRollBoot() {
  Fusion f;
  f.config().filter.expected_mag_norm = 0.0f;
  f.reset();

  std::mt19937 rng(17);
  std::normal_distribution<float> n01(0.0f, 1.0f);
  auto noise = [&](float s) { return Vec3(s * n01(rng), s * n01(rng), s * n01(rng)); };

  const float dt = 0.005f;
  const float deg = float(M_PI / 180.0);
  const float w = 2.0f * float(M_PI) / 8.0f;
  const Vec3 f_n(0.0f, 0.0f, -g_std);
  const Vec3 mag_n(20.0f, 0.0f, 45.0f);
  const Vec3 gyro_bias(0.47f * deg, 0.0f, 0.0f);
  auto attitude = [&](float t) {
    return Mat3(Eigen::AngleAxisf(1.0f, Vec3::UnitZ()) *
                Eigen::AngleAxisf(0.5f * deg * std::sin(w * t), Vec3::UnitX()));
  };

  for (int k = 0; k < 200 * 10; ++k) {
    const float t = k * dt;
    const Mat3 R_nb = attitude(t);
    const Vec3 rate(0.5f * deg * w * std::cos(w * t), 0.0f, 0.0f);
    f.setMagBody(R_nb.transpose() * mag_n + noise(0.4f), true);
    f.update(dt, rate + gyro_bias + noise(0.0016f), R_nb.transpose() * f_n + noise(0.015f));
    const auto s = f.snapshot();
    if (!s.initialized) continue;
    // A seed is taken when the handed-over bias is not the bootstrap's own,
    // which is still near zero.
    const Vec3 b = s.tvg.gyro_bias_b;
    const float seed_error = (b - gyro_bias).norm() / deg;
    const bool seeded = b.norm() > 0.05f * deg;
    const bool ok = !seeded || seed_error < 0.1f;
    std::printf("NLO mag turns slow roll boot seeded=%d seed error=%.3f deg/s %s\n",
                int(seeded), seed_error, ok ? "OK" : "FAIL");
    return ok;
  }
  std::printf("NLO mag turns slow roll boot never initialized FAIL\n");
  return false;
}

int main() {
  const bool a = run("still boot", true, {2.0f, 5.0f, 0.5f, 0.3f});
  const bool a_device = run("still boot, AtomS3R gyro", true, {2.0f, 5.0f, 0.5f, 0.3f}, true);
  const bool a_roll = runSlowRollBoot();
  const bool b = run("moving boot", false, {8.0f, 20.0f, 5.0f, 1e9f});
  bool c = true;
  for (float phase : {0.0f, 1.3f, 2.6f, 3.9f, 5.2f}) c = runHandled(phase) && c;
  return (a && a_device && a_roll && b && c) ? 0 : 1;
}
