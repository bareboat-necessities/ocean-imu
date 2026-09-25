// Device configuration of the NLO (no GNSS, 3D magnetometer, as the
// atomS3R_ins_nlo sketch runs it) on a bench history the wave simulator does
// not exercise: a gyro residual bias of 0.3 deg/s per axis, an accelerometer
// bias, and the body yawing. Checks both startups the adapter sees:
//
//   still boot : the device lies still while the bootstrap runs, so the
//                stillness-gated bias seed applies; then it is turned 90 deg
//                every 40 s, and afterwards spun continuously at 5 deg/s.
//   moving boot: it is already yawing at 3 deg/s at boot, so no seed is taken
//                and the bias has to be learned by the observer while it
//                turns. That takes minutes, so only boundedness is checked.
//
// Before the horizontal theta floor and the bias seed, the still-boot history
// reached 11 deg of tilt, 24 deg of heading error and 23 m of heave, and the
// moving-boot history 30 deg, 180 deg and 200 m, within ten minutes.
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
  float tilt_deg, hdg_deg, heave_m;
};

static bool run(const char* name, bool still_boot, Limits lim) {
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
  const Vec3 acc_bias(0.05f, -0.05f, 0.025f);
  const Vec3 gyro_bias(0.005f, -0.005f, 0.005f);

  Mat3 R_nb = Eigen::AngleAxisf(1.0f, Vec3::UnitZ()).toRotationMatrix();
  float max_tilt = 0.0f, max_hdg_err = 0.0f, max_heave = 0.0f;

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
    f.update(dt, Vec3(0.0f, 0.0f, yaw_rate) + gyro_bias + noise(0.0016f),
             R_nb.transpose() * f_n + acc_bias + noise(0.015f));

    if (t < 150.0f) continue;  // startup: attitude gains still ramping
    const auto s = f.snapshot();
    const float tilt = std::hypot(s.euler_rad.x(), s.euler_rad.y()) / deg;
    const float hdg_err = wrap180((s.euler_rad.z() - std::atan2(R_nb(1, 0), R_nb(0, 0))) / deg);
    max_tilt = std::max(max_tilt, tilt);
    max_hdg_err = std::max(max_hdg_err, std::fabs(hdg_err));
    max_heave = std::max(max_heave, std::fabs(s.disp_zu.z()));
  }

  const bool ok = std::isfinite(max_tilt) && max_tilt < lim.tilt_deg &&
                  max_hdg_err < lim.hdg_deg && max_heave < lim.heave_m;
  std::printf("NLO mag turns %-11s max|tilt|=%.2f deg max|hdg err|=%.2f deg max|heave|=%.3f m %s\n",
              name, max_tilt, max_hdg_err, max_heave, ok ? "OK" : "FAIL");
  return ok;
}

int main() {
  const bool a = run("still boot", true, {2.0f, 5.0f, 0.5f});
  const bool b = run("moving boot", false, {8.0f, 20.0f, 5.0f});
  return (a && b) ? 0 : 1;
}
