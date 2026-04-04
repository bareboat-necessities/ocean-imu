/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R Tilt-Compensated Compass + IMU Calibration Wizard (Mahony)
*/

#include <Arduino.h>
#include <M5Unified.h>

#define NO_BOSCH_API

// 1 = graphical compass by default, 0 = text UI by default
#ifndef COMPASS_UI_DEFAULT_GRAPHICS
#define COMPASS_UI_DEFAULT_GRAPHICS 1
#endif

// 0 = keep debug serial
// 1 = emit NMEA0183 (HDM + XDR + ROT) like pypilot
#ifndef COMPASS_SERIAL_NMEA
#define COMPASS_SERIAL_NMEA 1
#endif

#ifndef COMPASS_NMEA_TALKER
#define COMPASS_NMEA_TALKER "II"
#endif

#include "AtomS3R/AtomS3R_CompassAppBase.h"
#include "ahrs/Mahony_AHRS.h"

using namespace atoms3r_compass;

class MahonyBackend : public IAttitudeBackend {
 public:
  void reset() override {
    const float twoKp = 2.0f * 6.0f;
    const float twoKi = 2.0f * 0.0002f;
    mahony_AHRS_init(&mahony_, twoKp, twoKi);

    mahony_.q0 = 1.0f;
    mahony_.q1 = 0.0f;
    mahony_.q2 = 0.0f;
    mahony_.q3 = 0.0f;
    mahony_.integralFBx = mahony_.integralFBy = mahony_.integralFBz = 0.0f;
    valid_ = true;
  }

  void step(const CalibratedSample& s, AttitudeSolution& out) override {
    Vector3f a_att = -s.a_cal;
    const float an = a_att.norm();
    if (an > 1e-6f) a_att *= (ImuCalCfg::g_std / an);

    float pd = 0.0f, rd = 0.0f, yd = 0.0f;
    if (s.mag_ok && s.mag_fresh) {
      mahony_AHRS_update_mag(&mahony_, s.w_cal.x(), s.w_cal.y(), s.w_cal.z(), a_att.x(), a_att.y(), a_att.z(), s.m_unit.x(),
                             s.m_unit.y(), s.m_unit.z(), &pd, &rd, &yd, s.dt);
    } else {
      mahony_AHRS_update(&mahony_, s.w_cal.x(), s.w_cal.y(), s.w_cal.z(), a_att.x(), a_att.y(), a_att.z(), &pd, &rd, &yd, s.dt);
    }

    out = makeAttitudeFromQuat(mahony_.q1, mahony_.q2, mahony_.q3, mahony_.q0);
  }

  bool isValid() const override { return valid_; }

 private:
  Mahony_AHRS_Vars mahony_{};
  bool valid_ = false;
};

class MahonyCompassApp : public CompassAppBase {
 public:
  MahonyCompassApp()
      : CompassAppBase(std::make_unique<MahonyBackend>(), MagGateConfig{12, 0.02f, 10, 20}, "Mahony") {}
};

static MahonyCompassApp g_app;

void setup() { g_app.begin(); }
void loop() { g_app.tick(); }
