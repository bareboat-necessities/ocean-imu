/*
  Copyright 2026, Mikhail Grushinskiy

  AtomS3R Tilt-Compensated Compass + IMU Calibration Wizard (qMEKF)
*/

#include <Arduino.h>
#include <M5Unified.h>

#include <ArduinoOceanImu.h>

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
#include "ahrs/KalmanQMEKF.h"

using namespace atoms3r_compass;

class QmekfBackend : public IAttitudeBackend {
 public:
  void reset() override {
    const float g = ImuCalCfg::g_std;

    Vector3f sigma_a; sigma_a <<  0.06f * g,  0.06f * g,   0.06f * g;
    Vector3f sigma_g; sigma_g <<    0.0030f,    0.0030f,     0.0030f;
    Vector3f sigma_m; sigma_m <<     0.020f,     0.020f,      0.020f;

    if (mekf_) {
      mekf_->~QuaternionMEKF<float, true>();
      mekf_ = nullptr;
    }

    mekf_ = new (storage_) QuaternionMEKF<float, true>(sigma_a, sigma_g, sigma_m, 0.5f, 1e-2f, 1e-9f);
    inited_ = false;
  }

  void step(const CalibratedSample& s, AttitudeSolution& out) override {
    if (!inited_) {
      Vector3f a_init = s.a_cal;
      const float an0 = a_init.norm();
      if (an0 > 1e-6f) a_init *= (ImuCalCfg::g_std / an0);

      if (s.mag_ok)
        mekf_->initialize_from_acc_mag(a_init, s.m_unit);
      else
        mekf_->initialize_from_acc(a_init);

      inited_ = true;
    }

    mekf_->time_update(s.w_cal, s.dt);

    Vector3f a_att = s.a_cal;
    const float an = a_att.norm();
    if (an > 1e-6f) a_att *= (ImuCalCfg::g_std / an);
    mekf_->measurement_update_acc_only(a_att);

    if (s.mag_ok && s.mag_fresh) mekf_->measurement_update_mag_only(s.m_unit);

    const auto q = mekf_->quaternion();
    out = makeAttitudeFromQuat(q(0), q(1), q(2), q(3));
  }

  bool isValid() const override { return inited_; }

 private:
  alignas(QuaternionMEKF<float, true>) uint8_t storage_[sizeof(QuaternionMEKF<float, true>)];
  QuaternionMEKF<float, true>* mekf_ = nullptr;
  bool inited_ = false;
};

class QmekfCompassApp : public CompassAppBase {
 public:
  QmekfCompassApp() : CompassAppBase(std::make_unique<QmekfBackend>(), MagGateConfig{35, 0.001f, 20}, "QMEKF") {}
};

static QmekfCompassApp g_app;

void setup() { g_app.begin(); }
void loop() { g_app.tick(); }
