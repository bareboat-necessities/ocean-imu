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
    const float g = ImuCalCfg::g_std;  // noise below is specified in nominal g (unit conversion)

    Vector3f sigma_a; sigma_a <<  0.06f * g,  0.06f * g,   0.06f * g;
    Vector3f sigma_g; sigma_g <<    0.0030f,    0.0030f,     0.0030f;
    Vector3f sigma_m; sigma_m << kSigmaMag, kSigmaMag, kSigmaMag;

    if (mekf_) {
      mekf_->~QuaternionMEKF<float, true>();
      mekf_ = nullptr;
    }

    mekf_ = new (storage_) QuaternionMEKF<float, true>(sigma_a, sigma_g, sigma_m, 0.5f, 1e-2f, 1e-9f);
    inited_ = false;
  }

  void step(const CalibratedSample& s, AttitudeSolution& out) override {
    if (!inited_) {
      // A compass needs its magnetic reference: initializing from the
      // accelerometer alone would leave the filter's field reference unset.
      if (!s.mag_ok) {
        out = AttitudeSolution{};
        return;
      }
      Vector3f a_init = s.a_cal;
      const float an0 = a_init.norm();
      if (an0 > 1e-6f) a_init *= (ImuCalCfg::g_cal_local / an0);  // reference at the physical static norm
      mekf_->initialize_from_acc_mag(a_init, s.m_unit);
      inited_ = true;
    }

    mekf_->time_update(s.w_cal, s.dt);

    Vector3f a_att = s.a_cal;
    const float an = a_att.norm();
    if (an > 1e-6f) a_att *= (ImuCalCfg::g_cal_local / an);
    mekf_->measurement_update_acc_only(a_att);

    // Heading-only magnetic correction, as in a tilt-compensated compass: the
    // field's dip never pulls roll/pitch, which stay with the accelerometer.
    if (s.mag_ok && s.mag_fresh) mekf_->measurement_update_mag_heading(s.m_unit, kSigmaMag);

#if !COMPASS_SERIAL_NMEA
    if (++debug_count_ >= 200) {  // about once per second at 200 Hz
      debug_count_ = 0;
      const Vector3f b = mekf_->gyroscope_bias() * RAD_TO_DEG;
      Serial.printf("[QMEKF] gyro bias estimate [%+.3f %+.3f %+.3f] deg/s\n", (double)b.x(), (double)b.y(),
                    (double)b.z());
    }
#endif

    const auto q = mekf_->quaternion();
    out = makeAttitudeFromQuat(q(0), q(1), q(2), q(3));
  }

  bool isValid() const override { return inited_; }

 private:
  static constexpr float kSigmaMag = 0.020f;  // per component of the unit field vector
  int debug_count_ = 0;
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
