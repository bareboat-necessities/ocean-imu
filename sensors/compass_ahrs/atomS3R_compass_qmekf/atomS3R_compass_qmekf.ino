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
    // Full calibrated uT vector: normalization would remove a measurement
    // degree of freedom and invalidate a 3-DOF vector NIS test.
    Vector3f sigma_m; sigma_m <<     0.80f,     0.80f,      0.80f;

    if (mekf_) {
      mekf_->~QuaternionMEKF<float, true>();
      mekf_ = nullptr;
    }

    // sigma_g is a gyro white-noise density (rad/s/sqrt(Hz)); the gyro-bias
    // random-walk density 2e-7 rad^2/s^3 equals the former 1e-9 per 200 Hz step.
    mekf_ = new (storage_) QuaternionMEKF<float, true>(sigma_a, sigma_g, sigma_m, 0.5f, 1e-2f, 2e-7f);
    inited_ = false;mag_reference_ready_=false;mag_elapsed_s_=0;mag_reference_world_.setZero();
  }

  void configureNoise(const atoms3r_ical::RuntimeCals& runtime) override {
    if(mekf_ && runtime.configureProcessNoise(*mekf_))mekf_->set_Racc_std(runtime.noiseAccelStd());
  }

  void step(const CalibratedSample& s, AttitudeSolution& out) override {
    if (!inited_) {
      Vector3f a_init = s.a_cal;
      const float an0 = a_init.norm();
      if (an0 > 1e-6f) a_init *= (ImuCalCfg::g_cal_local / an0);  // reference at the physical static norm

      if (s.mag_ok)
        mekf_->initialize_from_acc_mag(a_init, s.m_cal);
      else
        mekf_->initialize_from_acc(a_init);

      inited_ = true;
      if(s.mag_ok)rememberMagReference_(s.m_cal);
    }

    mekf_->time_update(s.w_cal, s.dt);

    Vector3f a_att = s.a_cal;
    const float an = a_att.norm();
    if (an > 1e-6f) a_att *= (ImuCalCfg::g_cal_local / an);
    mekf_->measurement_update_acc_only(a_att);

    if(std::isfinite(s.dt) && s.dt>0)mag_elapsed_s_+=s.dt;
    if (s.mag_ok && s.mag_fresh) {
      // First magnetic acquisition establishes the yaw gauge. Once acquired,
      // the trusted reference must be used for NIS BEFORE any adaptation.
      if(!mag_reference_ready_) {
        mekf_->initialize_from_acc_mag(a_att,s.m_cal);
        rememberMagReference_(s.m_cal);
      }
      mekf_->measurement_update_mag_only(s.m_cal);
      const auto& health=mekf_->magneticHealth();
      if(health.accepted && health.health==ocean_imu::MagneticInnovationGate<float>::Health::HEALTHY) {
        const auto qc=mekf_->quaternion();
        const Vector3f mw=Eigen::Quaternionf(qc(3),qc(0),qc(1),qc(2))*s.m_cal;
        const Vector3f target(std::hypot(mw.x(),mw.y()),0,mw.z());
        // Retain gradual dip correction for boot noise, but an outlier or a
        // recovery packet cannot redefine its own predicted measurement.
        const float alpha=1.f-std::exp(-std::min(.2f,mag_elapsed_s_)/10.f);
        mag_reference_world_+=alpha*(target-mag_reference_world_);
        mekf_->set_mag_world_ref(mag_reference_world_);
      }
      mag_elapsed_s_=0;
    }

    const auto q = mekf_->quaternion();
    out = makeAttitudeFromQuat(q(0), q(1), q(2), q(3));
  }

  bool isValid() const override { return inited_; }
  const auto& magneticHealth() const {return mekf_->magneticHealth();}

 private:
  void rememberMagReference_(const Vector3f& measurement) {
    const auto qc=mekf_->quaternion();
    mag_reference_world_=Eigen::Quaternionf(qc(3),qc(0),qc(1),qc(2))*measurement;
    mag_reference_ready_=true;mag_elapsed_s_=0;
  }
  Vector3f mag_reference_world_=Vector3f::Zero();
  bool mag_reference_ready_=false;float mag_elapsed_s_=0;
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
