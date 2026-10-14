#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Shared raw-IMU spike rejection for the AtomS3R compass and marine INS
  sketches. Sporadic single-sample BMI270 accelerometer spikes are detected
  on |a| and, to catch sign flips that keep |a|, on the vertical (body z)
  component. A rejected read is replaced as a whole: the accelerometer AND the
  gyroscope vectors of the same observed window sample are substituted,
  because a corrupted accelerometer read makes the gyro read of the same
  instant suspect too.

  Apply it to the raw mapped sample before any other processing
  (calibration, fusion, gyro-bias averaging). Build with
  -DATOMS3R_IMU_SPIKE_FILTER=0 to disable it in every sketch.
*/

#include "avg/TimeAwareSpikeFilter3D.h"

#ifndef ATOMS3R_IMU_SPIKE_FILTER
  #define ATOMS3R_IMU_SPIKE_FILTER 1
#endif

namespace atoms3r_ical {

struct ImuSpikeCfg {
  static constexpr bool  enabled        = ATOMS3R_IMU_SPIKE_FILTER != 0;
  static constexpr int   window         = 5;       // moving window size
  static constexpr float jerk_mps3      = 1800.0f; // |a| derivative threshold
  static constexpr int   vert_axis      = 2;       // body-NED z also checked alone
  static constexpr float vert_jerk_mps3 = 1800.0f; // a_z derivative threshold
};

class ImuSpikeFilter {
 public:
  void reset() { filter_.reset(); }

  // a, w: raw accelerometer / gyro vectors of one IMU read (any 3-vector type
  // with x()/y()/z() lvalue accessors), overwritten in place on a spike.
  // Returns true when the read was flagged and replaced.
  template <class Vec3>
  bool apply(Vec3& a, Vec3& w, float dt) {
    if (!ImuSpikeCfg::enabled) return false;
    return filter_.filterWithDelta(a.x(), a.y(), a.z(), w.x(), w.y(), w.z(), dt);
  }

 private:
  TimeAwareSpikeFilter3D filter_{ImuSpikeCfg::window, ImuSpikeCfg::jerk_mps3,
                                 ImuSpikeCfg::vert_axis, ImuSpikeCfg::vert_jerk_mps3};
};

}  // namespace atoms3r_ical
