#pragma once

#include "AtomS3R/AtomS3R_SerialBuffer.h"
#include "AtomS3R/AtomS3R_OutputTiming.h"
#include "nmea/NmeaCompass.h"
#include "wave_dir/WaveDirectionReport.h"

namespace atoms3r_ical {

// Hold publication-time values while later IMU ticks send them in wire order.
class MarineNmeaBatch {
 public:
  struct Status {
    bool heave_good, mag_locked, gyro_bias_learning, acc_bias_estimating;
    float imu_hz, mag_hz;
  };

  void begin(float heave_m, float speed_mps, float rot_dpm, bool valid) {
    heave_m_ = heave_m;
    speed_mps_ = speed_mps;
    rot_dpm_ = rot_dpm;
    valid_ = valid;
    phase_ = Phase::Heave;
    wave_due_ = temp_due_ = status_due_ = false;
  }

  // Defer the two tilt sentences when the transport needs smaller bursts.
  void addTilt(float pitch_deg, float roll_deg) {
    pitch_deg_ = pitch_deg;
    roll_deg_ = roll_deg;
    phase_ = Phase::Tilt;
  }

  void beginCompass(float pitch_deg, float roll_deg, float rot_dpm, bool valid) {
    addTilt(pitch_deg, roll_deg);
    rot_dpm_ = rot_dpm;
    valid_ = valid;
    phase_ = Phase::CompassTilt;
  }

  void addWave(const wave_direction::WaveDirectionReport& wave, bool live) {
    wave_ = wave;
    live_ = live;
    wave_due_ = true;
  }

  void addTemperature(float temp_c) {
    temp_c_ = temp_c;
    temp_due_ = true;
  }

  void addStatus(const Status& status) {
    status_ = status;
    status_due_ = true;
  }

  bool pending() const { return phase_ != Phase::Idle; }
  void reset() {
    phase_ = Phase::Idle;
    output_.reset();
  }

  SerialBuffer<2 * nmea_detail::kMaxSentenceBytes>& output() { return output_; }

  template<class Port>
  bool readyForSnapshot(Port& port, const char* talker, uint32_t now, uint32_t& last) {
    const bool due = now - last >= atoms3r_timing::NMEA_SERIAL_MS;
    if (!due && !pending() && !output_.pending()) return false;
    if (!port) {
      reset();
      last = now;
      return false;
    }
    if (output_.pending()) {
      output_.drain(port);
      return false;
    }
    if (port.availableForWrite() <= 0) return false;
    if (sendNext(talker, output_)) {
      output_.drain(port);
      return false;
    }
    last = now;
    return true;
  }

  // Called once per sampling slot. Send at most two sentences per call.
  template<class Out = decltype(Serial)>
  bool sendNext(const char* talker, Out& out = Serial) {
    switch (phase_) {
      case Phase::Idle:
        return false;
      case Phase::CompassTilt:
      case Phase::Tilt:
        nmea_xdr_pitch_roll(talker, pitch_deg_, roll_deg_, out);
        phase_ = phase_ == Phase::CompassTilt ? Phase::Turn : Phase::Heave;
        break;
      case Phase::Heave:
        nmea_xdr_heave(talker, heave_m_, out);
        nmea_xdr_heave_speed(talker, speed_mps_, out);
        phase_ = wave_due_ ? Phase::WaveAngles : afterWave();
        break;
      case Phase::WaveAngles:
        nmea_xdr_wave_axis_rel(talker, wave_.axis_deg, wave_.axis_ok, out);
        nmea_xdr_wave_direction_rel(talker, wave_.dir_deg, wave_.dir_ok, out);
        phase_ = Phase::WaveQuality;
        break;
      case Phase::WaveQuality:
        nmea_txt_wave_direction_sign(talker, wave_.sign_raw, wave_.sign_polarity, wave_.sign_ok, out);
        nmea_txt_wave_direction_confidence(talker, wave_.conf_pct, live_, out);
        phase_ = afterWave();
        break;
      case Phase::TemperatureAndStatus:
        if (temp_due_) nmea_xdr_imu_temp(talker, temp_c_, out);
        if (status_due_) {
          nmea_txt_ins_status(talker, valid_, status_.heave_good, status_.mag_locked,
                              status_.gyro_bias_learning, status_.acc_bias_estimating,
                              status_.imu_hz, status_.mag_hz, out);
        }
        phase_ = Phase::Turn;
        break;
      case Phase::Turn:
        nmea_rot(talker, rot_dpm_, valid_, out);
        phase_ = Phase::Idle;
        break;
    }
    return true;
  }

 private:
  enum class Phase { Idle, CompassTilt, Tilt, Heave, WaveAngles, WaveQuality, TemperatureAndStatus, Turn };
  Phase afterWave() const {
    return temp_due_ || status_due_ ? Phase::TemperatureAndStatus : Phase::Turn;
  }

  SerialBuffer<2 * nmea_detail::kMaxSentenceBytes> output_;
  Phase phase_ = Phase::Idle;
  float pitch_deg_ = 0, roll_deg_ = 0;
  float heave_m_ = 0, speed_mps_ = 0, rot_dpm_ = 0, temp_c_ = NAN;
  bool valid_ = false, live_ = false;
  bool wave_due_ = false, temp_due_ = false, status_due_ = false;
  wave_direction::WaveDirectionReport wave_{};
  Status status_{};
};

}  // namespace atoms3r_ical
