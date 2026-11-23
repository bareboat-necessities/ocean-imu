#pragma once

/*
  Copyright (c) 2026 Mikhail Grushinskiy
*/

#ifdef ARDUINO
#include <Arduino.h>
#endif

#include <math.h>
#include "nmea/NmeaChecksum.h"

// Send NMEA sentence WITHOUT the "*hh" part (we append checksum + CRLF).
static inline void nmea_send(const char* s_no_checksum) {
  const int cs = nmea0183_checksum(s_no_checksum);
  Serial.printf("%s*%02X\r\n", s_no_checksum, cs);
}

static inline uint8_t nmeaChecksumBody_(const char* body) {
  uint8_t cs = 0;
  while (*body) {
    cs ^= static_cast<uint8_t>(*body++);
  }
  return cs;
}

namespace nmea_detail {

// The original HDM/ROT/XDR buffers allow 81 bytes including '$'. The
// WAV and status buffers allow 95 body bytes, plus the leading '$'.
class Sentence {
public:
  explicit Sentence(bool shortTalker = false)
      : end_(data_ + 1), limit_(data_ + (shortTalker ? 96 : 81)), wholeBody_(shortTalker) {
    data_[0] = '$';
  }

  void append(const char* text, size_t maxChars = SIZE_MAX) {
    while (maxChars && end_ < limit_ && *text) {
      *end_++ = *text++;
      --maxChars;
    }
  }

  void decimal(float value, unsigned digits) {
    static constexpr float scales[] = {1, 10, 100, 1000, 10000};
    const float magnitude = fabsf(value) * scales[digits];
    char buf[48];  // sign, all binary32 integer digits, fraction and terminator
    // Below 2^23 every half-integer is representable in binary32. Rounding
    // the product can't cross the half-integer without landing exactly on it.
    // Let libc resolve that ambiguous case from the original input.
    if (magnitude < 8388608.0f) {
      uint32_t rounded = static_cast<uint32_t>(magnitude);
      const float fraction = magnitude - static_cast<float>(rounded);
      if (fraction != 0.5f) {
        rounded += fraction > 0.5f;
        char* p = buf + sizeof(buf) - 1;
        *p = '\0';
        for (unsigned i = 0; i < digits; ++i) {
          *--p = static_cast<char>('0' + rounded % 10);
          rounded /= 10;
        }
        if (digits) *--p = '.';
        do {
          *--p = static_cast<char>('0' + rounded % 10);
          rounded /= 10;
        } while (rounded);
        if (std::signbit(value)) *--p = '-';
        append(p);
        return;
      }
    }
    snprintf(buf, sizeof(buf), "%.*f", static_cast<int>(digits), static_cast<double>(value));
    append(buf);
  }

  void send() {
    *end_ = '\0';
    const uint8_t checksum = wholeBody_ ? nmeaChecksumBody_(data_ + 1) : nmea0183_checksum(data_);
    static constexpr char hex[] = "0123456789ABCDEF";
    *end_++ = '*';
    *end_++ = hex[checksum >> 4];
    *end_++ = hex[checksum & 15];
    *end_++ = '\r';
    *end_++ = '\n';
    Serial.write(reinterpret_cast<const uint8_t*>(data_), end_ - data_);
  }

private:
  char data_[96 + 5];
  char* end_;
  char* limit_;
  bool wholeBody_;
};

inline void sendValue(const char* talker, const char* prefix, float value, unsigned precision,
                      const char* suffix, bool shortTalker = false) {
  Sentence sentence(shortTalker);
  sentence.append(talker, shortTalker ? 2 : SIZE_MAX);
  sentence.append(prefix);
  sentence.decimal(value, precision);
  sentence.append(suffix);
  sentence.send();
}

}  // namespace nmea_detail

static inline float wrap360f_(float deg) {
  while (deg < 0.0f) deg += 360.0f;
  while (deg >= 360.0f) deg -= 360.0f;
  return deg;
}

// $--HDM,xxx.x,M*hh  (magnetic heading)
static inline void nmea_hdm(const char* talker2, float heading_deg) {
  nmea_detail::sendValue(talker2, "HDM,", wrap360f_(heading_deg), 1, ",M");
}

// $--ROT,x.x,A*hh  (rate of turn, degrees per minute; A=valid, V=invalid)
static inline void nmea_rot(const char* talker2, float rot_deg_per_min, bool valid) {
  nmea_detail::sendValue(talker2, "ROT,", rot_deg_per_min, 1, valid ? ",A" : ",V");
}

// $--XDR,A,x.x,D,PTCH*hh  (angular displacement in degrees)
// $--XDR,A,x.x,D,ROLL*hh
static inline void nmea_xdr_pitch_roll(const char* talker2, float pitch_deg, float roll_deg) {
  nmea_detail::sendValue(talker2, "XDR,A,", pitch_deg, 1, ",D,PTCH");
  nmea_detail::sendValue(talker2, "XDR,A,", roll_deg, 1, ",D,ROLL");
}

// $--XDR,D,x.x,M,DRT1*hh  (heave / vertical displacement, meters)
static inline void nmea_xdr_heave(const char* talker2, float heave_m) {
  nmea_detail::sendValue(talker2, "XDR,D,", heave_m, 4, ",M,DRT1");
}

// $--XDR,F,x.x,H,FRT1*hh  (dominant wave frequency, Hz)
static inline void nmea_xdr_freq(const char* talker2, float freq_hz) {
  nmea_detail::sendValue(talker2, "XDR,F,", freq_hz, 4, ",H,FRT1");
}

// $--XDR,V,x.xxx,M,VHSPD*hh  (vertical heave speed, meters/second, up-positive)
static inline void nmea_xdr_heave_speed(const char* talker2, float heave_speed_mps) {
  nmea_detail::sendValue(talker2, "XDR,V,", heave_speed_mps, 3, ",M,VHSPD");
}

// $--XDR,C,x.x,C,IMUT*hh  (IMU chip temperature, degrees C)
static inline void nmea_xdr_imu_temp(const char* talker2, float temp_c) {
  if (!std::isfinite(temp_c)) return;
  nmea_detail::sendValue(talker2, "XDR,C,", temp_c, 1, ",C,IMUT");
}

static inline void nmeaPrintBody_(const char* body) {
  if (strlen(body) <= 95) {
    nmea_detail::Sentence sentence(true);
    sentence.append(body);
    sentence.send();
    return;
  }
  const uint8_t cs = nmeaChecksumBody_(body);
  Serial.printf("$%s*%02X\r\n", body, static_cast<unsigned>(cs));
}

static inline void nmea_xdr_wave_axis_rel(const char* talker, float wave_axis_deg, bool valid)
{
  if (!valid || !std::isfinite(wave_axis_deg)) {
    return;
  }
  nmea_detail::sendValue(talker, "XDR,A,", wave_axis_deg, 1, ",D,WAVAXIS", true);
}

static inline void nmea_xdr_wave_direction_rel(const char* talker, float wave_dir_deg, bool valid)
{
  if (!valid || !std::isfinite(wave_dir_deg)) {
    return;
  }
  nmea_detail::sendValue(talker, "XDR,A,", wave_dir_deg, 1, ",D,WAVDIR", true);
}

static inline void nmea_txt_wave_direction_sign(const char* talker, int raw_sign, int polarity, bool valid)
{
  if (!valid) {
    return;
  }
  char body[96];
  snprintf(body, sizeof(body), "%.2sTXT,01,01,00,WAVSGN=%d POL=%d", talker, raw_sign, polarity);
  nmeaPrintBody_(body);
}

static inline void nmea_txt_wave_direction_confidence(const char* talker, float confidence_pct, bool valid)
{
  if (!valid || !std::isfinite(confidence_pct)) {
    return;
  }
  if (confidence_pct < 0.0f) confidence_pct = 0.0f;
  if (confidence_pct > 100.0f) confidence_pct = 100.0f;
  nmea_detail::sendValue(talker, "TXT,01,01,00,WAVCONF=", confidence_pct, 0, "", true);
}

static inline void nmea_txt_ins_status(const char* talker,
                                       bool attitude_good,
                                       bool heave_good,
                                       bool mag_locked,
                                       bool gyro_bias_learning,
                                       bool acc_bias_estimating,
                                       float imu_hz,
                                       float mag_hz)
{
  nmea_detail::Sentence sentence(true);
  sentence.append(talker, 2);
  sentence.append("TXT,01,01,00,ATT=");
  sentence.append(attitude_good ? "Y HEV=" : "N HEV=");
  sentence.append(heave_good ? "Y MAG=" : "N MAG=");
  sentence.append(mag_locked ? "Y GB=" : "N GB=");
  sentence.append(gyro_bias_learning ? "Y AB=" : "N AB=");
  sentence.append(acc_bias_estimating ? "Y IHz=" : "N IHz=");
  sentence.decimal(imu_hz, 1);
  sentence.append(" MHz=");
  sentence.decimal(mag_hz, 1);
  sentence.send();
}
