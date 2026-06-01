#pragma once

#include <Arduino.h>

#include <math.h>
#include "nmea/NmeaChecksum.h"

// Send NMEA sentence WITHOUT the "*hh" part (we append checksum + CRLF).
static inline void nmea_send(const char* s_no_checksum) {
  const int cs = nmea0183_checksum(s_no_checksum);
  Serial.printf("%s*%02X\r\n", s_no_checksum, cs);
}

static inline float wrap360f_(float deg) {
  while (deg < 0.0f) deg += 360.0f;
  while (deg >= 360.0f) deg -= 360.0f;
  return deg;
}

// $--HDM,xxx.x,M*hh  (magnetic heading)
static inline void nmea_hdm(const char* talker2, float heading_deg) {
  char s[82];
  heading_deg = wrap360f_(heading_deg);
  snprintf(s, sizeof(s), "$%sHDM,%.1f,M", talker2, (double)heading_deg);
  nmea_send(s);
}

// $--ROT,x.x,A*hh  (rate of turn, degrees per minute; A=valid, V=invalid)
static inline void nmea_rot(const char* talker2, float rot_deg_per_min, bool valid) {
  char s[82];
  // pypilot typically sends ROT in deg/min; sign may need flip depending on your axis convention
  snprintf(s, sizeof(s), "$%sROT,%.1f,%c", talker2, (double)rot_deg_per_min, valid ? 'A' : 'V');
  nmea_send(s);
}

// $--XDR,A,x.x,D,PTCH*hh  (angular displacement in degrees)
// $--XDR,A,x.x,D,ROLL*hh
static inline void nmea_xdr_pitch_roll(const char* talker2, float pitch_deg, float roll_deg) {
  char s[82];
  // Keep it short to stay under 82 chars
  snprintf(s, sizeof(s), "$%sXDR,A,%.1f,D,PTCH", talker2, (double)pitch_deg);
  nmea_send(s);
  snprintf(s, sizeof(s), "$%sXDR,A,%.1f,D,ROLL", talker2, (double)roll_deg);
  nmea_send(s);
}

// $--XDR,D,x.x,M,DRT1*hh  (heave / vertical displacement, meters)
static inline void nmea_xdr_heave(const char* talker2, float heave_m) {
  char s[82];
  snprintf(s, sizeof(s), "$%sXDR,D,%.4f,M,DRT1", talker2, (double)heave_m);
  nmea_send(s);
}

// $--XDR,F,x.x,H,FRT1*hh  (dominant wave frequency, Hz)
static inline void nmea_xdr_freq(const char* talker2, float freq_hz) {
  char s[82];
  snprintf(s, sizeof(s), "$%sXDR,F,%.4f,H,FRT1", talker2, (double)freq_hz);
  nmea_send(s);
}

// $--XDR,V,x.xxx,M,VHSPD*hh  (vertical heave speed, meters/second, up-positive)
static inline void nmea_xdr_heave_speed(const char* talker2, float heave_speed_mps) {
  char s[82];
  snprintf(s, sizeof(s), "$%sXDR,V,%.3f,M,VHSPD", talker2, (double)heave_speed_mps);
  nmea_send(s);
}

// $--XDR,C,x.x,C,IMUT*hh  (IMU chip temperature, degrees C)
static inline void nmea_xdr_imu_temp(const char* talker2, float temp_c) {
  if (!std::isfinite(temp_c)) {
    return;
  }
  char s[82];
  snprintf(s, sizeof(s), "$%sXDR,C,%.1f,C,IMUT", talker2, (double)temp_c);
  nmea_send(s);
}

static inline uint8_t nmeaChecksumBody_(const char* body) {
  uint8_t cs = 0;
  while (*body) {
    cs ^= static_cast<uint8_t>(*body++);
  }
  return cs;
}

static inline void nmeaPrintBody_(const char* body) {
  const uint8_t cs = nmeaChecksumBody_(body);
  Serial.printf("$%s*%02X\r\n", body, static_cast<unsigned>(cs));
}

static inline void nmea_xdr_wave_axis_rel(const char* talker, float wave_axis_deg, bool valid)
{
  if (!valid || !std::isfinite(wave_axis_deg)) {
    return;
  }
  char body[96];
  snprintf(body, sizeof(body), "%.2sXDR,A,%.1f,D,WAVAXIS", talker, static_cast<double>(wave_axis_deg));
  nmeaPrintBody_(body);
}

static inline void nmea_xdr_wave_direction_rel(const char* talker, float wave_dir_deg, bool valid)
{
  if (!valid || !std::isfinite(wave_dir_deg)) {
    return;
  }
  char body[96];
  snprintf(body, sizeof(body), "%.2sXDR,A,%.1f,D,WAVDIR", talker, static_cast<double>(wave_dir_deg));
  nmeaPrintBody_(body);
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
  char body[96];
  if (confidence_pct < 0.0f) confidence_pct = 0.0f;
  if (confidence_pct > 100.0f) confidence_pct = 100.0f;
  snprintf(body, sizeof(body), "%.2sTXT,01,01,00,WAVCONF=%.0f", talker, static_cast<double>(confidence_pct));
  nmeaPrintBody_(body);
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
  char body[96];
  snprintf(body, sizeof(body),
           "%.2sTXT,01,01,00,ATT=%c HEV=%c MAG=%c GB=%c AB=%c IHz=%.1f MHz=%.1f",
           talker,
           attitude_good ? 'Y' : 'N',
           heave_good ? 'Y' : 'N',
           mag_locked ? 'Y' : 'N',
           gyro_bias_learning ? 'Y' : 'N',
           acc_bias_estimating ? 'Y' : 'N',
           static_cast<double>(imu_hz),
           static_cast<double>(mag_hz));
  nmeaPrintBody_(body);
}
