// Check each wire write including the legacy sentence-length and checksum rules.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <limits>
#include <random>
#include <string>

struct {
  std::deque<std::string> writes;
  void write(const uint8_t* p, size_t n) { writes.emplace_back(reinterpret_cast<const char*>(p), n); }
  void printf(const char* fmt, ...) {
    char line[256];
    va_list args;
    va_start(args, fmt);
    const int n = vsnprintf(line, sizeof(line), fmt, args);
    va_end(args);
    assert(n >= 0 && static_cast<size_t>(n) < sizeof(line));
    writes.emplace_back(line, n);
  }
} Serial;
#include "nmea/NmeaCompass.h"

void expect(bool wholeBody, const char* fmt, ...) {
  char line[103];
  va_list args;
  va_start(args, fmt);
  vsnprintf(line, wholeBody ? 97 : 82, fmt, args);
  va_end(args);
  unsigned cs = 0;
  for (const char* p = line + 1; *p; ++p) {
    if (!wholeBody && (*p == '*' || *p == '\r')) break;
    cs ^= static_cast<uint8_t>(*p);
  }
  snprintf(line + strlen(line), 7, "*%02X\r\n", cs);
  assert(!Serial.writes.empty() && Serial.writes.front() == line);
  Serial.writes.pop_front();
}

void checkValues(float value) {
  nmea_rot("II", value, true);
  expect(false, "$IIROT,%.1f,A", value);
  nmea_xdr_heave_speed("II", value);
  expect(false, "$IIXDR,V,%.3f,M,VHSPD", value);
  nmea_xdr_heave("II", value);
  expect(false, "$IIXDR,D,%.4f,M,DRT1", value);
  nmea_txt_wave_direction_confidence("II", value, true);
  if (std::isfinite(value)) expect(true, "$IITXT,01,01,00,WAVCONF=%.0f", std::clamp(value, 0.0f, 100.0f));
  assert(Serial.writes.empty());
}

int main() {
  // Exactly halfway, either side, signed zero, nonfinite and large inputs.
  for (float x : {0.0f, -0.0f, 2.5f, -2.5f, 127.5f, 0.25f, -31.25f, 0.12345f,
                  0.1235f, 83.35f, 838860.75f, std::numeric_limits<float>::max()}) {
    checkValues(x);
    checkValues(nextafterf(x, INFINITY));
    checkValues(nextafterf(x, -INFINITY));
  }
  checkValues(INFINITY);
  checkValues(-INFINITY);
  checkValues(NAN);
  checkValues(-NAN);
  std::mt19937 random(4801);
  for (int i = 0; i < 16000; ++i) {
    uint32_t bits = random();
    float value;
    std::memcpy(&value, &bits, sizeof(value));
    checkValues(value);
  }

  for (const char* talker : {"II", "ABCDEFG", "I*", "I\r", "",
                              "a very long talker which makes the HDM buffer truncate the sentence"
                              " before the measurement"}) {
    nmea_hdm(talker, -3.0f);
    expect(false, "$%sHDM,357.0,M", talker);
    nmea_xdr_wave_axis_rel(talker, 48.33f, true);
    expect(true, "$%.2sXDR,A,48.3,D,WAVAXIS", talker);
    nmea_xdr_pitch_roll(talker, -0.0f, -1.25f);
    expect(false, "$%sXDR,A,%.1f,D,PTCH", talker, -0.0);
    expect(false, "$%sXDR,A,%.1f,D,ROLL", talker, -1.25);
  }
  // %.2s accepts a two-byte talker with no terminating zero.
  const char boundedTalker[2] = {'I', 'N'};
  nmea_xdr_wave_direction_rel(boundedTalker, 42.0f, true);
  expect(true, "$INXDR,A,42.0,D,WAVDIR");
  nmea_xdr_freq("BB", 3.112f);
  expect(false, "$BBXDR,F,%.4f,H,FRT1", 3.112f);
  nmea_xdr_imu_temp("II", 35.1f);
  expect(false, "$IIXDR,C,35.1,C,IMUT");
  nmea_xdr_wave_direction_rel("II", 187.9f, true);
  expect(true, "$IIXDR,A,187.9,D,WAVDIR");
  nmea_txt_wave_direction_sign("II", -1, 1, true);
  expect(true, "$IITXT,01,01,00,WAVSGN=-1 POL=1");
  for (float v : {200.0f, NAN, std::numeric_limits<float>::max()}) {
    nmea_txt_ins_status("I*", true, false, true, false, true, v, -v);
    expect(true, "$I*TXT,01,01,00,ATT=Y HEV=N MAG=Y GB=N AB=Y IHz=%.1f MHz=%.1f", v, -v);
  }
  nmea_xdr_imu_temp("II", NAN);
  nmea_xdr_wave_direction_rel("II", NAN, true);
  nmea_xdr_wave_axis_rel("II", 5.0f, false);
  nmea_txt_wave_direction_sign("II", 1, -1, false);
  nmea_txt_wave_direction_confidence("II", 2.0f, false);
  assert(Serial.writes.empty());
}
