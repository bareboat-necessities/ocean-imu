// Run the firmware's diagnostic publisher with a finite byte ring and USB sink.
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

struct Print {
  virtual size_t write(uint8_t) = 0;
  virtual size_t write(const uint8_t*, size_t) = 0;
  template<class... Args> size_t printf(const char* fmt, Args... args) {
    char line[1024];
    int n = std::snprintf(line, sizeof(line), fmt, args...);
    assert(n >= 0 && static_cast<size_t>(n) < sizeof(line));
    return write(reinterpret_cast<const uint8_t*>(line), n);
  }
  size_t println() { return printf("\r\n"); }
};

struct : Print {
  int room = 256;
  size_t limit = SIZE_MAX;
  bool connected = true;
  std::string bytes;
  explicit operator bool() const { return connected; }
  int availableForWrite() const { return room; }
  size_t write(uint8_t c) override { return write(&c, 1); }
  size_t write(const uint8_t* p, size_t n) override {
    assert(connected && n <= static_cast<size_t>(room));
    n = std::min(n, limit);
    bytes.append(reinterpret_cast<const char*>(p), n);
    room -= static_cast<int>(n);
    return n;
  }
} Serial;

struct StaticRingbuffer_t {
  uint8_t* data;
  size_t capacity, head, used, held;
};
using RingbufHandle_t = StaticRingbuffer_t*;
enum { RINGBUF_TYPE_BYTEBUF, pdTRUE };
RingbufHandle_t xRingbufferCreateStatic(size_t cap, int, uint8_t* data, RingbufHandle_t q) {
  *q = {data, cap, 0, 0, 0};
  return q;
}
size_t xRingbufferGetCurFreeSize(RingbufHandle_t q) { return q->capacity - q->used; }
int xRingbufferSend(RingbufHandle_t q, const uint8_t* p, size_t n, int) {
  if (n > xRingbufferGetCurFreeSize(q)) return 0;
  for (size_t i = 0; i < n; ++i) q->data[(q->head + q->used++) % q->capacity] = p[i];
  return pdTRUE;
}
void* xRingbufferReceiveUpTo(RingbufHandle_t q, size_t* n, int, size_t limit) {
  assert(q->held == 0);
  *n = q->held = std::min({q->used, q->capacity - q->head, limit});
  return *n ? q->data + q->head : nullptr;
}
void vRingbufferReturnItem(RingbufHandle_t q, void* item) {
  assert(item == q->data + q->head && q->held > 0);
  q->head = (q->head + q->held) % q->capacity;
  q->used -= q->held;
  q->held = 0;
}

uint32_t micros() { return 0; }
void delayMicroseconds(uint32_t) {}
#include "AtomS3R/AtomS3R_DiagnosticSerial.h"
#include "AtomS3R/AtomS3R_ImuCommon.h"

template<size_t Capacity>
void check() {
  atoms3r_ical::DiagnosticSerial<Capacity> log;
  log.begin();
  Serial.bytes.clear();
  Serial.room = 4096;
  Serial.limit = SIZE_MAX;
  atoms3r_ical::ImuLogRow row;
  atoms3r_ical::printImuLogRow(row);
  const auto expected = Serial.bytes;
  Serial.bytes.clear();

  // Exercise wraparound of the application ring and short/zero HWCDC writes.
  for (unsigned i = 0; i < 48; ++i) {
    assert(log.readyForSample());
    atoms3r_ical::printImuLogRow(row, log);
    log.println("command reply");
    for (unsigned slot = 0; !log.readyForSample(); ++slot) {
      assert(slot < 200);
      Serial.room = 23;
      Serial.limit = slot % 2 ? 11 : 0;
      log.service();
    }
    assert(Serial.bytes == expected + "command reply\r\n");
    Serial.bytes.clear();
  }
  log.println("old session");
  Serial.limit = 1;
  log.service();
  Serial.room = 0;
  for (int slot = 0; slot < 2000; ++slot) {
    assert(!log.readyForSample());
    log.service();
  }
  assert(Serial.bytes == "o");
  Serial.connected = false;
  log.service();
  assert(log.println("disconnected") == 0);
  Serial.connected = true;
  Serial.limit = SIZE_MAX;
  Serial.room = 256;
  assert(log.readyForSample());
  log.println("fresh");
  log.service();
  assert(Serial.bytes == "ofresh\r\n");
}

int main() {
  check<1024>();
  check<8192>();
}
