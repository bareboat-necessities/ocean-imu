#pragma once
// Host stand-in for the Arduino core used by the INS sketch regression test.
// The simulated clock advances only through delay()/delayMicroseconds(), so a
// sketch loop tick lasts exactly its scheduled period.
#include <cmath>
#include <cstdarg>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

namespace ins_sketch_host {
inline uint64_t& clock_us() { static uint64_t now = 1000000; return now; }
}
inline uint32_t micros() { return uint32_t(ins_sketch_host::clock_us()); }
inline uint32_t millis() { return uint32_t(ins_sketch_host::clock_us() / 1000u); }
inline void delayMicroseconds(uint32_t us) { ins_sketch_host::clock_us() += us; }
inline void delay(uint32_t ms) { ins_sketch_host::clock_us() += uint64_t(ms) * 1000u; }
#define SET_LOOP_TASK_STACK_SIZE(bytes) static_assert((bytes) > 0, "loop stack size")

// Formats and discards output; formatting still checks the printf arguments.
struct Print {
  size_t printf(const char* fmt, ...) __attribute__((format(printf, 2, 3))) {
    char buffer[512];
    va_list args;
    va_start(args, fmt);
    const int n = std::vsnprintf(buffer, sizeof(buffer), fmt, args);
    va_end(args);
    return n > 0 ? size_t(n) : 0;
  }
  template <class... T> size_t print(T...) { return 0; }
  template <class... T> size_t println(T...) { return 0; }
  size_t write(const uint8_t*, size_t n) { return n; }
  size_t write(uint8_t) { return 1; }
  int availableForWrite() { return 4096; }
  void flush() {}
};
// Bytes written with write() go to capture() when a test installs one;
// printf() output is formatted and discarded as before.
struct HardwareSerial : Print {
  using Sink = void (*)(const char* data, size_t n);
  static Sink& capture() { static Sink sink = nullptr; return sink; }
  void begin(unsigned long) {}
  void setTxBufferSize(size_t) {}
  void setTxTimeoutMs(uint32_t) {}
  size_t write(const uint8_t* data, size_t n) {
    if (capture()) capture()(reinterpret_cast<const char*>(data), n);
    return n;
  }
  size_t write(uint8_t c) { return write(&c, 1); }
  // A test queues serial commands with input().
  static std::string& input() { static std::string pending; return pending; }
  int available() { return int(input().size()); }
  int read() {
    if (input().empty()) return -1;
    const int c = static_cast<unsigned char>(input().front());
    input().erase(0, 1);
    return c;
  }
  explicit operator bool() const { return true; }
};
inline HardwareSerial Serial;
