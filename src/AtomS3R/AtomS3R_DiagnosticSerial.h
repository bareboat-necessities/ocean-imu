#pragma once

#include <Arduino.h>
#include <freertos/ringbuf.h>
#include <algorithm>

namespace atoms3r_ical {

// Diagnostic reports can be produced in bursts by loopTask and the mag-fit
// worker. The RTOS byte ring serializes their writes. Only loopTask drains it.
template<size_t Capacity>
class DiagnosticSerial : public Print {
 public:
  void begin() { queue_ = xRingbufferCreateStatic(Capacity, RINGBUF_TYPE_BYTEBUF, storage_, &state_); }

  size_t write(uint8_t byte) override { return write(&byte, 1); }

  size_t write(const uint8_t* data, size_t count) override {
    if (!count || !Serial) return 0;
    return xRingbufferSend(queue_, data, count, 0) == pdTRUE ? count : 0;
  }

  using Print::println;
  size_t println(const char* text) { return printf("%s\r\n", text); }

  bool readyForSample() const {
    // During a reader stall keep at most the current telemetry row. Reserve
    // the rest of the queue for command replies and one-time reports.
    return !item_ && xRingbufferGetCurFreeSize(queue_) == Capacity && Serial && Serial.availableForWrite() > 0;
  }

  void service() {
    if (!Serial) {
      release();
      // A byte ring can wrap once. Discard the old connection's backlog.
      for (unsigned i = 0; i < 2 && receive(Capacity); ++i) release();
      return;
    }
    const int room = Serial.availableForWrite();
    if (room <= 0 || (!item_ && !receive(static_cast<size_t>(room)))) return;
    const size_t count = std::min(remaining_, static_cast<size_t>(room));
    const size_t written = Serial.write(cursor_, count);
    cursor_ += written;
    remaining_ -= written;
    if (!remaining_) release();
  }

 private:
  bool receive(size_t limit) {
    item_ = xRingbufferReceiveUpTo(queue_, &remaining_, 0, limit);
    cursor_ = static_cast<const uint8_t*>(item_);
    return item_ != nullptr;
  }

  void release() {
    if (item_) vRingbufferReturnItem(queue_, item_);
    item_ = nullptr;
  }

  alignas(4) uint8_t storage_[Capacity];
  StaticRingbuffer_t state_;
  RingbufHandle_t queue_ = nullptr;
  void* item_ = nullptr;
  const uint8_t* cursor_ = nullptr;
  size_t remaining_ = 0;
};

}  // namespace atoms3r_ical
