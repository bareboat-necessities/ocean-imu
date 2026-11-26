#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>

namespace atoms3r_ical {

// One publication (or NMEA phase). Keep its unsent suffix across sampling slots.
template<size_t Capacity>
class SerialBuffer {
 public:
  bool pending() const { return sent_ < size_; }
  void reset() { size_ = sent_ = 0; }

  size_t write(const uint8_t* bytes, size_t count) {
    if (count > Capacity - size_) return 0;
    std::memcpy(data_ + size_, bytes, count);
    size_ += count;
    return count;
  }

  template<class... Args>
  void printf(const char* format, Args... args) {
    const size_t room = sizeof(data_) - size_;
    const int n = std::snprintf(data_ + size_, room, format, args...);
    if (n > 0 && static_cast<size_t>(n) < room) size_ += static_cast<size_t>(n);
  }

  template<class Port>
  void drain(Port& port) {
    if (!pending()) return;
    const int room = port.availableForWrite();
    if (room <= 0) return;
    const size_t n = std::min(size_ - sent_, static_cast<size_t>(room));
    sent_ += port.write(reinterpret_cast<const uint8_t*>(data_ + sent_), n);
    if (!pending()) reset();
  }

  // Poll even when a partial line is the only work. The HWCDC connection
  // probe also re-arms TX after a reconnect, before the free-space check.
  template<class Port>
  bool ready(Port& port, bool due) {
    if (!due && !pending()) return false;
    if (!port) {
      reset();
      return false;
    }
    if (pending()) {
      drain(port);
      return false;
    }
    return due && port.availableForWrite() > 0;
  }

 private:
  char data_[Capacity + 1];
  size_t size_ = 0, sent_ = 0;
};

}  // namespace atoms3r_ical
