#pragma once

// Copyright 2026, Mikhail Grushinskiy
//
// Age of the last successfully acquired magnetometer reading (host-testable).
//
// The AtomS3R source repeats its last reading between new ones, as
// M5Unified's cache does. A stalled AUX stream or failed compensated reads
// would otherwise repeat one vector indefinitely and every caller would see
// it as current. Past MAX_AGE_MS without a new reading the source reports no
// magnetometer value instead. A still sensor still produces new readings
// (noise changes the low bits), so only a stalled stream reaches the limit.

#include <cstdint>

namespace atoms3r_ical {

struct MagAcquisition {
  static constexpr uint32_t MAX_AGE_MS = 200;  // 5 sample periods at 25 Hz

  bool have = false;
  uint32_t last_ms = 0;
  // Host observation metadata, NOT a BMM150 conversion timestamp. Sequence
  // identifies successful distinct cache observations, not physical DRDY.
  uint32_t sequence = 0, frame_us = 0, read_start_us = 0, read_end_us = 0;

  void acquired(uint32_t now_ms) { have = true; last_ms = now_ms; }
  void acquired(uint32_t now_ms, uint32_t frame, uint32_t start, uint32_t end) {
    acquired(now_ms);
    ++sequence;
    frame_us = frame; read_start_us = start; read_end_us = end;
  }
  bool fresh(uint32_t now_ms) const { return have && uint32_t(now_ms - last_ms) <= MAX_AGE_MS; }
};

}  // namespace atoms3r_ical
