#pragma once

// Copyright 2026, Mikhail Grushinskiy

#include <stdint.h>
#include <string.h>

namespace atoms3r_ical {

// Best-effort diagnostics from the wizard's UI thread. In particular, the
// burst of block records at the end of a pose must not wait for USB TX space.
// Submit a complete line in one write only when it fits, so a skipped record
// cannot leave its prefix attached to the next record. The caller is the only
// serial writer during pose capture/completion.
template <typename SerialT>
bool tryCalLogLine(SerialT& serial, const char* line) {
  if (!line || !serial) return false;
  const size_t len = strlen(line);
  char packet[322];  // largest AccelCalProcedure log buffer, plus CRLF
  if (len + 2 > sizeof(packet)) return false;
  const int room = serial.availableForWrite();
  if (room < 0 || (size_t)room < len + 2) return false;
  memcpy(packet, line, len);
  packet[len] = '\r';
  packet[len + 1] = '\n';
  return serial.write((const uint8_t*)packet, len + 2) == len + 2;
}

} // namespace atoms3r_ical
