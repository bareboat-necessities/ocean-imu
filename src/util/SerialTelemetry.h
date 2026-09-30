#pragma once

/* Copyright (c) 2026 Mikhail Grushinskiy */

#include <stddef.h>

namespace ocean_imu {
namespace telemetry {

// Complete records only. The ESP32 implementation copies the record into a
// bounded queue with zero queue-wait time. It never calls Serial on the producer
// task. False means dropped, not "retry until accepted". Call from the sketch's
// task after Serial.begin(), not from an ISR. Resources are created once on the
// first call; failure disables telemetry rather than falling back to blocking IO.
bool tryEnqueueSerialRecord(const char* record, size_t length);

}  // namespace telemetry
}  // namespace ocean_imu
