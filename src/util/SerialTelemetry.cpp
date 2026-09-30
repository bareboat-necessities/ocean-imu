/* Copyright (c) 2026 Mikhail Grushinskiy */

#if defined(ARDUINO_ARCH_ESP32)

#include "util/SerialTelemetry.h"
#include <Arduino.h>
#include <cstring>
#include "freertos/FreeRTOS.h"
#include "freertos/queue.h"
#include "freertos/task.h"

namespace ocean_imu {
namespace telemetry {
namespace {

constexpr size_t kRecordBytes = 128;
constexpr UBaseType_t kQueueRecords = 16;
constexpr uint32_t kMaxQueueAgeMs = 250;
constexpr uint32_t kWorkerStackBytes = 4096;

struct Record {
  uint32_t enqueued_ms;
  size_t length;
  char bytes[kRecordBytes];
};

bool fresh(const Record& record) {
  return static_cast<uint32_t>(millis() - record.enqueued_ms) <= kMaxQueueAgeMs;
}

void transmitTask(void* parameter) {
  const auto queue = static_cast<QueueHandle_t>(parameter);
  // A driver short write must not join a truncated sentence to the next one.
  // Keep an unfinished CRLF separator until it has actually been accepted.
  constexpr char separator[] = "\r\n";
  size_t separator_offset = 2;
  Record record{};
  for (;;) {
    if (xQueueReceive(queue, &record, portMAX_DELAY) != pdTRUE) continue;
    if (fresh(record)) {
      if (separator_offset < 2) {
        const size_t left = 2 - separator_offset;
        const int available = Serial.availableForWrite();
        if (available >= 0 && static_cast<size_t>(available) >= left && fresh(record)) {
          const size_t sent = Serial.write(
              reinterpret_cast<const uint8_t*>(separator + separator_offset), left);
          separator_offset += sent > left ? left : sent;
        }
      }
      if (separator_offset == 2 && fresh(record)) {
        // These driver calls may wait for a mutex or touch the USB peripheral.
        // They belong ONLY to this worker, never to the IMU/filter task.
        const int available = Serial.availableForWrite();
        if (available >= 0 && static_cast<size_t>(available) >= record.length && fresh(record)) {
          const size_t sent = Serial.write(
              reinterpret_cast<const uint8_t*>(record.bytes), record.length);
          if (sent > 0 && sent < record.length) separator_offset = 0;
        }
      }
    }
    // Pace bursts and let the idle task run, including with a non-reading host.
    // This is a worker delay; it does not delay the sensor task.
    vTaskDelay(1);
  }
}

QueueHandle_t createTransport() {
  QueueHandle_t queue = xQueueCreate(kQueueRecords, sizeof(Record));
  if (queue == nullptr) return nullptr;
  BaseType_t core = 0;
#if !CONFIG_FREERTOS_UNICORE
  core = xPortGetCoreID() == 0 ? 1 : 0;
#endif
  if (xTaskCreatePinnedToCore(transmitTask, "imu-serial-tx", kWorkerStackBytes,
                              queue, tskIDLE_PRIORITY + 1, nullptr, core) != pdPASS) {
    vQueueDelete(queue);
    return nullptr;
  }
  return queue;
}

}  // namespace

bool tryEnqueueSerialRecord(const char* record, size_t length) {
  if (record == nullptr || length == 0 || length > kRecordBytes) return false;
  // All shipping NMEA producers run on the sketch task. Initialization happens
  // once, before publishing a record; the worker does not re-enter this path.
  static QueueHandle_t queue = createTransport();
  if (queue == nullptr) return false;
  Record queued{};
  queued.enqueued_ms = millis();
  queued.length = length;
  std::memcpy(queued.bytes, record, length);
  return xQueueSend(queue, &queued, 0) == pdTRUE;
}

}  // namespace telemetry
}  // namespace ocean_imu

#endif  // ARDUINO_ARCH_ESP32
