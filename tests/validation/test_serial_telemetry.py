"""Compile the shipping ESP32 transport against adversarial Arduino/RTOS stubs.

No wall-clock or hardware timing claims: driver hooks hold the TX worker while
running the producer, and reject every attempt to touch Serial from that producer.
"""
import pathlib
import shutil
import subprocess
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
ARDUINO = r'''
#pragma once
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cmath>
uint32_t millis();
struct FakeSerial {
  int availableForWrite();
  size_t write(const uint8_t*, size_t);
};
extern FakeSerial Serial;
'''
RTOS = r'''
#pragma once
#include <cstddef>
#include <cstdint>
using BaseType_t = int;
using UBaseType_t = unsigned;
using TickType_t = uint32_t;
using QueueHandle_t = void*;
using TaskFunction_t = void (*)(void*);
constexpr int pdTRUE = 1, pdFALSE = 0, pdPASS = 1;
constexpr unsigned tskIDLE_PRIORITY = 0;
constexpr TickType_t portMAX_DELAY = UINT32_MAX;
QueueHandle_t xQueueCreate(UBaseType_t, UBaseType_t);
void vQueueDelete(QueueHandle_t);
BaseType_t xQueueSend(QueueHandle_t, const void*, TickType_t);
BaseType_t xQueueReceive(QueueHandle_t, void*, TickType_t);
BaseType_t xTaskCreatePinnedToCore(TaskFunction_t, const char*, uint32_t,
                                void*, UBaseType_t, void*, BaseType_t);
BaseType_t xPortGetCoreID();
void vTaskDelay(TickType_t);
'''
HARNESS = r'''
#include <cassert>
#include <cstring>
#include <deque>
#include <functional>
#include <iostream>
#include <string>
#include <vector>
#include "Arduino.h"
#include "freertos/FreeRTOS.h"
#include "nmea/NmeaCompass.h"
#include "nmea/NmeaXDR.h"
using ocean_imu::telemetry::tryEnqueueSerialRecord;

struct Queue {
  size_t capacity, item_size;
  std::deque<std::vector<char>> items;
};
struct Empty {};
static Queue* queue;
static TaskFunction_t worker;
static void* worker_arg;
static bool in_sensor = true, fail_queue = false, fail_task = false;
static int queue_creates = 0, queue_deletes = 0, task_creates = 0;
static int driver_calls = 0, delays = 0, capacity = 1024;
static uint32_t now_ms = 0;
static size_t write_limit = 9999;
static std::function<void()> available_hook, write_hook, delay_hook;
static std::string output;
FakeSerial Serial;
uint32_t millis() { return now_ms; }

QueueHandle_t xQueueCreate(UBaseType_t count, UBaseType_t size) {
  ++queue_creates;
  if (fail_queue) return nullptr;
  queue = new Queue{count, size, {}};
  return queue;
}
void vQueueDelete(QueueHandle_t q) {
  ++queue_deletes; delete static_cast<Queue*>(q); queue = nullptr;
}
BaseType_t xTaskCreatePinnedToCore(TaskFunction_t fn, const char*, uint32_t stack,
                                 void* arg, UBaseType_t priority, void*, BaseType_t core) {
  ++task_creates;
  assert(stack >= 4096 && priority == tskIDLE_PRIORITY + 1);
  assert(core == 0);  // sensor is on core 1; unicore also uses core 0
  if (fail_task) return pdFALSE;
  worker = fn; worker_arg = arg;
  return pdPASS;
}
BaseType_t xPortGetCoreID() { return 1; }
BaseType_t xQueueSend(QueueHandle_t handle, const void* data, TickType_t wait) {
  assert(in_sensor && wait == 0);  // never wait for the TX consumer
  auto& q = *static_cast<Queue*>(handle);
  if (q.items.size() == q.capacity) return pdFALSE;
  const auto* bytes = static_cast<const char*>(data);
  q.items.emplace_back(bytes, bytes + q.item_size);
  return pdTRUE;
}
BaseType_t xQueueReceive(QueueHandle_t handle, void* data, TickType_t wait) {
  assert(!in_sensor && wait == portMAX_DELAY);
  auto& q = *static_cast<Queue*>(handle);
  if (q.items.empty()) throw Empty{};  // stop this deterministic worker slice
  std::memcpy(data, q.items.front().data(), q.item_size);
  q.items.pop_front();
  return pdTRUE;
}
void vTaskDelay(TickType_t ticks) {
  assert(!in_sensor && ticks == 1); ++delays;
  if (delay_hook) delay_hook();
}
int FakeSerial::availableForWrite() {
  assert(!in_sensor); ++driver_calls;
  auto hook = std::move(available_hook); available_hook = {};
  if (hook) hook();
  return capacity;
}
size_t FakeSerial::write(const uint8_t* bytes, size_t size) {
  assert(!in_sensor); ++driver_calls;
  auto hook = std::move(write_hook); write_hook = {};
  if (hook) hook();
  const size_t sent = std::min(size, write_limit);
  output.append(reinterpret_cast<const char*>(bytes), sent);
  write_limit = 9999;
  return sent;
}
void drain() {
  assert(worker);
  in_sensor = false;
  try { worker(worker_arg); } catch (const Empty&) {}
  in_sensor = true;
}
void publishWhileDriverIsHeld() {
  const int before = driver_calls;
  in_sensor = true;
  int dropped = 0;
  for (int samples = 0; samples < 10000; ++samples) {
    constexpr char stress[] = "$IIHDM,10.0,M*00\r\n";
    if (!tryEnqueueSerialRecord(stress, sizeof(stress) - 1)) ++dropped;
    // Exercise both public NMEA formatting routes as well as the transport.
    nmea_hdm("II", 10.0f);
    nmea_txt_wave_direction_confidence("II", 80.0f, true);
  }
  assert(dropped > 0 && queue->items.size() == queue->capacity);
  assert(driver_calls == before);  // producer never entered the blocked driver
  now_ms += 300;                  // queued backlog should expire, not replay
  in_sensor = false;
}
std::string sentence(const char* prefix) {
  char s[128];
  std::snprintf(s, sizeof(s), "%s*%02X\r\n", prefix, nmea0183_checksum(prefix));
  return s;
}
int main(int argc, char** argv) {
  assert(argc == 2);
  const std::string mode = argv[1];
  if (mode == "queue-failure" || mode == "task-failure") {
    fail_queue = mode == "queue-failure";
    fail_task = mode == "task-failure";
    for (int i = 0; i < 1000; ++i) assert(!tryEnqueueSerialRecord("x\n", 2));
    assert(queue_creates == 1 && driver_calls == 0);
    assert(task_creates == (fail_queue ? 0 : 1));
    assert(queue_deletes == (fail_task ? 1 : 0));
  } else if (mode == "invalid") {
    char oversized[129]{};
    assert(!tryEnqueueSerialRecord(nullptr, 2));
    assert(!tryEnqueueSerialRecord("", 0));
    assert(!tryEnqueueSerialRecord(oversized, sizeof(oversized)));
    assert(queue_creates == 0 && driver_calls == 0);
  } else if (mode == "blocked-available" || mode == "blocked-write") {
    nmea_hdm("II", 12.0f);
    assert(driver_calls == 0);
    if (mode == "blocked-available") available_hook = publishWhileDriverIsHeld;
    else write_hook = publishWhileDriverIsHeld;
    drain();
    assert(output == (mode == "blocked-write" ? sentence("$IIHDM,12.0,M") : ""));
    output.clear();
    nmea_hdm("II", 42.0f); drain();
    assert(output == sentence("$IIHDM,42.0,M"));
  } else if (mode == "full-reconnect") {
    for (int i = 0; i < 10000; ++i) nmea_hdm("II", 12.0f);
    assert(driver_calls == 0 && queue->items.size() == 16);
    capacity = 0; drain();
    assert(output.empty());
    capacity = 1024;
    nmea_hdm("II", 43.0f); drain();
    assert(output == sentence("$IIHDM,43.0,M"));
  } else if (mode == "copy-format") {
    char record[] = "copy me\r\n";
    assert(tryEnqueueSerialRecord(record, sizeof(record) - 1));
    record[0] = '!';
    nmea_hdm("II", -1.0f);
    nmea_txt_wave_direction_confidence("II", 80.0f, true);
    gen_nmea0183_xdr("$BBXDR,D,%.1f,M,DRT1", 1.5f);
    drain();
    assert(output == std::string("copy me\r\n") + sentence("$IIHDM,359.0,M")
                     + sentence("$IITXT,01,01,00,WAVCONF=80")
                     + sentence("$BBXDR,D,1.5,M,DRT1"));
  } else if (mode == "short-write") {
    int stage = 0;
    delay_hook = [&] {
      in_sensor = true;
      if (++stage == 1) {
        assert(output == "$IIHD");
        write_limit = 1;  // accept only CR of the separator
        nmea_hdm("II", 13.0f);
      } else if (stage == 2) {
        assert(output == "$IIHD\r");
        nmea_hdm("II", 14.0f);
      } else {
        assert(stage == 3);
        assert(output == std::string("$IIHD\r\n") + sentence("$IIHDM,14.0,M"));
      }
      in_sensor = false;
    };
    write_limit = 5;
    nmea_hdm("II", 12.0f); drain();
    assert(stage == 3);
  } else if (mode == "max-record") {
    char record[128];
    std::memset(record, 'x', sizeof(record));
    const std::string expected(record, sizeof(record));
    assert(tryEnqueueSerialRecord(record, sizeof(record)));
    record[0] = '!'; drain();
    assert(output == expected);
  } else if (mode == "clock-wrap") {
    now_ms = UINT32_MAX - 20;
    nmea_hdm("II", 1.0f);
    now_ms += 40; drain();
    assert(output == sentence("$IIHDM,1.0,M"));
    output.clear();
    now_ms = UINT32_MAX - 20;
    nmea_hdm("II", 2.0f);
    now_ms += 300; drain();
    assert(output.empty());
  } else {
    assert(false);
  }
  if (queue) delete queue;
  std::cout << mode << " PASS\n";
}
'''


class SerialTelemetryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which("g++")
        if compiler is None:
            raise RuntimeError("g++ is required for the serial transport regression")
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        work = pathlib.Path(cls.tmp.name)
        (work / "freertos").mkdir()
        (work / "Arduino.h").write_text(ARDUINO)
        (work / "freertos/FreeRTOS.h").write_text(RTOS)
        for name in ("queue.h", "task.h"):
            (work / "freertos" / name).write_text('#include "FreeRTOS.h"\n')
        (work / "harness.cpp").write_text(HARNESS)
        cls.binaries = []
        for unicore in (0, 1):
            binary = work / f"regression-{unicore}"
            subprocess.run([
                compiler, "-std=c++17", "-Wall", "-Wextra", "-Werror", "-O2",
                "-DARDUINO_ARCH_ESP32", f"-DCONFIG_FREERTOS_UNICORE={unicore}",
                "-I", str(work), "-I", str(ROOT / "src"),
                str(work / "harness.cpp"), str(ROOT / "src/util/SerialTelemetry.cpp"),
                "-o", str(binary),
            ], check=True, capture_output=True, text=True, timeout=60)
            cls.binaries.append(binary)

    def test_transport(self):
        for binary in self.binaries:
            for mode in ("queue-failure", "task-failure", "invalid", "blocked-available",
                         "blocked-write", "full-reconnect", "copy-format", "short-write",
                         "max-record", "clock-wrap"):
                with self.subTest(core=binary.name, scenario=mode):
                    result = subprocess.run([str(binary), mode], check=True,
                                            capture_output=True, text=True, timeout=10)
                    self.assertIn("PASS", result.stdout)


if __name__ == "__main__":
    unittest.main()
