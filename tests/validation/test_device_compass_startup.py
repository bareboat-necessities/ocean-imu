#!/usr/bin/env python3
"""Run actual sketch compass code with the shipping first-sample tilt observer.

Host tests do not replace an AtomS3R bench test. They verify first-sample
availability, heading geometry, invalid-input handling, and the output wiring.
"""

import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SKETCHES = tuple(ROOT / "sensors/full_marine_ins" / name / (name + ".ino") for name in (
    "atomS3R_ins_kalman_ou2", "atomS3R_ins_kalman_ou3", "atomS3R_ins_carma3", "atomS3R_ins_tfg_carma3",
    "atomS3R_ins_carma_shared_modes", "atomS3R_ins_tfg_carma_shared_modes",
    "atomS3R_ins_tfg",
))


def function(source: str, name: str) -> str:
    """Extract a definition from the sketch, not a call or reimplementation."""
    match = re.search(
        r"^[ \t]*(?:static[ \t]+(?:inline[ \t]+)?)?(?:void|bool|float|Vector3f)[ \t]+"
        + re.escape(name) + r"\(", source, re.MULTILINE,
    )
    if match is None:
        raise AssertionError(f"Missing definition: {name}")
    brace = source.index("{", match.end())
    depth, end = 1, brace + 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    return source[match.start():end]


class DeviceCompassStartupTest(unittest.TestCase):
    def test_sampling_phase_and_overruns(self):
        # Exercise the same pacer and real timestamp wiring as every marine sketch.
        other_sketches = ("atomS3R_ins_nlo", "atomS3R_ins_pii_observer")
        for path in (*SKETCHES, *(ROOT / "sensors/full_marine_ins" / name / (name + ".ino")
                                  for name in other_sketches)):
            source = path.read_text()
            tick = function(source, "tick")
            update = function(source, "updateFilter_")
            self.assertLess(tick.index("pacer_.wait();"), tick.index("sample_us = micros();"))
            self.assertNotIn("waitForNextLoopTick_", tick)
            self.assertIn("navigation_imu_.read(", tick)
            self.assertIn("Input::pollButton();", tick)
            if path.stem in ("atomS3R_ins_kalman_ou2", "atomS3R_ins_kalman_ou3", "atomS3R_ins_tfg",
                             "atomS3R_ins_tfg_carma_shared_modes", *other_sketches):
                self.assertLess(tick.index("updateFilter_(sample"), tick.index("Input::pollButton();"))
                self.assertIn("AtomS3RNavigationImu::kReadyWindowUs", tick)
                self.assertIn("pacer_.followReadyPoll(sample_us, sample.sample_us)", tick)
                self.assertIn("navigation_imu_.configure(M5.Imu, true);", function(source, "reinitImu_"))
                self.assertIn("sample_dt_.reset();", function(source, "resetFusion_"))
                if path.stem == "atomS3R_ins_tfg_carma_shared_modes":
                    self.assertIn("timing_.host_steps += sample_dt_.usedHostTime();", update)
                    self.assertIn("timing_.missed_sensor_frames += sample_dt_.missedFrames();", update)
                elif path.stem in ("atomS3R_ins_tfg", *other_sketches):
                    self.assertIn("dt_, sample_dt_.usedHostTime(), sample_dt_.missedFrames());", tick)
                else:
                    self.assertIn("timing_.finish(sample.sample_us, got_sample, LOOP_PERIOD_US, dt_, sample_dt_.usedHostTime());", tick)
                self.assertIn("sample_dt_.advance(sample.sample_us, sample.bmi270_frame, dt_)", tick)
                self.assertNotIn("dt_ =", update)
                if path.stem == "atomS3R_ins_pii_observer":
                    for method in ("updateIMU", "updateIMUMag"):
                        self.assertRegex(update, rf"fusion_\.{method}\([^;]+dt_\);")
                else:
                    self.assertIn("fusion_.update(dt_,", update)
            else:
                self.assertIn("sample_dt_.update(s.sample_us)", update)
            self.assertIn("display_.waitForMenu();", function(source, "drawHomePending_"))
            if path.stem in other_sketches:
                self.assertIn("PollMagEveryFrame{}", tick)
            else:
                self.assertIn("mag_gate_, poll_ms", tick)

        cpp = r'''
#include <cassert>
#include <initializer_list>
#include "AtomS3R/AtomS3R_SamplePacer.h"
uint64_t clock_us;
uint64_t last_read;
uint64_t interrupt_at = UINT64_MAX;
uint32_t interrupt_us;
uint32_t micros() {
    if (clock_us == interrupt_at) clock_us += interrupt_us;
    last_read = clock_us++;
    return static_cast<uint32_t>(last_read);
}
int main() {
    constexpr uint32_t period = 5000;
    atoms3r_ical::SamplePacer pacer(period);
    const uint64_t start = UINT32_MAX-12000;
    clock_us = start;
    uint64_t last_sample = 0;
    for (unsigned i = 0; i < 1200; ++i) {
        pacer.wait();
        assert(last_read == start + i*uint64_t(period));
        // Button handling precedes the real acquisition timestamp.
        clock_us += 3;
        const uint64_t sample = clock_us;
        if (i != 61) {  // One poll returned no complete sensor frame.
            if (i) assert(sample-last_sample == (i == 62 ? 2*period : period));
            last_sample = sample;
        }
        // Variable accepted work, followed by 20 us of loopTask overhead.
        clock_us = last_read + 4300 + i%150 + 20;
    }
    const uint64_t nominal = last_read + period;
    interrupt_at = nominal-1;
    interrupt_us = 31;
    pacer.wait();
    assert(last_read == nominal+30);
    clock_us = last_read+4300;
    pacer.wait();
    assert(last_read == nominal+period); // 30 us of jitter did not accumulate.

    for (unsigned cost : {5020u, 6000u}) {
        for (int i = 0; i < 1000; ++i) {
            clock_us = last_read+cost;
            const auto ready = clock_us;
            pacer.wait();
            assert(last_read == ready); // No pad to a later slot under overload.
        }
    }

    clock_us = last_read+30000001;
    const auto resumed = clock_us;
    pacer.wait();
    assert(last_read == resumed);
    clock_us = last_read+100;
    pacer.wait();
    assert(last_read == resumed+period); // No catch-up poll after a long outage.

    const uint64_t stalled = last_read+period;
    interrupt_at = stalled-1;
    interrupt_us = 18001;
    clock_us = last_read+100;
    pacer.wait();
    assert(last_read == stalled+18000);
    clock_us = last_read+100;
    pacer.wait();
    assert(last_read == stalled+18000+period);

    // Independent host and sensor clocks. A bounded retry alone eventually
    // misses a slot as the sensor's ready edge drifts through the polling phase.
    for (bool follow : {false, true}) {
        atoms3r_ical::SamplePacer drift(period);
        clock_us = start;
        uint64_t ready_at = start, previous = start, max_gap = 0;
        unsigned empty_slots = 0;
        for (unsigned i = 0; i < 2000; ++i) {
            drift.wait();
            const uint64_t first = last_read;
            uint64_t poll = first;
            while (poll < ready_at && poll-first < 200) poll += 100;
            if (poll >= ready_at) {
                if (follow) drift.followReadyPoll(uint32_t(first), uint32_t(poll));
                if (poll-previous > max_gap) max_gap = poll-previous;
                previous = poll;
                do { ready_at += 5006; } while (ready_at <= poll);
            } else {
                ++empty_slots;
            }
            clock_us = poll + 4300;
        }
        if (follow) assert(empty_slots == 0 && max_gap == 5100);
        else assert(empty_slots == 3 && max_gap == 9800);
    }
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            path, exe = Path(tmp) / "pacer.cpp", Path(tmp) / "pacer"
            (Path(tmp) / "Arduino.h").write_text("#include <cstdint>\nuint32_t micros();\n")
            path.write_text(cpp)
            subprocess.run([*shlex.split(os.environ.get("CXX", "c++")), "-std=c++17",
                            "-Wall", "-Wextra", "-Werror", "-I", tmp, "-I", str(ROOT / "src"),
                            str(path), "-o", str(exe)], check=True)
            subprocess.run([str(exe)], check=True)

    def test_estimator_profile_counts_sensor_gaps(self):
        # A repeated frame after a gap must not count the previous gap twice.
        cpp = r'''
#include <cstdint>
#include "AtomS3R/AtomS3R_Bmi270FrameClock.h"
#include "AtomS3R/AtomS3R_LoopTiming.h"
uint32_t host;
uint32_t micros() { return host; }
Console Serial;
int main() {
    atoms3r_ical::LoopTiming timing;
    atoms3r_ical::Bmi270FrameClock clock;
    const uint32_t phases[] = {150,350,50};
    uint32_t frame = 0;
    float dt = 0;
    for (int i=0; i<4000; ++i) {
        if (i == 1001) frame += 4;
        else if (i != 1002) ++frame;
        host = frame*5000 + phases[i%3];
        const bool accepted = clock.advance(host, i<2000 ? int32_t(frame) : -1, dt);
        timing.begin();
        if (i<3000) timing.finish(host, accepted, 5000, dt, clock.usedHostTime(), clock.missedFrames());
        else timing.finish(host, accepted, 5000, dt, clock.usedHostTime());
    }
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "Arduino.h").write_text(
                "#pragma once\n#include <cstdio>\n#include <cstdint>\nuint32_t micros();\n"
                "struct Console { decltype(&std::printf) printf = &std::printf; };\nextern Console Serial;\n")
            source, exe = Path(tmp) / "profile.cpp", Path(tmp) / "profile"
            source.write_text(cpp)
            subprocess.run([*shlex.split(os.environ.get("CXX", "c++")), "-std=c++17", "-Wall",
                            "-Wextra", "-Werror", "-DOCEAN_IMU_PROFILE_LOOP=1", "-I", tmp,
                            "-I", str(ROOT / "src"), str(source), "-o", str(exe)], check=True)
            output = subprocess.check_output([str(exe)], text=True)
        self.assertEqual([line for line in output.splitlines() if line.startswith("[ESTIMATOR]")], [
            "[ESTIMATOR] dt_min_ms=5.000 dt_max_ms=5.000 host_steps=0 missed_sensor_frames=0",
            "[ESTIMATOR] dt_min_ms=5.000 dt_max_ms=20.000 host_steps=0 missed_sensor_frames=3",
            "[ESTIMATOR] dt_min_ms=4.700 dt_max_ms=5.200 host_steps=1000 missed_sensor_frames=0",
            "[ESTIMATOR] dt_min_ms=4.700 dt_max_ms=5.200 host_steps=1000",
        ])

    def test_psram_build_internal_mekf_allocation_and_fallback(self):
        cpp = r'''
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <memory>
#include <utility>
constexpr unsigned MALLOC_CAP_INTERNAL=1, MALLOC_CAP_8BIT=2;
bool allow_internal = true;
unsigned allocated = 0, freed = 0, destroyed = 0, ordinary_deletes = 0;
void* heap_caps_aligned_alloc(size_t align, size_t bytes, unsigned caps) {
    assert(caps == 3 && align == 64);
    if (!allow_internal) return nullptr;
    ++allocated;
    return std::aligned_alloc(align, bytes);
}
void heap_caps_free(void* p) { ++freed; std::free(p); }
struct alignas(64) Mekf {
    explicit Mekf(int value) { if (!value) throw value; }
    ~Mekf() { ++destroyed; }
    static void operator delete(void* p) { ++ordinary_deletes; ::operator delete(p, std::align_val_t(64)); }
};
#include "util/FilterStorage.h"
int main() {
    {
        auto owner = ocean_imu::makeFilter<Mekf>(1);
        assert(allocated == 1 && !freed && !destroyed);
        auto* first = owner.get();
        try { owner = ocean_imu::makeFilter<Mekf>(0); assert(false); } catch(int) {}
        assert(owner.get() == first && allocated == 2 && freed == 1 && !destroyed);
        allow_internal = false;
        owner = ocean_imu::makeFilter<Mekf>(1);
        assert(freed == 2 && destroyed == 1 && ordinary_deletes == 0);
    }
    assert(destroyed == 2 && ordinary_deletes == 1);
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "esp_heap_caps.h").touch()
            path, exe = Path(tmp) / "allocation.cpp", Path(tmp) / "allocation"
            path.write_text(cpp)
            subprocess.run([*shlex.split(os.environ.get("CXX", "c++")), "-std=c++17",
                             "-DARDUINO_ARCH_ESP32", "-DCONFIG_IDF_TARGET_ESP32S3", "-DBOARD_HAS_PSRAM",
                            "-I", tmp, "-I", str(ROOT / "src"), str(path), "-o", str(exe)], check=True)
            subprocess.run([str(exe)], check=True)

    def test_serial_recovery(self):
        here = Path(__file__).resolve().parent
        template = (here / "shared_tfg_serial_test.cpp.in").read_text()
        for name in ("atomS3R_ins_tfg", "atomS3R_ins_tfg_carma_shared_modes",
                     "atomS3R_ins_pii_observer", "atomS3R_ins_nlo"):
            sketch = (ROOT / "sensors/full_marine_ins" / name / f"{name}.ino").read_text()
            tfg = name in ("atomS3R_ins_tfg", "atomS3R_ins_tfg_carma_shared_modes")
            fixture = template.replace("@STREAM_SERIAL@", function(sketch, "streamSerial_"))
            fixture = fixture.replace("@WRITE_DEBUG_CHUNK@", function(sketch, "writeDebugChunk_") if tfg else "")
            with tempfile.TemporaryDirectory() as tmp:
                source, exe = Path(tmp) / "serial.cpp", Path(tmp) / "serial"
                source.write_text(fixture)
                for nmea, plotter in (((1, 0), (0, 0), (0, 1)) if tfg else ((1, 0),)):
                    with self.subTest(sketch=name, nmea=nmea, plotter=plotter):
                        subprocess.run([*shlex.split(os.environ.get("CXX", "c++")), "-std=c++17",
                                        "-Wall", "-Wextra", "-Werror", f"-DSEA_STATE_SERIAL_NMEA={nmea}",
                                        f"-DARDUINO_PLOTTER={plotter}", "-I", str(ROOT / "src"),
                                        str(source), "-o", str(exe)], check=True)
                        subprocess.run([str(exe)], check=True)

    def test_startup_compass_then_fused_live_output(self):
        for path in SKETCHES:
            with self.subTest(sketch=path.name):
                source = path.read_text()
                update = function(source, "updateFilter_")
                for consumer in ("streamSerial_", "updateUI_"):
                    self.assertIn("refreshOutput_();", function(source, consumer))
                self.assertIn("rot_.update(world_z, dt_);", update)
                self.assertIn("if (!have_output_) refreshOutput_();", update)
                update += function(source, "refreshOutput_")
                compass = function(source, "updateCompassHeading_")
                serial = function(source, "streamSerial_")
                self.assertIn("m_cal_ = runtime_.applyMag(s.m);", update)
                self.assertIn("updateCompassHeading_(", update)
                self.assertNotIn("heading_deg_ = heading_est_deg", update)
                self.assertNotIn("fusion_", compass)
                self.assertIn("heading_valid_ = heading_fused_ || heading_mag_ok_;", compass)
                self.assertIn("heading_deg_ = heading_fused_ ? ins::wrap360Deg(fused_heading_deg)", compass)
                self.assertIn("ins::magneticHeadingFromDownAndMagBody(", compass)
                self.assertIn('#include "util/MagneticHeading.h"', source)
                self.assertIn("if (heading_valid_) {\n      nmea_hdm", serial)
                self.assertIn("fusion_.isLive()", serial)
                self.assertRegex(update, r"WaveDirectionReport::from\([^;]*fusion_\.isLive\(\)\);")
                self.assertIn("runWizardFlow_(true)", function(source, "begin"))
                # A missing or cancelled calibration must not stop operation.
                self.assertIn("if (!have_blob_)", function(source, "begin"))
                self.assertNotIn("while (!have_blob_)", function(source, "begin"))
                self.assertNotIn("if (!have_blob_) return;", update)
                if "kalman_ou" in path.name or "carma3" in path.name or "carma_shared_modes" in path.name:
                    self.assertIn("fusion_.attitudeQuat()", update)
                    self.assertIn("startupProxyQuat()", update)
                    self.assertIn("startupProxyTiltQuat()", update)
                    self.assertIn("q_compass_tilt", update)
                    self.assertIn("updateCompassHeading_(q_compass_tilt, attitude_ok,", update)
                    self.assertIn("startupProxyInitialized()", update)
                    self.assertIn("const bool live = fusion_.isLive();", update)
                else:
                    self.assertIn("fusion_.startupTiltQuaternion(q_bw)", update)
                    self.assertIn("if (!fusion_.isLive()", update)
                    self.assertIn("heave_m_ = displacement_det_out_.wave_clean.z();", update)
                    self.assertIn("heave_raw_m_        = displacement_up_m_.z();", update)
                    self.assertIn("fusion_.update(dt_, w_cal_, a_cal_, 35.0f);", update)
                    self.assertIn("updateWaveDirection_(q_bw, 35.0f, dt_);", update)
                    self.assertNotIn("fusion_.update(dt_, w_cal_, a_cal_, tempC);", update)

    def test_kalman_compass_redraw_under_sustained_load(self):
        # Run the sketch publisher and display task with a deterministic RTOS/LCD.
        header = (ROOT / "src/AtomS3R/AtomS3R_MarineCompassDisplay.h").read_text()
        display = header[header.index("namespace atoms3r_ical {"):]
        for path in SKETCHES:
            with self.subTest(sketch=path.name), tempfile.TemporaryDirectory() as tmp:
                sketch = path.read_text()
                self.assertIn("display_.waitForMenu();", function(sketch, "drawHomePending_"))
                cpp = (Path(__file__).parent / "compass_display_test.cpp.in").read_text()
                cpp = cpp.replace("/* DISPLAY */", display)
                cpp = cpp.replace("/* UPDATE */", function(sketch, "updateUI_"))
                source, exe = Path(tmp) / "redraw.cpp", Path(tmp) / "redraw"
                source.write_text(cpp)
                subprocess.run([*shlex.split(os.environ.get("CXX", "c++")), "-std=c++17",
                                "-Wall", "-Wextra", "-Werror", "-I", str(ROOT / "src"),
                                str(source), "-o", str(exe)], check=True)
                subprocess.run([str(exe)], check=True)

    def test_bmi270_navigation_burst_matches_m5unified(self):
        header = (ROOT / "src/AtomS3R/AtomS3R_ImuCal.h").read_text()
        reader = header.split("class AtomS3RNavigationImu {", 1)[1].split("\n};", 1)[0]
        cpp = r'''
#include <cassert>
#include <cstdint>
#include <cstring>
#include <cmath>
#include "AtomS3R/AtomS3R_MagAcquisition.h"
#include <initializer_list>
#include "AtomS3R/AtomS3R_Bmm150Compensation.h"
#include "AtomS3R/AtomS3R_Bmi270FrameClock.h"
using namespace atoms3r_ical;
uint32_t clock_us = 10000;
uint32_t micros() { return clock_us; }
uint32_t millis() { return clock_us / 1000; }
namespace m5 {
enum class board_t { board_M5AtomS3R, other };
enum imu_t { imu_bmi270, other };
struct imu_3d_t { float x, y, z; };
struct imu_data_t { imu_3d_t accel, gyro, mag; };
struct IMU_Base {
    struct imu_convert_param_t { float accel_res, gyro_res, mag_res; };
    int16_t raw[6]{};
    uint8_t ready = 0, next_ready = 0;
    uint8_t aux[8]{0x80,1,0x10,2,0x38,3,0x40,4};
    uint8_t accel_odr = 9, gyro_odr = 9;
    uint32_t sensor_ticks = 0;
    int clock_bursts = 0;
    bool combined_ok = true;
    int status_reads = 0, data_reads = 0, combined_reads = 0, aux_reads = 0;
    bool status_ok = true, burst_ok = true;
    void getConvertParam(imu_convert_param_t* p) { *p = {8.f/32768, 2000.f/32768, 10.f*4912/32760}; }
    bool readRegister(uint8_t reg, uint8_t* dst, size_t len) {
        clock_us += 100;
        if (reg == 0x40) {
            assert(len == 3);
            dst[0] = 0xa0 | accel_odr; dst[1] = 2; dst[2] = 0xe0 | gyro_odr;
            return true;
        }
        if ((reg == 0x04 && len == 23) || (reg == 0x0c && len == 15)) {
            len -= 3;
            dst[len] = uint8_t(sensor_ticks);
            dst[len+1] = uint8_t(sensor_ticks >> 8);
            dst[len+2] = uint8_t(sensor_ticks >> 16);
            ++clock_bursts;
        }
        if (reg == 0x04) {
            if (len == 8) { ++aux_reads; std::memcpy(dst, aux, 8); return true; }
            assert(len == 20); ++combined_reads;
            if (!combined_ok) return false;
            std::memcpy(dst, aux, 8); std::memcpy(dst+8, raw, 12);
            return true;
        }
        if (reg == 0x1d) {
            assert(len == 1); ++status_reads;
            if (!status_ok) return false;
            *dst = ready; ready = next_ready; next_ready = 0;
            return true;
        }
        assert(reg == 0x0c && len == 12); ++data_reads;
        if (!burst_ok) return false;
        std::memcpy(dst, raw, len);
        return true;
    }
} dev;
struct IMU_Class {
    enum { sensor_mask_mag = 4 };
    imu_t type = imu_bmi270;
    IMU_Base* secondary = nullptr;
    int offsets[9]{}, updates = 0, data_reads = 0;
    imu_t getType() { return type; }
    IMU_Base* getImuInstancePtr(int i) { return i == 0 ? &dev : secondary; }
    int getOffsetData(int i) { return offsets[i]; }
    void getTemp(float* t) { *t = 25.f; }
    unsigned update() { ++updates; return 7; }
    imu_data_t getImuData() { ++data_reads; return {{1,2,3},{4,5,6},{7,8,9}}; }
};
}
struct {
    m5::IMU_Class Imu;
    m5::board_t board = m5::board_t::board_M5AtomS3R;
    m5::board_t getBoard() { return board; }
} M5;
constexpr unsigned ATOMS3R_IMU_MASK_ACCEL = 1, ATOMS3R_IMU_MASK_GYRO = 2, kImuMaskAccelGyro = 3;
struct Vector3f {
    float x, y, z;
    Vector3f(float x, float y, float z) : x(x), y(y), z(z) {}
};
Vector3f map_sensor_xyz_to_body_ned_(float x, float y, float z, float scale) {
    return {y*scale, x*scale, -z*scale};
}
// Leave inertial axes in M5 coordinates here to compare with its conversion.
m5::imu_3d_t map_acc_to_body_ned_(m5::imu_3d_t a) { return a; }
m5::imu_3d_t map_gyr_to_body_ned_(m5::imu_3d_t w) { return w; }
''' + function(header, "map_mag_to_body_uT_") + r'''
struct {
    bool compensated = true, have_last = false; m5::IMU_Base* bmi270 = &m5::dev;
    uint8_t last_aux[8]{}; m5::imu_3d_t last_m5{};
    Vector3f last_body{NAN,NAN,NAN}; MagAcquisition acquisition; Bmm150Trim trim;
} source;
auto& atoms3rMagSource() { return source; }
''' + function(header, "readMagBody_") + r'''
struct ImuSample {
    int32_t bmi270_frame = -1;
    uint32_t sample_us, mask;
    m5::imu_3d_t a, w;
    Vector3f m{NAN,NAN,NAN};
    float tempC;
    MagAcquisition mag_timing;
};
struct ImuTemperatureCache { float read(m5::IMU_Class&, uint32_t) { return 25.f; } };
struct Gate { uint32_t next; bool due(uint32_t ms) const { return ms >= next; } };
''' + function(header, "readImuMapped") + "\nclass AtomS3RNavigationImu {" + reader + "\n};\n" + r'''
// Reproduce the pinned M5Unified IMU_Class::getImuData 16.16 arithmetic with
// unsigned shifts/negation to define Xtensa wrap at the negative full scale.
float reference(int16_t sample, bool invert, float scale) {
    uint32_t fixed = uint32_t(sample) << 16;
    if (invert) fixed = 0u - fixed;
    return int32_t(fixed) * (scale / 65536.f);
}
int main() {
    auto& dev = m5::dev;
    AtomS3RNavigationImu reader;
    reader.configure(M5.Imu);
    m5::imu_data_t data{};
    for (int i = INT16_MIN; i <= INT16_MAX; ++i) {
        for (int j = 0; j < 6; ++j) dev.raw[j] = int16_t(i + 577*j);
        dev.ready = 0xe0;
        assert(reader.update(M5.Imu, data) == 7);
        for (int j = 0; j <= 1; ++j) {
            const auto& got = j ? data.gyro : data.accel;
            const int16_t* raw = dev.raw + 3*j;
            const float scale = (j ? 2000.f : 8.f) / 32768.f;
            assert(got.x == reference(raw[1], false, scale));
            assert(got.y == reference(raw[0], true, scale));
            assert(got.z == reference(raw[2], false, scale));
        }
    }
    for (unsigned mask = 0; mask <= 7; ++mask) {
        dev.ready = (mask&1 ? 0x80 : 0) | (mask&2 ? 0x40 : 0) | (mask&4 ? 0x20 : 0);
        assert(reader.update(M5.Imu, data) == mask);
        const int reads = dev.data_reads;
        assert(reader.update(M5.Imu, data) == 0 && dev.data_reads == reads);
    }
    // Inspect the complete transaction and compensated field on a due poll.
    source.trim.valid = true;
    source.trim.dig_z1 = source.trim.dig_xyz1 = 10000; source.trim.dig_z2 = -12500;
    ImuSample s; ImuTemperatureCache temperature; uint32_t poll_ms;
    auto poll = [&](uint32_t due) {
        dev.ready = 0xc0;
        return reader.read(M5.Imu, micros(), s, temperature, Gate{due}, poll_ms);
    };
    const uint32_t frame = micros();
    const int separate = dev.data_reads;
    assert(poll(0) && dev.combined_reads == 1 && dev.aux_reads == 0 && dev.data_reads == separate);
    assert(source.acquisition.sequence == 1 && source.acquisition.frame_us == frame);
    assert(source.acquisition.read_start_us == frame+100 && source.acquisition.read_end_us == frame+200);
    assert(s.a.x == data.accel.x && s.w.z == data.gyro.z);
    float expected[3];
    assert(bmm150Compensate(Bmm150Raw::fromData(dev.aux), source.trim, expected));
    assert(s.m.x == expected[1] && s.m.y == -expected[0] && s.m.z == expected[2]);
    dev.aux[0] ^= 1; // Only a status bit changed
    assert(poll(0) && source.acquisition.sequence == 1);
    assert(poll(UINT32_MAX) && dev.combined_reads == 2 && dev.aux_reads == 0 && dev.data_reads == separate+1);

    dev.combined_ok = false; dev.aux[0] ^= 8;
    assert(poll(0) && dev.combined_reads == 3 && dev.data_reads == separate+2 && dev.aux_reads == 1);
    assert(source.acquisition.sequence == 2);
    // Poll becomes due during an ordinary read on the millisecond boundary.
    clock_us = (millis()+1)*1000-50; dev.aux[0] ^= 8;
    assert(poll(millis()+1) && dev.aux_reads == 2 && dev.data_reads == separate+3);
    assert(source.acquisition.sequence == 3);

    // A just-late frame is accepted in this slot, even across micros() wrap.
    // These flags clear on read; both sensor halves must have arrived in this slot.
    const uint32_t start = UINT32_MAX-49;
    dev.combined_ok = true;
    clock_us = start;
    auto read_window = [&]() {
        return reader.read(M5.Imu, micros(), s, temperature, Gate{UINT32_MAX}, poll_ms,
                           AtomS3RNavigationImu::kReadyWindowUs);
    };
    dev.ready = 0x80; dev.next_ready = 0x40;
    const int before = dev.data_reads + dev.combined_reads;
    assert(read_window() && s.sample_us == uint32_t(start+100));
    assert(dev.data_reads + dev.combined_reads == before+1);
    // Repeat one delayed frame between two on-time releases, like a 200 Hz loop.
    for (uint32_t window : {0u, AtomS3RNavigationImu::kReadyWindowUs}) {
        uint32_t previous = 0;
        for (unsigned i = 0; i < 3; ++i) {
            clock_us = start + 5000*i;
            dev.ready = i == 1 ? 0 : 0xc0;
            dev.next_ready = i == 1 ? 0xc0 : 0;
            bool ok = reader.read(M5.Imu, micros(), s, temperature, Gate{UINT32_MAX}, poll_ms, window);
            assert(ok == (i != 1 || window != 0));
            if (!ok) continue;
            if (i == 1) assert(s.sample_us - previous == 5100);
            if (i == 2) assert(s.sample_us - previous == (window ? 4900u : 10000u));
            previous = s.sample_us;
        }
    }
    const int complete = dev.data_reads + dev.combined_reads;
    // A wholly missing or incomplete frame cannot leak into the next slot.
    const int statuses = dev.status_reads;
    dev.ready = 0x80;
    assert(!read_window() && dev.status_reads == statuses+3);
    clock_us += 4700;
    dev.ready = 0x40;
    assert(!read_window() && dev.data_reads + dev.combined_reads == complete);
    dev.ready = 0; dev.next_ready = 0xc0;
    const auto next_slot = micros();
    assert(read_window() && s.sample_us == next_slot+100);
    dev.status_ok = false;
    const int failed_status = dev.status_reads;
    assert(!read_window() && dev.status_reads == failed_status+1);

    // Exercise the frame clock with the shipping BMI270 reader.
    dev.status_ok = true;
    reader.configure(M5.Imu, true);
    Bmi270FrameClock clock;
    float dt = 0;
    auto sample = [&](uint32_t host, uint32_t ticks, bool aux) {
        dev.sensor_ticks = ticks & 0xffffff;
        dev.ready = 0xc0; clock_us = host;
        assert(reader.read(M5.Imu, host, s, temperature, Gate{aux ? 0u : UINT32_MAX}, poll_ms));
        return clock.advance(s.sample_us, s.bmi270_frame, dt);
    };
    auto near = [&](float expected) { assert(std::fabs(dt-expected) < 1e-7f); };
    const uint32_t host = UINT32_MAX-8000;
    // Cross both counter wraps. The sample phase within each burst varies.
    assert(sample(host, 0xffff00+5, true)); near(.005f);
    assert(sample(host+5200, 0xffff80+10, false)); near(.005f);
    assert(sample(host+9900, 2, true)); near(.005f);
    assert(!clock.usedHostTime());
    assert(!sample(host+10200, 9, false)); // a new ready flag for the same XYZ slot
    dev.combined_ok = false; // retry with an inertial + clock burst
    assert(sample(host+20000, 256+5, true)); near(.010f);
    assert(!sample(host+31000, 256+7, false)); // retained data after a host stall
    dev.combined_ok = true;
    assert(sample(host+35100, 640+9, false)); near(.015f);
    assert(s.bmi270_frame == 5 && dev.clock_bursts >= 6);
    assert(sample(host+40300, 0, false)); near(.0052f); // unexpected counter reset
    assert(clock.usedHostTime());
    assert(sample(host+45100, 128, false)); near(.005f);
    assert(!clock.usedHostTime());
    assert(sample(host+400045100, 256, false)); near(400.f); // ambiguous long hiatus

    assert(clock.usedHostTime());

    clock.reset();
    assert(!clock.usedHostTime());
    dev.gyro_odr = 8; // mixed ODR: ordinary bursts and measured host intervals
    reader.configure(M5.Imu, true);
    const int previous_clock_bursts = dev.clock_bursts;
    assert(sample(1000, 80, true)); near(.005f);
    assert(sample(6200, 140, false)); near(.0052f);
    assert(sample(10900, 260, false)); near(.0047f);
    assert(clock.usedHostTime() && s.bmi270_frame == -1 && dev.clock_bursts == previous_clock_bursts);

    dev.ready = 0xc0; dev.status_ok = false;
    const int reads = dev.data_reads;
    assert(reader.update(M5.Imu, data) == 0 && dev.data_reads == reads);
    dev.status_ok = true; dev.burst_ok = false;
    assert(reader.update(M5.Imu, data) == 0 && dev.data_reads == reads+1);
    assert(M5.Imu.updates == 0);
    auto fallback = [&]() {
        const int calls = M5.Imu.updates, data_calls = M5.Imu.data_reads;
        assert(reader.update(M5.Imu, data) == 7);
        assert(M5.Imu.updates == calls+1 && M5.Imu.data_reads == data_calls+1);
        assert(data.accel.x == 1 && data.gyro.y == 5 && data.mag.z == 9);
    };
    source.compensated = false;
    fallback(); // A source change after configure cannot reuse compensated reads.
    reader.configure(M5.Imu);
    fallback(); // Non-opt-in sketches continue to use M5Unified with raw calibration.
    reader.configure(M5.Imu, true);
    assert(!reader.frameClockEnabled()); // GYR_CONF is still 100 Hz
    fallback();

    dev.gyro_odr = 9;
    dev.burst_ok = true;
    for (bool trim_available : {true, false}) {
        // The first case models a saved calibration fitted in M5 raw units.
        // In the second case reading BMM150 factory trim failed at startup.
        source.bmi270 = trim_available ? &dev : nullptr;
        source.have_last = false;
        source.acquisition = {};
        reader.configure(M5.Imu, true);
        assert(reader.frameClockEnabled());
        const int m5_updates = M5.Imu.updates;
        const int mag_reads = dev.combined_reads;
        auto read_raw = [&](uint32_t host, uint32_t frame_id, uint8_t ready) {
            clock_us = host;
            dev.sensor_ticks = (frame_id << 7) | (frame_id % 127);
            dev.ready = ready;
            return read_window();
        };
        auto check_body = [&](float m5_x, float m5_y, float m5_z) {
            assert(s.m.x == m5_y * .1f && s.m.y == m5_x * .1f && s.m.z == -m5_z * .1f);
        };

        assert(read_raw(10000, 2, 0xc0));
        check_body(7, 8, 9); // Use M5's held field until the first AUX-ready event.
        assert(dev.combined_reads == mag_reads);
        const float mag_res = 10.f*4912/32760;
        for (int i = INT16_MIN; i <= INT16_MAX; ++i) {
            const int16_t words[] = {int16_t(i), int16_t(i+555), int16_t(i+999), 0};
            std::memcpy(dev.aux, words, sizeof(words));
            assert(read_raw(clock_us+5000, uint32_t(i-INT16_MIN+3), 0xe0));
            // BMI270_Class's raw X/Y are signed values after dropping two bits.
            // Z retains all but bit 0, including the INT16_MIN wrap on inversion.
            check_body(reference(words[0] >> 2, true, mag_res),
                       reference(words[1] >> 2, false, mag_res),
                       reference(int16_t(words[2] & 0xfffe), true, mag_res));
        }
        const auto held = s.m;
        const uint32_t seq = s.mag_timing.sequence;
        std::memset(dev.aux, 0, sizeof(dev.aux));
        assert(read_raw(clock_us+5000, 10, 0xc0)); // AUX byte change without DRDY
        assert(s.mag_timing.sequence == seq && s.m.x == held.x && s.m.y == held.y && s.m.z == held.z);
        dev.combined_ok = false;
        assert(read_raw(clock_us+5000, 11, 0xe0)); // Recover XYZ/time, hold the previous field
        assert(s.bmi270_frame == 11 && s.mag_timing.sequence == seq);
        assert(s.m.x == held.x && s.m.y == held.y && s.m.z == held.z);
        dev.combined_ok = true;
        assert(!read_raw(clock_us+5000, 12, 0x20)); // Mag-only update is held for next full frame
        assert(read_raw(clock_us+5000, 13, 0xc0));
        check_body(0, 0, 0);
        assert(s.mag_timing.sequence == seq+1);

        // Replay 1000 consecutive 200 Hz frames. Host delta cycles through
        // 5200, 4700 and 5100 us. Include a genuinely missed slot and a
        // duplicate, neither of which becomes an ordinary 5 ms update.
        clock.reset();
        const uint32_t phases[] = {150, 350, 50};
        constexpr uint32_t missing = 420;
        const uint32_t base = UINT32_MAX - 12000;
        for (uint32_t i = 0; i < 1000; ++i) {
            const uint32_t host_us = base + i*5000 + phases[i%3];
            const bool got = read_raw(host_us, i, i == missing ? 0 : (i%5 ? 0xc0 : 0xe0));
            assert(got == (i != missing));
            if (!got) continue;
            assert(s.sample_us == host_us && s.bmi270_frame == int32_t(i));
            assert(clock.advance(s.sample_us, s.bmi270_frame, dt));
            assert(dt == (i == missing+1 ? .010f : .005f) && !clock.usedHostTime());
            assert(clock.missedFrames() == (i == missing+1 ? 1u : 0u));
            if (i == 501) {
                assert(read_raw(host_us+300, i, 0xc0));
                assert(!clock.advance(s.sample_us, s.bmi270_frame, dt));
            }
        }
        assert(M5.Imu.updates == m5_updates);
    }

    source.bmi270 = &dev;
    // True fallbacks retain M5's offsets, devices and measured host intervals.
    M5.Imu.offsets[8] = 7; reader.configure(M5.Imu, true);
    assert(!reader.frameClockEnabled()); fallback();
    M5.Imu.offsets[8] = 0;
    M5.Imu.secondary = &dev; reader.configure(M5.Imu, true);
    assert(!reader.frameClockEnabled()); fallback();
    M5.Imu.secondary = nullptr;
    source.compensated = true;
    M5.board = m5::board_t::other; reader.configure(M5.Imu, true); fallback();
    M5.board = m5::board_t::board_M5AtomS3R;
    M5.Imu.type = m5::other; reader.configure(M5.Imu, true); fallback();
    M5.Imu.type = m5::imu_bmi270;
    M5.Imu.secondary = &dev; reader.configure(M5.Imu, true); fallback();
    M5.Imu.secondary = nullptr;
    M5.Imu.offsets[5] = 7; reader.configure(M5.Imu, true); fallback();
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            path, exe = Path(tmp) / "burst.cpp", Path(tmp) / "burst"
            path.write_text(cpp)
            subprocess.run([*shlex.split(os.environ.get("CXX", "c++")), "-std=c++17",
                            "-I", str(ROOT / "src"), str(path), "-o", str(exe)], check=True)
            subprocess.run([str(exe)], check=True)

    def test_navigation_aux_reads_and_wizard_cadence(self):
        header = (ROOT / "src/AtomS3R/AtomS3R_ImuCal.h").read_text()
        cache = header.split("struct ImuTemperatureCache {", 1)[1].split("\n};", 1)[0]
        gate = (ROOT / "src/util/MagneticHeading.h").read_text().split("class MagFreshGate {", 1)[1].split("\n};", 1)[0]
        cpp = r'''
#include <cassert>
#include <cmath>
#include <cstdint>
#include "AtomS3R/AtomS3R_MagAcquisition.h"
using atoms3r_ical::MagAcquisition;
using atoms3r_ical::Bmi270AuxBurst;
struct ImuSample { int a, w, m; float tempC; uint32_t mask, sample_us; MagAcquisition mag_timing; };
struct Data { int accel, gyro, mag; };
namespace m5 { using imu_data_t = Data; }
struct {
    struct {
        int reads = 0, temps = 0;
        bool temp_ok = true;
        Data getImuData() { ++reads; return {reads, 2*reads, 3*reads}; }
        void getTemp(float* t) { ++temps; if (temp_ok) *t = 25.f; }
    } Imu;
} M5;
uint32_t now_ms;
uint32_t millis() { return now_ms; }
constexpr unsigned ATOMS3R_IMU_MASK_GYRO = 2, kImuMaskAccelGyro = 3;
int map_acc_to_body_ned_(int a) { return a; }
int map_gyr_to_body_ned_(int w) { return w; }
struct { bool compensated = true; int last_body = 0; MagAcquisition acquisition; } source;
auto& atoms3rMagSource() { return source; }
int bursts = 0;
int readMagBody_(int m, uint32_t time, const Bmi270AuxBurst*) {
    ++bursts;
    source.acquisition.acquired(now_ms, time, time+1, time+2);
    return source.last_body = m;
}
''' + "struct ImuTemperatureCache {" + cache + "\n};\nclass MagFreshGate {" + gate + "\n};\n" + function(header, "readImuMapped") + r'''
int main() {
    ImuSample s{};
    ImuTemperatureCache t;
    MagFreshGate gate;
    int accepted_mag = 0;
    // Cross the micros() rollover during one second of healthy samples.
    uint32_t start = UINT32_MAX - 50000;
    for (int i=0; i<200; ++i) {
        now_ms = 100 + 5*i;
        assert(readImuMapped(M5.Imu, 3, start+5000u*i, s, false, &t, gate.due(now_ms)));
        accepted_mag += gate.update(true, now_ms);
        assert(s.a == i+1 && s.w == 2*(i+1) && s.tempC == 25.f);
        assert(s.mag_timing.sequence == unsigned(bursts));
    }
    assert(accepted_mag == 29 && bursts == 29 && M5.Imu.temps == 10);
    assert(!readImuMapped(M5.Imu, 2, 0, s) && M5.Imu.reads == 200);
    now_ms += MagAcquisition::MAX_AGE_MS+1;
    assert(readImuMapped(M5.Imu, 3, 1200000, s, false, &t, false) && bursts == 30);
    source.compensated = false;
    assert(readImuMapped(M5.Imu, 3, 1205000, s, false, &t, false) && bursts == 31);
    int temps = M5.Imu.temps;
    assert(readImuMapped(M5.Imu, 2, 1210000, s, true));
    assert(readImuMapped(M5.Imu, 3, 1215000, s));
    assert(M5.Imu.temps == temps+2 && bursts == 33);
    M5.Imu.temp_ok = false;
    assert(std::isnan(t.read(M5.Imu, 1500000)));
    M5.Imu.temp_ok = true;
    assert(t.read(M5.Imu, 1505000) == 25.f);
    const Data acquired{901, 902, 0};
    int reads = M5.Imu.reads;
    source.compensated = true;
    assert(readImuMapped(M5.Imu, 3, 1510000, s, false, &t, false, &acquired));
    assert(s.a == 901 && s.w == 902 && s.sample_us == 1510000 && s.tempC == 25.f);
    assert(M5.Imu.reads == reads && bursts == 33);
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            path, exe = Path(tmp) / "acquisition.cpp", Path(tmp) / "acquisition"
            path.write_text(cpp)
            subprocess.run([*shlex.split(os.environ.get("CXX", "c++")), "-std=c++17",
                            "-I", str(ROOT / "src"), str(path), "-o", str(exe)], check=True)
            subprocess.run([str(exe)], check=True)

    def test_first_sample_geometry_and_invalid_inputs(self):
        candidates = [Path(os.environ.get("EIGEN_INCLUDE_DIR", "/usr/include/eigen3")),
                      ROOT / "third_party/eigen", ROOT / "vendor/eigen"]
        eigen = next((p for p in candidates if (p / "Eigen/Dense").is_file()), None)
        self.assertIsNotNone(eigen, "Install Eigen or set EIGEN_INCLUDE_DIR")
        for path in SKETCHES:
            with self.subTest(sketch=path.name), tempfile.TemporaryDirectory() as tmp:
                source = path.read_text()
                cpp = r'''
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "tuner/VerticalAccelComplementary.h"
#include "util/AngleUtils.h"
#include "util/MagneticHeading.h"
#include "util/QuaternionUtils.h"
using Vector3f = Eigen::Vector3f;
namespace ins = ocean_imu::ins;
constexpr float RAD_TO_DEG = 57.29577951308232f;
constexpr float D2R = 1.0f / RAD_TO_DEG;
void require(bool ok, const char* message) {
    if (!ok) { std::cerr << message << '\n'; std::exit(1); }
}
struct Compass {
    Vector3f m_cal_ = Vector3f::Zero();
    bool mag_ok_ = false;
    bool heading_mag_ok_ = false;
    float heading_mag_deg_ = NAN;
    bool heading_valid_ = false;
    bool heading_fused_ = false;
    float heading_deg_ = NAN;
''' + function(source, "updateCompassHeading_") + r'''
};
int main() {
    constexpr float g = 9.80665f;
    const Vector3f field(20.0f, 0.0f, 43.0f);
    int cases = 0;
    float worst = 0.0f;
    for (int hdg = 0; hdg < 360; hdg += 5) {
        for (float roll : {-60.0f, -30.0f, 0.0f, 30.0f, 60.0f}) {
            for (float pitch : {-60.0f, -30.0f, 0.0f, 30.0f, 60.0f}) {
                const Eigen::Quaternionf truth =
                    Eigen::AngleAxisf(hdg * D2R, Vector3f::UnitZ()) *
                    Eigen::AngleAxisf(pitch * D2R, Vector3f::UnitY()) *
                    Eigen::AngleAxisf(roll * D2R, Vector3f::UnitX());
                const Vector3f acc = truth.conjugate() * Vector3f(0, 0, -g);
                VerticalAccelComplementary observer(0.2f, 0.02f);
                require(!observer.isInitialized(), "observer unexpectedly initialized");
                observer.update(0.005f, Vector3f::Zero(), acc, g);
                require(observer.isInitialized(), "first sample did not seed tilt");
                require(!observer.isReady(), "test accidentally waited for warmup");
                Compass c;
                c.mag_ok_ = true;
                c.m_cal_ = truth.conjugate() * field;
                c.updateCompassHeading_(observer.tiltQuaternion(), true);
                require(c.heading_valid_, "first sample magnetic heading unavailable");
                float error = std::abs(std::remainder(c.heading_deg_ - hdg, 360.0f));
                worst = std::max(worst, error);
                require(error < 0.02f, "first sample magnetic heading incorrect");
                const Eigen::Quaternionf wrong_yaw =
                    Eigen::AngleAxisf(1.9f, Vector3f::UnitZ()) * truth;
                c.updateCompassHeading_(wrong_yaw, true);
                error = std::abs(std::remainder(c.heading_deg_ - hdg, 360.0f));
                require(c.heading_valid_ && error < 0.02f, "filter yaw leaked into compass");
                ++cases;
            }
        }
    }
    Compass c;
    const Eigen::Quaternionf q = Eigen::Quaternionf::Identity();
    c.m_cal_ = field;
    c.updateCompassHeading_(q, true);
    require(!c.heading_valid_ && std::isnan(c.heading_deg_), "missing mag reported valid");
    c.mag_ok_ = true;
    c.updateCompassHeading_(q, false);
    require(!c.heading_valid_, "missing tilt reported valid");
    c.updateCompassHeading_(q, true);
    require(c.heading_valid_ && std::abs(c.heading_deg_) < 0.01f, "north is not valid zero");
    c.updateCompassHeading_(q, true, true, 137.0f);
    require(c.heading_valid_ && c.heading_fused_ && std::abs(c.heading_deg_ - 137.0f) < 0.01f,
            "ready filter yaw was replaced by direct magnetic heading");
    c.mag_ok_ = false;
    c.updateCompassHeading_(q, true, true, 138.0f);
    require(c.heading_valid_ && c.heading_fused_, "ready gyro-propagated yaw lost between mag samples");
    c.mag_ok_ = true;
    c.updateCompassHeading_(q, true, true, NAN);
    require(c.heading_valid_ && !c.heading_fused_ && std::abs(c.heading_deg_) < 0.01f,
            "nonfinite fused yaw did not fall back to measured heading");
    const Vector3f bad_fields[] = {Vector3f::Zero(), Vector3f(0, 0, 40),
                                  Vector3f(NAN, 0, 1), Vector3f(INFINITY, 0, 1)};
    for (const Vector3f& bad : bad_fields) {
        c.m_cal_ = bad;
        c.updateCompassHeading_(q, true);
        require(!c.heading_valid_ && std::isnan(c.heading_deg_), "bad mag reported valid");
    }
    float out = 123.0f;
    require(!ins::magneticHeadingFromDownAndMagBody(Vector3f::Zero(), field, out), "zero down accepted");
    require(std::isnan(out), "invalid geometry retained old heading");
    require(!ins::magneticHeadingFromDownAndMagBody(Vector3f(NAN, 0, 1), field, out), "NaN down accepted");
    require(!ins::magneticHeadingFromDownAndMagBody(Vector3f::UnitX(), field, out), "vertical bow accepted");
    std::cout << "first-sample compass: " << cases << " orientations PASS; max error="
              << worst << " deg; invalid inputs and yaw independence PASS\n";
}
'''
                cpp_path, exe_path = Path(tmp) / "compass.cpp", Path(tmp) / "compass"
                cpp_path.write_text(cpp)
                command = shlex.split(os.environ.get("CXX", "g++")) + [
                    "-std=c++17", "-O1", "-Wall", "-Wextra", "-Werror",
                    "-DEIGEN_NON_ARDUINO", "-I", str(eigen), "-I", str(ROOT / "src"),
                    str(cpp_path), "-o", str(exe_path),
                ]
                subprocess.run(command, check=True, timeout=120)
                subprocess.run([str(exe_path)], check=True, timeout=30)


if __name__ == "__main__":
    unittest.main()
