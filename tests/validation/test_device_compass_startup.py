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
    def test_startup_compass_then_fused_live_output(self):
        for path in SKETCHES:
            with self.subTest(sketch=path.name):
                source = path.read_text()
                update = function(source, "updateFilter_")
                if path.stem == "atomS3R_ins_tfg_carma_shared_modes":
                    for consumer in ("streamSerial_", "updateUI_graphics_", "updateUI_text_"):
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

    def test_carma_compass_redraw_under_sustained_load(self):
        # Exercise the sketch scheduler and the real row transfer with a fake LCD.
        ui = (ROOT / "src/AtomS3R/AtomS3R_CompassUI.h").read_text()
        for path in (p for p in SKETCHES if "carma" in p.name):
            with self.subTest(sketch=path.name), tempfile.TemporaryDirectory() as tmp:
                cpp = r'''
#include <algorithm>
#include <cassert>
#include <cstdint>
#define SEA_STATE_ENABLE_WIZARD 1
uint32_t clock_us, clock_ms;
uint32_t micros() { return clock_us; }
uint32_t millis() { return clock_ms; }
namespace lgfx { using swap565_t = uint16_t; }
uint16_t pixels[128 * 128];
int next_row, frames, transfers, max_rows;
struct {
    struct {
        void pushImage(int x, int y, int w, int rows, const uint16_t* data) {
            assert(x == 0 && w == 128 && y == next_row && data == pixels + y * w);
            assert(rows > 0 && rows <= max_rows && y + rows <= 128);
            next_row = (y + rows) % 128;
            frames += next_row == 0;
            ++transfers;
        }
    } Display;
} M5;
struct CompassUI {
    static constexpr uint32_t stripBudgetUs = 1200;
    bool _ok = true;
    int _w = 128, _h = 128, _nextRow = 128;
    struct { void* getBuffer() { return pixels; } } _frame;
''' + function(ui, "pushNextStrip") + r'''
};
struct App {
    CompassUI compass_ui_;
    bool use_graphics_ = true, compass_ui_ready_ = true;
    int tap_count_ = 0;
    uint32_t last_ui_ms_ = 0;
    static constexpr uint32_t LOOP_PERIOD_US = 5000, UI_REFRESH_MS = 100;
    void updateUI_graphics_() { assert(next_row == 0); compass_ui_._nextRow = 0; }
    void updateUI_text_() { assert(false); }
''' + function(path.read_text(), "updateUI_") + r'''
};
int main() {
    for (uint32_t used : {3000u, 4500u, 5100u}) {
        App app;
        next_row = transfers = frames = 0;
        max_rows = used == 3000 ? 16 : used == 4500 ? 5 : 1;
        for (unsigned tick = 0; tick < 300; ++tick) {
            const uint32_t start = UINT32_MAX - 2500u + tick * 5200u;
            clock_us = start + used;
            clock_ms = 100 + tick * 5u;
            const int previous = transfers;
            app.tap_count_ = tick == 10;
            app.updateUI_(start);
            assert(transfers - previous <= 1);
            if (app.tap_count_) assert(transfers == previous);
        }
        assert(frames >= 2);
    }
}
'''
                source = Path(tmp) / "redraw.cpp"
                exe = Path(tmp) / "redraw"
                source.write_text(cpp)
                subprocess.run([*shlex.split(os.environ.get("CXX", "c++")), "-std=c++17",
                                str(source), "-o", str(exe)], check=True)
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
struct ImuSample { int a, w, m; float tempC; uint32_t mask, sample_us; MagAcquisition mag_timing; };
struct Data { int accel, gyro, mag; };
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
int readMagBody_(int m, uint32_t time) {
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
