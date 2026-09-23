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
    "atomS3R_ins_kalman_ou2", "atomS3R_ins_kalman_ou3", "atomS3R_ins_tfg",
))


def function(source: str, name: str) -> str:
    """Extract a definition from the sketch, not a call or reimplementation."""
    match = re.search(
        r"^[ \t]*(?:static[ \t]+inline[ \t]+)?(?:void|bool|float|Vector3f)[ \t]+"
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
                self.assertIn("if (!have_blob_)", function(source, "begin"))
                if "kalman_ou" in path.name:
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
