#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "ou3_lever_arm_study.py"
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location("ou3_lever_arm_study", TOOL)
mod = importlib.util.module_from_spec(spec)
assert spec and spec.loader
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


class LeverArmKinematicsTests(unittest.TestCase):
    def test_centripetal_term(self):
        a = mod.lever_acceleration(
            (0.0, 0.0, 2.0),
            (0.0, 0.0, 0.0),
            (0.1, 0.0, 0.0),
        )
        self.assertAlmostEqual(a[0], -0.4, places=12)
        self.assertAlmostEqual(a[1], 0.0, places=12)
        self.assertAlmostEqual(a[2], 0.0, places=12)

    def test_tangential_term(self):
        a = mod.lever_acceleration(
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 3.0),
            (0.1, 0.0, 0.0),
        )
        self.assertAlmostEqual(a[0], 0.0, places=12)
        self.assertAlmostEqual(a[1], 0.3, places=12)
        self.assertAlmostEqual(a[2], 0.0, places=12)

    def test_an_arm_along_the_rate_contributes_nothing(self):
        """Why a purely vertical offset is the mild one under yaw."""
        a = mod.lever_acceleration(
            (0.0, 0.0, 2.0),
            (0.0, 0.0, 3.0),
            (0.0, 0.0, 0.3),
        )
        self.assertAlmostEqual(max(abs(v) for v in a), 0.0, places=12)


class LeverArmInvocationTests(unittest.TestCase):
    """The study drives the simulator's own lever-arm stage, not a rewrite."""

    def test_env_carries_the_offset_model_and_band(self):
        env = mod.lever_env((0.0, 1.0, 0.0), 0.30, "gyro", 15.0)
        self.assertEqual(env["W3D_IMU_LEVER_ARM_M"], "0,0.3,0")
        self.assertEqual(env["W3D_IMU_LEVER_ARM_MODEL"], "gyro")
        self.assertEqual(env["W3D_IMU_LEVER_ARM_CUTOFF_HZ"], "15")

    def test_every_arm_maps_to_a_model_the_simulator_accepts(self):
        source = (ROOT / "src" / "util" / "W3dSimCommon.cpp").read_text(
            encoding="utf-8"
        )
        for model in mod.MODES.values():
            with self.subTest(model=model):
                self.assertIn(f'text == "{model}"', source)

    def test_lever_arm_diagnostics_are_parsed(self):
        stdout = "\n".join(
            (
                "IMU_LEVER_ARM offset_m=[0 0.3 0] norm_m=0.3 model=gyro",
                "IMU_LEVER_ARM_RESULT norm_m=0.3 model=gyro samples=240000 "
                "installed_rms_mps2=0.37 residual_rms_mps2=0.05",
                "QUALITY_GATE: PASS=1 REASON=ok",
            )
        )
        got = mod.parse_lever_arm_result(stdout)
        self.assertAlmostEqual(got["installed_rms_mps2"], 0.37)
        self.assertAlmostEqual(got["residual_rms_mps2"], 0.05)

    def test_a_baseline_run_reports_no_injected_term(self):
        got = mod.parse_lever_arm_result("QUALITY_GATE: PASS=1 REASON=ok")
        self.assertEqual(got["installed_rms_mps2"], 0.0)
        self.assertEqual(got["residual_rms_mps2"], 0.0)


class LeverArmSummaryTests(unittest.TestCase):
    def rows(self) -> list[dict[str, object]]:
        def row(mode: str, axis: str, distance: float, disp: float, tilt: float):
            return {
                "mode": mode,
                "axis": axis,
                "distance_m": distance,
                "spectrum": "JONSWAP",
                "hs_m": 1.5,
                "disp_x_rms_m": disp,
                "disp_y_rms_m": disp,
                "disp_z_rms_m": disp,
                "disp_3d_rms_m": disp,
                "roll_rms_deg": tilt,
                "pitch_rms_deg": tilt,
                "yaw_rms_deg": 1.0,
                "accel_bias_3d_rms_mps2": 0.03,
                "gyro_bias_3d_rms_radps": 1e-4,
                "installed_rms_mps2": 0.0 if mode == "baseline" else 0.30,
                "residual_rms_mps2": {
                    "baseline": 0.0,
                    "unmodeled": 0.30,
                    "gyro": 0.05,
                    "exact": 0.0,
                }[mode],
            }

        return [
            row("baseline", "cg", 0.0, 0.100, 0.20),
            row("unmodeled", "y-fore-aft", 0.30, 0.200, 0.40),
            row("gyro", "y-fore-aft", 0.30, 0.125, 0.22),
            row("exact", "y-fore-aft", 0.30, 0.100, 0.20),
        ]

    def summaries(self) -> dict[str, dict[str, object]]:
        return {s["mode"]: s for s in mod.summarize(self.rows())}

    def test_ratios_are_relative_to_the_cg_baseline(self):
        got = self.summaries()
        self.assertAlmostEqual(got["unmodeled"]["disp_3d_ratio_to_baseline"], 2.0)
        self.assertAlmostEqual(got["exact"]["disp_3d_ratio_to_baseline"], 1.0)
        self.assertAlmostEqual(got["unmodeled"]["tilt_ratio_to_baseline"], 2.0)

    def test_excess_removed_is_measured_against_the_unmodeled_arm(self):
        got = self.summaries()
        # The exact arm returns the whole excess; the gyro arm 75% of it.
        self.assertAlmostEqual(got["exact"]["excess_removed_fraction"], 1.0)
        self.assertAlmostEqual(got["gyro"]["excess_removed_fraction"], 0.75)
        self.assertTrue(
            math.isnan(float(got["unmodeled"]["excess_removed_fraction"]))
        )

    def test_the_injected_term_is_carried_into_the_summary(self):
        got = self.summaries()
        self.assertAlmostEqual(got["unmodeled"]["installed_rms_mps2"], 0.30)
        self.assertAlmostEqual(got["unmodeled"]["residual_rms_mps2"], 0.30)
        self.assertAlmostEqual(got["exact"]["residual_rms_mps2"], 0.0)


class LeverArmCutoffSweepTests(unittest.TestCase):
    def test_sweep_is_pooled_per_cutoff_and_normalized(self):
        rows = [
            {
                "cutoff_hz": 2.0,
                "spectrum": "JONSWAP",
                "hs_m": 1.5,
                "disp_3d_rms_m": 0.2,
                "installed_rms_mps2": 0.3,
                "residual_rms_mps2": 0.15,
            },
            {
                "cutoff_hz": 15.0,
                "spectrum": "JONSWAP",
                "hs_m": 1.5,
                "disp_3d_rms_m": 0.11,
                "installed_rms_mps2": 0.3,
                "residual_rms_mps2": 0.03,
            },
        ]
        got = {s["cutoff_hz"]: s for s in mod.summarize_sweep(rows, 0.1)}
        self.assertAlmostEqual(got[2.0]["disp_3d_ratio_to_baseline"], 2.0)
        self.assertAlmostEqual(got[15.0]["disp_3d_ratio_to_baseline"], 1.1)
        self.assertAlmostEqual(got[2.0]["residual_fraction"], 0.5)
        self.assertAlmostEqual(got[15.0]["residual_fraction"], 0.1)

    def test_the_swept_band_brackets_the_deployed_one(self):
        self.assertIn(mod.DEFAULT_CUTOFF_HZ, mod.SWEEP_CUTOFFS_HZ)
        self.assertLess(min(mod.SWEEP_CUTOFFS_HZ), mod.DEFAULT_CUTOFF_HZ)
        self.assertGreater(max(mod.SWEEP_CUTOFFS_HZ), mod.DEFAULT_CUTOFF_HZ)

    def test_the_tool_and_the_simulator_agree_on_the_deployed_band(self):
        header = (ROOT / "src" / "util" / "W3dSimCommon.h").read_text(
            encoding="utf-8"
        )
        self.assertIn(
            f"float derivative_cutoff_hz = {mod.DEFAULT_CUTOFF_HZ:.1f}f;", header
        )


if __name__ == "__main__":
    unittest.main()
