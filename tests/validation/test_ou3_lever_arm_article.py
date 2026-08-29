#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "doc" / "kalman_ou_iii"
GENERATOR = ROOT / "tools" / "ou3_lever_arm_tex.py"

spec = importlib.util.spec_from_file_location("ou3_lever_arm_tex", GENERATOR)
mod = importlib.util.module_from_spec(spec)
assert spec and spec.loader
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


class LeverArmArticleContractTests(unittest.TestCase):
    def test_post_results_includes_lever_arm_study(self):
        text = (DOC / "w3d-post-results-investigations.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertIn(r"\input{w3d-imu-lever-arm-study.tex-part}", text)

    def test_section_states_rigid_body_model_and_oracle_boundary(self):
        text = (DOC / "w3d-imu-lever-arm-study.tex-part").read_text(
            encoding="utf-8"
        )
        for token in (
            r"\dot{\omega}\times r",
            r"\omega\times(\omega\times r)",
            "10;20;30",
            "unmodeled",
            "exact-model",
            "oracle bound",
            "w3d-imu-lever-arm-results.tex-part",
            "ou3_lever_arm_3d.svg",
            "ou3_lever_arm_tilt.svg",
        ):
            self.assertIn(token, text)

    def test_generated_fragment_reports_worst_direction_and_exact_bound(self):
        baseline = {
            "mode": "baseline",
            "axis": "cg",
            "distance_m": "0",
            "disp_3d_rms_m": "0.5",
            "max_tilt_rms_deg": "0.2",
            "disp_3d_ratio_to_baseline": "1",
            "tilt_ratio_to_baseline": "1",
        }
        rows = [baseline]
        for distance in (0.1, 0.2, 0.3):
            for i, axis in enumerate(("x-athwartships", "y-fore-aft", "z-vertical")):
                rows.append(
                    {
                        "mode": "unmodeled",
                        "axis": axis,
                        "distance_m": str(distance),
                        "disp_3d_ratio_to_baseline": str(1.0 + distance * (i + 1)),
                        "tilt_ratio_to_baseline": str(1.0 + distance * (3 - i)),
                    }
                )
                rows.append(
                    {
                        "mode": "exact",
                        "axis": axis,
                        "distance_m": str(distance),
                        "disp_3d_ratio_to_baseline": "1.001",
                        "tilt_ratio_to_baseline": "1.001",
                    }
                )
        text = mod.generate(rows)
        self.assertIn(r"\SI{30}{cm}", text)
        self.assertIn(r"\num{1.900}", text)
        self.assertIn(r"\num{1.001}", text)
        self.assertIn("vertical direction", text)
        self.assertIn(r"\label{tab:imu-lever-arm}", text)


if __name__ == "__main__":
    unittest.main()
