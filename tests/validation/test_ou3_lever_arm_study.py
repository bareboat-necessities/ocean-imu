#!/usr/bin/env python3
from __future__ import annotations

import csv
import importlib.util
import math
import sys
import tempfile
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

    def test_exact_rewrite_recovers_cg_acceleration(self):
        header = [
            "time",
            "disp_x",
            "disp_y",
            "disp_z",
            "vel_x",
            "vel_y",
            "vel_z",
            "acc_x",
            "acc_y",
            "acc_z",
            "acc_bx",
            "acc_by",
            "acc_bz",
            "gyro_x",
            "gyro_y",
            "gyro_z",
            "roll_deg",
            "pitch_deg",
            "yaw_deg",
            "q_wb_zu_w",
            "q_wb_zu_x",
            "q_wb_zu_y",
            "q_wb_zu_z",
        ]
        rows = []
        for i, t in enumerate((0.0, 0.005, 0.010)):
            row = {name: "0" for name in header}
            row.update(
                {
                    "time": str(t),
                    "acc_bx": "1.25",
                    "acc_by": "-0.5",
                    "acc_bz": "9.7",
                    "gyro_z": str(0.1 + 0.2 * i),
                    "q_wb_zu_w": "1",
                }
            )
            rows.append(row)

        with tempfile.TemporaryDirectory() as d:
            src = Path(d) / "source.csv"
            dst = Path(d) / "exact.csv"
            with src.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=header, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            mod.rewrite_record(src, dst, (0.3, 0.0, 0.0), exact=True)
            with dst.open(newline="") as f:
                got = list(csv.DictReader(f))

        for row in got:
            self.assertTrue(
                math.isclose(float(row["acc_bx"]), 1.25, rel_tol=0.0, abs_tol=2e-8)
            )
            self.assertTrue(
                math.isclose(float(row["acc_by"]), -0.5, rel_tol=0.0, abs_tol=2e-8)
            )
            self.assertTrue(
                math.isclose(float(row["acc_bz"]), 9.7, rel_tol=0.0, abs_tol=2e-8)
            )


if __name__ == "__main__":
    unittest.main()
