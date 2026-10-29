#!/usr/bin/env python3
"""Contract of the CARMA off-CoG study (tools/carma_lever_arm_study.py).

The study drives the simulator's installation stage and the estimators' own
corrections through environment variables; these checks keep the tool, the
simulator, the committed evidence and docs/carma-imu-lever-arm.md agreeing.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location("carma_lever_arm_study", ROOT / "tools" / "carma_lever_arm_study.py")
mod = importlib.util.module_from_spec(spec)
assert spec and spec.loader
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)

RESULTS = ROOT / "reports" / "results" / "carma_lever_arm_study"
DOC = ROOT / "docs" / "carma-imu-lever-arm.md"


class StudyToolTests(unittest.TestCase):
    def test_offsets_reach_the_simulator_through_its_frame_map(self):
        # src/ahrs/FrameConversions.h: v_ned = [v_zu.y, v_zu.x, -v_zu.z], an
        # involution, so the study's inverse is the same permutation.
        source = (ROOT / "src" / "ahrs" / "FrameConversions.h").read_text(encoding="utf-8")
        self.assertIn("return Vector3f(v.y(), v.x(), -v.z());", source)
        for v in ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0), (0.3, -0.2, 0.7)):
            zu = mod.ned_to_zu(v)
            self.assertEqual((zu[1], zu[0], -zu[2]), v)

    def test_arm_environments(self):
        def offset(env):
            return [float(v) for v in env["W3D_IMU_LEVER_ARM_M"].split(",")]

        env = mod.arm_env("wrapper", mod.AXES["z-up"], 1.0, 15.0)
        self.assertEqual(offset(env), [0.0, 0.0, 1.0])  # above the CG is z-up +1
        self.assertEqual(env["W3D_IMU_LEVER_ARM_MODEL"], "none")
        self.assertEqual(env["W3D_IMU_LEVER_ARM_FILTER"], "wrapper")
        env = mod.arm_env("exact", mod.AXES["x-forward"], 0.3, 15.0)
        self.assertEqual(offset(env), [0.0, 0.3, 0.0])  # forward is z-up +y
        self.assertNotIn("W3D_IMU_LEVER_ARM_FILTER", env)

    def test_every_arm_is_accepted_by_the_simulators(self):
        stage = (ROOT / "src" / "util" / "W3dSimCommon.cpp").read_text(encoding="utf-8")
        sim = (ROOT / "tests" / "kalman_carma" / "carma_lever_arm-sim.cpp").read_text(encoding="utf-8")
        for arm, (model, placement) in mod.ARMS.items():
            with self.subTest(arm=arm):
                self.assertIn(f'text == "{model}"', stage)
                if placement is not None:
                    self.assertIn(f'"{placement}"', sim)

    def test_the_sweep_case_is_in_the_matrix(self):
        self.assertIn(mod.SWEEP_AXIS, mod.AXES)
        self.assertIn(mod.SWEEP_DISTANCE_M, mod.DISTANCES_M)
        self.assertIn(mod.DEFAULT_CUTOFF_HZ, mod.SWEEP_CUTOFFS_HZ)


class CommittedEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = json.loads((RESULTS / "manifest.json").read_text(encoding="utf-8"))
        cls.summary = mod.ou3.read_csv(RESULTS / "lever_arm_summary.csv")
        cls.sweep = mod.ou3.read_csv(RESULTS / "lever_arm_cutoff_summary.csv")

    def test_manifest_hashes_match_the_files(self):
        self.assertEqual(self.manifest["mode"], "full")
        for name, digest in self.manifest["files"].items():
            with self.subTest(file=name):
                self.assertEqual(hashlib.sha256((RESULTS / name).read_bytes()).hexdigest(), digest)

    def test_the_full_matrix_is_present(self):
        cells = {(r["filter"], r["mode"], r["axis"], float(r["distance_m"])) for r in self.summary}
        for filt in mod.FILTERS:
            self.assertIn((filt, "baseline", "cg", 0.0), cells)
            for arm in mod.ARMS:
                for axis in mod.AXES:
                    for d in mod.DISTANCES_M:
                        self.assertIn((filt, arm, axis, d), cells)
        self.assertTrue(all(int(r["records"]) == len(mod.RECORDS) for r in self.summary))
        swept = {(r["filter"], float(r["cutoff_hz"])) for r in self.sweep}
        self.assertEqual(swept, {(f, c) for f in mod.FILTERS for c in mod.SWEEP_CUTOFFS_HZ})

    def test_the_oracle_returns_to_the_cg(self):
        # Removing exactly what was installed must leave the CG run: the
        # arithmetic of the installation and its removal, not a model, so the
        # tolerance is tight.
        for r in self.summary:
            if r["mode"] == "exact":
                with self.subTest(case=(r["filter"], r["axis"], r["distance_m"])):
                    self.assertAlmostEqual(float(r["disp_3d_ratio_to_baseline"]), 1.0, delta=0.002)
                    self.assertAlmostEqual(float(r["tilt_ratio_to_baseline"]), 1.0, delta=0.002)

    def test_the_report_is_rendered_from_the_committed_summaries(self):
        report = (RESULTS / "lever_arm_report.md").read_text(encoding="utf-8")
        rendered = mod.markdown_report(self.summary, self.sweep, float(self.manifest["scoring_window_sec"]),
                                       self.manifest["source_commit"], True)
        self.assertEqual(report, rendered)

    def test_the_document_quotes_the_committed_table(self):
        doc = DOC.read_text(encoding="utf-8")
        match = re.search(r"<!-- BEGIN carma_lever_arm_headline -->\n(.*?)<!-- END carma_lever_arm_headline -->",
                          doc, re.S)
        self.assertIsNotNone(match, "docs/carma-imu-lever-arm.md lost its generated table")
        self.assertEqual(match.group(1), mod.headline_markdown(self.summary, self.sweep))


if __name__ == "__main__":
    unittest.main()
