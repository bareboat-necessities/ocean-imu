#!/usr/bin/env python3
from __future__ import annotations

import csv
import importlib.util
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "doc" / "kalman_ou_iii"
RESULTS = ROOT / "reports" / "results" / "ou3_lever_arm_study"
RUNS = RESULTS / "lever_arm_runs.csv"
SUMMARY = RESULTS / "lever_arm_summary.csv"
CUTOFF_SUMMARY = RESULTS / "lever_arm_cutoff_summary.csv"
FRAGMENT = DOC / "w3d-imu-lever-arm-results.tex-part"
GENERATOR = ROOT / "tools" / "ou3_lever_arm_tex.py"
FIGURES = (
    "ou3_lever_arm_penalty.svg",
    "ou3_lever_arm_tilt.svg",
    "ou3_lever_arm_sea_state.svg",
    "ou3_lever_arm_mechanism.svg",
    "ou3_lever_arm_cutoff.svg",
)

spec = importlib.util.spec_from_file_location("ou3_lever_arm_tex", GENERATOR)
mod = importlib.util.module_from_spec(spec)
assert spec and spec.loader
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


def load(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


class LeverArmArticleContractTests(unittest.TestCase):
    def test_post_results_includes_lever_arm_study(self):
        text = (DOC / "w3d-post-results-investigations.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertIn(r"\input{w3d-imu-lever-arm-study.tex-part}", text)

    def test_section_states_the_model_the_stages_and_the_oracle_boundary(self):
        text = (DOC / "w3d-imu-lever-arm-study.tex-part").read_text(
            encoding="utf-8"
        )
        for token in (
            r"\dot{\omega}\times r",
            r"\omega\times(\omega\times r)",
            "10;20;30",
            "unmodeled",
            "exact-model",
            "gyro-model",
            "oracle bound",
            # The two stages are the whole point of running this in the
            # simulator instead of pre-processing records.
            "before sensor corruption",
            "immediately before fusion",
            "w3d-imu-lever-arm-results.tex-part",
            *FIGURES,
        ):
            with self.subTest(token=token):
                self.assertIn(token, text)

    def test_section_renders_every_figure(self):
        text = (DOC / "w3d-imu-lever-arm-study.tex-part").read_text(
            encoding="utf-8"
        )
        for name in FIGURES:
            stem = name.removesuffix(".svg")
            with self.subTest(name=name):
                self.assertIn(rf"\IfFileExists{{{name}}}", text)
                self.assertIn(rf"inkscapelatex=false]{{{stem}}}", text)
        for label in (
            "fig:imu-lever-arm-penalty",
            "fig:imu-lever-arm-tilt",
            "fig:imu-lever-arm-mechanism",
            "fig:imu-lever-arm-sea-state",
            "fig:imu-lever-arm-cutoff",
        ):
            with self.subTest(label=label):
                self.assertIn(rf"\label{{{label}}}", text)

    def test_figures_are_mirrored_from_generated_evidence(self):
        for name in FIGURES:
            with self.subTest(name=name):
                generated = RESULTS / name
                mirrored = DOC / name
                self.assertTrue(generated.exists(), generated)
                self.assertEqual(generated.read_bytes(), mirrored.read_bytes(), name)

    def test_generator_publishes_the_article_figures(self):
        tool = GENERATOR.read_text(encoding="utf-8")
        study = (ROOT / "tools" / "ou3_lever_arm_study.py").read_text(
            encoding="utf-8"
        )
        for token in (
            "write_ratio_plot",
            "write_sea_state_plot",
            "write_mechanism_plot",
            "write_cutoff_plot",
            *FIGURES,
            'rcParams["svg.hashsalt"]',
            'metadata={"Date": None}',
        ):
            with self.subTest(token=token):
                self.assertIn(token, study)
        self.assertIn("lever_arm_cutoff_summary.csv", tool)


class LeverArmCommittedEvidenceTests(unittest.TestCase):
    """The article fragment must be exactly what the committed study yields."""

    def test_fragment_is_regenerable_from_the_committed_summaries(self):
        expected = mod.generate(load(SUMMARY), load(CUTOFF_SUMMARY))
        self.assertEqual(FRAGMENT.read_text(encoding="utf-8"), expected)

    def test_no_row_opens_with_an_unbraced_bracket(self):
        """A row after ``\\\\`` that opens with ``[`` is read as its optional
        vertical-space argument, and LaTeX stops with "Missing number".  This
        took down the published build once; braces are the fix, and this is
        the check that keeps them there."""
        for path in (FRAGMENT, DOC / "w3d-imu-lever-arm-study.tex-part"):
            lines = path.read_text(encoding="utf-8").splitlines()
            for previous, current in zip(lines, lines[1:]):
                if not previous.rstrip().endswith(r"\\"):
                    continue
                with self.subTest(path=path.name, line=current):
                    self.assertFalse(
                        current.lstrip().startswith("["),
                        f"{path.name}: row after a line break opens with an "
                        f"unbraced bracket: {current.strip()!r}",
                    )

    def test_generator_braces_the_bracketed_unit_row(self):
        """The same guard on the generator, so a regenerated fragment cannot
        reintroduce it even before anything is committed."""
        text = mod.generate(load(SUMMARY), load(CUTOFF_SUMMARY))
        self.assertIn("{[cm]}", text)
        lines = text.splitlines()
        for previous, current in zip(lines, lines[1:]):
            if previous.rstrip().endswith(r"\\"):
                self.assertFalse(current.lstrip().startswith("["), current)

    def test_exact_model_returns_to_the_cg_baseline(self):
        """The section's central claim, checked against the numbers."""
        rows = load(SUMMARY)
        for row in rows:
            if row["mode"] != "exact":
                continue
            with self.subTest(axis=row["axis"], distance=row["distance_m"]):
                self.assertAlmostEqual(
                    float(row["disp_3d_ratio_to_baseline"]), 1.0, delta=0.01
                )
                self.assertAlmostEqual(
                    float(row["tilt_ratio_to_baseline"]), 1.0, delta=0.01
                )
                self.assertLess(float(row["residual_rms_mps2"]), 1e-3)

    def test_installed_force_scales_with_offset_and_unmodeled_force_remains(self):
        rows = [row for row in load(SUMMARY) if row["mode"] == "unmodeled"]
        for axis in {row["axis"] for row in rows}:
            scales = []
            for row in rows:
                if row["axis"] != axis:
                    continue
                installed = float(row["installed_rms_mps2"])
                self.assertAlmostEqual(float(row["residual_rms_mps2"]), installed)
                scales.append(installed / float(row["distance_m"]))
            self.assertLess(max(scales) - min(scales), 1e-5)

    def test_gyro_model_reduces_the_installed_force(self):
        for row in load(SUMMARY):
            if row["mode"] == "gyro":
                self.assertLess(float(row["residual_rms_mps2"]),
                                float(row["installed_rms_mps2"]))

    def test_the_swept_derivative_band_is_two_sided(self):
        rows = sorted(load(CUTOFF_SUMMARY), key=lambda row: float(row["cutoff_hz"]))
        ratios = [float(row["disp_3d_ratio_to_baseline"]) for row in rows]
        best = min(range(len(ratios)), key=lambda i: ratios[i])
        self.assertGreater(best, 0, "the narrowest band should not be the best")
        self.assertLess(best, len(ratios) - 1, "the widest band should not be the best")
        self.assertTrue(math.isfinite(ratios[best]))

    def test_the_deployed_band_is_inside_the_flat_basin(self):
        rows = {float(row["cutoff_hz"]): row for row in load(CUTOFF_SUMMARY)}
        deployed = rows[15.0]
        best = min(rows.values(), key=lambda row: float(row["disp_3d_ratio_to_baseline"]))
        self.assertLess(
            float(deployed["disp_3d_ratio_to_baseline"])
            - float(best["disp_3d_ratio_to_baseline"]),
            0.005,
        )


if __name__ == "__main__":
    unittest.main()
