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
    "ou3_lever_arm_calibration.svg",
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
            # The self-calibrating arm, and the two things that bound it.
            "sec:imu-lever-arm-selfcal",
            "M(\\omega,\\dot{\\omega})",
            "annihilates any $r$ parallel to the instantaneous rotation axis",
            "not an accuracy claim",
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
            "fig:imu-lever-arm-calibration",
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
            "write_calibration_plot",
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

    def test_the_unmodeled_penalty_is_real_and_grows_with_the_arm(self):
        rows = [row for row in load(SUMMARY) if row["mode"] == "unmodeled"]
        worst_disp = max(rows, key=lambda row: float(row["disp_3d_ratio_to_baseline"]))
        worst_tilt = max(rows, key=lambda row: float(row["tilt_ratio_to_baseline"]))
        self.assertGreater(float(worst_disp["disp_3d_ratio_to_baseline"]), 1.015)
        self.assertGreater(float(worst_tilt["tilt_ratio_to_baseline"]), 1.5)
        for axis in {row["axis"] for row in rows}:
            per_axis = sorted(
                (row for row in rows if row["axis"] == axis),
                key=lambda row: float(row["distance_m"]),
            )
            injected = [float(row["installed_rms_mps2"]) for row in per_axis]
            with self.subTest(axis=axis):
                # The injected term is linear in |r| by construction; the
                # scored penalty need not be monotone, the mechanism is.
                self.assertEqual(injected, sorted(injected))

    def test_the_gyro_model_recovers_most_of_the_penalty(self):
        rows = {
            (row["mode"], row["axis"], row["distance_m"]): row
            for row in load(SUMMARY)
        }
        checked = 0
        for (mode, axis, distance), row in rows.items():
            if mode != "gyro":
                continue
            reference = rows[("unmodeled", axis, distance)]
            if float(reference["disp_3d_ratio_to_baseline"]) < 1.015:
                continue  # no penalty worth recovering at this offset
            checked += 1
            with self.subTest(axis=axis, distance=distance):
                self.assertGreater(float(row["excess_removed_fraction"]), 0.5)
        self.assertGreater(checked, 0)

    def test_the_self_calibrating_arm_is_reported_on_its_calibration(self):
        """The arm is only publishable with the caveat attached to it.

        A self-calibrating lever arm that quotes a score without quoting how
        far the calibration actually got, and how far the covariance is from
        saying so, reads as a solved problem.  It is not one, and the section
        has to keep saying which.
        """
        rows = [row for row in load(SUMMARY) if row["mode"] == "estimated"]
        self.assertTrue(rows, "the study must run the estimated arm")
        for row in rows:
            with self.subTest(axis=row["axis"], distance=row["distance_m"]):
                # Every estimated case carries a calibration error and the
                # uncertainty the filter reported alongside it.
                self.assertTrue(math.isfinite(float(row["lever_estimate_err_m"])))
                self.assertTrue(math.isfinite(float(row["lever_sigma_max_m"])))

        text = FRAGMENT.read_text(encoding="utf-8")
        for token in ("recovers", "standard deviation"):
            with self.subTest(token=token):
                self.assertIn(token, text)

    def test_the_reported_uncertainty_understates_the_calibration_error(self):
        """The finding the section leads with, checked against the numbers.

        If a re-run ever made the covariance an honest accuracy statement, the
        section's warning would be wrong and would have to be rewritten -- so
        this fails rather than passing quietly.
        """
        rows = [
            row
            for row in load(SUMMARY)
            if row["mode"] == "estimated"
            and math.isclose(float(row["distance_m"]), 0.30, abs_tol=1e-9)
        ]
        self.assertTrue(rows)
        worst = max(rows, key=lambda row: float(row["lever_estimate_err_m"]))
        self.assertGreater(
            float(worst["lever_estimate_err_m"]),
            5.0 * float(worst["lever_sigma_max_m"]),
        )

    def test_the_two_channels_peak_on_opposite_sea_states(self):
        """The section says so in prose; the runs have to say it too."""
        runs = load(RUNS)

        def value(row: dict[str, str], field: str) -> float:
            if field == "tilt":
                return max(float(row["roll_rms_deg"]), float(row["pitch_rms_deg"]))
            return float(row["disp_3d_rms_m"])

        def ratio(axis: str, spectrum: str, hs: float, field: str) -> float:
            def pick(mode: str, want_axis: str, distance: float) -> dict[str, str]:
                return next(
                    row
                    for row in runs
                    if row["mode"] == mode
                    and row["axis"] == want_axis
                    and math.isclose(
                        float(row["distance_m"]), distance, abs_tol=1e-9
                    )
                    and row["spectrum"] == spectrum
                    and math.isclose(float(row["hs_m"]), hs, abs_tol=1e-9)
                )

            return value(pick("unmodeled", axis, 0.30), field) / value(
                pick("baseline", "cg", 0.0), field
            )

        for spectrum in ("JONSWAP", "PM-Stokes"):
            with self.subTest(spectrum=spectrum, channel="displacement"):
                self.assertGreater(
                    ratio("z-vertical", spectrum, 0.27, "disp"),
                    ratio("z-vertical", spectrum, 8.50, "disp"),
                )
            with self.subTest(spectrum=spectrum, channel="tilt"):
                self.assertGreater(
                    ratio("x-athwartships", spectrum, 8.50, "tilt"),
                    ratio("x-athwartships", spectrum, 0.27, "tilt"),
                )

    def test_attitude_degrades_before_displacement(self):
        """The tilt figure's claim, from the pooled summary."""
        rows = {
            (row["axis"], row["distance_m"]): row
            for row in load(SUMMARY)
            if row["mode"] == "unmodeled"
        }
        worst = rows[("x-athwartships", "0.3")]
        self.assertGreater(float(worst["tilt_ratio_to_baseline"]), 1.5)
        self.assertLess(float(worst["disp_3d_ratio_to_baseline"]), 1.05)

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
