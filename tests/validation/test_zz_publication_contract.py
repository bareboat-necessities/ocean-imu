"""Semantic publication contract for the main OU--III article.

The numerical/statistical contracts live in test_ou_validation.py and the
evidence-provenance tests. This file checks publication structure and claim
scope rather than pinning editorial prose word for word.
"""

import json
import re
import unittest
from pathlib import Path

import test_ou_validation as core


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"


def _flat(path: Path) -> str:
    return " ".join(path.read_text(encoding="utf-8").split())


def _macro_value(text: str, name: str) -> str:
    match = re.search(
        rf"\\providecommand\{{\\{re.escape(name)}\}}\{{([^}}]+)\}}", text
    )
    if match is None:
        raise AssertionError(f"generated macro not found: {name}")
    return match.group(1)


def _abstract_reports_committed_stationary_aggregate(self):
    with (self.RESULTS / "ou_validation_manifest.json").open(
        encoding="utf-8"
    ) as stream:
        aggregate = json.load(stream)["stationary_normalized_aggregate"]

    manuscript = (DOC / "kalman_ou-w3d.tex").read_text(encoding="utf-8")
    abstract = manuscript.split("\\begin{abstract}", 1)[1].split(
        "\\end{abstract}", 1
    )[0]
    macros = (DOC / "w3d-ou-validation-macros-generated.tex-part").read_text(
        encoding="utf-8"
    )
    expected = {
        "OUValidationOUIIINormalizedMean": aggregate["OU_III"]["mean"],
        "OUValidationOUIIINormalizedStd": aggregate["OU_III"]["std"],
        "OUValidationOUIINormalizedMean": aggregate["OU_II"]["mean"],
        "OUValidationOUIINormalizedStd": aggregate["OU_II"]["std"],
        "OUValidationNormalizedDifference": aggregate["OU_III_minus_OU_II"][
            "mean_paired_difference"
        ],
        "OUValidationNormalizedDifferenceLow": aggregate["OU_III_minus_OU_II"][
            "bootstrap_ci95_low"
        ],
        "OUValidationNormalizedDifferenceHigh": aggregate["OU_III_minus_OU_II"][
            "bootstrap_ci95_high"
        ],
    }
    for name, value in expected.items():
        self.assertIn(rf"\{name}", abstract)
        self.assertAlmostEqual(float(_macro_value(macros, name)), float(value), places=2)

    flat = " ".join(abstract.split())
    self.assertIn("paired ten-seed study", flat)
    self.assertIn("stationary JONSWAP primary endpoint", flat)
    self.assertIn("bootstrap 95\\%", flat)
    self.assertIn("current isotropic OU--III configuration", flat)
    self.assertIn("integral-regularization channel", flat)
    self.assertIn("evidence remains simulation based", flat)


def _fixed_reference_and_transition_limits_are_stated(self):
    protocol = self.read_flat("w3d-sim-charts.tex-part")
    fixed_points = self.read("w3d-ou-validation-tuning-points-generated.tex-part")
    self.assertIn("scenario-matched nondeployable reference", protocol)
    self.assertIn("not an error-minimizing oracle", protocol)
    self.assertIn("tab:ou_fixed_points", fixed_points)
    self.assertIn("FixedNominal", fixed_points)
    self.assertIn("FixedOracle", fixed_points)
    self.assertIn("kinematically closed transition", protocol)
    self.assertIn("Both endpoint spectra coexist", protocol)
    self.assertIn("need not have a unique", protocol)
    self.assertIn("pure-start, blend, and pure-end", protocol)
    self.assertIn("not a continuously evolving wave spectrum", protocol)


def _inference_is_qualified_rather_than_asserted(self):
    protocol = self.read_flat("w3d-sim-charts.tex-part")
    results = self.read_flat("w3d-baseline-comparison.tex-part")
    self.assertIn("declared primary endpoint", protocol)
    self.assertIn("sample mean and standard deviation", protocol)
    self.assertIn("percentile-bootstrap intervals", protocol)
    self.assertIn("exact paired sign-flip test", protocol)
    self.assertIn("Descriptive secondary comparisons", protocol)
    self.assertIn("Only OU--II is evaluated under the full paired inferential protocol", results)
    self.assertIn("frequency-domain double-integration", results)
    self.assertNotIn("statistically powered", protocol + results)


def _three_dimensional_and_channel_results_are_reported(self):
    protocol = self.read_flat("w3d-sim-charts.tex-part")
    results = self.read_flat("w3d-baseline-comparison.tex-part")
    figures = self.read_flat("w3d-results.tex-part")
    conclusion = self.read_flat("w3d-conclusion-summary.tex-part")

    self.assertIn("tab:ou_mc_axes", results)
    self.assertIn("lower paired mean 3-D displacement RMS", results)
    self.assertIn("all \\OUValidationThreeDScenarios", results)
    self.assertIn("paired bootstrap interval excludes zero in every", results)
    self.assertIn("par:channel-ablation", protocol)
    self.assertTrue(
        "integral-state regularization" in conclusion
        or "integral regularizer" in conclusion
    )
    for asset in (
        "w3d_ou3_jonswap_medium_xykin.pgf",
        "w3d_ou3_jonswap_medium_zkin.pgf",
        "w3d_ou3_jonswap_medium.pgf",
        "w3d_ou3_jonswap_medium_acc_bias.pgf",
        "w3d_ou3_jonswap_medium_gyro_bias.pgf",
    ):
        self.assertIn(asset, figures)


def _transition_and_secondary_ensembles_are_rescored(self):
    protocol = self.read_flat("w3d-sim-charts.tex-part")
    results = self.read_flat("w3d-baseline-comparison.tex-part")
    generated = self.read_flat("w3d-ou-validation-results-generated.tex-part")
    self.assertIn("tab:ou_transition_segments", results)
    self.assertIn("tab:ou_mc_pmstokes", generated)
    self.assertIn("OUValidationPMStokesDifference", results)
    self.assertIn("tab:ou_mc_direction", generated)
    self.assertIn("OUValidationDirectionAbsError", results)
    self.assertIn("crossing seas", protocol)


def _baseline_fairness_thresholds_and_hardware_limits_are_recorded(self):
    baseline = self.read_flat("w3d-baseline-comparison.tex-part")
    fusion = self.read("w3d-fus-methods.tex-part")
    results = self.read_flat("w3d-results.tex-part")
    self.assertIn("Parameters are frozen before evaluation", baseline)
    self.assertIn("reference displacement is never used to select gains", baseline)
    self.assertIn("same timestamped synthetic motion and inertial records", baseline)
    self.assertIn("tab:baseline-tuning-policy", baseline)
    self.assertIn("tab:implementation-gates", fusion)
    self.assertIn("functional source portability only", results)
    self.assertIn("Processor load, memory use, power, deadline margin", results)
    self.assertIn("were not measured", results)


def _contribution_is_framed_at_the_width_of_the_evidence(self):
    intro = self.read_flat("w3d-intro.tex-part")
    manuscript = self.read("kalman_ou-w3d.tex")
    title_match = re.search(r"\\title\{(.+?)\}", manuscript)
    self.assertIsNotNone(title_match)
    title = title_match.group(1)
    self.assertIn("OU-Regularized Quaternion MEKF", title)
    self.assertIn("Ornstein--Uhlenbeck process itself", intro)
    self.assertIn("local reference, not geodetic position", intro)
    self.assertIn("evidence in this paper is primarily synthetic", intro)
    self.assertIn("temperature is not swept", intro)
    self.assertIn("outside the Live-state ISS theorem", intro)



class ReorganizedPublicationContractTests(unittest.TestCase):
    def test_article_order_is_math_then_engineering_then_evidence(self):
        main = (DOC / "kalman_ou-w3d.tex").read_text(encoding="utf-8")
        ordered = [
            r"\input{w3d-state.tex-part}",
            r"\input{w3d-kalm-up.tex-part}",
            r"\input{w3d-iss-stability.tex-part}",
            r"\input{w3d-adaptation-motivation.tex-part}",
            r"\input{w3d-fus-methods.tex-part}",
            r"\input{w3d-init.tex-part}",
            r"\input{w3d-mag-hard-iron.tex-part}",
            r"\input{w3d-sim-charts.tex-part}",
            r"\input{w3d-results.tex-part}",
        ]
        positions = [main.index(item) for item in ordered]
        self.assertEqual(positions, sorted(positions))
        adaptation = (DOC / "w3d-adaptation-motivation.tex-part").read_text(
            encoding="utf-8"
        ).lstrip()
        self.assertTrue(
            adaptation.startswith(r"\input{w3d-regularization-theory.tex-part}")
        )
        pure = (DOC / "w3d-regularization-theory.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertIn(r"\input{w3d-rs-riccati-analysis.tex-part}", pure)
        self.assertIn(r"\input{w3d-rs-reduced-mse.tex-part}", pure)

    def test_intro_does_not_preempt_the_adaptation_derivation(self):
        intro = (DOC / "w3d-intro.tex-part").read_text(encoding="utf-8")
        self.assertNotIn(r"r_S\propto", intro)
        self.assertNotIn(r"\tau^{5/2}", intro)
        self.assertLess(intro.count(r"\begin{itemize}"), 2)

    def test_startup_describes_only_deployed_proxy_policy(self):
        startup = (DOC / "w3d-init.tex-part").read_text(encoding="utf-8")
        self.assertIn("measurement-only Mahony observer", startup)
        self.assertIn("starts directly in its normal", startup)
        for legacy in ("StagedMekf", "TunerWarm", "degraded warmup"):
            self.assertNotIn(legacy, startup)

    def test_methodology_contains_directional_3d_spectra(self):
        methodology = (DOC / "w3d-sim-charts.tex-part").read_text(encoding="utf-8")
        self.assertIn("spectrum_jonswap_medium_3d.pgf", methodology)
        self.assertIn("spectrum_pmstokes_medium_3d.pgf", methodology)
        self.assertIn("three-dimensional displacement", methodology)

    def test_adaptation_figure_uses_wave_band_schedule_not_carrier(self):
        draw = (
            REPO_ROOT / "plots" / "kalman_ou_iii" / "draw_plots.sh"
        ).read_text(encoding="utf-8")
        results = (DOC / "w3d-results.tex-part").read_text(encoding="utf-8")
        self.assertIn("wave_tuning_freq_hz", draw)
        self.assertIn("1.0 / (2.0 * tau_for_plot)", draw)
        self.assertIn("acceleration-band carrier", results)
        self.assertIn("intentionally omitted", results)
        self.assertNotIn("w3d-sim-charts.tex-part", draw)
        self.assertNotIn("article bibliography anchor", draw)

    def test_reduced_horizontal_analysis_keeps_scope_boundary(self):
        reduced = (DOC / "w3d-rs-reduced-mse.tex-part").read_text(encoding="utf-8")
        design = (DOC / "w3d-rs-design-interpretation.tex-part").read_text(
            encoding="utf-8"
        )
        results = (DOC / "w3d-results.tex-part").read_text(encoding="utf-8")
        self.assertIn(r"q_{\mathrm{eff}}^{1/14}", reduced)
        self.assertIn(r"q_{\mathrm{eff}}^{4/7}", reduced)
        self.assertIn("approximately $65$ and $102$", results)
        self.assertIn("not universal ocean constants", results)
        self.assertIn("isotropic", results)
        for legacy in ("1.138", "0.861", "0.36", "1.87"):
            self.assertNotIn(legacy, reduced + design + results)

    def test_future_work_preserves_major_limitations_and_companion_paths(self):
        conclusion = _flat(DOC / "w3d-conclusion-summary.tex-part")
        for phrase in (
            "IMU-to-center-of-gravity lever-arm",
            "GNSS",
            "response-amplitude operators",
            "Multimodal and multidirectional seas",
            "Wind-heel separation",
            "Temperature-dependent bias",
            "Magnetic calibration under heading excitation",
            "Independent physical validation",
        ):
            self.assertIn(phrase, conclusion)

    def test_bias_discretization_matches_the_stability_model(self):
        analytic = (DOC / "w3d-analytic-coeff.tex-part").read_text(encoding="utf-8")
        stability = (DOC / "w3d-iss-stability.tex-part").read_text(encoding="utf-8")
        self.assertIn(r"e^{-h/\tau_b}\Id", analytic)
        self.assertIn(r"1-e^{-2h/\tau_b}", analytic)
        self.assertIn(r"\phi_{b,k}=e^{-h_k/\tau_b}", stability)
        self.assertIn(r"1-e^{-2h_k/\tau_b}", stability)

    def test_stability_scope_remains_bounded_gap_and_no_reset(self):
        stability = (DOC / "w3d-iss-stability.tex-part").read_text(encoding="utf-8")
        self.assertIn("T_j:=t_j-t_{j-1}\\in[T_-,T_+]", stability)
        self.assertIn("need not be equally spaced", stability)
        self.assertIn(r"T_+\le T_{S,\max}+h_{\max}", stability)
        self.assertTrue("hard reset" in stability or "hard-reset" in stability)


