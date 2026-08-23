"""Stable semantic publication contract for the main OU--III article.

Numerical results, replay provenance, and evidence integrity are tested
elsewhere.  This module pins article structure and technical invariants while
allowing ordinary editorial changes and equivalent mathematical notation.
"""

import json
import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
RESULTS = REPO_ROOT / "reports" / "results" / "ou_validation"


def _read(name: str) -> str:
    return (DOC / name).read_text(encoding="utf-8")


def _macro_value(text: str, name: str) -> str:
    match = re.search(
        rf"\\providecommand\{{\\{re.escape(name)}\}}\{{([^}}]+)\}}", text
    )
    if match is None:
        raise AssertionError(f"generated macro not found: {name}")
    return match.group(1)


def _assert_any(test: unittest.TestCase, text: str, *alternatives: str) -> None:
    folded = text.casefold()
    test.assertTrue(
        any(term.casefold() in folded for term in alternatives),
        f"expected one of {alternatives!r}",
    )


# Compatibility entry points imported by test_ou_validation.py.
def _abstract_reports_committed_stationary_aggregate(self):
    manuscript = _read("kalman_ou-w3d.tex")
    abstract = manuscript.split(r"\begin{abstract}", 1)[1].split(
        r"\end{abstract}", 1
    )[0]
    macros = _read("w3d-ou-validation-macros-generated.tex-part")
    with (RESULTS / "ou_validation_manifest.json").open(encoding="utf-8") as stream:
        aggregate = json.load(stream)["stationary_normalized_aggregate"]

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


def _fixed_reference_and_transition_limits_are_stated(self):
    protocol = _read("w3d-sim-charts.tex-part")
    fixed_points = _read("w3d-ou-validation-tuning-points-generated.tex-part")
    roundtrip = _read("w3d-roundtrip-transition-ablation.tex-part")
    self.assertIn(r"\label{par:fixed-points}", protocol)
    self.assertIn("tab:ou_fixed_points", fixed_points)
    self.assertIn("FixedNominal", fixed_points)
    self.assertIn("FixedOracle", fixed_points)
    self.assertIn(r"\input{w3d-roundtrip-transition-scores-generated.tex-part}", roundtrip)
    self.assertIn("400--520", roundtrip)
    self.assertIn("800--920", roundtrip)
    self.assertNotIn("tab:ou_transition_segments", roundtrip)


def _inference_is_qualified_rather_than_asserted(self):
    protocol = _read("w3d-sim-charts.tex-part")
    _assert_any(self, protocol, "paired monte carlo", "paired")
    _assert_any(self, protocol, "bootstrap", "resampling")
    _assert_any(self, protocol, "sign-flip", "sign flip")
    _assert_any(self, protocol, "primary endpoint", "primary comparison")
    self.assertNotIn("statistically powered", protocol.casefold())


def _three_dimensional_and_channel_results_are_reported(self):
    protocol = _read("w3d-sim-charts.tex-part")
    results = _read("w3d-baseline-comparison.tex-part")
    generated = _read("w3d-ou-validation-results-generated.tex-part")
    self.assertIn("tab:ou_mc_axes", results)
    self.assertIn(r"\label{par:channel-ablation}", protocol)
    self.assertIn("tab:ou_mc_channels", generated)


def _transition_and_secondary_ensembles_are_rescored(self):
    results = _read("w3d-baseline-comparison.tex-part")
    generated = _read("w3d-ou-validation-results-generated.tex-part")
    roundtrip = _read("w3d-roundtrip-transition-ablation.tex-part")
    self.assertNotIn("tab:ou_transition_segments", results)
    self.assertIn("w3d-roundtrip-transition-scores-generated.tex-part", roundtrip)
    self.assertIn("rise and fall crossfade scores kept separate", roundtrip)
    self.assertIn("tab:ou_mc_pmstokes", generated)
    self.assertIn("tab:ou_mc_direction", generated)


def _contribution_is_framed_at_the_width_of_the_evidence(self):
    manuscript = _read("kalman_ou-w3d.tex")
    intro = _read("w3d-intro.tex-part")
    conclusion = _read("w3d-conclusion-summary.tex-part")
    scope = "\n".join((intro, conclusion))
    title_match = re.search(r"\\title\{(.+?)\}", manuscript)
    self.assertIsNotNone(title_match)
    title = title_match.group(1)
    self.assertIn("MEKF", title)
    _assert_any(self, title, "OU", "Ornstein")
    _assert_any(self, scope, "simulation", "synthetic")
    _assert_any(self, scope, "local reference", "geodetic", "global navigation")
    _assert_any(self, scope, "temperature", "thermal")


def _baseline_fairness_thresholds_and_hardware_limits_are_recorded(self):
    baseline = _read("w3d-baseline-comparison.tex-part")
    fusion = _read("w3d-fus-methods.tex-part")
    results = _read("w3d-results.tex-part")
    conclusion = _read("w3d-conclusion-summary.tex-part")
    self.assertIn("tab:baseline-tuning-policy", baseline)
    self.assertIn("tab:implementation-gates", fusion)
    for method in ("OU--III", "OU--II", "PII", "TVG--NLO"):
        self.assertIn(method, baseline)
    _assert_any(self, results + conclusion, "processor", "memory", "deadline")


class ReorganizedPublicationContractTests(unittest.TestCase):
    def test_article_order_follows_design_then_stability_then_evidence(self):
        main = _read("kalman_ou-w3d.tex")
        ordered = [
            r"\input{w3d-state.tex-part}",
            r"\input{w3d-kalm-up.tex-part}",
            r"\input{w3d-adaptation-motivation.tex-part}",
            r"\input{w3d-adaptation-cadence-interpretation.tex-part}",
            r"\input{w3d-reduced-mse-envelope.tex-part}",
            r"\input{w3d-reduced-mse-validation.tex-part}",
            r"\input{w3d-adaptation-deployed-law.tex-part}",
            r"\input{w3d-rs-anisotropy-design.tex-part}",
            r"\input{w3d-fus-methods.tex-part}",
            r"\input{w3d-init.tex-part}",
            r"\input{w3d-mag-hard-iron.tex-part}",
            r"\input{w3d-iss-stability.tex-part}",
            r"\input{w3d-semiglobal-stability.tex-part}",
            r"\input{w3d-sim-charts.tex-part}",
            r"\input{w3d-results.tex-part}",
            r"\input{w3d-post-results-investigations.tex-part}",
            r"\input{w3d-conclusion-summary.tex-part}",
        ]
        positions = [main.index(item) for item in ordered]
        self.assertEqual(positions, sorted(positions))
        adaptation = _read("w3d-adaptation-motivation.tex-part")
        self.assertIn(
            r"\section{Sea-State Adaptation and Integral-State Regularization}",
            adaptation,
        )
        self.assertIn(r"\label{sec:adaptation}", adaptation)
        self.assertIn(r"\input{w3d-adaptation-observables.tex-part}", adaptation)
        for retired_input in (
            r"\input{w3d-regularization-theory.tex-part}",
            r"\input{w3d-rs-scale-theorem.tex-part}",
            r"\input{w3d-jonswap-integrated-elevation-clarification.tex-part}",
            r"\input{w3d-rs-design-interpretation.tex-part}",
        ):
            self.assertNotIn(retired_input, adaptation)

    def test_adaptation_narrative_orders_observables_sigma_tau_and_regularizer(self):
        adaptation = _read("w3d-adaptation-motivation.tex-part")
        observables = _read("w3d-adaptation-observables.tex-part")
        deployed = _read("w3d-adaptation-deployed-law.tex-part")
        anisotropy = _read("w3d-rs-anisotropy-design.tex-part")
        fusion = _read("w3d-fus-methods.tex-part")
        self.assertLess(
            adaptation.index(r"\input{w3d-adaptation-observables.tex-part}"),
            adaptation.index(r"\subsection{Measured acceleration scale and OU prior amplitude}"),
        )
        self.assertLess(
            adaptation.index(r"\subsection{Measured acceleration scale and OU prior amplitude}"),
            adaptation.index(r"\subsection{Sea time scale and OU correlation time}"),
        )
        self.assertIn(r"\label{eq:sigma-noise-subtraction}", observables)
        self.assertIn(r"\label{eq:wave-band-period}", observables)
        self.assertIn(r"\label{eq:adapt-default-rs-law}", deployed)
        self.assertIn(r"\label{eq:adapt-rs-base}", deployed)
        self.assertIn(r"\label{eq:adapt-rs-axis-ratios}", anisotropy)
        self.assertIn("isotropic", anisotropy.casefold())
        self.assertIn(r"\label{eq:adapt-one-sample-staging}", fusion)
        self.assertNotIn("activated at approximately", fusion.casefold())

    def test_abstract_uses_generated_validation_values(self):
        _abstract_reports_committed_stationary_aggregate(self)

    def test_intro_does_not_preempt_the_adaptation_derivation(self):
        intro = _read("w3d-intro.tex-part")
        self.assertNotIn(r"r_S\propto", intro)
        self.assertNotIn(r"\tau^{5/2}", intro)

    def test_startup_keeps_the_deployed_proxy_architecture(self):
        startup = _read("w3d-init.tex-part")
        self.assertIn(r"\label{sec:startup-policy}", startup)
        self.assertIn(r"\label{sec:proxy-handoff}", startup)
        _assert_any(self, startup, "measurement-only", "measurement only")
        for term in ("Mahony", "MEKF", "Live", "ISS"):
            self.assertIn(term, startup)
        for legacy in ("StagedMekf", "TunerWarm", "degraded warmup"):
            self.assertNotIn(legacy, startup)

    def test_required_methodology_and_result_assets_remain_wired(self):
        methodology = _read("w3d-sim-charts.tex-part")
        results = _read("w3d-baseline-comparison.tex-part")
        generated = _read("w3d-ou-validation-results-generated.tex-part")
        figures = _read("w3d-results.tex-part")
        roundtrip = _read("w3d-roundtrip-transition-ablation.tex-part")
        for asset in (
            "spectrum_jonswap_medium_3d.pgf",
            "spectrum_pmstokes_medium_3d.pgf",
            "w3d_ou3_jonswap_medium_xykin.pgf",
            "w3d_ou3_jonswap_medium_zkin.pgf",
            "w3d_ou3_jonswap_medium.pgf",
            "w3d_ou3_jonswap_medium_acc_bias.pgf",
            "w3d_ou3_jonswap_medium_gyro_bias.pgf",
        ):
            self.assertIn(asset, methodology + figures)
        self.assertIn("tab:ou_mc_axes", results)
        self.assertIn("ou_rs_roundtrip_transition", roundtrip)
        self.assertIn("w3d-roundtrip-transition-scores-generated.tex-part", roundtrip)
        self.assertNotIn("ou_validation_transition", results + roundtrip)
        self.assertIn("tab:ou_mc_pmstokes", generated)
        self.assertIn("tab:ou_mc_direction", generated)

    def test_adaptation_figure_uses_deployed_wave_band_schedule(self):
        draw = (REPO_ROOT / "plots" / "kalman_ou_iii" / "draw_plots.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn("wave_tuning_freq_hz", draw)
        self.assertIn("1.0 / (2.0 * tau_for_plot)", draw)

    def test_adaptation_derivation_and_bias_discretization_remain_present(self):
        adaptation = _read("w3d-adaptation-motivation.tex-part")
        observables = _read("w3d-adaptation-observables.tex-part")
        deployed = _read("w3d-adaptation-deployed-law.tex-part")
        cadence = _read("w3d-adaptation-cadence-interpretation.tex-part")
        lti = _read("w3d-lti-discrete.tex-part")
        analytic = _read("w3d-analytic-coeff.tex-part")
        stability = _read("w3d-iss-stability.tex-part")
        combined_adaptation = "\n".join((adaptation, observables, cadence, deployed))
        for marker in (
            r"\label{eq:adapt-sigma-map}",
            r"\label{eq:adapt-Lambda}",
            r"\label{eq:adapt-qeff-strong}",
            r"\label{eq:adapt-G-mag}",
            r"\label{eq:adapt-G-error}",
            r"\label{eq:adapt-RS-general}",
            r"\label{eq:adapt-rs-strong-law}",
            r"\label{eq:adapt-rs-base}",
            r"\label{eq:adapt-cR-spectral-meaning}",
            r"\label{eq:adapt-default-rs-law}",
            r"\label{eq:sigma-noise-subtraction}",
            r"\label{eq:wave-band-period}",
        ):
            self.assertIn(marker, combined_adaptation)
        for marker in (
            r"\tau^{5/2}",
            r"q_{\mathrm{eff}}",
            r"\widehat\sigma_{a,B}",
            r"\sigma_{aw}",
        ):
            self.assertIn(marker, combined_adaptation)

        self.assertIn(r"\label{eq:ba-ou-phi}", lti)
        self.assertIn(r"\label{eq:ba-ou-Qd}", lti)
        self.assertIn(r"e^{-h/\tau_b}\mat I_3", lti)
        self.assertIn(r"1-e^{-2h/\tau_b}", lti)
        self.assertIn(r"\eqref{eq:ba-ou-phi}", analytic)
        self.assertIn(r"\eqref{eq:ba-ou-Qd}", analytic)
        self.assertIn(r"\label{eq:iss-bias-contraction}", stability)
        self.assertIn(r"\|\Phi_b(k,j)\|", stability)
        _assert_any(self, stability, r"\tau_{b,+}", r"\overline\tau_b")
        self.assertNotIn(r"1-e^{-2h_k/\tau_b}", stability)

    def test_reference_point_and_sim_scope_remain_explicit(self):
        protocol = _read("w3d-sim-charts.tex-part")
        _assert_any(self, protocol, "center of gravity", "CoG")
        _assert_any(self, protocol, "simulation", "synthetic")

    def test_generated_evidence_contracts_remain_in_use(self):
        _abstract_reports_committed_stationary_aggregate(self)
        _fixed_reference_and_transition_limits_are_stated(self)
        _inference_is_qualified_rather_than_asserted(self)
        _three_dimensional_and_channel_results_are_reported(self)
        _transition_and_secondary_ensembles_are_rescored(self)
        _contribution_is_framed_at_the_width_of_the_evidence(self)
        _baseline_fairness_thresholds_and_hardware_limits_are_recorded(self)


if __name__ == "__main__":
    unittest.main()
