"""Stable semantic publication contract for the main OU--III article.

This file intentionally avoids sentence-level editorial checks. Numerical
results, replay provenance, and evidence integrity are tested elsewhere. Here we
only pin structure and technical invariants that should survive ordinary prose
editing: document wiring, generated-result macros, required assets, equations,
deployed plotting logic, theorem scope markers, and removal of retired startup
modes.
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
    """Require a concept without requiring one particular sentence."""
    folded = text.casefold()
    test.assertTrue(
        any(term.casefold() in folded for term in alternatives),
        f"expected one of {alternatives!r}",
    )


# Compatibility entry points imported by test_ou_validation.py.  Keep these
# names stable, but deliberately make their contracts semantic rather than
# editorial: labels, generated tables/macros, broad scope concepts, and method
# families are protected; exact sentences are not.
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
    results = _read("w3d-baseline-comparison.tex-part")

    self.assertIn(r"\label{par:fixed-points}", protocol)
    self.assertIn("tab:ou_fixed_points", fixed_points)
    self.assertIn("FixedNominal", fixed_points)
    self.assertIn("FixedOracle", fixed_points)
    self.assertIn("tab:ou_transition_segments", results)


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
    self.assertIn("tab:ou_transition_segments", results)
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
    def test_article_order_is_math_then_engineering_then_evidence(self):
        main = _read("kalman_ou-w3d.tex")
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

        adaptation = _read("w3d-adaptation-motivation.tex-part")
        pure = _read("w3d-regularization-theory.tex-part")
        self.assertIn(r"\input{w3d-regularization-theory.tex-part}", adaptation)
        self.assertIn(r"\input{w3d-rs-riccati-analysis.tex-part}", pure)
        self.assertIn(r"\input{w3d-rs-reduced-mse.tex-part}", pure)

    def test_abstract_uses_generated_validation_values(self):
        _abstract_reports_committed_stationary_aggregate(self)

    def test_intro_does_not_preempt_the_adaptation_derivation(self):
        intro = _read("w3d-intro.tex-part")
        self.assertNotIn(r"r_S\propto", intro)
        self.assertNotIn(r"\tau^{5/2}", intro)

    def test_startup_keeps_the_deployed_proxy_architecture(self):
        startup = _read("w3d-init.tex-part")

        # These are architecture markers, not editorial phrases. The exact
        # description of handoff is free to change.
        self.assertIn(r"\label{sec:startup-policy}", startup)
        self.assertIn(r"\label{sec:proxy-handoff}", startup)
        _assert_any(self, startup, "measurement-only", "measurement only")
        self.assertIn("Mahony", startup)
        self.assertIn("MEKF", startup)
        self.assertIn("Live", startup)
        self.assertIn("ISS", startup)

        for legacy in ("StagedMekf", "TunerWarm", "degraded warmup"):
            self.assertNotIn(legacy, startup)

    def test_required_methodology_and_result_assets_remain_wired(self):
        methodology = _read("w3d-sim-charts.tex-part")
        results = _read("w3d-baseline-comparison.tex-part")
        generated = _read("w3d-ou-validation-results-generated.tex-part")
        figures = _read("w3d-results.tex-part")

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
        self.assertIn("tab:ou_transition_segments", results)
        self.assertIn("tab:ou_mc_pmstokes", generated)
        self.assertIn("tab:ou_mc_direction", generated)

    def test_adaptation_figure_uses_deployed_wave_band_schedule(self):
        draw = (
            REPO_ROOT / "plots" / "kalman_ou_iii" / "draw_plots.sh"
        ).read_text(encoding="utf-8")
        self.assertIn("wave_tuning_freq_hz", draw)
        self.assertIn("1.0 / (2.0 * tau_for_plot)", draw)

    def test_reduced_model_and_bias_discretization_equations_remain_present(self):
        reduced = _read("w3d-rs-reduced-mse.tex-part")
        analytic = _read("w3d-analytic-coeff.tex-part")
        stability = _read("w3d-iss-stability.tex-part")

        self.assertIn(r"q_{\mathrm{eff}}^{1/14}", reduced)
        self.assertIn(r"q_{\mathrm{eff}}^{4/7}", reduced)
        self.assertIn(r"e^{-h/\tau_b}\Id", analytic)
        self.assertIn(r"1-e^{-2h/\tau_b}", analytic)
        self.assertIn(r"\phi_{b,k}=e^{-h_k/\tau_b}", stability)
        self.assertIn(r"1-e^{-2h_k/\tau_b}", stability)

    def test_stability_contract_keeps_bounded_update_gap_and_reset_scope(self):
        stability = _read("w3d-iss-stability.tex-part")
        self.assertIn(r"T_j:=t_j-t_{j-1}\in[T_-,T_+]", stability)
        self.assertIn(r"T_+\le T_{S,\max}+h_{\max}", stability)
        self.assertIn("ISS", stability)
        _assert_any(self, stability, "reset", "re-lock", "reconfiguration")

    def test_scope_limitations_are_concepts_not_fixed_sentences(self):
        scope = "\n".join(
            (
                _read("w3d-intro.tex-part"),
                _read("w3d-results.tex-part"),
                _read("w3d-conclusion-summary.tex-part"),
            )
        )
        # Keep the important claim boundaries while allowing arbitrary prose.
        _assert_any(self, scope, "simulation", "synthetic")
        _assert_any(self, scope, "local reference", "geodetic", "global navigation")
        _assert_any(self, scope, "temperature", "thermal")
        _assert_any(self, scope, "processor", "memory", "deadline")

        conclusion = _read("w3d-conclusion-summary.tex-part")
        self.assertIn(r"\label{sec:future-work}", conclusion)
        self.assertIn(r"\begin{itemize}", conclusion)

    def test_retired_horizontal_tuning_constants_do_not_reappear(self):
        publication = "\n".join(
            (
                _read("w3d-rs-reduced-mse.tex-part"),
                _read("w3d-rs-design-interpretation.tex-part"),
                _read("w3d-results.tex-part"),
            )
        )
        for legacy in ("1.138", "0.861", "0.36", "1.87"):
            self.assertNotIn(legacy, publication)


if __name__ == "__main__":
    unittest.main()
