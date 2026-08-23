"""Publication contract for the OU--III stochastic concentration certificate."""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"


def _read(name: str) -> str:
    return (DOC / name).read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", text)


class OUIIIStochasticConcentrationContractTests(unittest.TestCase):
    def test_capture_keeps_mean_square_drift_without_markov_tail(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        self.assertIn(r"\label{thm:capture-stochastic-live}", capture)
        self.assertIn(r"\label{eq:capture-stochastic-drift}", capture)
        self.assertIn(r"\label{eq:capture-stochastic-ms}", capture)
        self.assertNotIn(r"\label{eq:capture-stochastic-finite-hp}", capture)
        self.assertNotIn("Markov's inequality", capture)
        flat = _flat(capture)
        self.assertIn("conditional Gaussian quadratic-form concentration", flat)
        self.assertIn("martingale Bernstein/Freedman", flat)

    def test_phase_f_has_conditional_gaussian_quadratic_tail(self):
        proof = _read("w3d-stability-widening-phase-f.tex-part")
        for marker in (
            r"\label{eq:widen-gaussian-conditional}",
            r"\label{eq:widen-gaussian-tail-constants}",
            r"\label{lem:widen-gaussian-tail}",
            r"\label{eq:widen-gaussian-tail}",
            r"\label{eq:widen-gaussian-tail-rate}",
            r"\label{eq:widen-stochastic-exit}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn(r"K\ell_+e^{-t_*}", proof)
        self.assertIn("No independence between physical samples or source words is required", flat)
        self.assertNotIn(r"\frac{s_2}{w_*^2}", proof)
        self.assertNotIn(r"\frac{s_4}{w_*^4}", proof)

    def test_phase_f_has_freedman_word_excursion_bound(self):
        proof = _read("w3d-stability-widening-phase-f.tex-part")
        for marker in (
            r"\label{eq:widen-stochastic-martingale-difference}",
            r"\label{eq:widen-stochastic-martingale-bounds}",
            r"\label{eq:widen-stochastic-drift-envelope}",
            r"\label{eq:widen-stochastic-martingale-recursion}",
            r"\label{eq:widen-stochastic-freedman-variance}",
            r"\label{thm:widen-stochastic-concentration}",
            r"\label{eq:widen-stochastic-excursion-margin}",
            r"\label{eq:widen-stochastic-freedman}",
            r"\label{eq:widen-stochastic-total-tail}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn("Bernstein/Freedman inequality", flat)
        self.assertIn("does not assume that gating preserves the conditional mean", flat)
        self.assertIn("exponential Gaussian and martingale concentration", flat)
        self.assertNotIn("Markov's inequality", proof)

    def test_total_tail_combines_excursion_and_localization_events(self):
        proof = _read("w3d-stability-widening-phase-f.tex-part")
        flat = _flat(proof)
        self.assertIn(r"+N\ell_+e^{-t_*}", proof)
        self.assertIn(r"V_W:=\frac{v_W}{1-\lambda_s^2}", proof)
        self.assertIn("raw and localized trajectories coincide", flat)
        self.assertIn("finite-horizon and high-probability", flat)

    def test_conclusion_no_longer_lists_concentration_as_future_work(self):
        conclusion = _read("w3d-conclusion-summary.tex-part")
        flat = _flat(conclusion)
        self.assertNotIn("Sharper stochastic concentration", conclusion)
        self.assertIn("Conditional Gaussian quadratic-form concentration", flat)
        self.assertIn("Bernstein/Freedman", flat)
        self.assertIn(r"$\overline\Sigma$", conclusion)
        self.assertIn(r"$(b_W,v_W)$", conclusion)


if __name__ == "__main__":
    unittest.main()
