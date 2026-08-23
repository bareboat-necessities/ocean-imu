"""Publication contract for the widened OU--III analytical stability chain."""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"


def _read(name: str) -> str:
    return (DOC / name).read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", text)


class OUIIIWidenedStabilityContractTests(unittest.TestCase):
    def test_main_article_wires_single_current_stability_chain(self):
        main = _read("kalman_ou-w3d.tex")
        ordered = (
            r"\input{w3d-iss-stability.tex-part}",
            r"\input{w3d-analytical-stability-widening.tex-part}",
            r"\input{w3d-stability-widening-phase-c.tex-part}",
            r"\input{w3d-stability-widening-source-path.tex-part}",
            r"\input{w3d-semiglobal-stability.tex-part}",
            r"\input{w3d-stability-widening-phase-d.tex-part}",
            r"\input{w3d-stability-widening-phase-e.tex-part}",
            r"\input{w3d-stability-widening-phase-f.tex-part}",
            r"\input{w3d-sim-charts.tex-part}",
        )
        positions = [main.index(item) for item in ordered]
        self.assertEqual(positions, sorted(positions))
        self.assertNotIn("w3d-block-local-iss", main)

    def test_semiglobal_file_wires_only_funnel_and_hybrid_components(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        prefix = (
            r"\input{w3d-finite-live-capture.tex-part}" "\n"
            r"\input{w3d-hybrid-stability.tex-part}"
        )
        self.assertTrue(proof.startswith(prefix))
        self.assertNotIn("w3d-block-local-iss", proof)
        for marker in (
            r"\ref{thm:iss-21-ues}",
            r"\ref{thm:finite-live-capture}",
            r"\ref{thm:hybrid-quotient-capture}",
            r"\ref{thm:capture-stochastic-live}",
        ):
            self.assertIn(marker, proof)

    def test_funnel_capture_has_no_mandatory_outer_radius(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        startup = _read("w3d-semiglobal-stability.tex-part")
        init = _read("w3d-init.tex-part")
        conclusion = _read("w3d-conclusion-summary.tex-part")
        for text in (capture, startup, init, conclusion):
            self.assertNotIn(r"R_C", text)
            self.assertNotIn(r"\vct R_H\prec", text)
        self.assertIn(r"\label{eq:capture-handoff-inclusion}", capture)
        self.assertIn(r"\label{eq:capture-funnel-recursion}", capture)
        self.assertIn(r"\mathcal H_{k_H}\subset\mathcal D_{k_H}", startup)
        self.assertIn("radius summary is not propagated", startup)

    def test_no_dead_certificate_names_are_referenced(self):
        compiled = "\n".join(
            _read(name)
            for name in (
                "w3d-iss-stability.tex-part",
                "w3d-analytical-stability-widening.tex-part",
                "w3d-stability-widening-phase-c.tex-part",
                "w3d-stability-widening-source-path.tex-part",
                "w3d-finite-live-capture.tex-part",
                "w3d-hybrid-stability.tex-part",
                "w3d-semiglobal-stability.tex-part",
                "w3d-stability-widening-phase-d.tex-part",
                "w3d-stability-widening-phase-e.tex-part",
                "w3d-stability-widening-phase-f.tex-part",
                "w3d-init.tex-part",
                "w3d-conclusion-summary.tex-part",
            )
        )
        for dead in (
            "thm:iss-block-local",
            "eq:iss-block-basin",
            "eq:widen-kappaN",
            "thm:widen-explicit-P",
            "widen-cylinder-iss",
            "widen-component-iss",
            "widen-direct-iss",
            "eq:capture-comparison-sequence",
            "eq:capture-comparison-envelope",
        ):
            self.assertNotIn(dead, compiled)

    def test_stability_prose_is_forward_only(self):
        compiled = "\n".join(
            _read(name)
            for name in (
                "w3d-iss-stability.tex-part",
                "w3d-analytical-stability-widening.tex-part",
                "w3d-stability-widening-phase-c.tex-part",
                "w3d-stability-widening-source-path.tex-part",
                "w3d-finite-live-capture.tex-part",
                "w3d-hybrid-stability.tex-part",
                "w3d-semiglobal-stability.tex-part",
            )
        ).casefold()
        for phrase in (
            "previous theorem",
            "old theorem",
            "improves the old",
            "previous proof",
            "former theorem",
            "we tried",
            "fallback certificate",
            "certificate fallback",
            "fallback route",
        ):
            self.assertNotIn(phrase, compiled)

    def test_stochastic_claim_remains_conditional(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        flat = _flat(capture)
        self.assertIn(r"\label{thm:capture-stochastic-live}", capture)
        self.assertIn(r"\label{eq:capture-stochastic-drift}", capture)
        self.assertIn("stochastic certificate condition", flat)
        self.assertIn("is not inferred", flat)

    def test_hybrid_and_quotient_results_are_in_conclusion(self):
        conclusion = _read("w3d-conclusion-summary.tex-part")
        flat = _flat(conclusion).casefold()
        self.assertIn(r"\ref{thm:hybrid-live-recovery}", conclusion)
        self.assertIn("yaw quotient", flat)
        self.assertIn("source-shaped handoff reachability", flat)


if __name__ == "__main__":
    unittest.main()
