"""Publication contract for the large-angle/funnel OU--III stability certificate."""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"


def _read(name: str) -> str:
    return (DOC / name).read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", text)


class OUIIILargeAngleFunnelContractTests(unittest.TestCase):
    def test_phase_c_uses_exact_group_energy_change(self):
        proof = _read("w3d-stability-widening-phase-c.tex-part")
        for marker in (
            r"\label{eq:widen-SO3-energy}",
            r"\label{eq:widen-exact-source-correction}",
            r"\label{eq:widen-exact-attitude-update}",
            r"\label{eq:widen-exact-energy-change}",
            r"\label{lem:widen-exact-group-identity}",
            r"\label{eq:widen-large-angle-sector}",
            r"\label{thm:widen-large-angle-dissipation}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn("exact for arbitrary finite", flat)
        self.assertIn("full $S\\to\\theta$", proof)
        self.assertIn("not by a small-angle Taylor radius", flat)
        self.assertIn(r"\theta_{*,i}<\pi", proof)

    def test_phase_c_has_no_radius_certificate_fallback(self):
        proof = _read("w3d-stability-widening-phase-c.tex-part")
        for dead in (
            "widen-cylinder-iss",
            "widen-cylinder-load",
            "widen-component-iss",
            "widen-direct-iss",
            "w3d-stability-widening-structured.tex-part",
        ):
            self.assertNotIn(dead, proof)
        self.assertNotIn("fallback", proof.casefold())

    def test_source_path_is_single_path_metric_route(self):
        proof = _read("w3d-stability-widening-source-path.tex-part")
        for marker in (
            r"\label{eq:widen-src-graph}",
            r"\label{eq:widen-src-word-family}",
            r"\label{eq:widen-path-lmi}",
            r"\label{thm:widen-path-ues}",
            r"\label{eq:widen-path-group-metric}",
            r"\label{eq:widen-path-direct-margin}",
            r"\label{thm:widen-path-direct-iss}",
            r"\label{eq:widen-path-contract}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn("sole quantitative normal-Live certificate path", flat)
        self.assertNotIn("fallback", flat)
        self.assertNotIn(r"\kappa_{N,m}^{\rm src}", proof)
        self.assertNotIn("special choice", flat)

    def test_capture_uses_funnel_levels_not_recursive_set_boxes(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        for marker in (
            r"\label{eq:capture-funnel-metric}",
            r"\label{eq:capture-funnel-section}",
            r"\label{eq:capture-funnel-word}",
            r"\label{eq:capture-funnel-recursion}",
            r"\label{eq:capture-inner-funnel}",
            r"\label{thm:finite-live-capture}",
            r"\label{eq:capture-certificate-pass}",
        ):
            self.assertIn(marker, capture)
        flat = _flat(capture).casefold()
        self.assertIn("propagates verified lyapunov funnel levels", flat)
        self.assertIn("positive one-sample lyapunov increments are allowed", flat)
        self.assertIn("are not propagated and are not certificate state variables", flat)
        self.assertNotIn("radius-only implementation", flat)
        self.assertNotIn("optional outer approximation", flat)
        self.assertNotIn("fallback", flat)

    def test_old_structured_fallback_file_is_removed(self):
        self.assertFalse((DOC / "w3d-stability-widening-structured.tex-part").exists())


if __name__ == "__main__":
    unittest.main()
