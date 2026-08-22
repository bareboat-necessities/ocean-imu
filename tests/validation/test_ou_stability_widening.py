"""Publication/source contract for the OU--III analytical stability widening."""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", text)


class OUIIIStabilityWideningContractTests(unittest.TestCase):
    def test_phase_a_is_wired_into_manuscript(self):
        main = _read(DOC / "kalman_ou-w3d.tex")
        self.assertIn(r"\input{w3d-analytical-stability-widening.tex-part}", main)

    def test_phase_a_has_explicit_vector_information_bound(self):
        proof = _read(DOC / "w3d-analytical-stability-widening.tex-part")
        for marker in (
            r"\label{eq:widen-vector-packet}",
            r"\label{eq:widen-vector-mu}",
            r"\label{lem:widen-vector-coercivity}",
            r"\label{eq:widen-gamma-lower}",
            r"\label{lem:widen-two-packet-bg}",
            r"\label{eq:widen-alpha6}",
        ):
            self.assertIn(marker, proof)
        self.assertIn(r"1-\sqrt{1-s_{fm}^2}", proof)
        self.assertIn(r"\alpha_6^{\rm an}", proof)

    def test_phase_a_keeps_general_pe_as_fallback(self):
        proof = _flat(_read(DOC / "w3d-analytical-stability-widening.tex-part"))
        self.assertIn("not a replacement for the general PE route", proof)
        self.assertIn("original Gramian condition", proof)
        self.assertIn("random-walk limit", proof)

    def test_phase_a_envelope_does_not_depend_on_adaptation_law(self):
        proof = _read(DOC / "w3d-analytical-stability-widening.tex-part")
        flat = _flat(proof)
        for marker in (
            r"\Pi^{\rm env}",
            r"\label{eq:widen-envelope-bounds}",
            r"\label{lem:widen-envelope-translation}",
            r"\label{thm:widen-envelope-ues}",
        ):
            self.assertIn(marker, proof)
        self.assertIn("No EMA recurrence", flat)
        self.assertIn("does not depend on the adaptation-law exponents", flat)
        self.assertIn("performance mechanism rather than a stability-critical feedback law", flat)


if __name__ == "__main__":
    unittest.main()
