"""Publication contract for the current, conditional OU-III stability study."""

import re
import unittest

from test_publication_references import STUDY, reachable_sources


class StabilityPublicationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.study = "\n".join(reachable_sources(STUDY).values())
        cls.flat = re.sub(r"\s+", " ", cls.study)

    def test_used_theorem_environments_are_declared(self):
        declared = set(re.findall(r"\\newtheorem\{([^}]+)\}", self.study))
        for name in ("lemma", "theorem", "proposition", "remark", "corollary"):
            if r"\begin{" + name + "}" in self.study:
                self.assertIn(name, declared, f"undefined LaTeX theorem environment: {name}")

    def test_simultaneous_physical_contracts_are_defined_in_the_study(self):
        self.assertIn(r"h\in{\cal M}_{\rm marine}\cap{\cal B}_{\rm IMU}\cap{\cal S}_{\rm mag}", self.flat)
        self.assertIn("The same \\(h\\) determines physical motion, delivered sensors", self.flat)
        self.assertIn("Complete physical stillness may persist arbitrarily long.", self.flat)
        self.assertIn("Physical history continues across proof-word boundaries.", self.flat)
        self.assertIn("e_a=b_{a,s}+b_{a,f},\\qquad e_g=b_{g,s}+b_{g,f}", self.flat)
        self.assertIn("The primitive history is carried through adjacent proof words and is never reset by the proof.",
                      self.flat)

    def test_only_applied_informative_magnetic_corrections_establish_service(self):
        self.assertIn("recurring accepted magnetic information", self.flat)
        self.assertIn("The proof uses the literal accepted-event chronology", self.flat)
        self.assertIn("does not replace it with continuous yaw observation", self.flat)

    def test_capture_and_the_tail_inherit_one_execution(self):
        self.assertIn("Generated filter quantities are functions of this history, not independent theorem inputs.",
                      self.flat)
        self.assertIn("advances the same mean and covariance history through the literal operation stream", self.flat)
        self.assertIn("No generated quantity is selected independently.", self.flat)
        for state in ("frontend and tuner state", "covariance", "accepted events", "scheduler state",
                      "estimator mean"):
            self.assertIn(state, self.flat)

    def test_finite_error_statement_preserves_metric_and_supply(self):
        for marker in (
            r"e_N=Me_0+b",
            r"G_\gamma=J_0-M^TJ_NM-\gamma J_0",
            r"\chi_\gamma=b^TJ_Nb+z^TG_\gamma^{-1}z",
            r"\sqrt V\le0.15",
        ):
            self.assertIn(marker, self.flat)
        self.assertIn("Homogeneous contraction alone is insufficient", self.flat)
        self.assertIn("at every required prefix", self.flat)

    def test_finite_error_statement_and_qualification_remain_conditional(self):
        self.assertIn("None of these labels promotes an end-to-end theorem.", self.flat)
        self.assertIn("--- OPEN]", self.flat)
        self.assertIn("regional practical stability is not yet claimed", self.flat)
        self.assertIn("The regional theorem is intentionally stated as OPEN", self.flat)
        self.assertIn("Numerical diagnostics are evidence only when their provenance identifies the shipping source",
                      self.flat)


if __name__ == "__main__":
    unittest.main()
