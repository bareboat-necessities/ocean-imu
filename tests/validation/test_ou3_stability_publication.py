"""Publication contract for the current, conditional OU-III stability study."""

import re
import unittest

from test_publication_references import STUDY, reachable_sources


class StabilityPublicationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.study = "\n".join(reachable_sources(STUDY).values())
        cls.flat = re.sub(r"\s+", " ", cls.study)

    def test_simultaneous_physical_contracts_are_defined_in_the_study(self):
        for marker in (
            r"\label{eq:intersection}", r"\label{eq:pointwise}",
            r"\label{eq:pac}", r"\label{eq:bias-cont}",
            r"\label{eq:bias-discrete}", r"\label{eq:mag-service}",
        ):
            self.assertIn(marker, self.study)
        self.assertIn("same-history requirement", self.flat)
        self.assertIn("Quiet water", self.flat)
        self.assertIn("Permanent displacement DC is inadmissible", self.flat)
        self.assertIn("not restarted at word boundaries", self.flat)

    def test_only_applied_informative_magnetic_corrections_establish_service(self):
        self.assertIn(r"\label{eq:G}", self.study)
        self.assertIn(r"\succeq\mu_M I_2", self.study)
        self.assertIn("whose shipping correction was actually applied", self.flat)
        self.assertIn("A maximum event gap alone is insufficient", self.flat)
        self.assertIn(r"\label{eq:unipotent}", self.study)
        self.assertIn("not a second stability theorem", self.flat)

    def test_capture_and_the_tail_inherit_one_execution(self):
        self.assertIn(r"T_c=T_c(h,x_0)<\infty", self.study)
        self.assertIn("No common finite startup deadline is assumed", self.flat)
        self.assertIn("without reseeding", self.flat)
        for state in (
            "full covariance", "physical bias histories", "frontend",
            "tuner state", "committed parameters", "magnetic reference state",
            "scheduler", "clocks", "physical source continuation",
        ):
            self.assertIn(state, self.flat)

    def test_finite_error_statement_and_qualification_remain_conditional(self):
        for marker in (
            r"\label{eq:diss}", r"0\le\rho<1",
            r"m\norm{e_j}^2\le V_j\le M\norm{e_j}^2",
            r"\label{tab:limits}",
        ):
            self.assertIn(marker, self.study.replace("&", ""))
        self.assertIn("every-prefix bounds and forward retention", self.flat)
        self.assertIn("This implication is conditional until", self.flat)
        self.assertIn("Finite simulations do not establish", self.flat)
        self.assertIn("assembled-history qualification remains open", self.flat)
        self.assertIn("End-to-end regional practical stability remains", self.flat)


if __name__ == "__main__":
    unittest.main()
