import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))

SPEC = importlib.util.spec_from_file_location(
    "ou3_neighborhood_radius_search", TOOLS / "ou3_neighborhood_radius_search.py"
)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MOD
SPEC.loader.exec_module(MOD)


class Ou3NeighborhoodRadiusSearchTests(unittest.TestCase):
    def base_case(self):
        return {
            "status": "FAIL_SAMPLED",
            "pass_sampled": False,
            "source_match_all": True,
            "measurement_acceptance_match_all": True,
            "theta_max_rad": 0.2,
            "W0": 1.0,
            "W1": 0.9,
        }

    def test_source_word_failure_has_priority(self):
        c = self.base_case()
        c["source_match_all"] = False
        c["measurement_acceptance_match_all"] = False
        self.assertEqual(MOD.classify_case_failure(c), "SOURCE_WORD_IDENTITY")

    def test_measurement_gate_failure_is_distinct(self):
        c = self.base_case()
        c["measurement_acceptance_match_all"] = False
        self.assertEqual(
            MOD.classify_case_failure(c), "MEASUREMENT_GATING_CONSISTENCY"
        )

    def test_chart_escape_is_distinct(self):
        c = self.base_case()
        c["theta_max_rad"] = 3.2
        self.assertEqual(MOD.classify_case_failure(c), "SO3_CHART_SAFETY")

    def test_loss_of_W_decrease_is_distinct(self):
        c = self.base_case()
        c["W1"] = 1.01
        self.assertEqual(MOD.classify_case_failure(c), "POSITIVE_W_DECREASE")

    def test_nonfinite_prefix_is_distinct(self):
        c = self.base_case()
        c["theta_max_rad"] = float("nan")
        self.assertEqual(MOD.classify_case_failure(c), "PREFIX_FINITE_SAFETY")

    def test_passed_case_has_no_failure(self):
        c = self.base_case()
        c["status"] = "PASS_SAMPLED"
        c["pass_sampled"] = True
        self.assertIsNone(MOD.classify_case_failure(c))

    def test_summary_keeps_sampled_qualification(self):
        c = self.base_case()
        c.update({"case": "x", "mode": "A", "direction": "theta_x"})
        report = {
            "status": "FAIL_OR_INCOMPLETE_SAMPLED",
            "case_count": 1,
            "valid_endpoint_case_count": 1,
            "cases": [c],
        }
        s = MOD.summarize_round(report, 4.0)
        self.assertFalse(s["pass_all_sampled"])
        self.assertEqual(s["failure_reason_counts"]["POSITIVE_W_DECREASE"], 1)


if __name__ == "__main__":
    unittest.main()
