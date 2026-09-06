from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_riccati_metric_p3 as mod  # noqa: E402


class Sea3P3PromotionStateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_complete_sea3_and_literal_word_are_mandatory(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(self.d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(self.d["complete_SEA3_source_consumed"])
        self.assertTrue(self.d["literal_full_word_assembler_consumed"])
        self.assertTrue(self.d["literal_full_word_assembler_validation_pass"])
        self.assertTrue(self.d["literal_shipping_event_order_pass"])
        self.assertTrue(self.d["reset_complete_literal_execution_consumed"])
        self.assertTrue(self.d["no_fallback_route_enabled"])

    def test_pass_requires_universal_certificate_chain_not_finite_materialization(self):
        self.assertTrue(self.d["P3_CANONICAL_PASS"])
        self.assertFalse(self.d["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"])
        self.assertFalse(self.d["finite_source_family_materialization_required"])
        d = deepcopy(self.d)
        d["UNIVERSAL_COMPLETE_SEA3_CERTIFICATE_CHAIN_CLOSED"] = False
        failures = mod.validate(d)
        self.assertTrue(any("UNIVERSAL_COMPLETE_SEA3_CERTIFICATE_CHAIN_CLOSED" in x for x in failures), failures)

    def test_full_matrix_flag_cannot_be_removed_under_pass(self):
        d = deepcopy(self.d)
        d["P3_FULL_MATRIX_COMPARISON_CLOSED"] = False
        failures = mod.validate(d)
        self.assertTrue(any("P3_FULL_MATRIX_COMPARISON_CLOSED" in x for x in failures), failures)

    def test_fallback_flag_invalidates_gate(self):
        d = deepcopy(self.d)
        key = next(iter(d["no_fallback_generators"]))
        d["no_fallback_generators"][key] = True
        failures = mod.validate(d)
        self.assertTrue(any("fallback" in x for x in failures), failures)

    def test_p4_must_follow_p3_exactly(self):
        d = deepcopy(self.d)
        d["P3_CANONICAL_PASS"] = False
        d["P4_MAY_CONSUME_P3"] = True
        failures = mod.validate(d)
        self.assertTrue(any("P4 promotion" in x for x in failures), failures)

    def test_physical_left_inclusion_cannot_be_smuggled_into_conditional_p3(self):
        d = deepcopy(self.d)
        d["global_physical_deployment_left_inclusion_closed_here"] = True
        failures = mod.validate(d)
        self.assertTrue(any("global_physical_deployment_left_inclusion_closed_here" in x for x in failures), failures)


if __name__ == "__main__":
    unittest.main()
