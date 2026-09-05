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

    def test_literal_full_word_assembler_is_mandatory(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertTrue(self.d["literal_full_word_assembler_consumed"])
        self.assertTrue(self.d["literal_full_word_assembler_validation_pass"])
        self.assertTrue(self.d["literal_shipping_event_order_pass"])
        self.assertTrue(self.d["literal_same_source_state_feeds_F_Q_TS_RS"])
        self.assertTrue(self.d["literal_every_valid_accelerometer_required"])
        self.assertTrue(self.d["literal_all_due_S_required"])
        self.assertTrue(self.d["literal_async_magnetometer_family_retained"])
        self.assertTrue(self.d["literal_aw_floor_event_family_retained"])

    def test_pass_cannot_be_flipped_without_full_HA_closure(self):
        d = deepcopy(self.d)
        d["P3_CANONICAL_PASS"] = True
        failures = mod.validate(d)
        self.assertTrue(any("inconsistent" in x or "not passed" in x for x in failures), failures)

    def test_closed_stage_cannot_be_claimed_while_canonical_pass_is_false(self):
        d = deepcopy(self.d)
        d["P3_FULL_WORD_ENCLOSED"] = True
        failures = mod.validate(d)
        self.assertTrue(any("inconsistent promotion flags" in x for x in failures), failures)


if __name__ == "__main__":
    unittest.main()
