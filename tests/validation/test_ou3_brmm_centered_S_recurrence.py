import copy
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_centered_S_recurrence as SREC


class CenteredSRecurrenceTests(unittest.TestCase):
    def test_materialized_source_condition_is_fail_closed(self):
        d = SREC.build()
        self.assertEqual(SREC.validate(d), [])
        self.assertTrue(d["fresh_common_origin_reduction_preserved"])
        self.assertTrue(d["condition_is_origin_invariant"])
        self.assertTrue(d["condition_bounds_every_handoff_centered_S"])
        self.assertTrue(d["condition_excludes_nonzero_constant_position_history"])
        self.assertIsNone(d["D_S_max_m_s"])
        self.assertFalse(d["D_S_numeric_qualification_closed"])
        self.assertFalse(d["COMPLETE_BRMM_INDEFINITE_S_QUALIFIED"])
        self.assertFalse(d["P4_PASS"])
        self.assertFalse(d["P5_MAY_START"])
        self.assertEqual(d["P3_delta"], 1e-18)

    def test_legacy_entry_ball_cannot_be_smuggled_back(self):
        d = SREC.build()
        d["legacy_300_m_s_fresh_entry_ball_used"] = True
        self.assertIn("legacy_300_m_s_fresh_entry_ball_used not false", SREC.validate(d))

    def test_numeric_source_bound_cannot_be_fabricated(self):
        d = SREC.build()
        d["D_S_max_m_s"] = 300.0
        self.assertIn("D_S_max must remain unfrozen until retention computes it", SREC.validate(d))

    def test_false_promotion_is_rejected(self):
        d = SREC.build()
        d["P4_PASS"] = True
        d["COMPLETE_BRMM_INDEFINITE_S_QUALIFIED"] = True
        f = SREC.validate(d)
        self.assertIn("P4_PASS not false", f)
        self.assertIn("COMPLETE_BRMM_INDEFINITE_S_QUALIFIED not false", f)

    def test_wordwise_rezero_is_rejected(self):
        d = SREC.build()
        d["wordwise_rezero_of_S_used"] = True
        self.assertIn("wordwise_rezero_of_S_used not false", SREC.validate(d))


if __name__ == "__main__":
    unittest.main()
