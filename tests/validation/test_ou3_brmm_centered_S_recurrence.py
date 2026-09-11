import copy
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_centered_S_recurrence as SREC


class CenteredSRecurrenceTests(unittest.TestCase):
    def test_materialized_source_condition_is_numeric_and_fail_closed_for_P4(self):
        d = SREC.build()
        self.assertEqual(SREC.validate(d), [])
        self.assertTrue(d["fresh_common_origin_reduction_preserved"])
        self.assertTrue(d["condition_is_origin_invariant"])
        self.assertTrue(d["condition_bounds_every_handoff_centered_S"])
        self.assertTrue(d["condition_excludes_nonzero_constant_position_history"])
        self.assertEqual(d["D_S_max_m_s"], 1100.0)
        self.assertTrue(d["COMPLETE_BRMM_INDEFINITE_S_THEOREM_CLOSED"])
        self.assertFalse(d["working_radius_sets_physical_D_S"])
        self.assertNotIn("D_S_retention_search_upper_m_s", d)
        self.assertTrue(d["D_S_numeric_qualification_closed"])
        self.assertTrue(d["COMPLETE_BRMM_INDEFINITE_S_QUALIFIED"])
        self.assertFalse(d["P4_PASS"])
        self.assertFalse(d["P5_MAY_START"])
        self.assertEqual(d["P3_delta"], 1e-18)

    def test_legacy_entry_ball_cannot_be_smuggled_back(self):
        d = SREC.build()
        d["legacy_300_m_s_fresh_entry_ball_used"] = True
        self.assertIn("legacy_300_m_s_fresh_entry_ball_used differs from derived centered-S recurrence", SREC.validate(d))

    def test_numeric_source_bound_cannot_be_fabricated(self):
        d = SREC.build()
        d["D_S_max_m_s"] = 300.0
        self.assertIn("D_S_max_m_s differs from derived centered-S recurrence", SREC.validate(d))

    def test_false_p4_promotion_is_rejected_even_with_numeric_source_bound(self):
        d = SREC.build()
        d["P4_PASS"] = True
        f = SREC.validate(d)
        self.assertIn("P4_PASS differs from derived centered-S recurrence", f)

    def test_wordwise_rezero_is_rejected(self):
        d = SREC.build()
        d["wordwise_rezero_of_S_used"] = True
        self.assertIn("wordwise_rezero_of_S_used differs from derived centered-S recurrence", SREC.validate(d))


if __name__ == "__main__":
    unittest.main()
