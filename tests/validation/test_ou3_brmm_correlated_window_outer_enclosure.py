import copy
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_correlated_window_outer_enclosure as O


class CorrelatedWindowOuterEnclosureTests(unittest.TestCase):
    def test_status_closes_left_inclusion_but_not_executor_or_p4(self):
        d = O.build()
        self.assertEqual(O.validate(d), [])
        self.assertTrue(d["left_inclusion_closed"])
        self.assertTrue(d["validated_correlated_outer_enclosure_closed"])
        self.assertFalse(d["executor_601_sample_relation_materialized_here"])
        self.assertFalse(d["P4_PASS"])
        self.assertEqual(d["P3_delta"], 1e-18)

    def test_same_history_correlations_are_retained(self):
        d = O.build()
        self.assertTrue(d["correlation_retained_across_samples"])
        self.assertTrue(d["correlation_retained_across_axes"])
        self.assertTrue(d["correlation_retained_between_translation_and_moments"])
        self.assertFalse(d["independent_sample_boxes_used"])
        self.assertFalse(d["independent_axis_boxes_used"])
        t = d["constraints"]["translation"]
        self.assertTrue(t["primitive_out_is_next_primitive_in"])
        self.assertTrue(t["one_live_centered_S_origin_for_all_samples"])
        m = d["constraints"]["acceleration_moments"]
        self.assertTrue(m["same_transition_witness_as_translation"])
        self.assertTrue(m["independent_J0_J1_J2_forbidden"])
        self.assertTrue(m["independent_axes_forbidden"])

    def test_detaching_moment_witness_rejected(self):
        d = O.build()
        d["constraints"]["acceleration_moments"]["same_transition_witness_as_translation"] = False
        self.assertIn("moment witness detached from translation", O.validate(d))

    def test_cartesianization_rejected(self):
        d = O.build()
        d["independent_sample_boxes_used"] = True
        d["independent_axis_boxes_used"] = True
        f = O.validate(d)
        self.assertIn("independent_sample_boxes_used not false", f)
        self.assertIn("independent_axis_boxes_used not false", f)

    def test_centered_origin_reset_rejected(self):
        d = O.build()
        d["constraints"]["translation"]["one_live_centered_S_origin_for_all_samples"] = False
        self.assertIn("centered-S origin continuity lost", O.validate(d))

    def test_false_promotion_rejected(self):
        d = O.build()
        d["P4_PASS"] = True
        d["executor_601_sample_relation_materialized_here"] = True
        f = O.validate(d)
        self.assertIn("P4_PASS not false", f)
        self.assertIn("executor_601_sample_relation_materialized_here not false", f)


if __name__ == "__main__":
    unittest.main()
