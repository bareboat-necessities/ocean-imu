import math
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from ou3_interval import Interval
import ou3_p5_sample1_v41_adaptive_directional_current_v42 as V42


class Sample1V41AdaptiveDirectionalCurrentV42Tests(unittest.TestCase):
    def test_binary_split_is_an_exact_closed_cover(self):
        x = Interval(-1.25, 2.75)
        parts = V42._split_interval(x)
        self.assertIsNotNone(parts)
        left, right = parts
        self.assertEqual(left.lo, x.lo)
        self.assertEqual(right.hi, x.hi)
        self.assertEqual(left.hi, right.lo)
        self.assertGreater(left.hi, x.lo)
        self.assertLess(right.lo, x.hi)

    def test_directional_split_selects_largest_dot_uncertainty(self):
        chart = {
            "cx": Interval(-1.0, 1.0),
            "cy": Interval(-2.0, 2.0),
            "cz": Interval(-3.0, 3.0),
            "cyz_norm_upper": 4.0,
        }
        vd = [Interval(2.0, 2.0), Interval(0.25, 0.25), Interval(0.1, 0.1)]
        # dot-width scores are 4, 1, and 0.6 respectively.
        self.assertEqual(V42._split_dimension(vd, chart), 0)

    def test_adaptive_cover_is_nonworsening_and_union_safe(self):
        chart = {
            "cx": Interval(-0.5, 0.5),
            "cy": Interval(-0.2, 0.2),
            "cz": Interval(-0.1, 0.1),
            "cyz_norm_upper": math.sqrt(0.05),
        }
        wd = Interval(0.35, 0.35)
        vd = [Interval(0.6, 0.6), Interval(0.15, 0.15), Interval(0.05, 0.05)]
        base = V42._leaf_product(
            q_parent=0.6, wd=wd, vd=vd, chart=chart,
            support_fn=V42.V18._support_product_scalar,
            qplus_fn=V42.V14._qplus_from_product_scalar,
        )
        self.assertIsNotNone(base)
        refined = V42._adaptive_product_cover(
            q_parent=0.6, wd=wd, vd=vd, chart=chart, max_depth=2,
            support_fn=V42.V18._support_product_scalar,
            qplus_fn=V42.V14._qplus_from_product_scalar,
        )
        self.assertFalse(refined["empty"])
        self.assertGreaterEqual(refined["leaf_evaluations"], 1)
        self.assertLessEqual(refined["q_upper"], base["q_upper"])
        self.assertGreaterEqual(refined["abs_W_lower"], base["abs_W_lower"])

    def test_contract_keeps_filter_limit_and_promotion_fail_closed(self):
        d = {
            "schema": V42.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_V41_ADAPTIVE_DIRECTIONAL_CURRENT_V42",
            "source_generated_not_trajectory_fit": True,
            "source_replay_used": False,
            "filter_changed": False,
            "V41_full_source_cell0_parent_retained": True,
            "V40_exact_Joseph_first_PSD_retained_through_V41": True,
            "adaptive_directional_current_subdivision_used": True,
            "adaptive_children_outside_parent_q_ball_discarded_only_by_norm_lower": True,
            "adaptive_union_uses_min_abs_W_and_max_q": True,
            "temporary_adaptive_helpers_restored_after_build": True,
            "adaptive_depth": 1,
            "adaptive_eligible_open_support_calls": 1,
            "deployed_correction_limit_rad": 6.0,
            "deployed_correction_limit_increased": False,
            "q_target": V42.Q_TARGET,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "P5_SAMPLE1_V41_ADAPTIVE_DIRECTIONAL_CURRENT_V42": "NOT_ESTABLISHED",
            "first_unclosed_q8_cell": {"post_sample1_cayley_norm_upper": 8.1},
            "failures": [],
        }
        self.assertEqual(V42.validate(d), [])


if __name__ == "__main__":
    unittest.main()
