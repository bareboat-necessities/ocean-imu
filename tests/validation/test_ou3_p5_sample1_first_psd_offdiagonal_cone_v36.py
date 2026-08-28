import pathlib
import sys
import unittest

TOOLS = pathlib.Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import ou3_p5_sample1_first_psd_offdiagonal_cone_v36 as V36


class Sample1FirstPSDOffdiagonalConeV36Tests(unittest.TestCase):
    def test_psd_offdiag_operator_bound_is_eps(self):
        eps = 2.5e-7
        got = V36._psd_offdiag_operator_upper(eps)
        self.assertGreaterEqual(got, eps)
        self.assertLessEqual(got, eps * (1.0 + 1e-12))
        self.assertLess(got, 2.0 * eps)

    def test_invalid_eps_rejected(self):
        with self.assertRaises(ValueError):
            V36._psd_offdiag_operator_upper(-1.0)

    def test_validation_keeps_promotion_guards(self):
        d = {
            "schema": V36.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_FIRST_PSD_OFFDIAGONAL_CONE_V36",
            "source_generated_not_trajectory_fit": True,
            "V12D_baseline_revalidated": True,
            "PSD_remainder_order_interval_0_E_epsI_used": True,
            "PSD_diagonal_absorbed_in_existing_t_Y_intervals": True,
            "PSD_principal_minor_offdiag_abs_le_eps_over_2": True,
            "PSD_zero_diagonal_remainder_operator_le_eps": True,
            "source_replay_used": False,
            "filter_changed": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "P5_SAMPLE1_FIRST_PSD_OFFDIAGONAL_CONE_V36": "PASS",
            "failures": [],
        }
        self.assertEqual(V36.validate(d), [])
        d["N_H_words_set_here"] = True
        self.assertIn("N_H_words_set_here is not false", V36.validate(d))


if __name__ == "__main__":
    unittest.main()
