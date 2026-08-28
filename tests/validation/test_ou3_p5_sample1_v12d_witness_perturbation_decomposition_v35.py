import pathlib
import sys
import unittest

TOOLS = pathlib.Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import ou3_p5_sample1_v12d_witness_perturbation_decomposition_v35 as V35


class Sample1V12DWitnessPerturbationDecompositionV35Tests(unittest.TestCase):
    def test_fraction(self):
        self.assertEqual(V35._frac(0.0, 0.0), 0.0)
        self.assertAlmostEqual(V35._frac(1.0, 4.0), 0.25)

    def test_validation_keeps_promotion_guards(self):
        d = {
            "schema": V35.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_V12D_WITNESS_PERTURBATION_DECOMPOSITION_V35",
            "source_generated_not_trajectory_fit": True,
            "V12D_parent_revalidated": True,
            "source_replay_used": False,
            "filter_changed": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "P5_SAMPLE1_V12D_WITNESS_PERTURBATION_DECOMPOSITION_V35": "PASS",
            "failures": [],
        }
        self.assertEqual(V35.validate(d), [])
        d["N_H_words_set_here"] = True
        self.assertIn("N_H_words_set_here is not false", V35.validate(d))


if __name__ == "__main__":
    unittest.main()
