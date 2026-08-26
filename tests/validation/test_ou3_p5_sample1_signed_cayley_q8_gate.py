from __future__ import annotations
import sys
from pathlib import Path
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_signed_cayley_q8_gate as G


class Sample1V13EToV14DGateTests(unittest.TestCase):
    def test_v13e_nonclosure_blocks_v14d_and_keeps_witness(self):
        witness = {"source_row": 3, "radial_upper_rad": 6.2}
        radial = {
            "P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13E": "NOT_ESTABLISHED",
            "evaluated_signed_subcells": 10,
            "above_6rad_subcells": 2,
            "unclosed_radial_subcells": 1,
            "max_radial_upper": 6.2,
            "minimum_radial_lower_above_6": 5.9,
            "first_unclosed_radial_subcell": witness,
            "worst_radial_subcell": witness,
        }
        with mock.patch.object(G.V13E, "build", return_value=radial), \
             mock.patch.object(G.V13E, "validate", return_value=[]), \
             mock.patch.object(G.V14D, "build") as q8_build:
            d = G.build()
        q8_build.assert_not_called()
        self.assertEqual(G.validate(d), [])
        self.assertEqual(d["P5_SAMPLE1_V13E_TO_V14D_GATE"], "NOT_ESTABLISHED")
        self.assertFalse(d["V14D_invoked"])
        self.assertEqual(d["V13E_first_unclosed_radial_subcell"], witness)
        self.assertFalse(d["N_H_words_set_here"])

    def test_v13e_pass_invokes_v14d_and_propagates_q8_witness(self):
        radial = {
            "P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13E": "PASS",
            "evaluated_signed_subcells": 10,
            "above_6rad_subcells": 0,
            "unclosed_radial_subcells": 0,
            "max_radial_upper": 5.5,
            "minimum_radial_lower_above_6": None,
            "first_unclosed_radial_subcell": None,
            "worst_radial_subcell": {"radial_upper_rad": 5.5},
        }
        witness = {"post_sample1_cayley_norm_upper": 9.1}
        q8 = {
            "P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14D": "NOT_ESTABLISHED",
            "evaluated_signed_cayley_cells": 10,
            "product_scalar_antipode_cells": 1,
            "unclosed_q8_cells": 1,
            "minimum_abs_product_scalar_lower": 0.0,
            "max_post_sample1_cayley_norm_upper": float("inf"),
            "first_unclosed_q8_cell": witness,
            "worst_q8_cell": witness,
        }
        with mock.patch.object(G.V13E, "build", return_value=radial), \
             mock.patch.object(G.V13E, "validate", return_value=[]), \
             mock.patch.object(G.V14D, "build", return_value=q8), \
             mock.patch.object(G.V14D, "validate", return_value=[]):
            d = G.build()
        self.assertEqual(G.validate(d), [])
        self.assertTrue(d["V14D_invoked"])
        self.assertEqual(d["V14D_first_unclosed_q8_cell"], witness)
        self.assertEqual(d["P5_SAMPLE1_V13E_TO_V14D_GATE"], "NOT_ESTABLISHED")

    def test_v14d_pass_advances_only_to_source_phase_lift(self):
        radial = {
            "P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13E": "PASS",
            "evaluated_signed_subcells": 10,
            "above_6rad_subcells": 0,
            "unclosed_radial_subcells": 0,
            "max_radial_upper": 5.5,
            "minimum_radial_lower_above_6": None,
            "first_unclosed_radial_subcell": None,
            "worst_radial_subcell": None,
        }
        q8 = {
            "P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14D": "PASS",
            "evaluated_signed_cayley_cells": 10,
            "product_scalar_antipode_cells": 0,
            "unclosed_q8_cells": 0,
            "minimum_abs_product_scalar_lower": 0.7,
            "max_post_sample1_cayley_norm_upper": 4.0,
            "first_unclosed_q8_cell": None,
            "worst_q8_cell": None,
        }
        with mock.patch.object(G.V13E, "build", return_value=radial), \
             mock.patch.object(G.V13E, "validate", return_value=[]), \
             mock.patch.object(G.V14D, "build", return_value=q8), \
             mock.patch.object(G.V14D, "validate", return_value=[]):
            d = G.build()
        self.assertEqual(G.validate(d), [])
        self.assertEqual(d["P5_SAMPLE1_V13E_TO_V14D_GATE"], "PASS")
        self.assertIn("LIFT_CLOSED_SAMPLE1_PREFIX", d["next_obligation"])
        self.assertFalse(d["q8_word_promoted_here"])
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])


if __name__ == "__main__":
    unittest.main()
