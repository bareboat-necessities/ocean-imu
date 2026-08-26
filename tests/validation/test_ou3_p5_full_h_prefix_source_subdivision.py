import sys
from pathlib import Path
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_full_h_prefix_source_subdivision as S


class Ou3P5FullHPrefixSourceSubdivisionTests(unittest.TestCase):
    def _fake_v3(self, *, closed=True, first_failure=None):
        return {
            "schema": S.V3.SCHEMA,
            "P5_FULL_H_PREFIX_MATRIX_CERTIFICATE": "PASS" if closed else "NOT_ESTABLISHED",
            "complete_q_le_8_prefix_family_closed": closed,
            "max_reached_cayley_norm_upper": 1.25 if closed else 7.5,
            "first_failure": first_failure,
            "full_18x18_covariance_propagated": True,
            "H_R_S_K_r_d_eff_recomputed_in_same_prefix_cell": True,
            "shipping_Joseph_update_used": True,
            "immediate_left_error_reset_congruence_used": True,
            "deployed_quaternion_composed_before_result_cayley": True,
        }

    def test_source_cover_has_cartesian_cardinality_and_tau_dependent_cadence(self):
        cells = S._source_children(2)
        self.assertEqual(len(cells), 8)
        taus = {tuple(c["tau_s"].as_list()) for c in cells}
        periods = {tuple(c["pseudo_period_s"].as_list()) for c in cells}
        self.assertEqual(len(taus), 2)
        self.assertGreaterEqual(len(periods), 2)

    def test_build_cell_keeps_selected_source_cell_disjoint(self):
        seen = {}

        def fake_build(_domain):
            seen["source"] = S._serialize_source(S.V1._source_cell())
            return self._fake_v3(closed=True)

        with mock.patch.object(S.V3, "build", side_effect=fake_build), \
             mock.patch.object(S.V3, "validate", return_value=[]):
            d = S.build_cell(source_pieces=2, source_cell_index=3)
        self.assertEqual(S.validate_cell(d), [])
        self.assertEqual(d["selected_source_cell"], seen["source"])
        self.assertFalse(d["source_partition_cells_hulled_together"])
        self.assertEqual(d["P5_FULL_H_PREFIX_SOURCE_CELL_CERTIFICATE"], "PASS")

    def test_nonclosure_is_fail_closed_with_same_witness(self):
        witness = {"sample": 4, "operation": "accelerometer", "reason": "cell too wide"}
        with mock.patch.object(S.V3, "build", return_value=self._fake_v3(closed=False, first_failure=witness)), \
             mock.patch.object(S.V3, "validate", return_value=[]):
            d = S.build_cell(source_pieces=1, source_cell_index=0)
        self.assertEqual(S.validate_cell(d), [])
        self.assertEqual(d["P5_FULL_H_PREFIX_SOURCE_CELL_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertEqual(d["first_failure"], witness)
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])

    def test_aggregate_requires_every_child_and_all_closed(self):
        children = []
        for i in range(8):
            d = self._fake_v3(closed=True)
            d.update({
                "schema": S.CHILD_SCHEMA,
                "qualification": "OU3_P5_FULL_H_PREFIX_SOURCE_SUBDIVISION_CHILD",
                "source_partition_is_finite_cover": True,
                "source_partition_cells_hulled_together": False,
                "source_partition_pieces_per_axis": 2,
                "source_partition_total_cells": 8,
                "source_partition_cell_index": i,
                "selected_source_cell": {
                    "dt_s": 0.005,
                    "tau_s": [1.0, 2.0],
                    "sigma_aw_mps2": [0.1, 0.2],
                    "R_S_filter_std": [0.1, 0.2],
                    "pseudo_period_s": [0.1, 0.2],
                },
                "v3_validation_failures": [],
                "same_full_18x18_Joseph_reset_quaternion_map_as_v3": True,
                "filter_changed": False,
                "whole_word_promoted_here": False,
                "N_H_words_set_here": False,
                "P5_FULL_H_PREFIX_SOURCE_CELL_CERTIFICATE": "PASS",
                "next_obligation": "AGGREGATE_ALL_SOURCE_PARTITION_CELLS_AND_CHECK_COMPLETE_Q8_CLOSURE",
            })
            children.append(d)
        agg = S.aggregate(children, source_pieces=2)
        self.assertEqual(S.validate_aggregate(agg), [])
        self.assertEqual(agg["P5_FULL_H_PREFIX_SOURCE_SUBDIVISION_CERTIFICATE"], "PASS")
        self.assertTrue(agg["complete_q_le_8_prefix_family_closed_over_source_partition"])
        self.assertFalse(agg["whole_word_promoted_here"])

        missing = S.aggregate(children[:-1], source_pieces=2)
        self.assertEqual(missing["P5_FULL_H_PREFIX_SOURCE_SUBDIVISION_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertTrue(any("missing source cells" in x for x in S.validate_aggregate(missing)))


if __name__ == "__main__":
    unittest.main()
