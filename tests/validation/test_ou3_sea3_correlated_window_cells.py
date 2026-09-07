#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

from ou3_interval import Interval
import ou3_sea3_correlated_window_cells as CELLS


class Sea3CorrelatedWindowCellsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = CELLS.build()

    def test_status_is_structural_and_nonpromoting(self):
        self.assertEqual(CELLS.validate(self.payload), [])
        self.assertEqual(self.payload["canonical_source"], CELLS.CANONICAL_SOURCE)
        self.assertEqual(self.payload["P3_delta_preserved"], 1.0e-18)
        self.assertTrue(self.payload["correlated_window_source_cell_data_model_available"])
        self.assertTrue(self.payload["all_coupled_parent_constraints_retained_on_children"])
        self.assertTrue(self.payload["derived_shipping_coordinates_are_not_refinement_coordinates"])
        self.assertFalse(self.payload["hard_realization_oracle_closed_here"])
        self.assertFalse(self.payload["joint_source_output_map_closed_here"])
        self.assertFalse(self.payload["source_uniform_absolute_bias_path_closed_here"])
        self.assertFalse(self.payload["source_uniform_complete_window_cover_closed_here"])
        self.assertFalse(self.payload["P4_promoted_here"])

    def test_binary_split_retains_common_witness_and_constraints(self):
        parent = CELLS._root_smoke_cell()
        left, right = CELLS.split_source_cell(parent, "lambda.H1_fraction")
        self.assertEqual(
            CELLS.validate_binary_split(parent, (left, right), "lambda.H1_fraction"),
            [],
        )
        for child in (left, right):
            self.assertEqual(child.parent_cell_id, parent.cell_id)
            self.assertEqual(child.source_identity, parent.source_identity)
            self.assertEqual(
                child.hard_realization_witness_id, parent.hard_realization_witness_id
            )
            self.assertEqual(
                child.joint_response_witness_id, parent.joint_response_witness_id
            )
            self.assertEqual(
                child.shared_constraint_tokens[: len(parent.shared_constraint_tokens)],
                parent.shared_constraint_tokens,
            )

    def test_derived_shipping_coordinate_cannot_be_split_or_declared_source(self):
        parent = CELLS._root_smoke_cell()
        with self.assertRaises(KeyError):
            CELLS.split_source_cell(parent, "tau_committed")
        with self.assertRaises(ValueError):
            CELLS.SourceBound(
                "tau_committed",
                Interval.outward_bounds(1.0, 2.0),
                "illegal independent tuner hull",
            )
        with self.assertRaises(ValueError):
            CELLS.SourceBound(
                "R_S.current",
                Interval.outward_bounds(1.0, 2.0),
                "illegal independent R_S hull",
            )

    def test_source_split_does_not_cartesianize_coupled_partition_constraints(self):
        parent = CELLS._root_smoke_cell()
        left, right = CELLS.split_source_cell(
            parent, "lambda.H1_fraction", split_value=0.4
        )
        required = {
            "H_s^2=sum_r H_r^2",
            "active partition peak-steepness",
            "lambda_{k+1} in coupled Rhat_lambda(lambda_k)",
            "one common continuum phase history with no reseed",
            "one common joint translational/rotational response witness",
            "B^601_SEA3 common-witness membership",
        }
        self.assertTrue(required.issubset(set(left.shared_constraint_tokens)))
        self.assertTrue(required.issubset(set(right.shared_constraint_tokens)))
        self.assertEqual(left.bound("lambda.H1_fraction").interval.lo, 0.0)
        self.assertEqual(right.bound("lambda.H1_fraction").interval.hi, 1.0)

    def test_selector_attachment_requires_exact_bias_prefix_family(self):
        smoke = self.payload["smoke"]
        self.assertEqual(
            smoke["selector_prefixes_attached"], smoke["selector_lineage_length"]
        )
        self.assertTrue(smoke["point_bias_attached_same_history_every_prefix"])
        self.assertFalse(smoke["point_bias_history_is_source_uniform"])


if __name__ == "__main__":
    unittest.main()
