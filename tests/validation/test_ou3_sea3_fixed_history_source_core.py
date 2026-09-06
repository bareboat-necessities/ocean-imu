#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS / "stability"))

import ou3_sea3_fixed_history_source_core as CORE


class CompleteSea3FixedHistorySourceCoreTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = CORE.build()

    def test_fixed_history_is_legal_complete_sea3_member(self) -> None:
        d = self.payload
        self.assertEqual(CORE.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(d["source_membership"]["membership_is_by_operator_construction_not_quadrature"])
        self.assertTrue(d["source_membership"]["same_driver_field_entire_window"])
        self.assertTrue(d["source_membership"]["same_driver_field_translation_and_rotation"])
        self.assertEqual(d["source_membership"]["driver_norm"], 1.0)
        self.assertTrue(d["SEA3_fixed_history"]["partition_peak_steepness_admissible"])
        self.assertTrue(d["fixed_response_member"]["inside_declared_continuum_RAO_family"])

    def test_quadrature_is_only_evaluation_not_source_discretization(self) -> None:
        d = self.payload
        q = d["quadrature_diagnostic"]
        self.assertFalse(q["source_modes_are_quadrature_nodes"])
        self.assertLess(q["convergence_relative_to_peak"], 5e-5)
        self.assertFalse(d["fixed_response_member"]["finite_RAO_grid_used"])
        self.assertFalse(d["finite_harmonic_source_used"])
        self.assertFalse(d["trajectory_replay_used"])
        self.assertFalse(d["independent_sample_boxes_used"])

    def test_it_stops_before_frontend_seed_or_proof_promotion(self) -> None:
        d = self.payload
        self.assertEqual(d["sample_count"], 601)
        self.assertEqual(len(d["source_core"]), 601)
        self.assertFalse(d["front_end_entry_derived_from_same_history"])
        self.assertFalse(d["live_covariance_seed_derived_from_same_history"])
        self.assertFalse(d["complete_executor_artifact_materialized"])
        self.assertFalse(d["P4_promoted"])
        self.assertFalse(d["P5_promoted"])


if __name__ == "__main__":
    unittest.main()
