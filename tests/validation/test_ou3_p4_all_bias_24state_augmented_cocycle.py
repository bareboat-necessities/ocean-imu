from __future__ import annotations

import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

from ou3_interval import Interval
import ou3_p4_all_bias_24state_augmented_cocycle as COCYCLE
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS

I = Interval.point


class AllBias24StateAugmentedCocycleTests(unittest.TestCase):
    def _selector_fixture(self):
        frontend = SELECTORS.FRONTEND._point_state()
        P0_H, P0_A, _ = SELECTORS._live_structured_point_covariance_fixture(
            frontend, COCYCLE.GRAPH.DEFAULT_DOMAIN
        )
        sample = SELECTORS._point_sample()
        endpoints, selectors, _ = SELECTORS.execute_with_prefix_selectors(
            frontend_entry=frontend,
            P0_H=P0_H,
            P0_A=P0_A,
            samples=[sample, sample],
            domain_path=COCYCLE.GRAPH.DEFAULT_DOMAIN,
            branch_limit=128,
        )
        return selectors, endpoints[0].source_cell_id

    def test_status_is_all_family_and_nonpromoting(self):
        d = COCYCLE.build()
        self.assertEqual(COCYCLE.validate(d), [])
        self.assertEqual(set(d["family_lifts"]), {"BIAS0", "BIAS1", "BIAS2"})
        self.assertTrue(d["source_reachable_selector_family_relation_consumed"])
        self.assertTrue(d["suffix_propagated_prediction_supply_map_materializer_available"])
        self.assertFalse(d["source_uniform_joint_augmented_LDLT_closed_here"])
        self.assertFalse(d["P4_PASS"])

    def test_two_prefix_cocycle_materializes_all_bias_families(self):
        selectors, endpoint = self._selector_fixture()
        state = [I(0.0)] * 21
        for family in ("BIAS0", "BIAS1", "BIAS2"):
            r = COCYCLE.materialize_A21_joint_cocycle(
                selectors, endpoint, family=family, initial_error_state=state
            )
            self.assertEqual(len(r.A_word), 24)
            self.assertEqual(len(r.A_word[0]), 24)
            self.assertEqual(len(r.B_word), 24)
            self.assertEqual(len(r.B_word[0]), 12)  # two predictions x six supply coordinates
            self.assertEqual(r.predictions, 2)
            self.assertGreaterEqual(r.measurements, 2)
            self.assertEqual(r.prediction_source_column_ranges, ((0, 6), (6, 12)))
            self.assertTrue(all(x != "not_applicable" for x in r.projection_branches))

    def test_unknown_bias_family_rejected(self):
        selectors, endpoint = self._selector_fixture()
        with self.assertRaises(ValueError):
            COCYCLE.materialize_A21_joint_cocycle(
                selectors, endpoint, family="BIASX", initial_error_state=[I(0.0)] * 21
            )

    def test_wrong_initial_dimension_rejected(self):
        selectors, endpoint = self._selector_fixture()
        with self.assertRaises(ValueError):
            COCYCLE.materialize_A21_joint_cocycle(
                selectors, endpoint, family="BIAS0", initial_error_state=[I(0.0)] * 18
            )


if __name__ == "__main__":
    unittest.main()
