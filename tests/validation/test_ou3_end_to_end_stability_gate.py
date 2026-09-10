#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_end_to_end_stability_gate as E2E


class EndToEndStabilityGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = E2E.build()

    def test_validates_and_stays_fail_closed(self):
        self.assertEqual(E2E.validate(self.d), [])
        self.assertFalse(self.d["END_TO_END_STABILITY_PASS"])
        self.assertFalse(self.d["END_TO_END_DEPLOYMENT_PASS"])
        self.assertFalse(self.d["P4_PASS"])
        self.assertFalse(self.d["P5_MAY_START"])
        self.assertFalse(self.d["P5_CAPTURE_PASS"])
        self.assertFalse(self.d["P4_INVARIANCE_PASS"])

    def test_theorem_is_the_end_to_end_statement_not_P4(self):
        self.assertIn("finite time", self.d["theorem"])
        self.assertIn("practical ISS", self.d["theorem"])
        self.assertIn("BIAS0/BIAS1/BIAS2", self.d["theorem"])

    def test_promotion_is_the_exact_conjunction(self):
        d = dict(self.d)
        d["P5_CAPTURE_PASS"] = True
        d["END_TO_END_STABILITY_PASS"] = True
        self.assertIn("P5 capture promotion is not the exact conjunction", E2E.validate(d))
        self.assertIn("end-to-end passed while P4 is false", E2E.validate(d))

    def test_every_open_capture_obligation_carries_a_failure_class(self):
        capture = self.d["P5_capture_obligations"]
        self.assertTrue(capture)
        open_names = [k for k, v in capture.items() if not v["established"]]
        self.assertTrue(open_names)
        for name in open_names:
            self.assertIn(capture[name]["class"], ("A", "B", "C", "D", "E"), name)
        self.assertEqual(set(self.d["remaining_P5_capture_obligations"]), set(open_names))

    def test_cold_to_tuner_warm_is_the_one_discharged_capture_step(self):
        row = self.d["P5_capture_obligations"]["cold_to_tuner_warm_finite_time"]
        self.assertTrue(row["established"])
        self.assertIsNone(row["class"])
        self.assertGreater(row["bound_s"], 0.0)

    def test_no_startup_to_live_teleport_is_permitted(self):
        cov = self.d["transition_coverage"]
        self.assertFalse(cov["startup_to_live_teleport_permitted"])
        self.assertTrue(cov["live_entry_pinning_from_deployed_code"])
        self.assertTrue(cov["linear_block_unpropagated_before_go_live"])
        self.assertTrue(cov["attitude_linear_cross_zeroed_at_entry"])
        self.assertEqual(tuple(cov["deployed_stages_enumerated"]), E2E.DEPLOYED_STAGES)
        self.assertEqual(tuple(cov["deployed_live_events_enumerated"]), E2E.DEPLOYED_LIVE_EVENTS)

    def test_basin_is_carried_as_a_maximisation(self):
        basin = self.d["certified_basin"]
        self.assertIn("maximise", basin["objective"])
        self.assertFalse(basin["entry_radii_reduced_for_proof_convenience"])
        self.assertTrue(basin["frontier_is_a_frozen_map_diagnostic"])
        self.assertGreater(basin["declared_box_chart_overshoot_factor_H18"], 1.0)
        self.assertGreater(basin["max_uniform_scale_of_declared_box"], 0.0)
        self.assertGreater(basin["max_integral_radius_with_others_declared_m_s"], 0.0)
        # The volume optimum admits a larger integral radius than the
        # others-at-declared frontier point, which is the whole reason the basin
        # is carried as a polytope rather than one number.
        self.assertGreater(basin["max_volume_radii_H18"]["integral_displacement"],
                           basin["max_integral_radius_with_others_declared_m_s"])

    def test_open_P4_blockers_are_mirrored_from_the_final_gate(self):
        self.assertTrue(self.d["remaining_P4_mathematical_blockers"])
        self.assertFalse(self.d["P4_invariance_obligations"]["P4_motion_practical_ISS"]["established"])
        self.assertTrue(
            self.d["P4_invariance_obligations"]["conditional_binary32_additive_ISS"]["established"])
        self.assertTrue(
            self.d["P4_invariance_obligations"]["all_bias_families_closed"]["established"])

    def test_deployment_pass_is_validated_as_its_own_conjunction(self):
        d = dict(self.d)
        d["END_TO_END_DEPLOYMENT_PASS"] = True
        failures = E2E.validate(d)
        self.assertIn("deployment pass set while the end-to-end theorem is open", failures)
        self.assertIn("deployment pass set while deployment blockers remain", failures)

    def test_mirrored_P4_bits_cannot_disagree_with_their_blocker_list(self):
        d = dict(self.d)
        d["P4_PASS"] = True
        self.assertIn("P4 pass set while mathematical blockers remain", E2E.validate(d))
        d = dict(self.d)
        d["P4_MOTION_PASS"] = True
        failures = E2E.validate(d)
        self.assertIn("P4 motion pass set while mathematical blockers remain", failures)
        self.assertIn("mirrored P4_MOTION_PASS disagrees with its own invariance row", failures)
        d = dict(self.d)
        d["P5_MAY_START"] = True
        self.assertIn("P5 may start while P4 is false", E2E.validate(d))

    def test_theorem_invariants_are_unchanged(self):
        self.assertEqual(self.d["P3_delta"], 1e-18)
        self.assertFalse(self.d["filter_changed"])
        self.assertFalse(self.d["quality_gates_changed"])
        self.assertFalse(self.d["declared_domain_shrunk"])


if __name__ == "__main__":
    unittest.main()
