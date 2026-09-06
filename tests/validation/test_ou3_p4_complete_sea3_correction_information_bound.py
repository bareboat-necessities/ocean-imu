"""Serialization guards; operation arithmetic is tested separately.

These tests do not run or replace the canonical P3 certificate build.
"""
from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "stability"))
import ou3_p4_complete_sea3_correction_information_bound as CORR


class CompleteSea3CorrectionInformationBoundTests(unittest.TestCase):
    def setUp(self):
        # Only serializer fixtures. Numeric candidate energy radii must be null.
        candidates = [{"attitude_angle_deg": a, "cayley_norm_upper": q}
                      for a, q in ((30.0, 0.54), (25.0, 0.45), (20.0, 0.36), (15.0, 0.27))]
        self.d = CORR._status(candidates, True)

    def test_open_bound_status_validates_without_numerical_promotion(self):
        self.assertEqual([], CORR.validate(self.d))
        self.assertFalse(self.d['reset_transport_correction_radius_source_closed'])
        self.assertFalse(self.d['candidate_metric_energy_balls_derived'])
        self.assertFalse(self.d['P4_promoted_here'])
        self.assertEqual(3, len(self.d['open_obligations']))

    def test_same_source_and_information_identity_are_retained(self):
        self.assertEqual('COMPLETE_SEA3_NORMAL_LIVE_WORD', self.d['canonical_source'])
        self.assertTrue(self.d['P3_frozen_not_modified'])
        self.assertTrue(self.d['all_due_S_updates_and_actual_RS_remain_in_complete_word'])
        self.assertTrue(self.d['all_valid_accelerometer_updates_remain_in_complete_word'])
        self.assertTrue(self.d['same_operation_correction_information_identity_valid'])
        self.assertIn('(Q_aw-I)*delta_a_w', self.d['storage_vector'])
        self.assertTrue(self.d['full_posterior_inverse_required_for_defect_cost'])

    def test_both_modes_keep_four_candidates_without_fake_balls(self):
        self.assertEqual([30.0, 25.0, 20.0, 15.0], self.d['candidate_angles_deg'])
        for mode, n in (('H', 18), ('A', 21)):
            self.assertEqual(n, self.d['modes'][mode]['dimension'])
            for row in self.d['modes'][mode]['candidate_cells']:
                self.assertIsNone(row['derived_metric_energy_radius_upper'])
                self.assertFalse(row['candidate_metric_energy_ball_certified'])

    def test_tampered_radius_or_closed_flag_is_rejected(self):
        for flag in ('candidate_metric_energy_balls_derived', 'reset_transport_correction_radius_source_closed',
                     'storage_is_original_physical_metric_isometry', 'P4_promoted_here'):
            bad = copy.deepcopy(self.d)
            bad[flag] = True
            self.assertTrue(CORR.validate(bad))
        bad = copy.deepcopy(self.d)
        bad['modes']['A']['candidate_cells'][0]['derived_metric_energy_radius_upper'] = 1e-20
        self.assertTrue(CORR.validate(bad))

    def test_conditional_p3_still_required(self):
        self.d['P3_conditional_complete_SEA3_consumed'] = False
        self.assertTrue(CORR.validate(self.d))


if __name__ == '__main__':
    unittest.main()
