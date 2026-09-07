#!/usr/bin/env python3
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

import ou3_p4_complete_sea3_reset_gauge_attachment as ATTACH


class CompleteSea3ResetGaugeAttachmentTests(unittest.TestCase):
    def test_frozen_p3_attaches_to_homogeneous_physical_zero_reset_gauge(self):
        d = ATTACH.build()
        self.assertEqual(ATTACH.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertEqual(d["P3_delta_consumed"], 1e-18)
        self.assertTrue(d["P3_frozen_not_modified"])
        self.assertTrue(d["P3_reset_margin_preserved_for_every_finite_injection"])
        self.assertTrue(d["P3_reset_congruence_metric_isometry"])
        self.assertTrue(d["arbitrary_finite_reset_energy_identity_enclosed"])
        self.assertTrue(d["homogeneous_zero_error_reset_G_is_identity"])
        self.assertTrue(d["homogeneous_physical_events_fix_zero"])
        self.assertTrue(d["homogeneous_physical_zero_error_tangent_is_I_minus_KH"])
        self.assertTrue(d["P3_margin_valid_in_zero_reset_congruent_representative"])
        self.assertTrue(d["zero_error_P3_to_physical_P4_tangent_attachment_closed"])

    def test_bridge_does_not_promote_finite_reset_or_forcing_obligations(self):
        d = ATTACH.build()
        self.assertFalse(d["nonzero_replay_reset_inserted_into_homogeneous_tangent"])
        self.assertFalse(d["reset_injection_added_as_new_SEA3_source_coordinate"])
        self.assertFalse(d["finite_error_nonlinear_reset_transport_closed_here"])
        self.assertFalse(d["nonzero_forcing_noise_reset_attachment_closed_here"])
        self.assertFalse(d["filter_changed"])
        self.assertFalse(d["declared_domain_changed"])
        self.assertFalse(d["source_family_replaced"])
        self.assertFalse(d["trajectory_replay_used"])
        self.assertFalse(d["P4_promoted_here"])
        self.assertFalse(d["P5_may_start_here"])


if __name__ == "__main__":
    unittest.main()
