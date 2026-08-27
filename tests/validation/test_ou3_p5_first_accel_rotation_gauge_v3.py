import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_first_accel_rotation_gauge_v3 as G

SPECTRAL_FALLBACK_OBLIGATION = (
    "rotation-gauged innovation still required loose spectral inverse fallback"
)


class Ou3P5FirstAccelRotationGaugeV3Tests(unittest.TestCase):
    """V3 is the rotation-gauged stage that runs to a result.

    It certifies the first-prefix source sparsity before the generic interval
    product, so the subnormal dust that stopped V2 no longer masquerades as
    reachable cross-axis covariance.  The stage still does not close: part of
    the source cover leaves the innovation box too wide for a fixed-pivot
    inverse and falls back on the deliberately loose S>=R enclosure.  That is
    reported as the open subdivision target, which is what the structured-gain
    and sample-1 stages downstream of here refine.
    """

    @classmethod
    def setUpClass(cls):
        cls.d = G.build(
            source_pieces=2,
            yaw_axis_face_pieces=4,
            force_magnitude_pieces=4,
        )

    def test_structural_zeros_are_source_certified_not_thresholded(self):
        d = self.d
        self.assertTrue(d["first_prefix_source_sparsity_certified_before_interval_product"])
        self.assertGreater(d["source_sparsity_cells_checked"], 0)
        self.assertFalse(d["structural_zero_canonicalization_uses_numeric_threshold"])
        self.assertFalse(d["cross_axis_interval_dust_treated_as_physical_covariance"])
        self.assertTrue(d["axis_equivalent_same_axis_intervals_hulled"])
        self.assertGreaterEqual(d["arithmetic_zero_dust_abs_upper_max"], 0.0)

    def test_certified_dust_stays_far_below_any_reachable_covariance(self):
        self.assertLess(self.d["arithmetic_zero_dust_abs_upper_max"], 1e-300)

    def test_source_and_deployed_range_contracts_remain_unchanged(self):
        d = self.d
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])
        self.assertTrue(d["rotation_gauge_sets_J_aw_to_identity"])
        self.assertTrue(d["rotation_gauge_sets_specific_force_direction_to_e3"])

    def test_source_sparsity_certificate_admits_the_whole_child_cover(self):
        d = self.d
        self.assertGreater(d["evaluated_child_count"], 0)
        self.assertEqual(
            d["evaluated_child_count"],
            d["source_phase_cell_count"] * d["yaw_axis_cell_count"] * d["force_magnitude_cell_count"],
        )
        self.assertEqual(
            d["fixed_pivot_inverse_count"] + d["spectral_fallback_inverse_count"],
            d["evaluated_child_count"],
        )

    def test_spectral_fallback_is_the_only_outstanding_stage_obligation(self):
        d = self.d
        self.assertGreater(d["fixed_pivot_inverse_count"], 0)
        self.assertGreater(d["spectral_fallback_inverse_count"], 0)
        self.assertEqual(G.validate(d), [SPECTRAL_FALLBACK_OBLIGATION])

    def test_stage_is_fail_closed_on_the_children_it_cannot_bound(self):
        d = self.d
        self.assertEqual(d["P5_FIRST_ACCEL_ROTATION_GAUGED_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertFalse(d["all_first_accelerometer_children_inside_validated_correction_range"])
        self.assertGreater(d["children_above_validated_correction_limit"], 0)
        self.assertGreater(d["max_first_accelerometer_correction_norm_upper_rad"], 6.0)
        self.assertEqual(
            d["next_obligation"],
            "REFINE_FIRST_ACCEL_ATTITUDE_COVARIANCE_AND_EFFECTIVE_AW_DIRECTION_COUPLING",
        )

    def test_nonclosure_carries_a_source_child_witness(self):
        witness = self.d["first_unclosed_child"]
        self.assertIsNotNone(witness)
        self.assertEqual(witness["inverse_backend"], "SPD_S_GE_R_SPECTRAL_ENTRY_ENCLOSURE")
        self.assertGreater(witness["correction_norm_upper_rad"], 6.0)
        for key in ("tau_s", "sigma_aw_mps2", "R_S_filter_std", "pseudo_period_s"):
            self.assertEqual(len(witness[key]), 2)

    def test_stage_never_promotes_complete_word(self):
        d = self.d
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])


if __name__ == "__main__":
    unittest.main()
