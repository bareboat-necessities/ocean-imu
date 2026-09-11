#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools" / "stability"
sys.path.insert(0, str(TOOLS))

import ou3_p4_source_uniform_estimator_event_attachment as ATTACH


class SourceUniformEstimatorEventAttachmentTest(unittest.TestCase):
    def test_relation_attachment_closes_without_promoting_p4(self) -> None:
        d = ATTACH.build()
        self.assertEqual(ATTACH.validate(d), [])
        self.assertTrue(d["source_uniform_estimator_owned_event_attachment_relation_closed"])
        self.assertTrue(d["kernel_authoritative_for_current_Riccati_slice_only"])
        self.assertTrue(d["JOINT_image_authoritative_for_next_frontend_state"])
        self.assertTrue(d["current_schedule_committed_before_post_measurement_adaptation"])
        self.assertTrue(d["every_joint_successor_retained"])
        self.assertTrue(d["all_kernel_frontend_children_share_current_post_Riccati_state"])
        self.assertTrue(d["trusted_event_local_P_H_R_imported_without_reconstruction"])
        self.assertTrue(d["exact_nonlinear_state_succession_materialized_between_literal_events"])
        self.assertTrue(d["same_finite_map_generates_state_and_Jacobian"])
        self.assertTrue(d["actual_applied_RS_owned_by_same_joint_image"])
        self.assertFalse(d["kernel_frontend_successor_used_to_select_next_theorem_state"])
        self.assertFalse(d["kernel_joint_child_frontend_equality_required"])
        self.assertFalse(d["favorable_successor_selected"])
        self.assertFalse(d["branch_ordinal_correspondence_assumed"])
        self.assertFalse(d["independent_P_H_R_K_reconstruction_used"])
        self.assertFalse(d["numeric_COMPLETE_BRMM_enumeration_used"])
        self.assertFalse(d["production_augmented_PrefixInput_assembled_here"])
        self.assertFalse(d["P4_MOTION_PASS"])
        self.assertFalse(d["P4_PASS"])
        self.assertFalse(d["P5_MAY_START"])

    def test_prior_frequency_and_physical_potential_survive_literal_attachment(self) -> None:
        from dataclasses import replace
        from unittest.mock import patch
        import ou3_p4_joint_brmm_frontend_transition as J
        import ou3_p4_brmm_same_signal_statistics as ST
        import ou3_brmm_complete_window_execution_kernel as K
        import ou3_brmm_physical_wave_source as W
        import ou3_brmm_frontend_state_step as F
        I = ATTACH.I
        zero = (I(0), I(0), I(0))
        js = J._smoke_state()
        names = ("accel_prev", "high_pass_1", "high_pass_1_prev", "high_pass_2",
                 "velocity", "elevation", "velocity_mean", "velocity_sq",
                 "elevation_mean", "elevation_sq", "raw_period_s")
        wpe = replace(js.frontend.wpe, **{k: I(0) for k in names},
                      elapsed_s=I(100), weight=I(1), log_period_s=None, usable_period=False)
        front = replace(js.frontend, wpe=wpe)
        js = replace(js, frontend=front, statistics=ST.from_shipping(wpe, front.tuner))
        wave = K.WavePrimitivePayload(W.certificate_id(W.spectral_certificate(())),
                 "p0", "p1", "one-Live-origin", zero, zero, zero, zero, zero, zero)
        sample = K.SampleCoordinates(
            gyro_measurement=K.MAHONY.Vec3(*zero), omega_body_corrected=zero,
            specific_force=K.MAHONY.Vec3(I(0), I(0), I(-9.80665)),
            f_cog_body=(I(0), I(0), I(-9.80665)), R_wb=ATTACH._identity(3),
            due_S=True, aw_floor_requested=False, wave_primitive=wave)
        branch = K.ExecutionBranch(front, ATTACH.WORD.initialize_word("H", ATTACH._identity(18)),
                  ATTACH.WORD.initialize_word("A", ATTACH._identity(21)), "prior-root")
        with patch.object(F, "advance", side_effect=AssertionError("obsolete usable-period guard called")):
            pairs = ATTACH.synchronize_sample(
                branch=branch, joint_state=js, sample=sample,
                state_in_H=[I(0)]*18, state_in_A=[I(0)]*21, radial_scale=I(0),
                true_bias=zero, bias_projection_limit=.4, tau_ba=I(1800),
                sample_index=0, next_cell_prefix="prior-test")
        self.assertTrue(pairs)
        for pair in pairs:
            for attached in pair:
                self.assertTrue(all(c.wave_primitive == wave for c in attached.cells))
                self.assertEqual(attached.image.frequency_hz, I(.2))
                bad = list(attached.cells)
                bad[0] = replace(bad[0], wave_primitive=None)
                self.assertTrue(ATTACH.BIND.validate_event_cells_against_selector(
                    attached.selector, bad, mode=attached.mode))
        # This executes a legitimate prior predecessor; it is not a finite
        # capture or complete 601-sample covariance/physical source cover.

    def test_smoke_retains_both_modes_and_literal_order(self) -> None:
        s = ATTACH._smoke()
        self.assertGreater(s["joint_successors_retained"], 0)
        self.assertTrue(s["all_H_cells_bound"])
        self.assertTrue(s["all_A_cells_bound"])
        self.assertTrue(s["all_event_orders_equal"])
        self.assertTrue(s["authoritative_next_frontend_is_joint"])


if __name__ == "__main__":
    unittest.main()
