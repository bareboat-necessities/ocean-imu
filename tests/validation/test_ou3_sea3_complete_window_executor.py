import sys
from pathlib import Path
import unittest

TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_sea3_complete_window_executor as EXEC
import ou3_sea3_hard_finite_window_source as SEA0


class Sea3CompleteWindowExecutorTest(unittest.TestCase):
    def test_status_is_wired_but_fail_closed_before_sea0_window(self):
        d = EXEC.build_status()
        self.assertEqual(EXEC.validate_status(d), [])
        self.assertTrue(d["only_canonical_SEA0_provider_artifact_accepted"])
        self.assertTrue(d["literal_word_shipping_parity_pass"])
        self.assertTrue(d["frontend_shipping_parity_pass"])
        self.assertTrue(d["prediction_primitives_ready"])
        self.assertFalse(d["raw_sample_array_accepted"])
        self.assertFalse(d["independent_F_Q_schedule_accepted"])
        self.assertFalse(d["independent_RS_schedule_accepted"])
        self.assertFalse(d["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"])
        self.assertFalse(d["FULL_H18_WORD_EXECUTED"])
        self.assertFalse(d["FULL_A21_WORD_EXECUTED"])
        self.assertFalse(d["FULL_H18_A21_LDLT_CLOSED"])
        self.assertFalse(d["P3_PROMOTED"])
        self.assertEqual(d["execution_contract"]["delta_lower_required"], 1.0e-18)

    def test_raw_or_self_declared_artifact_cannot_execute(self):
        candidate = {
            "schema": SEA0.SCHEMA,
            "qualification": SEA0.QUALIFICATION,
            "canonical_source": SEA0.CANONICAL_SOURCE,
            "window_horizon_s": SEA0.HORIZON_S,
            "sample_period_s": SEA0.DT_S,
            "complete_window_samples": SEA0.SAMPLES,
            "SEA3_parameter_domain_compact": True,
            "compact_transition_relation_is_theorem_domain": True,
            "phase_continuous": True,
            "same_xs_lambda_history_for_all_channels": True,
            "same_realization_drives_translation_rotation_frontend_tuner_geometry": True,
            "all_valid_accelerometer_samples_retained": True,
            "all_due_S_updates_retained": True,
            "actual_applied_per_axis_RS_retained": True,
            "asynchronous_vector_PE_events_retained": True,
            "covariance_floor_events_retained": True,
            "trajectory_replay_used": False,
            "gaussian_good_event_used": False,
            "spectral_moment_only_source_used": False,
            "arbitrary_bounded_input_source_used": False,
            "fixed_lambda_word_used": False,
            "independent_axis_boxes_used": False,
            "independent_SEA_RAO_cartesian_product_used": False,
            "finite_RAO_grid_used": False,
            "independent_tuner_schedule_used": False,
            "retired_P2_graph_used": False,
            "selected_four_S_word_used": False,
            "finite_window_representation": "equivalent_hard_finite_window_constraint",
            "provider_generated_source_family": True,
            "front_end_entry_witness_id": "entry",
            "live_covariance_seed_witness_id": "seed",
            "front_end_entry": {},
            "live_covariance_seed": {},
            "transitions": [],
        }
        failures = EXEC.validate_window_artifact(candidate)
        self.assertTrue(failures)
        self.assertIn(
            "validated SEA0 hard finite-window provider is not implemented",
            failures,
        )
        with self.assertRaises(ValueError):
            EXEC.execute_verified_window(candidate)


if __name__ == "__main__":
    unittest.main()
