import copy
import sys
from pathlib import Path
import unittest

TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_sea3_hard_finite_window_source as SEA0


class Sea3HardFiniteWindowSourceTest(unittest.TestCase):
    def _candidate(self):
        transitions = []
        zero3 = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        eye3 = [
            [[1.0, 1.0], [0.0, 0.0], [0.0, 0.0]],
            [[0.0, 0.0], [1.0, 1.0], [0.0, 0.0]],
            [[0.0, 0.0], [0.0, 0.0], [1.0, 1.0]],
        ]
        for k in range(SEA0.SAMPLES):
            tr = f"tr-{k}"
            resp = f"resp-{k}"
            transitions.append({
                "k": k,
                "source_identity": "sea3-root",
                "xs_in_id": f"xs-{k}",
                "xs_out_id": f"xs-{k+1}",
                "lambda_in_id": f"lambda-{k}",
                "lambda_out_id": f"lambda-{k+1}",
                "source_transition_witness_id": tr,
                "joint_response_witness_id": resp,
                "joint_physical_output": {
                    "source_transition_witness_id": tr,
                    "joint_response_witness_id": resp,
                    "gyro_measurement_interval": zero3,
                    "omega_body_corrected_interval": zero3,
                    "specific_force_body_interval": zero3,
                    "f_cog_body_interval": zero3,
                    "R_wb_interval": eye3,
                },
                "source_events": {
                    "source_transition_witness_id": tr,
                    "magnetometer_events_after_imu": [],
                    "aw_covariance_floor_requested": False,
                    "S_zero_due": False,
                },
            })
        return {
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
            "front_end_entry_witness_id": "frontend-entry",
            "live_covariance_seed_witness_id": "live-seed",
            "front_end_entry": {},
            "live_covariance_seed": {},
            "transitions": transitions,
        }

    def test_status_preserves_compact_sea3_and_fails_closed(self):
        d = SEA0.build()
        self.assertEqual(SEA0.validate_status(d), [])
        self.assertTrue(d["SEA3_parameter_domain_compact"])
        self.assertTrue(d["compact_transition_relation_is_theorem_domain"])
        ingredients = d["executable_provider_ingredients"]
        self.assertTrue(ingredients["machine_readable_R_lambda_closed"])
        self.assertFalse(ingredients["hard_shaping_state_or_excitation_bound_closed"])
        self.assertFalse(ingredients["joint_translational_rotational_shaping_closed"])
        rlambda = d["R_lambda_certificate"]
        self.assertTrue(rlambda["actual_rate_bounded_R_lambda_subset_Rhat"])
        self.assertFalse(rlambda["rate_constants_fitted_or_invented"])
        self.assertFalse(rlambda["fixed_lambda_word_used"])
        self.assertTrue(d["executor_payload_contract"]["raw_gyro_and_corrected_rate_are_distinct_coordinates"])
        self.assertFalse(d["executor_payload_contract"]["precomputed_aw_covariance_floor_increment_allowed"])
        self.assertFalse(d["provider_implementation_closed"])
        self.assertFalse(d["source_reachable_event_family_materialized"])
        self.assertFalse(d["P3_promoted"])

    def test_structurally_complete_candidate_still_cannot_self_promote(self):
        d = self._candidate()
        self.assertEqual(SEA0.validate_candidate_structure(d), [])
        failures = SEA0.validate_artifact(d)
        self.assertIn(
            "validated SEA0 hard finite-window provider is not implemented",
            failures,
        )

    def test_every_shortcut_is_rejected(self):
        base = self._candidate()
        for key in SEA0._FORBIDDEN_TRUE_FLAGS:
            with self.subTest(key=key):
                d = copy.deepcopy(base)
                d[key] = True
                failures = SEA0.validate_candidate_structure(d)
                self.assertTrue(any(key in x for x in failures), failures)

    def test_broken_xs_or_lambda_history_is_rejected(self):
        for field, expected in (
            ("xs_in_id", "x^s phase continuity"),
            ("lambda_in_id", "lambda transition continuity"),
        ):
            with self.subTest(field=field):
                d = self._candidate()
                d["transitions"][137][field] = "not-the-predecessor"
                failures = SEA0.validate_candidate_structure(d)
                self.assertTrue(any(expected in x for x in failures), failures)

    def test_detached_physical_or_event_witness_is_rejected(self):
        d = self._candidate()
        d["transitions"][9]["joint_physical_output"]["source_transition_witness_id"] = "other"
        failures = SEA0.validate_candidate_structure(d)
        self.assertTrue(any("physical output detached" in x for x in failures), failures)

        d = self._candidate()
        d["transitions"][9]["source_events"]["source_transition_witness_id"] = "other"
        failures = SEA0.validate_candidate_structure(d)
        self.assertTrue(any("source events detached" in x for x in failures), failures)

    def test_precomputed_covariance_floor_increment_is_rejected(self):
        d = self._candidate()
        d["transitions"][11]["source_events"]["aw_covariance_floor_increment"] = [[0.0]]
        failures = SEA0.validate_candidate_structure(d)
        self.assertTrue(any("illegally serializes" in x for x in failures), failures)

    def test_raw_gyro_and_corrected_rate_are_both_required(self):
        for key in ("gyro_measurement_interval", "omega_body_corrected_interval"):
            with self.subTest(key=key):
                d = self._candidate()
                del d["transitions"][3]["joint_physical_output"][key]
                failures = SEA0.validate_candidate_structure(d)
                self.assertTrue(any(key in x for x in failures), failures)

    def test_wrong_window_length_is_rejected(self):
        d = self._candidate()
        d["transitions"].pop()
        failures = SEA0.validate_candidate_structure(d)
        self.assertTrue(any("exactly 601" in x for x in failures), failures)


if __name__ == "__main__":
    unittest.main()