from pathlib import Path
import sys
import unittest
from unittest import mock

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import ou3_p4_canonical_gate as GATE
import ou3_p4_complete_word_dissipation as W

TINY = 2.0e-18


def _p3(margin_h=3.0e-18, margin_a=TINY):
    return {
        "qualification": "OU3_P3_CANONICAL_THEOREM_INTERFACE",
        "canonical_definition_frozen": True,
        "only_this_module_may_promote_P3_for_P4": True,
        "useful_gate": 1.0e-18,
        "P3_CANONICAL_PASS": True,
        "P4_MAY_CONSUME_P3": True,
        "worst_H_A_margin": min(margin_h, margin_a),
        "mode_margins": {"H": margin_h, "A": margin_a},
        "validation_pass": True,
        "validation_failures": [],
    }


def _envelope():
    return {
        "translation_covariance_upper_groups": {
            "v": 25.0, "p": 400.0, "S": 90000.0, "a_w": 8.66,
        },
        "H_bias_covariance_upper": {
            "theta_covariance_upper": 0.25,
            "gyro_bias_covariance_upper": 1.0e-4,
            "accel_bias_covariance_upper": None,
        },
        "A_bias_covariance_upper": {
            "theta_covariance_upper": 0.25,
            "gyro_bias_covariance_upper": 1.0e-4,
            "accel_bias_covariance_upper": 0.25,
        },
    }


def _metric():
    rows = [
        {
            "source_node": i,
            "boundary_history_envelope": _envelope(),
            "positive_phase_history_envelope": _envelope(),
        }
        for i in range(800)
    ]
    lower = {
        "theta": 1.0e-10, "b_g": 1.0e-14, "v": 1.0e-16,
        "p": 1.0e-20, "S": 1.0e-30, "a_w": 1.0e-12,
    }
    return {
        "qualification": "OU3_P4_CANONICAL_P3_SOURCE_PHASE_METRIC_ATTACHMENT",
        "same_history_P3_frontier_consumed": True,
        "independent_cartesian_tau_sigma_R_S_extrema_used": False,
        "finite_source_phase_classes": 800 * 26,
        "canonical_P3_candidate_numeric_pass_observed": True,
        "translation_post_acc_S_at_metric_boundary": True,
        "H_A_fresh_process_floor_at_same_endpoint_vector_packet": True,
        "magnetometer_translation_jacobian_zero_on_declared_branch": True,
        "zero_lever_arm_branch": True,
        "dormant_transparent_vibration_guard_branch": True,
        "frozen_clock": {"absorbing_hold_arbitrary_duration_covered": True},
        "endpoint_rows": rows,
        "global_source_phase_covariance_lower_group_diagonal": {
            "H": dict(lower),
            "A": dict(lower, b_a=1.0e-11),
        },
    }


def _timing():
    return {
        "qualification": "OU3_P4_SOURCE_COMPLETE_WORD_TIMING_DECOMPOSITION",
        "configured_dt_s": 0.005,
        "word_samples_upper": 635,
        "word_horizon_s": 3.1686363527551307,
        "S_firing_times_are_source_intervals_not_fixed_samples": True,
        "S_residual_exactly_linear_selector": True,
        "S_nonlinear_eta_identically_zero": True,
        "S_timing_consumed_by_linear_P3_translation_UCO": True,
        "nonlinear_timing_obligations_reduce_to_vector_measurements": True,
        "ready_for_source_complete_nonlinear_remainder_composition": True,
        "fixed_minimum_gap_S_schedule_is_source_complete": False,
    }


def _clock():
    return {
        "qualification": "OU3_P2_SAMPLE_CLOCK_COMMIT_REACHABILITY_REFINEMENT",
        "P2_SAMPLE_CLOCK_REFINEMENT_CERTIFICATE": "PASS",
        "source_graph_all_to_all": False,
        "frozen_clock_self_loop_included": True,
        "clock": {"floating_clock_stagnation_verified": True},
        "partition": {"states": 800},
        "transition_edges": 78633,
    }


def _audit():
    return {
        "qualification": "OU3_P4_SOURCE_CORRELATED_SIGNED_JOSEPH_FEASIBILITY_AUDIT",
        "same_history_P3_metric_consumed": True,
    }


def _build(p3=None, metric=None, timing=None, clock=None, audit=None):
    with mock.patch.object(W.METRIC, "validate", return_value=[]), \
         mock.patch.object(W.TIMING, "validate", return_value=[]), \
         mock.patch.object(W.CLOCK, "validate", return_value=[]), \
         mock.patch.object(W.AUDIT, "validate", return_value=[]):
        return W.build(
            _p3() if p3 is None else p3,
            _metric() if metric is None else metric,
            _timing() if timing is None else timing,
            _clock() if clock is None else clock,
            _audit() if audit is None else audit,
        )


class Ou3P4CompleteWordDissipationTests(unittest.TestCase):
    def test_complete_word_is_established_on_certified_prerequisites(self):
        d = _build()
        self.assertEqual([], W.validate(d))
        self.assertEqual([], d["unmet_theorem_obligations"])
        self.assertTrue(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"])
        self.assertFalse(d["P4_USABLE_CERTIFICATE_PROMOTED"])
        self.assertEqual(18, d["H_dimension"])
        self.assertEqual(21, d["A_dimension"])
        self.assertEqual(0.80, d["outer_angle_rad"])
        self.assertEqual(4 * 635, d["word_state_operations"])
        for mode in ("H", "A"):
            row = d["modes"][mode]
            self.assertEqual(1.0, row["metric_eigenvalue_lower_m_minus"])
            self.assertGreater(row["metric_eigenvalue_upper_m_plus"], 1.0)
            self.assertGreater(row["quadratic_defect_constant_upper"], 0.0)
            self.assertGreater(row["word_defect_gain_B_m_upper"], 0.0)
            self.assertGreater(row["inner_level_W_star_lower"], 0.0)
            self.assertTrue(row["funnel_consistent"])

    def test_strict_gap_is_carried_exactly_below_binary64_epsilon(self):
        d = _build()
        for mode, margin in (("H", 3.0e-18), ("A", TINY)):
            gap = d[f"one_minus_rho_{mode}_lower"]
            self.assertGreater(gap, 0.0)
            self.assertLessEqual(gap, margin / 2.0)
            # 1-gap is not representable below one, so the companion rounds up.
            self.assertEqual(1.0, d[f"rho_{mode}_upper"])
            self.assertTrue(d["modes"][mode]["strict_gap_below_binary64_epsilon"])
            self.assertEqual(gap, d[f"strict_dissipation_margin_{mode}_lower"])

    def test_inner_funnel_stays_inside_every_declared_operation_region(self):
        d = _build()
        for mode in ("H", "A"):
            row = d["modes"][mode]
            theta = row["inner_attitude_radius_theta_star_lower_rad"]
            self.assertGreater(theta, 0.0)
            self.assertLess(theta, 1.0)
            for gain in row["operation_defect_budget"]["accepted_correction_gain_upper"].values():
                self.assertLessEqual(gain * theta, W.CORRECTION_REGION_RADIUS)
        self.assertLessEqual(
            d["modes"]["A"]["inner_attitude_radius_theta_star_lower_rad"],
            d["active_accelerometer_bias_projection_interior_margin_mps2"],
        )

    def test_failing_canonical_p3_blocks_the_complete_word(self):
        p3 = _p3()
        p3["P3_CANONICAL_PASS"] = False
        p3["P4_MAY_CONSUME_P3"] = False
        d = _build(p3=p3)
        self.assertEqual([], W.validate(d))
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"])
        self.assertTrue(any("canonical P3" in x for x in d["unmet_theorem_obligations"]))

    def test_all_to_all_source_graph_blocks_the_complete_word(self):
        clock = _clock()
        clock["source_graph_all_to_all"] = True
        d = _build(clock=clock)
        self.assertEqual([], W.validate(d))
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"])
        self.assertTrue(
            any("sample-clock" in x for x in d["unmet_theorem_obligations"])
        )

    def test_margin_below_unchanged_gate_blocks_its_mode(self):
        d = _build(p3=_p3(margin_h=9.0e-19))
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"])
        self.assertNotIn("H", d["modes"])
        self.assertIn("A", d["modes"])
        self.assertTrue(any("below the unchanged gate" in x for x in d["failures"]))

    def test_canonical_gate_accepts_the_produced_candidate(self):
        d = _build()
        gate = GATE.build(
            _p3(),
            {
                "qualification": "OU3_P4_GLOBAL_CAYLEY_SECTOR_GEOMETRY",
                "source_generated_not_trajectory_fit": True,
                "declared_filter_entrance_covered": True,
                "chart_antipode_excluded": True,
                "usable_sector_geometry_pass": True,
                "pass": True,
                "trajectory_replay_used": False,
                "filter_changed": False,
                "full_18_21_state_Joseph_word_established_here": False,
                "signed_EKF_remainder_charged_here": False,
                "outer_angle_rad": 0.80,
                "validation_pass": True,
                "validation_failures": [],
            },
            {
                "qualification": "OU3_P4_GLOBAL_VECTOR_NONLINEAR_REMAINDER_SECTOR",
                "source_generated_not_trajectory_fit": True,
                "declared_filter_entrance_covered": True,
                "accelerometer_bias_cancels_exactly_from_eta": True,
                "penalties_are_homogeneous_quadratic_not_affine_beta": True,
                "measurement_covariance_isotropy_required": True,
                "usable_sector_remainder_primitive_pass": True,
                "pass": True,
                "trajectory_replay_used": False,
                "filter_changed": False,
                "complete_Joseph_word_established_here": False,
                "P4_USABLE_CERTIFICATE_PROMOTED": False,
                "outer_angle_rad": 0.80,
                "validation_pass": True,
                "validation_failures": [],
            },
            dict(
                _timing(),
                source_generated_not_trajectory_fit=True,
                trajectory_replay_used=False,
                filter_changed=False,
                old_terminal_192_201_cluster_required_for_promotion=False,
                P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE=False,
                validation_pass=True,
                validation_failures=[],
            ),
            {
                "path_graph_ready": True,
                "P2_SOURCE_PATH_CERTIFICATE": "PASS",
                "usable_P4_promoted": False,
                "partition": {"states": 800},
                "validation_pass": True,
                "validation_failures": [],
            },
            {
                "qualification": "OU3_P2_SOURCE_NODE_CELL_MATERIALIZATION",
                "source_only": True,
                "state_order_matches_P2_nested_loops": True,
                "trajectory_replay_used": False,
                "filter_changed": False,
                "source_graph_rebuilt_or_pruned_here": False,
                "P4_metric_attached_here": False,
                "P4_USABLE_CERTIFICATE_PROMOTED": False,
                "partition": {"states": 800},
                "validation_pass": True,
                "validation_failures": [],
            },
            dict(
                _clock(),
                source_only=True,
                same_physical_partition_as_P2=True,
                EMA_updated_every_valid_sample=True,
                EMA_composed_sample_by_sample=True,
                sample_varying_target_and_horizon_boxes_admitted=True,
                commit_only_stages_current_smoothed_candidate=True,
                pending_candidate_applied_before_next_sample=True,
                arbitrary_late_commit_jump_removed=True,
                trajectory_replay_used=False,
                filter_changed=False,
                P4_USABLE_CERTIFICATE_PROMOTED=False,
                finite_stage_gap_lower_samples=13,
                finite_stage_gap_upper_samples=26,
                validation_pass=True,
                validation_failures=[],
            ),
            d,
        )
        self.assertEqual([], GATE.validate(gate))
        self.assertEqual([], gate["P4_CANONICAL_FAIL_REASONS"])
        self.assertTrue(gate["P4_CANONICAL_PASS"])
        self.assertTrue(gate["P5_MAY_CONSUME_P4"])
        self.assertFalse(gate["P5_FINITE_CAPTURE_ESTABLISHED_HERE"])
        for mode in ("H", "A"):
            self.assertGreater(gate["one_minus_rho_lower"][mode], 0.0)


if __name__ == "__main__":
    unittest.main()
