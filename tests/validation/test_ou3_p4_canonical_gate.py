from pathlib import Path
import sys
import unittest

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import ou3_p4_canonical_gate as G


def validated(payload):
    out = dict(payload)
    out["validation_pass"] = True
    out["validation_failures"] = []
    return out


def fixtures():
    p3 = validated({
        "qualification": "OU3_P3_CANONICAL_THEOREM_INTERFACE",
        "canonical_definition_frozen": True,
        "only_this_module_may_promote_P3_for_P4": True,
        "useful_gate": 1.0e-18,
        "P3_CANONICAL_PASS": True,
        "P4_MAY_CONSUME_P3": True,
        "worst_H_A_margin": 2.0e-18,
        "mode_margins": {"H": 3.0e-18, "A": 2.0e-18},
    })
    cayley = validated({
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
    })
    remainder = validated({
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
    })
    timing = validated({
        "qualification": "OU3_P4_SOURCE_COMPLETE_WORD_TIMING_DECOMPOSITION",
        "source_generated_not_trajectory_fit": True,
        "S_firing_times_are_source_intervals_not_fixed_samples": True,
        "S_residual_exactly_linear_selector": True,
        "S_nonlinear_eta_identically_zero": True,
        "S_timing_consumed_by_linear_P3_translation_UCO": True,
        "nonlinear_timing_obligations_reduce_to_vector_measurements": True,
        "ready_for_source_complete_nonlinear_remainder_composition": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "fixed_minimum_gap_S_schedule_is_source_complete": False,
        "old_terminal_192_201_cluster_required_for_promotion": False,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "word_horizon_s": 3.1686363527551307,
    })
    path = validated({
        "path_graph_ready": True,
        "P2_SOURCE_PATH_CERTIFICATE": "PASS",
        "usable_P4_promoted": False,
        "partition": {"states": 800},
    })
    nodes = validated({
        "qualification": "OU3_P2_SOURCE_NODE_CELL_MATERIALIZATION",
        "source_only": True,
        "state_order_matches_P2_nested_loops": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "source_graph_rebuilt_or_pruned_here": False,
        "P4_metric_attached_here": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "partition": {"states": 800},
    })
    clock = validated({
        "qualification": "OU3_P2_SAMPLE_CLOCK_COMMIT_REACHABILITY_REFINEMENT",
        "P2_SAMPLE_CLOCK_REFINEMENT_CERTIFICATE": "PASS",
        "source_only": True,
        "same_physical_partition_as_P2": True,
        "EMA_updated_every_valid_sample": True,
        "EMA_composed_sample_by_sample": True,
        "sample_varying_target_and_horizon_boxes_admitted": True,
        "commit_only_stages_current_smoothed_candidate": True,
        "pending_candidate_applied_before_next_sample": True,
        "arbitrary_late_commit_jump_removed": True,
        "frozen_clock_self_loop_included": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "source_graph_all_to_all": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "partition": {"states": 800},
        "clock": {"floating_clock_stagnation_verified": True},
        "finite_stage_gap_lower_samples": 13,
        "finite_stage_gap_upper_samples": 26,
        "transition_edges": 12345,
    })
    candidate = {
        "qualification": "OU3_P4_COMPLETE_NONLINEAR_WORD_DISSIPATION_V1",
        "source_generated_not_trajectory_fit": True,
        "canonical_P3_artifact_consumed": True,
        "same_source_history_for_metric_and_nonlinear_word": True,
        "implemented_prediction_measurement_Joseph_order_covered": True,
        "source_complete_vector_packet_language_covered": True,
        "S_linear_timing_discharged_by_canonical_P3": True,
        "Cayley_exact_geometry_consumed": True,
        "homogeneous_vector_remainder_consumed": True,
        "finite_speed_sample_clock_graph_consumed": True,
        "frozen_clock_absorbing_hold_branch_covered": True,
        "zero_lever_arm_branch": True,
        "dormant_transparent_vibration_guard_branch": True,
        "full_H18_state_word_covered": True,
        "full_A21_state_word_covered": True,
        "signed_nonlinear_remainder_charged": True,
        "complete_word_generalized_Jacobian_or_equivalent_bound": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "H_dimension": 18,
        "A_dimension": 21,
        "outer_angle_rad": 0.80,
        "canonical_P3_worst_H_A_margin_consumed": 2.0e-18,
        "source_word_horizon_s": 3.1686363527551307,
        "sample_clock_transition_edges": 12345,
        "rho_H_upper": 0.92,
        "rho_A_upper": 0.95,
        "strict_dissipation_margin_H_lower": 0.07,
        "strict_dissipation_margin_A_lower": 0.04,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED": True,
    }
    return p3, cayley, remainder, timing, path, nodes, clock, candidate


def build_with(candidate_marker="fixture"):
    p3, cayley, remainder, timing, path, nodes, clock, candidate = fixtures()
    if candidate_marker is None:
        candidate = None
    return G.build(p3, cayley, remainder, timing, path, nodes, clock, candidate)


class Ou3P4CanonicalGateTests(unittest.TestCase):
    def test_synthetic_complete_word_candidate_can_pass_only_frozen_interface(self):
        d = build_with()
        self.assertEqual(G.validate(d), [])
        self.assertTrue(d["P4_CANONICAL_PASS"])
        self.assertTrue(d["P5_MAY_CONSUME_P4"])
        self.assertFalse(d["P5_FINITE_CAPTURE_ESTABLISHED_HERE"])
        self.assertEqual(d["rho_upper"], {"H": 0.92, "A": 0.95})
        self.assertEqual(d["required_dimensions"], {"H": 18, "A": 21})
        self.assertEqual(d["required_outer_angle_rad"], 0.80)
        self.assertEqual(d["P3_useful_gate"], 1.0e-18)

    def test_sector_and_source_primitives_cannot_substitute_for_complete_word(self):
        d = build_with(None)
        self.assertEqual(G.validate(d), [])
        self.assertFalse(d["P4_CANONICAL_PASS"])
        self.assertFalse(d["P5_MAY_CONSUME_P4"])
        self.assertTrue(any("candidate is missing" in x for x in d["P4_CANONICAL_FAIL_REASONS"]))

    def test_canonical_p3_must_actually_pass_unchanged_gate(self):
        p3, cayley, remainder, timing, path, nodes, clock, candidate = fixtures()
        p3["P3_CANONICAL_PASS"] = False
        p3["P4_MAY_CONSUME_P3"] = False
        p3["worst_H_A_margin"] = 9.0e-19
        p3["mode_margins"]["A"] = 9.0e-19
        candidate["canonical_P3_worst_H_A_margin_consumed"] = 9.0e-19
        d = G.build(p3, cayley, remainder, timing, path, nodes, clock, candidate)
        self.assertEqual(G.validate(d), [])
        self.assertFalse(d["P4_CANONICAL_PASS"])
        reasons = "\n".join(d["P4_CANONICAL_FAIL_REASONS"])
        self.assertIn("canonical P3 has not passed", reasons)
        self.assertIn("worst H/A margin is missing or below 1e-18", reasons)

    def test_rho_must_be_strictly_below_one_for_both_modes(self):
        p3, cayley, remainder, timing, path, nodes, clock, candidate = fixtures()
        candidate["rho_A_upper"] = 1.0
        d = G.build(p3, cayley, remainder, timing, path, nodes, clock, candidate)
        self.assertEqual(G.validate(d), [])
        self.assertFalse(d["P4_CANONICAL_PASS"])
        self.assertTrue(any("rho_A_upper" in x for x in d["P4_CANONICAL_FAIL_REASONS"]))

    def test_candidate_must_bind_exact_source_graph_and_p3_metric(self):
        p3, cayley, remainder, timing, path, nodes, clock, candidate = fixtures()
        candidate["sample_clock_transition_edges"] += 1
        candidate["canonical_P3_worst_H_A_margin_consumed"] *= 2.0
        d = G.build(p3, cayley, remainder, timing, path, nodes, clock, candidate)
        self.assertEqual(G.validate(d), [])
        self.assertFalse(d["P4_CANONICAL_PASS"])
        reasons = "\n".join(d["P4_CANONICAL_FAIL_REASONS"])
        self.assertIn("sample-clock edge family differs", reasons)
        self.assertIn("exact canonical P3 metric margin", reasons)

    def test_complete_word_cannot_silently_drop_active_accel_bias_mode(self):
        p3, cayley, remainder, timing, path, nodes, clock, candidate = fixtures()
        candidate["full_A21_state_word_covered"] = False
        candidate["A_dimension"] = 18
        d = G.build(p3, cayley, remainder, timing, path, nodes, clock, candidate)
        self.assertEqual(G.validate(d), [])
        self.assertFalse(d["P4_CANONICAL_PASS"])
        reasons = "\n".join(d["P4_CANONICAL_FAIL_REASONS"])
        self.assertIn("full_A21_state_word_covered", reasons)
        self.assertIn("A dimension is not 21", reasons)


if __name__ == "__main__":
    unittest.main()
