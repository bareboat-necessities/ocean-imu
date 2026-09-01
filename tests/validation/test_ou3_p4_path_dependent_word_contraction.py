from copy import deepcopy
from decimal import Decimal
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_path_dependent_word_contraction as ROUTE
from ou3_proof_module_state import preserve_module_bindings


class PathDependentWordContractionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with preserve_module_bindings():
            cls.d = ROUTE.build()
        cls.assertions = ROUTE.validate(cls.d)

    def _good_candidate(self):
        route = self.d
        g = route["source_graph"]
        return {
            "qualification": ROUTE.JACOBIAN_QUALIFICATION,
            "source_only": True,
            "trajectory_replay_used": False,
            "outward_validated": True,
            "exact_shipping_complete_word_map": True,
            "per_operation_sector_invariance_required": False,
            "N_times_global_defect_used": False,
            "P3_delta_used_as_nonlinear_radius": False,
            "proof_operating_domain_sha256": route["proof_operating_domain_sha256"],
            "certified_outer_angle_rad": route["outer_geometry_angle_rad"],
            "certified_outer_cayley_norm_upper": route["outer_geometry_cayley_norm_upper"],
            "full_state_domain_exact_match": True,
            "source_partition_state_count": g["partition_state_count"],
            "source_transition_edge_count": g["transition_edge_count"],
            "recurrent_state_count": g["recurrent_state_count"],
            "all_required_reachable_word_edges_checked": True,
            "sequential_reset_jacobians_included": True,
            "accepted_rejected_not_due_branches_covered": True,
            "modes": {
                "H": {
                    "dimension": 18,
                    "endpoint_metric_uses_actual_reachable_source_node": True,
                    "full_state_cross_terms_retained": True,
                    "max_whitened_generalized_jacobian_norm_upper": 0.97,
                },
                "A": {
                    "dimension": 21,
                    "endpoint_metric_uses_actual_reachable_source_node": True,
                    "full_state_cross_terms_retained": True,
                    "max_whitened_generalized_jacobian_norm_upper": 0.98,
                    "accelerometer_bias_projection_generalized_jacobian_included": True,
                },
            },
        }

    def test_route_validates_but_does_not_promote_p4(self):
        d = self.d
        self.assertEqual(self.assertions, [])
        self.assertTrue(d["route_contract_pass"])
        self.assertFalse(d["P4_COMPLETE_WORD_DIFFERENTIAL_CONTRACTION_ESTABLISHED_HERE"])
        self.assertFalse(d["P4_USABLE_CERTIFICATE_PROMOTED"])
        self.assertFalse(d["P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE"])

    def test_independent_route_preserves_45deg_and_exact_0p80rad_geometry(self):
        d = self.d
        self.assertEqual(d["P5_entrance_angle_deg"], 45.0)
        self.assertTrue(d["P5_45deg_entrance_preserved"])
        self.assertEqual(d["outer_geometry_angle_rad"], 0.80)
        self.assertLess(d["outer_geometry_cayley_norm_upper"], 1.0)
        self.assertFalse(d["declared_domain_changed"])
        self.assertFalse(d["filter_changed"])
        self.assertFalse(d["source_replay_used"])

    def test_candidate_is_bound_to_exact_full_state_and_outer_domain(self):
        d = self.d
        self.assertEqual(len(d["proof_operating_domain_sha256"]), 64)
        c = d["candidate_domain_contract"]
        self.assertTrue(c["proof_operating_domain_hash_exact_match_required"])
        self.assertTrue(c["full_state_domain_exact_match_required"])
        self.assertEqual(c["certified_outer_angle_rad_required"], 0.80)
        self.assertEqual(c["certified_outer_cayley_norm_must_cover"], d["outer_geometry_cayley_norm_upper"])

    def test_complete_word_is_atomic_and_old_scalar_routes_are_not_required(self):
        r = self.d["route_distinctions"]
        self.assertEqual(r["atomic_object"], "COMPLETE_SOURCE_WORD_RETURN_MAP")
        self.assertTrue(r["individual_packet_rank_may_be_singular"])
        self.assertFalse(r["per_packet_full_rank_required"])
        self.assertFalse(r["per_operation_sector_invariance_required"])
        self.assertFalse(r["N_times_global_lipschitz_defect_used"])
        self.assertFalse(r["translation_nontranslation_schur_split_required"])
        self.assertTrue(r["full_18_21_state_jacobian_test_required"])
        self.assertFalse(r["P3_delta_used_as_nonlinear_radius"])

    def test_metrics_are_source_node_correlated_without_arbitrary_rescaling(self):
        d = self.d
        for mode, dim in (("H", 18), ("A", 21)):
            m = d["modes"][mode]
            self.assertEqual(m["dimension"], dim)
            self.assertTrue(m["same_mode_global_scale_on_all_nodes"])
            self.assertTrue(m["full_attitude_linear_cross_terms_retained"])
            self.assertFalse(m["block_diagonal_metric_used"])
            self.assertFalse(m["translation_nontranslation_schur_split_required"])
            self.assertFalse(m["per_packet_full_rank_required"])
            self.assertFalse(m["per_operation_contraction_required"])
            self.assertIn("J_w(x)^T M_h J_w(x)", m["acceptance_matrix_inequality"])

    def test_tiny_p3_delta_is_visible_but_not_used_as_nonlinear_budget(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            delta = Decimal(m["P3_linear_origin_contraction_gap_diagnostic"]["delta"])
            gap = Decimal(m["P3_linear_origin_contraction_gap_diagnostic"]["one_minus_sqrt_one_minus_delta"])
            self.assertGreater(delta, Decimal(0))
            self.assertGreater(gap, Decimal(0))
            self.assertLess(gap, Decimal("1e-15"))
            self.assertFalse(m["P3_delta_used_as_nonlinear_radius"])
            self.assertFalse(m["P3_delta_used_as_whole_word_jacobian_acceptance_threshold"])

    def test_source_graph_scope_is_nonempty_and_source_faithful(self):
        g = self.d["source_graph"]
        self.assertTrue(g["path_graph_ready"])
        self.assertGreater(g["partition_state_count"], 0)
        self.assertGreater(g["transition_edge_count"], 0)
        self.assertGreater(g["recurrent_state_count"], 0)
        self.assertTrue(g["raw_tuner_sigma_subfloor_states_included"])
        self.assertFalse(g["RS_target_powf_tightening_used"])

    def test_future_jacobian_contract_accepts_only_strict_full_state_source_result(self):
        good = self._good_candidate()
        status = ROUTE._candidate_status(good, self.d)
        self.assertTrue(status["contract_accepted"])
        self.assertTrue(status["all_modes_strict"])
        self.assertEqual(status["reasons"], [])

        bad = deepcopy(good)
        bad["modes"]["H"]["max_whitened_generalized_jacobian_norm_upper"] = 1.0
        status = ROUTE._candidate_status(bad, self.d)
        self.assertFalse(status["contract_accepted"])
        self.assertTrue(any("H:" in x for x in status["reasons"]))

    def test_future_candidate_on_smaller_or_different_domain_is_rejected(self):
        bad = self._good_candidate()
        bad["certified_outer_angle_rad"] = 0.70
        status = ROUTE._candidate_status(bad, self.d)
        self.assertFalse(status["contract_accepted"])
        self.assertTrue(any("outer-angle" in x for x in status["reasons"]))

        bad = self._good_candidate()
        bad["proof_operating_domain_sha256"] = "0" * 64
        status = ROUTE._candidate_status(bad, self.d)
        self.assertFalse(status["contract_accepted"])
        self.assertTrue(any("domain hash" in x for x in status["reasons"]))

        bad = self._good_candidate()
        bad["full_state_domain_exact_match"] = False
        status = ROUTE._candidate_status(bad, self.d)
        self.assertFalse(status["contract_accepted"])
        self.assertTrue(any("full-state theorem-domain" in x for x in status["reasons"]))

    def test_boolean_gamma_is_rejected_not_interpreted_as_zero(self):
        bad = self._good_candidate()
        bad["modes"]["H"]["max_whitened_generalized_jacobian_norm_upper"] = False
        status = ROUTE._candidate_status(bad, self.d)
        self.assertFalse(status["contract_accepted"])
        self.assertTrue(any("H:" in x for x in status["reasons"]))

    def test_candidate_metadata_can_never_promote_p4_in_route_contract(self):
        self.assertFalse(self.d["P4_USABLE_CERTIFICATE_PROMOTED"])
        self.assertFalse(self.d["jacobian_candidate"]["provided"])

    def test_A_mode_projection_must_be_inside_generalized_jacobian(self):
        self.assertTrue(self.d["modes"]["A"]["A_mode_projection_generalized_jacobian_required"])
        self.assertTrue(self.d["finite_angle_acceptance"]["A_mode_bias_projection_requires_generalized_jacobian"])

    def test_outer_geometry_mutation_is_rejected(self):
        d = deepcopy(self.d)
        d["outer_geometry_angle_rad"] = 0.81
        failures = ROUTE.validate(d)
        self.assertTrue(any("exactly 0.80" in x for x in failures))

    def test_next_obligation_is_direct_full_word_jacobian_not_packet_budget(self):
        text = self.d["next_obligation"]
        self.assertIn("18/21-state", text)
        self.assertIn("generalized-Jacobian", text)
        self.assertIn("below one", text)
        self.assertIn("do not scalarize packet defects", text)


if __name__ == "__main__":
    unittest.main()
