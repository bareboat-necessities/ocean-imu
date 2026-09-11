import unittest
from tools.stability.ou3_alt_contraction import proof_plan as P

class PlanGuardTests(unittest.TestCase):
    def test_diagnostic_cannot_be_main_proof_task(self):
        with self.assertRaises(RuntimeError):
            P.require_theorem_task(obligation='Live word',evidence_kind='unreachable_perturbation',complete_physical_word=False,requested_phase='physical_map')
    def test_storage_is_hard_blocked_before_complete_word(self):
        with self.assertRaises(RuntimeError):
            P.require_theorem_task(obligation='common storage',evidence_kind='analytic_source',complete_physical_word=False,requested_phase='storage_search')
        with self.assertRaises(RuntimeError):P.assert_storage_search_allowed({})
    def test_six_assembly_flags_do_not_qualify_a_finite_master(self):
        s={k:True for k in ('same_history_complete_BRMM_word','physical_prediction_forcing_attached','physical_S_residual_attached','all_bias_families_attached','all_literal_branches_attached','H18_A21_edge_attached')}
        with self.assertRaisesRegex(RuntimeError,'finite-state storage blocked'):
            P.assert_finite_storage_master(s)
    def test_derivative_word_is_not_a_finite_state_identity(self):
        from tools.stability.ou3_alt_contraction import physical_word as W
        with self.assertRaisesRegex(RuntimeError,'not a Jacobian cocycle'):
            P.assert_finite_storage_master(W.finite_storage_readiness())
    def test_finite_guard_requires_every_branch_and_reference_graph(self):
        s=dict(map_representation='finite_physical_descriptor',finite_error_identity_for_every_event=True,
               physical_reference_forcing_retained=True,all_coefficient_product_graphs_retained=True,
               all_configured_branches_bound_to_finite_graph=True)
        P.assert_finite_storage_master(s)  # Process check only, never a proof certificate.
        for key in tuple(s):
            q=dict(s);q.pop(key)
            with self.assertRaises(RuntimeError):P.assert_finite_storage_master(q)
    def test_storage_guard_names_every_required_closure(self):
        s={k:True for k in ('same_history_complete_BRMM_word','physical_prediction_forcing_attached','physical_S_residual_attached','all_bias_families_attached','all_literal_branches_attached','H18_A21_edge_attached')}
        self.assertTrue(P.storage_search_allowed(s));P.assert_storage_search_allowed(s)

if __name__=='__main__':unittest.main()
