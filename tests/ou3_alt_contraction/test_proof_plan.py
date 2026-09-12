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
    def test_old_assembly_flags_do_not_qualify_a_finite_master(self):
        s={k:True for k in ('same_history_complete_BRMM_word','physical_prediction_forcing_attached','physical_S_residual_attached','all_bias_families_attached','all_literal_branches_attached','H18_A21_edge_attached')}
        with self.assertRaisesRegex(RuntimeError,'finite-state storage blocked'):
            P.assert_finite_storage_master(s)
    def test_derivative_word_is_not_a_finite_state_identity(self):
        from tools.stability.ou3_alt_contraction import physical_word as W
        with self.assertRaisesRegex(RuntimeError,'not a Jacobian cocycle'):
            P.assert_finite_storage_master(W.finite_storage_readiness())
    def test_legacy_phase1_ledger_cannot_unlock_storage(self):
        from tools.stability.ou3_alt_contraction import phase1_closure as L
        d=L.build()
        self.assertTrue(d['legacy_source_ancestry_ledger_closed'])
        self.assertFalse(d['finite_storage_master_closed'])
        self.assertFalse(d['storage_search_allowed'])
        self.assertTrue(d['common_joint24_storage_must_wait_for_finite_master'])
        self.assertIn('finite-state storage blocked',d['finite_storage_gate_error'])
        self.assertEqual(L.validate(d),[])
    def test_finite_guard_requires_every_branch_reference_and_zeroheel_scope(self):
        s=dict(map_representation='finite_physical_descriptor',finite_error_identity_for_every_event=True,
               physical_reference_forcing_retained=True,all_coefficient_product_graphs_retained=True,
               all_configured_branches_bound_to_finite_graph=True,zero_wind_heel_scope_enforced=True)
        P.assert_finite_storage_master(s)
        for key in tuple(s):
            q=dict(s);q.pop(key)
            with self.assertRaises(RuntimeError):P.assert_finite_storage_master(q)
    def test_storage_guard_names_every_required_closure(self):
        # Legacy guard remains available for old source-ledger callers, but it
        # is insufficient for finite storage and MUST NOT be used by Phase-1.
        s={k:True for k in ('same_history_complete_BRMM_word','physical_prediction_forcing_attached','physical_S_residual_attached','all_bias_families_attached','all_literal_branches_attached','H18_A21_edge_attached','zero_wind_heel_scope_enforced')}
        self.assertTrue(P.storage_search_allowed(s));P.assert_storage_search_allowed(s)
        with self.assertRaisesRegex(RuntimeError,'finite-state storage blocked'):
            P.assert_finite_storage_master(s)
    def test_windheel_scope_cannot_be_omitted(self):
        s={k:True for k in ('same_history_complete_BRMM_word','physical_prediction_forcing_attached','physical_S_residual_attached','all_bias_families_attached','all_literal_branches_attached','H18_A21_edge_attached')}
        self.assertFalse(P.storage_search_allowed(s))
        with self.assertRaisesRegex(RuntimeError,'zero_wind_heel_scope_enforced'):P.assert_storage_search_allowed(s)

if __name__=='__main__':unittest.main()
