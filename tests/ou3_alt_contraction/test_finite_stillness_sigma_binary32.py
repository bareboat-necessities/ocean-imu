"""Binary32 sigma-relevant StillnessAdapter projection regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as R
from tools.stability.ou3_alt_contraction import finite_stillness_sigma_binary32 as X
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as E


def expw(x):
    lo,hi=E.exp_minus_enclosure(x); return B.rn32((lo+hi)/2)


def energy_values(state,cfg,a):
    g=B.rn32(cfg.gravity); alpha=B.rn32(cfg.energy_alpha)
    an=B.div(a,g); inst=B.mul(an,an); decay=B.sub(B.rn32(1),alpha)
    return X._sum_products(decay,state.energy,alpha,inst)


class Tests(unittest.TestCase):
    def test_calm_sample_sets_actual_machine_still_predicate_and_time(self):
        cfg=R.Config(); s=X.State(); a=B.rn32(0); vals=energy_values(s,cfg,a); en=vals[0]
        h=B.rn32(F(1,200)); st=B.add(s.still_time,h); out=X.step(s,cfg,vertical_lp=a,dt=h,energy_successor=en,attenuation_exp=expw(st))
        self.assertTrue(out.state.is_still); self.assertEqual(out.state.still_time,st)
        self.assertEqual(out.attenuation,out.exp_result); self.assertIn(out.state.energy,out.energy_values)

    def test_moving_sample_resets_time_and_uses_literal_one_attenuation(self):
        cfg=R.Config(); s=X.State(B.rn32(0),B.rn32(1),True,7); a=B.rn32(20); vals=energy_values(s,cfg,a)
        # Every local outcome is far above the threshold for this input.
        en=max(vals); out=X.step(s,cfg,vertical_lp=a,dt=B.rn32(F(1,200)),energy_successor=en)
        self.assertFalse(out.state.is_still); self.assertEqual(out.state.still_time,0); self.assertEqual(out.attenuation,B.rn32(1)); self.assertIsNone(out.exp_result)

    def test_detached_energy_successor_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'outside local contraction set'):
            X.step(X.State(),R.Config(),vertical_lp=B.rn32(0),dt=B.rn32(F(1,200)),energy_successor=B.rn32(1),attenuation_exp=B.rn32(1))

    def test_full_sixty_second_attenuation_domain_has_rigorous_rne_witness(self):
        t=B.rn32(60); lo,hi=E.exp_minus_enclosure(t); e=B.rn32((lo+hi)/2)
        self.assertTrue(E._interval_hits_rne_cell(lo,hi,e))
        self.assertGreater(lo,0); self.assertLess(hi,F(1,10**20))

    def test_readiness_closes_projection_not_vertical_input_or_platform(self):
        r=X.readiness()
        for k in ('shipping_sigma_relevant_stillness_source_shape_matches','sigma_projection_independent_of_tracker_frequency_relaxation_state',
                  'energy_normalization_square_and_EMA_binary32_graph_materialized','actual_energy_successor_bound_to_all_local_contraction_choices',
                  'still_predicate_and_capped_time_derived_from_actual_machine_energy','attenuation_exp_bound_to_same_stored_machine_still_time_by_error_profile',
                  'moving_branch_attenuation_is_literal_one'):
            self.assertTrue(r[k])
        for k in ('target_exp_libm_correspondence_closed','upstream_vertical_LP_machine_production_closed',
                  'target_compiler_contraction_membership_closed','source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
