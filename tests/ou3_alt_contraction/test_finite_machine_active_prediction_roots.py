"""Machine-applied prediction-root provenance regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_machine_active_prediction_roots as X
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
import test_finite_source_bound_live_word as BASE
import test_finite_source_bound_prediction_word as PBASE


class Tests(unittest.TestCase):
    def test_same_active_parameters_reproduce_exact_source_bound_roots(self):
        s=BASE.root_state(); witness,segment,raw,_=PBASE.operands(s)
        physical=SOURCE.QualifiedPhysicalSegment(s.source.root,witness,segment)
        active=s.live.live.live.active; args=PBASE.root_args(); args.pop('temperature_c')
        m=X.build(s,physical,raw,active,mode='separate',**args)
        e=PBASE.X.build(s,physical,raw,**args)
        self.assertEqual(m.angular,e.angular); self.assertEqual(m.ou,e.ou)
        self.assertEqual(m.qaxis,e.qaxis); self.assertEqual(m.bias,e.bias)
        self.assertEqual(m.active_join.supply.tau,0)
        self.assertEqual(m.active_join.supply.Sigma_aw,((0,0,0),(0,0,0),(0,0,0)))

    def test_nonprediction_active_fields_remain_visible_without_changing_roots(self):
        s=BASE.root_state(); witness,segment,raw,_=PBASE.operands(s)
        physical=SOURCE.QualifiedPhysicalSegment(s.source.root,witness,segment)
        exact=s.live.live.live.active
        changed=ACTIVE.ActiveParameters(exact.tau,exact.Sigma_aw,exact.pseudo_period+F(1,100),
                                        tuple(tuple(x+F(1,100) if i==j else x for j,x in enumerate(row))
                                              for i,row in enumerate(exact.R_S)))
        args=PBASE.root_args(); args.pop('temperature_c')
        out=X.build(s,physical,raw,changed,mode='fma',**args)
        self.assertEqual(out.ou.tau,exact.tau); self.assertEqual(out.qaxis.sigma_aw,exact.Sigma_aw)
        self.assertEqual(out.active_join.supply.pseudo_period,F(1,100))
        self.assertNotEqual(out.active_join.supply.R_S,((0,0,0),(0,0,0),(0,0,0)))

    def test_machine_active_cannot_detach_from_carried_source_root(self):
        s=BASE.root_state(); witness,segment,raw,_=PBASE.operands(s)
        other=replace(s.source.root,history_id='different-machine-root')
        physical=SOURCE.QualifiedPhysicalSegment(other,witness,segment)
        args=PBASE.root_args(); args.pop('temperature_c')
        with self.assertRaisesRegex(ValueError,'next carried source transition'):
            X.build(s,physical,raw,s.live.live.live.active,mode='separate',**args)

    def test_readiness_keeps_event_injection_and_storage_open(self):
        r=X.readiness()
        self.assertTrue(r['machine_prediction_root_relation_attached'])
        self.assertTrue(r['small_general_Qaxis_branch_is_recomputed_from_machine_tau_not_shadow_tau'])
        for k in ('machine_root_effect_injected_into_joint24_event_relation','machine_pseudo_period_scheduler_effect_attached',
                  'machine_RS_measurement_effect_attached','OU_Qaxis_target_libm_and_Eigen_correspondence_closed',
                  'source_uniform_machine_coefficient_supply_bound_closed','source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
