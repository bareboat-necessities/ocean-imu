#!/usr/bin/env python3
from __future__ import annotations
import copy,sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_brmm_same_Dtheta_reset_binding as B
import ou3_p4_same_cell_correction_domain as C

class SameDthetaResetFailClosedTest(unittest.TestCase):
    def test_none_legacy_delta_remains_open_and_production_route_stays_structural(self):
        corr=C.build();self.assertEqual(C.validate(corr),[])
        sz,_=B._smokes();mode='H18';kind=sz.kind
        mutated=copy.deepcopy(corr)
        mutated['modes'][mode]['events'][kind]['same_cell_attitude_correction_norm_upper']=None
        legacy=B.bind_event(sz,mutated)
        self.assertIsNone(legacy['delta'])
        self.assertFalse(legacy['closed']);self.assertFalse(legacy['chart_safe'])
        self.assertFalse(legacy['scalar_delta_certified'])
        prod=B.bind_event_first_exit(sz)
        self.assertTrue(prod['closed'])
        self.assertFalse(prod['correction_domain_target_proved_here'])
        self.assertTrue(prod['correction_domain_target_must_be_proved_by_same_augmented_first_exit_master'])
        self.assertFalse(prod['rowwise_K_bound_used'])
        self.assertFalse(prod['scalar_covariance_residual_ceiling_used'])

    def test_module_remains_fail_closed_for_p4(self):
        d=B.build();self.assertEqual(B.validate(d),[])
        self.assertTrue(d['legacy_uncertified_scalar_delta_fails_closed_without_exception'])
        self.assertTrue(d['first_exit_same_Dtheta_reset_structure_closed'])
        self.assertFalse(d['production_correction_domain_target_closed_here'])
        self.assertFalse(d['P4_PASS']);self.assertFalse(d['P5_MAY_START'])

if __name__=='__main__':unittest.main()
