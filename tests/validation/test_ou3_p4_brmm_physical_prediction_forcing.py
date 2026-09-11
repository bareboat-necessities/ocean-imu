#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
from ou3_interval import Interval
import ou3_p4_brmm_physical_prediction_forcing as F

class PhysicalPredictionForcingTest(unittest.TestCase):
    def test_exact_subtraction_identity_and_fail_closed_status(self):
        d=F.build();self.assertEqual(F.validate(d),[])
        self.assertTrue(d['physical_source_forcing_matrix_materialized'])
        self.assertTrue(d['physical_moment_forcing_and_latent_acceleration_mismatch_combined_once'])
        self.assertTrue(d['point_identity_check']['passed'])
        self.assertFalse(d['source_uniform_quadratic_enclosure_of_full_15D_source_closed_here'])
        self.assertFalse(d['P4_PASS'])
    def test_zero_physical_OU_defect_when_witness_matches_shipping_OU(self):
        I=F.I;h=I(.005);tau=I(2.0)
        axis=F.SHIPPING.translation_axis_transition(tau,h);a0=(I(.3),I(-.2),I(.1));alpha=axis[3][3]
        a1=tuple(alpha*x for x in a0)
        J0=tuple(axis[0][3]*x for x in a0);J1=tuple(axis[1][3]*x for x in a0);J2=tuple(axis[2][3]*x for x in a0)
        q=F.witness(a0=a0,a1=a1,J0=J0,J1=J1,J2=J2);d=F.forcing(q,tau,h)
        for v in (*d.dv,*d.dp,*d.dS,*d.daw):self.assertLessEqual(v.abs_upper(),1e-14)
    def test_source_matrix_has_required_signs(self):
        I=F.I;M=F.source_matrix(I(2.0),I(.005));axis=F.SHIPPING.translation_axis_transition(I(2.0),I(.005))
        self.assertLess(M[0][0].hi,0.0);self.assertEqual(M[0][6].lo,1.0)
        self.assertLess(M[9][0].hi,0.0);self.assertEqual(M[9][3].lo,1.0)
        self.assertLessEqual(abs(M[0][0].lo+axis[0][3].hi),1e-12)
    def test_H_A_injections_do_not_touch_attitude_or_bias_rows(self):
        M=F.source_matrix(F.I(2.0),F.I(.005))
        for mode,n in [('H',18),('A',21)]:
            G=F.inject_error_state(mode,M);self.assertEqual((len(G),len(G[0])),(n,15))
            self.assertTrue(all(x.lo==0 and x.hi==0 for row in G[:6] for x in row))
            if mode=='A':self.assertTrue(all(x.lo==0 and x.hi==0 for row in G[18:21] for x in row))

if __name__=='__main__':unittest.main()
