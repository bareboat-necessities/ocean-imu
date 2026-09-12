"""Finite magnetic source / Rmag / measurement regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_runtime as X
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MR
import test_finite_core as FC


class Tests(unittest.TestCase):
    def sample(self,state=None):
        s=state or FC.root('H')
        model=X.Model((22,0,43),(F(1,5),F(1,4),F(1,3)))
        n=(F(1,100),F(-1,200),F(1,300))
        phys=X.SENSOR.q_rotate(s.reference.q_world_to_body,model.world_reference)
        raw=tuple(phys[i]+n[i] for i in range(3))
        return s,X.Sample(s.reference,raw,n,model)

    def test_packet_identity_and_Rmag_are_not_free(self):
        s,m=self.sample()
        self.assertEqual(m.observed_internal,tuple(X.SENSOR.q_rotate(s.reference.q_world_to_body,m.model.world_reference)[i]+m.residual_internal[i] for i in range(3)))
        self.assertEqual(m.model.covariance,((F(1,25),0,0),(0,F(1,16),0),(0,0,F(1,9))))
        with self.assertRaisesRegex(ValueError,'detached'):
            X.Sample(s.reference,(0,0,0),m.residual_internal,m.model)

    def test_sanity_rejection_precedes_LDLT_and_is_identity(self):
        s=FC.root('H'); model=X.Model((0,0,0),(1,1,1)); sample=X.Sample(s.reference,(0,0,0),(0,0,0),model)
        out=X.update(s,sample)
        self.assertFalse(out.attempted_measurement); self.assertEqual(out.state,s)
        with self.assertRaisesRegex(ValueError,'consumes no LDLT'):
            X.update(s,sample,ldlt=MR.SafeLDLT(True,None,1,F(1,10**7)))

    def test_safe_ldlt_rejection_still_represents_attempted_measurement(self):
        s,m=self.sample(); ldlt=MR.SafeLDLT(False,False,1,F(1,10**7))
        out=X.update(s,m,ldlt=ldlt)
        self.assertTrue(out.attempted_measurement); self.assertFalse(out.accepted)
        self.assertEqual(out.state,s)

    def test_accepted_mag_uses_same_packet_world_reference_and_covariance(self):
        s,m=self.sample(); ldlt=MR.SafeLDLT(True,None,1,F(1,10**7))
        out=X.update(s,m,ldlt=ldlt)
        self.assertTrue(out.accepted); self.assertIsNotNone(out.measurement.accepted_graph)
        self.assertEqual(out.measurement.accepted_graph.state,out.state)
        self.assertNotEqual(out.state,s)

    def test_detached_endpoint_and_free_ldlt_are_rejected(self):
        s,m=self.sample()
        # Use another internally valid finite state rooted at a different
        # physical history.  Do not manufacture an inconsistent State merely
        # to reach the endpoint guard: CORE.State itself must remain fail-closed.
        other=FC.root('H','BIAS2')
        with self.assertRaisesRegex(ValueError,'current physical endpoint'):
            X.update(other,m,ldlt=MR.SafeLDLT(True,None,1,F(1,10**7)))
        with self.assertRaisesRegex(TypeError,'requires safe-LDLT'):
            X.update(s,m)

    def test_readiness_is_fail_closed(self):
        r=X.readiness()
        self.assertTrue(r['same_mag_packet_supplies_observation_and_residual'])
        self.assertTrue(r['Rmag_diagonal_from_configured_sigma_only'])
        self.assertFalse(r['mag_world_reference_startup_ancestry_attached'])
        self.assertFalse(r['mag_noise_source_bound_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
