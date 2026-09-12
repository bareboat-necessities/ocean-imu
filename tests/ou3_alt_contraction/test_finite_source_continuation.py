from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_source_continuation as X
from tools.stability.ou3_alt_contraction import finite_physical_prediction as P


def kin(t, *, live=F(0), beta=(0,0,0)):
    return P.PhysicalKinematics(F(t),(1,0,0,0),(0,0,0),(0,0,0),(0,0,0),
                                (0,0,0),(0,0,0),beta,F(live))


def seg(k, *, phi=F(1), driver=(0,0,0), beta0=(0,0,0), beta1=None):
    t0=F(k-1,200); t1=F(k,200)
    if beta1 is None:
        beta1=tuple(F(phi)*F(beta0[i])+F(driver[i]) for i in range(3))
    return P.PhysicalSegment(kin(t0,beta=beta0),kin(t1,beta=beta1),
                             (0,0,0),(0,0,0),(0,0,0),F(phi),driver)


def root():
    return X.SourceRoot('hist','generator',0,'BIAS0',
                        'BIAS0:one-physical-history:phi-root-and-driver')


def wit(k,parent,child,pin,pout):
    return X.StepWitness(k,parent,child,pin,pout)


class Tests(unittest.TestCase):
    def test_consecutive_segments_preserve_one_source_and_bias_root(self):
        r=root(); c=X.begin(r)
        c=X.append(c,witness=wit(1,'root','c1','p0','p1'),segment=seg(1))
        c=X.append(c,witness=wit(2,'c1','c2','p1','p2'),segment=seg(2))
        self.assertEqual(c.next_ordinal,3); self.assertFalse(c.complete)
        self.assertEqual(c.steps[0].root,c.steps[1].root)
        self.assertEqual(c.steps[0].segment.after,c.steps[1].segment.before)

    def test_bias_family_phi_is_not_free_per_segment(self):
        r=root()
        with self.assertRaises(ValueError):
            X.QualifiedPhysicalSegment(r,wit(1,'root','c1','p0','p1'),seg(1,phi=F(1,2)))

    def test_source_cell_and_primitive_chains_are_hard(self):
        r=root(); c=X.begin(r)
        c=X.append(c,witness=wit(1,'root','c1','p0','p1'),segment=seg(1))
        with self.assertRaises(ValueError):
            X.append(c,witness=wit(2,'wrong','c2','p1','p2'),segment=seg(2))
        with self.assertRaises(ValueError):
            X.append(c,witness=wit(2,'c1','c2','wrong','p2'),segment=seg(2))

    def test_physical_endpoint_chain_cannot_jump(self):
        r=root(); q1=X.QualifiedPhysicalSegment(r,wit(1,'root','c1','p0','p1'),seg(1))
        # Same canonical clock, but change the next predecessor quaternion so it
        # is not the prior exact physical endpoint. Segment algebra remains valid.
        b=P.PhysicalKinematics(F(1,200),(0,1,0,0),(0,0,0),(0,0,0),(0,0,0),
                               (0,0,0),(0,0,0),(0,0,0),0)
        a=P.PhysicalKinematics(F(2,200),(0,1,0,0),(0,0,0),(0,0,0),(0,0,0),
                               (0,0,0),(0,0,0),(0,0,0),0)
        s2=P.PhysicalSegment(b,a,(0,0,0),(0,0,0),(0,0,0),1,(0,0,0))
        q2=X.QualifiedPhysicalSegment(r,wit(2,'c1','c2','p1','p2'),s2)
        with self.assertRaises(ValueError): X.Continuation(r,(q1,q2))

    def test_readiness_closes_ancestry_but_not_finite_master(self):
        d=X.readiness()
        self.assertTrue(d['correlated_COMPLETE_BRMM_left_inclusion_consumed'])
        self.assertTrue(d['bias_phi_driver_and_true_beta_hard_contracts_checked_per_segment'])
        self.assertTrue(d['source_cell_parent_child_and_primitive_continuity_checked'])
        self.assertFalse(d['finite_estimator_coefficients_bound_to_same_source_continuation'])
        self.assertFalse(d['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(d['storage_search_allowed']); self.assertFalse(d['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
