"""Magnetic theorem-source qualification regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_source_qualification as X
from tools.stability.ou3_alt_contraction import finite_mag_startup_source as S


def sample():
    p=S.PhysicalEndpoint(7,(1,0,0,0),'hist')
    m=S.Model((3,4,0),(F(3,5),F(4,5),0),'model')
    return S.make_sample(p,m,(0,F(3,4),F(4,3)),'m1')


def nw(v,n): return X.NormWitness(v,F(n))


class Tests(unittest.TestCase):
    def test_missing_theorem_envelope_fails_closed(self):
        s=sample()
        with self.assertRaisesRegex(ValueError,'undeclared'):
            X.qualify(s,None)
        with self.assertRaisesRegex(ValueError,'unqualified'):
            X.assert_source_qualified(None)

    def test_declared_envelope_qualifies_exact_same_source_sample(self):
        s=sample(); e=X.Envelope(5,1,F(25,12),'MAG-DET-1')
        q=X.qualify(s,e,
            field_norm=nw(s.model.world_field,5),
            hard_iron_norm=nw(s.model.hard_iron_body,1),
            residual_norm=nw(s.residual_body,F(25,12)))
        self.assertTrue(q.qualified); self.assertIs(q.sample,s); self.assertEqual(q.envelope.assumption_id,'MAG-DET-1')
        self.assertTrue(X.assert_source_qualified(q))

    def test_Rmag_or_detached_norm_cannot_substitute_for_source_bound(self):
        s=sample(); e=X.Envelope(5,1,F(25,12),'MAG-DET-1')
        with self.assertRaisesRegex(ValueError,'detached'):
            X.qualify(s,e,
                field_norm=nw((0,0,5),5),
                hard_iron_norm=nw(s.model.hard_iron_body,1),
                residual_norm=nw(s.residual_body,F(25,12)))

    def test_each_declared_bound_is_hard(self):
        s=sample()
        cases=(
          X.Envelope(F(49,10),1,F(25,12),'field'),
          X.Envelope(5,F(9,10),F(25,12),'hard'),
          X.Envelope(5,1,2,'noise'),
        )
        witnesses=dict(field_norm=nw(s.model.world_field,5),
                       hard_iron_norm=nw(s.model.hard_iron_body,1),
                       residual_norm=nw(s.residual_body,F(25,12)))
        for e in cases:
            with self.subTest(e=e.assumption_id), self.assertRaisesRegex(ValueError,'exceeds'):
                X.qualify(s,e,**witnesses)

    def test_readiness_blocks_complete_word_and_storage(self):
        r=X.readiness()
        self.assertTrue(r['missing_magnetic_envelope_fails_closed'])
        self.assertTrue(r['Rmag_is_not_used_as_deterministic_noise_bound'])
        self.assertFalse(r['canonical_theorem_currently_declares_magnetic_envelope'])
        self.assertFalse(r['source_uniform_magnetic_word_qualified'])
        self.assertFalse(r['complete_word_finite_identity'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_END_TO_END_PASS'])

if __name__=='__main__': unittest.main()
