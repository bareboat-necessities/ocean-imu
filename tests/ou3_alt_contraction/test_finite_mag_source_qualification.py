"""Magnetic theorem-source qualification regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_source_qualification as X
from tools.stability.ou3_alt_contraction import finite_mag_startup_source as S


def sample(field=(36,27,0), hard=(6,8,0), residual=(0,0,3)):
    p=S.PhysicalEndpoint(7,(1,0,0,0),'hist')
    m=S.Model(field,hard,'model')
    return S.make_sample(p,m,residual,'m1')


def nw(v,n): return X.NormWitness(v,F(n))
def hw(v,n): return X.HorizontalNormWitness(tuple(v[:2]),F(n))


def qualify(s,e=None,field_n=45,horiz_n=45,hard_n=10,resid_n=3):
    return X.qualify(s,e,
        field_norm=nw(s.model.world_field,field_n),
        horizontal_norm=hw(s.model.world_field,horiz_n),
        hard_iron_norm=nw(s.model.hard_iron_body,hard_n),
        residual_norm=nw(s.residual_body,resid_n))


class Tests(unittest.TestCase):
    def test_default_named_BMM150_envelope_qualifies_boundary_sample(self):
        s=sample(); q=qualify(s)
        self.assertTrue(q.qualified); self.assertIs(q.sample,s)
        self.assertEqual(q.envelope.assumption_id,'MAG-BMM150-DET-v1')
        self.assertEqual(q.envelope.world_field_norm_min,20)
        self.assertEqual(q.envelope.world_field_norm_max,75)
        self.assertEqual(q.envelope.world_field_horizontal_min,15)
        self.assertEqual(q.envelope.hard_iron_norm_max,10)
        self.assertEqual(q.envelope.residual_norm_max,3)
        self.assertTrue(X.assert_source_qualified(q))

    def test_wrong_named_assumption_fails_closed(self):
        s=sample(); e=X.Envelope(20,75,15,10,3,'OTHER')
        with self.assertRaisesRegex(ValueError,'MAG-BMM150'):
            qualify(s,e)

    def test_Rmag_or_detached_norm_cannot_substitute_for_source_bound(self):
        s=sample()
        with self.assertRaisesRegex(ValueError,'detached'):
            X.qualify(s,
                field_norm=nw((0,0,45),45), horizontal_norm=hw(s.model.world_field,45),
                hard_iron_norm=nw(s.model.hard_iron_body,10),
                residual_norm=nw(s.residual_body,3))
        with self.assertRaisesRegex(ValueError,'unqualified'):
            X.assert_source_qualified(None)

    def test_total_field_lower_and_upper_bounds_are_hard(self):
        low=sample(field=(12,16,0))  # norm 20, accepted boundary
        self.assertTrue(X.qualify(low,field_norm=nw(low.model.world_field,20),
                                  horizontal_norm=hw(low.model.world_field,20),
                                  hard_iron_norm=nw(low.model.hard_iron_body,10),
                                  residual_norm=nw(low.residual_body,3)).qualified)
        below=sample(field=(9,12,0)) # norm/horizontal 15, below total minimum
        with self.assertRaisesRegex(ValueError,'world field violates'):
            X.qualify(below,field_norm=nw(below.model.world_field,15),
                      horizontal_norm=hw(below.model.world_field,15),
                      hard_iron_norm=nw(below.model.hard_iron_body,10),
                      residual_norm=nw(below.residual_body,3))
        high=sample(field=(60,45,0)) # norm 75 accepted boundary
        self.assertTrue(X.qualify(high,field_norm=nw(high.model.world_field,75),
                                  horizontal_norm=hw(high.model.world_field,75),
                                  hard_iron_norm=nw(high.model.hard_iron_body,10),
                                  residual_norm=nw(high.residual_body,3)).qualified)
        over=sample(field=(64,48,0)) # norm 80
        with self.assertRaisesRegex(ValueError,'world field violates'):
            X.qualify(over,field_norm=nw(over.model.world_field,80),
                      horizontal_norm=hw(over.model.world_field,80),
                      hard_iron_norm=nw(over.model.hard_iron_body,10),
                      residual_norm=nw(over.residual_body,3))

    def test_horizontal_yaw_observability_bound_is_hard(self):
        s=sample(field=(9,12,16)) # total 21.93..., horizontal exactly 15; use 3-4-? exact total is sqrt481 not rational
        # Use a rational-Pythagorean total instead: (9,12,20) -> total 25, horizontal 15.
        s=sample(field=(9,12,20))
        self.assertTrue(X.qualify(s,field_norm=nw(s.model.world_field,25),
                                  horizontal_norm=hw(s.model.world_field,15),
                                  hard_iron_norm=nw(s.model.hard_iron_body,10),
                                  residual_norm=nw(s.residual_body,3)).qualified)
        bad=sample(field=(F(42,5),F(56,5),24)) # horizontal 14, total sqrt(772) not exact rational
        # Exact rational construction: (0,14,48) has total 50.
        bad=sample(field=(0,14,48))
        with self.assertRaisesRegex(ValueError,'yaw-observability'):
            X.qualify(bad,field_norm=nw(bad.model.world_field,50),
                      horizontal_norm=hw(bad.model.world_field,14),
                      hard_iron_norm=nw(bad.model.hard_iron_body,10),
                      residual_norm=nw(bad.residual_body,3))

    def test_hard_iron_and_residual_bounds_are_hard(self):
        hard=sample(hard=(0,0,11))
        with self.assertRaisesRegex(ValueError,'hard iron exceeds'):
            X.qualify(hard,field_norm=nw(hard.model.world_field,45),
                      horizontal_norm=hw(hard.model.world_field,45),
                      hard_iron_norm=nw(hard.model.hard_iron_body,11),
                      residual_norm=nw(hard.residual_body,3))
        noise=sample(residual=(0,0,4))
        with self.assertRaisesRegex(ValueError,'mag residual exceeds'):
            X.qualify(noise,field_norm=nw(noise.model.world_field,45),
                      horizontal_norm=hw(noise.model.world_field,45),
                      hard_iron_norm=nw(noise.model.hard_iron_body,10),
                      residual_norm=nw(noise.residual_body,4))

    def test_readiness_still_blocks_word_storage_but_not_source_spec(self):
        r=X.readiness()
        self.assertTrue(r['named_BMM150_deterministic_envelope_declared'])
        self.assertTrue(r['horizontal_yaw_observability_lower_bound_declared'])
        self.assertTrue(r['canonical_theorem_currently_declares_magnetic_envelope'])
        self.assertTrue(r['individual_magnetic_sample_can_be_source_qualified'])
        self.assertTrue(r['Rmag_is_not_used_as_deterministic_noise_bound'])
        self.assertFalse(r['source_uniform_magnetic_word_qualified'])
        self.assertFalse(r['complete_word_finite_identity'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_END_TO_END_PASS'])

if __name__=='__main__': unittest.main()
