"""Same-event magnetic forcing regressions; not contraction/storage evidence."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_source_bound_mag_forcing as X
import test_finite_source_bound_live_word as BASE


class Tests(unittest.TestCase):
    def test_source_owned_mag_event_exposes_derived_effective_residual(self):
        s=BASE.root_state()
        out=X.mag_step(s,**BASE.mag_kwargs(s))
        event=out.word.event.event
        if event.qualification is None:
            self.assertIsNone(out.forcing)
            self.assertIsNone(event.effective_residual)
            return
        f=out.forcing
        self.assertEqual(f.raw_residual_body,event.qualification.sample.residual_body)
        self.assertEqual(f.effective_residual,event.effective_residual)
        expected=tuple(f.rotated_reference_mismatch[i]+f.hard_iron_body[i]
                       +f.applied_bias_correction[i]+f.raw_residual_body[i]
                       for i in range(3))
        self.assertEqual(f.effective_residual,expected)
        self.assertEqual(out.state.source.next_ordinal,1)

    def test_mag_forcing_record_rejects_detached_effective_residual(self):
        with self.assertRaisesRegex(ValueError,'detached'):
            X.MagForcing((1,0,0),(0,1,0),(0,0,1),(1,1,1),(2,2,2))

    def test_supply_vector_retains_components_instead_of_collapsing_box(self):
        f=X.MagForcing((1,0,0),(0,1,0),(0,0,-1),(2,3,4),(3,4,3))
        self.assertEqual(f.supply_vector,
                         (F(2),F(3),F(4),F(1),F(0),F(0),F(0),F(1),F(0),F(0),F(0),F(-1)))

    def test_readiness_does_not_promote_magnetic_smallness_or_storage(self):
        r=X.readiness()
        self.assertTrue(r['magnetic_effective_residual_derived_from_same_qualified_event'])
        self.assertTrue(r['raw_magnetic_residual_retained_separately'])
        self.assertTrue(r['active_reference_mismatch_retained_separately'])
        self.assertTrue(r['hard_iron_and_applied_correction_retained_separately'])
        self.assertTrue(r['independent_free_magnetic_innovation_forcing_forbidden_at_this_entry'])
        self.assertFalse(r['magnetic_forcing_smallness_claimed'])
        self.assertFalse(r['deployment_roundoff_supply_attached'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
