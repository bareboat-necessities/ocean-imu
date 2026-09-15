from fractions import Fraction as F
import unittest
from tools.stability.ou3_alt_contraction import finite_admitted_machine_clock_qualified_interleaved_prefix as X
import test_finite_admitted_machine_joined_frontend_racc_interleaved_prefix as BASE
import test_finite_live_magnetic_word as MAG

class Tests(unittest.TestCase):
    def test_same_joined_event_gets_binary64_clock_qualification(self):
        base,kw,join,_,_,_=BASE.fixture(); s=X.begin(base)
        out=X.imu_step(s,separate_racc_accel_ldlt=MAG.REJECT,fma_racc_accel_ldlt=MAG.REJECT,**join,**kw)
        self.assertEqual(out.state.clock_steps,1)
        self.assertEqual(out.state.base.source_steps,base.source_steps+1)
        self.assertEqual(out.separate.due,out.fma.due)
        self.assertEqual(out.separate.deployed_cadence,X.CLOCK.ADAPT_FLOAT)

    def test_event_local_cadence_override_is_forbidden(self):
        base,kw,join,_,_,_=BASE.fixture(); s=X.begin(base)
        with self.assertRaisesRegex(TypeError,'carried by runtime'):
            X.imu_step(s,separate_racc_accel_ldlt=MAG.REJECT,fma_racc_accel_ldlt=MAG.REJECT,
                       aw_sync_adapt_every=F(1,5),**join,**kw)

    def test_runtime_and_deployment_cadence_are_same_compiled_value(self):
        base,_,_,_,_,_=BASE.fixture()
        exact,deployed=X._carried_cadence(base)
        self.assertEqual(exact,X.CLOCK.ADAPT_REAL)
        self.assertEqual(deployed,X.CLOCK.ADAPT_FLOAT)

    def test_complete_cannot_skip_600_clock_qualified_edges(self):
        base,_,_,_,_,_=BASE.fixture(); s=X.begin(base)
        with self.assertRaises((ValueError,TypeError)): X.complete(s)

    def test_readiness_does_not_unlock_storage(self):
        r=X.readiness()
        self.assertTrue(r['default_aw_sync_due_partition_binary64_closed'])
        self.assertTrue(r['complete_word_requires_clock_qualification_on_all_600_IMU_edges'])
        self.assertTrue(r['adapt_every_runtime_value_ancestry_closed_for_current_word'])
        self.assertTrue(r['canonical_5ms_source_dt_ancestry_closed_for_current_word'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])

if __name__=='__main__': unittest.main()
