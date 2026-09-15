import unittest
from fractions import Fraction as F
from tools.stability.ou3_alt_contraction import finite_startup_handoff_seed as SEED
from tools.stability.ou3_alt_contraction import finite_startup_handoff_init as INIT
from tools.stability.ou3_alt_contraction import finite_startup_live_entry as LIVE
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as TILT
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE
class Tests(unittest.TestCase):
    def test_ungauged_proxy_normalizes_and_enters_h18_with_free_yaw_covariance(self):
        seed=SEED.seed(V.State(q=(F(2),0,0,0),initialized=True),None);self.assertFalse(seed.gauged)
        P=tuple(tuple(F(1) if i==j else F(0) for j in range(21)) for i in range(21));x=(F(0),)*21
        h=INIT.initialize_from_seed_zero_heel(seed,x,P,scope=SCOPE.Scope(),q_norm=TILT.SqrtWitness(F(4),F(2)))
        self.assertFalse(h.gauged);self.assertEqual(h.q_seed,(F(1),0,0,0));self.assertEqual(h.yaw_sigma,SEED.YAW_SIGMA_FREE)
        active=ACTIVE.ActiveParameters(F(1),((F(1),0,0),(0,F(1),0),(0,0,F(1))),F(1),((F(1),0,0),(0,F(1),0),(0,0,F(1))))
        live=LIVE.enter_live(h,active,scope=SCOPE.Scope());self.assertFalse(live.gauged);self.assertEqual(live.startup_stage,'Live');self.assertFalse(live.acc_bias_updates_enabled)
    def test_ungauged_requires_same_quaternion_norm_witness(self):
        seed=SEED.seed(V.State(q=(F(2),0,0,0),initialized=True),None);P=tuple(tuple(F(1) if i==j else F(0) for j in range(21)) for i in range(21))
        with self.assertRaises((TypeError,ValueError)): INIT.initialize_from_seed_zero_heel(seed,(F(0),)*21,P,scope=SCOPE.Scope(),q_norm=TILT.SqrtWitness(F(1),F(1)))
if __name__=='__main__': unittest.main()
