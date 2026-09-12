"""Fresh H18 joint24 construction regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_fresh_joint24_entry as X
from tools.stability.ou3_alt_contraction import finite_startup_live_entry as LIVE
from tools.stability.ou3_alt_contraction import finite_startup_handoff_init as INIT
from tools.stability.ou3_alt_contraction import finite_startup_handoff_seed as SEED
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE


def active():
    sig=((F(1),0,0),(0,F(2),0),(0,0,F(3)))
    rs=((F(4),0,0),(0,F(5),0),(0,0,F(6)))
    return ACTIVE.ActiveParameters(F(11,10),sig,F(3,200),rs)


def entry():
    s=SEED.Result((1,0,0,0),SEED.TILT_SIGMA,SEED.YAW_SIGMA_GAUGED,False,True)
    # Distinct estimator coordinates make every truth-minus-estimate block visible.
    x=tuple([F(9),F(8),F(7)] + [F(i,10) for i in range(1,19)])
    P=[[F(1 if i==j else 0) for j in range(21)] for i in range(21)]
    h=INIT.initialize_from_gauged_seed_zero_heel(s,x,P,scope=SCOPE.certified_scope())
    return LIVE.enter_live(h,active(),scope=SCOPE.certified_scope())


def reference(time=0):
    return CORE.Reference(F(time),(1,0,0,0),
        (F(11,10),F(12,10),F(13,10)),
        (F(21,10),F(22,10),F(23,10)),
        (0,0,0),
        (F(31,10),F(32,10),F(33,10)),
        (F(41,100),F(42,100),F(43,100)),
        (F(51,100),F(52,100),F(53,100)),
        F(0),'fresh-history','fresh-bias','BIAS0')


class Tests(unittest.TestCase):
    def test_builds_exact_CORE_H_state_without_entry_box(self):
        e=entry(); r=reference(); s=X.build(e,r,scope=SCOPE.certified_scope())
        self.assertEqual(s.mode,'H'); self.assertEqual(s.reference,r); self.assertEqual(s.q_hat,e.q_hat)
        self.assertEqual(s.z[:3],(0,0,0))
        self.assertEqual(s.z[3:6],tuple(r.gyro_bias[i]-e.x[3+i] for i in range(3)))
        self.assertEqual(s.z[6:9],tuple(r.velocity[i]-e.x[6+i] for i in range(3)))
        self.assertEqual(s.z[9:12],tuple(r.position[i]-e.x[9+i] for i in range(3)))
        self.assertEqual(s.z[12:15],tuple(-e.x[12+i] for i in range(3)))
        self.assertEqual(s.z[15:18],tuple(r.acceleration[i]-e.x[15+i] for i in range(3)))
        self.assertEqual(s.z[18:21],tuple(r.beta[i]-e.x[18+i] for i in range(3)))
        self.assertEqual(s.z[21:24],r.beta)
        self.assertEqual(s.covariance,e.P)

    def test_internal_handoff_quaternion_is_world_to_body(self):
        e=entry(); self.assertEqual(e.q_hat,(1,0,0,0))
        s=X.build(e,reference(),scope=SCOPE.certified_scope())
        self.assertEqual(s.z[:3],(0,0,0))

    def test_reference_must_be_exact_live_origin(self):
        e=entry(); r=reference(F(1,200))
        with self.assertRaisesRegex(ValueError,'Live origin'):
            X.build(e,r,scope=SCOPE.certified_scope())

    def test_readiness_keeps_basin_and_binary32_open(self):
        r=X.readiness()
        self.assertTrue(r['fresh_entry_joint24_derived_from_truth_minus_estimator'])
        self.assertTrue(r['covariance_consistency_not_used_as_entry_assumption'])
        self.assertFalse(r['fresh_entry_inside_storage_basin_proved'])
        self.assertFalse(r['ALT_STARTUP_PASS']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
