import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))

from tools.stability.ou3_theorem.imu_bias import (
    BiasLimits,BiasSample,accel_prediction_error,audit_bias_trace,
    correction_error,estimator_prediction_factor,projected_error,projection_defect,
    successor_allowed,
)

LIMITS=BiasLimits(0.22516660498395405,0.001,0.02,1.0e-5)

class ImuBiasTests(unittest.TestCase):
    def test_successor_is_not_an_independent_box(self):
        self.assertFalse(successor_allowed((0,0,0),(0.1,0,0),0.005,LIMITS.B_a_mps2,LIMITS.D_a_mps3))
    def test_same_history_slow_drift_passes(self):
        s=[BiasSample("h",0.0,(0.01,0,0),(0.001,0,0)),
           BiasSample("h",1.0,(0.0105,0,0),(0.001005,0,0))]
        self.assertTrue(audit_bias_trace(s,LIMITS)["finite_prefix_pass"])
    def test_one_prediction_relation_covers_held_and_active(self):
        phi=math.exp(-0.1/5000.0)
        self.assertEqual(estimator_prediction_factor("H18",phi),1.0)
        self.assertEqual(estimator_prediction_factor("A21",phi),phi)
        held=accel_prediction_error("H18",phi,(0.02,0,0),(0.1,0,0),(0.001,0,0))
        active=accel_prediction_error("A21",phi,(0.02,0,0),(0.1,0,0),(0.001,0,0))
        self.assertAlmostEqual(held[0],0.021)
        self.assertAlmostEqual(active[0],phi*0.02+(1.0-phi)*0.1+0.001)
    def test_correction_and_projection_are_separate_operations(self):
        e_corr=correction_error((0.05,0,0),(0.01,0,0))
        self.assertEqual(e_corr,(0.04,0.0,0.0))
        defect=projection_defect((0.5,0,0),0.4)
        self.assertAlmostEqual(defect[0],0.1)
        e_plus=projected_error((0.2,0,0),(0.5,0,0),0.4)
        self.assertAlmostEqual(e_plus[0],-0.2)

if __name__=="__main__": unittest.main()
