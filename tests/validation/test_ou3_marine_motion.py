import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))

from tools.stability.ou3_theorem.marine_motion import (
    MarineContinuationCertificate,MarineLimits,MarineSample,audit_sampled_trace,
    constant_displacement_continuation_admitted,continuation_admitted,
    quiet_water_admitted,quiet_water_continuation_certificate,
)

LIMITS=MarineLimits(8.1,5.5,8.8,math.radians(35.0),1100.0)

class MarineMotionTests(unittest.TestCase):
    def test_permanent_nonzero_wave_displacement_is_rejected(self):
        self.assertFalse(constant_displacement_continuation_admitted((2.0,0.0,0.0)))
    def test_quiet_water_is_admitted(self):
        self.assertTrue(quiet_water_admitted())
    def test_same_history_and_bounded_potential_prefix(self):
        s=[
          MarineSample("h",0.0,(0,0,0),(1,0,0),(0,0,0),(0,0,0),(0,0,0)),
          MarineSample("h",1.0,(1,0,0),(1,0,0),(0,0,0),(0,0,0),(0.5,0,0)),
        ]
        r=audit_sampled_trace(s,LIMITS)
        self.assertTrue(r["finite_prefix_pass"],r["failures"])
        self.assertFalse(r["all_time_membership_certified_by_finite_trace"])
    def test_all_time_certificate_requires_attitude_rate_and_reference_acceleration(self):
        quiet=quiet_water_continuation_certificate("h")
        self.assertTrue(continuation_admitted(quiet,LIMITS))
        missing_attitude=MarineContinuationCertificate(
            "h",0.0,0.0,0.0,0.0,0.0,True,False,True,True,True
        )
        self.assertFalse(continuation_admitted(missing_attitude,LIMITS))
        hidden_reference_accel=MarineContinuationCertificate(
            "h",0.0,0.0,0.0,0.0,0.0,True,True,True,False,True
        )
        self.assertFalse(continuation_admitted(hidden_reference_accel,LIMITS))
        restarted_potential=MarineContinuationCertificate(
            "h",0.0,0.0,0.0,0.0,0.0,True,True,False,True,True
        )
        self.assertFalse(continuation_admitted(restarted_potential,LIMITS))

    def test_detached_coordinates_fail(self):
        s=[
          MarineSample("a",0.0,(0,0,0),(0,0,0),(0,0,0),(0,0,0),(0,0,0)),
          MarineSample("b",1.0,(0,0,0),(0,0,0),(0,0,0),(0,0,0),(0,0,0)),
        ]
        self.assertFalse(audit_sampled_trace(s,LIMITS)["finite_prefix_pass"])

if __name__=="__main__": unittest.main()
