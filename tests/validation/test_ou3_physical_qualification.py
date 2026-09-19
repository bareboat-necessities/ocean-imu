from __future__ import annotations
import sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))

from tools.stability.ou3_theorem.marine_motion import (
    MarineContinuationCertificate,MarineLimits,
)
from tools.stability.ou3_theorem.imu_bias import BiasContinuationCertificate,BiasLimits
from tools.stability.ou3_theorem.magnetic_service import MagneticServiceContinuationCertificate
from tools.stability.ou3_theorem.physical_qualification import (
    PhysicalQualification,qualification_admitted,
)

class PhysicalQualificationTests(unittest.TestCase):
    def test_same_history_all_time_composition(self):
        m=MarineContinuationCertificate("h",1,1,1,.1,1,True,True,True,True,True)
        b=BiasContinuationCertificate("h",.1,.01,.01,.001,True,True,True)
        g=MagneticServiceContinuationCertificate("h",1,1,True,True,True,True)
        q=PhysicalQualification(m,b,g,True,True)
        r=qualification_admitted(q,marine_limits=MarineLimits(2,2,2,.2,2),
            bias_limits=BiasLimits(.2,.02,.02,.002),
            magnetic_window_s=1,magnetic_information_floor=1)
        self.assertTrue(r["qualified"])

    def test_detached_history_fails(self):
        m=MarineContinuationCertificate("m",1,1,1,.1,1,True,True,True,True,True)
        b=BiasContinuationCertificate("b",.1,.01,.01,.001,True,True,True)
        g=MagneticServiceContinuationCertificate("g",1,1,True,True,True,True)
        q=PhysicalQualification(m,b,g,True,True)
        r=qualification_admitted(q,marine_limits=MarineLimits(2,2,2,.2,2),
            bias_limits=BiasLimits(.2,.02,.02,.002),
            magnetic_window_s=1,magnetic_information_floor=1)
        self.assertFalse(r["qualified"]);self.assertFalse(r["same_history"])

if __name__=="__main__":unittest.main()
