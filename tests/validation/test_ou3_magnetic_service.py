import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))

from tools.stability.ou3_theorem.magnetic_service import (
    MagneticEvent,audit_window,transported_whitened_rows,
)

class MagneticServiceTests(unittest.TestCase):
    def event(self,t,rows,*,applied=True,gauged=True,finite=True,saturated=False):
        return MagneticEvent("h",t,applied,gauged,finite,saturated,tuple(rows))

    def test_whitening_uses_actual_innovation_covariance(self):
        # H=diag(1,2), S=diag(4,9), Phi=E=I gives
        # G=diag(1/2,2/3). Replacing S by measurement R would change this.
        rows=transported_whitened_rows(
            ((1.0,0.0),(0.0,2.0)),
            ((4.0,0.0),(0.0,9.0)),
            ((1.0,0.0),(0.0,1.0)),
            ((1.0,0.0),(0.0,1.0)),
        )
        self.assertAlmostEqual(rows[0][0],0.5)
        self.assertAlmostEqual(rows[0][1],0.0)
        self.assertAlmostEqual(rows[1][0],0.0)
        self.assertAlmostEqual(rows[1][1],2.0/3.0)

    def test_complete_preceding_transition_is_consumed(self):
        rows=transported_whitened_rows(
            ((1.0,0.0),(0.0,1.0)),
            ((1.0,0.0),(0.0,1.0)),
            ((1.0,3.0),(0.0,1.0)),
            ((1.0,0.0),(0.0,1.0)),
        )
        self.assertEqual(rows,((1.0,3.0),(0.0,1.0)))

    def test_only_actually_applied_measurements_supply_information(self):
        events=[self.event(0.2,((4.0,0.0),),applied=False),
                self.event(0.4,((0.0,4.0),),finite=False),
                self.event(0.6,((4.0,0.0),)),
                self.event(0.8,((0.0,4.0),))]
        r=audit_window(events,0.0,1.0,1.0)
        self.assertEqual(r["actually_applied_informative_events"],2)
        self.assertTrue(r["information_pass"])

    def test_small_event_gap_is_not_enough_without_rank(self):
        r=audit_window([self.event(0.1*k,((2.0,0.0),)) for k in range(1,10)],0.0,1.0,1.0)
        self.assertLess(r["max_usable_event_gap_s_diagnostic_only"],0.2)
        self.assertEqual(r["information_min_eigenvalue"],0.0)
        self.assertFalse(r["information_pass"])

    def test_saturated_measurement_does_not_count(self):
        r=audit_window([self.event(0.5,((3.0,0.0),(0.0,3.0)),saturated=True)],0.0,1.0,1.0)
        self.assertEqual(r["actually_applied_informative_events"],0)
        self.assertFalse(r["information_pass"])

if __name__=="__main__": unittest.main()
