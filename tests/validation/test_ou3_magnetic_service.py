import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))

from tools.stability.ou3_theorem.magnetic_service import (
    MagneticEvent,audit_window,transported_whitened_heading_bias_rows,
)

class MagneticServiceTests(unittest.TestCase):
    def test_rows_are_built_from_actual_H_S_and_preceding_transition(self):
        H=((1.0,0.0),(0.0,2.0),(0.0,0.0))
        S=((4.0,0.0,0.0),(0.0,9.0,0.0),(0.0,0.0,1.0))
        PhiE=((1.0,0.0),(0.0,1.0))
        G=transported_whitened_heading_bias_rows(H,S,PhiE)
        self.assertAlmostEqual(G[0][0],0.5)
        self.assertAlmostEqual(G[1][1],2.0/3.0)
        self.assertEqual(G[2],(0.0,0.0))

    def event(self,t,rows,*,applied=True,gauged=True,finite=True,saturated=False):
        return MagneticEvent("h",t,applied,gauged,finite,saturated,tuple(rows))
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
