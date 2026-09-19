import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))

from tools.stability.ou3_theorem.magnetic_service import (
    MagneticEvent,MagneticServiceContinuationCertificate,audit_window,
    continuation_admitted,transported_whitened_rows,
)

class MagneticServiceTests(unittest.TestCase):
    def event(self,t,rows,*,applied=True,gauged=True,finite=True,saturated=False):
        return MagneticEvent("h",t,applied,gauged,finite,saturated,tuple(rows))

    def test_source_uniform_service_requires_all_time_certificate(self):
        cert=MagneticServiceContinuationCertificate(
            "h",1.0,1.0,True,True,True,True)
        self.assertTrue(continuation_admitted(
            cert,required_window_s=1.0,required_information_floor=1.0))
        bad=MagneticServiceContinuationCertificate(
            "h",1.0,1.0,True,True,True,False)
        self.assertFalse(continuation_admitted(
            bad,required_window_s=1.0,required_information_floor=1.0))

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

    def test_rectangular_error_transport_keeps_original_information_root(self):
        # A synthetic dimension-changing map checks algebra only; it does not
        # claim to be a source-qualified shipping release certificate.
        transition = [[float(i == j) for j in range(18)] for i in range(21)]
        transition[18][2] = 0.25
        transition[18][5] = 0.5
        injection = [[0.0, 0.0] for _ in range(18)]
        injection[2][0] = 1.0
        injection[5][1] = 0.02
        H = [[0.0]*21 for _ in range(3)]
        H[0][2], H[1][5], H[2][18] = 1.0, 1.0, 1.0
        S = ((1,0,0),(0,1,0),(0,0,1))
        rows = transported_whitened_rows(H, S, transition, injection)
        self.assertEqual(rows, ((1.0,0.0),(0.0,0.02),(0.25,0.01)))
        with self.assertRaises(ValueError):
            transported_whitened_rows(H, S, transition, injection + [[0,0]])

if __name__=="__main__": unittest.main()
