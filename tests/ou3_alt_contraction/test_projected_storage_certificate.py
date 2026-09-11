import unittest
import numpy as np

from ou3_interval import Interval, matrix_point
from tools.stability.ou3_alt_contraction import projected_storage_certificate as C


class ProjectedStorageCertificateTest(unittest.TestCase):
    def test_point_interval_family_certifies_neutral_supply_example(self):
        A=matrix_point([[0.5,0.7],[0.0,1.0]])
        M=np.array([[2.0,0.3],[0.3,1.0]])
        d=C.certify_projected_interval_family(A,M,0.95,(1,))
        self.assertTrue(d['metric_spd'])
        self.assertTrue(d['projected_strict'])
        self.assertGreater(d['worst_projected_ldlt_pivot_lower'],0.0)
        self.assertFalse(d['ordinary_eigenvalue_used_for_certificate'])

    def test_interval_family_is_certified_as_family_not_midpoint(self):
        # A11 in [0.45,0.55] is still uniformly contracting on ker(C).
        A=[[Interval.outward_bounds(0.45,0.55),Interval.outward_bounds(0.6,0.8)],
           [Interval.point(0.0),Interval.point(1.0)]]
        d=C.certify_projected_interval_family(A,np.eye(2),0.9,(1,))
        self.assertTrue(d['projected_strict'])
        self.assertTrue(d['interval_family_consumed'])

    def test_unsupplied_neutral_family_fails_outward_ldlt(self):
        A=matrix_point([[1.0,0.0],[0.0,1.0]])
        d=C.certify_projected_interval_family(A,np.eye(2),0.9,(1,))
        self.assertTrue(d['metric_spd'])
        self.assertFalse(d['projected_strict'])

    def test_non_spd_metric_fails_before_storage_claim(self):
        A=matrix_point([[0.5,0.0],[0.0,1.0]])
        d=C.certify_projected_interval_family(A,[[1.0,0.0],[0.0,-1.0]],0.9,(1,))
        self.assertFalse(d['metric_spd'])
        self.assertFalse(d['projected_strict'])

    def test_theorem_backend_remains_fail_closed_without_universal_family(self):
        d=C.build(); f=C.validate(d)
        self.assertEqual(f,[])
        self.assertTrue(d['outward_full_matrix_LDLT_terminal_gate'])
        self.assertFalse(d['source_uniform_endpoint_family_enclosure_closed'])
        self.assertFalse(d['common_M_source_uniform_projected_LDLT_closed'])
        self.assertFalse(d['ALT_LIVE_PASS'])


if __name__=='__main__':
    unittest.main()
