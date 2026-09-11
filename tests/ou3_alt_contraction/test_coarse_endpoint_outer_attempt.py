import unittest
from ou3_interval import matrix_point,matrix_mul
from tools.stability.ou3_alt_contraction import coarse_endpoint_outer_attempt as E


class CoarseEndpointOuterAttemptTest(unittest.TestCase):
    def test_interval_power_contains_literal_two_branch_products(self):
        # Small exact regression: hull({A,B})^3 contains every explicit product.
        A=matrix_point([[0.5,0.0],[0.0,1.0]]);B=matrix_point([[0.8,0.0],[0.0,1.0]])
        from tools.stability.ou3_alt_contraction.endpoint_family_induction import hull_matrices,matrix_encloses
        H=hull_matrices([A,B]);P=E.interval_matrix_power(H,3)
        for X in (A,B):
            for Y in (A,B):
                for Z in (A,B):
                    self.assertTrue(matrix_encloses(P,matrix_mul(Z,matrix_mul(Y,X))))

    def test_first_600_step_universal_outer_attempt(self):
        d=E.build();self.assertEqual(E.validate(d),[])
        print('ALT_COARSE_600_ENDPOINT',E.summary(d))
        self.assertTrue(d['binary_interval_power_encloses_all_length_600_IMU_core_products'])
        self.assertTrue(d['async_magnetometer_excluded_and_star_obligation_retained'])
        # Finiteness is diagnostic of representation quality, never theorem truth.
        self.assertFalse(d['async_magnetometer_nonexpansive_same_M_closed'])
        self.assertFalse(d['common_storage_closed'])
        self.assertFalse(d['ALT_LIVE_PASS'])


if __name__=='__main__':unittest.main()
