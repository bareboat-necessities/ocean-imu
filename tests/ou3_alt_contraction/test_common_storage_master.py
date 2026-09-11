import unittest
import numpy as np

from tools.stability.ou3_alt_contraction import common_storage_master as C


class CommonStorageMasterTest(unittest.TestCase):
    def test_joint24_selector_is_exactly_bias_reference_coordinates(self):
        S=C.joint24_neutral_selector(); Z=C.joint24_motion_injection()
        self.assertEqual(S.shape,(6,24)); self.assertEqual(Z.shape,(24,18))
        np.testing.assert_array_equal(S@Z,np.zeros((6,18)))
        np.testing.assert_array_equal(S[:,18:24],np.eye(6))
        np.testing.assert_array_equal(Z[:18,:],np.eye(18))
        np.testing.assert_array_equal(Z[18:24,:],np.zeros((6,18)))

    def test_neutral_coordinate_is_closed_by_finite_supply(self):
        # x1 contracts; x2 is exactly held.  Strict homogeneous contraction is
        # impossible for rho<1, but a bounded x2 supply closes the full SPD
        # storage inequality exactly as required for H18.
        A=np.diag([0.5,1.0]); M=np.eye(2); rho=0.9
        R=C.unsupplied_storage_form(A,M,rho,(1,))
        self.assertLess(float(np.linalg.eigvalsh(R)[-1]),0.0)
        d=C.coordinate_supply_completion(A,M,rho,(1,),margin=1e-12)
        self.assertTrue(d['feasible']); self.assertGreater(d['beta_scalar'],0.0)
        self.assertLess(d['completed_lambda_max'],0.0)

    def test_unsupplied_neutral_direction_cannot_be_hidden_by_supply(self):
        # x1 is neutral but only x2 is admitted as supply.  Finsler restriction
        # correctly fails; no Beta on x2 can manufacture contraction of x1.
        A=np.eye(2); M=np.eye(2); rho=0.9
        d=C.coordinate_supply_completion(A,M,rho,(1,),margin=1e-12)
        self.assertFalse(d['feasible'])
        self.assertGreater(d['projected_lambda_max'],0.0)

    def test_cross_coupling_uses_schur_completion_not_diagonal_handwave(self):
        A=np.array([[0.5,0.7],[0.0,1.0]])
        M=np.array([[2.0,0.3],[0.3,1.0]])
        d=C.coordinate_supply_completion(A,M,0.95,(1,),margin=1e-10)
        self.assertTrue(d['feasible'])
        self.assertLess(d['completed_lambda_max'],0.0)

    def test_bounded_additive_supply_multiplier(self):
        self.assertAlmostEqual(C.bounded_additive_supply_multiplier(0.5,0.75),3.0)
        self.assertEqual(C.bounded_additive_supply_multiplier(0.0,0.5),1.0)
        with self.assertRaises(ValueError): C.bounded_additive_supply_multiplier(0.9,0.8)

    def test_phase1_bound_reduction_remains_fail_closed(self):
        d=C.build(); f=C.validate(d)
        self.assertEqual(f,[])
        self.assertTrue(d['phase1_storage_search_allowed_consumed'])
        self.assertFalse(d['common_M_source_uniform_search_closed'])
        self.assertFalse(d['source_uniform_outward_projected_LDLT_closed'])
        self.assertFalse(d['ALT_LIVE_PASS'])


if __name__=='__main__':
    unittest.main()
