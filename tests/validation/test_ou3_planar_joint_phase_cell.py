import unittest
from tools.stability.ou3_theorem.planar_joint_phase_cell import calculate

class JointPhaseCellThresholdTest(unittest.TestCase):
    def test_positive_fail_closed_budgets(self):
        f={"result_type":"FINITE DIAGNOSTIC ONLY","parity_blocks":{
            "even":{"relative_Frobenius_gain":0.7426143368769409},
            "odd":{"relative_Frobenius_gain":0.8738212970667967}}}
        m={"all_time_frontend_tube_verified":False,
           "squared_norm_decrement_lower":"27/100000",
           "h_binary32_exact":"5368709/1073741824"}
        out=calculate(f,m,7.024764605642485)
        self.assertIsNone(out["coupling_product_budget"])
        self.assertIsNone(out["symmetric_coupling_norm_budget"])
        self.assertFalse(out["mahony_gain_used_as_mekf_mean_gain"])
        self.assertAlmostEqual(out["weyl_information_perturbation_budget_to_muM_1"],
                               6.024764605642485)
        self.assertTrue(out["generated_reverse_feedback_ports_zero"])
        self.assertEqual(out["remaining_bidirectional_loop"],"MEKF mean <-> covariance through state-dependent H and P-dependent K")
        self.assertFalse(out["joint_phase_cell_forward_invariant"])
        self.assertFalse(out["every_placed_window_magnetic_service_verified"])
        self.assertFalse(out["theorem_closed"])

if __name__=="__main__":unittest.main()
