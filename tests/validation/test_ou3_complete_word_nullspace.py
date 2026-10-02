import unittest
from tools.stability.ou3_theorem.complete_word_nullspace import compatibility_line_certificate
from tools.stability.ou3_theorem.nullspace_strata import audit_service_row_strata

class T(unittest.TestCase):
    def test_exact_line(self):
        blocks=[("S",[[1,0,0]]),("mag",[[0,1,-1]])]
        r=compatibility_line_certificate(blocks,[0,1,1])
        self.assertEqual(r["rank"],2); self.assertEqual(r["nullity"],1)
        self.assertTrue(r["kernel_equals_compatibility_line"])

    def test_hidden_second_kernel_fails_closed(self):
        blocks=[("S",[[1,0,0,0]]),("mag",[[0,1,-1,0]])]
        r=compatibility_line_certificate(blocks,[0,1,1,0])
        self.assertEqual(r["nullity"],2)
        self.assertFalse(r["kernel_equals_compatibility_line"])

    def test_service_row_stratum_uses_no_gain(self):
        st={"name":"regular","process_rows":[[1,0,0]],"S_rows":[],
            "mag_rows":[[0,1,-1]],"acc_rows":[],"gyro_rows":[],"ba_rows":[],
            "compatibility_line":[0,1,1]}
        r=audit_service_row_strata([st])
        self.assertTrue(r["all_strata_verified"])
        self.assertFalse(r["uses_K_or_innovation_covariance"])
        self.assertFalse(r["callback_pattern_enumeration_used"])

    def test_wrong_line_fails_closed(self):
        blocks=[("S",[[1,0]]),("mag",[[0,1]])]
        r=compatibility_line_certificate(blocks,[0,1])
        self.assertFalse(r["compatibility_line_annihilated"])
        self.assertFalse(r["kernel_equals_compatibility_line"])

if __name__=="__main__": unittest.main()
