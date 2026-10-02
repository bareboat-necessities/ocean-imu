import unittest
from tools.stability.ou3_theorem.complete_word_nullspace import compatibility_line_certificate, graph_kernel_certificate
from tools.stability.ou3_theorem.nullspace_strata import audit_service_row_strata, magnetic_service_kernel_certificate, invertible_transport_kernel_rule

class T(unittest.TestCase):
    def test_word_dependent_graph_kernel(self):
        # theta constraint leaves span((1,1)); arbitrary word-dependent A
        # then forces BA=-A theta without changing nullity.
        r=graph_kernel_certificate([[1,-1]],[[2,3],[5,7]])
        self.assertEqual(r["attitude_nullity"],1)
        self.assertEqual(r["graph_stack_nullity"],1)
        self.assertTrue(r["kernel_is_word_dependent_BA_graph"])
        self.assertFalse(r["fixed_global_compatibility_vector_required"])

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

    def test_continuous_magnetic_and_transport_cover(self):
        m=magnetic_service_kernel_certificate(1)
        self.assertTrue(m["continuous_service_row_cover_complete"])
        self.assertEqual(m["magnetic_service_common_kernel_dimension"],0)
        t=invertible_transport_kernel_rule()
        self.assertTrue(t["kernel_dimension_invariant_under_transport"])
        self.assertFalse(t["continuous_transport_subdivision_required"])

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
