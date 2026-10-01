import unittest
from tools.stability.ou3_theorem.local_defect_composition import certificate, verify_endpoint

class T(unittest.TestCase):
    def test_exact_local_composition(self):
        r=certificate()
        self.assertTrue(r["synthetic_endpoint_identity_exact"])
        self.assertFalse(r["endpoint_residual_used_as_input"])
        self.assertFalse(r["native_literal_boundary_export_complete"])

    def test_nontrivial_three_step(self):
        b=[
          {"A":[["1"]],"e_before":[["2"]],"e_after":[["5"]]},
          {"A":[["2"]],"e_before":[["5"]],"e_after":[["9"]]},
          {"A":[["3"]],"e_before":[["9"]],"e_after":[["28"]]}]
        self.assertTrue(verify_endpoint(b)["endpoint_identity_exact"])

if __name__=="__main__":
    unittest.main()
