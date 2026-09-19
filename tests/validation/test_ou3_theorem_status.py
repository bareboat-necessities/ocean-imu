import json,sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.theorem_status import status_report

class TheoremStatusTests(unittest.TestCase):
    def test_single_architecture_is_fail_closed(self):
        r=status_report()
        self.assertEqual(r["principal_assumptions"],["MARINE MOTION","IMU BIAS","MAGNETIC SERVICE"])
        self.assertFalse(r["parallel_no_magnetometer_stability_path"])
        self.assertFalse(r["theorem_closed"])
        self.assertFalse(r["regional_practical_stability_claimed"])
        self.assertIn("H18-to-A21 release",r["proof_path"])
    def test_full_state_claims_require_more_than_partial_lin_and_service(self):
        o=status_report()["obligations"]
        self.assertTrue(o["LIN_endpoint_matrix_path_action"])
        self.assertTrue(o["complete_word_covariance_energy_identity"])
        for key in ("full_state_magnetic_information_lifting",
                    "constructive_full_A21_mu_rho_enclosure",
                    "source_uniform_A21_linear_dissipativity",
                    "constructive_fixed_coordinate_A21_mu_enclosure"):
            self.assertFalse(o[key])
    def test_committed_status_matches_code(self):
        p=ROOT/"reports/results/ou3_stability/theorem-status.json"
        self.assertEqual(json.loads(p.read_text()),status_report())

if __name__=="__main__": unittest.main()
