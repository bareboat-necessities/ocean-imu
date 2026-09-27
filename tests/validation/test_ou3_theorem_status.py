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
        self.assertTrue(o["historical_AG_readout_covariance_implication"])
        self.assertFalse(o["uniform_historical_AG_readout_action"])
        self.assertTrue(o["source_covariance_projection_guard"])
        self.assertFalse(o["certified_tail_prefix_retention"])
        for key in ("full_state_magnetic_information_lifting",
                    "constructive_full_A21_mu_rho_enclosure",
                    "source_uniform_A21_linear_dissipativity",
                    "constructive_fixed_coordinate_A21_mu_enclosure"):
            self.assertFalse(o[key])
    def test_regime_and_pivot_algebra_does_not_promote_stability(self):
        r=status_report()
        self.assertEqual(r["physical_regimes"],["STILL","TRANSITION","MOVING"])
        self.assertFalse(r["shipping_mode_switch_enabled"])
        self.assertTrue(r["exact_rest_detector_obstruction"]["identical_rest_motion_IMU_histories"])
        o=r["obligations"]
        for key in ("stationary_observability_structure", "stationary_gyro_information_bound",
                    "finite_transition_prefix_composition", "conditional_two_group_six_column_bound",
                    "singular_floor_to_all_six_historical_pivots"):
            self.assertTrue(o[key])
        for key in ("stationary_compatible_class_practical_stability", "sound_runtime_regime_certification",
                    "finite_regime_detection_and_transition_retention", "recurring_regime_storage_budget",
                    "temporal_margins_to_six_historical_pivots", "full_21_covariance_upper"):
            self.assertFalse(o[key])

    def test_committed_status_matches_code(self):
        p=ROOT/"reports/results/ou3_stability/theorem-status.json"
        self.assertEqual(json.loads(p.read_text()),status_report())

if __name__=="__main__": unittest.main()
