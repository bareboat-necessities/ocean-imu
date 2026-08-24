import copy
import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
SPEC = importlib.util.spec_from_file_location(
    "ou3_stability_theorem_certificate",
    ROOT / "tools" / "ou3_stability_theorem_certificate.py",
)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)


class Ou3StabilityTheoremCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = MOD.build()

    def test_composed_certificate_passes_without_replay_evidence(self):
        d = self.payload
        self.assertEqual(d["schema"], 2)
        self.assertEqual(d["status"], "PASS_CONDITIONAL_LOCAL_ISS", d["failures"])
        self.assertEqual(MOD.validate(d), [])
        self.assertFalse(d["sampled_evidence_used"])
        self.assertFalse(d["trajectory_fit_used"])
        self.assertTrue(all(d["theorem_obligations"].values()))

    def test_both_fixed_dimensional_live_modes_are_ues(self):
        linear = self.payload["linearized_normal_live"]
        self.assertEqual(linear["H_dimension"], 18)
        self.assertEqual(linear["A_dimension"], 21)
        self.assertTrue(linear["bounded_stabilizing_Riccati"])
        self.assertTrue(linear["uniform_exponential_stability"])
        q = self.payload["quantitative_anchors"]
        self.assertGreater(q["H_prediction_Q_lambda_min_lower"], 0.0)
        self.assertGreater(q["A_prediction_Q_lambda_min_lower"], 0.0)
        self.assertGreater(q["translation_information_lower"], 0.0)
        self.assertGreater(q["vector_alpha_6_information_lower"], 0.0)

    def test_all_stable_tails_are_strict(self):
        q = self.payload["quantitative_anchors"]
        self.assertGreater(q["stable_aw_alpha_upper"], 0.0)
        self.assertLess(q["stable_aw_alpha_upper"], 1.0)
        ba = q["active_accelerometer_bias"]
        self.assertTrue(ba["pass"], ba.get("failure"))
        self.assertGreater(ba["alpha_interval"][0], 0.0)
        self.assertLess(ba["alpha_interval"][1], 1.0)

    def test_full_heading_packet_hypotheses_are_explicit_not_empirical(self):
        pe = self.payload["persistent_excitation_operating_envelope"]
        self.assertTrue(pe["persistent_excitation_is_theorem_hypothesis"])
        self.assertTrue(pe["accepted_accelerometer_packet_at_vector_times_required"])
        self.assertTrue(pe["accepted_magnetometer_consecutive_pair_required"])
        self.assertTrue(pe["measurement_gate_margin_required"])
        self.assertGreater(pe["specific_force_norm_lower_mps2"], 0.0)
        self.assertGreater(pe["magnetic_vector_norm_lower_uT"], 0.0)
        self.assertGreater(pe["vector_sine_separation_lower"], 0.0)
        self.assertLess(pe["vector_sine_separation_lower"], 1.0)
        self.assertEqual(len(pe["packet_gap_s"]), 2)
        self.assertGreater(pe["packet_gap_s"][0], 0.0)
        self.assertIn("NOT_INFERRED", pe["qualification"])

    def test_nonlinear_theorem_requires_branch_regular_gate_margin(self):
        nonlinear = self.payload["nonlinear_normal_live"]
        limits = self.payload["scope_limits"]
        self.assertTrue(nonlinear["branch_regular_source_word_required"])
        self.assertTrue(limits["measurement_gate_boundary_points_not_certified"])
        self.assertTrue(limits["consecutive_accepted_mag_pair_required_for_vector_uco"])
        self.assertTrue(limits["permanent_or_unbounded_measurement_rejection_not_certified"])

    def test_certificate_does_not_overclaim_numeric_deployment_basin(self):
        nonlinear = self.payload["nonlinear_normal_live"]
        limits = self.payload["scope_limits"]
        claims = self.payload["claim_separation"]
        self.assertTrue(nonlinear["local_iss"])
        self.assertTrue(nonlinear["nonzero_neighborhood_exists"])
        self.assertFalse(nonlinear["explicit_numeric_basin_radius_produced"])
        self.assertTrue(limits["full_heading_requires_persistent_excitation"])
        self.assertEqual(
            claims["numerical_source_complete_deployment_funnel"],
            "NOT_ESTABLISHED_BY_THIS_PRODUCER",
        )

    def test_validator_rejects_hidden_PE_branch_or_sampled_promotion(self):
        d = copy.deepcopy(self.payload)
        d["persistent_excitation_operating_envelope"][
            "persistent_excitation_is_theorem_hypothesis"
        ] = False
        self.assertTrue(any("PE" in x for x in MOD.validate(d)))

        d = copy.deepcopy(self.payload)
        d["persistent_excitation_operating_envelope"][
            "accepted_magnetometer_consecutive_pair_required"
        ] = False
        self.assertTrue(any("consecutive" in x for x in MOD.validate(d)))

        d = copy.deepcopy(self.payload)
        d["persistent_excitation_operating_envelope"]["measurement_gate_margin_required"] = False
        self.assertTrue(any("gate-margin" in x for x in MOD.validate(d)))

        d = copy.deepcopy(self.payload)
        d["nonlinear_normal_live"]["branch_regular_source_word_required"] = False
        self.assertTrue(any("branch-regular" in x for x in MOD.validate(d)))

        d = copy.deepcopy(self.payload)
        d["sampled_evidence_used"] = True
        self.assertTrue(any("sampled" in x for x in MOD.validate(d)))

    def test_periodic_aw_sync_is_exactly_nonexpansive(self):
        aw = self.payload["periodic_aw_covariance_sync"]
        self.assertTrue(aw["pass"])
        self.assertEqual(aw["proof_mode"], "PSD_NONEXPANSIVE")
        self.assertLessEqual(aw["jump_gain_upper"], 1.0)


if __name__ == "__main__":
    unittest.main()
