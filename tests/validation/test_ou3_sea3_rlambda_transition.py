import sys
from pathlib import Path
import unittest

TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_sea3_rlambda_transition as RL


class Sea3RLambdaTransitionTest(unittest.TestCase):
    def setUp(self):
        self.gravity = 9.80665
        self.hs_upper = 8.5
        self.gamma = [1.0, 7.0]

    @staticmethod
    def part(h, fp, gamma=1.0, beta=0.25, spread=0.5):
        nu = 1.0 if fp == float("inf") else fp / (1.0 + fp)
        return RL.CompactPartition(h, nu, gamma, beta, spread)

    def test_build_closes_only_machine_readable_outer_relation(self):
        d = RL.build()
        self.assertEqual(RL.validate(d), [])
        self.assertTrue(d["machine_readable_R_lambda_closed"])
        self.assertTrue(d["actual_rate_bounded_R_lambda_subset_Rhat"])
        self.assertFalse(d["rate_constants_fitted_or_invented"])
        self.assertFalse(d["hard_shaping_state_materialized_here"])
        self.assertFalse(d["joint_response_materialized_here"])
        self.assertFalse(d["P3_promoted"])

    def test_nonfixed_multimodal_transition_is_admitted(self):
        a = (
            self.part(2.0, 0.15, 1.0, 0.1, 0.2),
            self.part(1.5, 0.08, 3.3, 0.7, 0.6),
            self.part(0.0, float("inf"), 7.0, 1.0, 1.0),
        )
        b = (
            self.part(1.4, 0.12, 2.2, 0.9, 0.4),
            self.part(1.8, 0.07, 6.0, 0.2, 0.8),
            self.part(0.6, 0.05, 1.0, 0.5, 1.0),
        )
        self.assertTrue(RL.transition_admissible(
            a, b,
            gravity=self.gravity,
            Hs_upper_m=self.hs_upper,
            gamma_range=self.gamma,
        ))
        self.assertNotEqual(a, b)

    def test_exact_partition_energy_coupling_is_enforced(self):
        bad = (
            self.part(6.0, 0.05),
            self.part(6.0, 0.05),
            self.part(1.0, 0.05),
        )
        self.assertFalse(RL.lambda_member(
            bad,
            gravity=self.gravity,
            Hs_upper_m=self.hs_upper,
            gamma_range=self.gamma,
        ))

    def test_peak_steepness_rejects_incompatible_active_partition(self):
        bad = (
            self.part(3.0, 0.20),
            self.part(0.0, 0.1),
            self.part(0.0, 0.1),
        )
        self.assertFalse(RL.lambda_member(
            bad,
            gravity=self.gravity,
            Hs_upper_m=self.hs_upper,
            gamma_range=self.gamma,
        ))

    def test_gamma_interval_is_enforced(self):
        bad = (
            self.part(1.0, 0.1, 7.01),
            self.part(0.0, 0.1),
            self.part(0.0, 0.1),
        )
        self.assertFalse(RL.lambda_member(
            bad,
            gravity=self.gravity,
            Hs_upper_m=self.hs_upper,
            gamma_range=self.gamma,
        ))

    def test_compact_frequency_boundary_is_only_inactive_closure(self):
        inactive = (
            self.part(0.0, float("inf")),
            self.part(1.0, 0.1),
            self.part(0.0, 0.1),
        )
        active = (
            self.part(1.0, float("inf")),
            self.part(0.0, 0.1),
            self.part(0.0, 0.1),
        )
        self.assertTrue(RL.lambda_member(
            inactive,
            gravity=self.gravity,
            Hs_upper_m=self.hs_upper,
            gamma_range=self.gamma,
        ))
        self.assertFalse(RL.lambda_member(
            active,
            gravity=self.gravity,
            Hs_upper_m=self.hs_upper,
            gamma_range=self.gamma,
        ))


if __name__ == "__main__":
    unittest.main()