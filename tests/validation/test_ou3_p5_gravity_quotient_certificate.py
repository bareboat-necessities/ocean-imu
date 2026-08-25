import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_gravity_quotient_certificate as Q


class Ou3P5GravityQuotientCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = Q.build()

    def test_reduced_detectability_is_strict_and_source_declared(self):
        d = self.d
        self.assertEqual(Q.validate(d), [])
        self.assertEqual(d["P5_GRAVITY_QUOTIENT_REDUCED_DETECTABILITY_CERTIFICATE"], "PASS")
        info = d["reduced_attitude_bias_information"]
        self.assertTrue(info["strict"])
        self.assertTrue(math.isfinite(info["alpha_4_quotient_information_lower"]))
        self.assertGreater(info["alpha_4_quotient_information_lower"], 0.0)
        h = d["gravity_packet_hypothesis"]
        self.assertEqual(h["hypothesis_origin"], "DECLARED_DEPLOYMENT_THEOREM_ASSUMPTION_NOT_REPLAY_FIT")
        self.assertGreater(h["gravity_alignment_cosine_lower"], 0.0)
        self.assertEqual(h["accepted_accelerometer_packets_per_word_lower"], 2)

    def test_axial_bias_is_input_not_fake_gauge_or_removed_filter_state(self):
        d = self.d
        self.assertTrue(d["absolute_yaw_quotiented"])
        self.assertTrue(d["axial_gyro_bias_removed_from_strict_metric"])
        self.assertFalse(d["axial_gyro_bias_removed_from_filter_state"])
        self.assertIn("NEUTRAL_BOUNDED", d["axial_gyro_bias_role"])
        sm = d["schur_metric_policy"]
        self.assertFalse(sm["axial_bias_declared_physical_symmetry"])
        self.assertTrue(sm["source_information_cross_terms_retained"])

    def test_old_zero_dynamics_obstruction_is_consumed_not_ignored(self):
        d = self.d
        self.assertTrue(d["old_yaw_only_full_bias_metric_rejected"])
        self.assertTrue(d["old_yaw_only_counterexample_consumed"])
        self.assertEqual(d["old_yaw_only_counterexample_status"], "PASS")

    def test_translation_stays_in_same_complete_four_S_word(self):
        tr = self.d["translation_word"]
        self.assertEqual(tr["state_order"], ["v", "p", "S", "a_w"])
        self.assertTrue(tr["four_S_complete_chain_required"])
        self.assertTrue(tr["source_complete"])
        self.assertGreater(tr["Q_axis_lambda_min_lower"], 0.0)

    def test_complete_nonlinear_quotient_word_is_not_promoted_yet(self):
        d = self.d
        self.assertEqual(d["P5_UNGAUGED_TIMEOUT_QUOTIENT_WORD_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertIn("b_g_parallel", d["next_obligation"])
        text = (ROOT / "tools" / "ou3_p5_gravity_quotient_certificate.py").read_text(encoding="utf-8")
        self.assertNotIn("ou3_exact_replay", text)
        self.assertNotIn("numpy", text)


if __name__ == "__main__":
    unittest.main()
