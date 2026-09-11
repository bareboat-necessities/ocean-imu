"""Independent exact-arithmetic checks of sampled BRMM statistics and scope."""
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))


def module(path):
    spec = importlib.util.spec_from_file_location(Path(path).stem, ROOT/path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


audit = module("tests/kalman_ou_iii/ou3_brmm_audit.py")
contract = module("tools/stability/ou3_brmm_contract.py")


class BrmmTests(unittest.TestCase):
    def assert_encloses(self, lo, exact, hi):
        self.assertLessEqual(F(float(lo)), exact)
        self.assertLessEqual(exact, F(float(hi)))

    def test_every_window_against_exact_rational_statistics(self):
        a = np.random.default_rng(17).normal(size=(51, 3))
        w = audit.windows(a, sample_hz=4)
        for start in range(12):
            values = [[F(float(x)) for x in row] for row in a[start:start+40]]
            impulse = [sum(row[k] for row in values)/4 for k in range(3)]
            energy = sum(x*x for row in values for x in row)/4-sum(x*x for x in impulse)/10
            self.assert_encloses(w["energy_lo"][start], energy, w["energy_hi"][start])
            for k in range(3):
                self.assert_encloses(w["impulse_lo"][start, k], impulse[k], w["impulse_hi"][start, k])

    def test_every_positive_lobe_is_an_exact_disjoint_witness(self):
        a = np.zeros((105, 3))
        a[:, 0] = np.sin(2*np.pi*np.arange(105)/8)
        w = audit.windows(a, sample_hz=8)
        self.assertTrue(np.all(w["lobe_J_lower"] > 0))
        for start, index in enumerate(w["lobe_witness_index"]):
            s, t, u, v = w["lobe_candidate_sample_offsets"][index]
            self.assertLessEqual(t, u)
            first = sum(F(float(x)) for x in a[start+s:start+t, 0])/8
            second = sum(F(float(x)) for x in a[start+u:start+v, 0])/8
            bound = F(float(w["lobe_J_lower"][start]))**2
            self.assertGreaterEqual(first*first, bound)
            self.assertGreaterEqual(second*second, bound)
            self.assertLessEqual(first*second, -bound)

    def test_constant_acceleration_exposes_quiet_branch_loophole(self):
        # Small DC passes AC and 2 m/s window cap; the bounded primitive excludes it.
        a = np.tile([.1, 0., 0.], (200, 1))
        w = audit.windows(a, sample_hz=4)
        s = audit.stats(w, .03)
        self.assertEqual(s["Q_windows"], 161)
        self.assertEqual(s["Q_impulse_cap_2_violations"], 0)
        self.assertTrue(np.all(w["impulse_norm_lo"] > .99))
        d = contract.build()
        self.assertTrue(d["bounded_velocity_primitive_required"])
        self.assertTrue(d["nonzero_fixed_physical_acceleration_DC_prohibited_in_Q_and_O"])
        self.assertFalse(d["fixed_window_impulse_cap_alone_excludes_all_DC"])
        self.assertFalse(d["constant_measurement_bias_is_marine_motion"])

    def test_zero_motion_requires_no_global_excitation(self):
        w = audit.windows(np.zeros((40, 3)), sample_hz=4)
        self.assertEqual(audit.stats(w, .03)["Q_windows"], 1)
        self.assertEqual(w["lobe_witness_index"][0], -1)
        self.assert_encloses(w["energy_lo"][0], F(0), w["energy_hi"][0])
        self.assert_encloses(w["tv_lo"][0], F(0), w["tv_hi"][0])

    def test_long_period_one_sided_window_can_be_oscillatory(self):
        a = np.zeros((80, 3))
        a[:, 0] = np.cos(2*np.pi*.02*np.arange(80)/8)
        s = audit.stats(audit.windows(a, sample_hz=8), .03)
        self.assertEqual(s["O_windows"], 1)
        self.assertEqual(s["O_lobe_search_unresolved"], 1)
        self.assertFalse(s["lobe_search_failure_proves_nonexistence"])

    def test_total_variation_internal_jumps_exact(self):
        a = np.zeros((43, 3))
        a[:, 0] = (np.arange(43) % 3)/4
        w = audit.windows(a, sample_hz=4)
        for start in range(4):
            exact = sum(abs(F(float(y))-F(float(x))) for x, y in zip(
                a[start:start+39, 0], a[start+1:start+40, 0]))
            self.assert_encloses(w["tv_lo"][start], exact, w["tv_hi"][start])

    def test_invalid_samples_and_clock_rejected(self):
        for a in (np.zeros((39, 3)), np.full((40, 3), np.nan), np.zeros((40, 2))):
            with self.assertRaises(ValueError):
                audit.windows(a, sample_hz=4)
        with self.assertRaises(ValueError):
            audit.windows(np.zeros((50, 3)), sample_hz=5)

    def test_qualified_primitives_do_not_promote_remaining_source_obligations(self):
        d = contract.build()
        self.assertEqual(contract.validate(d), [])
        self.assertFalse(d["spectral_membership_required"])
        self.assertEqual(d["P3_delta"], 1e-18)
        physical = d["physical_constants"]
        primitives = d["uniform_primitive_qualification"]["uniform_physical_primitives"]
        self.assertTrue(d["BRMM_PRIMITIVE_QUALIFICATION_PASS"])
        self.assertEqual(physical["V_m"], primitives["V_m_norm_upper_mps"])
        self.assertEqual(physical["P_m"], primitives["P_m_norm_upper_m"])
        self.assertEqual(physical["A_m"], primitives["acceleration_norm_upper_mps2"])
        self.assertGreater(physical["V_m"], 0)
        self.assertGreater(physical["P_m"], 0)
        self.assertTrue(all(physical[key] is None for key in
                            ("T_R", "E_q", "V_R", "E_min", "J_min", "chi", "V_a", "S_m")))
        self.assertFalse(d["sampled_audit_sets_physical_theorem_constants"])
        for key in ("V_m", "P_m", "A_m", "S_m"):
            changed = {**physical, key: 1.0}
            self.assertTrue(contract.validate({**d, "physical_constants": changed}))
        for gate in ("BRMM_SOURCE_UNIFORM_PASS", "source_declaration_is_P3_certificate", "BRMM_P4_MOTION_PASS",
                     "BRMM_P5_MOTION_MAY_START"):
            self.assertTrue(contract.validate({**d, gate: True}))


if __name__ == "__main__":
    unittest.main()
