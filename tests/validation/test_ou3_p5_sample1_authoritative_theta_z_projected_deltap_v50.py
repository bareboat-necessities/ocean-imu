import math
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from ou3_interval import Interval
import ou3_p5_sample1_authoritative_theta_z_projected_deltap_v50 as V50


class Sample1AuthoritativeThetaZProjectedDeltaPV50Tests(unittest.TestCase):
    def test_authoritative_witness_schema_and_target_are_frozen(self):
        self.assertEqual(V50.WITNESS, (0, 0, 23))
        self.assertEqual(V50.SCHEMA, 5000)
        self.assertEqual(V50.Q_TARGET, 8.0)

    def test_nominal_first_psd_attitude_component_matrix_uses_sparse_transport(self):
        d = V50._nominal_first_psd_attitude_component_matrix(
            beta_upper=0.5, offdiag_upper=2.0)
        self.assertGreaterEqual(d[0][1], 0.5)
        self.assertGreaterEqual(d[0][2], 1.0)
        self.assertGreaterEqual(d[1][2], 1.0)
        self.assertEqual(d[0][0], 0.0)
        self.assertEqual(d[1][1], 0.0)
        self.assertEqual(d[2][2], 0.0)
        self.assertEqual(d[2][0], d[0][2])
        self.assertEqual(d[2][1], d[1][2])

    def test_abs_congruence_preserves_identity_component_matrix(self):
        z = Interval.point(0.0)
        o = Interval.point(1.0)
        L = [[o, z, z], [z, o, z], [z, z, o]]
        M = [[0.0, 2.0, 3.0], [2.0, 0.0, 4.0], [3.0, 4.0, 0.0]]
        out = V50._abs_congruence_upper(L, M)
        for i in range(3):
            for j in range(3):
                self.assertGreaterEqual(out[i][j], M[i][j])
                self.assertLessEqual(out[i][j], M[i][j] + 1.0e-12)

    def test_theta_z_projection_uses_requested_component_combinations(self):
        d = V50._theta_z_projected_upper(
            dP_zx=1.0, dP_zy=2.0, dP_zz=3.0,
            fy_abs=4.0, fz_abs=5.0)
        self.assertGreaterEqual(
            d["minus_fz_DeltaP_zy_plus_fy_DeltaP_zz_abs_upper"], 22.0)
        self.assertGreaterEqual(d["fz_DeltaP_zx_abs_upper"], 5.0)
        self.assertGreaterEqual(d["minus_fy_DeltaP_zx_abs_upper"], 4.0)
        self.assertGreaterEqual(
            d["theta_z_projected_DeltaP_Htheta_componentwise_upper"],
            math.sqrt(22.0 ** 2 + 5.0 ** 2 + 4.0 ** 2))

    def test_invalid_component_inputs_fail_closed(self):
        with self.assertRaises(ValueError):
            V50._nominal_first_psd_attitude_component_matrix(
                beta_upper=-1.0, offdiag_upper=1.0)
        with self.assertRaises(ValueError):
            V50._theta_z_projected_upper(
                dP_zx=-1.0, dP_zy=1.0, dP_zz=1.0,
                fy_abs=1.0, fz_abs=1.0)


if __name__ == "__main__":
    unittest.main()
