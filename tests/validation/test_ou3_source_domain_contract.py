import importlib.util
import math
from pathlib import Path
import struct
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location(
    "ou3_source_domain_contract", ROOT / "tools" / "ou3_source_domain_contract.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def f32(value):
    return struct.unpack(">f", struct.pack(">f", float(value)))[0]


def approximate_qd(h, tau, sigma2, n=4096):
    dr = h / n
    out = [[0.0] * 4 for _ in range(4)]
    qc = 2.0 * sigma2 / tau
    for k in range(n):
        r = (k + 0.5) * dr
        x = r / tau
        a = math.exp(-x)
        g = (
            tau * (1.0 - a),
            tau * tau * (x + math.expm1(-x)),
            tau ** 3 * (0.5 * x * x - x - math.expm1(-x)),
            a,
        )
        for i in range(4):
            for j in range(4):
                out[i][j] += qc * g[i] * g[j] * dr
    return out


class SourceDomainContractTests(unittest.TestCase):
    def test_contract_uses_shipping_clamps_and_keeps_theorem_unpromoted(self):
        d = mod.build(mod.DEFAULT_HEADER)
        self.assertEqual(d["schema"], 2)
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertTrue(d["source_complete_parameter_domain"])
        self.assertFalse(d["validated_arithmetic"])
        self.assertFalse(d["outward_rounded"])
        self.assertEqual(
            d["implementation_scalar_semantics"]["type"], "IEEE754_BINARY32"
        )
        self.assertEqual(
            d["continuous_parameters"]["tau_aw_s"], [f32(0.02), f32(12.0)]
        )
        self.assertEqual(
            d["continuous_parameters"]["sigma_aw_mps2"], [f32(0.05), f32(6.0)]
        )
        self.assertEqual(set(d["discrete_source_branches"]["mode"]), {"H", "A"})

    def test_constexpr_arithmetic_rounds_as_binary32_after_each_operation(self):
        text = """
        constexpr float A = 0.1f;
        constexpr float B = A + A;
        constexpr float C = B + A;
        constexpr float DT = 1.0f / 200.0f;
        """
        a = f32(0.1)
        b = f32(a + a)
        c = f32(b + a)
        dt = f32(f32(1.0) / f32(200.0))
        self.assertEqual(mod.parse_const(text, "A"), a)
        self.assertEqual(mod.parse_const(text, "B"), b)
        self.assertEqual(mod.parse_const(text, "C"), c)
        self.assertEqual(mod.parse_const(text, "DT"), dt)
        self.assertNotEqual(dt, 0.005)

    def test_validated_parameter_box_outwardly_contains_every_source_endpoint(self):
        d = mod.build(mod.DEFAULT_HEADER)
        box = d["validated_parameter_box"]
        self.assertTrue(box["validated_arithmetic"])
        self.assertTrue(box["outward_rounded"])
        self.assertEqual(box["theorem_promotion"], "NOT_ESTABLISHED")
        self.assertFalse(box["continuous_word_enclosed"])
        self.assertFalse(box["nonlinear_word_enclosed"])

        for name, source_bounds in d["continuous_parameters"].items():
            lo, hi = box["continuous_parameters"][name]
            self.assertLess(lo, source_bounds[0])
            self.assertGreater(hi, source_bounds[1])
            self.assertEqual(lo, math.nextafter(source_bounds[0], -math.inf))
            self.assertEqual(hi, math.nextafter(source_bounds[1], math.inf))

        for name, source_value in d["timing_constants_s"].items():
            lo, hi = box["timing_constants_s"][name]
            self.assertLess(lo, source_value)
            self.assertGreater(hi, source_value)
            self.assertEqual(lo, math.nextafter(source_value, -math.inf))
            self.assertEqual(hi, math.nextafter(source_value, math.inf))

    def test_deployment_step_domain_distinguishes_type_bound_from_supported_guard(self):
        d = mod.build(mod.DEFAULT_HEADER)
        step = d["accepted_update_step_domain_s"]
        self.assertGreater(step["lower_closed"], 0.0)
        self.assertTrue(math.isfinite(step["upper_closed"]))
        self.assertGreater(step["upper_closed"], 1.0e30)
        self.assertTrue(step["type_level_finite_upper_bound"])
        self.assertFalse(step["operational_safety_upper_guard"])
        self.assertFalse(step["proof_usable_supported_upper_bound"])
        self.assertIn(
            "NO_PROOF_USABLE_ACCEPTED_DT_GUARD",
            d["validated_ou_primitive_backend"]["theorem_promotion"],
        )
        self.assertFalse(
            d["validated_ou_primitive_backend"]
            ["process_covariance_shipping_float_path_enclosed"]
        )

    def test_rational_taylor_ou_primitive_box_contains_direct_evaluations(self):
        box = mod.validated_ou_primitives((0.001, 0.25), (0.02, 12.0))
        self.assertTrue(box["validated_arithmetic"])
        self.assertTrue(box["outward_rounded"])
        self.assertLessEqual(box["alpha"][0], box["alpha"][1])
        self.assertLessEqual(box["phi_pa_s2"][0], box["phi_pa_s2"][1])
        self.assertLessEqual(box["phi_Sa_s3"][0], box["phi_Sa_s3"][1])

        for h in (0.001, 0.005, 0.05, 0.25):
            for tau in (0.02, 0.1, 1.0, 12.0):
                x = h / tau
                alpha = math.exp(-x)
                phi_pa = tau * tau * (x + math.expm1(-x))
                phi_sa = tau ** 3 * (0.5 * x * x - x - math.expm1(-x))
                self.assertLessEqual(box["alpha"][0], alpha)
                self.assertGreaterEqual(box["alpha"][1], alpha)
                self.assertLessEqual(box["phi_pa_s2"][0], phi_pa)
                self.assertGreaterEqual(box["phi_pa_s2"][1], phi_pa)
                self.assertLessEqual(box["phi_Sa_s3"][0], phi_sa)
                self.assertGreaterEqual(box["phi_Sa_s3"][1], phi_sa)

    def test_validated_axis_transition_contains_deployed_formula_grid(self):
        box = mod.validated_phi_axis4((0.001, 0.25), (0.02, 12.0))
        self.assertTrue(box["validated_arithmetic"])
        self.assertTrue(box["outward_rounded"])
        self.assertEqual(box["state_order"], ["v", "p", "S", "a_w"])
        M = box["Phi_interval"]

        def contains(i, j, value):
            self.assertLessEqual(M[i][j][0], value, (i, j, value, M[i][j]))
            self.assertGreaterEqual(M[i][j][1], value, (i, j, value, M[i][j]))

        for h in (0.001, 0.005, 0.05, 0.25):
            for tau in (0.02, 0.1, 1.0, 12.0):
                x = h / tau
                alpha = math.exp(-x)
                phi_va = tau * (1.0 - alpha)
                phi_pa = tau * tau * (x + math.expm1(-x))
                phi_sa = tau ** 3 * (0.5 * x * x - x - math.expm1(-x))
                exact = (
                    (1.0, 0.0, 0.0, phi_va),
                    (h, 1.0, 0.0, phi_pa),
                    (0.5 * h * h, h, 1.0, phi_sa),
                    (0.0, 0.0, 0.0, alpha),
                )
                for i in range(4):
                    for j in range(4):
                        contains(i, j, exact[i][j])

    def test_validated_process_covariance_contains_midpoint_quadrature_samples(self):
        box = mod.validated_qd_axis4_kernel(
            (0.001, 0.05), (0.02, 2.0), (0.0025, 4.0), cells=12
        )
        self.assertTrue(box["validated_arithmetic"])
        self.assertTrue(box["mathematical_integral_enclosed"])
        self.assertFalse(box["shipping_binary32_closed_form_enclosed"])
        self.assertFalse(box["shipping_psd_cleanup_enclosed"])
        Q = box["Qd_interval"]

        for h, tau, sigma2 in (
            (0.001, 0.02, 0.0025),
            (0.005, 0.1, 0.25),
            (0.02, 1.0, 1.0),
            (0.05, 2.0, 4.0),
        ):
            q = approximate_qd(h, tau, sigma2)
            for i in range(4):
                for j in range(4):
                    self.assertLessEqual(Q[i][j][0], q[i][j] + 1e-15)
                    self.assertGreaterEqual(Q[i][j][1] + 1e-15, q[i][j])

    def test_validated_covariance_prediction_contains_direct_point_prediction(self):
        P0 = [
            [[0.5, 0.5], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
            [[0.0, 0.0], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0]],
            [[0.0, 0.0], [0.0, 0.0], [2.0, 2.0], [0.0, 0.0]],
            [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.1, 0.1]],
        ]
        h, tau, sigma2 = 0.005, 1.0, 0.25
        out = mod.validated_covariance_predict_axis4(
            P0, (h, h), (tau, tau), (sigma2, sigma2), cells=12
        )
        self.assertTrue(out["validated_arithmetic"])
        self.assertTrue(out["mathematical_prediction_enclosed"])
        self.assertFalse(out["shipping_covariance_prediction_enclosed"])
        Pbox = out["P_predicted_interval"]

        x = h / tau
        a = math.exp(-x)
        Phi = [
            [1.0, 0.0, 0.0, tau * (1.0 - a)],
            [h, 1.0, 0.0, tau * tau * (x + math.expm1(-x))],
            [0.5 * h * h, h, 1.0,
             tau ** 3 * (0.5 * x * x - x - math.expm1(-x))],
            [0.0, 0.0, 0.0, a],
        ]
        Ppoint = [[P0[i][j][0] for j in range(4)] for i in range(4)]
        AP = [[sum(Phi[i][k] * Ppoint[k][j] for k in range(4))
               for j in range(4)] for i in range(4)]
        pred = [[sum(AP[i][k] * Phi[j][k] for k in range(4))
                 for j in range(4)] for i in range(4)]
        q = approximate_qd(h, tau, sigma2)
        for i in range(4):
            for j in range(4):
                value = pred[i][j] + q[i][j]
                self.assertLessEqual(Pbox[i][j][0], value + 1e-15)
                self.assertGreaterEqual(Pbox[i][j][1] + 1e-15, value)

    def test_contract_names_every_hybrid_transition_required_for_deployment(self):
        d = mod.build(mod.DEFAULT_HEADER)
        self.assertEqual(
            set(d["hybrid_obligations"]),
            {
                "startup_handoff",
                "held_to_active",
                "magnetic_lock",
                "magnetic_regauge_refinement",
                "tilt_reset",
                "tilt_relock",
                "cooldown_reentry",
                "periodic_aw_covariance_sync",
            },
        )
        self.assertEqual(
            d["periodic_aw_covariance_sync_proof"]["required_mode"],
            "PSD_NONEXPANSIVE",
        )


if __name__ == "__main__":
    unittest.main()
