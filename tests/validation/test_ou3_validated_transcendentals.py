from decimal import Decimal, getcontext
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_validated_transcendentals as mod

getcontext().prec = 90


def dec(x: float) -> Decimal:
    return Decimal.from_float(float(x))


class ValidatedTranscendentalTests(unittest.TestCase):
    def test_exp_and_expm1_enclose_high_precision_reference(self):
        for x in (-0.5, -0.25, -0.01, -0.0004, 0.0, 0.01, 0.5):
            with self.subTest(x=x):
                exact_exp = dec(x).exp()
                exact_em1 = exact_exp - Decimal(1)
                E = mod.exp_point(x)
                M = mod.expm1_point(x)
                self.assertLessEqual(dec(E.lo), exact_exp)
                self.assertGreaterEqual(dec(E.hi), exact_exp)
                self.assertLessEqual(dec(M.lo), exact_em1)
                self.assertGreaterEqual(dec(M.hi), exact_em1)

    def test_ou_positive_kernels_avoid_expm1_cancellation(self):
        for x in (0.0004, 0.01, 0.25, 0.5):
            with self.subTest(x=x):
                X = dec(x)
                exp_neg = (-X).exp()
                exact_pa = X + exp_neg - Decimal(1)
                exact_Sa = X * X / Decimal(2) - X + Decimal(1) - exp_neg
                pa = mod.ou_phi_pa_kernel_point(x)
                Sa = mod.ou_phi_Sa_kernel_point(x)
                self.assertLessEqual(dec(pa.lo), exact_pa)
                self.assertGreaterEqual(dec(pa.hi), exact_pa)
                self.assertLessEqual(dec(Sa.lo), exact_Sa)
                self.assertGreaterEqual(dec(Sa.hi), exact_Sa)
                self.assertGreaterEqual(pa.lo, 0.0)
                self.assertGreaterEqual(Sa.lo, 0.0)

    def test_interval_endpoints_use_monotonicity(self):
        from ou3_interval import Interval
        X = Interval(-0.25, -0.01)
        E = mod.exp_interval(X)
        self.assertLessEqual(dec(E.lo), dec(X.lo).exp())
        self.assertGreaterEqual(dec(E.hi), dec(X.hi).exp())
        K = mod.ou_phi_pa_kernel_interval(Interval(0.01, 0.25))
        exact_lo = dec(0.01) + (-dec(0.01)).exp() - Decimal(1)
        exact_hi = dec(0.25) + (-dec(0.25)).exp() - Decimal(1)
        self.assertLessEqual(dec(K.lo), exact_lo)
        self.assertGreaterEqual(dec(K.hi), exact_hi)

    def test_audited_range_is_explicit(self):
        with self.assertRaises(ValueError):
            mod.exp_point(0.5000001)
        with self.assertRaises(ValueError):
            mod.expm1_point(-0.5000001)

    def test_proof_module_does_not_call_libm_exp(self):
        text = (ROOT / "tools" / "ou3_validated_transcendentals.py").read_text(encoding="utf-8")
        self.assertNotIn("math.exp(", text)
        self.assertNotIn("math.expm1(", text)


if __name__ == "__main__":
    unittest.main()
