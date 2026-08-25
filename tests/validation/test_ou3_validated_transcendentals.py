from decimal import Decimal, getcontext
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_validated_transcendentals as mod

getcontext().prec = 100


def dec(x: float) -> Decimal:
    return Decimal.from_float(float(x))


def dec_sin(x: Decimal, terms: int = 90) -> Decimal:
    total = Decimal(0)
    term = x
    total += term
    x2 = x * x
    for n in range(1, terms):
        term *= -x2 / Decimal((2 * n) * (2 * n + 1))
        total += term
    return total


def dec_cos(x: Decimal, terms: int = 90) -> Decimal:
    total = Decimal(1)
    term = Decimal(1)
    x2 = x * x
    for n in range(1, terms):
        term *= -x2 / Decimal((2 * n - 1) * (2 * n))
        total += term
    return total


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

    def test_sin_cos_enclose_high_precision_reference_on_entire_p4_range(self):
        for x in (-4.0, -3.0, -1.5, -0.01, 0.0, 0.01, 1.5, 3.0, 4.0):
            with self.subTest(x=x):
                X = dec(x)
                s = dec_sin(X)
                c = dec_cos(X)
                S = mod.sin_point(x)
                C = mod.cos_point(x)
                self.assertLessEqual(dec(S.lo), s)
                self.assertGreaterEqual(dec(S.hi), s)
                self.assertLessEqual(dec(C.lo), c)
                self.assertGreaterEqual(dec(C.hi), c)

    def test_sinc_cosc_are_regular_at_zero_and_enclose_reference(self):
        self.assertTrue(mod.sinc_point(0.0).contains(1.0))
        self.assertTrue(mod.cosc_point(0.0).contains(0.5))
        for x in (1.0e-8, 0.01, 0.5, 1.5, 3.0):
            with self.subTest(x=x):
                X = dec(x)
                sinc = dec_sin(X) / X
                cosc = (Decimal(1) - dec_cos(X)) / (X * X)
                S = mod.sinc_point(x)
                C = mod.cosc_point(x)
                self.assertLessEqual(dec(S.lo), sinc)
                self.assertGreaterEqual(dec(S.hi), sinc)
                self.assertLessEqual(dec(C.lo), cosc)
                self.assertGreaterEqual(dec(C.hi), cosc)

    def test_audited_ranges_are_explicit(self):
        with self.assertRaises(ValueError):
            mod.exp_point(0.5000001)
        with self.assertRaises(ValueError):
            mod.expm1_point(-0.5000001)
        with self.assertRaises(ValueError):
            mod.sin_point(4.0000001)
        with self.assertRaises(ValueError):
            mod.cosc_point(-4.0000001)

    def test_proof_module_does_not_call_libm_transcendentals(self):
        text = (ROOT / "tools" / "ou3_validated_transcendentals.py").read_text(encoding="utf-8")
        for dead in ("math.exp(", "math.expm1(", "math.sin(", "math.cos("):
            self.assertNotIn(dead, text)


if __name__ == "__main__":
    unittest.main()
