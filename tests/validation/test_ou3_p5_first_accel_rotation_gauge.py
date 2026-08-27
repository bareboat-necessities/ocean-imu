import json
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval, matrix_add, matrix_mul, matrix_transpose
import ou3_p5_first_accel_rotation_gauge as G
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3

# Captured at import, before any producer installs a replacement axis backend on
# the module, so these tests exercise the V1 gate exactly as V1 wrote it no
# matter which later stage has already run in the same process.
V1_SCALAR_AXIS_STRUCTURE = G._scalar_axis_structure


def _first_prefix_source_matrices():
    FULL3._install_backend()
    domain_path = Path(G.DEFAULT_DOMAIN).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    src, _phase = G._source_phase_children(1)[0]
    P0 = FULL._initial_covariance(src, domain_path)
    F, Q, _Rstep = FULL._transition_and_Q(src, domain)
    return P0, F, Q


def _first_prefix_covariance():
    """Predicted covariance of the first source/phase child, as build() forms it."""
    P0, F, Q = _first_prefix_source_matrices()
    return FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, P0), matrix_transpose(F)), Q))


def _isotropic_first_prefix(pss: Interval, psa: Interval, paw: Interval):
    P = FULL._zero(18, 18)
    for ax in range(3):
        P[12 + ax][12 + ax] = pss
        P[12 + ax][15 + ax] = psa
        P[15 + ax][12 + ax] = psa
        P[15 + ax][15 + ax] = paw
    return P


class Ou3P5FirstAccelRotationGaugeSourceStructureTests(unittest.TestCase):
    """The gauge's axis isotropy is a property of the source matrices."""

    @classmethod
    def setUpClass(cls):
        cls.P0, cls.F, cls.Q = _first_prefix_source_matrices()

    def test_first_prefix_source_matrices_are_exactly_axis_isotropic(self):
        for M in (self.P0, self.F, self.Q):
            for base in ((12, 12), (12, 15), (15, 15)):
                first = M[base[0]][base[1]]
                for ax in range(3):
                    self.assertEqual(M[base[0] + ax][base[1] + ax], first)

    def test_first_prefix_source_matrices_have_no_cross_axis_linear_entries(self):
        for M in (self.P0, self.F, self.Q):
            for ai in range(3):
                for aj in range(3):
                    if ai == aj:
                        continue
                    for a, b in ((12 + ai, 12 + aj), (12 + ai, 15 + aj), (15 + ai, 15 + aj)):
                        self.assertEqual(M[a][b], Interval.point(0.0))


class Ou3P5FirstAccelRotationGaugeTests(unittest.TestCase):
    """V1 is retained as the historical stage that stops at its own gate.

    Its structural requirements are correct but it compares the three
    *interval representations* of numerically identical axis expressions, so
    the generic 18x18 product's last-ulp widening stops it before it can
    evaluate a single correction cell.  V2 relaxes that to an axis hull and V3
    certifies the source sparsity; both are exercised by their own tests.
    """

    @classmethod
    def setUpClass(cls):
        cls.Pp = _first_prefix_covariance()

    def test_gate_stops_on_bit_identical_axis_interval_representations(self):
        with self.assertRaisesRegex(RuntimeError, "lost axis symmetry"):
            V1_SCALAR_AXIS_STRUCTURE(self.Pp)

    def test_the_rejected_axis_representations_differ_only_by_arithmetic_noise(self):
        for base in ((12, 12), (12, 15), (15, 15)):
            axes = [self.Pp[base[0] + ax][base[1] + ax] for ax in range(3)]
            scale = max(abs(x.lo) for x in axes) + max(abs(x.hi) for x in axes)
            self.assertGreater(scale, 0.0)
            for attr in ("lo", "hi"):
                values = [getattr(x, attr) for x in axes]
                self.assertLessEqual((max(values) - min(values)) / scale, 1e-12)

    def test_gate_accepts_an_exactly_isotropic_first_prefix(self):
        pss = Interval(1.0, 2.0)
        psa = Interval(-0.5, 0.5)
        paw = Interval(3.0, 4.0)
        self.assertEqual(
            V1_SCALAR_AXIS_STRUCTURE(_isotropic_first_prefix(pss, psa, paw)),
            (pss, psa, paw),
        )

    def test_gate_rejects_reachable_cross_axis_covariance(self):
        P = _isotropic_first_prefix(Interval(1.0, 2.0), Interval(-0.5, 0.5), Interval(3.0, 4.0))
        P[12][13] = Interval(0.25, 0.25)
        with self.assertRaisesRegex(RuntimeError, "gained cross-axis terms"):
            V1_SCALAR_AXIS_STRUCTURE(P)

    def test_gate_rejects_nonzero_theta_aw_covariance(self):
        P = _isotropic_first_prefix(Interval(1.0, 2.0), Interval(-0.5, 0.5), Interval(3.0, 4.0))
        P[0][15] = Interval(0.0, 1e-9)
        with self.assertRaisesRegex(RuntimeError, "theta/a_w covariance is not exactly zero"):
            V1_SCALAR_AXIS_STRUCTURE(P)


class Ou3P5FirstAccelRotationGaugeValidateTests(unittest.TestCase):
    def test_validate_reports_a_surviving_spectral_fallback_as_an_obligation(self):
        message = "rotation-gauged innovation still required loose spectral inverse fallback"
        self.assertIn(message, G.validate({"spectral_fallback_inverse_count": 6160}))
        self.assertNotIn(message, G.validate({"spectral_fallback_inverse_count": 0}))

    def test_validate_keeps_the_stage_out_of_the_promotion_path(self):
        d = {
            "whole_word_promoted_here": True,
            "N_H_words_set_here": True,
            "deployed_correction_limit_increased": True,
            "deployed_correction_limit_rad": 8.0,
        }
        failures = G.validate(d)
        self.assertIn("whole_word_promoted_here is not false", failures)
        self.assertIn("N_H_words_set_here is not false", failures)
        self.assertIn("deployed_correction_limit_increased is not false", failures)
        self.assertIn("deployed correction range changed", failures)


if __name__ == "__main__":
    unittest.main()
