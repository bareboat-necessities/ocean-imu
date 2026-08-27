import json
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval, matrix_add, matrix_mul, matrix_transpose
import ou3_p5_first_accel_rotation_gauge as V1G
import ou3_p5_first_accel_rotation_gauge_v2 as G
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3

# Captured at import so the V1 gate stays available for comparison whatever
# else has already installed a backend on the shared module.
V1_SCALAR_AXIS_STRUCTURE = V1G._scalar_axis_structure
HULLED = G._scalar_axis_structure_hulled


def _first_prefix_covariance():
    FULL3._install_backend()
    domain_path = Path(G.DEFAULT_DOMAIN).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    src, _phase = V1G._source_phase_children(1)[0]
    P0 = FULL._initial_covariance(src, domain_path)
    F, Q, _Rstep = FULL._transition_and_Q(src, domain)
    return FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, P0), matrix_transpose(F)), Q))


def _axis_split_first_prefix():
    """Isotropic first prefix whose axis representations differ by one ulp."""
    P = FULL._zero(18, 18)
    for base, value in (((12, 12), 2500.01), ((12, 15), 1e-9), ((15, 15), 0.3)):
        for ax in range(3):
            lo = FULL.down(value) if ax else value
            hi = FULL.up(value) if ax == 2 else value
            P[base[0] + ax][base[1] + ax] = Interval(lo, hi)
            if base == (12, 15):
                P[base[1] + ax][base[0] + ax] = Interval(lo, hi)
    return P


class Ou3P5FirstAccelRotationGaugeV2Tests(unittest.TestCase):
    """V2 relaxes V1's representation gate but keeps its structural gates.

    V2 therefore clears the last-ulp axis widening that stopped V1 and then
    stops on the cross-axis subnormal dust the generic 18x18 interval product
    leaves around the exact structural zero.  V3 is the stage that certifies
    that dust against the source sparsity and runs to a result.
    """

    @classmethod
    def setUpClass(cls):
        cls.Pp = _first_prefix_covariance()

    def test_hull_clears_the_axis_representation_gate_that_stopped_v1(self):
        P = _axis_split_first_prefix()
        with self.assertRaisesRegex(RuntimeError, "lost axis symmetry"):
            V1_SCALAR_AXIS_STRUCTURE(P)
        pss, psa, paw = HULLED(P)
        for value, base in ((pss, (12, 12)), (psa, (12, 15)), (paw, (15, 15))):
            for ax in range(3):
                axis = P[base[0] + ax][base[1] + ax]
                self.assertLessEqual(value.lo, axis.lo)
                self.assertGreaterEqual(value.hi, axis.hi)

    def test_hull_still_stops_on_first_prefix_cross_axis_interval_dust(self):
        with self.assertRaisesRegex(RuntimeError, "gained cross-axis terms"):
            HULLED(self.Pp)

    def test_the_rejected_cross_axis_entries_are_subnormal_arithmetic_dust(self):
        worst = 0.0
        for ai in range(3):
            for aj in range(3):
                if ai == aj:
                    continue
                for a, b in ((12 + ai, 12 + aj), (12 + ai, 15 + aj), (15 + ai, 15 + aj)):
                    z = self.Pp[a][b]
                    self.assertLessEqual(z.lo, 0.0)
                    self.assertGreaterEqual(z.hi, 0.0)
                    worst = max(worst, z.abs_upper())
        self.assertLess(worst, 1e-300)

    def test_hull_rejects_reachable_cross_axis_covariance(self):
        P = _axis_split_first_prefix()
        P[12][13] = Interval(0.25, 0.25)
        with self.assertRaisesRegex(RuntimeError, "gained cross-axis terms"):
            HULLED(P)

    def test_hull_rejects_nonzero_theta_aw_covariance(self):
        P = _axis_split_first_prefix()
        P[0][15] = Interval(0.0, 1e-9)
        with self.assertRaisesRegex(RuntimeError, "theta/a_w covariance is not exactly zero"):
            HULLED(P)


class Ou3P5FirstAccelRotationGaugeV2ContractTests(unittest.TestCase):
    def test_v2_declares_its_own_schema_over_the_v1_stage(self):
        self.assertEqual(G.SCHEMA, 2)
        self.assertNotEqual(G.SCHEMA, V1G.SCHEMA)
        self.assertEqual(G.DEFAULT_DOMAIN, V1G.DEFAULT_DOMAIN)

    def test_validate_rejects_a_weakened_structural_contract(self):
        d = {
            "schema": G.SCHEMA,
            "axis_isotropic_source_intervals_hulled_across_equivalent_axes": False,
            "bit_identical_axis_interval_endpoints_required": True,
            "cross_axis_covariance_still_required_exact_zero": False,
            "theta_aw_covariance_still_required_exact_zero": False,
        }
        failures = G.validate(d)
        self.assertIn("axis-equivalent interval hull is not active", failures)
        self.assertIn("bit-identical interval endpoint gate remains active", failures)
        self.assertIn("cross-axis zero structure was weakened", failures)
        self.assertIn("theta/a_w zero structure was weakened", failures)

    def test_validate_inherits_the_v1_stage_gates(self):
        message = "rotation-gauged innovation still required loose spectral inverse fallback"
        self.assertIn(message, G.validate({"spectral_fallback_inverse_count": 6160}))
        self.assertIn("schema mismatch", G.validate({"schema": V1G.SCHEMA}))


if __name__ == "__main__":
    unittest.main()
