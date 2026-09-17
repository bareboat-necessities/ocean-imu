"""Every literal joint24 event path must compose, not just the ones a word used.

The joint24/joint27 lift differentiates an 18- or 21-state slice against more
independent variables than the slice has entries.  A residual that read its AD
derivative width from the slice length instead of the lift width raised
``AD derivative dimensions differ`` for the H18 magnetometer and both A21
measurement events, so no word containing them could ever be composed.
"""
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / 'tools/stability')]

from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_source_cover_contract as COVER
from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import physical_lineage as LINEAGE

I = Interval.point
MEASUREMENT_KINDS = ('accelerometer', 'magnetometer', 'S_zero')


def eye(n):
    return [[I(1 if i == j else 0) for j in range(n)] for i in range(n)]


def primitive():
    zero = (I(0), I(0), I(0))
    return KERNEL.WavePrimitivePayload('gen', 'p0', 'p1', 'live',
                                       zero, zero, zero, zero, zero, zero)


def cell(mode, kind):
    n = 18 if mode == 'H' else 21
    fields = dict(
        source_token='c1:e0', predecessor_token='c1', mode=mode, sample_index=0,
        event_ordinal=0, kind=kind, state=[I(0)] * n, P=eye(n), dt_s=I(0.005),
        tau_applied_s=I(1.1), sigma_aw_mps2=I(0.5), pseudo_elapsed_s=I(0.1),
        radial_scale=I(1), estimator_source_token='c1',
        estimator_predecessor_token='root', estimator_generated_coefficients=True,
        wave_primitive=primitive(), R=eye(3))
    if kind == 'accelerometer':
        fields.update(f_hat=[I(0), I(0), I(-9.80665)], R_hat=eye(3))
    if kind == 'magnetometer':
        fields.update(m_body=[I(20), I(0), I(40)])
    if kind == 'S_zero':
        fields.update(R_provenance=EVENTS.ACTUAL_RS_PROVENANCE)
    if mode == 'A':
        fields.update(true_bias=[I(0)] * 3, bias_projection_limit=0.4)
    return COVER.SourceCoverCell(**fields)


def lineage(mode, kind):
    selector = SimpleNamespace(
        prefix_length=1, parent_source_cell_id='root', source_cell_id='c1',
        sample_coordinates=SimpleNamespace(omega_body_corrected=(I(0),) * 3))
    return SimpleNamespace(mode=mode, cells=(cell(mode, kind),), selector=selector)


def compose(mode, kind):
    contract = next(c for c in BIAS.contracts() if c.name == 'BIAS0')
    return LINEAGE.compose_attached_sample(
        lineage(mode, kind), bias_contract=contract,
        held_bias_error=[Interval(-0.75, 0.75)] * 3 if mode == 'H' else None,
        true_bias=[Interval(-0.35, 0.35)] * 3 if mode == 'H' else None,
        tau_ba=I(300.0) if mode == 'A' else None)


class LineageEventPathTests(unittest.TestCase):
    def test_every_mode_and_measurement_kind_composes(self):
        for mode in ('H', 'A'):
            for kind in MEASUREMENT_KINDS:
                with self.subTest(mode=mode, kind=kind):
                    composed = compose(mode, kind)
                    J = composed['J_joint24']
                    self.assertEqual((len(J), len(J[0])), (24, 24))

    def test_magnetometer_event_is_sensitive_to_attitude(self):
        for mode in ('H', 'A'):
            with self.subTest(mode=mode):
                J = compose(mode, 'magnetometer')['J_joint24']
                touched = any(not (J[i][j].lo == 0 and J[i][j].hi == 0)
                              for i in range(3) for j in range(3) if i != j)
                self.assertTrue(touched)


class ResidualDerivativeWidthTests(unittest.TestCase):
    def test_residuals_accept_a_slice_narrower_than_its_lift(self):
        lifted = EVENTS.AD.independent_vector([I(0)] * 24, n=24, offset=0)
        z = lifted[:18]
        self.assertEqual(len(z), 18)
        self.assertEqual(EVENTS._derivative_dimension(z), 24)
        for residual in (EVENTS.residual_magnetometer(z, [I(20), I(0), I(40)]),
                         EVENTS.residual_accelerometer(z, [I(0), I(0), I(-9.80665)], eye(3))):
            self.assertEqual(len(residual), 3)
            self.assertEqual(residual[0].n, 24)

    def test_unlifted_state_keeps_its_own_width(self):
        z = EVENTS._state_ad([I(0)] * 21)
        self.assertEqual(EVENTS._derivative_dimension(z), 21)
        # The A21 accelerometer residual adds the estimator bias row, and that
        # choice follows the state layout rather than the derivative width.
        y = EVENTS.residual_accelerometer(z, [I(0), I(0), I(-9.80665)], eye(3))
        self.assertEqual(y[0].n, 21)

    def test_empty_slice_fails_closed(self):
        with self.assertRaises(ValueError):
            EVENTS._derivative_dimension([])


if __name__ == '__main__':
    unittest.main()
