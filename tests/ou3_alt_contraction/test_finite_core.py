"""Exact algebra regressions, not BRMM admission or stability diagnostics."""
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.stability.ou3_alt_contraction import finite_core as C
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


def root(mode='A', family='BIAS1'):
    ref = C.Reference(time=F(7, 200), live_origin=0,
                      q_world_to_body=(1, 0, 0, 0), acceleration=(F(1, 10), 0, 0),
                      velocity=(F(1, 20), 0, 0), position=(F(1, 30), 0, 0),
                      centered_S=(F(1, 1000), 0, 0), beta=(F(1, 100), 0, 0),
                      gyro_bias=(F(1, 10000), 0, 0), history_id='unqualified-physical-history',
                      bias_root='one-bias-root', bias_family=family)
    z = [F(0)]*24
    z[18], z[21:24] = F(1, 200), list(ref.beta)
    u = [F((i % 3)-1, 100) if mode == 'A' or i < 18 else F(0) for i in range(21)]
    cov = M.plus(M.scaled(M.eye(21), F(1, 100)), [[x*y for y in u] for x in u])
    if mode == 'H':
        for i in range(18, 21):
            cov[i][i] = F(1, 62500)
    return C.State(mode, tuple(z), cov, (1, 0, 0, 0), ref)


def physical_successor(state):
    b = state.reference
    h = F(1, 200)
    J0, J1, J2 = (F(1, 1800), 0, 0), (F(1, 750000), 0, 0), (F(1, 450000000), 0, 0)
    pt, wb = F(999, 1000), (F(1, 100000), 0, 0)
    after = replace(b, time=b.time+h, q_world_to_body=(1, F(1, 3000), F(-1, 4000), 0),
                    acceleration=(F(1, 9), 0, 0),
                    velocity=tuple(b.velocity[i]+J0[i] for i in range(3)),
                    position=tuple(b.position[i]+h*b.velocity[i]+J1[i] for i in range(3)),
                    centered_S=tuple(b.centered_S[i]+h*b.position[i]+h*h*b.velocity[i]/2+J2[i] for i in range(3)),
                    beta=tuple(pt*b.beta[i]+wb[i] for i in range(3)),
                    gyro_bias=tuple(b.gyro_bias[i]+F(1, 1000000) for i in range(3)))
    segment = C.PHYSICAL.PhysicalSegment(b, after, J0, J1, J2, pt, wb)
    # Coefficients exercise an identity, NOT source membership or tuning values.
    coeff = (F(1, 201), F(1, 81000), F(1, 48100000), F(99, 100), h)
    return segment, dict(gyro_body=[0, 0, 0], axis_coefficients=[coeff]*3,
                         phi_hat=F(999, 1000), F21=M.eye(21), Q21=M.zeros(21, 21))


class FiniteCoreTests(unittest.TestCase):
    def test_physical_prediction_keeps_mode_and_unqualified_bias_labels(self):
        for mode in ('H', 'A'):
            for family in ('BIAS0', 'BIAS1', 'BIAS2'):
                with self.subTest(mode=mode, family=family):
                    before = root(mode, family)
                    segment, kw = physical_successor(before)
                    after = segment.after
                    out = C.prediction(before, segment, **kw)
                    self.assertEqual(out.reference, after)
                    self.assertEqual(out.z[:3], C.cayley(P.quat_mul(after.q_world_to_body, P.quat_conj(out.q_hat))))
                    self.assertNotEqual(out.z[:3], before.z[:3])
                    self.assertEqual(out.z[21:24], after.beta)
                    for i in range(3):
                        self.assertEqual(after.gyro_bias[i]-out.z[3+i], before.reference.gyro_bias[i]-before.z[3+i])
                        factor = 1 if mode == 'H' else kw['phi_hat']
                        self.assertEqual(after.beta[i]-out.z[18+i], factor*(before.reference.beta[i]-before.z[18+i]))

    def test_prediction_matches_direct_truth_minus_nominal_translation(self):
        s = root()
        segment, kw = physical_successor(s)
        after = segment.after
        out = C.prediction(s, segment, **kw)
        pv, pp, ps, alpha, h = kw['axis_coefficients'][0]
        old = s.reference
        for i in range(3):
            vh, ph, Sh, ah = old.velocity[i]-s.z[6+i], old.position[i]-s.z[9+i], old.centered_S[i]-s.z[12+i], old.acceleration[i]-s.z[15+i]
            hats = (vh+pv*ah, ph+h*vh+pp*ah, Sh+h*ph+h*h*vh/2+ps*ah, alpha*ah)
            for j, field in enumerate(('velocity', 'position', 'centered_S', 'acceleration')):
                self.assertEqual(getattr(after, field)[i]-out.z[6+3*j+i], hats[j])

    def test_origin_root_primitive_and_shared_bias_mutations_rejected(self):
        s = root()
        segment, kw = physical_successor(s)
        for field, value in (
            ('history_id', 'other'), ('live_origin', F(-1)), ('bias_root', 'fresh'),
            ('bias_family', 'BIAS2'), ('time', segment.after.time+F(1, 200)),
            ('velocity', (0, 0, 0)), ('position', (0, 0, 0)),
            ('centered_S', (0, 0, 0)), ('beta', (0, 0, 0))):
            with self.assertRaises(ValueError):
                changed = replace(segment, after=replace(segment.after, **{field: value}))
                C.prediction(s, changed, **kw)
        with self.assertRaises(ValueError):
            changed = replace(segment, before=replace(segment.before, gyro_bias=(0, 0, 0)))
            C.prediction(s, changed, **kw)
        with self.assertRaises(TypeError):
            C.prediction(s, segment, q15=[0]*15, **kw)  # No detached q15 argument.
        with self.assertRaises(ValueError):
            replace(s.reference, time=s.reference.live_origin)  # Fresh S must be zero.

    def test_state_requires_21_covariance_and_same_physical_attitude_and_beta(self):
        s = root()
        with self.assertRaises(ValueError):
            replace(s, covariance=M.eye(18))
        z = list(s.z); z[21] += 1
        with self.assertRaises(ValueError):
            replace(s, z=tuple(z))
        with self.assertRaises(ValueError):
            replace(s, q_hat=(1, F(1, 10), 0, 0))
        for q in ((0, 0, 0, 0), (0, 1, 0, 0)):
            with self.assertRaises(ValueError):
                C.cayley(q)
        with self.assertRaises(TypeError):
            replace(s, q_hat=(1.0, 0, 0, 0))
        nonsymmetric = M.eye(21); nonsymmetric[0][3] = F(1, 10)
        with self.assertRaises(ValueError):
            C.covariance_reset(nonsymmetric, [0, 0, 0])

    def test_measurements_match_direct_nominal_update_and_full_covariance_reset(self):
        for mode in ('H', 'A'):
            s = root(mode)
            # Nonzero finite physical attitude, not a tangent perturbation.
            ref = replace(s.reference, q_world_to_body=(1, F(1, 10000), F(-1, 20000), 0))
            z = list(s.z); z[:3] = C.cayley(ref.q_world_to_body)
            s = replace(s, z=tuple(z), reference=ref)
            for kind, packet in (
                ('S_zero', {}),
                ('accelerometer', dict(observed=[F(11, 100), 0, F(-196133, 20000)])),
                ('magnetometer', dict(observed=[22, F(1, 10000), 43], magnetic_reference=[22, 0, 43])),
            ):
                with self.subTest(mode=mode, kind=kind):
                    result = C.measurement(s, kind, R=M.eye(3), **packet)
                    out, d, K, S, N = result.state, result.correction, result.gain, result.innovation, result.numerator
                    self.assertEqual(out.reference, s.reference)
                    self.assertEqual(out.z[3:21], tuple(s.z[i]-d[i] for i in range(3, 21)))
                    self.assertEqual(M.mm(K, S), list(map(list, N)))
                    Pj = M.plus(M.plus(M.plus(s.covariance, M.mm(K, M.transpose(N)), -1),
                                        M.mm(N, M.transpose(K)), -1), M.mm(M.mm(K, S), M.transpose(K)))
                    reset = M.eye(21)
                    for i, row in enumerate(M.plus(M.eye(3), M.scaled(M.skew(d[:3]), F(1, 2)))):
                        reset[i][:3] = row
                    dense = M.mm(M.mm(reset, Pj), M.transpose(reset))
                    self.assertEqual(out.covariance, tuple(map(tuple, dense)))
                    if mode == 'H':
                        self.assertEqual(d[18:21], (0, 0, 0))
                        self.assertEqual(tuple(row[18:] for row in out.covariance), tuple(row[18:] for row in s.covariance))

    def test_retry_shift_changes_both_gain_and_Joseph_not_only_final_gain(self):
        s = root('H')
        a = C.measurement(s, 'S_zero', R=M.eye(3))
        b = C.measurement(s, 'S_zero', R=M.eye(3), innovation_shift=F(1, 10000))
        self.assertEqual(M.plus(b.innovation, a.innovation, -1), M.scaled(M.eye(3), F(1, 10000)))
        self.assertNotEqual(a.gain, b.gain)
        self.assertNotEqual(a.state.covariance, b.state.covariance)
        with self.assertRaises(ValueError):
            C.measurement(s, 'S_zero', R=M.eye(3), innovation_shift=-1)
        with self.assertRaises(ValueError):
            C.solve3(M.zeros(3, 3), [0, 0, 0])
        # Nonsingular indefinite is NOT incorrectly advertised as an SPD gate.
        self.assertEqual(C.solve3([[1, 0, 0], [0, -1, 0], [0, 0, 1]], [1, 1, 1]), [1, -1, 1])

    def test_active_projection_keeps_beta_and_does_not_project_covariance(self):
        s = root()
        ref = replace(s.reference, centered_S=(0, 0, 0))
        z = list(s.z); z[12] = F(12, 5); z[18] = ref.beta[0]-F(1, 5)
        cov = M.eye(21); cov[12][18] = cov[18][12] = F(1, 2)
        s = replace(s, z=tuple(z), covariance=cov, reference=ref)
        with self.assertRaises(ValueError):
            C.measurement(s, 'S_zero', R=M.eye(3), alpha=1)
        out = C.measurement(s, 'S_zero', R=M.eye(3), alpha=F(1, 2)).state
        self.assertEqual(ref.beta[0]-out.z[18], F(2, 5))
        self.assertEqual(out.z[21:24], s.z[21:24])
        self.assertEqual(out.covariance[18][18], F(7, 8))

    def test_finite_accepted_composition_retains_covariance_and_every_prefix(self):
        s = root('H')
        packets = [dict(kind='S_zero', R=M.eye(3)),
                   dict(kind='accelerometer', R=M.eye(3), observed=[F(11, 100), 0, F(-196133, 20000)]),
                   dict(kind='magnetometer', R=M.eye(3), observed=[22, 0, 43], magnetic_reference=[22, 0, 43])]
        whole = C.compose_accepted(s, packets)
        prefix = C.compose_accepted(s, packets[:1])
        suffix = C.compose_accepted(prefix[-1], packets[1:])
        self.assertEqual(whole, prefix+suffix[1:])
        self.assertEqual(C.compose_accepted(s, []), (s,))
        self.assertTrue(all(x.reference is s.reference for x in whole))
        actual = C.measurement(whole[1], **packets[1])
        detached = C.measurement(replace(whole[1], covariance=s.covariance), **packets[1])
        self.assertNotEqual(actual.gain, detached.gain)
        self.assertNotEqual(actual.state.covariance, detached.state.covariance)

    def test_S_event_at_zero_error_is_not_homogeneous(self):
        s = root('H')
        z = [F(0)]*24; z[21:24] = s.reference.beta
        s = replace(s, z=tuple(z))
        result = C.measurement(s, 'S_zero', R=M.eye(3))
        self.assertNotEqual(result.state.z[12:15], s.z[12:15])
        with self.assertRaises(ValueError):
            C.measurement(s, 'S_zero', R=M.eye(3), observed=[0, 0, 0])

    def test_unimplemented_branches_are_not_silently_replaced_by_identity(self):
        s = root()
        with self.assertRaises(ValueError):
            C.measurement(s, 'rejected_after_nonfinite_correction', R=M.eye(3))
        with self.assertRaises(ValueError):
            C.measurement(s, 'magnetometer', R=M.eye(3), observed=[10000, 10000, 10000], magnetic_reference=[22, 0, 43])
        # There is no complete-source/storage promotion method on this core.
        self.assertFalse(hasattr(C, 'ALT_LIVE_PASS'))
        self.assertFalse(hasattr(C, 'storage_search_allowed'))


if __name__ == '__main__':
    unittest.main()
