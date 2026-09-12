"""All-variable moment identities and negative source-binding regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest
from unittest.mock import patch

from test_finite_measurement_graph import Poly
from tools.stability.ou3_alt_contraction import finite_brmm_moment_prefix as X
from tools.stability.ou3_alt_contraction import finite_source_continuation as S
from tools.stability.ou3_alt_contraction import finite_physical_prediction as P
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as RAW
from tools.stability.ou3_alt_contraction import bias_families as BIAS
import test_finite_source_continuation as BASE


def mm(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))


def plus(a, b):
    return tuple(tuple(x+y for x, y in zip(ar, br)) for ar, br in zip(a, b))


def tr(a):
    return tuple(zip(*a))


def segment(before=None, *, h=F(1, 200), J=None, phi=BASE.VALID_PHI):
    before = BASE.kin(0) if before is None else before
    J = (X.ZERO, X.ZERO, X.ZERO) if J is None else J
    after = replace(before, time=before.time+h,
                    velocity=tuple(before.velocity[i]+J[0][i] for i in range(3)),
                    position=tuple(before.position[i]+h*before.velocity[i]+J[1][i] for i in range(3)),
                    centered_S=tuple(before.centered_S[i]+h*before.position[i]+h*h*before.velocity[i]/2+J[2][i]
                                     for i in range(3)),
                    beta=tuple(phi*x for x in before.beta))
    return P.PhysicalSegment(before, after, *J, phi, X.ZERO)


class Tests(unittest.TestCase):
    def test_gram_concatenation_identity_all_duration_coefficients(self):
        with patch.object(Poly, 'N', 2):
            a, b = Poly.var(0), Poly.var(1)
            lhs = X.gram(a+b)
            rhs = plus(mm(mm(X.shift(b), X.gram(a)), tr(X.shift(b))), X.gram(b))
            self.assertEqual(lhs, rhs)
            self.assertEqual(mm(X.shift(a), X.shift(b)), X.shift(a+b))

    def test_concatenation_associativity_all_29_indeterminates(self):
        with patch.object(Poly, 'N', 29):
            v = [Poly.var(i) for i in range(29)]
            Ja, Jb, Jc = [tuple(tuple(v[start+3*j+i] for i in range(3)) for j in range(3))
                          for start in (0, 9, 18)]
            b, c = v[27:]
            self.assertEqual(X.concatenate(X.concatenate(Ja, Jb, b), Jc, c),
                             X.concatenate(Ja, X.concatenate(Jb, Jc, c), b+c))

    def test_exact_declared_caps_are_physical_not_filter_error_bounds(self):
        b = X.BOUNDS
        self.assertEqual((b.acceleration, b.velocity, b.position, b.centered_S),
                         (F(44, 5), F(11, 2), F(81, 10), F(1100)))
        self.assertEqual(b.angular_rate_upper, F(11, 18))
        for field, cap in (('acceleration', b.acceleration), ('velocity', b.velocity),
                           ('position', b.position), ('centered_S', b.centered_S)):
            # A later endpoint permits a nonzero centered S.
            q = replace(BASE.kin(1), **{field: (cap, 0, 0)})
            X.check_endpoint(q)
            with self.assertRaisesRegex(ValueError, 'vector cap'):
                X.check_endpoint(replace(q, **{field: (cap, F(1, 10**9), 0)}))

    def test_detached_moment_sign_corner_is_rejected_despite_scalar_caps(self):
        h = F(1, 200); A = X.BOUNDS.acceleration
        J = ((A*h, 0, 0), (-A*h*h/2, 0, 0), (A*h*h*h/6, 0, 0))
        q = segment(J=J)
        X.check_endpoint(q.before); X.check_endpoint(q.after)
        self.assertEqual(X.energy(h, J), 193*h*A*A)
        with self.assertRaisesRegex(ValueError, 'coupled three-axis'):
            S.QualifiedPhysicalSegment(BASE.root(), BASE.wit(1, 'root', 'c1', 'p0', 'p1'), q)

    def test_three_axes_share_one_acceleration_budget(self):
        h = F(1, 200); A = X.BOUNDS.acceleration
        J = ((A*h, A*h, 0), (A*h*h/2, A*h*h/2, 0), (A*h*h*h/6, A*h*h*h/6, 0))
        q = segment(J=J)
        self.assertEqual(X.energy(h, J), 2*h*A*A)
        with self.assertRaisesRegex(ValueError, 'coupled three-axis'):
            X.check_segment(q)

    def test_prefix_has_exact_moments_and_derived_energy_without_restarting_S(self):
        h = F(1, 200); acceleration = (F(3), F(-2), F(1))
        J = tuple(tuple(a*h**(j+1)/F((1, 2, 6)[j]) for a in acceleration) for j in range(3))
        q1 = segment(before=replace(BASE.kin(0), acceleration=acceleration), J=J)
        q2 = segment(q1.after, J=J)
        prefix = X.append(X.append(X.Prefix(), q1), q2)
        H = 2*h
        self.assertEqual(prefix.moments,
                         tuple(tuple(a*H**(j+1)/F((1, 2, 6)[j]) for a in acceleration) for j in range(3)))
        self.assertEqual(prefix.minimum_energy, H*14)
        self.assertEqual(prefix.segment_energy_sum, H*14)
        self.assertEqual(q2.after.centered_S, prefix.J2)
        with self.assertRaisesRegex(ValueError, 'energy budget'):
            replace(prefix, segment_energy_sum=H*14-F(1, 10**9))

    def test_piecewise_change_has_projection_loss_not_independent_prefix_supply(self):
        h = F(1, 200)
        J = ((h, 0, 0), (h*h/2, 0, 0), (h*h*h/6, 0, 0))
        q1 = segment(J=J); q2 = segment(q1.after, J=tuple(tuple(-x for x in row) for row in J))
        prefix = X.append(X.append(X.Prefix(), q1), q2)
        self.assertEqual(prefix.segment_energy_sum, 2*h)
        self.assertLess(prefix.minimum_energy, prefix.segment_energy_sum)
        self.assertEqual(prefix.J0, X.ZERO)
        self.assertNotEqual(prefix.J1, X.ZERO)

    def test_large_algebraically_valid_acceleration_is_not_source_qualified(self):
        q = segment(before=replace(BASE.kin(0), acceleration=(100, 0, 0)))
        with self.assertRaisesRegex(ValueError, 'acceleration vector cap'):
            S.QualifiedPhysicalSegment(BASE.root(), BASE.wit(1, 'root', 'c1', 'p0', 'p1'), q)

    def test_rotation_chord_is_scale_and_sign_invariant_and_rejects_jump(self):
        q = segment()
        X.check_segment(replace(q, after=replace(q.after, q_world_to_body=(-7, 0, 0, 0))))
        with self.assertRaisesRegex(ValueError, 'rotation chord'):
            X.check_segment(replace(q, after=replace(q.after, q_world_to_body=(0, 1, 0, 0))))

    def test_sampled_physical_rate_not_bias_corrected_estimator_rate_is_bounded(self):
        q = BASE.qseg(BASE.root()); raw = BASE.raw_for(q)
        changed = replace(raw, omega_sample_internal=(1, 0, 0), raw_gyro_body=(1, 0, 0))
        with self.assertRaisesRegex(ValueError, 'physical sampled angular-rate'):
            S.qualify_raw_imu(q, BASE.sensor_root(q.root), changed, 'packet')
        # A sensor disturbance is NOT an actual physical angular rate. It is
        # retained as forcing and not silently clipped to the physical cap.
        noisy = replace(raw, gyro_residual_internal=(1, 0, 0), raw_gyro_body=(1, 0, 0))
        S.qualify_raw_imu(q, BASE.sensor_root(q.root), noisy, 'noisy-packet')

    def test_all_three_BIAS_families_keep_one_actual_factor_not_just_token(self):
        for c in BIAS.contracts():
            with self.subTest(family=c.name):
                root = S.SourceRoot('hist', 'generator', 0, c.name, c.parameter_token)
                lo, hi = F.from_float(c.phi_true.lo), F.from_float(c.phi_true.hi)
                phi = (lo+hi)/2
                q1 = segment(phi=phi); q2 = segment(q1.after, phi=phi)
                cont = S.append(S.begin(root), witness=BASE.wit(1, 'root', 'c1', 'p0', 'p1'), segment=q1)
                out = S.append(cont, witness=BASE.wit(2, 'c1', 'c2', 'p1', 'p2'), segment=q2)
                self.assertEqual(out.moment_prefix.duration, F(1, 100))
                with self.assertRaisesRegex(ValueError, 'BIAS factor changed'):
                    S.append(cont, witness=BASE.wit(2, 'c1', 'c2', 'p1', 'p2'), segment=segment(q1.after, phi=hi))
                # Preserve the original outward endpoints exactly, including
                # BIAS2's non-relaxing endpoint phi=1, without decimal rounding.
                for edge in (lo, hi):
                    S.QualifiedPhysicalSegment(root, BASE.wit(1, 'root', 'c1', 'p0', 'p1'), segment(phi=edge))
        self.assertEqual(next(c for c in BIAS.contracts() if c.name=='BIAS2').phi_true.hi, 1)

    def test_true_bias_component_bound_not_only_looser_norm_bound(self):
        c = BASE.BIAS0
        b = (F.from_float(c.true_bias_component_bound)+F.from_float(c.true_bias_norm_bound))/2
        q = segment(before=BASE.kin(0, beta=(b, 0, 0)))
        with self.assertRaisesRegex(ValueError, 'family component contract'):
            S.QualifiedPhysicalSegment(BASE.root(), BASE.wit(1, 'root', 'c1', 'p0', 'p1'), q)

    def test_continuation_prefix_summary_is_derived_not_caller_supplied(self):
        q = segment(); cont = S.append(S.begin(BASE.root()), witness=BASE.wit(1, 'root', 'c1', 'p0', 'p1'), segment=q)
        self.assertEqual(cont.moment_prefix.J0, q.J0)
        with self.assertRaisesRegex((TypeError, ValueError), 'init=False'):
            replace(cont, moment_prefix=X.Prefix())


if __name__ == '__main__':
    unittest.main()
