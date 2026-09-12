"""Universal coefficient checks for finite physical prediction; not source admission."""
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(Path(__file__).parent)]
from test_finite_measurement_graph import Poly
from tools.stability.ou3_alt_contraction import finite_physical_prediction as G
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M


def symbolic(n):
    return [Poly.var(i) for i in range(n)]


def endpoint(**kwargs):
    values = dict(time=F(0), live_origin=F(0), q_world_to_body=(1, 0, 0, 0),
                  velocity=(0, 0, 0), position=(0, 0, 0), centered_S=(0, 0, 0),
                  acceleration=(0, 0, 0), gyro_bias=(0, 0, 0), beta=(F(1, 20), 0, 0))
    values.update(kwargs)
    return G.PhysicalKinematics(**values)


class PhysicalPredictionTests(unittest.TestCase):
    def test_finite_attitude_identity_all_fifteen_indeterminates(self):
        # q_true,0 = (2,c) q_hat,0.  This proof retains an arbitrary physical
        # increment rather than replacing it by the estimator's OU/gyro step.
        with patch.object(Poly, 'N', 15):
            variables = symbolic(15)
            c, qh, qt_step, qh_step = variables[:3], variables[3:7], variables[7:11], variables[11:15]
            qt0 = G.quaternion_product([2, *c], qh)
            qt1 = G.quaternion_product(qt_step, qt0)
            qh1 = G.quaternion_product(qh_step, qh)
            relative = G.quaternion_product(qt1, G.conjugate(qh1))
            W, numerator = G.attitude_polynomials(c, qt_step, qh_step)
            scale = M.dot(qh, qh)
            self.assertEqual(relative[0], scale*W)
            self.assertEqual([2*v for v in relative[1:]], [scale*v for v in numerator])

    def test_physical_step_defect_order_all_twelve_indeterminates(self):
        with patch.object(Poly, 'N', 12):
            x = symbolic(12)
            D, sample, physical = x[:4], x[4:8], x[8:12]
            # If D=Q_phys conjugate(Q_sample), D Q_sample is exactly a
            # nonzero projective multiple of Q_phys, not Q_sample Q_phys.
            defect = G.quaternion_product(physical, G.conjugate(sample))
            self.assertEqual(G.physical_step_from_defect(defect, sample),
                             [M.dot(sample, sample)*v for v in physical])
            self.assertEqual(G.physical_step_from_defect(D, sample), G.quaternion_product(D, sample))

    def test_translation_finite_identity_all_seventeen_indeterminates(self):
        with patch.object(Poly, 'N', 17):
            x = symbolic(17)
            v, p, S, a, ev, ep, eS, ea, h, va, pa, Sa, alpha, a1, J0, J1, J2 = x
            truth_next = [v+J0, p+h*v+J1, S+h*p+F(1, 2)*h*h*v+J2, a1]
            estimate_next = [(v-ev)+va*(a-ea),
                             (p-ep)+h*(v-ev)+pa*(a-ea),
                             (S-eS)+h*(p-ep)+F(1, 2)*h*h*(v-ev)+Sa*(a-ea),
                             alpha*(a-ea)]
            actual = [x-y for x, y in zip(truth_next, estimate_next)]
            expected = G.linear_prediction_polynomials([ev, ep, eS, ea], a, a1, J0, J1, J2,
                h=h, va=va, pa=pa, Sa=Sa, alpha=alpha)
            self.assertEqual(actual, expected)

    def test_nonzero_physical_attitude_anchor_is_retained(self):
        physical = (1, F(1, 400), F(-1, 700), F(1, 900))
        actual = G.attitude((0, 0, 0), physical, (1, 0, 0, 0))
        self.assertEqual(actual, [2*x for x in physical[1:]])
        self.assertNotEqual(actual, [0, 0, 0])
        with self.assertRaisesRegex(ValueError, 'pole'):
            G.attitude((0, 0, 0), (0, 1, 0, 0), (1, 0, 0, 0))

    def test_same_driver_and_gyro_bias_increment_for_held_and_active(self):
        before = endpoint()
        # BIAS families are not inferred from these algebraic cases. The
        # physical recurrence is checked separately for each factor/root.
        for physical_factor, driver in ((F(1), (0, 0, 0)),
                                        (F(999, 1000), (F(1, 10**6), 0, 0)),
                                        (F(19, 20), (F(-1, 10000), 0, 0))):
            beta_next = tuple(physical_factor*b+u for b, u in zip(before.beta, driver))
            after = endpoint(time=F(1, 200), beta=beta_next,
                             gyro_bias=(F(1, 100000), 0, 0))
            segment = G.PhysicalSegment(before, after, (0, 0, 0), (0, 0, 0), (0, 0, 0), physical_factor, driver)
            z = [F(0)]*24; z[18:21] = [F(1, 25), 0, 0]; z[21:24] = before.beta
            coeff = [(F(1, 200), F(1, 80000), F(1, 48000000), F(999, 1000), F(1, 200))]*3
            for active in (False, True):
                ph = F(999, 1000) if active else F(1)
                out = G.prediction(z, segment, nominal_step=(1, 0, 0, 0),
                                   axis_coefficients=coeff, active_bias=active, phi_hat=ph)
                self.assertEqual(out[3:6], list(after.gyro_bias))
                self.assertEqual(out[21:24], list(after.beta))
                self.assertEqual(out[18:21], [ph*z[18+i]+after.beta[i]-ph*before.beta[i] for i in range(3)])

    def test_reject_detached_moments_bias_time_and_origin(self):
        before = endpoint(position=(F(1, 10), 0, 0))
        after = endpoint(time=F(1, 200), position=before.position,
                         centered_S=(F(1, 2000), 0, 0))
        segment = G.PhysicalSegment(before, after, (0, 0, 0), (0, 0, 0), (0, 0, 0), F(1), (0, 0, 0))
        with self.assertRaisesRegex(ValueError, 'one history'):
            replace(segment, J1=(F(1), 0, 0))
        with self.assertRaisesRegex(ValueError, 'one history'):
            replace(segment, bias_driver=(F(1), 0, 0))
        with self.assertRaisesRegex(ValueError, 'one persistent Live origin'):
            replace(segment, after=replace(after, live_origin=F(-1)))
        with self.assertRaisesRegex(ValueError, 'fresh centered physical S'):
            endpoint(centered_S=(1, 0, 0))
        z = [F(0)]*24
        with self.assertRaisesRegex(ValueError, 'physical predecessor'):
            G.prediction(z, segment, nominal_step=(1, 0, 0, 0),
                         axis_coefficients=[(0, 0, 0, 1, F(1, 200))]*3, active_bias=False)

    def test_finite_word_prediction_composes_without_S_or_position_restart(self):
        # Algebraic constant-acceleration segment, NOT an admitted all-time
        # physical source. Only checks concatenation of the finite identity.
        h = F(1, 200); b = endpoint(acceleration=(F(1, 10), 0, 0))
        z = [F(0)]*24; z[15:18] = b.acceleration; z[18:21] = b.beta; z[21:24] = b.beta
        initial_origin = b.live_origin
        for _ in range(600):
            J0 = tuple(h*a for a in b.acceleration)
            J1 = tuple(h*h*a/2 for a in b.acceleration)
            J2 = tuple(h*h*h*a/6 for a in b.acceleration)
            a = replace(b, time=b.time+h,
                        velocity=tuple(b.velocity[i]+J0[i] for i in range(3)),
                        position=tuple(b.position[i]+h*b.velocity[i]+J1[i] for i in range(3)),
                        centered_S=tuple(b.centered_S[i]+h*b.position[i]+h*h*b.velocity[i]/2+J2[i] for i in range(3)))
            segment = G.PhysicalSegment(b, a, J0, J1, J2, F(1), (0, 0, 0))
            z = G.prediction(z, segment, nominal_step=(1, 0, 0, 0),
                            axis_coefficients=[(h, h*h/2, h*h*h/6, F(1), h)]*3,
                            active_bias=False)
            self.assertEqual(z[6:18], list(a.velocity+a.position+a.centered_S+a.acceleration))
            self.assertEqual(a.live_origin, initial_origin)
            b = a
        self.assertEqual(z[12], F(9, 20))


if __name__ == '__main__':
    unittest.main()
