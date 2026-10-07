"""Exact rational/dual-number regression; no secants or sampled orbit."""
from fractions import Fraction as F
import unittest
from tools.stability.ou3_theorem.matrix_certificates import add, transpose
from tools.stability.ou3_theorem.planar_linked_riccati_mean import (
    certificate, product, gain_differential, gain_finite_remainder,
    innovation_differential, optimal_covariance_differential, joseph_differential,
    reset_differential,
)


class Dual:
    def __init__(self, v, d=0):
        self.v, self.d = F(v), F(d)

    @staticmethod
    def lift(x):
        return x if isinstance(x, Dual) else Dual(x)

    def __add__(self, other):
        b = self.lift(other)
        return Dual(self.v+b.v, self.d+b.d)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, -self.d)

    def __sub__(self, other):
        return self + (-self.lift(other))

    def __rsub__(self, other):
        return self.lift(other) + (-self)

    def __mul__(self, other):
        b = self.lift(other)
        return Dual(self.v*b.v, self.d*b.v+self.v*b.d)

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = self.lift(other)
        return Dual(self.v/b.v, (self.d*b.v-self.v*b.d)/b.v**2)

    def __rtruediv__(self, other):
        return self.lift(other)/self


def inverse2(a):
    det = a[0][0]*a[1][1]-a[0][1]*a[1][0]
    return [[a[1][1]/det, -a[0][1]/det], [-a[1][0]/det, a[0][0]/det]]


def jet(a, da):
    return [[Dual(x, dx) for x, dx in zip(r, dr)] for r, dr in zip(a, da)]


def deriv(a):
    return [[x.d for x in r] for r in a]


def joseph(p, b, s):
    k = product(b, inverse2(s))
    return add(add(add(p, product(k, transpose(b)), -1),
                   product(b, transpose(k)), -1), product(k, s, transpose(k)))


class LinkedRiccatiMeanTest(unittest.TestCase):
    def setUp(self):
        self.p = [[F(2), F(1, 3), F(1, 5)], [F(1, 3), F(3), F(1, 7)], [F(1, 5), F(1, 7), F(4)]]
        self.h = [[F(1), F(2), F(1)], [F(-1), F(1), F(2)]]
        self.r = [[F(2), F(1, 7)], [F(1, 7), F(3)]]
        self.dp = [[F(1, 10), F(1, 30), F(-1, 20)], [F(1, 30), F(-1, 20), F(1, 40)], [F(-1, 20), F(1, 40), F(1, 10)]]
        self.dh = [[F(1, 30), F(-1, 10), F(1, 40)], [F(1, 20), F(1, 60), F(-1, 30)]]
        self.dr = [[F(1, 20), F(1, 50)], [F(1, 50), F(-1, 30)]]

    def test_linked_derivatives_with_noise_against_exact_dual_algebra(self):
        p, h, r, dp, dh, dr = self.p, self.h, self.r, self.dp, self.dh, self.dr
        s = add(product(h, p, transpose(h)), r)
        si = inverse2(s)
        k = product(p, transpose(h), si)
        pj, hj, rj = jet(p, dp), jet(h, dh), jet(r, dr)
        sj = add(product(hj, pj, transpose(hj)), rj)
        bj = product(pj, transpose(hj))
        self.assertEqual(gain_differential(p, h, k, si, dp, dh, dr), deriv(product(bj, inverse2(sj))))
        self.assertEqual(optimal_covariance_differential(p, h, k, dp, dh, dr), deriv(joseph(pj, bj, sj)))
        p1, h1, r1 = add(p, dp), add(h, dh), add(r, dr)
        si1 = inverse2(add(product(h1, p1, transpose(h1)), r1))
        delta_k = add(product(p1, transpose(h1), si1), k, -1)
        linear = gain_differential(p, h, k, si, dp, dh, dr)
        self.assertEqual(add(delta_k, linear, -1), gain_finite_remainder(p, h, k, si, si1, dp, dh, dr))

    def test_held_ba_distinct_innovation_and_gain_rows(self):
        p, h, r, dp, dh, dr = self.p, self.h, self.r, self.dp, self.dh, self.dr
        hg, dhg = [row[:2]+[F(0)] for row in h], [row[:2]+[F(0)] for row in dh]
        mask = [[F(1), F(0), F(0)], [F(0), F(1), F(0)], [F(0), F(0), F(0)]]
        b = product(mask, p, transpose(hg))
        db = product(mask, add(product(dp, transpose(hg)), product(p, transpose(dhg))))
        s = add(product(h, p, transpose(h)), r)
        ds = innovation_differential(p, h, dp, dh, dr)
        k = product(b, inverse2(s))
        dk = product(add(db, product(k, ds), -1), inverse2(s))
        pj, hj, rj = jet(p, dp), jet(h, dh), jet(r, dr)
        bj = product(mask, pj, transpose(jet(hg, dhg)))
        sj = add(product(hj, pj, transpose(hj)), rj)
        self.assertEqual(joseph_differential(k, b, s, dp, dk, db, ds), deriv(joseph(pj, bj, sj)))

    def test_reset_and_nonpromotion(self):
        g, dg = self.p, self.dp
        c, dc = self.p, self.dp
        self.assertEqual(reset_differential(g, c, dg, dc),
                         deriv(product(jet(g, dg), jet(c, dc), transpose(jet(g, dg)))))
        cert = certificate()
        self.assertFalse(cert['finite_word_uniform_bound_closed'])
        self.assertIn('dK*r+K*dr', cert['mean_differential'])


if __name__ == '__main__':
    unittest.main()
