from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_p4_complete_brmm_signed_information_ledger as LEDGER


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def tr(A):
    return [list(x) for x in zip(*A)]


def add(A, B):
    return [[a + b for a, b in zip(ar, br)] for ar, br in zip(A, B)]


def sub(A, B):
    return [[a - b for a, b in zip(ar, br)] for ar, br in zip(A, B)]


def mul(A, B):
    return [[sum((x * y for x, y in zip(row, col)), F(0)) for col in tr(B)] for row in A]


def mv(A, x):
    return [sum((a * b for a, b in zip(row, x)), F(0)) for row in A]


def inv(A):
    n = len(A)
    aug = [list(row) + unit for row, unit in zip(A, eye(n))]
    for j in range(n):
        p = next(i for i in range(j, n) if aug[i][j] != 0)
        aug[j], aug[p] = aug[p], aug[j]
        pivot = aug[j][j]
        aug[j] = [x / pivot for x in aug[j]]
        for i in range(n):
            if i == j:
                continue
            q = aug[i][j]
            aug[i] = [x - q * y for x, y in zip(aug[i], aug[j])]
    return [row[n:] for row in aug]


def quad(A, x):
    return sum((a * b for a, b in zip(x, mv(A, x))), F(0))


def skew(v):
    x, y, z = v
    return [[F(0), -z, y], [z, F(0), -x], [-y, x, F(0)]]


def fixture(n):
    P = [[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        P[i][i] = F(2) + F(i + 1, 20)
    # Dense low-rank PSD addition keeps full cross-covariance in the identity.
    u = [F((i % 7) - 3, 100) for i in range(n)]
    for i in range(n):
        for j in range(n):
            P[i][j] += u[i] * u[j]

    H = [[F(0) for _ in range(n)] for _ in range(3)]
    H[0][0], H[0][3], H[0][15] = F(1), F(1, 3), F(2, 5)
    H[1][1], H[1][4], H[1][16] = F(1), F(-1, 4), F(1, 7)
    H[2][2], H[2][5], H[2][17] = F(1), F(1, 5), F(-1, 6)
    if n == 21:
        H[0][18], H[1][19], H[2][20] = F(1), F(1), F(1)

    R = [[F(0) for _ in range(3)] for _ in range(3)]
    R[0][0], R[1][1], R[2][2] = F(3, 2), F(5, 4), F(7, 5)
    Pinv = inv(P)
    S = add(mul(mul(H, P), tr(H)), R)
    Sinv = inv(S)
    K = mul(mul(P, tr(H)), Sinv)
    A = sub(eye(n), mul(K, H))
    Pj = add(mul(mul(A, P), tr(A)), mul(mul(K, R), tr(K)))
    Jj = inv(Pj)
    Rinv = inv(R)

    e = [F((i % 9) - 4, 70 + i) for i in range(n)]
    eta = [F(1, 11), F(-1, 13), F(1, 17)]
    return P, Pinv, H, R, Rinv, S, Sinv, K, Pj, Jj, e, eta


class SignedInformationLedgerTests(unittest.TestCase):
    def test_exact_joseph_identity_h18_a21(self):
        for n in (18, 21):
            with self.subTest(n=n):
                _, Pinv, H, _, Rinv, _, Sinv, K, _, Jj, e, eta = fixture(n)
                row = LEDGER.joseph_signed_energy_terms(
                    Pinv, Jj, H, Rinv, Sinv, K, e, eta
                )
                self.assertEqual(row["identity_residual"], F(0))
                self.assertEqual(row["direct_delta"], row["signed_delta"])
                self.assertGreater(row["measurement_information"], 0)
                self.assertGreater(row["nonlinear_residual_energy"], 0)

    def test_exact_linear_S_row_has_zero_nonlinear_charge(self):
        for n in (18, 21):
            with self.subTest(n=n):
                _, Pinv, H, _, Rinv, _, Sinv, K, _, Jj, e, _ = fixture(n)
                row = LEDGER.joseph_signed_energy_terms(
                    Pinv, Jj, H, Rinv, Sinv, K, e, [F(0), F(0), F(0)]
                )
                self.assertEqual(row["nonlinear_residual_energy"], F(0))
                self.assertEqual(row["direct_delta"], -row["measurement_information"])
                self.assertLess(row["direct_delta"], 0)

    def test_exact_reset_congruence_and_combined_row(self):
        for n in (18, 21):
            with self.subTest(n=n):
                _, Pinv, H, _, Rinv, _, Sinv, K, Pj, Jj, e, eta = fixture(n)
                joseph = LEDGER.joseph_signed_energy_terms(
                    Pinv, Jj, H, Rinv, Sinv, K, e, eta
                )
                t = joseph["tangent_posterior"]

                d = [F(1, 20), F(-1, 25), F(1, 30)]
                G = eye(n)
                G3 = add(eye(3), [[F(1, 2) * x for x in row] for row in skew(d)])
                for i in range(3):
                    G[i][:3] = G3[i]
                Ginv = inv(G)
                rho = [F(0) for _ in range(n)]
                rho[0], rho[1], rho[2] = F(1, 1000), F(-1, 1200), F(1, 1500)

                reset = LEDGER.reset_signed_energy_terms(Jj, Ginv, t, rho)
                self.assertEqual(reset["identity_residual"], F(0))
                self.assertEqual(reset["direct_delta"], reset["signed_delta"])

                Pr = mul(mul(G, Pj), tr(G))
                er = [a + b for a, b in zip(mv(G, t), rho)]
                direct_total = quad(inv(Pr), er) - quad(Pinv, e)
                signed_total = joseph["signed_delta"] + reset["signed_delta"]
                self.assertEqual(direct_total, signed_total)

    def test_prediction_and_psd_floor_do_not_increase_moving_energy(self):
        # Exact small fixture for the two non-Joseph event classes used by the
        # signed ledger.  It is algebra only, not a source word.
        n = 18
        P = [[F(0) for _ in range(n)] for _ in range(n)]
        Q = [[F(0) for _ in range(n)] for _ in range(n)]
        D = [[F(0) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            P[i][i] = F(2 + (i % 3), 1)
            Q[i][i] = F(1, 10 + i)
            D[i][i] = F(1, 20 + i)
        Fm = eye(n)
        Fm[6][15], Fm[9][6], Fm[12][9] = F(1, 20), F(1, 30), F(1, 40)
        e = [F((i % 5) - 2, 40 + i) for i in range(n)]
        ep = mv(Fm, e)
        Pp = add(mul(mul(Fm, P), tr(Fm)), Q)
        self.assertLessEqual(quad(inv(Pp), ep), quad(inv(P), e))
        Pf = add(Pp, D)
        self.assertLessEqual(quad(inv(Pf), ep), quad(inv(Pp), ep))


if __name__ == "__main__":
    unittest.main()
