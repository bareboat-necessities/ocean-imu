"""Numerical equivalence and conditioning of the finite diagnostic reader."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.stability.ou3_theorem.information_ratio_source_diagnostic import (
    psd_root, replay, square_root_replay)


class InformationRatioReplayTest(unittest.TestCase):
    def test_literal_chronology_matches_joseph_replay(self):
        rng = np.random.default_rng(914)
        n = 21
        a = rng.normal(size=(n, n))
        prior = np.eye(n)+a@a.T
        ops = []
        for _ in range(12):
            f = np.eye(n)+rng.normal(size=(n, n))*.01
            q = rng.normal(size=(n, 4))*.001
            sync = rng.normal(size=(n, 3))*.0001
            h = rng.normal(size=(3, n))
            reset = np.eye(n)+rng.normal(size=(n, n))*.002
            ops.extend([('P', f, q@q.T), ('Y', sync@sync.T),
                        ('C', h, np.eye(3)*.03, 'A'), ('G', reset)])
        np.testing.assert_allclose(square_root_replay(prior, ops),
                                   replay(prior, ops)[0], rtol=1e-10, atol=1e-12)

    def test_diffuse_correlated_prior_matches_information_solution(self):
        rng = np.random.default_rng(82)
        basis, _ = np.linalg.qr(rng.normal(size=(21, 21)))
        prior = (basis*np.geomspace(1, 1e10, 21))@basis.T
        ops = [('C', np.eye(21)[i:i+3], np.eye(3)*.02, 'M')
               for i in range(0, 21, 3)]
        expected = np.linalg.inv(np.linalg.inv(prior)+50*np.eye(21))
        np.testing.assert_allclose(square_root_replay(prior, ops), expected,
                                   rtol=1e-5, atol=1e-10)

    def test_indefinite_covariance_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'indefinite exported covariance'):
            psd_root(np.diag([1.0, -1e-5]))


if __name__ == '__main__':
    unittest.main()
