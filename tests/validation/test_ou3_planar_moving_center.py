import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from tools.stability.ou3_theorem.planar_moving_center import diagnose, physical_arc_coordinates, riccati_center, verify_report
from tools.stability.ou3_theorem.planar_parity import EVEN, ODD


def pack(P):
    return np.r_[P[np.ix_(EVEN, EVEN)].ravel(), P[np.ix_(ODD, ODD)].ravel()]


class MovingCenterTests(unittest.TestCase):
    def test_curved_physical_arc_extrema_preserve_one_parameter(self):
        rng = np.random.default_rng(42)
        A = rng.normal(size=(21, 21)); P = A @ A.T+np.eye(21)
        k = 40123; g = 9.80665
        out = physical_arc_coordinates(P, k)
        psi = .02*math.sin(math.pi*(k*.005)/10)
        r = np.zeros(21); r[:3] = [-math.cos(psi), 0., -math.sin(psi)]; r[19] = g
        L = np.linalg.cholesky(P); u = np.linalg.solve(L, r); u /= np.linalg.norm(u)
        for beta in np.linspace(-2*math.atan(1/200), 2*math.atan(1/200), 101):
            e = beta*r
            e[18:21] += g*(math.cos(beta)-1)*np.array([-math.sin(psi), 0., math.cos(psi)])
            e[19] += g*(math.sin(beta)-beta)
            e = np.linalg.solve(L, e)
            self.assertLessEqual(abs(u @ e), out['central_physical_chart_gauge_amplitude']+1e-15)
            self.assertLessEqual(np.linalg.norm(e-u*(u @ e)), out['central_physical_chart_transverse_curvature']+1e-15)

    def test_linked_gain_not_carried_independent_gain(self):
        P = np.diag(np.linspace(.5, 1.5, 21))
        H = np.zeros((3, 21)); H[:, :3] = np.eye(3)
        R = np.eye(3)*.3
        expected = np.linalg.inv(np.linalg.inv(P)+H.T @ np.linalg.solve(R, H))
        np.testing.assert_allclose(riccati_center(P, H, R), expected, rtol=1e-14, atol=1e-14)

    def test_center_is_inherited_and_gap_is_rejected(self):
        I = np.eye(21)
        events = []
        for k in (1, 2):
            # Observed samples are deliberately inconsistent with the reference
            # map, so reseeding would incorrectly erase the cumulative defect.
            a = np.zeros(912)
            a[:225] = pack(k*I); a[225:450] = pack(I); a[450:675] = pack(I)
            a[675:900] = pack((k+1)*I); a[910] = .1
            b = np.zeros(293); b[:225] = pack((k+2)*I); b[282] = .1
            events.extend(((1, k, a), (7, k, b)))
        with tempfile.NamedTemporaryFile() as f:
            with patch('tools.stability.ou3_theorem.planar_moving_center.records', return_value=iter(events)):
                report = diagnose(Path(f.name), word_samples=1)
            self.assertAlmostEqual(report['words'][1]['inherited_center_defect'], np.sqrt(21)/4)
            self.assertFalse(report['center_reseeded_at_word_boundaries'])
            gap = [(kind, k if k == 1 else 3, a) for kind, k, a in events]
            with patch('tools.stability.ou3_theorem.planar_moving_center.records', return_value=iter(gap)):
                with self.assertRaisesRegex(ValueError, 'sample gap'):
                    diagnose(Path(f.name), word_samples=1)

    def test_source_bound_report_cannot_promote_or_hide_coefficient_oracle(self):
        root = Path(__file__).resolve().parents[2]/'reports/results/ou3_stability'
        c = json.loads((root/'planar-moving-center-diagnostic.json').read_text())
        stream_hash = json.loads((root/'planar-service-stream-diagnostic.json').read_text())['stream_sha256']
        self.assertTrue(verify_report(c, stream_hash))
        for key, value in (('center_is_autonomous_future_construction', True),
                           ('all_time_magnetic_service_verified', True),
                           ('center_uses_observed_future_H_and_G', False),
                           ('stream_sha256', hashlib.sha256(b'other stream').hexdigest())):
            bad = copy.deepcopy(c); bad[key] = value
            with self.assertRaises(ValueError):
                verify_report(bad, stream_hash)


if __name__ == '__main__':
    unittest.main()
