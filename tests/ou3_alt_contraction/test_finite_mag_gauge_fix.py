"""MagAutoTuner gauge-fix boundary regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_gauge_fix as X


class Tests(unittest.TestCase):
    def test_gauge_fix_uses_same_mean_and_zeroes_horizontal_y(self):
        out=X.gauge_fix((3,4,12),mag_norm_min=1,sqrt=X.HorizontalSqrt(25,5))
        self.assertEqual(out.world_reference,(5,0,12)); self.assertTrue(out.ready)

    def test_sqrt_witness_cannot_detach_from_mean(self):
        with self.assertRaisesRegex(ValueError,'detached'):
            X.gauge_fix((3,4,12),mag_norm_min=1,sqrt=X.HorizontalSqrt(16,4))

    def test_ready_gate_is_strict_norm_threshold(self):
        at=X.gauge_fix((3,4,0),mag_norm_min=5,sqrt=X.HorizontalSqrt(25,5))
        above=X.gauge_fix((3,4,1),mag_norm_min=5,sqrt=X.HorizontalSqrt(25,5))
        self.assertFalse(at.ready); self.assertTrue(above.ready)

    def test_readiness_keeps_mean_and_binary32_open(self):
        r=X.readiness()
        self.assertTrue(r['MagAutoTuner_horizontal_gauge_fix_materialized'])
        self.assertFalse(r['accepted_window_mean_runtime_ancestry_attached'])
        self.assertFalse(r['horizontal_sqrt_binary32_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
