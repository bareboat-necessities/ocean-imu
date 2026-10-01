import csv, math, tempfile, unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
from tools.stability.ou3_theorem.bmi270_lf_qualification import design_fir, response, analyze

class LFQualificationTests(unittest.TestCase):
    def test_explicit_diagnostic_filter(self):
        h=design_fir(20,2,1)
        self.assertAlmostEqual(sum(h),1,places=12)
        self.assertGreater(response(h,20,.1),.99)
        with self.assertRaises(ValueError):design_fir(20,10,1)

    def test_zero_finite_record_cannot_qualify_device_or_motion_margin(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'capture.csv'
            fields=['t_s','ax_mps2','ay_mps2','az_mps2','gx_rad_s','gy_rad_s','gz_rad_s',
                    'ax_ref_mps2','ay_ref_mps2','az_ref_mps2','gx_ref_rad_s','gy_ref_rad_s','gz_ref_rad_s']
            with p.open('w',newline='') as f:
                w=csv.writer(f);w.writerow(fields)
                for i in range(300):w.writerow([i/20]+[0]*12)
            r=analyze(p,.01,.001,cutoff_hz=2,transition_hz=1)
            self.assertFalse(r['qualified']);self.assertFalse(r['all_time_continuation_certified'])
            self.assertFalse(r['slow_fast_decomposition_certified'])
            self.assertFalse(r['motion_threshold_inferred'])
            self.assertNotIn('margin_rad',r)
            self.assertGreaterEqual(r['eps_a_LF_mps2'],.01)
            for bad in (math.nan,math.inf,-1):
                with self.assertRaises(ValueError):analyze(p,bad,0,cutoff_hz=2,transition_hz=1)

if __name__=='__main__':unittest.main()
