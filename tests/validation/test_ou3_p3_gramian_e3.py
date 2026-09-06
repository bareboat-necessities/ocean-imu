#!/usr/bin/env python3
import math,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; TOOLS=ROOT/"tools"; sys.path.insert(0,str(TOOLS)); sys.path.insert(0,str(TOOLS/"stability"))
import ou3_p3_gramian_e3 as E3
class GramianE3Tests(unittest.TestCase):
 def test_e3_bound_is_strict_and_sharpen_probe_preserves_stronger_result(self):
  upper=[1.,2.,3.,4.]; sharp0=E3.sharpen(upper,.5,1e-12,1e-3); probe={"horizon_s":.5,"unit_gramian_det_lower":1e-12,"q_c_min_lower":1e-3,"relative_process_floor_lower":sharp0["relative_process_floor_lower"]/2}; sharp=E3.sharpen_probe(probe,upper); self.assertEqual(E3.validate(sharp),[]); self.assertGreater(sharp["relative_process_floor_lower"],probe["relative_process_floor_lower"]); self.assertGreater(sharp["improvement_over_trace3"],1.)
 def test_common_covariance_scaling_still_has_inverse_effect(self):
  a=E3.sharpen([1.,2.,3.,4.],.5,1e-12,1e-3); b=E3.sharpen([10.,20.,30.,40.],.5,1e-12,1e-3); ratio=a["relative_process_floor_lower"]/b["relative_process_floor_lower"]; self.assertGreater(ratio,9.999999999); self.assertLess(ratio,10.000000001)
 def test_no_numerical_eigendecomposition_is_used(self):
  d=E3.sharpen([1.,1.,1.,1.],.5,1e-12,1e-3); self.assertFalse(d["numerical_eigendecomposition_used"]); self.assertTrue(math.isfinite(d["relative_process_floor_lower"])); self.assertGreater(d["relative_process_floor_lower"],0.)
if __name__=="__main__": unittest.main()
