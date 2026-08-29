from pathlib import Path
import sys,unittest
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools'))
import ou3_p4_support_aware_metric_certificate as C
class T(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.d=C.build()
 def test_pass(self):self.assertEqual(C.validate(self.d),[]);self.assertEqual(self.d['P4_SUPPORT_AWARE_METRIC_CERTIFICATE'],'PASS')
 def test_monotone(self):
  for mode in ('H','A'):
   m=self.d['modes'][mode];self.assertGreaterEqual(m['certified_level_W'],m['certified_level_W_before_support']);self.assertLessEqual(m['attitude_support_metric_lambda_max_upper'],m['global_metric_lambda_max_upper_previous'])
if __name__=='__main__':unittest.main()