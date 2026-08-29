from pathlib import Path
import sys,unittest
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools'))
import ou3_p4_direct_word_contraction_certificate as C
class T(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.d=C.build()
 def test_pass(self):self.assertEqual(C.validate(self.d),[]);self.assertEqual(self.d['P4_DIRECT_STRICT_WORD_CONTRACTION_CERTIFICATE'],'PASS')
 def test_monotone(self):
  for mode in ('H','A'):
   m=self.d['modes'][mode];self.assertGreaterEqual(m['certified_level_W'],m['certified_level_W_before_direct']);self.assertGreaterEqual(m['direct_W_factor_lower'],1.0);self.assertGreater(m['direct_strict_endpoint_gap_lower'],0.0)
if __name__=='__main__':unittest.main()