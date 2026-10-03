import unittest,tempfile,struct,numpy as np
from pathlib import Path
from tools.stability.ou3_theorem.planar_phase_profile import read_profile,candidate_radii
class PhaseProfileTests(unittest.TestCase):
 def test_reader(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"p";a=np.zeros((8000,227),dtype="<f4")
   p.write_bytes(b"OU3PRF1\0"+struct.pack("<I",8000)+a.tobytes())
   r=read_profile(p);self.assertEqual(len(r),8000);self.assertEqual(max(candidate_radii(r)["even"]),0)
if __name__=="__main__":unittest.main()
