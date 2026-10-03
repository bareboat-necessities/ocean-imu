import unittest, numpy as np
from tools.stability.ou3_theorem.planar_information_ceilings import acc_ceiling,S_ceiling,mag_ceiling,certificate
class PlanarInformationCeilingsTests(unittest.TestCase):
 def test_psd(self):
  for pair in (acc_ceiling(),S_ceiling(),mag_ceiling()):
   for J in pair:self.assertGreaterEqual(np.linalg.eigvalsh((J+J.T)/2).min(),-1e-10)
 def test_fail_closed_bindings(self):
  c=certificate();self.assertFalse(c["all_time_f_bound_verified"]);self.assertFalse(c["literal_profile_binding_verified"])
if __name__=="__main__":unittest.main()
