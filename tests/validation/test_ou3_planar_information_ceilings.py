import unittest, numpy as np
from tools.stability.ou3_theorem.planar_information_ceilings import acc_ceiling,S_ceiling,mag_ceiling,certificate
class PlanarInformationCeilingsTests(unittest.TestCase):
 def test_psd(self):
  for pair in (acc_ceiling(),S_ceiling(),mag_ceiling()):
   for J in pair:self.assertGreaterEqual(np.linalg.eigvalsh((J+J.T)/2).min(),-1e-10)
 def test_actual_rotated_rows_are_dominated(self):
  from tools.stability.ou3_theorem.planar_information_ceilings import skew
  from tools.stability.ou3_theorem.planar_parity import EVEN,ODD
  for f in (np.array([0.,0.,30.]),np.array([18.,0.,24.]),np.array([0.,0.,0.])):
   for angle in (0.,.3,-1.):
    c,s=np.cos(angle),np.sin(angle);R=np.array([[c,0,s],[0,1,0],[-s,0,c]])
    H=np.zeros((3,21));H[:,:3]=-skew(f);H[:,15:18]=R;H[:,18:21]=np.eye(3)
    for idx,J in zip((EVEN,ODD),acc_ceiling()):
     actual=(H.T@H/.04)[np.ix_(idx,idx)]
     self.assertGreaterEqual(np.linalg.eigvalsh(J-actual).min(),-1e-9)
 def test_old_rank_three_ceiling_has_false_null_direction(self):
  H=np.zeros((3,21));H[1,0]=-30;H[:,15:18]=np.eye(3);H[:,18:21]=np.eye(3)
  x=np.zeros(21);x[0]=1;x[15]=-30
  self.assertEqual(np.linalg.norm(30*x[:3]+x[15:18]+x[18:21]),0)
  self.assertGreater(np.linalg.norm(H@x),40)
 def test_S_floor_includes_actual_axis_factors(self):
  from tools.stability.ou3_theorem.planar_parity import EVEN,ODD
  E,O=S_ceiling()
  self.assertAlmostEqual(E[EVEN.index(12),EVEN.index(12)],1/.108**2)
  self.assertAlmostEqual(O[ODD.index(13),ODD.index(13)],1/.075**2)
 def test_fail_closed_bindings(self):
  c=certificate();self.assertFalse(c["all_time_f_bound_verified"]);self.assertFalse(c["literal_profile_binding_verified"])
if __name__=="__main__":unittest.main()
