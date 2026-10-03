# ruff: noqa: F401, F811
import unittest,json
from pathlib import Path
from tools.stability.ou3_theorem.source_uniform_root import root_cells
from tools.stability.ou3_theorem.history_interpolation import parse_knots,lipschitz_outer,acceleration_outer_from_velocity_knots,sample_history
class HistoryInterpolationTests(unittest.TestCase):
 def test_slow_outer_uses_all_knots(self):
  c=root_cells(60.)[0];k=parse_knots(c);x=lipschitz_outer(k[("slow_gyro",0)],15,1e-5,.02)
  self.assertLessEqual(x.hi,.02);self.assertGreaterEqual(x.lo,-.02)
 def test_sample_history_carries_every_gravity_axis_from_shared_knots(self):
  from tools.stability.ou3_theorem.reachable_history_enclosure import HistoryCell,HistoryInterval
  C=json.loads((Path(__file__).resolve().parents[2]/"tools/stability/ou3_theorem/constants.json").read_text())
  root=root_cells(60.)[0]
  gravity=(.6,.8,0.)
  coords=[]
  for name,x in root.coordinates:
   if name.startswith("gravity_dir_"):
    value=gravity[int(name.split("_")[2])]
    x=HistoryInterval(value,value)
   coords.append((name,x))
  cell=HistoryCell(tuple(coords),root.prefix_token)
  knots=parse_knots(cell)
  for row in sample_history(cell,[0.,.5,60.],C):
   for axis,value in enumerate(gravity):
    with self.subTest(time=row["t"],axis=axis):
     actual=row[f"gravity_dir_{axis}"]
     expected=lipschitz_outer(knots[("gravity_dir",axis)],row["t"],
                             C["marine_motion"]["Omega_max_rad_s"],1.)
     self.assertEqual((actual.lo,actual.hi),(expected.lo,expected.hi))
     self.assertLessEqual(actual.lo,value)
     self.assertGreaterEqual(actual.hi,value)
     if row["t"] in (0.,60.):self.assertLess(actual.hi-actual.lo,1e-12)
 def test_generated_outputs_absent_from_root(self):
  names={k for k,_ in root_cells(60.)[0].coordinates}
  for x in ("tau","sigma_aw","R_S","covariance","gain"):self.assertNotIn(x,names)
 def test_acceleration_uses_signed_velocity_secant(self):
  from tools.stability.ou3_theorem.reachable_history_enclosure import HistoryCell,HistoryInterval
  import tools.stability.ou3_theorem.history_interpolation as hi
  import json
  C=json.loads((Path(__file__).resolve().parents[2]/"tools/stability/ou3_theorem/constants.json").read_text())
  coords=[]
  for t,v in ((0.,0.),(1.,1.)):
   for a in range(3):
    x=v if a==0 else 0.
    coords.append((f"physical_v_{a}_{t}",HistoryInterval(x,x)))
  cell=HistoryCell(tuple(coords),"v")
  z=acceleration_outer_from_velocity_knots(cell,0.,0,C)
  self.assertLessEqual(z.lo,1.);self.assertGreaterEqual(z.hi,1.)
if __name__=="__main__":unittest.main()
