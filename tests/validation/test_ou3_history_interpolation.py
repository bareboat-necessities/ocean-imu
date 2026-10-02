# ruff: noqa: F401, F811
import unittest,json
from pathlib import Path
from tools.stability.ou3_theorem.source_uniform_root import root_cells
from tools.stability.ou3_theorem.history_interpolation import parse_knots,lipschitz_outer,acceleration_outer_from_velocity_knots
class HistoryInterpolationTests(unittest.TestCase):
 def test_slow_outer_uses_all_knots(self):
  c=root_cells(60.)[0];k=parse_knots(c);x=lipschitz_outer(k[("slow_gyro",0)],15,1e-5,.02)
  self.assertLessEqual(x.hi,.02);self.assertGreaterEqual(x.lo,-.02)
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
