import unittest,json
from pathlib import Path
from tools.stability.ou3_theorem.source_uniform_root import root_cells
from tools.stability.ou3_theorem.history_interpolation import parse_knots,lipschitz_outer
class HistoryInterpolationTests(unittest.TestCase):
 def test_slow_outer_uses_all_knots(self):
  c=root_cells(60.)[0];k=parse_knots(c);x=lipschitz_outer(k[("slow_gyro",0)],15,1e-5,.02)
  self.assertLessEqual(x.hi,.02);self.assertGreaterEqual(x.lo,-.02)
 def test_generated_outputs_absent_from_root(self):
  names={k for k,_ in root_cells(60.)[0].coordinates}
  for x in ("tau","sigma_aw","R_S","covariance","gain"):self.assertNotIn(x,names)
if __name__=="__main__":unittest.main()
