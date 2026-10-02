import unittest
from tools.stability.ou3_theorem.temporal_source_domain import Symbol,build_temporal_domain,nonlinear_side_constraints
class TemporalSourceDomainTests(unittest.TestCase):
 def test_slow_rate_and_fast_windows_present(self):
  s=[]
  for k,t in enumerate((0.,1.,2.)):
   for a in range(3):
    s += [Symbol(f"sa{k}{a}","slow_accel",a,t),Symbol(f"fa{k}{a}","fast_accel",a,t)]
  d=build_temporal_domain(tuple(s),2);A,b=d.matrices()
  self.assertGreater(A.shape[0],len(s));self.assertEqual(A.shape[1],len(s))
 def test_marine_and_mag_side_constraints_attached(self):
  z=nonlinear_side_constraints(())
  self.assertEqual(z["marine"]["attitude_span"]["T"],30)
  self.assertEqual(z["magnetic_service"]["T"],1)
if __name__=="__main__":unittest.main()
