import unittest,math
from tools.stability.ou3_theorem.causal_tuner_interval import *
class CausalTunerIntervalTests(unittest.TestCase):
 def point(self,x):return I(x,x)
 def test_point_band_matches_formula(self):
  b=BandBox(I(0,0),I(0,0),I(0,0),I(0,0),I(0,0))
  y=b.step(I(1,1),I(.005,.005),I(.2,.2))
  self.assertLess(y.hi-y.lo,1e-12)
 def test_variance_nonnegative(self):
  v=VarianceBox(I(0,0),I(0,0),I(0,0),I(0,0))
  for _ in range(10): z=v.update(I(.005,.005),I(-.2,.3),I(.2,.2))
  self.assertGreaterEqual(z.lo,0)
 def test_staging_is_one_sample_late(self):
  t=ShippingTunerBox(I(1.1,1.1),I(.1,.1),I(.5,.5))
  old=t.commit();t.stage(I(.005,.005),I(.2,.2),I(.04,.04),I(.01,.01),1.38,.9,I(.01,.01),I(.01,.01))
  self.assertEqual((t.tau.lo,t.sigma.lo,t.RS.lo),(old[0].lo,old[1].lo,old[2].lo))
  new=t.commit();self.assertNotEqual(new[0].lo,old[0].lo)
 def test_generated_RS_not_independent(self):
  r=spectral_mse_RS(I(1,1),I(.2,.2),I(.02,.02))
  self.assertGreater(r.lo,0);self.assertLess(r.hi-r.lo,1e-12)
 def test_mahony_fails_closed_on_zero_acc_norm(self):
  m=MahonyBox((I(1,1),I(0,0),I(0,0),I(0,0)),(I(0,0),)*3,True)
  with self.assertRaises(ArithmeticError):m.step(I(.005,.005),(I(0,0),)*3,(I(0,0),)*3,.2,.02)
if __name__=="__main__":unittest.main()
