import unittest,math
from tools.stability.ou3_theorem.causal_tuner_interval import *
class CausalTunerIntervalTests(unittest.TestCase):
 def point(self,x):return I(x,x)
 def test_mahony_quaternion_norm_dependency_does_not_create_zero(self):
  z=I(-1,1);m=MahonyBox((z,z,z,z),(I(0,0),I(0,0),I(0,0)),True,1.0)
  # Deliberately decorrelated component box contains q=0, while the coupled
  # invariant states the represented quaternion is normalized.
  out=m.step(I(.005,.005),(I(0,0),I(0,0),I(0,0)),
             (I(0,0),I(0,0),I(9.8,9.8)),.2,.02,9.8)
  self.assertTrue(math.isfinite(out.lo));self.assertEqual(m.q_norm_lower,1.0)
 def test_interval_division_contains_corner_quotients(self):
  z=I(2,4)/I(2,3)
  self.assertLessEqual(z.lo,2/3);self.assertGreaterEqual(z.hi,2)
  with self.assertRaises(ArithmeticError): _=I(1,2)/I(-1,1)
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
 def test_wave_period_prior_is_generated_until_usable(self):
  w=initial_wave_period_box()
  self.assertEqual(w.frequency(.2),(I(.2,.2)))
 def test_wave_period_point_state_advances(self):
  w=initial_wave_period_box()
  # Early leak transient has no gate ambiguity and must remain generated.
  for _ in range(10): w.step(I(.005,.005),I(.1,.1))
  self.assertGreater(w.elapsed.lo,0)
  self.assertFalse(w.usable)
 def test_closed_chain_has_no_frequency_input(self):
  import inspect
  sig=inspect.signature(ClosedCausalAdaptationBox.step)
  self.assertNotIn("wave_frequency",sig.parameters)
if __name__=="__main__":unittest.main()
