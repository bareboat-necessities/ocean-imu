import math, unittest
from tools.stability.ou3_alt_contraction import finite_bias_history_restriction as X

class Tests(unittest.TestCase):
    def test_all_three_admitted_source_families_restrict_one_history(self):
        rs=X.restrictions(); self.assertEqual(tuple(r.name for r in rs),('BIAS0','BIAS1','BIAS2'))
        self.assertTrue(all(r.one_history and r.source_admission_pass for r in rs))
        self.assertTrue(all(not r.hardware_admission_pass for r in rs))
        self.assertTrue(all(r.dt_s==X.DT for r in rs))

    def test_bias0_and_bias1_use_one_fixed_physical_root_not_per_sample_phi(self):
        r0,r1,_=X.restrictions()
        for r in (r0,r1):
            self.assertIsNone(r.canonical_phi)
            self.assertLess(r.phi_lo,r.phi_hi)
            self.assertIn('fixes phi=exp(-dt/tau)',r.phi_rule)

    def test_bias2_can_be_represented_without_any_relaxation_assumption(self):
        *_,r=X.restrictions(); self.assertEqual(r.canonical_phi,1.0); self.assertEqual(r.phi_hi,1.0)
        import ou3_p4_bias2_family as B2
        d=B2.build(); variation=float(d['variation_rate_abs_upper_mps3'])*X.DT
        self.assertLessEqual(math.nextafter(variation,math.inf),r.driver_component_upper)
        self.assertIn('w=beta_k-beta_{k-1}',r.phi_rule)

    def test_driver_is_a_derived_same_history_quantity(self):
        beta0=(.1,-.2,.3); beta1=(.11,-.19,.29); phi=.9997
        w=tuple(beta1[i]-phi*beta0[i] for i in range(3))
        self.assertEqual(tuple(phi*beta0[i]+w[i] for i in range(3)),beta1)

    def test_rational_executable_graph_is_not_promoted_to_all_real_histories(self):
        d=X.build(); self.assertEqual(X.validate(d),[])
        self.assertFalse(d['rational_Fraction_graph_represents_every_real_bias_history'])
        self.assertFalse(d['irrational_exp_phi_scalar_extension_closed'])
        self.assertFalse(d['binary32_exp_and_runtime_phi_correspondence_closed'])
        self.assertFalse(d['bias_history_attached_to_every_finite_shipping_branch'])
        self.assertFalse(d['storage_search_allowed']); self.assertFalse(d['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
