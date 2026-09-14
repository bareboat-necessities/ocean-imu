"""Adaptive-band corner/exp binary32 coefficient regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_band_coefficients_binary32 as X
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as E
import test_finite_source_bound_live_word as BASE


def cfg(): return BASE.root_state().runtime.band_cfg

def expw(x):
    lo,hi=E.exp_minus_enclosure(x); return B.rn32((lo+hi)/2)


def active_inputs(f=F(1,5),dt=F(1,200)):
    c=cfg(); ff=B.rn32(f); h=B.rn32(dt)
    # Reproduce ordinary corner graph only to obtain the two exp arguments used
    # to construct valid machine witnesses; produce() independently rechecks it.
    guard=B.div(X.NYQUIST_FACTOR,h); upper=min(B.rn32(c.max_hz),guard)
    low=max(B.rn32(c.min_hz),B.mul(B.rn32(c.low_ratio),ff)); low=min(low,B.div(upper,X.SPACING))
    high=min(upper,B.mul(B.rn32(c.high_ratio),ff)); high=max(high,B.mul(low,X.SPACING)); high=min(high,upper)
    xl=B.mul(B.mul(X.TWO_PI,low),h); xh=B.mul(B.mul(X.TWO_PI,high),h)
    return ff,h,expw(xl),expw(xh)


class Tests(unittest.TestCase):
    def test_active_shipping_corner_graph_binds_two_distinct_exp_calls(self):
        f,h,el,eh=active_inputs(); out=X.produce(cfg(),f_ref=f,dt=h,exp_low=el,exp_high=eh)
        self.assertTrue(out.active); self.assertGreater(out.high,out.low)
        self.assertEqual(out.q_low,B.sub(B.rn32(1),B.sub(B.rn32(1),el)))
        self.assertEqual(out.q_high,B.sub(B.rn32(1),B.sub(B.rn32(1),eh)))
        self.assertNotEqual(out.x_low,out.x_high)

    def test_detached_exp_result_fails_closed(self):
        f,h,el,eh=active_inputs()
        with self.assertRaisesRegex(ValueError,'low-corner exp witness detached'):
            X.produce(cfg(),f_ref=f,dt=h,exp_low=B.rn32(F(1,2)),exp_high=eh)

    def test_inactive_upper_limit_consumes_no_exp_witnesses(self):
        c=cfg(); f=B.rn32(F(1,5)); h=B.rn32(100)
        out=X.produce(c,f_ref=f,dt=h)
        self.assertFalse(out.active); self.assertIsNone(out.q_low); self.assertIsNone(out.q_high)
        with self.assertRaisesRegex(ValueError,'consumes no exp witnesses'):
            X.produce(c,f_ref=f,dt=h,exp_low=B.rn32(1),exp_high=B.rn32(1))

    def test_readiness_closes_graph_not_platform_or_reference_frequency(self):
        r=X.readiness()
        for k in ('shipping_band_corner_and_decay_source_shape_matches','corner_clamp_and_nyquist_guard_binary32_graph_materialized',
                  'two_pi_corner_dt_arguments_binary32_source_order_materialized','low_and_high_exp_calls_retained_distinct',
                  'exp_results_bound_to_same_arguments_by_tight_real_enclosure_and_RNE_cell','source_two_subtraction_alpha_q_rounding_materialized'):
            self.assertTrue(r[k])
        for k in ('target_exp_libm_correspondence_closed','reference_frequency_machine_production_closed',
                  'persistent_band_machine_history_composed','source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
