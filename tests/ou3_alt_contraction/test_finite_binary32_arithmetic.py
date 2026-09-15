"""Complete finite-lattice regressions at the source/API and arithmetic joins."""
from fractions import Fraction as F
import random
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_binary32_mahony as BITS
from tools.stability.ou3_alt_contraction import finite_machine_accel_guard_binary32 as GUARD
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as STILL_CONFIG
from tools.stability.ou3_alt_contraction import finite_stillness_sigma_binary32 as STILL
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as EXP
from tools.stability.ou3_alt_contraction import finite_scheduler_nextafter_binary32 as NEXT

MIN_SUB=F(1,1<<149)
MIN_NORMAL=F(1,1<<126)


class Tests(unittest.TestCase):
    def test_signed_gradual_underflow_ties(self):
        for sign in (-1,1):
            for index in (0,1,2,3,17,127,(1<<23)-2,(1<<23)-1):
                mid=sign*F(2*index+1,2)*MIN_SUB
                winner=sign*(index+(index%2))*MIN_SUB
                self.assertEqual(B.rn32(mid),winner)
                self.assertEqual(B.rn32(mid-sign*MIN_SUB/16),sign*index*MIN_SUB)
                self.assertEqual(B.rn32(mid+sign*MIN_SUB/16),sign*(index+1)*MIN_SUB)
        self.assertEqual(B.rn32(-MIN_SUB/4),0)  # Rational projection merges signed zeros.
        self.assertEqual(BITS.round_bits(-MIN_SUB/4),BITS.SIGN)

    def test_normal_subnormal_boundary_and_exact_representation(self):
        self.assertEqual(B.rn32(MIN_NORMAL-MIN_SUB/2),MIN_NORMAL)
        self.assertEqual(B.rn32(MIN_NORMAL-3*MIN_SUB/4),MIN_NORMAL-MIN_SUB)
        for x in (0,MIN_SUB,7*MIN_SUB,MIN_NORMAL-MIN_SUB,MIN_NORMAL,-MIN_SUB):
            self.assertTrue(B.is_binary32(x))
        self.assertFalse(B.is_binary32(MIN_SUB/2))
        self.assertFalse(B.is_binary32(F(3,2)*MIN_SUB))

    def test_cross_kernel_and_independent_rounding_cells(self):
        rng=random.Random(528)
        values=[0,MIN_SUB/2,MIN_NORMAL-MIN_SUB/2]
        for _ in range(300):
            bits=rng.randrange(BITS.MAX_FINITE+1)
            x=BITS.value(bits)
            previous=BITS.value(bits-1) if bits else -MIN_SUB
            following=BITS.value(bits+1) if bits<BITS.MAX_FINITE else BITS.pow2(128)
            values.extend((x,(x+previous)/2,(x+following)/2))
        # Largest-finite upper midpoint overflows and has its own explicit test.
        for x in values:
            for sign in (-1,1):
                real=sign*x
                encoded=BITS.round_bits(real)
                rounded=B.rn32(real)
                self.assertEqual(rounded,BITS.value(encoded))
                self.assertTrue(BITS.rounding_cell_contains(real,encoded))
                self.assertLessEqual(abs(rounded-real),abs(real)*F(1,1<<24)+MIN_SUB/2)

    def test_add_mul_div_fma_and_ema_keep_underflow_branches(self):
        self.assertEqual(B.sub(MIN_NORMAL,MIN_NORMAL-MIN_SUB),MIN_SUB)
        self.assertEqual(B.mul(MIN_NORMAL,F(1,1<<23)),MIN_SUB)
        self.assertEqual(B.mul(MIN_SUB,F(1,2)),0)
        self.assertEqual(B.div(3*MIN_SUB,2),2*MIN_SUB)
        self.assertEqual(B.fma(MIN_SUB,F(1,2),MIN_SUB),2*MIN_SUB)
        # Separate rounding loses a different half-ulp than one fused rounding.
        self.assertEqual(B.ema(3*MIN_SUB,0,F(1,2),contracted=False),MIN_SUB)
        self.assertEqual(B.ema(3*MIN_SUB,0,F(1,2),contracted=True),2*MIN_SUB)

    def test_cell_and_predecessor_helpers_share_the_complete_lattice(self):
        for bits in (1,2,3,17,(1<<23)-1,1<<23,(1<<23)+1):
            value=BITS.value(bits)
            low,high=EXP._positive_rne_cell(value)
            self.assertEqual(low,(value+BITS.value(bits-1))/2)
            self.assertEqual(high,(value+BITS.value(bits+1))/2)
            self.assertEqual(NEXT.predecessor_positive(value),BITS.value(bits-1))
            for endpoint in (low,high):
                self.assertEqual(EXP._interval_hits_rne_cell(endpoint,endpoint,value),not bits%2)
            self.assertTrue(EXP._interval_hits_rne_cell(low,high,value))

    def test_source_api_and_guard_accept_tiny_finite_components(self):
        exact=(MIN_SUB,-MIN_SUB,MIN_SUB/4)
        stored=(MIN_SUB,-MIN_SUB,F(0))
        self.assertEqual(GUARD.api_vec(exact,stored,'gyro'),stored)
        self.assertLess(sum(x*x for x in exact),F(1,50)**2)
        acceleration=(F(0),F(0),B.rn32(F(196133,20000)))
        result=GUARD.step(GUARD.State(),GUARD.Config(),raw_gyro=stored,
            raw_acc=acceleration,dt=B.rn32(F(1,200)))
        self.assertEqual(result.raw_gyro,stored)
        self.assertEqual(result.conditioned_acc,acceleration)

    def test_normal_stillness_input_may_have_subnormal_squared_energy(self):
        cfg=STILL_CONFIG.Config()
        tiny_normal=F(1,1<<62)
        h=B.rn32(F(1,200))
        normalized=B.div(tiny_normal,B.rn32(cfg.gravity))
        instantaneous=B.mul(normalized,normalized)
        self.assertGreater(instantaneous,0)
        self.assertLess(instantaneous,MIN_NORMAL)
        energy=B.mul(B.rn32(cfg.energy_alpha),instantaneous)
        lo,hi=EXP.exp_minus_enclosure(h)
        out=STILL.step(STILL.State(),cfg,vertical_lp=tiny_normal,dt=h,
                       energy_successor=energy,attenuation_exp=B.rn32((lo+hi)/2))
        self.assertTrue(out.state.is_still)
        self.assertEqual(out.inst_energy,instantaneous)
        self.assertEqual(out.state.energy,energy)

    def test_overflow_and_division_by_zero_still_fail_closed(self):
        maximum=BITS.value(BITS.MAX_FINITE)
        threshold=(maximum+BITS.pow2(128))/2
        self.assertEqual(B.rn32(threshold-1),maximum)
        for x in (threshold,-threshold,BITS.pow2(128)):
            with self.assertRaisesRegex(ValueError,'overflow'): B.rn32(x)
        with self.assertRaises(ZeroDivisionError): B.div(MIN_SUB,0)
        readiness=B.readiness()
        self.assertTrue(readiness['finite_binary32_gradual_underflow_RNE_exact'])
        self.assertFalse(readiness['target_gradual_underflow_correspondence_qualified'])
        self.assertFalse(readiness['subnormal_nan_inf_scope_closed'])


if __name__=='__main__': unittest.main()
