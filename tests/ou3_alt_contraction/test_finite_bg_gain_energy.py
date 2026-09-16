from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import os
import tempfile
import unittest
from tools.stability.ou3_alt_contraction import finite_bg_gain_energy as X


def diagonal(x): return tuple(tuple(F(x) if i==j else F(0) for j in range(3)) for i in range(3))


def fixture(posterior=F(1,2)):
    p=diagonal(1);zero=(F(0),)*3
    innovation=(F(1),F(-2),F(1,2));mean=tuple(x/2 for x in innovation)
    return (X.Row(0,0,p,zero),X.Row(1,1,p,zero),X.Row(1,2,p,zero),
            X.Row(1,12,p,zero,diagonal(F(1,2)),diagonal(1),innovation),
            X.Row(1,13,p,mean),X.Row(1,14,diagonal(posterior),mean),
            X.Row(1,100,diagonal(posterior),mean))


class Tests(unittest.TestCase):
    def test_exact_energy_telescopes_without_PSD_nonincrease(self):
        out=X.account(fixture(F(9,8)),steps=1)
        self.assertEqual(out['gain_energy'],F(3,4))
        self.assertEqual(out['Joseph_positive_defect_sum'],F(9,8))
        self.assertEqual(out['telescoping_energy_budget'],F(3,4))
        self.assertFalse(out['floating_Joseph_PSD_nonincrease_assumed'])
        self.assertTrue(out['same_history_telescoping_and_Cauchy_bound_closed'])
        self.assertFalse(out['source_uniform_defect_and_innovation_bounds_closed'])

    def test_hidden_suffix_state_change_or_noise_replacement_rejected(self):
        trace=list(fixture());trace[-1]=replace(trace[-1],bias=(F(9),F(0),F(0)))
        with self.assertRaisesRegex(ValueError,'unaccounted'):X.account(trace,steps=1)
        trace=list(fixture());trace[3]=replace(trace[3],noise=diagonal(0))
        with self.assertRaisesRegex(ValueError,'positive diagonal'):X.account(trace,steps=1)
        trace=list(fixture());trace[3]=replace(trace[3],covariance=diagonal(2))
        with self.assertRaisesRegex(ValueError,'detached'):X.account(trace,steps=1)

    def test_injection_roundoff_is_retained(self):
        trace=list(fixture());extra=F(1,2**23)
        for i in (4,5,6):trace[i]=replace(trace[i],bias=(trace[i].bias[0]+extra,*trace[i].bias[1:]))
        out=X.account(trace,steps=1)
        self.assertEqual(out['bias_norm_additive_upper'],extra)
        self.assertEqual(out['innovation_energy'],F(21,4))


@unittest.skipUnless(os.environ.get('OU3_ALT_REQUIRE_NATIVE')=='1','explicit native feasibility requested')
class NativeTests(unittest.TestCase):
    def test_actual_startup_H_and_A_words_account_for_all_bg_events(self):
        with tempfile.TemporaryDirectory() as td:
            for mode,scenario in (('H',0),('A',2)):
                r=X.run_native(Path(td),mode=mode,scenario=scenario)
                self.assertEqual(r['steps'],600)
                self.assertGreaterEqual(r['gain_events'],600)
                self.assertTrue(r['native']['guard_dormant_in_this_trace'])
                self.assertGreater(r['final_bg_trace'],0)
                self.assertLess(r['gain_energy'],F(1,100000))
                self.assertFalse(r['all_event_arithmetic_source_uniform'])


if __name__=='__main__':unittest.main()
