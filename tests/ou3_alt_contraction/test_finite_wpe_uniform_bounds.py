from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_wpe_uniform_bounds as U
from tools.stability.ou3_alt_contraction import finite_wpe_machine_binary32 as W


def config(): return W.SHADOW.WPEConfig(U.LAMBDA,4,F(1,20),20,180)
def rounded_exp(x):
    lo,hi=U.exp_interval(F(x)); return U.B.rn32((lo+hi)/2)


class Tests(unittest.TestCase):
    def test_exact_induction_covers_every_alpha_and_nonproducing_branches(self):
        r=U.build(); self.assertEqual(U.validate(r),[])
        self.assertTrue(all(x>0 for x in r['induction_margins'].values()))
        self.assertGreater(r['log_induction_margin'],0)
        self.assertGreater(r['raw_period_derived_lower'],U.RAW_LO)
        self.assertLess(r['raw_period_derived_upper'],U.RAW_HI)
        self.assertTrue(r['no_positive_variance_or_period_production_assumed'])
        self.assertFalse(r['upstream_observer_totality_or_startup_reachability_closed'])
        self.assertFalse(r['target_libm_and_compiler_correspondence_closed'])
        r['state_bounds']['hp1']=F(1)
        self.assertTrue(U.validate(r))
        self.assertEqual(U.build()['state_bounds']['hp1'],U.HP1)

    def test_eager_frequency_exp_is_checked_even_before_usable_latch(self):
        from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as WF
        from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as ST
        s=W.initial(config(),bounded_profile=True)
        logs=replace(s.logs,separate=replace(s.logs.separate,log_period=F(1)),
                     fma=replace(s.logs.fma,log_period=F(1)))
        s=replace(s,logs=logs)
        exact=W.SHADOW.WPEState(log_period=F(1))
        good=WF.FrequencyExp(F(-1),rounded_exp(-1))
        kw=dict(logs=logs,fma_getter=good,shadow_frequency=None,
                stats_cfg=ST.StatsConfig(4,F(3,10),60,F(1,20),5),
                exact_min_hz=F(3,100),exact_max_hz=F(6,5))
        out=WF.machine_frequencies(exact,s,separate_getter=good,**kw)
        self.assertEqual(out[0].external.exact_shadow_frequency,WF.PRIOR_EXACT)
        legacy=WF.getters(WF.bind_log_state(exact,F(1)),period_exp=rounded_exp(1),frequency_exp=rounded_exp(-1))
        WF.machine_frequencies(exact,s,separate_getter=legacy,**kw)
        for bad in (WF.FrequencyExp(F(-1),None),WF.FrequencyExp(F(-1),F(1))):
            with self.assertRaisesRegex(ValueError,'exp witness detached'):
                WF.machine_frequencies(exact,s,separate_getter=bad,**kw)

    def test_transcendental_relations_reject_detached_and_nonfinite_witnesses(self):
        for x in (F(-48),F(-1),F(0),F(1),F(48)):
            U.check_exp(x,rounded_exp(x))
            with self.assertRaisesRegex(ValueError,'exp witness detached'):
                U.check_exp(x,None)
        for x in (U.RAW_LO,F(1,2),F(1),F(2),U.RAW_HI):
            lo,hi=U.log_interval(x); y=U.B.rn32((lo+hi)/2)
            U.check_log(x,y)
            with self.assertRaisesRegex(ValueError,'log witness detached'):
                U.check_log(x,U.B.add(y,1))

    def test_first_valid_log_really_comes_from_machine_moments(self):
        # Component induction predecessor, not asserted to be a reachable
        # physical startup fixture. All bounded predecessors must be covered.
        from test_finite_admitted_wpe_machine_clock_interleaved_prefix import general_witness
        s=W.initial(config(),bounded_profile=True)
        moment=W.MOM.State(elapsed=U.B.rn32(30),weight=F(1),velocity_sq=F(9),elevation_sq=F(1))
        s=replace(s,separate=moment,fma=moment)
        ws=[]
        for mode in ('separate','fma'):
            w=general_witness(s,F(0),U.DT,mode)
            raw=w.raw_log.raw_period; lo,hi=U.log_interval(raw)
            lr=U.B.rn32((lo+hi)/2)
            w=replace(w,raw_log=W.RawLogBinding(raw,lr),log=W.LOG.InitWitness(lr),
                      usable=W.USABLE.PeriodWitness(lr,rounded_exp(lr)))
            ws.append(w)
        out,*_=W.step(s,dt=U.DT,vertical_accel=0,separate=ws[0],fma=ws[1])
        self.assertTrue(out.bounded_profile)
        U.check_state(out)
        bad=replace(ws[0],raw_log=replace(ws[0].raw_log,log_raw=F(0)),log=W.LOG.InitWitness(F(0)))
        with self.assertRaisesRegex(ValueError,'log witness detached'):
            W.step(s,dt=U.DT,vertical_accel=0,separate=bad,fma=ws[1])

    def test_bounds_are_required_by_the_qualified_runtime_state(self):
        s=W.initial(config(),bounded_profile=True)
        with self.assertRaisesRegex(ValueError,'source input exceeds'):
            W.step(s,dt=U.DT,vertical_accel=F(33),separate=W.ModeWitnesses({}),fma=W.ModeWitnesses({}))
        with self.assertRaisesRegex(ValueError,'uniform hp1 bound'):
            replace(s,separate=replace(s.separate,hp1=2*U.HP1))
        with self.assertRaisesRegex(ValueError,'default construction constants'):
            W.initial(replace(config(),min_horizon_sec=F(1)),bounded_profile=True)


if __name__=='__main__': unittest.main()
