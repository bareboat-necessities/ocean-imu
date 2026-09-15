import unittest
from dataclasses import replace
from tools.stability.ou3_alt_contraction import finite_startup_wpe_machine_history as X
from tools.stability.ou3_alt_contraction import finite_wpe_machine_binary32 as WM
from tools.stability.ou3_alt_contraction import finite_wpe_moment_binary32 as M
import test_finite_startup_joined_machine_history as BASE


def expw(arg):
    lo,hi=M.EXP.exp_minus_enclosure(arg); return M.B.rn32((lo+hi)/2)


def mode_witness(machine_state,vertical,dt):
    cfg=machine_state.cfg; mom=machine_state.separate
    d=expw(M.B.mul(cfg.lambda_,M.B.rn32(dt)))
    gain=M.B.div(M.B.sub(M.ONE,d),cfg.lambda_)
    hp1=M.B.mul(d,M.B.sub(M.B.add(mom.hp1,vertical),mom.accel_prev))
    hp2=M.B.mul(d,M.B.sub(M.B.add(mom.hp2,hp1),mom.hp1_prev))
    vel=M._sum_products(d,mom.velocity,gain,hp2)[0]
    elev=M._sum_products(d,mom.elevation,gain,vel)[0]
    return WM.ModeWitnesses(dict(decay_exp=d,velocity_successor=vel,elevation_successor=elev))


class Tests(unittest.TestCase):
    def test_lower_frequency_branch_cannot_ignore_machine_latches(self):
        b=BASE.initial(); wpe=WM.initial(b.runtime.wpe_cfg)
        for changed in ('separate_usable','fma_usable'):
            with self.assertRaisesRegex(ValueError,'literal reset'):
                X.State(b,replace(wpe,**{changed:True}))


    def test_two_cold_samples_attach_full_wpe_to_same_mahony_vertical(self):
        b=BASE.initial(); s=X.State(b,WM.initial(b.runtime.wpe_cfg))
        for _ in range(2):
            raw,kw=BASE.cold_operands(s.base)
            guard=BASE.X.GUARD.step(s.base.guard,s.base.guard_cfg,dt=kw['machine_dt'],raw_gyro=kw['machine_gyro_body'],raw_acc=kw['machine_acc_body'],**kw['guard_witnesses'])
            sep,_=BASE.VSBASE.source_witness(s.base.separate_source,guard,s.base.runtime,kw['dt'])
            sw=mode_witness(s.wpe,sep.band_input,kw['dt']); fw=WM.ModeWitnesses(dict(sw.moment))
            out=X.step(s,raw,separate_wpe=sw,fma_wpe=fw,**kw)
            self.assertEqual(out.state.wpe.logs,out.state.base.base.lower.wpe)
            self.assertEqual(out.separate_wpe.vertical_accel,out.lower.separate_source.band_input)
            s=out.state
        self.assertEqual(s.wpe.separate.samples,2)
        self.assertEqual(s.wpe.fma.samples,2)

    def test_independent_first_period_and_takeover_survive_startup_product(self):
        from test_finite_admitted_wpe_machine_clock_interleaved_prefix import general_witness
        from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as WF
        b=BASE.initial(); s=X.State(b,WM.initial(b.runtime.wpe_cfg))
        raw,kw=BASE.cold_operands(s.base)
        guard=BASE.X.GUARD.step(s.base.guard,s.base.guard_cfg,dt=kw['machine_dt'],raw_gyro=kw['machine_gyro_body'],raw_acc=kw['machine_acc_body'],**kw['guard_witnesses'])
        src,_=BASE.VSBASE.source_witness(s.base.separate_source,guard,s.base.runtime,kw['dt'])
        w=mode_witness(s.wpe,src.band_input,kw['dt'])
        s=X.step(s,raw,separate_wpe=w,fma_wpe=w,**kw).state
        # Isolate the composed implication with distinct moment histories.
        # This component predecessor is not asserted to arise from that one sample.
        sm=replace(s.wpe.separate,elapsed=M.B.rn32(4),weight=M.ONE,
                   velocity_sq=M.B.rn32(9),elevation_sq=M.ONE)
        s=replace(s,wpe=replace(s.wpe,separate=sm))
        raw,kw=BASE.cold_operands(s.base)
        guard=BASE.X.GUARD.step(s.base.guard,s.base.guard_cfg,dt=kw['machine_dt'],raw_gyro=kw['machine_gyro_body'],raw_acc=kw['machine_acc_body'],**kw['guard_witnesses'])
        src,_=BASE.VSBASE.source_witness(s.base.separate_source,guard,s.base.runtime,kw['dt'])
        sw=general_witness(s.wpe,src.band_input,kw['dt'],'separate')
        fw=general_witness(s.wpe,src.band_input,kw['dt'],'fma')
        out=X.step(s,raw,separate_wpe=sw,fma_wpe=fw,**kw)
        self.assertIsNone(out.state.base.base.lower.frontend.tuner.wpe.log_period)
        self.assertIsNotNone(out.state.wpe.logs.separate.log_period)
        self.assertIsNone(out.state.wpe.logs.fma.log_period)
        self.assertTrue(out.state.wpe.separate_usable)
        self.assertFalse(out.state.wpe.fma_usable)
        self.assertEqual(out.state.wpe.logs,out.state.base.base.lower.wpe)
        held,_=X.boundary(out.state)
        self.assertIs(held.wpe,out.state.wpe)
        # Read the next sample's frequencies directly from those retained histories.
        sf,ff=WF.machine_frequencies(held.base.base.lower.frontend.tuner.wpe,held.wpe,logs=held.wpe.logs,
            separate_getter=WF.FrequencyExp(-held.wpe.logs.separate.log_period,M.B.rn32(1)),
            fma_getter=None,shadow_frequency=None,stats_cfg=held.base.runtime.stats_cfg,
            exact_min_hz=held.base.runtime.candidate_cfg.min_freq,exact_max_hz=held.base.runtime.candidate_cfg.max_freq)
        self.assertEqual((sf.external.branch,ff.external.branch),('wpe','prior'))

    def test_log_witness_cannot_be_supplied_outside_full_wpe_product(self):
        b=BASE.initial(); s=X.State(b,WM.initial(b.runtime.wpe_cfg)); raw,kw=BASE.cold_operands(b)
        with self.assertRaisesRegex(TypeError,'owned by the full WPE'):
            X.step(s,raw,separate_wpe=WM.ModeWitnesses({}),fma_wpe=WM.ModeWitnesses({}),
                   separate_log_witness=None,**kw)

    def test_readiness_closes_startup_vertical_ancestry_not_branch_robustness(self):
        r=X.readiness()
        self.assertTrue(r['startup_frontend_vertical_ancestry_attached'])
        self.assertTrue(r['goLive_preserves_full_WPE_moment_log_machine_history_by_identity'])
        self.assertFalse(r['machine_vs_exact_WPE_period_branch_agreement_required'])
        self.assertTrue(r['independent_machine_WPE_production_and_frequency_branches_composed'])
        self.assertFalse(r['Live_600_step_WPE_machine_history_attached'])
        self.assertFalse(r['storage_search_allowed'])

if __name__=='__main__': unittest.main()
