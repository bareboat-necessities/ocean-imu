import unittest
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

    def test_log_witness_cannot_be_supplied_outside_full_wpe_product(self):
        b=BASE.initial(); s=X.State(b,WM.initial(b.runtime.wpe_cfg)); raw,kw=BASE.cold_operands(b)
        with self.assertRaisesRegex(TypeError,'owned by the full WPE'):
            X.step(s,raw,separate_wpe=WM.ModeWitnesses({}),fma_wpe=WM.ModeWitnesses({}),
                   separate_log_witness=None,**kw)

    def test_readiness_closes_startup_vertical_ancestry_not_branch_robustness(self):
        r=X.readiness()
        self.assertTrue(r['startup_frontend_vertical_ancestry_attached'])
        self.assertTrue(r['goLive_preserves_full_WPE_moment_log_machine_history_by_identity'])
        self.assertFalse(r['machine_vs_exact_WPE_period_branch_robustness_closed'])
        self.assertFalse(r['Live_600_step_WPE_machine_history_attached'])
        self.assertFalse(r['storage_search_allowed'])

if __name__=='__main__': unittest.main()
