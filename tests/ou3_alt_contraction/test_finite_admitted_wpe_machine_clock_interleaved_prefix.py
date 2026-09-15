from fractions import Fraction as F
import unittest
from tools.stability.ou3_alt_contraction import finite_admitted_wpe_machine_clock_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_admitted_machine_clock_qualified_interleaved_prefix as CLOCK
from tools.stability.ou3_alt_contraction import finite_admitted_machine_joined_frontend_racc_interleaved_prefix as JOIN
from tools.stability.ou3_alt_contraction import finite_wpe_machine_binary32 as WM
from tools.stability.ou3_alt_contraction import finite_wpe_moment_binary32 as M
import test_finite_admitted_machine_joined_frontend_racc_interleaved_prefix as BASE
import test_finite_live_magnetic_word as MAG


def expw(arg):
    lo,hi=M.EXP.exp_minus_enclosure(arg); return M.B.rn32((lo+hi)/2)


def witness(wpe,vertical,dt):
    cfg=wpe.cfg; mom=wpe.separate; h=M.B.rn32(dt)
    d=expw(M.B.mul(cfg.lambda_,h)); gain=M.B.div(M.B.sub(M.ONE,d),cfg.lambda_)
    hp1=M.B.mul(d,M.B.sub(M.B.add(mom.hp1,vertical),mom.accel_prev))
    hp2=M.B.mul(d,M.B.sub(M.B.add(mom.hp2,hp1),mom.hp1_prev))
    vel=M._sum_products(d,mom.velocity,gain,hp2)[0]; elev=M._sum_products(d,mom.elevation,gain,vel)[0]
    period=WM.HorizonPeriodWitness(wpe.logs.separate.log_period,M.B.rn32(1))
    horizon=min(max(M.B.mul(cfg.moment_horizon_periods,period.period_exp),cfg.min_horizon_sec),cfg.max_horizon_sec)
    md=expw(M.B.div(h,horizon)); a=M.B.sub(M.ONE,md)
    ms=M.MomentSuccessors(M._ema(mom.weight,M.ONE,a)[0],M._ema(mom.velocity_mean,vel,a)[0],
        M._ema(mom.velocity_sq,M.B.mul(vel,vel),a)[0],M._ema(mom.elevation_mean,elev,a)[0],
        M._ema(mom.elevation_sq,M.B.mul(elev,elev),a)[0])
    vm=M.B.div(ms.velocity_mean,ms.weight); em=M.B.div(ms.elevation_mean,ms.weight)
    vs=M.B.div(ms.velocity_sq,ms.weight); es=M.B.div(ms.elevation_sq,ms.weight)
    vv=max(M.ZERO,M.B.sub(vs,M.B.mul(vm,vm))); ev=max(M.ZERO,M.B.sub(es,M.B.mul(em,em)))
    return WM.ModeWitnesses(dict(decay_exp=d,velocity_successor=vel,elevation_successor=elev,
        moment_decay_exp=md,moment_successors=ms,velocity_var_successor=vv,elevation_var_successor=ev),horizon=period)


class Tests(unittest.TestCase):
    def test_one_live_edge_advances_same_full_wpe_history(self):
        base,kw,join,_,_,_=BASE.fixture(); clock=CLOCK.begin(base)
        mt=JOIN._mtune_state(base.base)
        cfg=WM.MOM.Config.from_shadow(JOIN._runtime(base.base).wpe_cfg)
        mom=M.State(elapsed=M.B.rn32(3),samples=0)
        wpe=WM.State(cfg,mom,mom,mt.base.wpe); s=X.begin(clock,wpe)
        probe=JOIN.imu_step(base,separate_racc_accel_ldlt=MAG.REJECT,fma_racc_accel_ldlt=MAG.REJECT,**join,**kw)
        sw=witness(wpe,probe.separate_source.band_input,join['machine_dt'])
        fw=WM.ModeWitnesses(dict(sw.moment),horizon=WM.HorizonPeriodWitness(wpe.logs.fma.log_period,M.B.rn32(1)))
        out=X.imu_step(s,separate_wpe=sw,fma_wpe=fw,separate_racc_accel_ldlt=MAG.REJECT,fma_racc_accel_ldlt=MAG.REJECT,**join,**kw)
        self.assertEqual(out.state.wpe_steps,1)
        self.assertEqual(out.state.wpe.logs,X._logs(out.state.base))
        self.assertEqual(out.separate_wpe.vertical_accel,out.lower.lower.separate_source.band_input)

    def test_complete_cannot_skip_wpe_edges(self):
        base,_,_,_,_,_=BASE.fixture(); clock=CLOCK.begin(base); mt=JOIN._mtune_state(base.base)
        cfg=WM.MOM.Config.from_shadow(JOIN._runtime(base.base).wpe_cfg); mom=M.State(elapsed=M.B.rn32(3),samples=0)
        s=X.begin(clock,WM.State(cfg,mom,mom,mt.base.wpe))
        with self.assertRaises((ValueError,TypeError)): X.complete(s)

    def test_readiness_attaches_all_600_edges_but_not_source_uniform_arithmetic(self):
        r=X.readiness()
        self.assertTrue(r['Live_600_step_WPE_machine_history_attached'])
        self.assertTrue(r['complete_word_requires_full_WPE_machine_step_on_all_600_IMU_edges'])
        self.assertFalse(r['machine_vs_exact_WPE_period_branch_robustness_closed'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])

if __name__=='__main__': unittest.main()
