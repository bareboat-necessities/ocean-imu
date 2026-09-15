from fractions import Fraction as F
import unittest
from tools.stability.ou3_alt_contraction import finite_wpe_machine_binary32 as X
from tools.stability.ou3_alt_contraction import finite_wpe_moment_binary32 as M
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as L
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W


def q(x): return M.B.rn32(F(x))
CFG=W.WPEConfig(1,4,F(1,20),1,180); DT=q(F(1,200))

def expw(arg):
    lo,hi=M.EXP.exp_minus_enclosure(arg); return M.B.rn32((lo+hi)/2)

def pre():
    leak=expw(M.B.mul(q(1),DT))
    return dict(decay_exp=leak,velocity_successor=q(0),elevation_successor=q(0))

class Tests(unittest.TestCase):
    def test_reset_product_advances_same_vertical_input_on_both_modes(self):
        s=X.initial(CFG); w=X.ModeWitnesses(pre())
        n,sr,fr,logs=X.step(s,dt=DT,vertical_accel=q(0),separate=w,fma=w)
        self.assertEqual((sr.vertical_accel,fr.vertical_accel),(q(0),q(0)))
        self.assertEqual(n.logs.samples,1)
        self.assertEqual((sr.branch,fr.branch),('pre-moment-start','pre-moment-start'))

    def test_valid_raw_period_must_bind_to_same_log_input(self):
        s=X.initial(CFG)
        m=M.State(elapsed=q(4),weight=q(1),velocity_sq=q(9),elevation_sq=q(1),samples=10)
        logs=L.State(L.Track(),L.Track(),10); s=X.State(s.cfg,m,m,logs)
        horizon=q(24); md=expw(M.B.div(DT,horizon)); a=M.B.sub(M.ONE,md)
        ms=M.MomentSuccessors(q(1),q(0),M._ema(q(9),q(0),a)[0],q(0),M._ema(q(1),q(0),a)[0])
        vvar=M.B.div(ms.velocity_sq,ms.weight); evar=M.B.div(ms.elevation_sq,ms.weight)
        omega=M.B.sub(M.B.div(vvar,evar),M.B.mul(s.cfg.lambda_,s.cfg.lambda_)); root=M.SQRT.sqrt32(omega)
        moment=dict(pre(),moment_decay_exp=md,moment_successors=ms,velocity_var_successor=vvar,elevation_var_successor=evar,sqrt_omega=root)
        raw=M.B.div(M.TWO_PI,root); lr=q(2)
        w=X.ModeWitnesses(moment,raw_log=X.RawLogBinding(raw,lr),log=L.InitWitness(lr),
            usable=X.USABLE.PeriodWitness(lr,q(8)))
        n,sr,fr,_=X.step(s,dt=DT,vertical_accel=q(0),separate=w,fma=w)
        self.assertEqual((sr.raw_period,fr.raw_period),(raw,raw))
        self.assertEqual((n.logs.separate.log_period,n.logs.fma.log_period),(lr,lr))
        self.assertFalse(n.separate_usable)
        self.assertFalse(n.fma_usable)
        from dataclasses import replace
        different=replace(w,usable=X.USABLE.PeriodWitness(lr,q(1)))
        mixed,_,_,_=X.step(s,dt=DT,vertical_accel=q(0),separate=different,fma=w)
        self.assertTrue(mixed.separate_usable)
        self.assertFalse(mixed.fma_usable)
        bad=X.ModeWitnesses(moment,raw_log=X.RawLogBinding(q(3),lr),log=L.InitWitness(lr))
        with self.assertRaisesRegex(ValueError,'raw period'):
            X.step(s,dt=DT,vertical_accel=q(0),separate=bad,fma=w)

    def test_readiness_attaches_persistent_raw_period_to_log_but_not_libm_or_storage(self):
        r=X.readiness()
        self.assertTrue(r['dual_compiler_persistent_moment_history_attached'])
        self.assertTrue(r['WPE_raw_period_binary32_production_attached_to_log_history'])
        self.assertFalse(r['target_exp_log_sqrt_libm_correspondence_closed'])
        self.assertFalse(r['startup_frontend_vertical_ancestry_attached'])
        self.assertFalse(r['storage_search_allowed'])

if __name__=='__main__': unittest.main()
