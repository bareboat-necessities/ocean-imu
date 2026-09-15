from fractions import Fraction as F
import unittest
from tools.stability.ou3_alt_contraction import finite_wpe_moment_binary32 as X
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W


def q(x): return X.B.rn32(F(x))
CFG=X.Config.from_shadow(W.WPEConfig(1,4,F(1,20),1,180))
DT=q(F(1,200))


def expw(arg):
    lo,hi=X.EXP.exp_minus_enclosure(arg)
    return X.B.rn32((lo+hi)/2)


def base_step(state, x=0, **kw):
    leakarg=X.B.mul(CFG.lambda_,DT); d=expw(leakarg)
    return X.step(state,CFG,dt=DT,vertical_accel=q(x),decay_exp=d,
                  velocity_successor=q(0),elevation_successor=q(0),**kw)


class Tests(unittest.TestCase):
    def test_pre_moment_start_is_literal_and_consumes_no_moment_witness(self):
        out=base_step(X.State())
        self.assertEqual(out.branch,'pre-moment-start')
        self.assertEqual(out.state.samples,1)
        with self.assertRaisesRegex(ValueError,'pre-moment-start'):
            base_step(X.State(),canonical_period=q(6))

    def test_moment_successors_are_bound_to_same_machine_recurrence(self):
        s=X.State(elapsed=q(3))
        h=q(4*6); mdarg=X.B.div(DT,h); md=expw(mdarg); a=X.B.sub(X.ONE,md)
        ms=X.MomentSuccessors(a,q(0),q(0),q(0),q(0))
        out=base_step(s,moment_decay_exp=md,moment_successors=ms)
        self.assertEqual(out.branch,'insufficient-weight')
        bad=X.MomentSuccessors(q(F(1,2)),q(0),q(0),q(0),q(0))
        with self.assertRaisesRegex(ValueError,'weight successor'):
            base_step(s,moment_decay_exp=md,moment_successors=bad)

    def test_valid_raw_period_is_derived_not_supplied(self):
        s=X.State(elapsed=q(4),weight=q(1),velocity_sq=q(9),elevation_sq=q(1))
        horizon=q(24); mdarg=X.B.div(DT,horizon); md=expw(mdarg); a=X.B.sub(X.ONE,md)
        ms=X.MomentSuccessors(q(1),q(0),X._ema(q(9),q(0),a)[0],q(0),X._ema(q(1),q(0),a)[0])
        vm=q(0); em=q(0); vsecond=X.B.div(ms.velocity_sq,ms.weight); esecond=X.B.div(ms.elevation_sq,ms.weight)
        vvar=max(X.ZERO,X.B.sub(vsecond,X.B.mul(vm,vm))); evar=max(X.ZERO,X.B.sub(esecond,X.B.mul(em,em)))
        omega=X.B.sub(X.B.div(vvar,evar),X.B.mul(CFG.lambda_,CFG.lambda_))
        root=X.SQRT.sqrt32(omega)
        out=base_step(s,moment_decay_exp=md,moment_successors=ms,velocity_var_successor=vvar,elevation_var_successor=evar,sqrt_omega=root)
        self.assertEqual(out.branch,'valid-period')
        self.assertEqual(out.raw_period,X.B.div(X.TWO_PI,root))
        with self.assertRaises(TypeError):
            X.step(s,CFG,dt=DT,vertical_accel=q(0),decay_exp=expw(X.B.mul(CFG.lambda_,DT)),
                   velocity_successor=q(0),elevation_successor=q(0),moment_decay_exp=md,moment_successors=ms,
                   velocity_var_successor=vvar,elevation_var_successor=evar)

    def test_readiness_closes_raw_topology_not_deployment_libm_or_storage(self):
        r=X.readiness()
        self.assertTrue(r['raw_period_has_no_independent_input_port'])
        self.assertTrue(r['raw_period_derived_from_same_machine_moments_and_sqrt_result'])
        self.assertFalse(r['target_exp_sqrt_libm_correspondence_closed'])
        self.assertFalse(r['dual_compiler_persistent_moment_history_attached'])
        self.assertFalse(r['storage_search_allowed'])

if __name__=='__main__': unittest.main()
