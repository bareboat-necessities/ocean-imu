"""Literal construction/ungauged handoff; not a startup capture certificate."""
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_machine_startup_core as X
from tools.stability.ou3_alt_contraction import finite_machine_core_continuation as CARRY
import test_finite_startup_joined_machine_history as START
import test_finite_ungauged_live_word as UNGAUGED


def produced_startup():
    state=START.initial(); raw,kwargs=START.cold_operands(state)
    return START.X.step(state,raw,**kwargs).state


class Tests(unittest.TestCase):
    def test_literal_construction_is_sealed_and_preserved_by_startup(self):
        state=START.initial(); root=state.core_construction
        with self.assertRaisesRegex(TypeError,'only by literal_reset'): X.Construction()
        with self.assertRaises(FrozenInstanceError): root.Pb0=F(1)
        self.assertEqual(root.Pb0,X.B.rn32(F(1,10**6)))
        self.assertEqual(root.ba_variance,X.B.mul(X.B.rn32(F(1,250)),X.B.rn32(F(1,250))))
        raw,kw=START.cold_operands(state)
        state=START.X.step(state,raw,**kw).state
        self.assertIs(state.core_construction,root)
        state,_=START.X.boundary(state)
        self.assertIs(state.core_construction,root)

    def test_actual_private_proxy_produces_complete_default_ungauged_CORE(self):
        state=produced_startup()
        reference=replace(state.last_raw.physical,time=state.base.lower.frontend.tuner.time,
                          live_origin=state.base.lower.frontend.tuner.time,centered_S=(0,0,0))
        go,_,_,_,_=START.admitted_fixture()
        for mode in ('separate','fma'):
            active=getattr(go.lower,mode+'_active')
            result=X.ungauged(state,reference,active,mode=mode)
            self.assertIs(result.source,getattr(state,mode+'_source'))
            self.assertIs(result.construction,state.core_construction)
            self.assertEqual(result.state.q_hat,(1,0,0,0))
            self.assertEqual(result.state.mode,'H')
            self.assertEqual(result.state.z[3:6],reference.gyro_bias)
            self.assertEqual(result.state.z[6:9],reference.velocity)
            self.assertEqual(result.state.z[9:12],reference.position)
            self.assertEqual(result.state.z[15:18],reference.acceleration)
            self.assertEqual(result.state.z[18:21],reference.beta)
            for k in range(3):
                self.assertEqual(result.state.covariance[3+k][3+k],state.core_construction.Pb0)
                self.assertEqual(result.state.covariance[6+k][6+k],1)
                self.assertEqual(result.state.covariance[9+k][9+k],400)
                self.assertEqual(result.state.covariance[12+k][12+k],2500)
                self.assertEqual(result.state.covariance[18+k][18+k],state.core_construction.ba_variance)
            self.assertEqual(tuple(row[15:18] for row in result.state.covariance[15:18]),active.Sigma_aw)
            self.assertFalse(result.covariance_fallback)
            self.assertTrue(all(X.B.is_binary32(op.result) for op in result.operations))
            with self.assertRaisesRegex(TypeError,'executed from startup'):
                X.Result(result.construction,result.source,mode,result.state,(),False)
        with self.assertRaises(TypeError):
            X.ungauged(state,reference,go.lower.separate_active,mode='separate',q_hat=(1,0,0,0))
        with self.assertRaises(TypeError):
            X.ungauged(state,reference,go.lower.separate_active,mode='separate',covariance=((1,),))

    def test_covariance_fallback_is_literal_and_has_no_final_axis_input(self):
        # Exercise the finite covariance helper's fallback on a degenerate
        # nonunit input; this is not a reachable normalized startup claim.
        for mode in ('separate','fma'):
            a=X.Arithmetic(mode)
            cov,fallback=X._attitude_covariance(a,(0,F(1,2),F(1,2),0),
                X.B.rn32(F(35,1000)),X.B.rn32(F(15708,10000)))
            self.assertTrue(fallback)
            var=X.B.mul(X.B.rn32(F(15708,10000)),X.B.rn32(F(15708,10000)))
            self.assertEqual(cov,tuple(tuple(var if i==j else 0 for j in range(3)) for i in range(3)))

    def test_no_seed_for_uninitialized_or_other_physical_startup(self):
        state=START.initial(); runtime=UNGAUGED.root(south=False)
        active=runtime.live.live.active; reference=runtime.live.live.mekf.reference
        with self.assertRaisesRegex(ValueError,'actual initialized'):
            X.ungauged(state,reference,active,mode='separate')
        with self.assertRaisesRegex(ValueError,'physical/bias lineage'):
            X.ungauged(produced_startup(),reference,active,mode='separate')

    def test_readiness_keeps_machine_yaw_and_target_schedule_open(self):
        r=X.readiness()
        self.assertTrue(r['ungauged_machine_proxy_to_full_CORE_producer_available'])
        self.assertFalse(r['caller_replacement_quaternion_or_covariance_port'])
        self.assertFalse(r['gauged_machine_pending_yaw_producer_available'])
        self.assertFalse(r['target_Eigen_compiler_operation_schedule_qualified'])
        self.assertFalse(r['source_uniform_startup_CORE_qualified'])
        self.assertFalse(CARRY.readiness()['machine_startup_quaternion_and_remaining_covariance_producer_attached'])


if __name__=='__main__': unittest.main()
