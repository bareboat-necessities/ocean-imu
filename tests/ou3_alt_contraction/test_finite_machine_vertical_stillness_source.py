"""Machine Mahony -> tracker LPF -> sigma-stillness source regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_machine_vertical_stillness_source as X
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_binary32_mahony as MAH
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as GUARD
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as SCFG
from tools.stability.ou3_alt_contraction import finite_stillness_sigma_binary32 as STILL
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as SIG
import test_finite_core as FC


def expw(x):
    lo,hi=SIG.exp_minus_enclosure(F(x)); return B.rn32((lo+hi)/2)


def quiet_guarded():
    g=B.rn32(F(196133,20000))
    phys=replace(FC.root('A').reference,acceleration=(0,0,0),beta=(0,0,0),gyro_bias=(0,0,0),q_world_to_body=(1,0,0,0))
    raw=SENSOR.RawImuSample(phys,(0,0,0),(0,0,0),(0,0,0),(0,0,0),(0,0,-g),(0,0,g))
    return SENSOR.guarded_sample(raw,GUARD.State(),GUARD.Config(),dt=B.rn32(F(1,200)))


def vertical_cfg():
    return V.Config(B.rn32(F(1,5)),B.rn32(F(1,50)),B.rn32(F(196133,20000)),B.rn32(20))


def still_cfg():
    return SCFG.Config(B.rn32(F(196133,20000)),B.rn32(F(1,20)),B.rn32(F(8,10000)),B.rn32(2),B.rn32(1),B.rn32(F(1,5)))


# Frozen dataclass; one module-level instance keeps it out of the argument
# defaults (B008) without changing the value any caller sees.
STILL_INIT=STILL.State()


def still_operands(vertical_lp,dt,cfg,state=STILL_INIT):
    g=B.rn32(cfg.gravity); a=B.rn32(cfg.energy_alpha); threshold=B.rn32(cfg.energy_thresh)
    an=B.div(vertical_lp,g); inst=B.mul(an,an); decay=B.sub(B.rn32(1),a)
    vals=STILL._sum_products(decay,state.energy,a,inst); en=vals[0]
    if en < threshold:
        st=min(B.add(state.still_time,dt),B.rn32(60)); return en,expw(st)
    return en,None


class Tests(unittest.TestCase):
    def test_first_tracker_LPF_sample_is_exact_seed_and_consumes_no_exp(self):
        s=X.LPFState(cutoff_hz=B.rn32(6)); x=B.rn32(F(3,7)); h=B.rn32(F(1,200))
        out=X.lpf_step(s,x=x,dt=h)
        self.assertEqual(out.state.value,x); self.assertTrue(out.state.initialized)
        self.assertIsNone(out.alpha); self.assertEqual(out.successor_values,(x,))
        with self.assertRaisesRegex(ValueError,'consumes no exp'):
            X.lpf_step(s,x=x,dt=h,alpha_exp=B.rn32(F(4,5)))

    def test_initialized_LPF_exp_and_contraction_are_same_expression_owned(self):
        s=X.LPFState(B.rn32(F(1,4)),B.rn32(6),True,1); x=B.rn32(F(1,2)); h=B.rn32(F(1,200))
        mag=B.mul(B.mul(B.mul(X.TWO,X.PI_F),s.cutoff_hz),h); alpha=expw(mag)
        om=B.sub(X.ONE,alpha); ax=B.mul(om,x); ap=B.mul(alpha,s.value)
        vals=tuple(sorted(set((B.add(ax,ap),B.fma(om,x,ap),B.fma(alpha,s.value,ax)))))
        out=X.lpf_step(s,x=x,dt=h,alpha_exp=alpha,successor=vals[-1])
        self.assertEqual(out.exp_argument,mag); self.assertEqual(out.successor_values,vals)
        with self.assertRaisesRegex(ValueError,'SAME rounded argument'):
            X.lpf_step(s,x=x,dt=h,alpha_exp=B.rn32(F(1,2)),successor=vals[-1])

    def test_complete_machine_source_shares_one_Mahony_vertical_between_band_and_LPF(self):
        guarded=quiet_guarded(); h=B.rn32(F(1,200)); vc=vertical_cfg(); sc=still_cfg()
        base=X.begin(V.State(initialized=True),0,False,STILL.State(),cutoff_hz=B.rn32(6))
        mah=MAH.step_initialized(base.vertical,vc,dt=h,gyro=guarded.raw_gyro_body,acc=guarded.conditioned_accel_body)
        lp=mah.vertical.vertical_accel
        en,atten=still_operands(lp,h,sc)
        out=X.step(base,guarded,vc,sc,dt=h,still_energy_successor=en,still_attenuation_exp=atten)
        self.assertEqual(out.band_input,out.mahony.vertical.vertical_accel)
        self.assertEqual(out.lpf.input,out.band_input)
        self.assertEqual(out.stillness.vertical_lp,out.lpf.state.value)
        self.assertEqual(out.state.samples,1)

    def test_detached_stillness_input_and_band_input_are_rejected_by_bindings(self):
        # Public binders fail closed on anything that is not the exact qualified source/result type.
        class Env: pass
        with self.assertRaises(TypeError): X.require_frontend_input(Env(),Env())
        with self.assertRaises(TypeError): X.require_sigma_stillness(Env(),Env())

    def test_source_shape_and_promotion_boundary(self):
        r=X.readiness()
        for k in ('shipping_vertical_LPF_stillness_source_shape_matches','initialized_private_Mahony_binary32_graph_consumed',
                  'same_machine_Mahony_vertical_drives_band_and_tracker_LPF','FreqInputLPF_binary32_recurrence_materialized',
                  'FreqInputLPF_exp_bound_to_same_rounded_argument_error_profile','same_stored_LPF_output_drives_binary32_stillness_projection',
                  'sigma_stillness_operands_bindable_to_same_machine_successor','machine_band_input_bindable_to_same_Mahony_successor'):
            self.assertTrue(r[k])
        for k in ('target_libm_and_compiler_profile_correspondence_closed','guard_binary32_machine_history_attached',
                  'startup_machine_history_attached','Live_600_step_machine_history_attached','source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
