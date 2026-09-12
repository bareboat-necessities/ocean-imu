"""Shared private-vertical frontend composition regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_frontend_runtime as X
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as B
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as S


def vertical(v=2): return V.Result(V.State(initialized=True,up=F(v)),F(v))
def bcfg(): return B.BandConfig(F(1,2),4,F(1,100),6,F(3,100),F(6,5),F(1,5))
def scfg(): return B.StatsConfig(4,F(3,10),60,F(1,20),5)


class Tests(unittest.TestCase):
    def test_first_tracker_lpf_sample_is_exact_shared_vertical_output(self):
        out=X.tracker_lpf_step(X.LPFState(),vertical(2),decay=X.LPFDecayWitness(F(1,2)))
        self.assertTrue(out.state.initialized); self.assertEqual(out.output,2)

    def test_initialized_tracker_lpf_uses_same_vertical_output(self):
        out=X.tracker_lpf_step(X.LPFState(4,True),vertical(2),decay=X.LPFDecayWitness(F(1,2)))
        self.assertEqual(out.output,3)

    def test_wpe_entry_owns_vertical_sample(self):
        out=X.wpe_step_from_vertical(W.WPEState(),W.WPEConfig(1,4,F(1,2),1,180),vertical(2),dt=F(1,10),decay=W.ExpWitness(F(1,2)))
        self.assertEqual(out.state.accel_prev,2)
        with self.assertRaises(TypeError):
            X.wpe_step_from_vertical(W.WPEState(),W.WPEConfig(1,4,F(1,2),1,180),vertical(2),dt=F(1,10),decay=W.ExpWitness(F(1,2)),vertical_accel=9)

    def test_sigma_band_entry_owns_same_vertical_sample(self):
        wp=W.UpdateResult(W.WPEState(log_period=1,usable_period=True),True,F(1,2),2,None,None)
        out=X.band_step_from_vertical(B.BandState(),B.StatsState(),wp,vertical(2),dt=F(1,10),band_cfg=bcfg(),stats_cfg=scfg(),band_decay=B.BandDecayWitness(F(1,2),F(1,2)),variance_decay=B.VarianceDecayWitness(F(1,2)),bench_noise_sigma=2,noise_sqrt=B.NoiseSqrtWitness(F(1,4)))
        self.assertEqual(out.filtered_accel,F(1,2))
        with self.assertRaises(TypeError):
            X.band_step_from_vertical(B.BandState(),B.StatsState(),wp,vertical(2),vertical_accel=9,dt=F(1,10),band_cfg=bcfg(),stats_cfg=scfg(),band_decay=B.BandDecayWitness(F(1,2),F(1,2)),variance_decay=B.VarianceDecayWitness(F(1,2)),bench_noise_sigma=2,noise_sqrt=B.NoiseSqrtWitness(F(1,4)))

    def test_stillness_owns_LPF_and_tracker_outputs(self):
        lpf=X.LPFResult(X.LPFState(0,True),0); tracker=X.TrackerOutputWitness(F(2,5))
        out=X.stillness_step_from_vertical(S.State(),S.Config(gravity=1),lpf,tracker,dt=1,attenuation=S.AttenuationWitness(F(1,2)))
        self.assertTrue(out.state.last_is_still); self.assertEqual(out.tracker_frequency_out,F(2,5))
        with self.assertRaises(TypeError):
            X.stillness_step_from_vertical(S.State(),S.Config(gravity=1),lpf,tracker,dt=1,a_vert_up_lp=9,attenuation=S.AttenuationWitness(F(1,2)))

    def test_readiness_fail_closed_at_tracker_and_sensor_ancestry(self):
        r=X.readiness(); self.assertTrue(r['one_private_vertical_output_shared_by_WPE_and_sigma_band']); self.assertTrue(r['same_private_vertical_output_drives_tracker_input_LPF']); self.assertFalse(r['tracker_output_algorithm_attached']); self.assertFalse(r['private_vertical_seed_and_fast_inv_sqrt_attached']); self.assertFalse(r['raw_sensor_BRMM_disturbance_relation_attached']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
