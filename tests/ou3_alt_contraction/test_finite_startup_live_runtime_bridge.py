"""TunerReady/goLive -> first Live-prefix state regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_startup_live_runtime_bridge as X
from tools.stability.ou3_alt_contraction import finite_startup_handoff_seed as SEED
from tools.stability.ou3_alt_contraction import finite_startup_handoff_init as INIT
from tools.stability.ou3_alt_contraction import finite_startup_live_entry as LIVE
from tools.stability.ou3_alt_contraction import finite_fresh_joint24_entry as FRESH
from tools.stability.ou3_alt_contraction import finite_tuner_commit as COMMIT
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_prefix as FRONT
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as GUARD
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_post_prediction as POST
from tools.stability.ou3_alt_contraction import finite_periodic_aw_sync as AWSYNC
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE
import test_finite_tuner_frontend_prefix as TP


def cfg():
    return COMMIT.CommitConfig(F(3,2),F(1,100),2,F(3,200),False,False,
                               F(2,5),5,F(4,5),F(6,5),F(7,5))


def objects(pending=True):
    tuner=TP.state(stage='TunerReady',stage_time=F(7,10),pending=pending,time=0)
    front=FRONT.State(GUARD.State(),tuner)
    c=COMMIT.commit(tuner.tune,cfg(),pending=True,live=True,
                    band_noise_floor_sigma=0,sync_covariance=True)
    active=ACTIVE.ActiveParameters.from_commit(c)
    seed=SEED.Result((1,0,0,0),SEED.TILT_SIGMA,SEED.YAW_SIGMA_GAUGED,False,True)
    x=tuple(F(0) for _ in range(21)); P=[[F(1 if i==j else 0) for j in range(21)] for i in range(21)]
    hand=INIT.initialize_from_gauged_seed_zero_heel(seed,x,P,scope=SCOPE.certified_scope())
    entry=LIVE.enter_live(hand,active,scope=SCOPE.certified_scope())
    ref=CORE.Reference(0,(1,0,0,0),(0,0,0),(0,0,0),(0,0,0),(0,0,0),
                       (0,0,0),(0,0,0),0,'bridge-history','bridge-bias','BIAS0')
    fresh=FRESH.build(entry,ref,scope=SCOPE.certified_scope())
    scheduler=POST.Scheduler(active.pseudo_period,0,0)
    return front,entry,fresh,active,scheduler


class Tests(unittest.TestCase):
    def test_goLive_preserves_frontend_memory_and_enters_first_live_state(self):
        front,entry,fresh,active,scheduler=objects(True)
        out=X.bridge(entry,fresh,front,scope=SCOPE.certified_scope(),commit_cfg=cfg(),
                     bench_noise_sigma=0,noise_sqrt=BAND.NoiseSqrtWitness(0),
                     scheduler=scheduler,racc=RACC.State(),aw_sync=AWSYNC.State())
        self.assertEqual(out.state.mekf,fresh); self.assertEqual(out.state.active,active)
        self.assertEqual(out.frontend_live.guard,front.guard)
        self.assertEqual(out.frontend_live.tuner.stage,'Live'); self.assertEqual(out.frontend_live.tuner.stage_time,0)
        self.assertTrue(out.frontend_live.tuner.pending)  # goLive does not consume online pending
        self.assertEqual(out.frontend_live.tuner.wpe,front.tuner.wpe)
        self.assertEqual(out.frontend_live.tuner.band,front.tuner.band)
        self.assertEqual(out.frontend_live.tuner.vertical,front.tuner.vertical)

    def test_goLive_commit_must_equal_fresh_entry_active_parameters(self):
        front,entry,fresh,active,scheduler=objects(False)
        bad=ACTIVE.ActiveParameters(active.tau+1,active.Sigma_aw,active.pseudo_period,active.R_S)
        from dataclasses import replace
        detached=replace(entry,active=bad)
        with self.assertRaisesRegex(ValueError,'detached from carried TunerReady'):
            X.bridge(detached,fresh,front,scope=SCOPE.certified_scope(),commit_cfg=cfg(),
                     bench_noise_sigma=0,noise_sqrt=BAND.NoiseSqrtWitness(0),
                     scheduler=scheduler,racc=RACC.State(),aw_sync=AWSYNC.State())

    def test_startup_cannot_supply_periodic_aw_floor_pending(self):
        front,entry,fresh,active,scheduler=objects(False)
        aw=AWSYNC.State(True,0,active.Sigma_aw)
        with self.assertRaisesRegex(ValueError,'startup cannot enter Live'):
            X.bridge(entry,fresh,front,scope=SCOPE.certified_scope(),commit_cfg=cfg(),
                     bench_noise_sigma=0,noise_sqrt=BAND.NoiseSqrtWitness(0),
                     scheduler=scheduler,racc=RACC.State(),aw_sync=aw)

    def test_readiness_keeps_first_sample_and_precision_open(self):
        r=X.readiness()
        self.assertTrue(r['fresh_H18_CORE_state_connected_to_first_Live_prefix_shape'])
        self.assertTrue(r['online_pending_bit_preserved_across_goLive'])
        self.assertFalse(r['first_Live_sample_executed_from_fresh_entry'])
        self.assertFalse(r['ALT_STARTUP_PASS']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
