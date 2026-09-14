"""Bounded IMU ISS-history regressions; not contraction evidence."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_imu_disturbance as D
from tools.stability.ou3_alt_contraction import finite_source_bound_imu_forcing as FORCE
import test_finite_source_bound_live_word as BASE


class Tests(unittest.TestCase):
    def test_symbolic_bound_accepts_exact_executed_supply_coordinates(self):
        s=BASE.root_state()
        h=D.BoundedHistory(s.sensor_root,'temperature-history',F(1000))
        f=FORCE.ImuForcing((1,0,0),(2,0,0),(3,0,0),(-1,0,0),0)
        r=D.bind(h,sensor_root=s.sensor_root,ordinal=1,forcing=f)
        self.assertEqual(r.forcing.supply_vector,(F(1),0,0,-1,0,0,3,0,0))
        self.assertEqual(r.history.supply_norm_upper,1000)

    def test_bound_is_not_Racc_or_a_fixed_sigma_multiple(self):
        s=BASE.root_state()
        h=D.BoundedHistory(s.sensor_root,'temperature-history',F(1,10))
        f=FORCE.ImuForcing((1,0,0),(0,0,0),(0,0,0),(0,0,0),0)
        with self.assertRaisesRegex(ValueError,'exceeds carried theorem bound'):
            D.bind(h,sensor_root=s.sensor_root,ordinal=1,forcing=f)
        self.assertTrue(D.readiness()['Racc_covariance_not_used_as_pathwise_noise_bound'])
        self.assertTrue(D.readiness()['symbolic_supply_bound_not_fixed_sigma_multiple'])

    def test_sensor_root_cannot_restart(self):
        s=BASE.root_state()
        h=D.BoundedHistory(s.sensor_root,'temperature-history',10)
        other=replace(s.sensor_root,gyro_residual_history_id='other-gyro-history')
        f=FORCE.ImuForcing((0,0,0),(0,0,0),(0,0,0),(0,0,0),0)
        with self.assertRaisesRegex(ValueError,'detached from carried sensor residual histories'):
            D.bind(h,sensor_root=other,ordinal=1,forcing=f)

    def test_ordinal_and_bound_are_finite_word_theorem_data(self):
        s=BASE.root_state(); h=D.BoundedHistory(s.sensor_root,'temperature-history',0)
        f=FORCE.ImuForcing((0,0,0),(0,0,0),(0,0,0),(0,0,0),0)
        D.bind(h,sensor_root=s.sensor_root,ordinal=600,forcing=f)
        with self.assertRaisesRegex(ValueError,'one of the 600'):
            D.bind(h,sensor_root=s.sensor_root,ordinal=601,forcing=f)

    def test_readiness_keeps_storage_closed(self):
        r=D.readiness()
        self.assertTrue(r['arbitrary_bounded_IMU_ISS_history_quantifier_available'])
        self.assertFalse(r['finite_horizon_probability_used_to_prune_disturbances'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
