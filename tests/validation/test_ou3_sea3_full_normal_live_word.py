from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT/'tools'))
import ou3_sea3_full_normal_live_word_reset as mod
from ou3_interval import Interval,matrix_identity,matrix_point

class Sea3LiteralFullNormalLiveWordTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.d=mod.build(); cls.failures=mod.validate(cls.d)
    def test_shipping_order_and_reset_parity(self):
        self.assertEqual(self.failures,[]); self.assertTrue(all(self.d['shipping_event_order_parity'].values()))
        self.assertTrue(self.d['shipping_reset_source_parity_pass']); self.assertTrue(self.d['literal_reset_execution_complete'])
        self.assertTrue(self.d['reset_injection_supplied_by_same_source_word']); self.assertFalse(self.d['reset_small_angle_bound_required'])
    def test_literal_word_keeps_every_event_family(self):
        self.assertGreaterEqual(self.d['imu_samples_upper'],600); self.assertGreaterEqual(self.d['guaranteed_S_updates_lower_over_word'],4)
        self.assertTrue(self.d['every_valid_imu_sample_requires_prediction']); self.assertTrue(self.d['every_valid_imu_sample_requires_accelerometer_Joseph'])
        self.assertTrue(self.d['S_scheduler_is_executed_not_replaced_by_selected_four']); self.assertTrue(self.d['magnetometer_is_asynchronous_external_event_family'])
    def test_accelerometer_full_jacobian(self):
        for mode in ('H','A'):
            H=mod.H_accelerometer(mode,[Interval.point(1),Interval.point(2),Interval.point(3)],matrix_identity(3)); self.assertEqual(len(H[0]),mod.state_dimension(mode))
            for i in range(3): self.assertEqual(H[i][mod.OFF_AW+i],Interval.point(1))
            if mode=='A':
                for i in range(3): self.assertEqual(H[i][mod.OFF_BA+i],Interval.point(1))
    def test_one_literal_sample_executes_resets_between_measurements(self):
        n=mod.state_dimension('H'); P0=matrix_point([[2.0 if i==j else 0.0 for j in range(n)] for i in range(n)]); w=mod.initialize_word('H',P0)
        Faa=matrix_identity(6); Qaa=matrix_point([[0.01 if i==j else 0 for j in range(6)] for i in range(6)]); Fll=matrix_identity(12); Qll=matrix_point([[0.01 if i==j else 0 for j in range(12)] for i in range(12)]); F,Q=mod.pack_prediction('H',Faa,Qaa,Fll,Qll)
        d=[Interval.point(0.01),Interval.point(-0.005),Interval.point(0.002)]
        mod.apply_imu_sample(w,F=F,Q=Q,f_cog_body=[Interval.point(0),Interval.point(0),Interval.point(-9.80665)],R_wb=matrix_identity(3),Racc=mod.diagonal_R([0.2]*3),due_S=True,rs_std_xyz=[Interval.point(0.72),Interval.point(0.72),Interval.point(1)],Delta_aw=matrix_point([[0.001 if i==j else 0 for j in range(3)] for i in range(3)]),S_reset_dtheta=d,acc_reset_dtheta=d)
        self.assertEqual(w.event_log,['prediction','aw_floor','S_zero','left_reset_S','accelerometer','left_reset_acc']); self.assertTrue(mod.BACKEND.decomposition_identity_enclosed(w.riccati))
    def test_reset_self_tests_cover_magnetometer_too(self):
        for mode in ('H','A'):
            s=self.d[f'{mode}_reset_execution_self_test']; self.assertTrue(s['decomposition_identity_enclosed']); self.assertTrue(s['all_three_measurements_followed_immediately_by_reset'])
    def test_reduced_routes_still_cannot_promote(self):
        r=self.d['no_reduced_promotion_routes']; self.assertTrue(all(v is False for v in r.values())); self.assertFalse(self.d['P3_CANONICAL_PASS'])

if __name__=='__main__': unittest.main()
