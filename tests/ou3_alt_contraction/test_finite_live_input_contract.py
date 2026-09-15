"""Input admission precedes execution; scalar-prefix bounds retain startup memory."""
from dataclasses import replace
from fractions import Fraction as F
import unittest
from unittest.mock import patch

from tools.stability.ou3_alt_contraction import finite_live_input_contract as X
from tools.stability.ou3_alt_contraction import finite_live_disturbance_overflow as OLD
from tools.stability.ou3_alt_contraction import finite_admitted_source_imu_word as IMU
from tools.stability.ou3_alt_contraction import finite_admitted_runtime_word as RUNTIME
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V


def packet(gyro=(0,0,0),accel=(0,0,-10)):
    raw=OLD.quiet_packet()
    gyro=tuple(map(F,gyro)); accel=tuple(map(F,accel))
    residual=tuple(accel[i]+(OLD.GRAVITY if i==2 else 0) for i in range(3))
    return replace(raw,gyro_residual_internal=gyro,raw_gyro_body=gyro,
                   accel_residual_internal=residual,raw_accel_body=accel)


class Tests(unittest.TestCase):
    def test_hardware_scale_and_literal_conversion_have_positive_margins(self):
        report=X.build()
        self.assertEqual(X.validate(report),[])
        self.assertEqual(report['raw_gyro_component_upper'],35)
        self.assertEqual(report['raw_accel_component_upper'],160)
        self.assertGreater(report['conversion']['gyro_conversion_margin_rad_s'],F(9,100))
        self.assertGreater(report['conversion']['accel_conversion_margin_mps2'],3)
        self.assertTrue(report['startup_residual_and_temporal_contracts_unchanged'])
        self.assertFalse(report['finite_raw_input_cap_is_full_deployment_certificate'])

    def test_component_endpoints_and_api_rounding_remain_source_bound(self):
        raw=packet(gyro=(35,-35,F(1,10)),accel=(-160,160,F(-1,3)))
        admitted=X.check_packet(raw)
        self.assertIs(admitted.raw,raw)
        self.assertEqual(admitted.machine_gyro,(F(35),F(-35),X.B.rn32(F(1,10))))
        self.assertEqual(admitted.machine_accel,tuple(X.B.rn32(v) for v in raw.raw_accel_body))
        with self.assertRaisesRegex(ValueError,'detached from same source API rounding'):
            replace(admitted,machine_gyro=(F(35),F(-35),F(0)))
        # The domain is on actual API values. Preserve the complete rounding
        # cell instead of pruning a legal real-shadow half-ULP sliver.
        rounded=X.check_packet(packet(
            gyro=(X.source_cell_upper(35),0,0),
            accel=(X.source_cell_upper(160),0,0)))
        self.assertEqual(rounded.machine_gyro[0],35)
        self.assertEqual(rounded.machine_accel[0],160)
        for sign in (-1,1):
            with self.assertRaisesRegex(ValueError,'gyro exceeds'):
                X.check_packet(packet(gyro=(sign*(35+F(1,1<<18)),0,0)))
            with self.assertRaisesRegex(ValueError,'accelerometer exceeds'):
                X.check_packet(packet(accel=(sign*(160+F(1,1<<16)),0,0)))

    def test_old_overflow_packet_is_rejected_before_either_theorem_event(self):
        raw=OLD.quiet_packet(pulse=True)
        with self.assertRaisesRegex(ValueError,'gyro exceeds'):
            X.check_packet(raw)
        with patch.object(IMU.PRED,'imu_step') as execute:
            with self.assertRaisesRegex(ValueError,'gyro exceeds'):
                IMU.imu_step(None,restricted=None,bias_restricted=None,witness=None,
                             raw=raw,packet_id='pulse',temperature_c=35)
            execute.assert_not_called()
        with patch.object(RUNTIME.IMU,'imu_step') as execute:
            with self.assertRaisesRegex(ValueError,'gyro exceeds'):
                RUNTIME.imu_step(None,restricted=None,bias_restricted=None,witness=None,
                                 raw=raw,packet_id='pulse')
            execute.assert_not_called()

    def test_actual_scalar_extremes_stay_inside_proved_prefix_bounds(self):
        certificate=X.prefix_certificate()
        cfg=V.Config(X.B.rn32(F(1,5)),X.B.rn32(F(1,50)),X.B.rn32(X.START.G),20)
        # Operation correspondence at the allowed extremes; it is not a
        # hand-installed source trajectory or evidence for capture/contraction.
        for acc in ((160,-160,160),(0,0,0),(F(1,1<<149),0,0)):
            state=V.State(q=(F(1,2),)*4,integral=(F(3),F(-3),F(3)),initialized=True)
            for gyro in ((35,-35,35),(-35,35,-35)):
                out=X.MAH.step_initialized(state,cfg,dt=X.CLOCK.DT_FLOAT,
                    gyro=gyro,acc=acc)
                for operation in out.operations: operation.validate()
                self.assertLess(sum(q*q for q in out.vertical.state.q),X.PREFIX.Q2)
                self.assertLess(max(abs(i) for i in out.vertical.state.integral),certificate['integral_component_upper'])
                self.assertLess(abs(out.vertical.state.up),certificate['vertical_abs_upper'])

    def test_whole_prefix_scalar_induction_does_not_reset_or_promote_MEKF(self):
        c=X.prefix_certificate()
        self.assertIsNone(c['max_steps'])
        self.assertTrue(c['startup_integral_not_reset_at_Live'])
        self.assertLess(c['accel_norm_sum_upper'],100000)
        self.assertEqual(c['integral_component_upper'],4096)
        self.assertLess(c['feedback_corrected_rate_component_upper'],4132)
        self.assertLess(c['Euler_norm_sum_upper'],5000)
        self.assertLess(c['vertical_abs_upper'],322)
        self.assertTrue(c['initialized_scalar_Mahony_prefix_totality_closed'])
        self.assertTrue(c['all_initialized_finite_prefixes_totality_closed'])
        self.assertEqual(c['finite_horizon_refinement']['max_steps'],30602)
        self.assertLess(c['finite_horizon_refinement']['integral_component_upper'],F(41,10))
        self.assertFalse(c['startup_seed_and_reachability_inferred'])
        self.assertFalse(c['target_compiler_and_libm_qualified'])
        self.assertFalse(c['complete_MEKF_covariance_solver_word_totality_closed'])

    def test_integral_and_elapsed_lattice_barriers_are_actual_RNE_fixed_points(self):
        c=X.prefix_certificate()
        barrier=c['integral_component_upper']; delta=c['integral_increment_upper']
        self.assertGreater(c['integral_barrier_rounding_margin'],F(1,10000))
        for sign in (-1,1):
            self.assertEqual(X.B.add(sign*barrier,sign*delta),sign*barrier)
        self.assertEqual(X.B.add(c['elapsed_float_upper'],X.CLOCK.DT_FLOAT),c['elapsed_float_upper'])
        cfg=V.Config(X.B.rn32(F(1,5)),X.B.rn32(F(1,50)),X.B.rn32(X.START.G),20)
        state=V.State(q=(F(1,2),)*4,integral=(barrier,-barrier,barrier),
                      initialized=True,elapsed=c['elapsed_float_upper'])
        out=X.MAH.step_initialized(state,cfg,dt=X.CLOCK.DT_FLOAT,
                                  gyro=(35,-35,35),acc=(-160,160,160))
        for operation in out.operations: operation.validate()
        self.assertLessEqual(max(abs(i) for i in out.vertical.state.integral),barrier)
        self.assertEqual(out.vertical.state.elapsed,c['elapsed_float_upper'])
        self.assertLess(sum(q*q for q in out.vertical.state.q),X.PREFIX.Q2)
        self.assertFalse(c['indefinite_clock_progress_inferred'])

    def test_default_app_temperature_zero_does_not_restrict_general_API(self):
        c=X.default_application_temperature_certificate()
        self.assertEqual(c['application_default_temperature_c'],35)
        self.assertEqual(c['application_default_thermal_model_component'],0)
        self.assertTrue(c['shipping_three_argument_call_thermal_term_zero'])
        self.assertFalse(c['general_four_argument_API_temperature_restricted'])
        self.assertFalse(c['general_thermal_history_totality_qualified'])


if __name__=='__main__': unittest.main()
