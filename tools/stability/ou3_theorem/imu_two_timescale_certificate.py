"""Exact scalar consequences of the two-timescale model; never hardware evidence."""
from fractions import Fraction as F
import json
from pathlib import Path


def certificate():
    c=json.loads(Path(__file__).with_name('constants.json').read_text(),parse_float=F)
    b=c['imu_bias']
    g,alpha,nu=F('9.80665'),F(1,100),F(1,2)
    da,dg=g*alpha*nu,alpha*nu**2
    # Rational enclosure of pi. Every reported lower cap is rounded DOWN by
    # choosing the upper half-period in the adverse terms, not a sampled max.
    pi_lo=F('3.14159265358979323846264338327950288')
    pi_hi=F('3.14159265358979323846264338327950289')
    lo,hi=pi_lo/nu,pi_hi/nu
    hmax=F(c['sensor_model']['sample_period_max_s'])
    ia=(g/nu)*(2*alpha-F(2,9)*alpha**3) # sin x >= x-x^3/6, positive half-wave
    ig=2*alpha
    def cap_lower(integral,slow_amp,slow_rate,signal_rate,sampled):
        quadrature=hi*signal_rate*hmax if sampled else F(0)
        ds=min(2*slow_amp,slow_rate*(hi+(hmax if sampled else 0)))
        return max(F(0),integral-quadrature-hi*ds/2)
    slow_a=cap_lower(ia,b['B_a_s_mps2'],b['D_a_s_mps3'],da,False)
    slow_g=cap_lower(ig,b['B_g_s_rad_s'],b['D_g_s_rad_s2'],dg,False)
    sampled_a=cap_lower(ia,b['B_a_s_mps2'],b['D_a_s_mps3'],da,True)
    sampled_g=cap_lower(ig,b['B_g_s_rad_s'],b['D_g_s_rad_s2'],dg,True)
    assert da>b['D_a_s_mps3'] and dg>b['D_g_s_rad_s2']
    assert 0<sampled_a<slow_a and 0<sampled_g<slow_g
    return {
        'qualification':'OU3_IMU_TWO_TIMESCALE_ALGEBRA_V1',
        'model':'e_a=b_a_s+b_a_f; e_g=b_g_s+b_g_f',
        'temporal_scope':'every placed window of calibrated delivered-sample hold; one carried decomposition',
        'fast_accel_numerical_qualification':'OPEN',
        'fast_gyro_numerical_qualification':'OPEN',
        'six_amplitude_rate_limits':'inherited declared candidates, not newly measured device guarantees',
        'raw_two_epoch_rule':'min(2Bs,Ds*h)+min(Bf,C/min(H,dt_i))+min(Bf,C/min(H,dt_j)); complete distinct cells, global class supremum',
        'raw_two_epoch_independent_fast_box_admission':False,
        'window_mean_difference_bound':'min(2Bs,Ds*h)+2 min(K*(L),K*(h))/L for continuous slow averages; add slow sampling defect for held samples',
        'matrix_weighted_fast_endpoint_term_retained':True,
        'pair_history_differences_require_two_reachable_error_histories':True,
        'oscillatory_witness':{
            'phi':'(1/100)*sin(t/2)',
            'accel_slow_required_rate_mps3':str(da),
            'gyro_slow_required_rate_rad_s2':str(dg),
            'accel_rate_to_candidate_limit_ratio':str(da/b['D_a_s_mps3']),
            'gyro_rate_to_candidate_limit_ratio':str(dg/b['D_g_s_rad_s2']),
            'all_slow_admissible_under_candidate_rates':False,
            'half_period_s':'2*pi',
            'half_period_rational_bracket':[str(lo),str(hi)],
            'continuous_necessary_accel_K_lower_mps':str(slow_a),
            'continuous_necessary_gyro_K_lower_rad':str(slow_g),
            'delivered_sample_necessary_accel_K_lower_mps':str(sampled_a),
            'delivered_sample_necessary_gyro_K_lower_rad':str(sampled_g),
            'bound_scope':'necessary for ANY slow/fast split, not fitted fast-model constants',
            'current_model_admissibility':'OPEN',
            'excluded_by_merely_renaming_noise_fast':False,
            'historical_V3_counterexample_preserved':True,
            'historical_BA_storage_lower_bound_transferred_to_new_split':False,
        },
        'marine_and_magnetic':{
            'assumptions_changed':False,
            'stronger_EXCITED_MOVING_added':False,
            'T_E_and_theta_E_numerically_qualified':False,
            'strict_general_gauge_breaking_margin_proved':False,
            'all_slow_sin_cubed_family_remains':'regime-certificate.json; phi=.001*sin(t/40)^3, T_E=80*pi, theta_E=.002, fast=0',
            'scope':'symbolic existing MARINE permits this family; not a counterexample for every fixed T_E/theta_E',
        },
        'source_uniform_finite_error_supply_closed':False,
        'runtime_calibration_tuner_and_gates_changed':False,
        'theorem_closed':False,
    }


if __name__=='__main__':
    print(json.dumps(certificate(),indent=2,sort_keys=True))
