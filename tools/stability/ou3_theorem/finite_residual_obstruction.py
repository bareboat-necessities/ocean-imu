"""HISTORICAL all-time obstruction under the V3 norm-only residual contract.

The witness is a physical truth producing the existing quiet packet history.
It is NOT inferred from homogeneous zero action. No runtime/noise setting is
modified, and no float32 all-time theorem is asserted.
"""
from __future__ import annotations
from fractions import Fraction as F
import json
from pathlib import Path
from .sampled_capture_obstruction import service_certificate

# Frozen old-domain limits: historical witness, NOT controlling IMU BIAS.
LEGACY_ACCEL_RESIDUAL = F(".3")
LEGACY_GYRO_RESIDUAL = F(".02")
LEGACY_ACCEL_BIAS = F(".22516660498395405")


def certificate():
    c = json.loads(Path(__file__).with_name('constants.json').read_text(), parse_float=F)
    g, g_model = F('9.80665'), F('9.8066501617431640625')
    alpha, nu, ba = F(1, 100), F(1, 2), F(1, 100)
    # R_bw=Rx(alpha*sin(nu*t)), B=75 ex, p=v=a=S=0.
    # ba=ba*ex; na=g*R_bw' ez-g_model*ez-ba*ex; ng=-phi'*ex.
    # Exact packets: f=-g_model*ez, gyro=0, mag=75*ex.
    accel_upper = g*alpha + ba + abs(g-g_model)
    gyro_upper = alpha*nu
    assert accel_upper < LEGACY_ACCEL_RESIDUAL
    assert gyro_upper < LEGACY_GYRO_RESIDUAL
    assert ba < LEGACY_ACCEL_BIAS
    assert gyro_upper < c['marine_motion']['Omega_max_rad_s']
    # Period=4*pi<88/7<60 and span=2*alpha>pi/180, using pi<22/7.
    assert 4*F(22, 7) < 60
    assert 2*alpha > F(22, 7)/180
    # Heading/axial-BG root has z projection cos(phi_root). Quiet axis
    # groups decouple, so retain its positive Gram and discard the other.
    quiet = service_certificate()
    service = F(quiet['actual_innovation_service_lower'])*(1-alpha**2/2)**2
    assert service > c['magnetic_service']['mu_M']
    r_lower = ba/F(1, 40)  # P_ba,ba <= (1/40)^2 I, actual bhat_a=0.
    assert r_lower > F(3, 20)
    return {
        'qualification': 'OU3_NORM_ONLY_FINITE_RESIDUAL_OBSTRUCTION_V1',
        'physical_roll': '(1/100)*sin(t/2)',
        'example_excited_window_s': 60,
        'example_excited_span_rad': 'pi/180',
        'all_time_gravity_span_rad': str(2*alpha),
        'physical_accel_bias_mps2': str(ba),
        'physical_bias_rates': 'zero',
        'translation_and_jerk_and_primitive': 'zero',
        'accel_residual_norm_upper_mps2': str(accel_upper),
        'gyro_residual_norm_upper_rad_s': str(gyro_upper),
        'actual_magnetic_service_lower': str(service),
        'packets': 'gyro=0; accel=-g_model*ez; mag=75*ex',
        'coupled_tuner_history': 'identical quiet history by deterministic causality',
        'physical_sqrt_V_lower': str(r_lower),
        'inner_radius': '3/20',
        'finite_residual_inner_retention_refuted_on_this_class': True,
        'scope': 'HISTORICAL V3 norm-only deterministic residuals, inherited real-arithmetic quiet source profile',
        'current_two_timescale_admissibility': 'OPEN; see imu-two-timescale-certificate.json',
        'current_model_storage_lower_bound_inferred': False,
        'all_time_float32_verified': False,
        'local_homogeneous_theorem_refuted': False,
        'end_to_end_stability_theorem_closed': False,
    }


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2, sort_keys=True))
