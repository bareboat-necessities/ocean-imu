"""Source-uniform branch certificate for shipping prediction quaternions.

Under the declared regional Normal-Live domain, every nominal and true-shadow
5 ms prediction increment is strictly below quat_from_delta_theta's 0.01 rad
branch threshold.  Therefore the finite prediction graph needs only the exact
polynomial (w,k*d) branch in real arithmetic; the trig branch cannot occur
inside this regional theorem domain.

This closes only branch selection. Binary32 FMA/normalization roundoff remains
a separate finite-precision obligation.
"""
from __future__ import annotations
import json,math
from pathlib import Path

REPO=Path(__file__).resolve().parents[3]
OPERATING=REPO/'tools/stability/ou3_proof_operating_domain.json'
CLOSURE=REPO/'tools/stability/ou3_p4_closure_domain.json'
CORE=REPO/'src/kalman_ou_common/KalmanOUCoreMath.h'
QUALIFICATION='OU3_ALT_PREDICTION_QUATERNION_SMALL_BRANCH_V1'
THRESHOLD_RAD=1.0e-2


def upward(x): return math.nextafter(float(x),math.inf)


def small_branch_parameters(dtheta2):
    """Real-arithmetic polynomial branch in terms of ||dtheta||^2."""
    u=float(dtheta2)
    if not (math.isfinite(u) and 0<=u<THRESHOLD_RAD**2):
        raise ValueError('increment outside strict polynomial quaternion branch')
    u2=u*u
    return 1-u/8+u2/384, .5-u/48+u2/3840


def build():
    op=json.loads(OPERATING.read_text());cl=json.loads(CLOSURE.read_text());src=CORE.read_text()
    h=upward(op['configured_runtime']['imu_dt_s'])
    omega=upward(math.radians(op['normal_live']['body_rate_norm_upper_deg_s']))
    bg=upward(cl['hard_entry_search']['base_coordinate_radii']['gyro_bias_norm_rad_s'])
    nominal=upward(omega*h)
    shadow=upward((omega+bg)*h)
    threshold_source=('if (theta < T(1e-2))' in src)
    source_prediction_call=('quat_from_delta_theta((-last_gyr_bias_corrected * Ts).eval())' in (REPO/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h').read_text())
    closed=bool(threshold_source and source_prediction_call and nominal<THRESHOLD_RAD and shadow<THRESHOLD_RAD)
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'runtime_dt_s_upper':h,'body_rate_norm_upper_rad_s':omega,'regional_gyro_bias_error_norm_upper_rad_s':bg,
      'nominal_step_angle_norm_upper_rad':nominal,'shadow_step_angle_norm_upper_rad':shadow,
      'shipping_quaternion_branch_threshold_rad':THRESHOLD_RAD,
      'shadow_margin_to_threshold_rad':math.nextafter(THRESHOLD_RAD-shadow,-math.inf),
      'shipping_threshold_source_parity':threshold_source,'shipping_prediction_call_source_parity':source_prediction_call,
      'nominal_prediction_always_polynomial_branch':bool(nominal<THRESHOLD_RAD),
      'shadow_prediction_always_polynomial_branch':bool(shadow<THRESHOLD_RAD),
      'source_uniform_prediction_quaternion_branch_closed':closed,
      'trigonometric_prediction_branch_reachable_in_regional_domain':False if closed else None,
      'binary32_polynomial_and_normalization_roundoff_closed':False,
      'ALT_LIVE_PASS':False,
      'next_obligation':'substitute the polynomial step-quaternion hard graph into every finite prediction descriptor, retaining omega, e_bg and h ancestry; then attach q15 and bias-root source relations and covariance successors',
    }


def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('shipping_threshold_source_parity','shipping_prediction_call_source_parity','nominal_prediction_always_polynomial_branch','shadow_prediction_always_polynomial_branch','source_uniform_prediction_quaternion_branch_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    if d.get('trigonometric_prediction_branch_reachable_in_regional_domain') is not False:f.append('trig branch not excluded')
    for k in ('binary32_polynomial_and_normalization_roundoff_closed','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if not 0<d.get('shadow_margin_to_threshold_rad',-1):f.append('no strict branch margin')
    return f
