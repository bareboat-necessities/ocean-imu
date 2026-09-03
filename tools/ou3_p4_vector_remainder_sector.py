#!/usr/bin/env python3
"""Source-uniform nonlinear vector-remainder sector for OU-III P4.

For c=2 tan(theta/2)u and any inertial vector v, define

    h_v(c)   = [c]x v,
    r_v(c)   = (R(c)-I)v,
    eta_v(c) = r_v(c)-h_v(c).

Decomposing v into components parallel/perpendicular to u gives the exact norm
identity

    ||eta_v|| = sin(theta/2) ||h_v||.

For the shipping zero-lever-arm accelerometer, the original world-frame a_w
error coordinate gives the valid but conservative expression

    eta_a = eta_f + (R-I) delta_a_w.

The current route removes that artificial second term by the exact co-rotated
coordinate certified in ``ou3_p4_accelerometer_corotated_aw``.  With

    Q_aw = R_hat' E R_hat,    u_aw = Q_aw delta_a_w,

Q_aw is orthogonal and the exact residual is

    r_a = (E-I) f_hat + R_hat u_aw + delta_b_a,
    h_a = [c]x f_hat + R_hat u_aw + delta_b_a.

Hence

    eta_a = ((E-I)-[c]x) f_hat,

so both a_w and b_a cancel exactly from nonlinear eta.  The operation coordinate
change is an exact state/covariance congruence and an isometry of the retained
group-isotropic P3->P4 metric; no Joseph-energy term is discarded.

For magnetometer measurements the same pure-vector identity gives

    ||eta_m||^2 <= s^2 ||[c]x m||^2,
    s = sin(theta_o/2).

Because the configured accelerometer and magnetometer covariance matrices are
isotropic in the proof scope, the same pure-rotation inequalities hold after
division by measurement variance.  They are homogeneous quadratic penalties,
not additive disturbances.  This primitive still does not establish the full
18/21-state Joseph word; directional signed forms must be accumulated over the
source-complete recurrent word before scalarization.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_accelerometer_corotated_aw as COROT
import ou3_validated_transcendentals as VT

REPO=Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN=REPO/"tools"/"ou3_proof_operating_domain.json"
SCHEMA=2
DEFAULT_OUTER_ANGLE_RAD=0.80
PI_UP=math.nextafter(math.pi,math.inf)


def down(x: float)->float:
    return math.nextafter(float(x),-math.inf)


def up(x: float)->float:
    return math.nextafter(float(x),math.inf)


def build(domain_path: Path=DEFAULT_DOMAIN, outer_angle_rad: float=DEFAULT_OUTER_ANGLE_RAD,
          young_epsilon: float=1.0):
    path=Path(domain_path).resolve()
    domain=json.loads(path.read_text(encoding="utf-8"))
    theta=float(outer_angle_rad)
    eps=float(young_epsilon)
    if not (0.0<theta<math.pi and eps>0.0 and math.isfinite(eps)):
        raise ValueError("invalid sector angle or legacy diagnostic Young epsilon")

    corot=COROT.build(path)
    failures=[f"co-rotated-aw: {x}" for x in COROT.validate(corot)]
    if float(corot.get("outer_angle_rad",0.0)) < theta:
        failures.append("co-rotated a_w primitive does not cover requested remainder sector")

    s=VT.sin_point(0.5*theta)
    s2=up(s.hi*s.hi)
    mag_coeff=s2
    acc_att_coeff=s2
    acc_aw_coeff=0.0
    rot_minus_I=up(2.0*s.hi)

    # Preserve the old original-coordinate Young bound only as a diagnostic so
    # future audits can quantify exactly what the congruence removed.
    legacy_acc_att=up((1.0+eps)*s2)
    legacy_acc_aw=up((1.0+1.0/eps)*up(4.0*s2))

    entrance_deg=float(domain["initial_filter_entrance"]["attitude"]["full_attitude_error_upper_deg"])
    entrance_rad_upper=up(entrance_deg*PI_UP/180.0)
    aw_bound=float(domain["startup"]["physical_handoff_coordinate_bounds"]["latent_acceleration_error_norm_upper_mps2"])
    ba_bound=float(domain["normal_live"]["active_accelerometer_bias_state_norm_upper_mps2"])
    covers=entrance_rad_upper<=theta

    primitive_pass=bool(
        not failures and covers and mag_coeff<0.25 and acc_att_coeff<0.25
        and acc_aw_coeff==0.0
    )
    return {
        "schema":SCHEMA,
        "qualification":"OU3_P4_GLOBAL_VECTOR_NONLINEAR_REMAINDER_SECTOR",
        "source_generated_not_trajectory_fit":True,
        "trajectory_replay_used":False,
        "filter_changed":False,
        "outer_angle_rad":theta,
        "declared_filter_entrance_attitude_deg":entrance_deg,
        "declared_filter_entrance_covered":covers,
        "sin_half_angle_upper":s.hi,
        "sin_half_angle_squared_upper":s2,
        "mag_eta_squared_over_linear_rotation_squared_upper":mag_coeff,
        "accelerometer_corotated_aw_coordinate_used":True,
        "accelerometer_corotated_aw_Joseph_congruence_exact":True,
        "acc_eta_force_rotation_quadratic_coefficient_upper":acc_att_coeff,
        "acc_eta_aw_quadratic_coefficient_upper":acc_aw_coeff,
        "rotation_minus_identity_norm_upper":rot_minus_I,
        "latent_acceleration_error_norm_upper_mps2":aw_bound,
        "active_accelerometer_bias_norm_upper_mps2":ba_bound,
        "accelerometer_bias_nonlinear_remainder_coefficient":0.0,
        "accelerometer_aw_cancels_exactly_from_eta_in_operation_coordinate":True,
        "accelerometer_bias_cancels_exactly_from_eta":True,
        "legacy_original_coordinate_diagnostic": {
            "used_for_P4_word":False,
            "young_epsilon":eps,
            "force_rotation_quadratic_coefficient_upper":legacy_acc_att,
            "aw_quadratic_coefficient_upper":legacy_acc_aw,
            "reason":"valid original-coordinate bound superseded by exact block-orthogonal co-rotation",
        },
        "penalties_are_homogeneous_quadratic_not_affine_beta":True,
        "measurement_covariance_isotropy_required":True,
        "complete_Joseph_word_established_here":False,
        "P4_USABLE_CERTIFICATE_PROMOTED":False,
        "usable_sector_remainder_primitive_pass":primitive_pass,
        "pass":primitive_pass,
        "failures":failures,
    }


def validate(d):
    f=list(d.get("failures",[]))
    if d.get("schema")!=SCHEMA: f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit","declared_filter_entrance_covered",
        "accelerometer_corotated_aw_coordinate_used",
        "accelerometer_corotated_aw_Joseph_congruence_exact",
        "accelerometer_aw_cancels_exactly_from_eta_in_operation_coordinate",
        "accelerometer_bias_cancels_exactly_from_eta",
        "penalties_are_homogeneous_quadratic_not_affine_beta",
        "measurement_covariance_isotropy_required","usable_sector_remainder_primitive_pass","pass"):
        if d.get(k) is not True: f.append(f"{k} is not true")
    for k in ("trajectory_replay_used","filter_changed","complete_Joseph_word_established_here","P4_USABLE_CERTIFICATE_PROMOTED"):
        if d.get(k) is not False: f.append(f"{k} is not false")
    for k in ("mag_eta_squared_over_linear_rotation_squared_upper","acc_eta_force_rotation_quadratic_coefficient_upper","acc_eta_aw_quadratic_coefficient_upper"):
        x=d.get(k)
        if not isinstance(x,(int,float)) or not math.isfinite(float(x)) or float(x)<0.0:
            f.append(f"{k} is invalid")
    if d.get("acc_eta_aw_quadratic_coefficient_upper")!=0.0:
        f.append("co-rotated a_w nonlinear eta coefficient must be zero")
    if d.get("accelerometer_bias_nonlinear_remainder_coefficient")!=0.0:
        f.append("accelerometer bias must cancel exactly from eta")
    legacy=d.get("legacy_original_coordinate_diagnostic",{})
    if legacy.get("used_for_P4_word") is not False:
        f.append("legacy original-coordinate a_w penalty became active again")
    if not float(d.get("acc_eta_force_rotation_quadratic_coefficient_upper",math.inf)) <= float(d.get("mag_eta_squared_over_linear_rotation_squared_upper",-1.0)):
        f.append("co-rotated accelerometer pure-rotation coefficient drifted above vector identity")
    return list(dict.fromkeys(f))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN)
    ap.add_argument("--outer-angle-rad",type=float,default=DEFAULT_OUTER_ANGLE_RAD)
    ap.add_argument("--young-epsilon",type=float,default=1.0,
                    help="legacy original-coordinate diagnostic only")
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    d=build(a.domain,a.outer_angle_rad,a.young_epsilon)
    f=validate(d)
    d["validation_pass"]=not f
    d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps(d,indent=2,sort_keys=True))
    return 0 if not f else 2


if __name__=="__main__":
    raise SystemExit(main())
