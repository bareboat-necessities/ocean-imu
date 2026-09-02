#!/usr/bin/env python3
"""Global Cayley-chart geometry certificate for a usable OU-III P4 sector.

The deployed attitude error coordinate is

    c = 2 tan(theta/2) u,

and the exact inverse-Cayley rotation is

    R(c) = I + 4/(4+||c||^2)[c]x + 2/(4+||c||^2)[c]x^2.

For an infinitesimal coordinate perturbation dc, the corresponding spatial
rotation perturbation is

    dphi = A(c) dc,
    A(c) = 4/(4+||c||^2) (I + [c]x/2).

Hence A has one singular value 1/(1+||c||^2/4) along c and two singular values
1/sqrt(1+||c||^2/4) transverse to c.  On theta<=theta_o this gives

    sigma_min(A) >= cos(theta_o/2)^2.

For any collection of accepted inertial vectors v_i with positive scalar
weights w_i, the exact residual Jacobian is J_i=-[R v_i]x A. Therefore

    sum w_i J_i'J_i
      = A' R (sum w_i [v_i]x'[v_i]x) R' A
      >= sigma_min(A)^2 * lambda_min(G_v) I.

This is a source-independent nonlinear geometry statement: body-frame gauge,
vector magnitudes and vector separation enter only through the same vector
Gramian G_v already certified by the PE producer. It does NOT by itself prove
the complete Joseph word, because the implemented EKF uses a linearized H and
the signed nonlinear remainder still has to be charged. Its purpose is to
establish that a 45-degree-scale P4 sector is not lost to Cayley conditioning.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
DEFAULT_OUTER_ANGLE_RAD = 0.80


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def build(domain_path: Path=DEFAULT_DOMAIN, outer_angle_rad: float=DEFAULT_OUTER_ANGLE_RAD):
    domain=json.loads(Path(domain_path).resolve().read_text(encoding="utf-8"))
    theta=float(outer_angle_rad)
    if not (0.0 < theta < math.pi):
        raise ValueError("outer Cayley angle must lie in (0,pi)")
    half=I(0.5*theta)
    s=VT.sin_interval(half)
    c=VT.cos_interval(half)
    if c.lo <= 0.0:
        raise RuntimeError("outer Cayley sector reaches chart antipode")
    q=(I(2.0)*s/c)
    q_upper=q.hi
    q2=up(q_upper*q_upper)
    chart_min=down(4.0/up(4.0+q2))
    chart_max=up(1.0/math.sqrt(down(1.0+q.lo*q.lo/4.0)))
    info_factor=down(chart_min*chart_min)

    # Exact operator identities at the sector boundary. These are diagnostics
    # for the later signed-Joseph remainder calculation, not a contraction claim.
    rotation_minus_identity=up(2.0*s.hi)
    chart_minus_identity=up(s.hi)  # ||A-I||_2 = sin(theta/2)

    entrance_deg=float(domain["initial_filter_entrance"]["attitude"]["full_attitude_error_upper_deg"])
    entrance_rad=up(entrance_deg*math.pi/180.0)
    covers_entrance=entrance_rad <= theta
    usable=bool(
        covers_entrance and chart_min >= 0.80 and info_factor >= 0.64
        and q_upper < 1.0
    )
    return {
        "schema":SCHEMA,
        "qualification":"OU3_P4_GLOBAL_CAYLEY_SECTOR_GEOMETRY",
        "source_generated_not_trajectory_fit":True,
        "trajectory_replay_used":False,
        "filter_changed":False,
        "outer_angle_rad":theta,
        "outer_angle_deg":theta*180.0/math.pi,
        "declared_filter_entrance_attitude_deg":entrance_deg,
        "declared_filter_entrance_covered":covers_entrance,
        "cayley_radius_upper":q_upper,
        "chart_antipode_excluded":q_upper < math.inf and theta < math.pi,
        "cayley_spatial_differential_formula":"A=4/(4+||c||^2)*(I+[c]x/2)",
        "chart_sigma_min_lower":chart_min,
        "chart_sigma_max_upper":chart_max,
        "exact_vector_information_retention_factor_lower":info_factor,
        "vector_information_relation":"G_exact(c)>=retention_factor*lambda_min(G_vector)*I",
        "rotation_minus_identity_norm_upper":rotation_minus_identity,
        "chart_differential_minus_identity_norm_upper":chart_minus_identity,
        "ordinary_float_eigensolver_used":False,
        "full_18_21_state_Joseph_word_established_here":False,
        "signed_EKF_remainder_charged_here":False,
        "usable_sector_geometry_pass":usable,
        "pass":usable,
    }


def validate(d):
    f=[]
    for k in ("source_generated_not_trajectory_fit","declared_filter_entrance_covered","chart_antipode_excluded","usable_sector_geometry_pass","pass"):
        if d.get(k) is not True: f.append(f"{k} is not true")
    for k in ("trajectory_replay_used","filter_changed","ordinary_float_eigensolver_used","full_18_21_state_Joseph_word_established_here","signed_EKF_remainder_charged_here"):
        if d.get(k) is not False: f.append(f"{k} is not false")
    for k in ("chart_sigma_min_lower","chart_sigma_max_upper","exact_vector_information_retention_factor_lower"):
        x=d.get(k)
        if not isinstance(x,(int,float)) or not math.isfinite(float(x)) or float(x)<=0.0:
            f.append(f"{k} is not finite positive")
    if float(d.get("chart_sigma_min_lower",0.0)) < 0.80:
        f.append("Cayley sector conditioning is below usable 0.80 floor")
    if float(d.get("exact_vector_information_retention_factor_lower",0.0)) < 0.64:
        f.append("exact vector information retention is below usable 0.64 floor")
    return f


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN)
    ap.add_argument("--outer-angle-rad",type=float,default=DEFAULT_OUTER_ANGLE_RAD)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    d=build(a.domain,a.outer_angle_rad)
    f=validate(d)
    d["validation_pass"]=not f
    d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps(d,indent=2,sort_keys=True))
    return 0 if not f else 2


if __name__=="__main__":
    raise SystemExit(main())
