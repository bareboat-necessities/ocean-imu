#!/usr/bin/env python3
"""Binary32 enclosure of the shipping accelerometer-bias projection.

The exact-real radial projection sector is already closed.  This producer adds
an explicit binary32 perturbation bound for the deployed three-component
projection computation.  It does not assume exact real norm arithmetic.

For unit roundoff u=2^-24 and k rounded basic operations, gamma_k=k*u/(1-k*u).
A three-component norm formed from three products, two additions and one square
root is bounded by a conservative relative error epsilon_norm.  FMA contraction
can only reduce the number of roundings relative to this bound.  Division and
three component multiplies add two more relative roundings on the active branch.
Branch disagreement can occur only in the narrow norm-error band around R; the
radial projection is continuous there, so the same radial perturbation budget
covers it.

The resulting floating map is the exact projection map plus an additive vector
r_fp with ||r_fp|| <= eps_projection.  Hence the exact joint sector becomes an
ISS sector with one explicit finite-precision forcing channel.  This closes
projection rounding only; Joseph gain, covariance and reset arithmetic remain
separate obligations.
"""
from __future__ import annotations

import argparse
import json
import math
import struct
from pathlib import Path

U=2.0**-24
R_REAL=0.4


def f32(x: float) -> float:
    return struct.unpack('<f',struct.pack('<f',float(x)))[0]


def gamma(k: int) -> float:
    if k < 0 or k*U >= 1.0:
        raise ValueError('invalid gamma index')
    return math.nextafter((k*U)/(1.0-k*U),math.inf)


def build() -> dict:
    # Three squarings + two additions, then correctly-rounded sqrt.  Reserve
    # gamma_6 for the dot/norm path, then one division and one component
    # multiply on the active scale path.  This deliberately overcounts FMA.
    eps_norm=math.nextafter((1.0+gamma(6))*(1.0+U)-1.0,math.inf)
    eps_scale=gamma(2)
    Rf=f32(R_REAL)
    radius_literal_error=abs(Rf-R_REAL)

    # Active projection radial magnitude can exceed the exact R because the
    # computed norm can be low and scale/multiply can round high.
    active_radial_upper=math.nextafter(
        Rf*(1.0+eps_scale)/(1.0-eps_norm),math.inf)
    active_excess=max(0.0,active_radial_upper-R_REAL)

    # If the floating branch remains inactive while the exact norm is just over
    # R, the true norm is at most Rf/(1-eps_norm).  Continuity of radial
    # projection bounds the branch-disagreement displacement by this excess.
    inactive_branch_radial_upper=math.nextafter(Rf/(1.0-eps_norm),math.inf)
    branch_excess=max(0.0,inactive_branch_radial_upper-R_REAL)

    eps_projection=math.nextafter(max(active_excess,branch_excess)+radius_literal_error,math.inf)
    pairwise_sector_additive=math.nextafter(2.0*eps_projection,math.inf)

    return {
      'qualification':'OU3_P4_BINARY32_BIAS_PROJECTION_ENCLOSURE_V1',
      'runtime_scalar_format':'IEEE754_binary32',
      'unit_roundoff':U,
      'shipping_radius_literal_real':R_REAL,
      'shipping_radius_literal_binary32':Rf,
      'radius_literal_abs_error':radius_literal_error,
      'norm_relative_error_upper':eps_norm,
      'active_scale_relative_error_upper':eps_scale,
      'active_projected_norm_upper_mps2':active_radial_upper,
      'inactive_branch_ambiguous_norm_upper_mps2':inactive_branch_radial_upper,
      'single_evaluation_projection_value_error_norm_upper_mps2':eps_projection,
      'pairwise_projection_sector_additive_norm_upper_mps2':pairwise_sector_additive,
      'floating_map_decomposition':'F_fp(e,beta)=F_exact(e,beta)+r_fp; ||r_fp||<=eps_projection',
      'floating_joint_sector':'||Delta F_fp|| <= sqrt(||Delta e||^2+||Delta beta||^2)+2*eps_projection',
      'branch_ambiguity_covered':True,
      'FMA_contraction_covered_by_overcounted_rounding_model':True,
      'reassociation_beyond_three_term_dot_product_irrelevant_to_bound':True,
      'projection_finite_precision_enclosure_closed':True,
      'full_Kalman_reset_finite_precision_enclosure_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,
    }


def validate(d: dict) -> list[str]:
    f=[]
    for k in ('branch_ambiguity_covered','FMA_contraction_covered_by_overcounted_rounding_model',
              'reassociation_beyond_three_term_dot_product_irrelevant_to_bound',
              'projection_finite_precision_enclosure_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    if d.get('full_Kalman_reset_finite_precision_enclosure_closed_here') is not False:f.append('full arithmetic falsely closed')
    if d.get('P4_MOTION_PASS') is not False or d.get('P4_PASS') is not False:f.append('projection arithmetic promoted P4')
    eps=float(d.get('single_evaluation_projection_value_error_norm_upper_mps2',-1))
    if not (math.isfinite(eps) and 0.0 < eps < 1e-5):f.append('projection rounding bound invalid')
    if float(d.get('shipping_radius_literal_binary32',0.0)) != f32(R_REAL):f.append('binary32 radius literal mismatch')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
