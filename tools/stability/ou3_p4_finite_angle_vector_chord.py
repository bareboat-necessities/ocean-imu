#!/usr/bin/env python3
"""Exact finite-angle vector-residual sector in deployed Cayley coordinates.

For c=2*tan(theta/2)*u and any body/reference vector m,

    ||(E(c)-I)m||^2
      = 4 sin^2(theta/2) ||u x m||^2
      = ||c x m||^2 / (1+||c||^2/4).

Thus on ||c||<=c_max,

    alpha(c_max) ||[m]x c||^2
      <= ||(E(c)-I)m||^2
      <= ||[m]x c||^2,
    alpha(c_max)=1/(1+c_max^2/4).

This is an exact global sector on the finite attitude chart, not a Taylor
remainder.  It removes the catastrophic standalone-eta scalarization for pure
vector orientation residuals.  With isotropic positive measurement weights,
the same factor applies to any sum of non-collinear vector information terms.
The complete P4 proof still has to retain accelerometer a_w/b_a coupling,
source-time propagation, Kalman gain mismatch and reset terms jointly.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_hard_entry_set as ENTRY

QUALIFICATION='OU3_P4_FINITE_ANGLE_CAYLEY_VECTOR_CHORD_SECTOR_V1'


def chord_factor(c_norm_upper: float) -> float:
    q=float(c_norm_upper)
    if not math.isfinite(q) or q<0.0: raise ValueError('finite nonnegative Cayley radius required')
    return math.nextafter(1.0/(1.0+0.25*q*q),-math.inf)


def build() -> dict:
    e=ENTRY.build(); ef=ENTRY.validate(e)
    if ef: raise RuntimeError('hard entry prerequisite failed: '+repr(ef))
    q=float(e['coordinate_radii']['attitude_cayley_norm'])
    alpha=chord_factor(q)
    theta=2.0*math.atan(q/2.0)
    return {
      'qualification':QUALIFICATION,
      'coordinate':'c=2*tan(theta/2)*u',
      'identity':'||(E(c)-I)m||^2=||c x m||^2/(1+||c||^2/4)',
      'full_declared_entry_cayley_radius':q,
      'full_declared_entry_angle_deg':math.degrees(theta),
      'finite_angle_linear_information_retention_factor_lower':alpha,
      'upper_factor':1.0,
      'valid_for_every_vector_m':True,
      'valid_for_every_c_in_full_declared_entry_ball':True,
      'Taylor_remainder_used':False,
      'finite_harmonic_or_replay_fit_used':False,
      'isotropic_weighted_sum_inherits_factor':True,
      'pure_vector_orientation_residual_standalone_eta_budget_required':False,
      'accelerometer_aw_ba_joint_coupling_closed_here':False,
      'same_history_prediction_and_reset_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
    }


def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('valid_for_every_vector_m','valid_for_every_c_in_full_declared_entry_ball','isotropic_weighted_sum_inherits_factor'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('Taylor_remainder_used','finite_harmonic_or_replay_fit_used','pure_vector_orientation_residual_standalone_eta_budget_required','accelerometer_aw_ba_joint_coupling_closed_here','same_history_prediction_and_reset_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    a=float(d.get('finite_angle_linear_information_retention_factor_lower',0))
    if not (0.85<a<1.0):f.append('unexpected full-entry chord factor')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,indent=2,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
