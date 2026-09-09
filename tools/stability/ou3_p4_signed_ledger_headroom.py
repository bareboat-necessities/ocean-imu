#!/usr/bin/env python3
"""Non-promoting headroom diagnostic for the signed Joseph/reset P4 route.

This does not replace the full same-history S-procedure.  It asks a simpler
question first: on the FULL declared P4 entry cell, how large are the exact
Cayley residual-sector coefficients relative to the already-certified H18
complete-word directional measurement information?

The comparison is intentionally pessimistic: it uses R^-1 residual energy and
the complete-word linear information lower.  Failure does not disprove P4; it
means scalarizing the nonlinear residual sector against the weakest complete-
word information direction loses too much correlation.  Success would be a
useful sufficient headroom indication before assembling the full joint LDLT.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_brmm_h18_information_composition as INFO
import ou3_p4_complete_brmm_residual_sector as RESIDUAL
import ou3_p4_hard_entry_set as ENTRY

QUALIFICATION='OU3_P4_SIGNED_LEDGER_FULL_ENTRY_HEADROOM_V1'


def build() -> dict:
    info=INFO.build(); residual=RESIDUAL.build(); entry=ENTRY.build()
    bad={'info':INFO.validate(info),'residual':RESIDUAL.validate(residual),'entry':ENTRY.validate(entry)}
    bad={k:v for k,v in bad.items() if v}
    if bad: raise RuntimeError(f'headroom prerequisites failed: {bad}')
    q=float(entry['coordinate_radii']['attitude_cayley_norm'])
    # Configured Racc=0.2^2, Rmag=0.3^2; these are the same theorem runtime
    # values consumed by the canonical P3/H18 information chain.
    racc_inv=1.0/(0.2**2); rmag_inv=1.0/(0.3**2)
    acc=RESIDUAL.sector_quadratic_diagonal(q,13.80665,racc_inv,include_aw=True)
    mag=RESIDUAL.sector_quadratic_diagonal(q,200.0,rmag_inv,include_aw=False)
    tri=info['triangular_information_composition']
    d=float(tri['D_H18_lambda_min_lower'])
    alpha=float(info['eta6_information_lower'])
    directional={k:float(v) for k,v in info['directional_translation_information_lower'].items()}
    # Two PE vector occurrences are the minimum declared witness.  This is a
    # diagnostic scalarization, NOT the final same-history ledger accounting.
    att_eta_two=2.0*(3.0*acc[0]+3.0*mag[0])
    aw_eta_two=2.0*(3.0*acc[3])
    return {
      'qualification':QUALIFICATION,
      'non_promoting':True,
      'full_declared_entry_cayley_radius':q,
      'H18_complete_word_information_lambda_min_lower':d,
      'H18_eta6_information_lower':alpha,
      'directional_translation_information_lower':directional,
      'accelerometer_residual_sector_diagonal':acc,
      'magnetometer_residual_sector_diagonal':mag,
      'two_required_PE_occurrences_attitude_eta_trace_coefficient':att_eta_two,
      'two_required_PE_occurrences_aw_eta_trace_coefficient':aw_eta_two,
      'attitude_eta_to_eta6_information_ratio':math.nextafter(att_eta_two/alpha,math.inf),
      'aw_eta_to_aw_directional_information_ratio':math.nextafter(aw_eta_two/directional['g^3*a_w'],math.inf),
      'weakest_scalar_information_ratio':math.nextafter(max(att_eta_two,aw_eta_two)/d,math.inf),
      'scalar_full_entry_residual_domination_indicated':max(att_eta_two,aw_eta_two)<d,
      'failure_interpretation':'scalarized residual-vs-weakest-information failure means retain event/directional correlation; it is not a P4 counterexample',
      'same_history_joint_LDLT_still_required':True,
      'reset_cross_defect_not_included_here':True,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
    }


def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('non_promoting','same_history_joint_LDLT_still_required','reset_cross_defect_not_included_here'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' falsely promoted')
    for k in ('attitude_eta_to_eta6_information_ratio','aw_eta_to_aw_directional_information_ratio','weakest_scalar_information_ratio'):
        if not (math.isfinite(float(d.get(k,math.nan))) and float(d[k])>=0):f.append(k+' invalid')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,indent=2,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
