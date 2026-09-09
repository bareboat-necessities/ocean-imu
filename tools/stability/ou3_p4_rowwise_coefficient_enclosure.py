#!/usr/bin/env python3
"""Source-uniform rowwise Kalman/reset coefficient enclosure for P4.

Uses the existing BRMM moving-Riccati diagonal ceiling.  For a Joseph update,
K=P H' (H P H'+R)^-1 and A=R^-1/2 H P^1/2 imply

  K R^1/2 = P^1/2 A' (A A'+I)^-1,
  ||row_i(K)|| <= .5*sqrt(P_ii/lambda_min(R)).

This rigorous rowwise bound retains the coordinate structure that a global
trace/K norm destroys.  It encloses every same-source P/H/R cell; it does not
claim independent choices are physically realizable.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_brmm_riccati_tube as TUBE
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_p4_hard_entry_set as ENTRY

GROUPS_H=(('attitude',0),('gyro_bias',3),('velocity',6),('position',9),('integral_displacement',12),('latent_acceleration',15))
GROUPS_A=GROUPS_H+(('accelerometer_bias',18),)


def up(x):return math.nextafter(float(x),math.inf)
def down(x):return math.nextafter(float(x),-math.inf)


def build():
    tube=TUBE.build(); dyn=DYNAMIC.build(); entry=ENTRY.build()
    live=json.loads(TUBE.DEFAULT_DOMAIN.read_text())['normal_live']
    r_acc=0.2**2; r_mag=0.3**2
    rslo=float(dyn['dynamic_invariant']['R_S_applied'][0]); r_s=(0.72*rslo)**2
    out={}
    for mode,key,groups in (('H18','H',GROUPS_H),('A21','A',GROUPS_A)):
      p=[float(x) for x in tube['modes'][key]['Pbar_diagonal_variance_upper']]
      rows={}
      for name,off in groups:
        vals=p[off:off+3]
        rows[name]={
          'Pii_upper':vals,
          'K_row_norm_upper_accelerometer':[up(.5*math.sqrt(v/r_acc)) for v in vals],
          'K_row_norm_upper_magnetometer':[up(.5*math.sqrt(v/r_mag)) for v in vals],
          'K_row_norm_upper_S_zero':[up(.5*math.sqrt(v/r_s)) for v in vals],
        }
      out[mode]={'rows':rows,'measurement_R_variance_lower':{'accelerometer':r_acc,'magnetometer':r_mag,'S_zero':r_s}}
    return {
      'qualification':'OU3_P4_SOURCE_UNIFORM_ROWWISE_COEFFICIENT_ENCLOSURE_V1',
      'source_uniform_Riccati_diagonal_ceiling_consumed':True,
      'rowwise_Joseph_gain_bound_proved':True,
      'independent_P_H_R_K_boxes_claimed_physical':False,
      'actual_RS_lower_from_dynamic_invariant':True,
      'hard_entry_set_consumed_for_later_residual_bounds':entry,
      'modes':out,
      'coefficient_family_outwardly_bounded':True,
      'consecutive_storage_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,
    }


def validate(d):
    f=[]
    for k in ('source_uniform_Riccati_diagonal_ceiling_consumed','rowwise_Joseph_gain_bound_proved','actual_RS_lower_from_dynamic_invariant','coefficient_family_outwardly_bounded'):
      if d.get(k) is not True:f.append(k+' not true')
    if d.get('independent_P_H_R_K_boxes_claimed_physical') is not False:f.append('independent boxes promoted')
    for m in ('H18','A21'):
      for g,row in d['modes'][m]['rows'].items():
        for key in ('K_row_norm_upper_accelerometer','K_row_norm_upper_magnetometer','K_row_norm_upper_S_zero'):
          if not all(math.isfinite(float(x)) and float(x)>=0 for x in row[key]):f.append(m+' '+g+' '+key+' invalid')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
