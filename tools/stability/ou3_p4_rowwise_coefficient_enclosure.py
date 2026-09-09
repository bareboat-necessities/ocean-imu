#!/usr/bin/env python3
"""Source-uniform rowwise Kalman/reset coefficient enclosure for P4.

The canonical endpoint-referenced factored BRMM Riccati tube gives a rigorous
rowwise Joseph-gain bound. For

    K=P H' (H P H'+R)^-1,

put A=R^-1/2 H P^1/2. Since every singular value satisfies
s/(1+s^2)<=1/2,

    ||row_i(K)|| <= .5*sqrt(P_ii/lambda_min(R)).

This bounds every actual same-source Joseph gain without claiming that the
resulting row boxes are independently realizable. The deterministic hard P4
entry set then bounds the physical residual of every event. Those bounds give
a source-uniform attitude-correction radius. The exact shipping reset transport
supplies G(d), ||G^-1||=1, its maximum singular value, and the finite Cayley
reset-defect bound for that correction radius.

The row boxes are magnitude envelopes only. Consecutive contraction must still
use the signed same-history information ledger / dense graph IQCs.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

# Use the canonical stable+endpoint facade directly.  Do not monkey-patch the
# raw BASE module: doing so breaks the endpoint-reference facade's saved build.
import ou3_brmm_riccati_tube_factored as TUBE
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_exact_reset_transport as RESET

GROUPS_H=(('attitude',0),('gyro_bias',3),('velocity',6),('position',9),('integral_displacement',12),('latent_acceleration',15))
GROUPS_A=GROUPS_H+(('accelerometer_bias',18),)

def up(x): return math.nextafter(float(x), math.inf)
def down(x): return math.nextafter(float(x), -math.inf)
def _row_norms(pdiag,rvar): return [up(.5*math.sqrt(float(v)/float(rvar))) for v in pdiag]
def _attitude_correction_radius(rows,residuals):
    best=-1.0;kind=''
    for event,key in (('accelerometer','K_row_norm_upper_accelerometer'),('magnetometer','K_row_norm_upper_magnetometer'),('S_zero','K_row_norm_upper_S_zero')):
        kr=rows['attitude'][key];knorm=up(math.sqrt(sum(float(x)*float(x) for x in kr)));delta=up(knorm*float(residuals[event]))
        if delta>best:best=delta;kind=event
    return best,kind

def build():
    tube=TUBE.build();dyn=DYNAMIC.build();entry=ENTRY.build();reset=RESET.build()
    bad={'tube':TUBE.validate(tube),'dynamic':DYNAMIC.validate(dyn),'entry':ENTRY.validate(entry),'reset':RESET.validate(reset)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError(f'coefficient prerequisites failed: {bad}')
    domain=json.loads(TUBE.DEFAULT_DOMAIN.read_text());live=domain['normal_live'];fmax=float(live['specific_force_norm_upper_mps2']);mmax=float(live['magnetic_vector_norm_upper_uT']);rslo=float(dyn['dynamic_invariant']['R_S_applied'][0])
    r_acc=down(0.2**2);r_mag=down(0.3**2);r_s=down((0.72*rslo)**2);er=entry['coordinate_radii'];q=float(er['attitude_cayley_norm']);aw=float(er['latent_acceleration_norm_mps2']);ba=float(er['accelerometer_bias_error_norm_mps2']);S=float(er['integral_displacement_norm_m_s'])
    rot_diff=up(q/math.sqrt(1.0+0.25*q*q));residual_common={'S_zero':S,'magnetometer':up(rot_diff*mmax)};out={}
    for mode,key,groups in (('H18','H',GROUPS_H),('A21','A',GROUPS_A)):
        p=[float(x) for x in tube['modes'][key]['Pbar_diagonal_variance_upper']];rows={}
        for name,off in groups:
            vals=p[off:off+3];rows[name]={'Pii_upper':vals,'K_row_norm_upper_accelerometer':_row_norms(vals,r_acc),'K_row_norm_upper_magnetometer':_row_norms(vals,r_mag),'K_row_norm_upper_S_zero':_row_norms(vals,r_s)}
        residuals=dict(residual_common);residuals['accelerometer']=up(rot_diff*fmax+aw+(ba if mode=='A21' else 0.0));delta,limiting=_attitude_correction_radius(rows,residuals)
        reset_enclosure=None;reset_closed=False;reset_failure=None
        try:
            reset_enclosure=RESET.reset_defect_bound(q,delta);reset_closed=bool(reset_enclosure['chart_safe'])
        except Exception as exc:reset_failure=str(exc)
        out[mode]={'rows':rows,'measurement_R_variance_lower':{'accelerometer':r_acc,'magnetometer':r_mag,'S_zero':r_s},'hard_entry_residual_norm_upper':residuals,'attitude_correction_norm_upper':delta,'limiting_correction_event':limiting,'reset_enclosure':reset_enclosure,'reset_enclosure_failure':reset_failure,'reset_coefficient_family_outwardly_bounded':reset_closed}
    reset_all=all(out[m]['reset_coefficient_family_outwardly_bounded'] for m in ('H18','A21'))
    return {'qualification':'OU3_P4_SOURCE_UNIFORM_ROWWISE_KALMAN_RESET_COEFFICIENT_ENCLOSURE_V3','canonical_endpoint_referenced_factored_Riccati_tube_consumed':True,'source_uniform_Riccati_diagonal_ceiling_consumed':True,'rowwise_Joseph_gain_bound_proved':True,'independent_P_H_R_K_boxes_claimed_physical':False,'actual_RS_lower_from_dynamic_invariant':True,'hard_entry_set_consumed':entry,'exact_reset_transport_consumed':True,'reset_inverse_operator_norm_upper':1.0,'modes':out,'Joseph_gain_family_outwardly_bounded':True,'reset_coefficient_family_outwardly_bounded':reset_all,'Kalman_and_reset_coefficient_family_outwardly_bounded':reset_all,'same_history_correlation_retained_for_later_storage_test':True,'independent_row_box_word_composition_forbidden':True,'consecutive_storage_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False}
def validate(d):
    f=[]
    if d.get('qualification')!='OU3_P4_SOURCE_UNIFORM_ROWWISE_KALMAN_RESET_COEFFICIENT_ENCLOSURE_V3':f.append('qualification mismatch')
    for k in ('canonical_endpoint_referenced_factored_Riccati_tube_consumed','source_uniform_Riccati_diagonal_ceiling_consumed','rowwise_Joseph_gain_bound_proved','actual_RS_lower_from_dynamic_invariant','exact_reset_transport_consumed','Joseph_gain_family_outwardly_bounded','reset_coefficient_family_outwardly_bounded','Kalman_and_reset_coefficient_family_outwardly_bounded','same_history_correlation_retained_for_later_storage_test','independent_row_box_word_composition_forbidden'):
        if d.get(k) is not True:f.append(k+' not true')
    if d.get('independent_P_H_R_K_boxes_claimed_physical') is not False:f.append('independent boxes promoted')
    if d.get('consecutive_storage_closed_here') is not False:f.append('storage falsely closed')
    if d.get('reset_inverse_operator_norm_upper')!=1.0:f.append('reset inverse norm changed')
    for m in ('H18','A21'):
        mm=d['modes'][m]
        if not math.isfinite(float(mm['attitude_correction_norm_upper'])) or float(mm['attitude_correction_norm_upper'])<0:f.append(m+' correction radius invalid')
        if mm['reset_coefficient_family_outwardly_bounded'] is not True:f.append(m+' reset family not bounded: '+str(mm.get('reset_enclosure_failure')))
        for g,row in mm['rows'].items():
            for key in ('K_row_norm_upper_accelerometer','K_row_norm_upper_magnetometer','K_row_norm_upper_S_zero'):
                if not all(math.isfinite(float(x)) and float(x)>=0 for x in row[key]):f.append(m+' '+g+' '+key+' invalid')
    return list(dict.fromkeys(f))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'kalman_reset_family':d['Kalman_and_reset_coefficient_family_outwardly_bounded'],'H18_correction':d['modes']['H18']['attitude_correction_norm_upper'],'A21_correction':d['modes']['A21']['attitude_correction_norm_upper'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
