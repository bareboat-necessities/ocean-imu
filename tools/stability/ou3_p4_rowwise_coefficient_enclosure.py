#!/usr/bin/env python3
"""Source-uniform rowwise Joseph-gain magnitude enclosure for P4.

The endpoint-referenced factored Riccati tube gives

  ||row_i(K)|| <= .5*sqrt(P_ii/lambda_min(R)).

for every actual Joseph gain. These row bounds are valid magnitude envelopes
for finite-precision accounting. They are deliberately NOT composed through a
word and are NOT used to choose a reset correction domain: doing so destroys
the same-history K/H correlation (notably for S=0) and produces uselessly huge
attitude-correction radii.

Reset G(d) has ||G^-1||=1 for every finite d, but its finite Cayley graph sector
requires a bound on d proved from the same augmented D=E_theta*K*Y graph. That
obligation lives in ou3_p4_affine_hard_tube_iqc.py / ou3_p4_reset_graph_iqc.py.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_brmm_riccati_tube_factored as TUBE
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_p4_hard_entry_set as ENTRY

GROUPS_H=(('attitude',0),('gyro_bias',3),('velocity',6),('position',9),('integral_displacement',12),('latent_acceleration',15))
GROUPS_A=GROUPS_H+(('accelerometer_bias',18),)
def up(x):return math.nextafter(float(x),math.inf)
def down(x):return math.nextafter(float(x),-math.inf)
def _row_norms(pdiag,rvar):return [up(.5*math.sqrt(float(v)/float(rvar))) for v in pdiag]
def _diagnostic_correction_radius(rows,residuals):
    best=-1.;kind=''
    for event,key in (('accelerometer','K_row_norm_upper_accelerometer'),('magnetometer','K_row_norm_upper_magnetometer'),('S_zero','K_row_norm_upper_S_zero')):
        kr=rows['attitude'][key];delta=up(math.sqrt(sum(float(x)*float(x) for x in kr))*float(residuals[event]))
        if delta>best:best=delta;kind=event
    return best,kind

def build():
    tube=TUBE.build();dyn=DYNAMIC.build();entry=ENTRY.build();bad={'tube':TUBE.validate(tube),'dynamic':DYNAMIC.validate(dyn),'entry':ENTRY.validate(entry)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('coefficient prerequisites failed: '+repr(bad))
    domain=json.loads(TUBE.DEFAULT_DOMAIN.read_text());live=domain['normal_live'];fmax=float(live['specific_force_norm_upper_mps2']);mmax=float(live['magnetic_vector_norm_upper_uT']);rslo=float(dyn['dynamic_invariant']['R_S_applied'][0])
    r_acc=down(.2**2);r_mag=down(.3**2);r_s=down((.72*rslo)**2);er=entry['coordinate_radii'];q=float(er['attitude_cayley_norm']);aw=float(er['latent_acceleration_norm_mps2']);ba=float(er['accelerometer_bias_error_norm_mps2']);S=float(er['integral_displacement_norm_m_s']);rot_diff=up(q/math.sqrt(1+.25*q*q));common={'S_zero':S,'magnetometer':up(rot_diff*mmax)};out={}
    for mode,key,groups in (('H18','H',GROUPS_H),('A21','A',GROUPS_A)):
        p=[float(x) for x in tube['modes'][key]['Pbar_diagonal_variance_upper']];rows={}
        for name,off in groups:
            vals=p[off:off+3];rows[name]={'Pii_upper':vals,'K_row_norm_upper_accelerometer':_row_norms(vals,r_acc),'K_row_norm_upper_magnetometer':_row_norms(vals,r_mag),'K_row_norm_upper_S_zero':_row_norms(vals,r_s)}
        residuals=dict(common);residuals['accelerometer']=up(rot_diff*fmax+aw+(ba if mode=='A21' else 0.0));delta,limiting=_diagnostic_correction_radius(rows,residuals)
        out[mode]={'rows':rows,'measurement_R_variance_lower':{'accelerometer':r_acc,'magnetometer':r_mag,'S_zero':r_s},'hard_entry_residual_norm_upper_diagnostic':residuals,'independent_rowbox_attitude_correction_norm_upper_diagnostic':delta,'independent_rowbox_limiting_event_diagnostic':limiting}
    return {'qualification':'OU3_P4_SOURCE_UNIFORM_ROWWISE_JOSEPH_GAIN_ENCLOSURE_V4','canonical_endpoint_referenced_factored_Riccati_tube_consumed':True,'source_uniform_Riccati_diagonal_ceiling_consumed':True,'rowwise_Joseph_gain_bound_proved':True,'Joseph_gain_family_outwardly_bounded':True,'actual_RS_lower_from_dynamic_invariant':True,'hard_entry_set_consumed':entry,'independent_P_H_R_K_boxes_claimed_physical':False,'rowwise_K_used_for_finite_precision_magnitude_only':True,'rowwise_K_reset_correction_domain_forbidden':True,'reset_inverse_operator_norm_upper_for_any_finite_correction':1.0,'reset_coefficient_family_requires_same_graph_correction_domain':True,'Kalman_and_reset_coefficient_family_outwardly_bounded':False,'same_history_correlation_retained_for_storage_test':True,'independent_row_box_word_composition_forbidden':True,'modes':out,'consecutive_storage_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False}
def validate(d):
    f=[]
    if d.get('qualification')!='OU3_P4_SOURCE_UNIFORM_ROWWISE_JOSEPH_GAIN_ENCLOSURE_V4':f.append('qualification mismatch')
    for k in ('canonical_endpoint_referenced_factored_Riccati_tube_consumed','source_uniform_Riccati_diagonal_ceiling_consumed','rowwise_Joseph_gain_bound_proved','Joseph_gain_family_outwardly_bounded','actual_RS_lower_from_dynamic_invariant','rowwise_K_used_for_finite_precision_magnitude_only','rowwise_K_reset_correction_domain_forbidden','reset_coefficient_family_requires_same_graph_correction_domain','same_history_correlation_retained_for_storage_test','independent_row_box_word_composition_forbidden'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_P_H_R_K_boxes_claimed_physical','Kalman_and_reset_coefficient_family_outwardly_bounded','consecutive_storage_closed_here','P4_MOTION_PASS','P4_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('reset_inverse_operator_norm_upper_for_any_finite_correction')!=1.0:f.append('reset inverse norm changed')
    for mode,m in d['modes'].items():
        x=float(m['independent_rowbox_attitude_correction_norm_upper_diagnostic'])
        if not(math.isfinite(x) and x>=0):f.append(mode+' diagnostic correction invalid')
        for g,row in m['rows'].items():
            for key in ('K_row_norm_upper_accelerometer','K_row_norm_upper_magnetometer','K_row_norm_upper_S_zero'):
                if not all(math.isfinite(float(x)) and float(x)>=0 for x in row[key]):f.append(mode+' '+g+' '+key+' invalid')
    return list(dict.fromkeys(f))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'Joseph_gain_family':d['Joseph_gain_family_outwardly_bounded'],'H18_rowbox_correction_diagnostic':d['modes']['H18']['independent_rowbox_attitude_correction_norm_upper_diagnostic'],'A21_rowbox_correction_diagnostic':d['modes']['A21']['independent_rowbox_attitude_correction_norm_upper_diagnostic'],'reset_domain_same_graph_required':d['reset_coefficient_family_requires_same_graph_correction_domain'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
