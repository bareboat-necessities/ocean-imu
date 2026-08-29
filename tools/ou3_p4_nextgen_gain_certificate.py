#!/usr/bin/env python3
from __future__ import annotations
import argparse,copy,json,math
from pathlib import Path
import ou3_p4_nonlinear_word_certificate as LEGACY
import ou3_p4_nextgen_directional_certificate as P4D
REPO=Path(__file__).resolve().parents[1]; DEFAULT_DOMAIN=REPO/'tools'/'ou3_proof_operating_domain.json'
def _refine_mode(mode,base):
 m=copy.deepcopy(base); smax=float(base['Sigma_lambda_max_upper']); meas=base['measurement_bounds']; r={'S_zero':float(meas['S_zero_variance_lower']),'accelerometer':float(meas['acc_variance_lower']),'magnetometer':float(meas['mag_variance_lower'])}; Kg=float(base['full_gain_norm_upper']); K={k:LEGACY.sqrt_up(LEGACY.div_up(smax,v)) for k,v in r.items()}
 if any(v>Kg for v in K.values()): raise RuntimeError(f'{mode}: class-local K exceeds global')
 h=base['directional_measurement_operator_norm_upper']; L={k:LEGACY.mul_up(K[k],float(h[k])) for k in K}; q=float(base['correction_quadratic_bound']['design_error_norm_radius']); Cva=float(base['vector_residual_quadratic_constant_acc_upper']); Cvm=float(base['vector_residual_quadratic_constant_mag_upper']); Cia=LEGACY.mul_up(K['accelerometer'],Cva); Cim=LEGACY.mul_up(K['magnetometer'],Cvm)
 cs=LEGACY._composition_quadratic_constant(L['S_zero'],0.,q); ca=LEGACY._composition_quadratic_constant(L['accelerometer'],Cia,q); cm=LEGACY._composition_quadratic_constant(L['magnetometer'],Cim,q); C={'prediction':float(base['directional_operation_quadratic_defect_constants_upper']['prediction']),'S_zero_accepted':float(cs['full_state_quadratic_defect_constant_upper']),'accelerometer_accepted':float(ca['full_state_quadratic_defect_constant_upper']),'magnetometer_accepted':float(cm['full_state_quadratic_defect_constant_upper'])}; s=LEGACY.add_up(LEGACY.add_up(C['prediction'],C['S_zero_accepted']),LEGACY.add_up(C['accelerometer_accepted'],C['magnetometer_accepted'])); prev=float(base['directional_operation_defect_sum_per_sample_upper']);
 if s>prev: raise RuntimeError(f'{mode}: defect sum regressed')
 n=int(base['word_samples_upper']); mmin=float(base['metric_lambda_min_lower']); mmax=float(base['metric_lambda_max_upper']); B=LEGACY.div_up(LEGACY.mul_up(LEGACY.mul_up(LEGACY.mul_up(LEGACY.PREFIX_BOOTSTRAP_W_FACTOR,float(n)),LEGACY.sqrt_up(mmax)),s),mmin); Bprev=float(base['transported_word_defect_B_upper']);
 if B>Bprev: raise RuntimeError(f'{mode}: B regressed')
 gap=P4D._endpoint_sqrt_gap_lower(float(base['P3_word_endpoint_delta_lower'])); sw=min(LEGACY.div_down(gap,B),LEGACY.div_down(1.,B)); W=LEGACY.mul_down(sw,sw); Wprev=float(base['certified_level_W']); swprev=float(base['certified_level_sqrt_W']);
 if W<Wprev: raise RuntimeError(f'{mode}: W regressed')
 qp=LEGACY.mul_up(2.,LEGACY.sqrt_up(LEGACY.div_up(W,mmin)))
 if not qp<q:
  sw=swprev; W=Wprev; qp=float(base['prefix_canonical_error_norm_upper'])
 corr={'S_zero':LEGACY.mul_up(L['S_zero'],qp),'accelerometer':LEGACY.add_up(LEGACY.mul_up(L['accelerometer'],qp),LEGACY.mul_up(Cia,qp*qp)),'magnetometer':LEGACY.add_up(LEGACY.mul_up(L['magnetometer'],qp),LEGACY.mul_up(Cim,qp*qp))}; cp=max(corr.values())
 if not cp<1e-2: raise RuntimeError(f'{mode}: quaternion branch failed')
 proj=copy.deepcopy(base.get('active_bias_projection'))
 if mode=='A' and proj:
  margin=float(proj['interior_margin_lower_mps2']); proj['certified_error_norm_prefix_upper']=qp; proj['projection_surface_reached_in_certified_funnel']=not(qp<margin)
  if not qp<margin: raise RuntimeError('A: projection reached')
 m.update({'measurement_specific_R_lambda_min_lower':r,'measurement_specific_gain_norm_upper':K,'global_gain_norm_upper_previous':Kg,'measurement_specific_gain_bounds_monotone':True,'gain_refined_linear_correction_gain_upper':L,'gain_refined_operation_quadratic_defect_constants_upper':C,'gain_refined_operation_defect_sum_per_sample_upper':s,'directional_operation_defect_sum_per_sample_upper_previous':prev,'gain_refined_defect_sum_monotone':True,'transported_word_defect_B_upper_previous_gain_stage':Bprev,'transported_word_defect_B_upper':B,'gain_stage_B_reduction_factor_lower':LEGACY.div_down(Bprev,B),'certified_level_W_previous_gain_stage':Wprev,'certified_level_sqrt_W_previous_gain_stage':swprev,'certified_level_W':W,'certified_level_sqrt_W':sw,'gain_stage_W_widening_factor_lower':LEGACY.div_down(W,Wprev),'total_W_widening_factor_vs_legacy_lower':LEGACY.div_down(W,float(base['certified_level_W_legacy'])),'prefix_canonical_error_norm_upper':qp,'accepted_correction_norm_prefix_upper':cp,'accepted_correction_norms_by_class_upper':corr,'active_bias_projection':proj,'nextgen_measurement_specific_gain_transport':True,'exact_nonlinear_word_pass':True}); return m
def build(domain_path=DEFAULT_DOMAIN):
 p=Path(domain_path).resolve(); prev=P4D.build(p); failures=[f'P4D: {x}' for x in P4D.validate(prev)]; modes={}
 if not failures:
  for mode in ('H','A'):
   try:modes[mode]=_refine_mode(mode,prev['modes'][mode])
   except Exception as e: failures.append(f'{mode}: {e}')
 out=copy.deepcopy(prev); out['modes']=modes; out['nextgen_measurement_specific_gain_refinement']=True; ok=not failures and len(modes)==2; out['P4_NEXTGEN_GAIN_WORD_CERTIFICATE']='PASS' if ok else 'FAIL'; out['P4_NEXTGEN_DIRECTIONAL_WORD_CERTIFICATE']='PASS' if ok else 'FAIL'; out['P4_NEXTGEN_WIDENED_WORD_CERTIFICATE']='PASS' if ok else 'FAIL'; out['P4_EXACT_NONLINEAR_WORD_CERTIFICATE']='PASS' if ok else 'FAIL'; out['failures']=failures; return out
def validate(d):
 f=list(d.get('failures',[]));
 for mode in ('H','A'):
  m=d.get('modes',{}).get(mode,{})
  if not m.get('nextgen_measurement_specific_gain_transport'): f.append(f'{mode}: missing gain refinement')
  elif float(m['certified_level_W'])<float(m['certified_level_W_previous_gain_stage']): f.append(f'{mode}: W regressed')
 return f
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args(); d=build(a.domain); f=validate(d); d['validation_pass']=not f; d['validation_failures']=f; a.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({'status':d['P4_NEXTGEN_GAIN_WORD_CERTIFICATE'],'failures':f},indent=2)); return 0 if not f else 2
if __name__=='__main__': raise SystemExit(main())