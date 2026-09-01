#!/usr/bin/env python3
from __future__ import annotations
import argparse,copy,json,math
from pathlib import Path
import ou3_p4_nonlinear_word_certificate as LEGACY
import ou3_p4_nextgen_widened_certificate as P4W
REPO=Path(__file__).resolve().parents[1]; DEFAULT_DOMAIN=REPO/'tools'/'ou3_proof_operating_domain.json'
def _sqrt_one_minus_up(x): return LEGACY.sqrt_up(math.nextafter(1.0-float(x),math.inf))
def _endpoint_sqrt_gap_lower(delta): return LEGACY.div_down(LEGACY.mul_down(.5,delta),LEGACY.add_up(_sqrt_one_minus_up(.5*delta),_sqrt_one_minus_up(delta)))
def _refine_mode(mode,base,dom):
 m=copy.deepcopy(base); live=dom['normal_live']; K=float(base['full_gain_norm_upper']); q=float(base['correction_quadratic_bound']['design_error_norm_radius']); f=float(live['specific_force_norm_upper_mps2']); mag=float(live['magnetic_vector_norm_upper_uT']); Hg=float(base['measurement_linear_operator_norm_upper']); H={'S_zero':1.0,'accelerometer':LEGACY.add_up(LEGACY.mul_up(2.,f),2.),'magnetometer':LEGACY.add_up(LEGACY.mul_up(2.,mag),2.)}
 if any(v>Hg for v in H.values()): raise RuntimeError(f'{mode}: directional H exceeds global')
 L={k:LEGACY.mul_up(K,v) for k,v in H.items()}; Cva=float(base['vector_residual_quadratic_constant_acc_upper']); Cvm=float(base['vector_residual_quadratic_constant_mag_upper']); Cia=LEGACY.mul_up(K,Cva); Cim=LEGACY.mul_up(K,Cvm)
 cs=LEGACY._composition_quadratic_constant(L['S_zero'],0.,q); ca=LEGACY._composition_quadratic_constant(L['accelerometer'],Cia,q); cm=LEGACY._composition_quadratic_constant(L['magnetometer'],Cim,q); C={'prediction':float(base['operation_specific_quadratic_defect_constants_upper']['prediction']),'S_zero_accepted':float(cs['full_state_quadratic_defect_constant_upper']),'accelerometer_accepted':float(ca['full_state_quadratic_defect_constant_upper']),'magnetometer_accepted':float(cm['full_state_quadratic_defect_constant_upper'])}; s=LEGACY.add_up(LEGACY.add_up(C['prediction'],C['S_zero_accepted']),LEGACY.add_up(C['accelerometer_accepted'],C['magnetometer_accepted'])); prev=float(base['operation_specific_defect_sum_per_sample_upper']);
 if s>prev: raise RuntimeError(f'{mode}: directional sum regressed')
 n=int(base['word_samples_upper']); mmin=float(base['metric_lambda_min_lower']); mmax=float(base['metric_lambda_max_upper']); B=LEGACY.div_up(LEGACY.mul_up(LEGACY.mul_up(LEGACY.mul_up(LEGACY.PREFIX_BOOTSTRAP_W_FACTOR,float(n)),LEGACY.sqrt_up(mmax)),s),mmin); Bprev=float(base['transported_word_defect_B_upper']);
 if B>Bprev: raise RuntimeError(f'{mode}: B regressed')
 delta=float(base['P3_word_endpoint_delta_lower']); gap=_endpoint_sqrt_gap_lower(delta); sw=min(LEGACY.div_down(gap,B),LEGACY.div_down(1.,B)); W=LEGACY.mul_down(sw,sw); Wprev=float(base['certified_level_W']);
 if W<Wprev: raise RuntimeError(f'{mode}: W regressed')
 qp=LEGACY.mul_up(2.,LEGACY.sqrt_up(LEGACY.div_up(W,mmin))); corr={'S_zero':LEGACY.mul_up(L['S_zero'],qp),'accelerometer':LEGACY.add_up(LEGACY.mul_up(L['accelerometer'],qp),LEGACY.mul_up(Cia,qp*qp)),'magnetometer':LEGACY.add_up(LEGACY.mul_up(L['magnetometer'],qp),LEGACY.mul_up(Cim,qp*qp))}; cp=max(corr.values())
 if not qp<q or not cp<1e-2: raise RuntimeError(f'{mode}: chart safety failed')
 proj=copy.deepcopy(base.get('active_bias_projection'))
 if mode=='A' and proj:
  margin=float(proj['interior_margin_lower_mps2']); proj['certified_error_norm_prefix_upper']=qp; proj['projection_surface_reached_in_certified_funnel']=not(qp<margin)
  if not qp<margin: raise RuntimeError('A: projection reached')
 m.update({'directional_measurement_operator_norm_upper':{'global_previous':Hg,**H},'directional_linear_correction_gain_upper':L,'directional_operation_quadratic_defect_constants_upper':C,'directional_operation_defect_sum_per_sample_upper':s,'operation_specific_defect_sum_per_sample_upper_previous':prev,'directional_defect_sum_monotone':True,'transported_word_defect_B_upper_previous':Bprev,'transported_word_defect_B_upper':B,'directional_B_reduction_factor_lower':LEGACY.div_down(Bprev,B),'exact_endpoint_sqrt_gap_lower':gap,'legacy_delta_over_8_budget_replaced':True,'prefix_bootstrap_B_sqrt_W_upper':LEGACY.mul_up(B,sw),'certified_level_W_previous_nextgen':Wprev,'certified_level_sqrt_W_previous_nextgen':float(base['certified_level_sqrt_W']),'certified_level_W':W,'certified_level_sqrt_W':sw,'secondgen_W_widening_factor_lower':LEGACY.div_down(W,Wprev),'total_W_widening_factor_vs_legacy_lower':LEGACY.div_down(W,float(base['certified_level_W_legacy'])),'prefix_canonical_error_norm_upper':qp,'accepted_correction_norm_prefix_upper':cp,'accepted_correction_norms_by_class_upper':corr,'active_bias_projection':proj,'nextgen_directional_operator_transport':True,'nextgen_exact_endpoint_budget':True,'exact_nonlinear_word_pass':True}); return m
def build(domain_path=DEFAULT_DOMAIN):
 p=Path(domain_path).resolve(); dom=json.loads(p.read_text()); prev=P4W.build(p); failures=[f'P4W: {x}' for x in P4W.validate(prev)]; modes={}
 if not failures:
  for mode in ('H','A'):
   try:modes[mode]=_refine_mode(mode,prev['modes'][mode],dom)
   except Exception as e: failures.append(f'{mode}: {e}')
 out=copy.deepcopy(prev); out['modes']=modes; out['nextgen_directional_operator_refinement']=True; out['nextgen_exact_endpoint_budget_refinement']=True; ok=not failures and len(modes)==2; out['P4_NEXTGEN_DIRECTIONAL_WORD_CERTIFICATE']='PASS' if ok else 'FAIL'; out['P4_NEXTGEN_WIDENED_WORD_CERTIFICATE']='PASS' if ok else 'FAIL'; out['P4_EXACT_NONLINEAR_WORD_CERTIFICATE']='PASS' if ok else 'FAIL'; out['failures']=failures; return out
def validate(d):
 f=list(d.get('failures',[]));
 for mode in ('H','A'):
  m=d.get('modes',{}).get(mode,{})
  if not m.get('nextgen_directional_operator_transport'): f.append(f'{mode}: missing directional')
  elif float(m['certified_level_W'])<float(m['certified_level_W_previous_nextgen']):f.append(f'{mode}: W regressed')
 return f
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args(); d=build(a.domain); f=validate(d); d['validation_pass']=not f; d['validation_failures']=f; a.output.write_text(json.dumps(d,indent=2,sort_keys=True)); return 0 if not f else 2
if __name__=='__main__': raise SystemExit(main())