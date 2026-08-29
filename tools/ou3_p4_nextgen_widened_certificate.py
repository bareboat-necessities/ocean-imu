#!/usr/bin/env python3
"""Next-generation widened OU-III P4 exact nonlinear word certificate.

Theorem-preserving operation-class refinement: prediction, S=0, accelerometer,
and magnetometer nonlinear defects are bounded separately before source-word
transport. Rejected/not-due branches remain covered by the per-sample upper
bound. No replay or source-language weakening is used.
"""
from __future__ import annotations
import argparse, copy, json, math
from pathlib import Path
import ou3_p4_nonlinear_word_certificate as LEGACY
REPO=Path(__file__).resolve().parents[1]; DEFAULT_DOMAIN=REPO/'tools'/'ou3_proof_operating_domain.json'

def _refine_mode(mode,base,domain):
 m=copy.deepcopy(base); live=domain['normal_live']; K=float(base['full_gain_norm_upper']); L=float(base['correction_quadratic_bound']['linear_correction_gain_L']); q=float(base['correction_quadratic_bound']['design_error_norm_radius']); f=float(live['specific_force_norm_upper_mps2']); mag=float(live['magnetic_vector_norm_upper_uT'])
 Cva=LEGACY.add_up(LEGACY.mul_up(LEGACY.ROTATION_REMAINDER_COEFF,f),1.5); Cvm=LEGACY.mul_up(LEGACY.ROTATION_REMAINDER_COEFF,mag); Cia=LEGACY.mul_up(K,Cva); Cim=LEGACY.mul_up(K,Cvm)
 cs=LEGACY._composition_quadratic_constant(L,0.0,q); ca=LEGACY._composition_quadratic_constant(L,Cia,q); cm=LEGACY._composition_quadratic_constant(L,Cim,q)
 cp=float(base['prediction_quadratic_bound']['full_state_quadratic_defect_constant_upper']); vals={'prediction':cp,'S_zero_accepted':float(cs['full_state_quadratic_defect_constant_upper']),'accelerometer_accepted':float(ca['full_state_quadratic_defect_constant_upper']),'magnetometer_accepted':float(cm['full_state_quadratic_defect_constant_upper'])}
 s=LEGACY.add_up(LEGACY.add_up(vals['prediction'],vals['S_zero_accepted']),LEGACY.add_up(vals['accelerometer_accepted'],vals['magnetometer_accepted'])); old=LEGACY.mul_up(4.0,float(base['uniform_operation_quadratic_defect_constant_upper']));
 if s>old: raise RuntimeError(f'{mode}: defect sum exceeds legacy budget')
 n=int(base['word_samples_upper']); mmin=float(base['metric_lambda_min_lower']); mmax=float(base['metric_lambda_max_upper']); delta=float(base['P3_word_endpoint_delta_lower']); B=LEGACY.div_up(LEGACY.mul_up(LEGACY.mul_up(LEGACY.mul_up(LEGACY.PREFIX_BOOTSTRAP_W_FACTOR,float(n)),LEGACY.sqrt_up(mmax)),s),mmin)
 oldB=float(base['transported_word_defect_B_upper']);
 if B>oldB: raise RuntimeError(f'{mode}: B regressed')
 sw=LEGACY.div_down(delta,LEGACY.mul_up(8.0,B)); W=LEGACY.mul_down(sw,sw); oldW=float(base['certified_level_W']);
 if W<oldW: raise RuntimeError(f'{mode}: W regressed')
 qp=LEGACY.mul_up(2.0,LEGACY.sqrt_up(LEGACY.div_up(W,mmin))); corr=LEGACY.add_up(LEGACY.mul_up(L,qp),LEGACY.mul_up(max(Cia,Cim),LEGACY.mul_up(qp,qp)))
 if not qp<q or not corr<1e-2: raise RuntimeError(f'{mode}: chart safety failed')
 proj=copy.deepcopy(base.get('active_bias_projection'))
 if mode=='A' and proj is not None:
  margin=float(proj['interior_margin_lower_mps2']); proj['certified_error_norm_prefix_upper']=qp; proj['projection_surface_reached_in_certified_funnel']=not(qp<margin)
  if not qp<margin: raise RuntimeError('A: projection reached')
 m.update({'vector_residual_quadratic_constant_acc_upper':Cva,'vector_residual_quadratic_constant_mag_upper':Cvm,'operation_specific_quadratic_defect_constants_upper':vals,'operation_specific_defect_sum_per_sample_upper':s,'legacy_four_operation_max_defect_per_sample_upper':old,'operation_specific_sum_no_larger_than_legacy_max_budget':True,'transported_word_defect_B_upper_legacy':oldB,'transported_word_defect_B_upper':B,'certified_level_W_legacy':oldW,'certified_level_sqrt_W_legacy':float(base['certified_level_sqrt_W']),'certified_level_W':W,'certified_level_sqrt_W':sw,'certified_level_W_widening_factor_lower':LEGACY.div_down(W,oldW),'certified_level_sqrt_W_widening_factor_lower':LEGACY.div_down(sw,float(base['certified_level_sqrt_W'])),'prefix_canonical_error_norm_upper':qp,'accepted_correction_norm_prefix_upper':corr,'active_bias_projection':proj,'nextgen_operation_specific_defect_transport':True,'exact_nonlinear_word_pass':True}); return m

def build(domain_path=DEFAULT_DOMAIN):
 p=Path(domain_path).resolve(); dom=json.loads(p.read_text()); legacy=LEGACY.build(p); failures=[f'legacy P4: {x}' for x in LEGACY.validate(legacy)]; modes={}
 if not failures:
  for mode in ('H','A'):
   try:modes[mode]=_refine_mode(mode,legacy['modes'][mode],dom)
   except Exception as e: failures.append(f'{mode}: {e}')
 out=copy.deepcopy(legacy); out['modes']=modes; out['nextgen_refinement_source_only']=True; out['legacy_P4_retained_as_baseline']=True; ok=not failures and len(modes)==2; out['P4_NEXTGEN_WIDENED_WORD_CERTIFICATE']='PASS' if ok else 'FAIL'; out['P4_EXACT_NONLINEAR_WORD_CERTIFICATE']='PASS' if ok else 'FAIL'; out['failures']=failures; return out

def validate(d):
 f=list(d.get('failures',[]));
 for mode in ('H','A'):
  m=d.get('modes',{}).get(mode,{})
  if not m.get('nextgen_operation_specific_defect_transport'): f.append(f'{mode}: missing refinement')
  elif float(m['certified_level_W'])<float(m['certified_level_W_legacy']): f.append(f'{mode}: W regressed')
 return f

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args(); d=build(a.domain); f=validate(d); d['validation_pass']=not f; d['validation_failures']=f; a.output.write_text(json.dumps(d,indent=2,sort_keys=True)); return 0 if not f else 2
if __name__=='__main__': raise SystemExit(main())