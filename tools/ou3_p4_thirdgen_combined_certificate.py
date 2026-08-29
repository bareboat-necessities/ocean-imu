#!/usr/bin/env python3
"""Third-generation combined OU-III P4 widening certificate.

Stacks the successful #431 operation-specific, directional-H, exact endpoint,
and measurement-specific gain refinements with the successful #432 proof-radius
continuation and exact prefix bootstrap.  Every q candidate rebuilds the
measurement-class Cayley/quaternion defect constants using the class-local gain
bounds, then solves the exact prefix-design inequality.  No replay, trajectory
sampling, source-language weakening, or filter change is used.
"""
from __future__ import annotations
import argparse,copy,json,math
from pathlib import Path
import ou3_p4_nonlinear_word_certificate as LEGACY
import ou3_p4_nextgen_gain_certificate as P4G
import ou3_p4_nextgen_directional_certificate as P4D
REPO=Path(__file__).resolve().parents[1]; DEFAULT_DOMAIN=REPO/'tools'/'ou3_proof_operating_domain.json'
def sqrt_down(x): return LEGACY.GROUP.sqrt_point(float(x)).lo
def _positive_root(B,a):
 disc=LEGACY.add_up(1.,LEGACY.mul_up(4.,LEGACY.mul_up(B,a))); return LEGACY.div_down(LEGACY.mul_down(2.,a),LEGACY.add_up(1.,LEGACY.sqrt_up(disc)))
def _grid(q0):
 vals={float(q0)}; q=max(float(q0)/64.,1e-12)
 for _ in range(30): vals.add(q); q*=2.
 return sorted(vals)
def _candidate(mode,base,q):
 mmin=float(base['metric_lambda_min_lower']); mmax=float(base['metric_lambda_max_upper']); n=int(base['word_samples_upper']); gap=P4D._endpoint_sqrt_gap_lower(float(base['P3_word_endpoint_delta_lower'])); K=base['measurement_specific_gain_norm_upper']; H=base['directional_measurement_operator_norm_upper']; Cva=float(base['vector_residual_quadratic_constant_acc_upper']); Cvm=float(base['vector_residual_quadratic_constant_mag_upper']); L={'S_zero':LEGACY.mul_up(float(K['S_zero']),float(H['S_zero'])),'accelerometer':LEGACY.mul_up(float(K['accelerometer']),float(H['accelerometer'])),'magnetometer':LEGACY.mul_up(float(K['magnetometer']),float(H['magnetometer']))}; Cia=LEGACY.mul_up(float(K['accelerometer']),Cva); Cim=LEGACY.mul_up(float(K['magnetometer']),Cvm)
 cs=LEGACY._composition_quadratic_constant(L['S_zero'],0.,q); ca=LEGACY._composition_quadratic_constant(L['accelerometer'],Cia,q); cm=LEGACY._composition_quadratic_constant(L['magnetometer'],Cim,q)
 dt=0.005; pred0=base['prediction_quadratic_bound']; Lpred=float(pred0['linear_correction_gain_L']); pred=LEGACY._composition_quadratic_constant(Lpred,0.,q)
 C={'prediction':float(pred['full_state_quadratic_defect_constant_upper']),'S_zero_accepted':float(cs['full_state_quadratic_defect_constant_upper']),'accelerometer_accepted':float(ca['full_state_quadratic_defect_constant_upper']),'magnetometer_accepted':float(cm['full_state_quadratic_defect_constant_upper'])}; s=LEGACY.add_up(LEGACY.add_up(C['prediction'],C['S_zero_accepted']),LEGACY.add_up(C['accelerometer_accepted'],C['magnetometer_accepted'])); B=LEGACY.div_up(LEGACY.mul_up(LEGACY.mul_up(LEGACY.mul_up(LEGACY.PREFIX_BOOTSTRAP_W_FACTOR,float(n)),LEGACY.sqrt_up(mmax)),s),mmin)
 if not math.isfinite(B) or B<=0: raise RuntimeError('invalid B')
 sm=sqrt_down(mmin); caps={'endpoint':LEGACY.div_down(gap,B),'bootstrap':LEGACY.div_down(1.,B),'design_radius':_positive_root(B,LEGACY.mul_down(q,sm)),'cayley_chart':_positive_root(B,LEGACY.mul_down(float(base['cayley_norm_limit']),sm))}
 proj=copy.deepcopy(base.get('active_bias_projection'))
 if mode=='A' and proj: caps['bias_projection']=_positive_root(B,LEGACY.mul_down(float(proj['interior_margin_lower_mps2']),sm))
 sw=min(caps.values()); W=LEGACY.mul_down(sw,sw); pf=LEGACY.add_up(1.,LEGACY.mul_up(B,sw)); qp=LEGACY.div_up(LEGACY.mul_up(pf,sw),sm)
 if not qp<=q or not qp<float(base['cayley_norm_limit']): raise RuntimeError('prefix outside design/chart')
 corr={'S_zero':LEGACY.mul_up(L['S_zero'],qp),'accelerometer':LEGACY.add_up(LEGACY.mul_up(L['accelerometer'],qp),LEGACY.mul_up(Cia,qp*qp)),'magnetometer':LEGACY.add_up(LEGACY.mul_up(L['magnetometer'],qp),LEGACY.mul_up(Cim,qp*qp))}; cp=max(corr.values())
 if not cp<1e-2: raise RuntimeError('quaternion branch')
 if mode=='A' and proj:
  margin=float(proj['interior_margin_lower_mps2']); proj['certified_error_norm_prefix_upper']=qp; proj['projection_surface_reached_in_certified_funnel']=not(qp<margin)
  if not qp<margin: raise RuntimeError('projection')
 if LEGACY.mul_up(B,sw)>gap or LEGACY.mul_up(B,sw)>1.: raise RuntimeError('endpoint/bootstrap')
 return {'q':q,'B':B,'W':W,'sqrtW':sw,'prefix_factor':pf,'qprefix':qp,'caps':caps,'active_cap':min(caps,key=caps.get),'C':C,'sum':s,'corr':corr,'corr_prefix':cp,'projection':proj}
def _refine_mode(mode,base):
 q0=float(base['correction_quadratic_bound']['design_error_norm_radius']); rows=[]
 for q in _grid(q0):
  try: rows.append(_candidate(mode,base,q))
  except Exception: pass
 if not rows: raise RuntimeError(f'{mode}: no combined radius candidate')
 best=max(rows,key=lambda x:x['W']); before=float(base['certified_level_W'])
 if best['W']<before: raise RuntimeError(f'{mode}: combined route regressed #431')
 m=copy.deepcopy(base); m.update({'thirdgen_combined_radius_continuation':True,'thirdgen_candidate_count_certified':len(rows),'thirdgen_selected_design_norm':best['q'],'thirdgen_selected_active_cap':best['active_cap'],'thirdgen_exact_prefix_factor_upper':best['prefix_factor'],'thirdgen_operation_defect_constants_upper':best['C'],'thirdgen_operation_defect_sum_upper':best['sum'],'transported_word_defect_B_upper_before_combined':float(base['transported_word_defect_B_upper']),'transported_word_defect_B_upper':best['B'],'certified_level_W_before_combined':before,'certified_level_sqrt_W_before_combined':float(base['certified_level_sqrt_W']),'certified_level_W':best['W'],'certified_level_sqrt_W':best['sqrtW'],'thirdgen_W_widening_factor_vs_431_lower':LEGACY.div_down(best['W'],before),'thirdgen_total_W_widening_factor_vs_legacy_lower':LEGACY.div_down(best['W'],float(base['certified_level_W_legacy'])),'prefix_W_factor_upper':LEGACY.mul_up(best['prefix_factor'],best['prefix_factor']),'prefix_canonical_error_norm_upper':best['qprefix'],'accepted_correction_norm_prefix_upper':best['corr_prefix'],'accepted_correction_norms_by_class_upper':best['corr'],'active_bias_projection':best['projection'],'thirdgen_candidates':[{'q':r['q'],'W':r['W'],'B':r['B'],'qprefix':r['qprefix'],'active_cap':r['active_cap']} for r in rows],'exact_nonlinear_word_pass':True}); return m
def build(domain_path=DEFAULT_DOMAIN):
 p=Path(domain_path).resolve(); prev=P4G.build(p); failures=[f'#431 gain P4: {x}' for x in P4G.validate(prev)]; modes={}
 if not failures:
  for mode in ('H','A'):
   try:modes[mode]=_refine_mode(mode,prev['modes'][mode])
   except Exception as e: failures.append(f'{mode}: {e}')
 out=copy.deepcopy(prev); out['qualification']='VALIDATED_THIRDGEN_COMBINED_RADIUS_GAIN_CAYLEY_SOURCE_WORD_CERTIFICATE'; out['claim']='P4_THIRDGEN_COMBINED_WIDENED_H_A_WORD_DISSIPATION'; out['modes']=modes; out['thirdgen_stacks_431_and_432_refinements']=True; out['thirdgen_source_only']=True; ok=not failures and len(modes)==2; out['P4_THIRDGEN_COMBINED_WORD_CERTIFICATE']='PASS' if ok else 'FAIL'; out['P4_EXACT_NONLINEAR_WORD_CERTIFICATE']='PASS' if ok else 'FAIL'; out['failures']=failures; return out
def validate(d):
 f=list(d.get('failures',[]));
 if not d.get('thirdgen_stacks_431_and_432_refinements'): f.append('combined stack marker missing')
 for mode in ('H','A'):
  m=d.get('modes',{}).get(mode,{})
  if not m.get('thirdgen_combined_radius_continuation'): f.append(f'{mode}: combined radius missing'); continue
  if float(m['certified_level_W'])<float(m['certified_level_W_before_combined']): f.append(f'{mode}: regressed #431')
  if not float(m['accepted_correction_norm_prefix_upper'])<1e-2: f.append(f'{mode}: quaternion safety')
  if mode=='A' and m['active_bias_projection']['projection_surface_reached_in_certified_funnel'] is not False: f.append('A: projection')
 if not f and d.get('P4_THIRDGEN_COMBINED_WORD_CERTIFICATE')!='PASS': f.append('status not PASS')
 return f
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args(); d=build(a.domain); f=validate(d); d['validation_pass']=not f; d['validation_failures']=f; a.output.write_text(json.dumps(d,indent=2,sort_keys=True)); print(json.dumps({'status':d['P4_THIRDGEN_COMBINED_WORD_CERTIFICATE'],'numerical':{m:{'W_legacy':d.get('modes',{}).get(m,{}).get('certified_level_W_legacy'),'W_431':d.get('modes',{}).get(m,{}).get('certified_level_W_before_combined'),'W_combined':d.get('modes',{}).get(m,{}).get('certified_level_W'),'factor_vs_431':d.get('modes',{}).get(m,{}).get('thirdgen_W_widening_factor_vs_431_lower'),'factor_vs_legacy':d.get('modes',{}).get(m,{}).get('thirdgen_total_W_widening_factor_vs_legacy_lower'),'q':d.get('modes',{}).get(m,{}).get('thirdgen_selected_design_norm'),'active_cap':d.get('modes',{}).get(m,{}).get('thirdgen_selected_active_cap')} for m in ('H','A')},'failures':f},indent=2)); return 0 if not f else 2
if __name__=='__main__': raise SystemExit(main())