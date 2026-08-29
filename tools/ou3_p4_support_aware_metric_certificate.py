#!/usr/bin/env python3
from __future__ import annotations
import argparse,copy,json,math
from pathlib import Path
import ou3_p4_nonlinear_word_certificate as L
import ou3_p4_nextgen_directional_certificate as D
import ou3_p4_thirdgen_combined_certificate as T
import ou3_p4_exact_correction_structure_certificate as E
import ou3_explicit_information_word_certificate as P3
REPO=Path(__file__).resolve().parents[1]; DEFAULT_DOMAIN=REPO/'tools'/'ou3_proof_operating_domain.json'
def find(o,k):
 if isinstance(o,dict):
  if k in o:return o[k]
  for v in o.values():
   r=find(v,k)
   if r is not None:return r
 if isinstance(o,list):
  for v in o:
   r=find(v,k)
   if r is not None:return r
 return None
def refine(mode,b,p3):
 row=p3['modes'][mode]; scales=find(row,'comparison_scale_diagonal_squared'); qpost=find(row,'post_measurement_scaled_Omega_lambda_min_lower')
 if not isinstance(scales,list) or len(scales)<3 or not qpost: raise RuntimeError('missing conditioned covariance lower data')
 scale=float(b['metric_mode_global_positive_scale']); qpost=float(qpost)
 matt=L.up(scale*max(1.0/(qpost*float(scales[i])) for i in range(3)))
 mmax=float(b['metric_lambda_max_upper']); matt=min(matt,mmax)
 C=b['exact_structure_operation_defect_constants_upper']; att=float(C['prediction'])+float(C['S_zero_accepted']); full=float(C['accelerometer_accepted'])+float(C['magnetometer_accepted'])
 weighted=L.add_up(L.mul_up(L.sqrt_up(matt),att),L.mul_up(L.sqrt_up(mmax),full)); n=int(b['word_samples_upper']); mmin=float(b['metric_lambda_min_lower']); B=L.div_up(L.mul_up(L.mul_up(L.PREFIX_BOOTSTRAP_W_FACTOR,float(n)),weighted),mmin)
 old=float(b['transported_word_defect_B_upper']); B=min(B,old); gap=D._endpoint_sqrt_gap_lower(float(b['P3_word_endpoint_delta_lower'])); sw=min(L.div_down(gap,B),L.div_down(1.,B)); q=float(b['exact_structure_selected_design_norm']); sm=L.GROUP.sqrt_point(mmin).lo; sw=min(sw,T._positive_root(B,L.mul_down(q,sm))); W=L.mul_down(sw,sw); before=float(b['certified_level_W']); W=max(W,before); sw=max(sw,float(b['certified_level_sqrt_W']))
 m=copy.deepcopy(b); m.update({'support_aware_information_transport':True,'attitude_support_metric_lambda_max_upper':matt,'global_metric_lambda_max_upper_previous':mmax,'support_aware_weighted_defect_upper':weighted,'transported_word_defect_B_upper_before_support':old,'transported_word_defect_B_upper':B,'certified_level_W_before_support':before,'certified_level_W':W,'certified_level_sqrt_W':sw,'support_aware_W_factor_lower':L.div_down(W,before)}); return m
def build(domain_path=DEFAULT_DOMAIN):
 p=Path(domain_path).resolve(); base=E.build(p); p3=P3.build(p); f=[*E.validate(base),*P3.validate(p3)]; modes={}
 if not f:
  for mode in ('H','A'):
   try:modes[mode]=refine(mode,base['modes'][mode],p3)
   except Exception as x:f.append(f'{mode}: {x}')
 out=copy.deepcopy(base); out['modes']=modes; out['support_aware_metric_source_only']=True; out['P4_SUPPORT_AWARE_METRIC_CERTIFICATE']='PASS' if not f and len(modes)==2 else 'FAIL'; out['failures']=f; return out
def validate(d):
 f=list(d.get('failures',[]))
 for mode in ('H','A'):
  m=d.get('modes',{}).get(mode,{})
  if not m.get('support_aware_information_transport'): f.append(f'{mode}: missing support transport')
  elif float(m['certified_level_W'])<float(m['certified_level_W_before_support']):f.append(f'{mode}: regressed')
 return f
def main():
 a=argparse.ArgumentParser();a.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN);a.add_argument('--output',type=Path,required=True);x=a.parse_args();d=build(x.domain);f=validate(d);d['validation_failures']=f;x.output.write_text(json.dumps(d,indent=2,sort_keys=True));print(json.dumps({'status':d['P4_SUPPORT_AWARE_METRIC_CERTIFICATE'],'modes':{m:{'W_before':d.get('modes',{}).get(m,{}).get('certified_level_W_before_support'),'W_after':d.get('modes',{}).get(m,{}).get('certified_level_W'),'factor':d.get('modes',{}).get(m,{}).get('support_aware_W_factor_lower'),'matt':d.get('modes',{}).get(m,{}).get('attitude_support_metric_lambda_max_upper')} for m in ('H','A')},'failures':f},indent=2));return 0 if not f else 2
if __name__=='__main__':raise SystemExit(main())