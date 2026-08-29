#!/usr/bin/env python3
"""Direct strict-contraction search for the validated nonlinear P4 word envelope.

Unlike the preceding P4 layers, this experiment does not insist on the fixed
advertised decrease delta/2.  For each source-safe proof radius q it asks the
minimal theorem question directly: is the exact-word enclosure strictly
contractive?

The homogeneous endpoint has sqrt(W) gain <=sqrt(1-delta).  The outward-rounded
nonlinear word defect at that radius is B(q) W.  Strict contraction therefore
follows from

  sqrt(1-delta) + B(q) sqrt(W) < 1.

The positive gap is evaluated cancellation-free as

  1-sqrt(1-delta) = delta/(1+sqrt(1-delta)).

We retain a small outward-safe interior factor below the strict boundary, check
the exact prefix/design/Cayley/projection/quaternion constraints, and select the
largest certified W over the proof-radius grid.  This is a direct return-map
contraction test of the validated exact-operation envelope; it is source-only
and is intentionally separate from the stronger fixed delta/2 decrease claim.
"""
from __future__ import annotations
import argparse,copy,json,math
from pathlib import Path
import ou3_p4_nonlinear_word_certificate as L
import ou3_p4_thirdgen_combined_certificate as T
import ou3_p4_exact_correction_structure_certificate as E
REPO=Path(__file__).resolve().parents[1];DEFAULT_DOMAIN=REPO/'tools'/'ou3_proof_operating_domain.json'
def strict_gap(delta):
 s=L.sqrt_up(math.nextafter(1.-float(delta),math.inf));return L.div_down(float(delta),L.add_up(1.,s))
def candidate(mode,base,domain,q):
 r=E._candidate(mode,base,domain,q);B=float(r['B']);gap=strict_gap(float(base['P3_word_endpoint_delta_lower'])); sm=L.GROUP.sqrt_point(float(base['metric_lambda_min_lower'])).lo
 caps={'strict_endpoint':L.mul_down(0.999999999999, L.div_down(gap,B)),'bootstrap':L.div_down(1.,B),'design_radius':T._positive_root(B,L.mul_down(q,sm)),'cayley_chart':T._positive_root(B,L.mul_down(float(base['cayley_norm_limit']),sm))}
 proj=base.get('active_bias_projection')
 if mode=='A' and proj:caps['bias_projection']=T._positive_root(B,L.mul_down(float(proj['interior_margin_lower_mps2']),sm))
 sw=min(caps.values());W=L.mul_down(sw,sw);pf=L.add_up(1.,L.mul_up(B,sw));qp=L.div_up(L.mul_up(pf,sw),sm)
 if not qp<=q or not qp<float(base['cayley_norm_limit']):raise RuntimeError('prefix')
 if not L.mul_up(B,sw)<gap:raise RuntimeError('not strict')
 return {'q':q,'B':B,'sqrtW':sw,'W':W,'qprefix':qp,'active_cap':min(caps,key=caps.get),'strict_gap':gap}
def refine(mode,base,domain):
 q0=float(base['thirdgen_selected_design_norm']);rows=[]
 for q in T._grid(q0):
  try:rows.append(candidate(mode,base,domain,q))
  except Exception:pass
 if not rows:raise RuntimeError('no direct contraction cell')
 best=max(rows,key=lambda x:x['W']);before=float(base['certified_level_W']);m=copy.deepcopy(base);m.update({'direct_strict_word_contraction_search':True,'direct_contraction_candidates':rows,'direct_selected_design_norm':best['q'],'direct_selected_active_cap':best['active_cap'],'direct_strict_endpoint_gap_lower':best['strict_gap'],'certified_level_W_before_direct':before,'certified_level_W':max(before,best['W']),'certified_level_sqrt_W':max(float(base['certified_level_sqrt_W']),best['sqrtW']),'prefix_canonical_error_norm_upper':best['qprefix'],'direct_W_factor_lower':L.div_down(max(before,best['W']),before),'direct_claim':'W_end<W0 with positive strict gap; fixed delta/2 decrease is not claimed at the enlarged level'});return m
def build(domain_path=DEFAULT_DOMAIN):
 p=Path(domain_path).resolve();domain=json.loads(p.read_text());base=E.build(p);f=list(E.validate(base));modes={}
 if not f:
  for mode in ('H','A'):
   try:modes[mode]=refine(mode,base['modes'][mode],domain)
   except Exception as x:f.append(f'{mode}: {x}')
 out=copy.deepcopy(base);out['modes']=modes;out['direct_word_contraction_source_only']=True;out['P4_DIRECT_STRICT_WORD_CONTRACTION_CERTIFICATE']='PASS' if not f and len(modes)==2 else 'FAIL';out['failures']=f;return out
def validate(d):
 f=list(d.get('failures',[]))
 for mode in ('H','A'):
  m=d.get('modes',{}).get(mode,{})
  if not m.get('direct_strict_word_contraction_search'):f.append(f'{mode}: missing direct search')
  elif float(m['certified_level_W'])<float(m['certified_level_W_before_direct']):f.append(f'{mode}: regressed')
 return f
def main():
 a=argparse.ArgumentParser();a.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN);a.add_argument('--output',type=Path,required=True);x=a.parse_args();d=build(x.domain);f=validate(d);d['validation_failures']=f;x.output.write_text(json.dumps(d,indent=2,sort_keys=True));print(json.dumps({'status':d['P4_DIRECT_STRICT_WORD_CONTRACTION_CERTIFICATE'],'modes':{m:{'W_before':d.get('modes',{}).get(m,{}).get('certified_level_W_before_direct'),'W_after':d.get('modes',{}).get(m,{}).get('certified_level_W'),'factor':d.get('modes',{}).get(m,{}).get('direct_W_factor_lower'),'q':d.get('modes',{}).get(m,{}).get('direct_selected_design_norm'),'cap':d.get('modes',{}).get(m,{}).get('direct_selected_active_cap')} for m in ('H','A')},'failures':f},indent=2));return 0 if not f else 2
if __name__=='__main__':raise SystemExit(main())