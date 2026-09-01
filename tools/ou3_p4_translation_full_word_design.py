#!/usr/bin/env python3
"""Fast source-derived design probe for complete-word P4 translation dissipation.

This is deliberately NOT a theorem producer.  It reconstructs the exact P3
limiting translation cell directly, then uses ordinary numpy arithmetic to ask
whether longer complete words are worth validating.  No replay data are used.
The validated theorem route lives in ou3_p4_translation_full_word_interval.py.
"""
from __future__ import annotations
import argparse,itertools,json,math
from pathlib import Path
import numpy as np
from ou3_interval import Interval
import ou3_source_reachable_matrix_p3 as P3BASE
import ou3_p4_worst_translation_cell as WORST

REPO=Path(__file__).resolve().parents[1]; DEFAULT_DOMAIN=REPO/'tools'/'ou3_proof_operating_domain.json'
HORIZONS_S=(1.0,2.0,4.0,8.0)

def _mid(A): return np.asarray([[(x.lo+x.hi)*0.5 for x in r] for r in A],float)
def _transition(tau,h):
 x=h/tau; a=math.exp(-x); em1=math.expm1(-x); pva=-tau*em1
 if abs(x)<1e-2:
  x2=x*x;x3=x2*x;x4=x3*x;x5=x4*x; ppa=tau*tau*(.5*x2-x3/6+x4/24); psa=tau**3*(x3/6-x4/24+x5/120)
 else: ppa=tau*tau*(x+em1); psa=tau**3*(.5*x*x-x-em1)
 return np.array([[1,0,0,pva],[h,1,0,ppa],[.5*h*h,h,1,psa],[0,0,0,a]],float)
def _Q(tau,sigma,h):
 x=h/tau; q=_mid(P3BASE.qbar_integrated_ou(Interval.outward_bounds(x,x))); D=np.diag([sigma*tau,sigma*tau*tau,sigma*tau**3,sigma]); return .5*((D@q@D.T)+(D@q@D.T).T)
def _upd(P,i,R):
 c=P[:,i].copy(); d=float(P[i,i]+R); out=P-np.outer(c,c)/d; return .5*(out+out.T)
def _one(tau,sigma,rs,h,T,upper,ra):
 A=_transition(tau,h);Q=_Q(tau,sigma,h);P=np.zeros((4,4));n=max(1,int(math.ceil(T/h)))
 for _ in range(n): P=A@P@A.T+Q;P=_upd(P,2,rs*rs);P=_upd(P,3,ra)
 D=np.diag(1/np.sqrt(upper));G=D@P@D;delta=float(np.linalg.eigvalsh(.5*(G+G.T))[0])
 return {'tau_s':tau,'sigma_aw_mps2':sigma,'R_S_std':rs,'horizon_s':T,'steps':n,'translation_complete_word_generalized_margin_design':delta}
def _pts(a,b): return (a,math.sqrt(a*b),b)
def build(domain_path=DEFAULT_DOMAIN):
 p=Path(domain_path).resolve();modes={};fail=[]
 for mode in ('H','A'):
  try:
   c=WORST.build_cell(mode,p);s=WORST.serializable(c);h=float(c['sched']['dt_s']);x=c['x'];tb=(h/x.hi,h/x.lo);sb=c['sigma'].as_list();rb=c['rs'].as_list();u=list(map(float,c['row']['Sigma_diagonal_upper']));upper=np.array([u[6],u[9],u[12],u[15]],float);ra=float(c['vector']['configured_measurement_bounds']['acc_measurement_std_mps2'])**2;rows=[]
   for T in HORIZONS_S:
    for tau,sigma,rs in itertools.product(_pts(*tb),_pts(*sb),_pts(*rb)): rows.append(_one(tau,sigma,rs,h,T,upper,ra))
   hs={}
   for T in HORIZONS_S:
    rr=[r for r in rows if r['horizon_s']==T];w=min(rr,key=lambda r:r['translation_complete_word_generalized_margin_design']);b=max(rr,key=lambda r:r['translation_complete_word_generalized_margin_design']);old=float(c['row']['direct_translation_generalized_margin_lower']);hs[str(T)]={'worst_grid_point':w,'best_grid_point':b,'old_single_seed_translation_margin_lower':old,'design_worst_to_old_margin_ratio':w['translation_complete_word_generalized_margin_design']/old}
   modes[mode]={'source_cell':s,'tau_s_derived_from_x_cell':list(tb),'translation_directional_covariance_upper':upper.tolist(),'grid_points_per_horizon':27,'horizons':hs}
  except Exception as e: fail.append(f'{mode}: {e}')
 return {'qualification':'OU3_P4_TRANSLATION_FULL_WORD_DESIGN_PROBE','source_parameters_from_theorem_domain':True,'trajectory_replay_used':False,'ordinary_floating_point_design_only':True,'validated_for_theorem_promotion':False,'P4_USABLE_CERTIFICATE_STATUS':'NOT_ESTABLISHED','modes':modes,'failures':fail}
def validate(d):
 f=list(d.get('failures',[]))
 if d.get('trajectory_replay_used') is not False or d.get('ordinary_floating_point_design_only') is not True or d.get('validated_for_theorem_promotion') is not False:f.append('design qualification invalid')
 for mode in ('H','A'):
  for T in HORIZONS_S:
   q=d.get('modes',{}).get(mode,{}).get('horizons',{}).get(str(T),{}).get('worst_grid_point',{}).get('translation_complete_word_generalized_margin_design')
   if not isinstance(q,(int,float)) or not math.isfinite(q) or q<=0:f.append(f'{mode}: horizon {T} invalid')
 if d.get('P4_USABLE_CERTIFICATE_STATUS')!='NOT_ESTABLISHED':f.append('design probe promoted P4')
 return f
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build(a.domain);f=validate(d);d['validation_failures']=f;a.output.write_text(json.dumps(d,indent=2,sort_keys=True));print(json.dumps({'modes':{m:{h:d.get('modes',{}).get(m,{}).get('horizons',{}).get(h,{}).get('design_worst_to_old_margin_ratio') for h in map(str,HORIZONS_S)} for m in ('H','A')},'failures':f},indent=2));return 0 if not f else 2
if __name__=='__main__':raise SystemExit(main())
