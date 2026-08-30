#!/usr/bin/env python3
"""Validated complete-word translation dissipation on the old limiting P3 cell.

The old P4 radius is destroyed by a translation margin obtained from one
covariance seed. Here we accumulate a certified covariance component through a
whole configurable word instead of discarding every later process contribution.

Use the P3-conditioned coordinates
 D=diag(sigma_min*h, sigma_min*h^2, sigma_min*h^3, sigma_min).
Every prediction injects at least rho_trans I in these coordinates.

For corrections, do not form a nearly singular artificial posterior. If L is a
lower covariance component and J=H'R^-1H is the maximum admissible measurement
information on this translation axis, then J <= beta I and

 (L^-1+J)^-1 >= (L^-1+beta I)^-1
                  >= L/(1+beta*lambda_max(L)).

Thus every possible S/accelerometer correction can be covered by a scalar
Loewner shrink of L.  For the orthogonal S and a_w coordinates,
beta=max(1/R_S_norm,1/R_acc_norm), not their sum.  Rejected/not-due updates are
less informative and are automatically covered.  This keeps a robust SPD lower
matrix while still accounting for the worst measurement information at every
sample.

All promoted arithmetic is outward-rounded. Tau is adaptively subdivided when
interval dependency prevents a positive Loewner lower. This file proves only
the previously limiting translation source cell and cannot by itself promote
complete P4.
"""
from __future__ import annotations
import argparse,json,math,re
from pathlib import Path
from ou3_interval import Interval,matrix_mul,matrix_transpose,symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_validated_transcendentals as VT
import ou3_p4_worst_translation_cell as WORST

REPO=Path(__file__).resolve().parents[1]; DEFAULT_DOMAIN=REPO/'tools'/'ou3_proof_operating_domain.json'; WRAPPER=REPO/'src'/'kalman_ou_iii'/'SeaStateFusionFilter_OU_III.h'
DEFAULT_HORIZON_S=1.0; MAX_TAU_SPLIT_DEPTH=14

def down(x):return math.nextafter(float(x),-math.inf)
def up(x):return math.nextafter(float(x),math.inf)
def I(x):return Interval.outward_bounds(float(x),float(x))
def _pm(A):return [[I(float(x)) for x in r] for r in A]

def _center_radius_lower(A):
 A=matrix_symmetric_hull(A);n=len(A);C=[[0.0]*n for _ in range(n)];rows=[0.0]*n
 for i in range(n):
  for j in range(n):
   a=A[i][j];c=min(max(.5*a.lo+.5*a.hi,a.lo),a.hi);C[i][j]=c;r=up(max(abs(a.lo-c),abs(a.hi-c)));rows[i]=up(rows[i]+r)
 rad=max(rows);L=[r[:] for r in C]
 for i in range(n):L[i][i]=down(L[i][i]-rad)
 if not symmetric_positive_definite_ldlt(_pm(L))[0]:raise RuntimeError(f'Loewner lower extraction lost SPD (conditioned radius={rad:.3e})')
 return L,rad

def _F(tau,h):
 x=I(h)/tau
 if x.hi>=1e-2:raise RuntimeError('limiting cell left shipping small-x branch')
 a=VT.exp_interval(-x);em1=VT.expm1_interval(-x);pva=-(tau*em1);x2=x*x;x3=x2*x;x4=x3*x;x5=x4*x
 ppa=tau*tau*(I(.5)*x2-I(1/6)*x3+I(1/24)*x4);psa=tau*tau*tau*(I(1/6)*x3-I(1/24)*x4+I(1/120)*x5);z=Interval.point(0);o=Interval.point(1)
 return [[o,z,z,pva/I(h)],[o,o,z,ppa/I(h*h)],[I(.5),o,o,psa/I(h*h*h)],[z,z,z,a]]

def _predict(L,F,rho):
 M=matrix_mul(matrix_mul(F,_pm(L)),matrix_transpose(F))
 for i in range(4):M[i][i]=M[i][i]+I(rho)
 return _center_radius_lower(M)

def _row_norm_upper(L):
 best=0.0
 for r in L:
  s=0.0
  for x in r:s=up(s+abs(float(x)))
  best=max(best,s)
 return up(best)

def _measurement_information_shrink(L,beta):
 lam=_row_norm_upper(L)
 f=(I(1)/(I(1)+I(beta)*I(lam))).lo
 if not 0<f<=1:raise RuntimeError('invalid correction-information shrink')
 M=[[I(f)*I(L[i][j]) for j in range(4)] for i in range(4)]
 return _center_radius_lower(M)

def _spd_delta(L,upper,d):
 A=_pm(L);q=I(d)
 for i in range(4):A[i][i]=A[i][i]-q*I(upper[i])
 return symmetric_positive_definite_ldlt(matrix_symmetric_hull(A))[0]
def _delta(L,upper):
 lo=0.0;trial=1e-36
 while trial<1 and _spd_delta(L,upper,trial):lo=trial;trial*=10
 if lo==0:return 0.0
 hi=min(1.0,trial)
 for _ in range(64):
  mid=math.sqrt(lo*hi)
  if _spd_delta(L,upper,mid):lo=mid
  else:hi=mid
 return down(lo)
def _member(t,n):
 m=re.search(rf'float\s+{re.escape(n)}\s*=\s*([0-9.eE+-]+)f\s*;',t)
 if not m:raise RuntimeError(f'cannot source-bind {n}')
 return float(m.group(1))

def _prop(tau,h,n,rho,beta,depth=0):
 try:
  F=_F(tau,h);L=[[0.0]*4 for _ in range(4)];mr=0.0
  for k in range(n):
   if k==0:L=[[rho if i==j else 0.0 for j in range(4)] for i in range(4)]
   else:L,r=_predict(L,F,rho);mr=max(mr,r)
   L,r=_measurement_information_shrink(L,beta);mr=max(mr,r)
  return [(tau,L,mr,depth)]
 except RuntimeError:
  if depth>=MAX_TAU_SPLIT_DEPTH:raise
  mid=math.sqrt(tau.lo*tau.hi);a=Interval.outward_bounds(tau.lo,mid);b=Interval.outward_bounds(mid,tau.hi)
  return _prop(a,h,n,rho,beta,depth+1)+_prop(b,h,n,rho,beta,depth+1)

def _mode(mode,p,horizon_s):
 c=WORST.build_cell(mode,p);s=WORST.serializable(c);row=c['row'];h=float(c['sched']['dt_s']);x=c['x'];tau=Interval.outward_bounds(h/x.hi,h/x.lo);sigma=float(c['sigma'].lo);rho=float(c['rho_translation_lower'])
 scale2=[(sigma*h)**2,(sigma*h*h)**2,(sigma*h*h*h)**2,sigma*sigma];u=list(map(float,row['Sigma_diagonal_upper']));physical=[u[6],u[9],u[12],u[15]];upper=[(I(physical[i])/I(scale2[i])).hi for i in range(4)]
 text=WRAPPER.read_text();rh=min(_member(text,'R_S_x_factor_'),_member(text,'R_S_y_factor_'),1.0);rs=(I(rh)*I(float(c['rs'].lo))).lo;rsvar=I(rs).square().lo;acc=float(c['vector']['configured_measurement_bounds']['acc_measurement_std_mps2']);accvar=I(acc).square().lo
 rS=(I(rsvar)/I(scale2[2])).lo;rA=(I(accvar)/I(scale2[3])).lo;invS=(I(1)/I(rS)).hi;invA=(I(1)/I(rA)).hi;beta=up(max(invS,invA));n=int(math.ceil(horizon_s/h));leaves=_prop(tau,h,n,rho,beta);cert=[]
 for t,L,rad,dep in leaves:
  d=_delta(L,upper)
  if d<=0 or not _spd_delta(L,upper,d):raise RuntimeError(f'nonpositive endpoint margin on tau leaf {t.as_list()}')
  cert.append({'tau_s':t.as_list(),'delta_lower':d,'max_conditioned_radius_removed':rad,'split_depth':dep})
 w=min(cert,key=lambda q:q['delta_lower']);old=float(row['direct_translation_generalized_margin_lower'])
 return {'source_cell':s,'conditioned_coordinates':'D^-1[v,p,S,a_w]','tau_interval_s':tau.as_list(),'tau_leaf_count':len(cert),'max_tau_split_depth_used':max(q['split_depth'] for q in cert),'steps':n,'horizon_s':horizon_s,'process_injection_lower_conditioned':rho,'maximum_measurement_information_beta_conditioned':beta,'artificial_S_variance_conditioned':rS,'artificial_acc_aw_variance_conditioned':rA,'translation_covariance_upper_conditioned':upper,'tau_leaf_certificates':cert,'complete_word_translation_margin_lower':w['delta_lower'],'limiting_tau_leaf':w,'old_single_seed_translation_margin_lower':old,'margin_widening_factor_lower':down(w['delta_lower']/old),'interval_ldlt_endpoint_recertified':True}
def build(domain_path=DEFAULT_DOMAIN,horizon_s=DEFAULT_HORIZON_S):
 horizon_s=float(horizon_s)
 if not math.isfinite(horizon_s) or horizon_s<=0:raise ValueError('horizon_s must be finite positive')
 p=Path(domain_path).resolve();m={};f=[]
 for mode in ('H','A'):
  try:m[mode]=_mode(mode,p,horizon_s)
  except Exception as e:f.append(f'{mode}: {e}')
 return {'qualification':'OU3_P4_VALIDATED_WORST_CELL_COMPLETE_WORD_TRANSLATION_DISSIPATION','source_only':True,'trajectory_replay_used':False,'outward_rounded':True,'horizon_s':horizon_s,'modes':m,'P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS':'PASS' if not f and len(m)==2 else 'NOT_ESTABLISHED','P4_USABLE_CERTIFICATE_STATUS':'NOT_ESTABLISHED','remaining_obligation':'extend complete-word propagation to every reachable source cell/edge and attitude-bias blocks, then validate exact nonlinear return map','failures':f}
def validate(d):
 f=list(d.get('failures',[]))
 if d.get('source_only') is not True or d.get('trajectory_replay_used') is not False or d.get('outward_rounded') is not True:f.append('qualification flags invalid')
 if not float(d.get('horizon_s',0))>0:f.append('invalid horizon')
 for mode in ('H','A'):
  m=d.get('modes',{}).get(mode,{})
  if not float(m.get('complete_word_translation_margin_lower',0))>0:f.append(f'{mode}: no complete-word translation margin')
  if m.get('interval_ldlt_endpoint_recertified') is not True:f.append(f'{mode}: endpoint not recertified')
  if not float(m.get('margin_widening_factor_lower',0))>1:f.append(f'{mode}: complete word did not widen seed')
 if d.get('P4_USABLE_CERTIFICATE_STATUS')!='NOT_ESTABLISHED':f.append('partial result prematurely promoted P4')
 return f
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN);ap.add_argument('--horizon-s',type=float,default=DEFAULT_HORIZON_S);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build(a.domain,a.horizon_s);f=validate(d);d['validation_failures']=f;a.output.write_text(json.dumps(d,indent=2,sort_keys=True));print(json.dumps({'translation_status':d['P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS'],'horizon_s':d['horizon_s'],'modes':{x:{'delta':d.get('modes',{}).get(x,{}).get('complete_word_translation_margin_lower'),'factor':d.get('modes',{}).get(x,{}).get('margin_widening_factor_lower'),'tau_leaves':d.get('modes',{}).get(x,{}).get('tau_leaf_count')} for x in ('H','A')},'failures':f},indent=2));return 0 if not f else 2
if __name__=='__main__':raise SystemExit(main())
