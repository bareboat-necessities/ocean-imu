#!/usr/bin/env python3
"""Validated complete-word translation dissipation on the old limiting P3 cell.

The old P4 radius is destroyed by a translation margin obtained from one
covariance seed. Here we accumulate a certified covariance component through a
whole configurable word instead of discarding every later process contribution.

Use the P3-conditioned coordinates
 D=diag(sigma_min*h, sigma_min*h^2, sigma_min*h^3, sigma_min).
Every prediction injects at least rho_trans I in these coordinates.

Corrections retain their actual translation directions.  If ``L`` is a certified
covariance lower and a possible scalar measurement contributes information
``J <= beta e_q e_q'``, monotonicity of the information update gives

    posterior(P,J) >= posterior(L,beta e_q e_q')
                    = L - beta L e_q e_q' L /(1+beta e_q' L e_q).

For the rank-one correction the inputs are deterministic binary64 certificate
numbers.  The Woodbury expression is therefore evaluated exactly as rational
arithmetic first, then converted to binary64 with one rigorous matrix rounding
allowance.  This avoids the catastrophic interval cancellation that occurs when
``L - beta*L[:,q]*L[q,:]/den`` is evaluated entry by entry with wide intervals.
The conversion subtracts an infinity-norm bound for the exact rational-to-float
rounding residual from the diagonal, which preserves a certified Loewner lower.

Prediction also preserves matrix structure.  Writing the interval transition as
``F = Fc + E`` gives

    F L F' = Fc L Fc' + Fc L E' + E L Fc' + E L E'.

Because ``L`` is positive definite, the final term is PSD and may be discarded
for a lower bound.  The symmetric cross term is enclosed by an outward-rounded
infinity-norm bound ``gamma`` assembled from elementwise radii of ``E`` and an
interval enclosure of ``Fc L``.  Hence

    F L F' >= Fc L Fc' - gamma I.

This avoids the severe dependency inflation of multiplying the full interval
matrix ``F L F'`` after directional measurement updates while remaining a
source-uniform Loewner lower.

Thus the S=0 pseudo can only remove information in the S direction and the
translation part of accelerometer corrections can only remove information in
the a_w direction.  We still apply *both* possible corrections at every IMU
sample, which is more informative than shipping scheduling/rejection and
therefore conservative.  This removes the previous artificial isotropic shrink
without relying on favorable timing.

All promoted arithmetic is outward-rounded. Tau is adaptively subdivided when
interval dependency prevents a positive Loewner lower. This file proves only
the previously limiting translation source cell and cannot by itself promote
complete P4.
"""
from __future__ import annotations
import argparse,json,math,re
from fractions import Fraction
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

def _down_fraction(q):
 f=float(q)
 if not math.isfinite(f):raise OverflowError('exact rational does not fit binary64')
 if Fraction.from_float(f)>q:f=math.nextafter(f,-math.inf)
 return f

def _up_fraction(q):
 f=float(q)
 if not math.isfinite(f):raise OverflowError('exact rational does not fit binary64')
 if Fraction.from_float(f)<q:f=math.nextafter(f,math.inf)
 return f

def _exact_rational_loewner_lower(Q):
 """Convert an exact symmetric rational matrix to a deterministic Loewner lower."""
 n=len(Q);C=[[float(Q[i][j]) for j in range(n)] for i in range(n)];rows=[]
 for i in range(n):
  s=Fraction(0)
  for j in range(n):s+=abs(Q[i][j]-Fraction.from_float(C[i][j]))
  rows.append(s)
 rad=_up_fraction(max(rows,default=Fraction(0)))
 L=[r[:] for r in C]
 for i in range(n):L[i][i]=down(L[i][i]-rad)
 if not symmetric_positive_definite_ldlt(_pm(L))[0]:raise RuntimeError(f'exact-rational Loewner lower lost SPD (rounding radius={rad:.3e})')
 return L,rad

def _F(tau,h):
 x=I(h)/tau
 if x.hi>=1e-2:raise RuntimeError('limiting cell left shipping small-x branch')
 a=VT.exp_interval(-x);em1=VT.expm1_interval(-x);pva=-(tau*em1);x2=x*x;x3=x2*x;x4=x3*x;x5=x4*x
 ppa=tau*tau*(I(.5)*x2-I(1/6)*x3+I(1/24)*x4);psa=tau*tau*tau*(I(1/6)*x3-I(1/24)*x4+I(1/120)*x5);z=Interval.point(0);o=Interval.point(1)
 return [[o,z,z,pva/I(h)],[o,o,z,ppa/I(h*h)],[I(.5),o,o,psa/I(h*h*h)],[z,z,z,a]]

def _midrad(A):
 C=[];R=[]
 for row in A:
  cr=[];rr=[]
  for a in row:
   c=min(max(.5*a.lo+.5*a.hi,a.lo),a.hi)
   r=up(max(abs(a.lo-c),abs(a.hi-c)))
   cr.append(c);rr.append(r)
  C.append(cr);R.append(rr)
 return C,R

def _abs_upper(a):return up(max(abs(a.lo),abs(a.hi)))

def _predict(L,F,rho):
 """Structured Loewner lower for F L F' + rho I."""
 Fc,R=_midrad(F)
 B=matrix_mul(_pm(Fc),_pm(L))
 # Cross term C = (Fc L)E' + E(L Fc').  It is symmetric.  Bound
 # ||C||_2 <= ||C||_inf by an outward-rounded row-sum enclosure.
 gamma=0.0
 for i in range(4):
  rowsum=0.0
  for j in range(4):
   cij=0.0
   for k in range(4):
    cij=up(cij+up(_abs_upper(B[i][k])*R[j][k]))
    cij=up(cij+up(R[i][k]*_abs_upper(B[j][k])))
   rowsum=up(rowsum+cij)
  gamma=max(gamma,rowsum)
 M=matrix_mul(B,matrix_transpose(_pm(Fc)))
 for i in range(4):M[i][i]=M[i][i]+I(rho)-I(gamma)
 Lout,rad=_center_radius_lower(M)
 return Lout,up(rad+gamma)

def _rank1_information_update_lower(L,beta,q):
 """Validated lower for (L^-1 + beta e_q e_q')^-1 without interval cancellation."""
 if not (0<=q<4 and beta>=0 and math.isfinite(beta)):raise RuntimeError('invalid rank-one information update')
 if beta==0:return L,0.0
 A=[[Fraction.from_float(float(L[i][j])) for j in range(4)] for i in range(4)]
 b=Fraction.from_float(float(beta)); den=Fraction(1)+b*A[q][q]
 if den<=0:raise RuntimeError('rank-one information denominator lost positivity')
 M=[[A[i][j]-b*A[i][q]*A[q][j]/den for j in range(4)] for i in range(4)]
 return _exact_rational_loewner_lower(M)

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

def _prop(tau,h,n,rho,betaS,betaA,depth=0):
 try:
  F=_F(tau,h);L=[[0.0]*4 for _ in range(4)];mr=0.0
  for k in range(n):
   if k==0:L=[[rho if i==j else 0.0 for j in range(4)] for i in range(4)]
   else:L,r=_predict(L,F,rho);mr=max(mr,r)
   L,r=_rank1_information_update_lower(L,betaS,2);mr=max(mr,r)
   L,r=_rank1_information_update_lower(L,betaA,3);mr=max(mr,r)
  return [(tau,L,mr,depth)]
 except RuntimeError:
  if depth>=MAX_TAU_SPLIT_DEPTH:raise
  mid=math.sqrt(tau.lo*tau.hi);a=Interval.outward_bounds(tau.lo,mid);b=Interval.outward_bounds(mid,tau.hi)
  return _prop(a,h,n,rho,betaS,betaA,depth+1)+_prop(b,h,n,rho,betaS,betaA,depth+1)

def _mode(mode,p,horizon_s):
 c=WORST.build_cell(mode,p);s=WORST.serializable(c);row=c['row'];h=float(c['sched']['dt_s']);x=c['x'];tau=Interval.outward_bounds(h/x.hi,h/x.lo);sigma=float(c['sigma'].lo);rho=float(c['rho_translation_lower'])
 scale2=[(sigma*h)**2,(sigma*h*h)**2,(sigma*h*h*h)**2,sigma*sigma];u=list(map(float,row['Sigma_diagonal_upper']));physical=[u[6],u[9],u[12],u[15]];upper=[(I(physical[i])/I(scale2[i])).hi for i in range(4)]
 text=WRAPPER.read_text();rh=min(_member(text,'R_S_x_factor_'),_member(text,'R_S_y_factor_'),1.0);rs=(I(rh)*I(float(c['rs'].lo))).lo;rsvar=I(rs).square().lo;acc=float(c['vector']['configured_measurement_bounds']['acc_measurement_std_mps2']);accvar=I(acc).square().lo
 rS=(I(rsvar)/I(scale2[2])).lo;rA=(I(accvar)/I(scale2[3])).lo;betaS=(I(1)/I(rS)).hi;betaA=(I(1)/I(rA)).hi;n=int(math.ceil(horizon_s/h));leaves=_prop(tau,h,n,rho,betaS,betaA);cert=[]
 for t,L,rad,dep in leaves:
  d=_delta(L,upper)
  if d<=0 or not _spd_delta(L,upper,d):raise RuntimeError(f'nonpositive endpoint margin on tau leaf {t.as_list()}')
  cert.append({'tau_s':t.as_list(),'delta_lower':d,'max_conditioned_radius_removed':rad,'split_depth':dep})
 w=min(cert,key=lambda q:q['delta_lower']);old=float(row['direct_translation_generalized_margin_lower'])
 return {'source_cell':s,'conditioned_coordinates':'D^-1[v,p,S,a_w]','tau_interval_s':tau.as_list(),'tau_leaf_count':len(cert),'max_tau_split_depth_used':max(q['split_depth'] for q in cert),'steps':n,'horizon_s':horizon_s,'process_injection_lower_conditioned':rho,'S_measurement_information_beta_conditioned':betaS,'accelerometer_aw_information_beta_conditioned':betaA,'measurement_information_geometry':'rank_one_S_and_aw_each_sample_exact_rational','prediction_enclosure':'midpoint_plus_symmetric_cross_term_loewner','corrections_allowed_every_sample_for_lower_bound':True,'artificial_S_variance_conditioned':rS,'artificial_acc_aw_variance_conditioned':rA,'translation_covariance_upper_conditioned':upper,'tau_leaf_certificates':cert,'complete_word_translation_margin_lower':w['delta_lower'],'limiting_tau_leaf':w,'old_single_seed_translation_margin_lower':old,'margin_widening_factor_lower':down(w['delta_lower']/old),'interval_ldlt_endpoint_recertified':True}
def build(domain_path=DEFAULT_DOMAIN,horizon_s=DEFAULT_HORIZON_S):
 horizon_s=float(horizon_s)
 if not math.isfinite(horizon_s) or horizon_s<=0:raise ValueError('horizon_s must be finite positive')
 p=Path(domain_path).resolve();m={};f=[]
 for mode in ('H','A'):
  try:m[mode]=_mode(mode,p,horizon_s)
  except Exception as e:f.append(f'{mode}: {e}')
 return {'qualification':'OU3_P4_VALIDATED_WORST_CELL_COMPLETE_WORD_TRANSLATION_DISSIPATION','source_only':True,'trajectory_replay_used':False,'outward_rounded':True,'horizon_s':horizon_s,'modes':m,'P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS':'PASS' if not f and len(m)==2 else 'NOT_ESTABLISHED','P4_USABLE_CERTIFICATE_STATUS':'NOT_ESTABLISHED','remaining_obligation':'extend directional complete-word propagation to every reachable source cell/edge and attitude-bias blocks, then validate exact nonlinear return map','failures':f}
def validate(d):
 f=list(d.get('failures',[]))
 if d.get('source_only') is not True or d.get('trajectory_replay_used') is not False or d.get('outward_rounded') is not True:f.append('qualification flags invalid')
 if not float(d.get('horizon_s',0))>0:f.append('invalid horizon')
 for mode in ('H','A'):
  m=d.get('modes',{}).get(mode,{})
  if not float(m.get('complete_word_translation_margin_lower',0))>0:f.append(f'{mode}: no complete-word translation margin')
  if m.get('interval_ldlt_endpoint_recertified') is not True:f.append(f'{mode}: endpoint not recertified')
  if m.get('measurement_information_geometry')!='rank_one_S_and_aw_each_sample_exact_rational':f.append(f'{mode}: directional measurement geometry missing')
  if m.get('prediction_enclosure')!='midpoint_plus_symmetric_cross_term_loewner':f.append(f'{mode}: structured prediction enclosure missing')
  if m.get('corrections_allowed_every_sample_for_lower_bound') is not True:f.append(f'{mode}: lower no longer covers maximum correction frequency')
  if not float(m.get('margin_widening_factor_lower',0))>1:f.append(f'{mode}: complete word did not widen seed')
 if d.get('P4_USABLE_CERTIFICATE_STATUS')!='NOT_ESTABLISHED':f.append('partial result prematurely promoted P4')
 return f
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN);ap.add_argument('--horizon-s',type=float,default=DEFAULT_HORIZON_S);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build(a.domain,a.horizon_s);f=validate(d);d['validation_failures']=f;a.output.write_text(json.dumps(d,indent=2,sort_keys=True));print(json.dumps({'translation_status':d['P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS'],'horizon_s':d['horizon_s'],'modes':{x:{'delta':d.get('modes',{}).get(x,{}).get('complete_word_translation_margin_lower'),'factor':d.get('modes',{}).get(x,{}).get('margin_widening_factor_lower'),'tau_leaves':d.get('modes',{}).get(x,{}).get('tau_leaf_count'),'geometry':d.get('modes',{}).get(x,{}).get('measurement_information_geometry'),'prediction':d.get('modes',{}).get(x,{}).get('prediction_enclosure')} for x in ('H','A')},'failures':f},indent=2));return 0 if not f else 2
if __name__=='__main__':raise SystemExit(main())
