# ruff: noqa: F401, F811
"""Closed-loop release operation stream from one causal history leaf.

This is proof-side mirror algebra, not a second estimator.  It uses the literal
shipping linear OU mean chronology and the same covariance-derived gain for
mean correction and reset.  Interval promotion is fail-closed where a history
cell does not yet separate a generated branch.
"""
from __future__ import annotations
import math,numpy as np
from .golive_release_seed import covariance_interval
from .interval_riccati_21 import (IMat,predict_covariance,verified_gain_interval_psd,
 joseph_covariance,shipping_prediction_intervals,shipping_acc_update_intervals,
 shipping_integral_update_intervals,selector_h)
from .source_uniform_release_chronology import Scheduler,pseudo_period
from .same_history_geometry_residual import skew
from .interval_kr import interval_matvec
N=21;G=9.80665

def point_imat(a):
 a=np.asarray(a,float);return IMat(tuple(tuple(float(x) for x in r) for r in a),
  tuple(tuple(0. for _ in r) for r in a))

def mid(a):return np.asarray(a.mid,float)

def ou_coeffs(h,tau):
 x=h/max(tau,1e-7);alpha=min(math.exp(-x),1.);pva=-tau*math.expm1(-x)
 if abs(x)<1e-2:
  x2=x*x;x3=x2*x;x4=x3*x;x5=x4*x
  ppa=tau*tau*(.5*x2-x3/6+x4/24)
  psa=tau**3*(x3/6-x4/24+x5/120)
 else:
  em1=math.expm1(-x);ppa=tau*tau*(x+em1);psa=tau**3*(.5*x*x-x-em1)
 return pva,ppa,psa,alpha

def predict_mean(x,h,tau,tau_ba=120.):
 y=np.asarray(x,float).copy();pva,ppa,psa,alpha=ou_coeffs(h,tau)
 v=x[6:9];p=x[9:12];S=x[12:15];aw=x[15:18]
 y[6:9]=v+pva*aw;y[9:12]=p+h*v+ppa*aw
 y[12:15]=S+h*p+.5*h*h*v+psa*aw;y[15:18]=alpha*aw
 y[18:21]=math.exp(-h/tau_ba)*x[18:21]
 return y,alpha

def exact_world_acc_H(x):
 """World-coordinate row orthogonally equivalent to literal body row.

 BA column is identity in the shipping local error coordinates; AW is identity
 after the common R factor is removed.
 """
 f=x[15:18]-np.array([0.,0.,G]);H=np.zeros((3,N))
 H[:,:3]=-skew(f);H[:,15:18]=np.eye(3);H[:,18:21]=np.eye(3)
 return H

def correction(x,P,H,R,residual,sensor):
 Hi=point_imat(H);Ri=point_imat(R)
 K,cert=verified_gain_interval_psd(P,Hi,Ri);Km=mid(K)
 dx=Km@np.asarray(residual,float);y=x+dx
 Pn=joseph_covariance(P,K,Hi,Ri)
 return y,Pn,{"kind":"correction","sensor":sensor,"H":Hi,"R":Ri,"K":K,
  "residual":np.asarray(residual,float),"dx":dx,"inverse_certificate":cert}

def reset_after_correction(x,P,dx):
 """Literal first-order covariance reset; nominal error-state theta is injected."""
 d=np.asarray(dx[:3],float);Gm=np.eye(N);Gm[:3,:3]+=0.5*skew(d)
 Gi=point_imat(Gm);Pn=__import__(
  "tools.stability.ou3_theorem.interval_riccati_21",fromlist=["matmul","transpose"])
 mm=Pn.matmul(Pn.matmul(Gi,P),Pn.transpose(Gi))
 y=x.copy();y[:3]=0.
 return y,mm,{"kind":"reset","G":Gi,"dtheta":d}

def generate_point_stream(history_payload,dt=.005):
 """Run one point-valued same-history leaf.

 Source-uniform interval cells are accepted only when every physical/tuner
 quantity used here is point-valued.  Otherwise caller must split the shared
 history cell; midpoint substitution is forbidden.
 """
 samples=history_payload["samples"];adapt=history_payload["adaptation_trace"]
 if len(samples)!=len(adapt):raise ValueError("history/adaptation length mismatch")
 x=np.zeros(N);P=covariance_interval();sch=Scheduler(0.,.015);events=[]
 for k,(row,tune) in enumerate(zip(samples,adapt)):
  def pointI(q,name):
   if q.lo!=q.hi:raise ArithmeticError(name+" history cell not point-resolved")
   return q.lo
  tau=pointI(tune["tau"],"tau");sigma=pointI(tune["sigma_aw"],"sigma_aw")
  h=dt;x,alpha=predict_mean(x,h,tau)
  # Covariance predictor retains literal ranges but pins this generated tau/sigma.
  F,Q=shipping_prediction_intervals(dt_min=h,dt_max=h,tau_min=tau,tau_max=tau,
    omega_max=.6108652381980153,tau_bacc=120.,gyro_white_density=.00135,
    gyro_bias_rw_density=math.sqrt(1e-11),aw_sigma_max=max(sigma,1e-6),
    accel_bias_drive_density=5e-4)
  P=predict_covariance(P,F,Q);events.append({"sample":k,"kind":"prediction","F":F,"Q":Q,
    "phi_aw":alpha,"mean":x.copy()})
  period=pseudo_period(tau);sch.retarget(period)
  if sch.advance(h):
   rs=pointI(tune["R_S"],"R_S");H=mid(selector_h(12));R=np.eye(3)*(rs*rs)
   r=-x[12:15];x,P,e=correction(x,P,H,R,r,"S");e["sample"]=k;events.append(e)
   x,P,e=reset_after_correction(x,P,e["dx"]);e["sample"]=k;events.append(e)
  # Physical acceleration and delivered sensor error must be point-resolved.
  pa=np.array([pointI(history_payload["physical"]["acceleration"][k].value_component[a]
      if False else row[f"physical_v_{a}"],"unused") for a in range(3)]) if False else None
  # literal_history_leaf_propagator exports the acceleration box separately below.
  avec=history_payload["physical_acceleration"][k]
  a=np.array([pointI(q,"physical acceleration") for q in avec])
  slow=np.array([pointI(row[f"slow_accel_{j}"],"slow accel") for j in range(3)])
  fast=np.array([pointI(history_payload["fast_accel"][k][j],"fast accel") for j in range(3)])
  sensor=slow+fast
  H=exact_world_acc_H(x);r=a-x[15:18]+sensor-x[18:21]
  # Racc uses nominal value here; a promoted interval leaf must carry the
  # vibration-conditioned R from the same front-end trace.
  R=np.eye(3)*(.2**2)
  x,P,e=correction(x,P,H,R,r,"acc");e["sample"]=k;events.append(e)
  x,P,e=reset_after_correction(x,P,e["dx"]);e["sample"]=k;events.append(e)
 return {"mean":x,"P":P,"events":events,"samples":len(samples),
         "same_history_closed_loop":True,"magnetic_callback_pattern_enumerated":False}

def interval_correction_product(P,H,R,r_mid,r_rad):
 """Covariance-derived interval K and SAME-cell K*r correction."""
 K,cert=verified_gain_interval_psd(P,H,R)
 dm,dr=interval_matvec(K,r_mid,r_rad)
 return {"K":K,"dx_mid":dm,"dx_rad":dr,
         "dtheta_mid":dm[:3],"dtheta_rad":dr[:3],
         "inverse_certificate":cert,"same_Kr_cell":True}

def racc_imat(entry):
 return IMat(tuple(tuple(float(x) for x in row) for row in entry["mid"]),
             tuple(tuple(float(x) for x in row) for row in entry["rad"]))

def _ival(q):
 return .5*(q.lo+q.hi),.5*(q.hi-q.lo)

def interval_world_acc_H(xm,xr):
 """H cell from the SAME AW mean interval."""
 fm=xm[15:18]-np.array([0.,0.,G]);fr=xr[15:18]
 Hm=np.zeros((3,N));Hr=np.zeros((3,N))
 Hm[:,:3]=-skew(fm)
 Hr[:,:3]=np.array([[0,fr[2],fr[1]],[fr[2],0,fr[0]],[fr[1],fr[0],0.]])
 Hm[:,15:18]=np.eye(3);Hm[:,18:21]=np.eye(3)
 return IMat(tuple(tuple(float(x) for x in row) for row in Hm),
             tuple(tuple(float(x) for x in row) for row in Hr))

def interval_predict_mean(xm,xr,h,tau_lo,tau_hi,tau_ba=120.):
 """Hull literal OU coefficients over tau endpoints; monotone enclosure by hull."""
 vals=[ou_coeffs(h,t) for t in (tau_lo,tau_hi)]
 coeff=[(min(z[j] for z in vals),max(z[j] for z in vals)) for j in range(4)]
 def mul_interval(cm,cr,lo,hi):
  cc=.5*(lo+hi);rr=.5*(hi-lo)
  return cc*cm,abs(cc)*cr+rr*np.abs(cm)+rr*cr
 ym=xm.copy();yr=xr.copy()
 for sl,j in ((slice(6,9),0),(slice(9,12),1),(slice(12,15),2),(slice(15,18),3)):
  pass
 vm,vr=xm[6:9],xr[6:9];pm,pr=xm[9:12],xr[9:12];sm,sr=xm[12:15],xr[12:15];am,ar=xm[15:18],xr[15:18]
 avm,avr=mul_interval(am,ar,*coeff[0]);apm,apr=mul_interval(am,ar,*coeff[1]);asm,asr=mul_interval(am,ar,*coeff[2]);aam,aar=mul_interval(am,ar,*coeff[3])
 ym[6:9]=vm+avm;yr[6:9]=vr+avr
 ym[9:12]=pm+h*vm+apm;yr[9:12]=pr+h*vr+apr
 ym[12:15]=sm+h*pm+.5*h*h*vm+asm;yr[12:15]=sr+h*pr+.5*h*h*vr+asr
 ym[15:18]=aam;yr[15:18]=aar
 b=math.exp(-h/tau_ba);ym[18:21]=b*xm[18:21];yr[18:21]=b*xr[18:21]
 return ym,np.nextafter(yr,np.inf),coeff[3]

def interval_reset_G(dm,dr):
 Gm=np.eye(N);Gr=np.zeros((N,N));Gm[:3,:3]+=.5*skew(dm)
 Gr[:3,:3]=.5*np.array([[0,dr[2],dr[1]],[dr[2],0,dr[0]],[dr[1],dr[0],0.]])
 return IMat(tuple(tuple(float(x) for x in row) for row in Gm),
             tuple(tuple(float(x) for x in row) for row in Gr))

def _dominant_radius(xr,rres,K=None,dtheta=None,guard=None):
 vals={"mean_AW":float(np.max(xr[15:18])),"mean_BA":float(np.max(xr[18:21])),
       "residual":float(np.max(rres))}
 if K is not None:vals["K"]=float(max(max(row) for row in K.rad))
 if dtheta is not None:vals["reset"]=float(np.max(dtheta))
 if guard is not None:vals["guard"]=float(guard)
 return max(vals,key=vals.get),vals

def generate_interval_stream(history_payload,dt=.005,max_samples=None):
 """Closed-loop midpoint-radius mean/covariance stream for one HistoryCell."""
 samples=history_payload["samples"];adapt=history_payload["adaptation_trace"]
 n=len(samples) if max_samples is None else min(len(samples),max_samples)
 xm=np.zeros(N);xr=np.zeros(N);P=covariance_interval();sch=Scheduler(0.,.015);events=[]
 for k in range(n):
  row,tune=samples[k],adapt[k]
  try:
   xm,xr,phi=interval_predict_mean(xm,xr,dt,tune["tau"].lo,tune["tau"].hi)
   F,Q=shipping_prediction_intervals(dt_min=dt,dt_max=dt,tau_min=tune["tau"].lo,tau_max=tune["tau"].hi,
    omega_max=.6108652381980153,tau_bacc=120.,gyro_white_density=.00135,
    gyro_bias_rw_density=math.sqrt(1e-11),aw_sigma_max=max(tune["sigma_aw"].hi,1e-6),
    accel_bias_drive_density=5e-4)
   P=predict_covariance(P,F,Q);events.append({"sample":k,"kind":"prediction","F":F,"Q":Q})
   plo,pup=pseudo_period(tune["tau"].lo),pseudo_period(tune["tau"].hi)
   if plo!=pup:raise ArithmeticError("S scheduler period interval unresolved")
   sch.retarget(plo)
   if sch.advance(dt):
    H=selector_h(12);rslo,rshi=tune["R_S"].lo,tune["R_S"].hi
    R=IMat(tuple(tuple((.5*(rslo*rslo+rshi*rshi) if i==j else 0.) for j in range(3)) for i in range(3)),
           tuple(tuple((.5*(rshi*rshi-rslo*rslo) if i==j else 0.) for j in range(3)) for i in range(3)))
    rm,rr=-xm[12:15],xr[12:15];q=interval_correction_product(P,H,R,rm,rr)
    xm+=q["dx_mid"];xr=np.nextafter(xr+q["dx_rad"],np.inf);P=joseph_covariance(P,q["K"],H,R)
    G=interval_reset_G(q["dtheta_mid"],q["dtheta_rad"])
    from .interval_riccati_21 import matmul,transpose
    P=matmul(matmul(G,P),transpose(G));xm[:3]=0.;xr[:3]=0.
    events.append({"sample":k,"kind":"S","H":H,"R":R,**q,"G":G})
   H=interval_world_acc_H(xm,xr)
   am=[];ar=[];em=[];er=[]
   for j in range(3):
    a0,a1=_ival(history_payload["physical_acceleration"][k][j]);am.append(a0);ar.append(a1)
    # conditioned delivered acceleration minus physical ideal acceleration is the
    # same-history sensor source after the guard.
    c0,c1=_ival(history_payload["conditioned_accel"][k][j]);em.append(c0-a0);er.append(c1+a1)
   am=np.array(am);ar=np.array(ar);em=np.array(em);er=np.array(er)
   rm=am-xm[15:18]+em-xm[18:21];rr=ar+xr[15:18]+er+xr[18:21]
   R=racc_imat(history_payload["Racc_interval"][k])
   q=interval_correction_product(P,H,R,rm,rr)
   xm+=q["dx_mid"];xr=np.nextafter(xr+q["dx_rad"],np.inf);P=joseph_covariance(P,q["K"],H,R)
   G=interval_reset_G(q["dtheta_mid"],q["dtheta_rad"])
   from .interval_riccati_21 import matmul,transpose
   P=matmul(matmul(G,P),transpose(G));xm[:3]=0.;xr[:3]=0.
   dom,parts=_dominant_radius(xr,rr,q["K"],q["dtheta_rad"],
      history_payload["Racc_interval"][k]["rad"][0][0])
   events.append({"sample":k,"kind":"acc","H":H,"R":R,**q,"G":G,"dominant_radius":dom,"radius_parts":parts})
   if not(np.all(np.isfinite(xm)) and np.all(np.isfinite(xr))):
    return {"verified":False,"first_failure_sample":k,"reason":"nonfinite mean enclosure",
            "dominant_radius":dom,"radius_parts":parts,"events":events}
  except (ArithmeticError,ValueError,OverflowError) as e:
   parts={"mean_AW":float(np.max(xr[15:18])),"mean_BA":float(np.max(xr[18:21]))}
   return {"verified":False,"first_failure_sample":k,"reason":str(e),
           "dominant_radius":max(parts,key=parts.get),"radius_parts":parts,"events":events}
 return {"verified":True,"samples":n,"mean_mid":xm,"mean_rad":xr,"P":P,"events":events,
         "magnetic_callback_pattern_enumerated":False}
