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
