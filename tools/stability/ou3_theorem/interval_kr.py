"""Joint interval K*r products for release corrections."""
from __future__ import annotations
import numpy as np
from .interval_riccati_21 import IMat

def vector_interval(mid,rad):
 m=np.asarray(mid,float);r=np.asarray(rad,float)
 if m.shape!=(3,) or r.shape!=(3,) or np.any(r<0):raise ValueError("3-vector interval")
 return m,r

def interval_matvec(K:IMat,r_mid,r_rad):
 """Outward algebraic enclosure of one K*r using the SAME K and r cells."""
 rm,rr=vector_interval(r_mid,r_rad)
 n,m=K.shape
 if m!=3:raise ValueError("K must have three innovation columns")
 outm=np.zeros(n);outr=np.zeros(n)
 for i in range(n):
  for j in range(3):
   km,kr=K.mid[i][j],K.rad[i][j]
   outm[i]+=km*rm[j]
   outr[i]+=abs(km)*rr[j]+abs(rm[j])*kr+kr*rr[j]
 return outm,np.nextafter(np.nextafter(outr,np.inf),np.inf)

def interval_residual_acc(a_mid,a_rad,aw_mid,aw_rad,sensor_mid,sensor_rad,ba_mid,ba_rad):
 m=np.asarray(a_mid)-np.asarray(aw_mid)+np.asarray(sensor_mid)-np.asarray(ba_mid)
 r=np.asarray(a_rad)+np.asarray(aw_rad)+np.asarray(sensor_rad)+np.asarray(ba_rad)
 return m,r

def interval_residual_S(S_mid,S_rad):
 return -np.asarray(S_mid,float),np.asarray(S_rad,float)

def attitude_reset_interval(dx_mid,dx_rad):
 return np.asarray(dx_mid[:3],float),np.asarray(dx_rad[:3],float)

def audit(K,r_mid,r_rad):
 dm,dr=interval_matvec(K,r_mid,r_rad)
 return {"dx_mid":dm,"dx_rad":dr,"dtheta_mid":dm[:3],"dtheta_rad":dr[:3],
         "same_Kr_cell":True,"independent_gain_residual_suprema":False}
