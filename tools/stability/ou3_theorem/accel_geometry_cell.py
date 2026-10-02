"""Dependency-preserving accelerometer Jacobian cell constructor.

Parameterizes H_acc by one force vector and one rotation cell instead of
independent H entries. This is the first refinement required by the release
innovation certificate.
"""
from __future__ import annotations
import numpy as np
from .interval_riccati_21 import IMat,N

def accel_h_from_force_rotation(f_mid,f_rad,R_mid,R_rad):
 f=np.asarray(f_mid,float);fr=np.asarray(f_rad,float);Rm=np.asarray(R_mid,float);Rr=np.asarray(R_rad,float)
 if f.shape!=(3,) or fr.shape!=(3,) or Rm.shape!=(3,3) or Rr.shape!=(3,3):raise ValueError("geometry shapes")
 # J_att=-[f]x
 hm=np.zeros((3,N));hr=np.zeros((3,N))
 skew=np.array([[0,-f[2],f[1]],[f[2],0,-f[0]],[-f[1],f[0],0.]])
 sr=np.array([[0,fr[2],fr[1]],[fr[2],0,fr[0]],[fr[1],fr[0],0.]])
 hm[:,:3]=-skew;hr[:,:3]=sr
 hm[:,15:18]=Rm;hr[:,15:18]=Rr
 for i in range(3):hm[i,18+i]=1.
 return IMat(tuple(tuple(float(x) for x in row) for row in hm),
             tuple(tuple(float(x) for x in row) for row in hr))

def rotation_cell_from_small_angle(theta_mid,theta_rad):
 """First-order SO(3) cell with rigorous quadratic remainder for |theta|<.2."""
 t=np.asarray(theta_mid,float);r=np.asarray(theta_rad,float);rho=np.linalg.norm(np.abs(t)+r)
 if rho>=.2:raise ArithmeticError("attitude cell too wide for local rotation enclosure")
 def skew(x):return np.array([[0,-x[2],x[1]],[x[2],0,-x[0]],[-x[1],x[0],0.]])
 Rm=np.eye(3)+skew(t);rem=.5*rho*rho/(1-rho/3)
 Rr=np.abs(skew(r))+rem*np.ones((3,3))
 return Rm,Rr
