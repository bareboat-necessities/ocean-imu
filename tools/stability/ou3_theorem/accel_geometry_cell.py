"""Dependency-preserving accelerometer Jacobian cell constructor.

Parameterizes H_acc by one force vector and one rotation cell instead of
independent H entries. This is the first refinement required by the release
innovation certificate.
"""
from __future__ import annotations
import numpy as np
import math
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

def rotation_cell_from_ball(radius):
 """Rigorous SO(3) entry enclosure for ||theta||<=radius via Rodrigues.

 R=I+sinc(r)[theta]x+(1-cos r)/r^2 [theta]x^2.
 For r<=radius, off-diagonal first term <=sin(radius), while every entry of
 [u]x^2=uu'-I has magnitude <=1 and coefficient <=1-cos(radius).
 Center at I; use common entry radius sin(rho)+(1-cos(rho)).
 This intentionally retains the norm-ball premise instead of a component cube.
 """
 rho=float(radius)
 if not(0<=rho<math.pi/2):raise ValueError("rotation ball radius")
 rad=math.sin(rho)+(1-math.cos(rho))
 return np.eye(3),np.full((3,3),rad)

def force_cell_from_ball(radius):
 """Component enclosure plus norm metadata for ||f||<=radius."""
 r=float(radius)
 if not(r>0):raise ValueError("force radius")
 return np.zeros(3),np.full(3,r),{"euclidean_radius":r,"component_cube_is_outer_only":True}

def rodrigues(theta):
 t=np.asarray(theta,float);a=np.linalg.norm(t)
 if a==0:return np.eye(3)
 K=np.array([[0,-t[2],t[1]],[t[2],0,-t[0]],[-t[1],t[0],0.]])/a
 return np.eye(3)+math.sin(a)*K+(1-math.cos(a))*(K@K)

def rotation_cell_centered(center,radius):
 """R=R(center) R(delta), ||delta||<=radius; interval around R(center).

 Multiplication by an orthogonal center preserves spectral perturbation norm.
 Entrywise radius uses ||R(delta)-I||2 <= 2 sin(radius/2).
 """
 Rc=rodrigues(center);eps=2*math.sin(float(radius)/2)
 return Rc,np.full((3,3),eps)
