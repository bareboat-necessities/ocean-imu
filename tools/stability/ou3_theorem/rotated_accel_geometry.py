# ruff: noqa: F401, F811
"""Joint accelerometer geometry in a rotated innovation frame.

Left multiplication of innovation by an orthogonal Q preserves the Kalman
correction when R is transformed consistently.  Choosing Q=R_wb' combines the
AW block exactly and carries force into body coordinates before enclosure.
"""
from __future__ import annotations
import numpy as np
from .interval_riccati_21 import IMat,N
from .rank_loss_interval_factor import exact

def skew(x):
 x=np.asarray(x,float);return np.array([[0,-x[2],x[1]],[x[2],0,-x[0]],[-x[1],x[0],0.]])

def transformed_point_h(Rwb,f_world):
 """Q H with Q=Rwb'; H=[-[f]x,0,0,0,Rwb,I].

 Identity R' [f]x = [R'f]x R' gives transformed attitude block
 -[f_b]x R', AW=I, BA=R'.
 """
 R=np.asarray(Rwb,float);f=np.asarray(f_world,float);Rt=R.T;fb=Rt@f
 H=np.zeros((3,N));H[:,:3]=-skew(fb)@Rt;H[:,15:18]=np.eye(3);H[:,18:21]=Rt
 return H,fb

def transformed_h_cell(Rmid,Rrad,fbody_mid,fbody_rad):
 """Dependency-preserving transformed H cell.

 The physical force is represented in the SAME local body frame. AW block is
 exact I. BA and attitude right factor share the same R' interval.
 """
 R=np.asarray(Rmid,float);Rr=np.asarray(Rrad,float);fb=np.asarray(fbody_mid,float);fr=np.asarray(fbody_rad,float)
 Rt=R.T;Rtr=Rr.T
 Hm=np.zeros((3,N));Hr=np.zeros((3,N))
 Sm=skew(fb);Sr=np.array([[0,fr[2],fr[1]],[fr[2],0,fr[0]],[fr[1],fr[0],0.]])
 Hm[:,:3]=-Sm@Rt
 # product radius |Sm|Rr + Sr|Rt| + Sr Rr
 Hr[:,:3]=np.abs(Sm)@Rtr+Sr@np.abs(Rt)+Sr@Rtr
 Hm[:,15:18]=np.eye(3)       # exact cancellation R'R
 Hm[:,18:21]=Rt;Hr[:,18:21]=Rtr
 return IMat(tuple(map(tuple,Hm)),tuple(map(tuple,Hr)))

def rotate_measurement_covariance_isotropic(Racc:IMat):
 """Shipping proof Racc enclosure is isotropic diagonal; rotation is invariant."""
 # Fail closed if interval representation contains off-diagonal/nonuniform center.
 m=np.asarray(Racc.mid,float);r=np.asarray(Racc.rad,float)
 if np.max(np.abs(m-np.eye(3)*m[0,0]))>1e-12:return None
 if np.max(np.abs(r-np.eye(3)*r[0,0]))>1e-12:return None
 return Racc
