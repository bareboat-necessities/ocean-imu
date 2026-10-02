"""Kernel-line enclosure from the literal joint attitude/BA graph.

For a cell, attitude null vector a and BA graph A_ba define r=(a,-A_ba a)
in the joint six coordinates, embedded into the 21-state root.
"""
from __future__ import annotations
import numpy as np

def joint_line(attitude_line,A_ba):
 a=np.asarray(attitude_line,float).reshape(3);G=np.asarray(A_ba,float)
 if G.shape!=(3,3):raise ValueError("3x3 BA graph")
 r=np.zeros(21);r[:3]=a;r[18:21]=-(G@a)
 if np.linalg.norm(r)==0:raise ValueError("nonzero line")
 return r

def line_angle_radius(mid,component_radius):
 r=np.asarray(mid,float);rad=np.asarray(component_radius,float)
 n=np.linalg.norm(r);eps=np.linalg.norm(rad)
 if eps>=n:return {"verified":False,"angle_rad":float("inf")}
 # Any vector in component box differs by <=eps; normalized direction angle <= asin(eps/n).
 return {"verified":True,"angle_rad":float(np.arcsin(min(1.,eps/n))),
         "euclidean_radius":float(eps),"mid_norm":float(n)}
