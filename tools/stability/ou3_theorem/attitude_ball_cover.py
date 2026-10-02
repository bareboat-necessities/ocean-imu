"""Certified lattice cover of the captured rotation-vector ball."""
from __future__ import annotations
from dataclasses import dataclass
import math,itertools,numpy as np
CAP=math.radians(6.9)
BETA=math.radians(1.10)  # strict margin below 1.15147686 deg threshold

@dataclass(frozen=True)
class AttitudeCell:
 center:np.ndarray;radius:float

def lattice_cover(cap=CAP,beta=BETA):
 """Cubic lattice with spacing h=2 beta/sqrt(3).

 Every point in R3 lies within sqrt(3)h/2=beta of a lattice center.
 Retain centers whose beta-ball intersects captured cap ball.
 """
 h=2*beta/math.sqrt(3);m=math.ceil((cap+beta)/h);cells=[]
 for ijk in itertools.product(range(-m,m+1),repeat=3):
  c=h*np.asarray(ijk,float)
  if np.linalg.norm(c)<=cap+beta+1e-15:cells.append(AttitudeCell(c,beta))
 return cells,{"spacing_rad":h,"cover_radius_rad":math.sqrt(3)*h/2,
               "cover_radius_deg":math.degrees(math.sqrt(3)*h/2),
               "cell_count":len(cells),"captured_radius_deg":math.degrees(cap),
               "analytic_full_space_cover":True}

def verify_cover_metadata():
 cells,c=lattice_cover()
 return {**c,"verified":c["cover_radius_rad"]<=BETA+1e-15,
         "beta_threshold_margin_deg":1.1514768600463867-c["cover_radius_deg"]}
