"""Same-history axis mismatch geometry for structured Riccati propagation."""
from __future__ import annotations
import math,numpy as np
from .changing_axis_inj import TwoAxis,rank_of_basis

def projector_difference_norm(angle):
 """||nn'-mm'||2 = sin(angle) for unit axes."""
 return abs(math.sin(angle))

def skew_difference_norm(angle):
 """||[n]x-[m]x||2 = ||n-m|| = 2 sin(angle/2)."""
 return 2*abs(math.sin(angle/2))

def rebase_bound(coeff_n,coeff_j,angle):
 """A_n(a,b,c) represented around A_m(a,b,c) plus exact mismatch norm."""
 return abs(coeff_n)*projector_difference_norm(angle)+abs(coeff_j)*skew_difference_norm(angle)

def minimal_two_axis_rank(angle_deg):
 n=np.array([0.,0.,1.]);a=math.radians(angle_deg);m=np.array([math.sin(a),0,math.cos(a)])
 return rank_of_basis(TwoAxis(n,m).basis())

def diagnostic():
 return [{"angle_deg":x,"projector_delta":projector_difference_norm(math.radians(x)),
          "skew_delta":skew_difference_norm(math.radians(x)),
          "two_axis_rank":minimal_two_axis_rank(x)}
         for x in (0,.1,1,2,5,10,30,90)]
