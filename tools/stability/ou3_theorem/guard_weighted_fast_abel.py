"""Integrated Abel bound for weighted guard FAST correction.

For Q = h sum_k c_k w_k ((L-I)u)_k, collect coefficients q_j of u_j.
With primitive U_j=h sum_{i<=j}u_i, summation by parts gives
|Q| <= C*(endpoint + total variation of q), with no C/h penalty.
This is the useful global application of the FAST signed-primitive premise.
"""
from __future__ import annotations
import math

def lp_impulse(alpha,n):
 return [(1-alpha)*(alpha**j) for j in range(n)]

def integrated_coefficients(alpha,weights,reader=None):
 n=len(weights);reader=[1.]*n if reader is None else list(reader)
 if len(reader)!=n:raise ValueError("reader length")
 l=lp_impulse(alpha,n);q=[0.]*n
 # output k: w_k[(L u)_k-u_k], then weighted by reader c_k.
 for k in range(n):
  cw=float(reader[k])*float(weights[k])
  q[k]-=cw
  for j in range(k+1):
   q[j]+=cw*l[k-j]
 return q

def abel_charge(q):
 if not q:return 0.
 return abs(q[-1])+sum(abs(q[j+1]-q[j]) for j in range(len(q)-1))

def integrated_primitive_bound(alpha,dt,C,B,H,weights,reader=None):
 q=integrated_coefficients(alpha,weights,reader)
 nH=max(1,int(math.floor(H/dt)))
 # If the whole functional lies inside H, the placed-window primitive applies
 # directly. Longer functionals must be tiled; do not silently reset history.
 if len(q)>nH:
  raise ArithmeticError("weighted guard functional exceeds FAST primitive horizon; tile with carried boundary primitive")
 charge=abel_charge(q)
 return {"bound":C*charge,"charge":charge,"coefficients":q,
         "uses_signed_primitive_cap":True,"pointwise_claim":False,
         "weight_variation":sum(abs(weights[k+1]-weights[k]) for k in range(len(weights)-1))}

def constant_weight(alpha,dt,C,H,n,w):
 return integrated_primitive_bound(alpha,dt,C,0.,H,[w]*n)
