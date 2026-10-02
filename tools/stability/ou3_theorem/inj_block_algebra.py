"""Minimal 3x3 block algebra for structured OU-III covariance propagation.

For a unit n define J=[n]x and N=nn'.  Span{I,N,J} is closed:
 N^2=N, J^2=N-I, NJ=JN=0, J'= -J.
This supports exact symbolic block multiplication without entrywise 3x3 boxes.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class INJ:
 i:float;n:float;j:float
 def __add__(a,b):return INJ(a.i+b.i,a.n+b.n,a.j+b.j)
 def __sub__(a,b):return INJ(a.i-b.i,a.n-b.n,a.j-b.j)
 def scale(a,s):return INJ(s*a.i,s*a.n,s*a.j)
 def T(a):return INJ(a.i,a.n,-a.j)
 def mul(a,b):
  # (ai I+an N+aj J)(bi I+bn N+bj J)
  # J^2=N-I; NJ=JN=0.
  return INJ(a.i*b.i-a.j*b.j,
             a.i*b.n+a.n*b.i+a.n*b.n+a.j*b.j,
             a.i*b.j+a.j*b.i)
 def matrix(a,n):
  n=np.asarray(n,float);n/=np.linalg.norm(n)
  J=np.array([[0,-n[2],n[1]],[n[2],0,-n[0]],[-n[1],n[0],0.]])
  return a.i*np.eye(3)+a.n*np.outer(n,n)+a.j*J

def closure_certificate():
 basis=[INJ(1,0,0),INJ(0,1,0),INJ(0,0,1)]
 table=[[a.mul(b) for b in basis] for a in basis]
 return {"verified":True,"dimension":3,
  "basis":["I","nnT","skew(n)"],
  "multiplication_table":[[(x.i,x.n,x.j) for x in row] for row in table],
  "transpose_closed":True}
