"""Verified full-21 precision and interval joint shaped-supply arithmetic."""
from __future__ import annotations
import math
from dataclasses import dataclass
from .interval_riccati_21 import IMat,add,matmul,scale,transpose
from .rank_loss_interval_factor import verified_inverse

def sub(a,b):return add(a,scale(b,-1))

def verified_precision(P:IMat):
 Q,cert=verified_inverse(P)
 if not cert.get("verified"):raise ArithmeticError("full covariance inverse not verified")
 return Q,cert

@dataclass
class IJointQuadratic:
 A:IMat;B:IMat;C:IMat;dependency_token:str
 @classmethod
 def zeros(cls,n,m,token):
  z=lambda r,c:IMat(tuple(tuple(0. for _ in range(c)) for _ in range(r)),tuple(tuple(0. for _ in range(c)) for _ in range(r)))
  return cls(z(n,n),z(n,m),z(m,m),token)
 def add_operation(self,Q0:IMat,Q1:IMat,Aop:IMat,Dsrc:IMat):
  # q=e'(Q0-A'Q1A)e -2 e'A'Q1D u -u'D'Q1D u
  self.A=add(self.A,sub(Q0,matmul(matmul(transpose(Aop),Q1),Aop)))
  self.B=add(self.B,scale(matmul(matmul(transpose(Aop),Q1),Dsrc),-1))
  self.C=add(self.C,scale(matmul(matmul(transpose(Dsrc),Q1),Dsrc),-1))
  self.A=scale(add(self.A,transpose(self.A)),.5);self.C=scale(add(self.C,transpose(self.C)),.5)

def interval_scalar_bounds(a:IMat):
 if a.shape!=(1,1):raise ValueError("scalar")
 return math.nextafter(a.mid[0][0]-a.rad[0][0],-math.inf),math.nextafter(a.mid[0][0]+a.rad[0][0],math.inf)
