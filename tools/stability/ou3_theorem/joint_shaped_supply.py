"""Joint quadratic-form accumulator for one reachable history cell."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

@dataclass
class JointQuadratic:
    # q(e,u)=e'Ae+2e'Bu+u'Cu; same u symbols feed loss and source.
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    dependency_token: str
    @classmethod
    def zeros(cls,n,m,token):
        return cls(np.zeros((n,n)),np.zeros((n,m)),np.zeros((m,m)),token)
    def add_operation(self,Q0,Q1,Aop,Dsrc):
        Q0=np.asarray(Q0,float);Q1=np.asarray(Q1,float);Aop=np.asarray(Aop,float);Dsrc=np.asarray(Dsrc,float)
        self.A += Q0-Aop.T@Q1@Aop
        self.B += -Aop.T@Q1@Dsrc
        self.C += -Dsrc.T@Q1@Dsrc
        self.A=(self.A+self.A.T)/2;self.C=(self.C+self.C.T)/2
    def matrix(self):
        return np.block([[self.A,self.B],[self.B.T,self.C]])

def linked_source_elimination(q:JointQuadratic, source_gram:np.ndarray):
    """Eliminate one SAME-history source ellipsoid u'R u<=1.

    Returns the generalized Schur object; callers must enclose dependency-aware
    operands before promoting it.
    """
    R=np.asarray(source_gram,float)
    if R.shape!=q.C.shape: raise ValueError("source metric shape")
    # S-procedure family A - B (lambda R + C)^-1 B'; lambda is retained for
    # optimization by the constructive cell solver, not fixed here.
    return {"A":q.A,"B":q.B,"C":q.C,"R":R,
            "dependency_token":q.dependency_token,
            "source_uniform_verified":False}
