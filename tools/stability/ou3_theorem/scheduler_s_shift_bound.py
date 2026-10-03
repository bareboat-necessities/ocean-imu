"""Analytical local S-shift commutator identity and norm bound.

R(P)=P-PH'(HPH'+R)^-1HP.  Compare R(FPF'+Q) with F R(P) F'+Q.
The difference is exactly
  F P H' A^-1 H P F' - C H' B^-1 H C,
where C=FPF'+Q, A=HPH'+R, B=HCH'+R,
when H is the same coordinate selector before/after the prediction.
For OU-III S=0, H selects S and F commutes with that selector only up to the
literal LIN chain, so the exact formula is evaluated with literal F.
"""
from __future__ import annotations
import numpy as np

def ric(P,H,R):
 S=H@P@H.T+R
 return P-P@H.T@np.linalg.solve(S,H@P)

def defect(P,F,Q,H,R):
 C=F@P@F.T+Q
 return ric(C,H,R)-(F@ric(P,H,R)@F.T+Q)

def identity_rhs(P,F,Q,H,R):
 C=F@P@F.T+Q
 A=H@P@H.T+R;B=H@C@H.T+R
 return F@P@H.T@np.linalg.solve(A,H@P)@F.T-C@H.T@np.linalg.solve(B,H@C)

def norm_bound(Plo,Phi,F,Q,H,R):
 """Rigorous but possibly loose spectral bound over Plo<=P<=Phi.

 Uses ||PH'|| <= ||Phi H'|| and innovation inverse <= 1/lambda_min(R).
 Bounds each rank-3 removal term separately. This is an interval-free outer
 bound; a linked difference bound can sharpen it if needed.
 """
 C_hi=F@Phi@F.T+Q
 rmin=float(np.linalg.eigvalsh(R).min())
 t1=np.linalg.norm(F@Phi@H.T,2)**2/rmin
 t2=np.linalg.norm(C_hi@H.T,2)**2/rmin
 return float(t1+t2)

def certificate():
 rng=np.random.default_rng(5);A=rng.normal(size=(7,7));P=A@A.T+.2*np.eye(7)
 F=np.eye(7)+.01*rng.normal(size=(7,7));Z=rng.normal(size=(7,7));Q=1e-4*Z@Z.T
 H=np.zeros((2,7));H[0,3]=1;H[1,4]=1;R=.1*np.eye(2)
 D=defect(P,F,Q,H,R);E=identity_rhs(P,F,Q,H,R)
 return {"qualification":"OU3_S_SHIFT_COMMUTATOR_IDENTITY_V1",
         "identity_defect":float(np.linalg.norm(D-E)),
         "defect_symmetric":float(np.linalg.norm(D-D.T)),
         "two_term_outer_bound":norm_bound(.9*P,1.1*P,F,Q,H,R),
         "phase_uniform_bound_verified":False,"theorem_closed":False}
if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
