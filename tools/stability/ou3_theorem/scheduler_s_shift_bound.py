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
 """Valid two-term spectral bound, evaluated in ordinary floating arithmetic.

 P<=Phi does NOT imply ||P H'||<=||Phi H'||. Instead PSD block Cauchy
 gives ||F P H'||^2 <= ||F Phi F'|| ||H Phi H'||. Each removal is
 also bounded by its complete prior covariance. Take the smaller of these
 two analytical ceilings; no directed-rounding or all-time claim is made.
 """
 Plo,Phi,F,Q,H,R=map(np.asarray,(Plo,Phi,F,Q,H,R))
 sym=lambda a:(a+a.T)/2
 if np.linalg.eigvalsh(sym(Plo)).min() < -1e-12:
  raise ValueError("PSD lower covariance required")
 if np.linalg.eigvalsh(sym(Phi-Plo)).min() < -1e-12:
  raise ValueError("ordered covariance faces required")
 if np.linalg.eigvalsh(sym(Q)).min() < -1e-12:
  raise ValueError("PSD process covariance required")
 rmin=float(np.linalg.eigvalsh(sym(R)).min())
 if rmin<=0:raise ValueError("positive measurement noise required")
 prior=F@Phi@F.T;C_hi=prior+Q
 norm=lambda a:float(np.linalg.norm(a,2))
 t1=min(norm(prior),norm(prior)*norm(H@Phi@H.T)/rmin)
 t2=min(norm(C_hi),norm(C_hi)*norm(H@C_hi@H.T)/rmin)
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
