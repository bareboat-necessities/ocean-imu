"""Stream anisotropic lower factors through a conservative literal event word.

This is a feasibility computation using matrix information ceilings and literal
default process factors. It is not promoted until the all-time f<=30 and exact
tuner/profile bindings are certified. It intentionally retains the full 21x21
matrix lower factor, then reports parity generalized prediction eta.
"""
from __future__ import annotations
import numpy as np, math
from .planar_anisotropic_factors import full_root_lower
from .planar_information_ceilings import acc_ceiling,S_ceiling,mag_ceiling
from .planar_parity import EVEN,ODD
from .lower_factor_transport import corrected_lower

def embed(Je,Jo):
 J=np.zeros((21,21));J[np.ix_(EVEN,EVEN)]=Je;J[np.ix_(ODD,ODD)]=Jo;return J

def process(dt=.005,tau=10.,sigma=.05,tau_b=5000.,q_ba=2.5e-7,sg=.00135,sbg=1e-5):
 # Conservative simplified literal block process for feasibility only.
 F=np.eye(21);Q=np.zeros((21,21))
 # theta-bg one-axis exact constant-rate transition/process.
 for a in range(3):
  th=a;bg=3+a;F[th,bg]=dt
  Q[th,th]=sg*sg*dt+sbg*sbg*dt**3/3;Q[th,bg]=Q[bg,th]=sbg*sbg*dt**2/2;Q[bg,bg]=sbg*sbg*dt
 # LIN exact transition is omitted here => not a certificate; use identity with
 # positive AW drive as feasibility. This makes the diagnostic explicitly non-promoting.
 phi=math.exp(-dt/tau)
 for a in range(3):
  aw=15+a;F[aw,aw]=phi;Q[aw,aw]=sigma*sigma*(1-phi*phi)
 # BA
 phib=math.exp(-dt/tau_b)
 for a in range(3):
  ba=18+a;F[ba,ba]=phib;Q[ba,ba]=q_ba*tau_b*.5*(1-phib*phib)
 return F,Q

def run():
 P=full_root_lower();ae,ao=acc_ceiling();se,so=S_ceiling();me,mo=mag_ceiling()
 Ja,Js,Jm=embed(ae,ao),embed(se,so),embed(me,mo)
 etas=[]; min_eig=1e99
 for k in range(200):
  F,Q=process();C=F@P@F.T
  # generalized eta via whitening current matrix lower bound
  L=np.linalg.cholesky(C);Li=np.linalg.inv(L);W=Li@Q@Li.T
  etas.append(max(0.,float(np.linalg.eigvalsh((W+W.T)/2).max())))
  P=C+Q
  # conservative event word: S at maximum plausible ~ every 24 samples,
  # accel every sample, mag every 8. More corrections shrink lower P.
  if k%24==0:P=np.linalg.inv(np.linalg.inv(P)+Js)
  P=np.linalg.inv(np.linalg.inv(P)+Ja)
  if (k+1)%8==0:P=np.linalg.inv(np.linalg.inv(P)+Jm)
  min_eig=min(min_eig,float(np.linalg.eigvalsh(P).min()))
 retention=float(np.prod([1/(1+x) for x in etas]))
 return {"qualification":"OU3_PLANAR_PREFIX_FACTOR_FEASIBILITY_V1",
         "result_type":"non-promoting feasibility; LIN F/Q and exact scheduler chronology not yet literal",
         "eta_max":max(etas),"retention_product":retention,"posterior_min_eigenvalue":min_eig,
         "information_ceilings":"literal-structure conservative matrices",
         "literal_prediction_sequence_verified":False,"theorem_closed":False}
if __name__=="__main__":
 import json;print(json.dumps(run(),indent=2,sort_keys=True))
