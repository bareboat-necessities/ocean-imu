"""Exact-rational one-second planar attitude/BG oracle certificate.

This proves a self-containing 2x2 covariance cell and magnetic incremental
information > I for an intentionally more-informative direct-attitude oracle.
It does NOT by itself prove that the literal accel+S nuisance chronology is
dominated by this oracle or that |f_hat|<=30 for all time.
"""
from __future__ import annotations
from fractions import Fraction as F
import json

def _mm(a,b):
 return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
def _tr(a): return [list(x) for x in zip(*a)]

def certificate():
 dt=F(1,200);B=F(75);ha=F(30);Rm=F(16,25);Ra=F(1,25)
 qg=F(27,20000)**2;qbg=F(1,10_000_000_000)
 P=[[F(1,2000),F(0)],[F(0),F(11,10_000_000)]];P0=[r[:] for r in P]
 X=[[F(1),F(0)],[F(0),F(1,50)]];I=[[F(0),F(0)],[F(0),F(0)]]
 def meas(P,X,h,R,collect):
  S=h*h*P[0][0]+R;y=[h*X[0][j] for j in range(2)]
  if collect:
   for i in range(2):
    for j in range(2): I[i][j]+=y[i]*y[j]/S
  K=[P[i][0]*h/S for i in range(2)]
  A=[[1-K[0]*h,F(0)],[-K[1]*h,F(1)]]
  Xn=_mm(A,X)
  Pn=[[P[i][j]-(P[i][0]*h)*(h*P[0][j])/S for j in range(2)] for i in range(2)]
  return Pn,Xn
 for k in range(200):
  A=[[F(1),dt],[F(0),F(1)]]
  Q=[[qg*dt+qbg*dt**3/F(3),qbg*dt**2/F(2)],[qbg*dt**2/F(2),qbg*dt]]
  Pp=_mm(_mm(A,P),_tr(A));P=[[Pp[i][j]+Q[i][j] for j in range(2)] for i in range(2)]
  X=_mm(A,X);P,X=meas(P,X,ha,Ra,False)
  if (k+1)%8==0:P,X=meas(P,X,B,Rm,True)
 D=[[P0[i][j]-P[i][j] for j in range(2)] for i in range(2)]
 detD=D[0][0]*D[1][1]-D[0][1]*D[1][0]
 J=[[I[i][j]-F(1 if i==j else 0) for j in range(2)] for i in range(2)]
 detJ=J[0][0]*J[1][1]-J[0][1]*J[1][0]
 invariant=(D[0][0]>0 and D[1][1]>0 and detD>0)
 service=(J[0][0]>0 and J[1][1]>0 and detJ>0)
 return {
  "qualification":"OU3_PLANAR_ORACLE_RATIONAL_V1",
  "arithmetic":"exact fractions; no floating comparison used for proof signs",
  "oracle_parameters":{"dt":"1/200","acc_attitude_gain":30,"R_acc":"1/25","B":75,"R_mag":"16/25",
                       "P_theta0":"1/2000","P_bg0":"11/10000000","gyro_density":"27/20000","bg_rw":"1/10000000000"},
  "one_second_covariance_self_inclusion_exact":invariant,
  "magnetic_information_minus_identity_spd_exact":service,
  "covariance_difference_diag_decimal":[float(D[0][0]),float(D[1][1])],
  "covariance_difference_det_decimal":float(detD),
  "information_matrix_decimal":[[float(x) for x in r] for r in I],
  "information_minus_identity_det_decimal":float(detJ),
  "literal_shipping_service_proved":False,
  "open_dependencies":["prove literal nonmag chronology is no more suppressive than oracle for magnetic incremental information",
                       "prove all-time literal predicted-force/attitude-information bound corresponding to oracle gain <=30",
                       "carry reset-coordinate effect in the comparison",
                       "apply comparison to both literal +/- service pairs and every placed root"],
  "theorem_closed":False,
 }

if __name__=="__main__": print(json.dumps(certificate(),indent=2,sort_keys=True))
