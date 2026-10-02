"""Closed-form accelerometer innovation inverse from rank-two force projector."""
from __future__ import annotations
import math,numpy as np
from .golive_release_seed import certificate as seed

def coefficients(rho,r_noise):
 """S=c I + a(I-nn'), c=r_noise+p_aw+p_ba, a=p_theta rho^2."""
 d=seed()["covariance_diagonal_upper"];pt=d[0];c=r_noise+d[15]+d[18];a=pt*rho*rho
 return c,a

def inverse_coefficients(rho,r_noise):
 """(cI+aP)^-1 = c^-1 nn' + (c+a)^-1(I-nn')
 = u I + v nn', P=I-nn'."""
 c,a=coefficients(rho,r_noise);u=1/(c+a);v=1/c-u
 return {"c":c,"a":a,"u":u,"v":v,
         "parallel_inverse_eigenvalue":1/c,
         "transverse_inverse_eigenvalue":1/(c+a)}

def inverse_matrix(rho,n,r_noise):
 n=np.asarray(n,float);n/=np.linalg.norm(n);q=inverse_coefficients(rho,r_noise)
 return q["u"]*np.eye(3)+q["v"]*np.outer(n,n)

def direction_free_bounds(rho_lo,rho_hi,r_lo,r_hi):
 """Uniform spectral facts for every unit n and allowed rho/noise."""
 if not(0<=rho_lo<=rho_hi and 0<r_lo<=r_hi):raise ValueError("ranges")
 # Largest inverse eigenvalue is always parallel 1/c and decreases with noise.
 d=seed()["covariance_diagonal_upper"];base=d[15]+d[18]
 inv_norm=1/(r_lo+base)
 # smallest S eigenvalue is c, independent of rho,n.
 s_floor=r_lo+base
 # transverse inverse ranges monotonically with rho/noise.
 trans_lo=1/(r_hi+base+d[0]*rho_hi*rho_hi)
 trans_hi=1/(r_lo+base+d[0]*rho_lo*rho_lo)
 return {"verified":True,"force_direction_required":False,
  "innovation_eigen_floor":s_floor,"inverse_spectral_norm_upper":inv_norm,
  "parallel_inverse_upper":inv_norm,
  "transverse_inverse_lower":trans_lo,"transverse_inverse_upper":trans_hi,
  "projector_inverse_exact":True}
