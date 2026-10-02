"""MEKF reset split for the changing-axis INJ covariance proof.

Shipping reset is G=I+0.5[dtheta]x on attitude coordinates.  It is not a
common rotation.  Relative to the current force axis n, split d=d_parallel*n+
d_perp.  The parallel part stays in INJ; the perpendicular part is a carried
matrix defect and must not be discarded.
"""
from __future__ import annotations
import numpy as np
from .inj_block_algebra import INJ

def split_reset_vector(dtheta,n):
 d=np.asarray(dtheta,float);n=np.asarray(n,float);n=n/np.linalg.norm(n)
 dp=float(d@n);perp=d-dp*n
 return dp,perp

def right_reset_inj(block:INJ,dtheta,n):
 """For an AW-theta block: P_aw,th^+ = P_aw,th^- G^T."""
 dp,perp=split_reset_vector(dtheta,n)
 # G^T=I-.5[d]x. Parallel contribution is I-.5*dp*J.
 structured=block.mul(INJ(1.,0.,-.5*dp))
 # Exact residual matrix from the perpendicular component.
 B=block.matrix(n)
 x,y,z=perp
 X=np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])
 defect=-.5*B@X
 return {"structured":structured,"defect":defect,
         "defect_norm":float(np.linalg.norm(defect,2)),
         "d_parallel":dp,"d_perp_norm":float(np.linalg.norm(perp)),
         "exact_split":True}

def left_reset_inj(block:INJ,dtheta,n):
 """For a theta-X block: P_th,X^+ = G P_th,X^-."""
 dp,perp=split_reset_vector(dtheta,n)
 structured=INJ(1.,0.,.5*dp).mul(block)
 B=block.matrix(n)
 x,y,z=perp
 X=np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])
 defect=.5*X@B
 return {"structured":structured,"defect":defect,
         "defect_norm":float(np.linalg.norm(defect,2)),
         "d_parallel":dp,"d_perp_norm":float(np.linalg.norm(perp)),
         "exact_split":True}

def reset_defect_bound(block:INJ,dtheta,n):
 """Sharp submultiplicative bound used only when matrix residual is not carried."""
 from .structured_gain_export import inj_spectral_norm
 _,perp=split_reset_vector(dtheta,n)
 return .5*inj_spectral_norm(block)*float(np.linalg.norm(perp))
