"""Diagnostics for exact common-R goLive innovation Gram."""
from __future__ import annotations
import math,numpy as np
from .blockwise_accel_gram import force_cone_gram_bounds,innovation_inverse_from_gram
from .release_interval_propagation import one_sample_boxes
RHO=18.60665;AX=np.array([0.,0.,1.])
def residual(alpha_deg,drho):
 _,_,_,R,*_=one_sample_boxes()
 G,Gr,m=force_cone_gram_bounds(max(0,RHO-drho),RHO,AX,math.radians(alpha_deg))
 return {"alpha_deg":alpha_deg,"drho":drho,**m,**innovation_inverse_from_gram(G,Gr,R)}
def diagnostic():
 rows=[residual(a,0) for a in [60,30,15,7.5,3,1,.5,.1,.01,0]]
 radial=[residual(0,d) for d in [18.60665,10,5,1,.1,.01,0]]
 return {"angle":rows,"radial":radial,"attitude_cells_required":False}
