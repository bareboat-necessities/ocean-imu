"""Exact parity decomposition for the planar MOVING OU-III history.

The state order is [theta,bg,v,p,S,aw,ba], each xyz.  Under the common
Ry(-psi) record, literal prediction and acc/mag/S rows preserve two coordinate
parities.  This module records the lossless permutation and the corresponding
sum/difference service-probe split.
"""
from __future__ import annotations
import numpy as np

# XZ-residual parity: theta_y,bg_y plus x/z LIN/AW/BA.
EVEN=(1,4, 6,8, 9,11, 12,14, 15,17, 18,20)
# Y-residual parity: theta_xz,bg_xz plus y LIN/AW/BA.
ODD=(0,2,3,5, 7,10,13,16,19)
assert sorted(EVEN+ODD)==list(range(21))

def permutation():
    p=np.zeros((21,21))
    for i,j in enumerate(EVEN+ODD): p[i,j]=1.
    return p

def service_probe_transform():
    """Orthogonal +/- transform on [h+,bg+,h-,bg-].

    Difference probes carry the y component (EVEN); sums carry the xz component
    (ODD). Scaling is orthogonal, so information eigenvalues are unchanged.
    """
    s=2**-.5
    return np.array([[s,0,-s,0],[0,s,0,-s],[s,0,s,0],[0,s,0,s]])

def off_parity_norm(a):
    a=np.asarray(a,float)
    return float(np.linalg.norm(a[np.ix_(EVEN,ODD)]))

def split_information(info):
    t=service_probe_transform(); z=t@np.asarray(info,float)@t.T
    return z[:2,:2],z[2:,2:],float(np.linalg.norm(z[:2,2:]))

def certificate():
    return {"qualification":"OU3_PLANAR_PARITY_DECOMPOSITION_V1",
            "even_indices":list(EVEN),"odd_indices":list(ODD),
            "block_dimensions":[len(EVEN),len(ODD)],
            "auxiliary_transformed_blocks":[2,2],
            "literal_service_criterion":"two separate 2x2 heading/BG Gramians for d_plus and d_minus",
            "lossless":True,
            "basis_change_orthogonal":True,
            "claim":"literal planar real-arithmetic factors preserve parity; interval enclosure may propagate the two blocks independently",
            "float32_transfer_open":True}

if __name__=="__main__":
 import json;print(json.dumps(certificate(),indent=2,sort_keys=True))
