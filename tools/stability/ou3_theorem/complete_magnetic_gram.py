"""Complete two-word Gram reduction under restricted magnetic service.

Arbitrary PSD magnetic nuisance can cancel protected service curvature along
the compatibility cone.  The source-uniform comparison must therefore be
made with the nonmagnetic joint Gram restricted to that cone.
"""
from __future__ import annotations
import math
from .interval_riccati_21 import IMat,matmul,transpose
from .rank_loss_interval_factor import generalized_ratio_lower

def scalar_magnetic_cancellation(h:float,n:float):
    """Construct rank-one service row cancelling one selected (h,n) direction."""
    if n==0:
        return {"cancellable":False,"minimum_service_action":h*h}
    a=1.0
    b=-h/n
    # Add an orthogonal protected row in the real 2-D service block as needed.
    return {"cancellable":True,"row":[a,b],
            "selected_direction_action":(a*h+b*n)**2}

def restricted_complete_ratio(R:IMat,Q_nonmag:IMat,N_terminal:IMat):
    """Conservative lower action/terminal ratio on compatibility coordinates."""
    qc=matmul(matmul(transpose(R),Q_nonmag),R)
    nc=matmul(matmul(transpose(R),N_terminal),R)
    z=generalized_ratio_lower(qc,nc)
    z["qualification"]="OU3_COMPLETE_MAG_CANCEL_REDUCTION_V1"
    return z

def current_contract():
    return {"magnetic_only_floor_required":False,
            "magnetic_relative_process_floor_required":False,
            "complete_compatibility_ratio_verified":False,
            "reason":"moving/rank-changing compatibility basis R_W not yet source-uniformly enclosed"}
