"""Exact necessity result for loss of informative heading service."""
from __future__ import annotations
import math

def centre_block(T_s: float) -> tuple[tuple[float,float],tuple[float,float]]:
    if not (math.isfinite(T_s) and T_s>0.0): raise ValueError("positive finite interval required")
    return ((1.0,T_s),(0.0,1.0))

def repeated_state(theta0: float,axial_bias_error: float,T_s: float,n: int) -> tuple[float,float]:
    if n<0: raise ValueError("nonnegative repeat count required")
    if not all(math.isfinite(v) for v in (theta0,axial_bias_error,T_s)) or T_s<=0.0:
        raise ValueError("finite state and positive interval required")
    return theta0+n*T_s*axial_bias_error,axial_bias_error

def obstruction_report(T_s: float=3.0) -> dict:
    A=centre_block(T_s)
    return {"block":[list(A[0]),list(A[1])],"eigenvalues":[1.0,1.0],"spectral_radius":1.0,
            "strict_full_state_contraction_possible_without_heading_information":False,
            "linear_heading_growth_for_nonzero_axial_bias_error":True,
            "role":"necessity/limitation result only"}
