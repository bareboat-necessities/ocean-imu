"""Exact scalar kernel-regularized reader identities."""

from __future__ import annotations
import math

def regularized_scalar(J:float,c:float,B:float):
    if J<0 or c<=0 or B<0: raise ValueError("J,B nonnegative; c positive")
    m=1.0/c
    post_var=1.0/(J+m)
    excess=B*post_var
    return {"posterior_variance":post_var,
            "terminal_excess":excess,
            "upper_cB":c*B,
            "root_residual_factor":m/(J+m),
            "observation_reader_factor":1.0/(J+m)}

def augmented_rank_one_leverage(c:float,B:float|None):
    """c*t'(Pi+c tt')^dagger*t. B=t'Pi^dagger t; None means B=+infinity limit."""
    if c<=0: raise ValueError("c positive")
    if B is None or math.isinf(B):
        return {"verified":True,"leverage_upper":1.0,"limit":True}
    if B<0: raise ValueError("B nonnegative")
    lev=c*B/(1.0+c*B)
    return {"verified":True,"leverage":lev,"leverage_upper":1.0,"limit":False}

def current_contract():
    return {"qualification":"OU3_SCALAR_KERNEL_REG_READER_V1",
            "rank_change_J_to_zero_bounded":True,
            "exact_source_controllability_required":False,
            "current_kernel_augmented_leverage_upper":1.0,
            "remaining":"current terminal rank-one direction versus next-word compatibility line"}
