"""Scalar terminal compatibility controllability energy.

For one compatibility root direction r, compute the known-root source energy
needed to reproduce its terminal component transverse to the next kernel while
remaining data-null.  This is the exact Pi^dagger quadratic and avoids a
full covariance upper bound.
"""
from __future__ import annotations
import math

def scalar_terminal_bound(E_perp:float,E_kernel:float,c_current:float,c_next:float):
    if min(E_perp,E_kernel)<0 or c_current<=0 or c_next<=0:
        raise ValueError("nonnegative energies and positive ceilings required")
    B=E_perp+E_kernel/c_next
    return {"B_term":B,"C_c_upper":c_current*B}

def source_range_status(in_range:bool,energy_upper:float|None=None):
    if not in_range:
        return {"verified":False,"energy_upper":math.inf,
                "reason":"terminal transverse compatibility image outside data-null source range"}
    if energy_upper is None or energy_upper<0 or not math.isfinite(energy_upper):
        return {"verified":False,"energy_upper":math.inf,
                "reason":"finite data-null source energy not certified"}
    return {"verified":True,"energy_upper":energy_upper}

def current_contract():
    return {"qualification":"OU3_TERMINAL_COMPAT_CONTROLLABILITY_V1",
            "data_null_source_range_inclusion_verified":False,
            "E_perp_source_uniform_upper":math.inf,
            "E_kernel_source_uniform_upper":math.inf,
            "B_term_source_uniform_upper":math.inf,
            "next_obligation":"bound scalar data-null source controllability energy and 1-D deterministic kernel transport"}
