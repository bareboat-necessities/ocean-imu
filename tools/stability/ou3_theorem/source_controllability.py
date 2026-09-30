"""One-vector data-null source controllability conditions.

These helpers encode the exact logical status: root compatibility does not
imply source-range inclusion.  Exact source mimic is a sufficient route, not
a prerequisite for the kernel-regularized reader.
"""
from __future__ import annotations
import math

def source_mimic_status(range_inclusion:bool,terminal_norm:float,
                        directional_modulus_lower:float|None=None):
    if not range_inclusion:
        return {"verified":False,"E_perp_upper":math.inf,
                "reason":"root data-nullity does not imply data-null source reachability"}
    if directional_modulus_lower is None or directional_modulus_lower<=0:
        return {"verified":False,"E_perp_upper":math.inf,
                "reason":"positive directional source-controllability modulus not certified"}
    return {"verified":True,
            "E_perp_upper":terminal_norm**2/directional_modulus_lower**2}

def obligations():
    return {"range_inclusion":False,"finite_E_perp":False,
            "finite_E_kernel":False,
            "exact_source_mimic_is_only_sufficient":True,
            "preferred_route":"kernel-regularized scalar compatibility reader"}
