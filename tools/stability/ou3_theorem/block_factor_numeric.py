"""Constructive real-arithmetic A21 block-factor normalization.

Consumes the promoted LIN path-energy factor and the retained AG/BA source
factors.  This certificate is deliberately separate from float32 whole-word
supply, which remains an additive practical-stability obligation.
"""
from __future__ import annotations
import json
from pathlib import Path
from .lin_path_certificate import certificate as lin_certificate

MU_N=2.04e-3
ELL_AG=4.999974644e-4
ELL_BA=5.618273739e-4

def certificate() -> dict:
    lin=lin_certificate()
    required=(
        lin["verified"],
        lin["deployment_noise_floor_binding_verified"],
        lin["shipping_interleaved_corrections_verified"],
        lin["aw_sync_psd_inflation_verified"],
    )
    if not all(required):
        return {"verified":False,"mu_cov":0.0,"rho0":1.0,
                "failure":"LIN shipping factor not promoted"}
    ell_lin=float(lin["factor_floor"])
    # Fixed-coordinate J >= mu_N I and P >= D D' imply
    # D' J D >= mu_N D'D.  The scalar floor is conservative but valid in
    # the declared proof coordinates; root coercivity gamma_root=1 comes from
    # the additive block-diagonal process injection.
    p_factor=min(ELL_AG*ELL_AG,ell_lin*ell_lin,ELL_BA*ELL_BA)
    mu_cov=MU_N*p_factor
    rho0=1.0/(1.0+mu_cov)
    return {
        "qualification":"OU3_A21_BLOCK_FACTOR_NUMERIC_REAL_ARITHMETIC_V1",
        "verified":mu_cov>0.0 and rho0<1.0,
        "mu_N":MU_N,
        "ell_AG":ELL_AG,"ell_LIN":ell_lin,"ell_BA":ELL_BA,
        "root_metric_gamma":1.0,
        "factor_covariance_floor":p_factor,
        "mu_cov":mu_cov,"rho0":rho0,
        "sqrt_rho0_margin":1.0-rho0**0.5,
        "cross_block_role":"finite prefix magnitude/totality; not root coercivity",
        "float32_whole_word_supply_closed":False,
        "theorem_closed":False,
    }

def main() -> int:
    r=certificate()
    print(json.dumps(r,indent=2,sort_keys=True))
    return 0 if r["verified"] else 1

if __name__=="__main__":
    raise SystemExit(main())
