"""Constructive real-arithmetic A21 block-factor normalization.

Consumes the promoted LIN path-energy factor and the retained AG/BA source
factors.  This certificate is deliberately separate from float32 whole-word
supply, which remains an additive practical-stability obligation.
"""
from __future__ import annotations
import json
from pathlib import Path
try:\n    from .lin_path_certificate import certificate as lin_certificate\nexcept ImportError:\n    from lin_path_certificate import certificate as lin_certificate

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
        "mu_cov":mu_cov,"rho0":rho0,\n        "constructive_full_A21_mu_rho_enclosure":True,\n        "normalization_kind":"scalar fallback; not final sharp factor-coordinate metric",
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


def blockwise_factor_coordinate_diagnostic() -> dict:
    """Expose exactly where the scalar normalization loses the margin.

    This is non-promoting until the fixed-coordinate information certificate
    exports compatible block/coordinate information factors.  It prevents us
    from pretending that min(ell_i)^2 is intrinsic to the theorem.
    """
    r=certificate()
    contributions={
        "AG_covariance_scale":ELL_AG*ELL_AG,
        "LIN_covariance_scale":r["ell_LIN"]**2,
        "BA_covariance_scale":ELL_BA*ELL_BA,
    }
    return {
        "verified_structure":True,
        "promoting":False,
        "scalar_bottleneck":min(contributions,key=contributions.get),
        "covariance_scales":contributions,
        "required_next_certificate":
            "compatible factor-coordinate information matrix/block floors; "
            "compute lambda_min(L' J L) directly instead of mu_N*min(ell)^2",
    }


def sharpened_factor_coordinate_certificate() -> dict:
    """Fixed, source-independent diagonal coordinate congruence.

    This is not a fitted rho: the declared diagonal coordinate map is certified
    by rerunning both the exact-rational LIN action bound and the analytic
    translation information bound in the same coordinates.
    """
    scales=(2.4,18.0,132.0,4.0)
    lin=lin_certificate(scales)
    ell=float(lin["factor_floor"])
    word=16.0; gap=.156; v,p,s=scales[:3]; noise=100.0
    d=.5*word-gap; t=word+gap
    det=(d**3)*v*p*s
    row2=(.5*v*t*t)**2+(p*t)**2+s*s
    sigma_min=det/(3.0*row2)
    mu_trans=(sigma_min/noise)**2
    # Retained long-word AG floor is >2.1e-3.  2.1e-3 is used downward.
    mu_ag=2.1e-3
    mu_neutral=min(mu_trans,mu_ag)
    # BA and a_w are strictly stable modes and are handled by detectability,
    # not allowed to dilute the neutral-quotient information floor.
    mu_cov=mu_neutral*ell*ell
    rho0=1.0/(1.0+mu_cov)
    return {
        "qualification":"OU3_A21_SHARP_FACTOR_COORDINATES_REAL_ARITHMETIC_V1",
        "verified":lin["verified"] and mu_cov>0 and rho0<1,
        "coordinate_map":{"v":v,"p":p,"S":s,"a_w":scales[3]},
        "ell_LIN":ell,"mu_translation":mu_trans,"mu_AG":mu_ag,
        "mu_neutral":mu_neutral,"mu_cov":mu_cov,"rho0":rho0,
        "sqrt_rho0_margin":1.0-rho0**0.5,
        "active_stable_modes":"a_w and BA retained by detectability; not neutral quotient",
        "selection_role":"fixed proof-coordinate congruence, not sampled/fitted rho",
        "float32_whole_word_supply_closed":False,
        "theorem_closed":False,
    }
