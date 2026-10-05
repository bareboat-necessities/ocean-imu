"""Linked correction calculus, with noise, held-gain and reset variations.

These identities concern the real operation lift on a fixed regular branch.
They do not differentiate a rounded binary32 program, assume a future orbit,
or supply a uniform cell. See the manuscript for the exact finite remainder.
"""
from __future__ import annotations

from .matrix_certificates import add, identity, matmul, transpose


def product(*xs):
    out = xs[0]
    for x in xs[1:]:
        out = matmul(out, x)
    return out


def innovation_differential(p, h, dp, dh, dr):
    return add(add(add(product(dh, p, transpose(h)),
                       product(h, dp, transpose(h))),
                   product(h, p, transpose(dh))), dr)


def gain_differential(p, h, k, sinv, dp, dh, dr):
    """Unmasked gain; sinv is the SAME operation's inverse innovation."""
    ds = innovation_differential(p, h, dp, dh, dr)
    db = add(product(dp, transpose(h)), product(p, transpose(dh)))
    return product(add(db, product(k, ds), -1), sinv)


def optimal_covariance_differential(p, h, k, dp, dh, dr):
    """Unmasked regular branch only; the held-BA branch uses joseph below."""
    a = add(identity(len(p)), product(k, h), -1)
    c = product(a, p)
    return add(add(add(product(a, dp, transpose(a)),
                       product(k, dh, c), -1),
                   product(c, transpose(dh), transpose(k)), -1),
               product(k, dr, transpose(k)))


def joseph_differential(k, b, s, dp, dk, db, ds):
    """Literal P-K B'-B K'+K S K', including a held-BA-masked B.

    Here B is the literal PCt operand, not an independently chosen matrix.
    If B=M P H_gain', dB=M(dP H_gain'+P dH_gain') on a fixed hold branch.
    Held accelerometer H_gain omits the BA columns, while H_S retains them.
    The implementation's LDLT diagonal bump, when taken, is part of S and dS.
    """
    out = dp
    for term in (product(dk, transpose(b)), product(k, transpose(db)),
                 product(db, transpose(k)), product(b, transpose(dk))):
        out = add(out, term, -1)
    for term in (product(dk, s, transpose(k)), product(k, ds, transpose(k)),
                 product(k, s, transpose(dk))):
        out = add(out, term)
    return out


def gain_finite_remainder(p, h, k, sinv, new_sinv, dp, dh, dr):
    """Exact DeltaK-DK for finite linked increments, not a Taylor estimate."""
    s1 = innovation_differential(p, h, dp, dh, dr)
    s2 = add(add(add(product(dh, dp, transpose(h)),
                       product(h, dp, transpose(dh))),
                   product(dh, p, transpose(dh))),
               product(dh, dp, transpose(dh)))
    dk = gain_differential(p, h, k, sinv, dp, dh, dr)
    numerator = add(add(product(dp, transpose(dh)), product(k, s2), -1),
                    product(dk, add(s1, s2)), -1)
    return product(numerator, new_sinv)


def reset_differential(g, c, dg, dc):
    return add(add(product(g, dc, transpose(g)), product(dg, c, transpose(g))),
               product(g, c, transpose(dg)))


def certificate():
    return {
        "qualification": "OU3_PLANAR_LINKED_RICCATI_MEAN_DIFFERENTIAL_V2",
        "result_type": "PROVED — analytical",
        "scope": "local real-operation identities on a fixed regular branch; no uniform orbit bound",
        "innovation_differential": "dS=dH P H'+H dP H'+H P dH'+dR",
        "gain_differential": "dK=(dP H'+P dH'-K dS) S^-1",
        "covariance_differential": "dC=A dP A'-K dH C-C dH' K'+K dR K'; A=I-KH; C=AP",
        "held_ba_rule": "use B=M P H_gain', dB=M(dP H_gain'+P dH_gain'); held accelerometer H_gain omits BA, H_S retains BA; do not use the unmasked optimal identity",
        "mean_differential": "d(x+)=d(x)+dK*r+K*dr; dr is differentiated from the same measurement model",
        "finite_gain_remainder": "DeltaK-DK=(DeltaP DeltaH'-K S2-DK DeltaS)(S+DeltaS)^-1",
        "S2": "DeltaH DeltaP H'+H DeltaP DeltaH'+DeltaH P DeltaH'+DeltaH DeltaP DeltaH'",
        "reset_differential": "d(G C G')=G dC G'+dG C G'+G C dG'; dG_theta=0.5[d(injected theta)]x",
        "noise_is_same_history": True,
        "independent_extrema_used": False,
        "structures_preserved": ["actual P,H,R,S,K,r", "held-BA PCt mask", "Joseph", "reset"],
        "relaxations_introduced": ["real operation lift; arithmetic charge separate", "fixed branch for differentials"],
        "finite_word_uniform_bound_closed": False,
        "theorem_closed": False,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
