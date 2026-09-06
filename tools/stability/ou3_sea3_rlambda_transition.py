#!/usr/bin/env python3
"""Machine-readable conservative transition relation for compact SEA3 lambda.

The SEA3 theorem declares a compact, rate-bounded relation

    lambda_{k+1} in R_lambda(lambda_k)

but the manuscript does not currently give numerical per-component rate
constants.  Inventing such constants would narrow the theorem domain.  For the
canonical proof we therefore use the stronger, conservative outer relation

    Rhat_lambda(lambda) = Lambda_SEA3,

where *both* endpoints must satisfy the full coupled compact SEA3 parameter
constraints.  Every theorem-admissible rate-bounded transition is contained in
this relation, so any certificate proved over Rhat_lambda automatically covers
the unknown tighter R_lambda.  Missing rate data can only reduce conservatism;
it is not needed for soundness.

This is not an independent rectangular parameter box.  The three partition
heights retain the exact total-energy coupling and each active partition keeps
the H/Tp peak-steepness constraint.  Frequency and spreading are represented
in compactified coordinates rather than by invented finite caps:

    nu = f_p / (1 + f_p) in [0,1],
    chi = s / (1 + s)   in [0,1].

The closure points nu=0 and chi=1 are conservative compactification limits.
For an active partition nu=1 is rejected by the steepness relation; inactive
H=0 slots may carry arbitrary compact parameter memory without changing the
sea.  Mean direction is represented by a turn coordinate beta_turn in [0,1]
with the endpoints identified physically.

This module closes only the machine-readable R_lambda ingredient of SEA0.  It
does not generate x^s, a response trajectory, tuner history, or a P3 word.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Sequence

import ou3_sea3_physical_admissibility as PHYSICAL
import ou3_sea3_spectral_moment_bridge as MOMENT

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_COUPLED_COMPACT_RLAMBDA_OUTER_RELATION_V1"
M_MAX = 3


@dataclass(frozen=True)
class CompactPartition:
    H_m: float
    peak_frequency_compact: float
    gamma: float
    mean_direction_turn: float
    spreading_compact: float


def _finite_unit(x: float) -> bool:
    y = float(x)
    return math.isfinite(y) and 0.0 <= y <= 1.0


def peak_frequency_hz_from_compact(nu: float) -> float:
    """Decode nu=f/(1+f); nu=1 is the compact f=+infinity boundary."""
    x = float(nu)
    if not _finite_unit(x):
        raise ValueError("peak-frequency compact coordinate must lie in [0,1]")
    if x == 1.0:
        return math.inf
    return x / (1.0 - x)


def spreading_from_compact(chi: float) -> float:
    """Decode chi=s/(1+s); chi=1 is the compact delta-spreading limit."""
    x = float(chi)
    if not _finite_unit(x):
        raise ValueError("spreading compact coordinate must lie in [0,1]")
    if x == 1.0:
        return math.inf
    return x / (1.0 - x)


def _partition_member(p: CompactPartition, gravity: float, gamma_range: Sequence[float]) -> bool:
    h = float(p.H_m)
    if not (math.isfinite(h) and h >= 0.0):
        return False
    if not _finite_unit(p.peak_frequency_compact):
        return False
    if not _finite_unit(p.mean_direction_turn):
        return False
    if not _finite_unit(p.spreading_compact):
        return False
    if len(gamma_range) != 2:
        return False
    g0, g1 = map(float, gamma_range)
    gamma = float(p.gamma)
    if not (math.isfinite(gamma) and g0 <= gamma <= g1):
        return False

    # Inactive slots carry compact parameter memory but inject zero sea energy.
    if h == 0.0:
        return True

    fp = peak_frequency_hz_from_compact(p.peak_frequency_compact)
    if not (math.isfinite(fp) and fp > 0.0):
        # nu=0 (infinite period) is retained only as a closure point; an active
        # finite-energy JONSWAP partition has positive finite peak frequency.
        return False
    tp = 1.0 / fp
    return PHYSICAL.partition_admissible(h, tp, gravity)


def lambda_member(
    partitions: Sequence[CompactPartition],
    *,
    gravity: float,
    Hs_upper_m: float,
    gamma_range: Sequence[float],
) -> bool:
    """Membership in the coupled compact SEA3 endpoint set Lambda_SEA3."""
    if len(partitions) != M_MAX:
        return False
    if not (math.isfinite(gravity) and gravity > 0.0):
        return False
    if not (math.isfinite(Hs_upper_m) and Hs_upper_m > 0.0):
        return False
    if not all(_partition_member(p, gravity, gamma_range) for p in partitions):
        return False
    energy = sum(float(p.H_m) ** 2 for p in partitions)
    return energy <= Hs_upper_m * Hs_upper_m


def transition_admissible(
    lambda_in: Sequence[CompactPartition],
    lambda_out: Sequence[CompactPartition],
    *,
    gravity: float,
    Hs_upper_m: float,
    gamma_range: Sequence[float],
) -> bool:
    """Conservative Rhat relation: both endpoints belong to coupled SEA3."""
    return lambda_member(
        lambda_in,
        gravity=gravity,
        Hs_upper_m=Hs_upper_m,
        gamma_range=gamma_range,
    ) and lambda_member(
        lambda_out,
        gravity=gravity,
        Hs_upper_m=Hs_upper_m,
        gamma_range=gamma_range,
    )


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    physical = PHYSICAL.build(Path(domain_path).resolve())
    pf = PHYSICAL.validate(physical)
    moment = MOMENT.build()
    mf = MOMENT.validate(moment)
    if pf or mf:
        raise RuntimeError(f"R_lambda prerequisites failed: physical={pf}, moment={mf}")

    gamma_range = list(map(float, moment["sea_family"]["declared_gamma_interval"]))
    hs_upper = float(physical["repository_total_Hs_upper_m"])
    gravity = float(physical["gravity_mps2"])

    # Nontrivial coupled smoke witnesses.  The first uses two active partitions
    # whose individual heights are below the total cap but whose vector energy
    # is checked jointly.  The second moves all compact coordinates while
    # remaining inside the same coupled theorem endpoint set.
    lam0 = (
        CompactPartition(3.0, 0.15 / 1.15, 1.0, 0.10, 0.25),
        CompactPartition(2.0, 0.1 / 1.1, 3.3, 0.65, 0.50),
        CompactPartition(0.0, 1.0, 7.0, 1.0, 1.0),
    )
    lam1 = (
        CompactPartition(2.5, 0.18 / 1.18, 2.0, 0.95, 0.35),
        CompactPartition(2.4, 0.09 / 1.09, 6.5, 0.15, 0.80),
        CompactPartition(0.5, 0.07 / 1.07, 1.0, 0.45, 1.0),
    )
    smoke_in = lambda_member(
        lam0, gravity=gravity, Hs_upper_m=hs_upper, gamma_range=gamma_range
    )
    smoke_out = lambda_member(
        lam1, gravity=gravity, Hs_upper_m=hs_upper, gamma_range=gamma_range
    )
    smoke_transition = transition_admissible(
        lam0,
        lam1,
        gravity=gravity,
        Hs_upper_m=hs_upper,
        gamma_range=gamma_range,
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generator": False,
        "trajectory_replay_used": False,
        "rate_constants_fitted_or_invented": False,
        "exact_unknown_rate_relation_claimed": False,
        "relation_role": "sound conservative outer relation for theorem R_lambda",
        "actual_rate_bounded_R_lambda_subset_Rhat": True,
        "future_rate_limits_only_tighten_relation": True,
        "machine_readable_R_lambda_closed": True,
        "fixed_lambda_word_used": False,
        "independent_component_rectangle_used": False,
        "coupled_partition_energy_retained": True,
        "coupled_peak_steepness_retained": True,
        "three_fixed_partition_slots_retained": True,
        "compact_coordinates": {
            "peak_frequency": "nu=f_p/(1+f_p) in [0,1]",
            "mean_direction": "beta_turn in [0,1] with endpoints physically identified",
            "spreading": "chi=s/(1+s) in [0,1]",
            "gamma_interval": gamma_range,
            "partition_height": "H_r>=0 with sum H_r^2 <= Hs_max^2",
            "Hs_upper_m": hs_upper,
            "active_partition_constraint": physical["three_partition_contract"][
                "active_partition_constraint"
            ],
        },
        "outer_relation": {
            "definition": "Rhat_lambda(lambda)=Lambda_SEA3 for every lambda in Lambda_SEA3",
            "both_endpoints_must_satisfy_full_coupled_membership": True,
            "sample_to_sample_parameter_jump_restricted_by_unknown_rate_numbers": False,
            "soundness_reason": (
                "the theorem's declared rate-bounded R_lambda is a subset of the full coupled compact endpoint relation; proving the latter covers every admissible rate-bounded history"
            ),
        },
        "smoke": {
            "lambda_in_member": smoke_in,
            "lambda_out_member": smoke_out,
            "transition_member": smoke_transition,
        },
        "hard_shaping_state_materialized_here": False,
        "joint_response_materialized_here": False,
        "P3_promoted": False,
        "next_obligation": (
            "propagate a hard phase-continuous x^s/equivalent finite-window constraint and the coupled translational/rotational response over this complete conservative lambda relation"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "actual_rate_bounded_R_lambda_subset_Rhat",
        "future_rate_limits_only_tighten_relation",
        "machine_readable_R_lambda_closed",
        "coupled_partition_energy_retained",
        "coupled_peak_steepness_retained",
        "three_fixed_partition_slots_retained",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_generator",
        "trajectory_replay_used",
        "rate_constants_fitted_or_invented",
        "exact_unknown_rate_relation_claimed",
        "fixed_lambda_word_used",
        "independent_component_rectangle_used",
        "hard_shaping_state_materialized_here",
        "joint_response_materialized_here",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    relation = d.get("outer_relation", {})
    if relation.get("both_endpoints_must_satisfy_full_coupled_membership") is not True:
        f.append("outer relation lost coupled endpoint membership")
    if relation.get("sample_to_sample_parameter_jump_restricted_by_unknown_rate_numbers") is not False:
        f.append("outer relation invented missing numerical rate limits")
    smoke = d.get("smoke", {})
    for key in ("lambda_in_member", "lambda_out_member", "transition_member"):
        if smoke.get(key) is not True:
            f.append(f"R_lambda smoke failed {key}")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "machine_readable_R_lambda_closed": d["machine_readable_R_lambda_closed"],
        "smoke": d["smoke"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())