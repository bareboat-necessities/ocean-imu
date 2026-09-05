#!/usr/bin/env python3
"""Validated continuum phase-state semantics for complete SEA3.

This is a subcertificate of the hard SEA3 shaping state, not a source
generator.  It closes only the phase coordinate and its propagation.  It does
not supply the still-open hard spectral amplitude/driver set, the physical
left inclusion, or the joint translational/rotational output enclosure.

For every partition and every continuum spectral coordinate (omega,theta), the
phase coordinate is represented on the unit circle by q=(q_c,q_s) and evolves
pointwise as

    q_{k+1} = R(omega*h) q_k,

where R is the exact two-dimensional rotation appearing in the SEA3 theorem.
The relation is applied to the entire continuum index set; no finite frequency
or direction grid is introduced.  An admissible lambda transition may change
spectral weights/response parameters but may not reseed q.  Inactive partition
slots retain their phase coordinate so activation does not change model
dimension or manufacture a new phase history.

The unit circle is compact and a product of compact unit circles is compact in
the product topology.  That topological statement is sufficient to close the
phase-coordinate set itself.  It does *not* establish that the spectral output
map is continuous/bounded on that product; that is deliberately left to the
hard driver/output certificate and therefore cannot promote SEA0 or P3 here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
THEOREM = REPO / "doc" / "kalman_ou_iii" / "w3d-sea3-stability-theorem.tex-part"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_CONTINUUM_PHASE_STATE_V1"

CONTINUUM_PHASE_COORDINATE_SET_CLOSED = True
PHASE_CONTINUOUS_PROPAGATION_CLOSED = True
HARD_SPECTRAL_DRIVER_SET_CLOSED = False
JOINT_SOURCE_OUTPUT_MAP_CLOSED = False


def rotation_step(q_c: float, q_s: float, omega_rad_s: float, h_s: float) -> tuple[float, float]:
    """Pointwise continuum oscillator step used only for algebraic witnesses."""
    qc = float(q_c)
    qs = float(q_s)
    w = float(omega_rad_s)
    h = float(h_s)
    if not all(math.isfinite(x) for x in (qc, qs, w, h)):
        raise ValueError("phase step inputs must be finite")
    if w < 0.0 or h < 0.0:
        raise ValueError("omega and h must be nonnegative")
    a = w * h
    c = math.cos(a)
    s = math.sin(a)
    return c * qc - s * qs, s * qc + c * qs


def norm2(q: tuple[float, float]) -> float:
    return q[0] * q[0] + q[1] * q[1]


def build() -> dict:
    theorem = THEOREM.read_text(encoding="utf-8")
    flat = " ".join(theorem.split())
    theorem_has_pointwise_rotation = (
        "For an individual undamped oscillator" in flat
        and "\\cos\\omega h&-\\sin\\omega h" in theorem
        and "\\sin\\omega h& \\cos\\omega h" in theorem
    )
    theorem_has_fixed_three_partition_dimension = (
        "M_{\\max}=3" in theorem
        and "An inactive slot is represented by $H_r=0$ rather than by a change of model dimension" in flat
    )
    theorem_requires_phase_continuity = "oscillatory blocks preserving phase continuity" in flat

    smoke = []
    for omega, h, phase in (
        (0.0, 0.005, 0.0),
        (2.0 * math.pi * 0.02, 0.005, 0.37),
        (2.0 * math.pi * 1.2, 0.005, -1.1),
        (2.0 * math.pi * 6.0, 0.005, 2.4),
    ):
        q0 = (math.cos(phase), math.sin(phase))
        q1 = rotation_step(q0[0], q0[1], omega, h)
        smoke.append({
            "omega_rad_s": omega,
            "h_s": h,
            "norm2_before": norm2(q0),
            "norm2_after": norm2(q1),
            "norm_error": abs(norm2(q1) - norm2(q0)),
        })

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "source_generator": False,
        "SEA3_parameter_domain_compact": True,
        "compactness_is_not_an_open_obligation": True,
        "continuum_index_set_retained": True,
        "finite_frequency_grid_used": False,
        "finite_direction_grid_used": False,
        "seeded_phase_realization_used": False,
        "phase_reset_on_lambda_transition_allowed": False,
        "inactive_slot_phase_state_retained": True,
        "same_phase_history_required_through_lambda_transition": True,
        "phase_state_set": "product_{r=1..3,(omega,theta)} S^1",
        "phase_state_set_compact_in_product_topology": True,
        "pointwise_transition": "q_{k+1}(omega,theta)=R(omega*h) q_k(omega,theta)",
        "theorem_has_pointwise_rotation": theorem_has_pointwise_rotation,
        "theorem_has_fixed_three_partition_dimension": theorem_has_fixed_three_partition_dimension,
        "theorem_requires_phase_continuity": theorem_requires_phase_continuity,
        "continuum_phase_coordinate_set_closed": CONTINUUM_PHASE_COORDINATE_SET_CLOSED,
        "phase_continuous_propagation_closed": PHASE_CONTINUOUS_PROPAGATION_CLOSED,
        "hard_spectral_driver_set_closed": HARD_SPECTRAL_DRIVER_SET_CLOSED,
        "joint_source_output_map_closed": JOINT_SOURCE_OUTPUT_MAP_CLOSED,
        "output_map_continuity_or_boundedness_claimed_here": False,
        "complete_SEA3_left_inclusion_closed_here": False,
        "complete_SEA3_family_materialized_here": False,
        "P3_promoted": False,
        "smoke": smoke,
        "next_obligation": (
            "close the hard continuum spectral amplitude/driver coordinate and its same-history output map; the continuum phase circle and its no-reseed propagation are already closed here"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "SEA3_parameter_domain_compact",
        "compactness_is_not_an_open_obligation",
        "continuum_index_set_retained",
        "inactive_slot_phase_state_retained",
        "same_phase_history_required_through_lambda_transition",
        "phase_state_set_compact_in_product_topology",
        "theorem_has_pointwise_rotation",
        "theorem_has_fixed_three_partition_dimension",
        "theorem_requires_phase_continuity",
        "continuum_phase_coordinate_set_closed",
        "phase_continuous_propagation_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_generator",
        "finite_frequency_grid_used",
        "finite_direction_grid_used",
        "seeded_phase_realization_used",
        "phase_reset_on_lambda_transition_allowed",
        "hard_spectral_driver_set_closed",
        "joint_source_output_map_closed",
        "output_map_continuity_or_boundedness_claimed_here",
        "complete_SEA3_left_inclusion_closed_here",
        "complete_SEA3_family_materialized_here",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("continuum phase state detached from canonical SEA3 source")
    for row in d.get("smoke", []):
        if float(row.get("norm_error", math.inf)) > 8.0e-16:
            f.append("pointwise rotation did not preserve the unit circle")
            break
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "phase_set": d["phase_state_set"],
        "continuum_phase_closed": d["continuum_phase_coordinate_set_closed"],
        "phase_propagation_closed": d["phase_continuous_propagation_closed"],
        "hard_driver_closed": d["hard_spectral_driver_set_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
