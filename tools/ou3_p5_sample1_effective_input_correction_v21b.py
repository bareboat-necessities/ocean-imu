#!/usr/bin/env python3
"""V21B: bind the V21 first-q8 witness audit to the live V14D chart semantics.

V21 selected the correct V18B witness indices, but it evaluated its copied
current-chart helper outside V14D.build.  The live V18B route evaluates that
helper while V14D has installed its source-faithful radial-sinc normalized
shipping quaternion.  Outside that context V21 silently used V14's older
quaternion primitive and reproduced the pre-V17 current radius q=0.659377...
instead of the authoritative V18B witness q=0.641523....

V21B changes no proof algebra.  It runs the complete V21 focused diagnostic
while the same V14D radial-sinc quaternion primitive is installed, restores the
module afterward, and requires the deterministic current radius to reproduce
the already measured V18B first-witness value.  The V21 correction intersection
remains a diagnostic parent: no q8, sample-1, word, or P5 promotion occurs here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Callable, TypeVar

import ou3_p5_sample1_effective_input_correction_v21 as V21
import ou3_p5_sample1_signed_cayley_q8_v14d as V14D

DEFAULT_DOMAIN = V21.DEFAULT_DOMAIN
SCHEMA = 2101
Q_TARGET = V21.Q_TARGET
V18B_FIRST_WITNESS_CURRENT_Q = 0.6415230535178351
_T = TypeVar("_T")


def _with_v14d_quaternion(fn: Callable[[], _T]) -> _T:
    """Execute fn under the exact V14D quaternion primitive and restore it."""
    original = V21.V14._normalized_shipping_quaternion
    V21.V14._normalized_shipping_quaternion = V14D.radial_sinc_normalized_shipping_quaternion
    try:
        return fn()
    finally:
        V21.V14._normalized_shipping_quaternion = original


def _matches_reference(q: float) -> bool:
    q = float(q)
    ref = V18B_FIRST_WITNESS_CURRENT_Q
    if not math.isfinite(q):
        return False
    tol = 64.0 * math.ulp(max(1.0, abs(ref)))
    return abs(q - ref) <= tol


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    path = Path(domain_path).resolve()

    def run_parent():
        return V21.build(
            path,
            source_pieces=source_pieces,
            source_cell_index=source_cell_index,
            p_pieces=p_pieces,
            tangent_pieces=tangent_pieces,
            axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces,
        )

    core = _with_v14d_quaternion(run_parent)
    inherited = V21.validate(core)
    q = float(core.get("sample1_current_cayley_norm_upper", math.inf))
    matches = _matches_reference(q)
    failures = list(inherited)
    if not matches:
        failures.append(
            "V21B current chart does not reproduce authoritative V18B first-witness q")

    parent_status = core.get("P5_SAMPLE1_EFFECTIVE_INPUT_CORRECTION_WITNESS_V21")
    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_V14D_BOUND_EFFECTIVE_INPUT_WITNESS_V21B",
        "V21_effective_input_parent_retained": True,
        "V14D_radial_sinc_quaternion_installed_for_V21_audit": True,
        "V14D_quaternion_restored_after_audit": True,
        "authoritative_V18B_first_witness_current_q_reference": V18B_FIRST_WITNESS_CURRENT_Q,
        "current_q_matches_authoritative_V18B_reference": matches,
        "P5_SAMPLE1_V14D_BOUND_EFFECTIVE_INPUT_WITNESS_V21B": (
            "PASS" if parent_status == "PASS" and not failures else "NOT_ESTABLISHED"
        ),
        "next_obligation": core.get(
            "next_obligation",
            "DERIVE_SOURCE_CORRELATED_SAMPLE1_AW_ERROR_COMPONENTS_AT_FIRST_Q8_WITNESS"),
        "failures": list(dict.fromkeys(failures)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_V14D_BOUND_EFFECTIVE_INPUT_WITNESS_V21B":
        f.append("qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "V21_effective_input_parent_retained",
        "V14D_radial_sinc_quaternion_installed_for_V21_audit",
        "V14D_quaternion_restored_after_audit",
        "current_q_matches_authoritative_V18B_reference",
        "V12D_PSD_S_perturbation_retained",
        "V10_one_plus_two_gain_retained",
        "exact_accelerometer_effective_input_lemma_used",
        "current_cayley_and_sample1_correction_jointly_mapped",
        "V13E_signed_subcell_intersected_not_replaced",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "q8_composed_here", "q8_word_promoted_here", "whole_word_promoted_here",
        "N_H_words_set_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    q = d.get("sample1_current_cayley_norm_upper")
    if not isinstance(q, (int, float)) or not _matches_reference(float(q)):
        f.append("authoritative first-witness current q mismatch")
    st = d.get("P5_SAMPLE1_V14D_BOUND_EFFECTIVE_INPUT_WITNESS_V21B")
    if st not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V21B status")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--residual-x-pieces", type=int, default=6)
    ap.add_argument("--parallel-pieces", type=int, default=6)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(
        x.domain,
        source_pieces=x.source_pieces,
        source_cell_index=x.source_cell_index,
        p_pieces=x.p_pieces,
        tangent_pieces=x.tangent_pieces,
        axial_pieces=x.axial_pieces,
        residual_x_pieces=x.residual_x_pieces,
        parallel_pieces=x.parallel_pieces,
    )
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_V14D_BOUND_EFFECTIVE_INPUT_WITNESS_V21B"],
        "q_current": d.get("sample1_current_cayley_norm_upper"),
        "q_reference": d["authoritative_V18B_first_witness_current_q_reference"],
        "baseline_radial_upper": d.get("baseline_correction_radial_upper_rad"),
        "refined_radial_upper": d.get("refined_correction_radial_upper_rad"),
        "baseline_geodesic_q": d.get("baseline_geodesic_q_upper"),
        "refined_geodesic_q": d.get("refined_geodesic_q_upper"),
        "refined_product_W": d.get("refined_product_abs_W_lower"),
        "refined_product_q": d.get("refined_product_q_upper"),
        "closed_q8": d.get("first_witness_closed_inside_q8"),
        "post_prediction_aw": d.get("post_prediction_physical_aw_error_norm_upper_mps2"),
        "bias": d.get("accelerometer_bias_error_norm_upper_mps2"),
        "attitude_defect": d.get("effective_attitude_defect_norm_upper_mps2"),
        "latent_cross": d.get("effective_latent_cross_norm_upper_mps2"),
        "nuisance": d.get("nominal_effective_residual_nuisance_norm_upper_mps2"),
        "V12D_correction_perturbation": d.get("V12D_correction_perturbation_norm_upper_rad"),
        "source_incompatible": d.get("source_subcell_incompatible"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
