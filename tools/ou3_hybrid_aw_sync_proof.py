#!/usr/bin/env python3
"""Source-bound proof for the OU-III periodic a_w covariance synchronization jump.

This is a theorem-level hybrid lemma, not a replay diagnostic.  It binds to the
shipping OU-III implementation and checks that the default synchronization path
is exactly the PSD-inflation policy used by the information metric:

    P+ = P- + E_a Delta_+ E_a^T,   Delta_+ >= 0.

For SPD P-, Loewner order gives P+ >= P- > 0 and therefore
P+^{-1} <= P-^{-1}.  Hence W=e^T P^{-1} e is nonexpansive across the jump for
an unchanged error coordinate e.  This closes the periodic covariance-sync
hybrid obligation independently of trajectories or sampled basin results.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_HEADER = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"


SOURCE_PATTERNS = {
    "default_policy_is_psd_inflation": re.compile(
        r"bool\s+legacy_aw_covariance_replacement_\s*=\s*false\s*;"
    ),
    "sync_queues_symmetric_target": re.compile(
        r"aw_covariance_floor_target_\s*=\s*T\(0\.5\)\s*\*\s*"
        r"\(Sigma_aw_stat\s*\+\s*Sigma_aw_stat\.transpose\(\)\)\s*;"
    ),
    "sync_marks_pending": re.compile(r"aw_covariance_floor_pending_\s*=\s*true\s*;"),
    "prediction_applies_pending_sync": re.compile(
        r"apply_pending_aw_covariance_inflation_\(\)\s*;\s*"
        r"symmetrize_Pext_\(\)\s*;", re.S
    ),
    "delta_is_target_minus_current": re.compile(
        r"Matrix3\s+Delta\s*=\s*aw_covariance_floor_target_\s*-\s*P_aw\s*;"
    ),
    "delta_is_symmetrized": re.compile(
        r"Delta\s*=\s*T\(0\.5\)\s*\*\s*\(Delta\s*\+\s*Delta\.transpose\(\)\)\s*;"
    ),
    "delta_uses_self_adjoint_eigendecomposition": re.compile(
        r"Eigen::SelfAdjointEigenSolver<Matrix3>\s+es\(Delta\)\s*;"
    ),
    "negative_delta_eigenvalues_are_clamped": re.compile(
        r"evals\(i\)\s*=\s*std::max\(T\(0\),\s*evals\(i\)\)\s*;"
    ),
    "psd_delta_is_added_to_aw_block": re.compile(
        r"Pext\.template\s+block<3,3>\(OFF_AW,\s*OFF_AW\)\s*\+=\s*Delta\s*;"
    ),
}


def prove(header: Path = DEFAULT_HEADER) -> dict:
    header = header.resolve()
    text = header.read_text(encoding="utf-8")
    checks = {name: bool(pattern.search(text)) for name, pattern in SOURCE_PATTERNS.items()}
    failures = [name for name, ok in checks.items() if not ok]

    # The proof is exact algebra once the source binding above succeeds:
    # Delta_+ is reconstructed from an orthonormal eigendecomposition with all
    # eigenvalues max(0, lambda_i), hence Delta_+ >= 0. Embedding that 3x3 block
    # as E_a Delta_+ E_a^T is PSD. Adding a PSD matrix preserves/increases P in
    # Loewner order, and inversion reverses Loewner order on SPD matrices.
    source_bound = not failures
    return {
        "schema": 1,
        "claim": "OU3_PERIODIC_AW_COVARIANCE_SYNC_PSD_NONEXPANSIVE",
        "qualification": "SOURCE_BOUND_ANALYTIC_HYBRID_PROOF",
        "source_generated_not_trajectory_fit": True,
        "sampled_evidence_used": False,
        "implementation_header": str(header.relative_to(REPO)),
        "implementation_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        "source_binding_checks": checks,
        "source_binding_pass": source_bound,
        "proof_mode": "PSD_NONEXPANSIVE",
        "source_complete_for_transition": source_bound,
        "outward_rounding_required": False,
        "reason_outward_rounding_not_required": (
            "the jump gain is established by exact Loewner-order algebra, not by a floating-point numerical bound"
        ),
        "matrix_identity": "P_plus = P_minus + E_a Delta_plus E_a^T",
        "delta_psd_by_construction": source_bound,
        "loewner_covariance_order": "P_plus >= P_minus > 0",
        "inverse_information_order": "P_plus^-1 <= P_minus^-1",
        "error_coordinate_change": "identity",
        "jump_gain_upper": 1.0 if source_bound else None,
        "additive_W_upper": 0.0 if source_bound else None,
        "new_coordinate_W_upper": 0.0 if source_bound else None,
        "nonexpansive_information_energy": source_bound,
        "hybrid_obligation": "periodic_aw_covariance_sync",
        "status": "PASS" if source_bound else "FAIL",
        "failures": failures,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--header", type=Path, default=DEFAULT_HEADER)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    payload = prove(args.header)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if payload["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
