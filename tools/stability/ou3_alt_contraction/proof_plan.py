"""Executable anti-dead-end phase guards for the ALT proof workflow.

These guards are process controls, not mathematical evidence.  They prevent a
metric/high-precision/storage task from being accidentally treated as useful
proof work before the physical source-uniform word exists.
"""
from __future__ import annotations

DIAGNOSTIC_ONLY = {"replay", "finite_seed", "unreachable_perturbation", "captured_trace"}


def require_theorem_task(*, obligation: str, evidence_kind: str, complete_physical_word: bool,
                         requested_phase: str) -> None:
    if not isinstance(obligation,str) or not obligation.strip():
        raise RuntimeError("ALT task has no named theorem obligation")
    if evidence_kind in DIAGNOSTIC_ONLY:
        raise RuntimeError("diagnostic-only evidence cannot be the main ALT proof task")
    if requested_phase in {"storage_search","high_precision_master","interval_refinement"} and not complete_physical_word:
        raise RuntimeError("complete source-uniform physical word is required before storage/high-precision/refinement")


def storage_search_allowed(status: dict) -> bool:
    required=("same_history_complete_BRMM_word","physical_prediction_forcing_attached",
              "physical_S_residual_attached","all_bias_families_attached",
              "all_literal_branches_attached","H18_A21_edge_attached")
    return all(status.get(k) is True for k in required)


def assert_storage_search_allowed(status: dict) -> None:
    if not storage_search_allowed(status):
        missing=[k for k in ("same_history_complete_BRMM_word","physical_prediction_forcing_attached",
              "physical_S_residual_attached","all_bias_families_attached",
              "all_literal_branches_attached","H18_A21_edge_attached") if status.get(k) is not True]
        raise RuntimeError("ALT storage search blocked; incomplete theorem word: "+", ".join(missing))


def assert_finite_storage_master(status: dict) -> None:
    """Reject a Jacobian/metadata assembly as input to a finite-state rho search.

    This is an additional process guard, not mathematical evidence. A finite
    identity and its source coverage must be proved by the supplying builder.
    In particular the older six assembly flags cannot imply these obligations.
    """
    required = ('finite_error_identity_for_every_event',
                'physical_reference_forcing_retained',
                'all_coefficient_product_graphs_retained',
                'all_configured_branches_bound_to_finite_graph')
    missing = [k for k in required if status.get(k) is not True]
    if status.get('map_representation') != 'finite_physical_descriptor':
        missing.insert(0, 'finite_physical_descriptor (not a Jacobian cocycle)')
    if missing:
        raise RuntimeError('ALT finite-state storage blocked: '+', '.join(missing))
