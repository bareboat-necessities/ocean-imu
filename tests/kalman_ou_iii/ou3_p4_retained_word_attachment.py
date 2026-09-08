"""Audit attachment of a captured point word; never certify a source by its ID.

The frame identities below apply to finite states and arbitrary finite maps.
They do not identify the captured noisy nominal path with a homogeneous
physical path. That identification requires source/nominal data absent from
OU3PHY1. No replacement Riccati recursion or new witness search is performed.
"""
from __future__ import annotations

import numpy as np

import ou3_p4_physical_finite_map_feasibility as BASE
import ou3_p4_physical_finite_map_feasibility_fast as FAST
import ou3_p4_complete_brmm_signed_information_ledger as SIGNED


def reset_matrix(dtheta, n):
    G = np.eye(n)
    G[:3, :3] += 0.5 * FAST._skew(dtheta)
    return G


def transport_cell(event, T, n):
    """Remove the captured reset by transporting ALL coordinates of one cell.

    z=T e, Pz=T P T^T. At a reset T_next=T G^-1. The physical
    nonlinear map must separately become T_next f(T^-1 z).
    """
    measurement = event["type"] not in (BASE.EV_PRED, BASE.EV_FLOOR)
    G = reset_matrix(event["dtheta"], n) if measurement else np.eye(n)
    Tnext = np.linalg.solve(G.T, T.T).T
    Pin = event["Pbefore_shipping"][:n, :n]
    Pout = event["Pafter_shipping"][:n, :n]
    result = {
        "Tnext": Tnext,
        "Pbefore": T @ Pin @ T.T,
        "Pafter": Tnext @ Pout @ Tnext.T,
    }
    if measurement:
        H = event["H"][:, :n]
        result["H"] = np.linalg.solve(T.T, H.T).T
        _, K, _ = BASE._joseph(Pin, H, event["R"])
        result["K_joseph"] = T @ K
        result["R"] = event["R"].copy()
    elif event["type"] == BASE.EV_PRED:
        result["F"] = np.linalg.solve(T.T, (Tnext @ event["linear_shipping"][:n, :n]).T).T
        result["Q"] = Tnext @ event["Q"][:n, :n] @ Tnext.T
    else:
        # Transport the actual captured PSD increment, never rerun the floor
        # policy on a covariance in different coordinates.
        result["floor_increment"] = T @ (Pout - Pin) @ T.T
    return result


def finite_map_pullback(state, Tbefore, Tafter, finite_map):
    return Tafter @ finite_map(np.linalg.solve(Tbefore, state))


def finite_reset_energy(before, after, K, y, Pj):
    """Attach the exact physical output to its own correction reset/storage.

    This is one event's identity, not a new whole-word endpoint ratio: later
    P/K cells cannot be reused after changing this event's covariance reset.
    """
    d = K @ y
    G = reset_matrix(d[:3], len(before))
    Gi = np.linalg.inv(G)
    q = before - d
    defect = after - G @ q  # includes A21 projection when it is active
    Jj = BASE._spd_inverse(Pj)
    terms = SIGNED.reset_signed_energy_terms(Jj.tolist(), Gi.tolist(), q.tolist(), defect.tolist())
    pulled = Gi @ after
    attached = float(pulled @ Jj @ pulled)
    direct = attached - float(q @ Jj @ q)
    return {
        "reset_defect_norm": float(np.linalg.norm(defect)),
        "reset_cross_term": float(terms["reset_cross_term"]),
        "reset_defect_energy": float(terms["reset_defect_energy"]),
        "signed_identity_residual": direct - float(terms["signed_delta"]),
        "relative_signed_identity_residual": abs(direct - float(terms["signed_delta"])) / max(1.0, attached),
        "attached_minus_unreset_storage": attached - float(after @ Jj @ after),
    }


def audit_captured_transport(raw, linear):
    n = raw["mode_dim"]
    T = np.eye(n)
    maxima = {key: 0.0 for key in ("F", "Q", "H", "Pafter")}
    first_changed_H = None
    reset_count = 0
    for index, (event, lp) in enumerate(zip(raw["events"], linear["path"], strict=True)):
        transformed = transport_cell(event, T, n)
        for key, original in (("F", event["linear_shipping"][:n, :n]),
                              ("Q", event["Q"][:n, :n]), ("H", event["H"][:, :n])):
            if key in transformed:
                difference = float(np.linalg.norm(transformed[key] - original))
                maxima[key] = max(maxima[key], difference)
                if key == "H" and difference > 1e-12 and first_changed_H is None:
                    first_changed_H = {"event_index": index, "event": event["name"],
                                       "time": event["time"], "required_H_transport_norm": difference}
        maxima["Pafter"] = max(maxima["Pafter"], float(np.linalg.norm(
            transformed["Pafter"] - lp["Pafter"]) / max(1.0, np.linalg.norm(transformed["Pafter"]))))
        reset_count += int(bool(np.linalg.norm(event["dtheta"]) > 0))
        T = transformed["Tnext"]
    return {
        "captured_nonzero_reset_count": reset_count,
        "max_required_transport_difference": maxima,
        "first_measurement_requiring_changed_H": first_changed_H,
        "captured_frame_endpoint": T.tolist(),
        "dropping_resets_with_unchanged_F_Q_H_is_not_a_gauge_attachment": first_changed_H is not None,
        "finite_coordinate_pullback_identity_available": True,
        "captured_nominal_to_homogeneous_physical_map_identified": False,
        "canonical_finite_word_storage_attached": False,
    }


def source_evidence_inventory(payload):
    """Identify unavailable proof inputs from the decoded capture schema.

    Presence would only make validation possible, never prove membership.
    The actual BRMM membership predicate is X^s_BRMM plus its joint output
    map; parameter bounds or a PSD alone do not define that predicate.
    """
    required = (
        "joint_shaping_root", "hard_driver_history", "continuum_phase_history",
        "partition_parameter_history", "joint_response_witness",
        "frontend_tuner_root", "raw_measurement_history", "nominal_state_history",
        "true_bias_root_and_history", "homogeneous_forcing_decomposition",
    )
    missing = [key for key in required if key not in payload]
    return {
        "capture_schema": "OU3PHY1",
        "decoded_event_fields": sorted(payload["events"][0]),
        "missing_common_realization_inputs": missing,
        "membership_decision": "UNDETERMINED" if missing else "REQUIRES_JOINT_SOURCE_VALIDATION",
        "missing_evidence_is_not_proof_of_nonmembership": True,
        "spectral_label_or_norm_caps_used_as_membership": False,
        "zero_bias_root_substituted_for_captured_truth": False,
        "complete_BRMM_membership_verified": False,
    }
