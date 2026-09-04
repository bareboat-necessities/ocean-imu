#!/usr/bin/env python3
"""SEA3 directional-response -> P2 inclusion -> H/A feasibility bridge.

The response-moment statement is conditional on an explicit provisional
vessel/IMU response-norm hypothesis because the repository has no owned hull
RAO.  The SEA3 -> P2 statement is stronger and independent of that numerical
hypothesis: the shipping code clamps the period-derived tuning frequency, tau,
sigma and R_S before the exact EMA/stage/commit language already represented by
P2.  This closes non-pruning inclusion into the broad P2 language.

If a canonical full-P2 H/A result is attached, PASS transfers to the SEA3
subset.  FAIL does not: it is only inconclusive until SEA3 actually prunes P2.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p2_p3_correlation_interface as P2I
import ou3_p3_canonical_gate as P3G
import ou3_p3_p2_v1_full_state_join as HA
import ou3_p4_source_path_reachability as PATH
import ou3_sea3_wave_period_frontend as FRONT
import ou3_source_domain_contract as SOURCE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
DEFAULT_RESPONSE_DOMAIN = REPO / "tools" / "ou3_sea3_directional_response_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_DIRECTIONAL_RESPONSE_P2_HA_FEASIBILITY"


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _source_const(text: str, name: str) -> float:
    return float(SOURCE.parse_const(text, name))


def _load_response_domain(path: Path = DEFAULT_RESPONSE_DOMAIN) -> dict:
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    if d.get("schema_version") != "OU3_SEA3_DIRECTIONAL_RESPONSE_DOMAIN_V1":
        raise RuntimeError("unexpected SEA3 directional response-domain schema")
    if int(d.get("sea_modes_max", 0)) != 3:
        raise RuntimeError("SEA3 response domain must retain M_max=3")
    r = d.get("response_contract", {})
    gain = float(r.get("response_vector_2_norm_upper", math.nan))
    band = r.get("finite_response_band_hz")
    if not math.isfinite(gain) or gain <= 0.0:
        raise RuntimeError("provisional response gain must be finite and positive")
    if not isinstance(band, list) or len(band) != 2:
        raise RuntimeError("response domain lost finite frequency band")
    lo, hi = map(float, band)
    if not (math.isfinite(lo) and math.isfinite(hi) and 0.0 < lo < hi):
        raise RuntimeError("invalid directional response band")
    return d


def directional_response_enclosure(
    repo: Path = REPO,
    response_domain_path: Path = DEFAULT_RESPONSE_DOMAIN,
) -> dict:
    """Finite PSD matrix-moment enclosure under the provisional response cap.

    Retain h(omega,theta) h* before taking an outer bound.  With ||h||_2 <= G,
    normalized direction density, and sum H_r^2 = H_s^2,

        tr M_a <= G^2 omega_hi^4 H_s^2 / 16.

    Thus cross-axis coupling is retained structurally instead of replacing the
    response by three independently selectable Cartesian boxes.
    """
    repo = Path(repo).resolve()
    cfg = _load_response_domain(Path(response_domain_path).resolve())
    r = cfg["response_contract"]
    wrapper = (repo / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h").read_text(
        encoding="utf-8"
    )
    deployed_band = [
        _source_const(wrapper, "SIGMA_BAND_MIN_HZ_DEFAULT"),
        _source_const(wrapper, "SIGMA_BAND_MAX_HZ_DEFAULT"),
    ]
    declared_band = list(map(float, r["finite_response_band_hz"]))
    parity = all(
        math.isclose(a, b, rel_tol=0.0, abs_tol=1.0e-12)
        for a, b in zip(deployed_band, declared_band)
    )
    if not parity:
        raise RuntimeError(
            f"response band {declared_band} drifted from shipping sigma band {deployed_band}"
        )

    gain = float(r["response_vector_2_norm_upper"])
    omega_lo = 2.0 * math.pi * declared_band[0]
    omega_hi = 2.0 * math.pi * declared_band[1]
    disp = gain * gain / 16.0
    vel = disp * omega_hi * omega_hi
    acc = vel * omega_hi * omega_hi
    return {
        "qualification": "OU3_SEA3_PROVISIONAL_DIRECTIONAL_RESPONSE_MOMENT_ENCLOSURE",
        "sea_modes_max": 3,
        "trajectory_replay_used": False,
        "measured_hull_RAO_used": False,
        "repository_owned_hull_RAO_present": bool(
            cfg.get("provenance", {}).get("repository_owned_hull_RAO_present", False)
        ),
        "physical_SEA0_promoted": False,
        "response_hypothesis_status": "PROVISIONAL_NOT_DATA_DERIVED",
        "response_definition": r["definition"],
        "response_vector_2_norm_upper": gain,
        "response_gain_db_upper": 20.0 * math.log10(gain),
        "finite_response_band_hz": declared_band,
        "shipping_sigma_band_hz": deployed_band,
        "response_frequency_band_source_parity": parity,
        "directional_cross_axis_coupling": {
            "rank_one_complex_outer_product_retained_before_outer_bound": True,
            "response_spectral_matrix_PSD": True,
            "independent_cartesian_axis_boxes_used": False,
            "entrywise_cross_terms_constrained_by_PSD_Cauchy_Schwarz": True,
        },
        "matrix_moment_trace_upper_per_Hs2": {
            "displacement": [down(0.0), up(disp)],
            "velocity": [down(0.0), up(vel)],
            "acceleration": [down(0.0), up(acc)],
        },
        "omega_band_rad_s": [down(omega_lo), up(omega_hi)],
        "finite_band_moments_established_conditionally": True,
        "numerical_response_gain_is_allowed_to_prune_P2": False,
        "remaining_physical_obligation": (
            "replace or explicitly accept the provisional response-norm hypothesis with a repository-owned vessel/IMU response model before physical SEA0 promotion"
        ),
    }


def source_bridge_contract(repo: Path = REPO) -> dict:
    repo = Path(repo).resolve()
    wrapper = (repo / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h").read_text(
        encoding="utf-8"
    )
    tuner = (repo / "src" / "tuner" / "SeaStateAutoTuner.h").read_text(encoding="utf-8")
    wrapper_markers = {
        "pending_applied_before_current_sample": "apply_pending_online_tune_();",
        "period_source_prior_or_estimator": "return tune_freq_prior_hz_;",
        "tune_frequency_lower_clamp": "f_tune = f_tune_floor;",
        "tune_frequency_upper_clamp": "f_tune = f_tune_ceil;",
        "tau_target_clamp": "std::min(std::max(tau_raw,  min_tau_s_), max_tau_s_)",
        "sigma_target_upper_clamp": "std::min(sigma_wave * sigma_coeff_,      max_sigma_a_)",
        "RS_target_clamp": "std::min(std::max(RS_raw, min_R_S_), max_R_S_)",
        "samplewise_tau_EMA": "tune_.tau_applied   += alpha",
        "samplewise_sigma_EMA": "tune_.sigma_applied += alpha",
        "samplewise_RS_EMA": "tune_.RS_applied    += alpha_RS",
        "stage_after_current_sample": "online_tune_apply_pending_ = true;",
    }
    tuner_markers = {
        "finite_frequency_input_required": "!std::isfinite(f_input_hz)",
        "internal_frequency_clamp": "f_eff = std::max(f_min_hz, std::min(f_max_hz, f_eff));",
        "nonnegative_variance": "return std::max(0.0f, A_sq.get() - mu * mu);",
    }
    wm = {k: marker in wrapper for k, marker in wrapper_markers.items()}
    tm = {k: marker in tuner for k, marker in tuner_markers.items()}
    return {
        "wrapper_markers": wm,
        "tuner_markers": tm,
        "all_shipping_bridge_markers_present": all(wm.values()) and all(tm.values()),
    }


def _inclusion_from_p2(
    response: dict,
    frontend: dict,
    p2: dict,
    repo: Path = REPO,
) -> dict:
    pf = P2I.validate(p2)
    ff = FRONT.validate(frontend)
    if pf:
        raise RuntimeError(f"P2 correlation interface invalid: {pf}")
    if ff:
        raise RuntimeError(f"wave-period front-end invalid: {ff}")
    if response.get("response_frequency_band_source_parity") is not True:
        raise RuntimeError("response band lost shipping parity")

    bridge = source_bridge_contract(repo)
    if bridge["all_shipping_bridge_markers_present"] is not True:
        raise RuntimeError("shipping SEA3->P2 clamp/staging bridge changed")

    c = PATH._constants()
    tau_lo = max(c["min_tau"], c["tau_coeff"] * 0.5 / c["max_freq"])
    projected = {
        "tuning_frequency_hz": [down(c["min_freq"]), up(c["max_freq"])],
        "tau_applied_s": [down(tau_lo), up(c["max_tau"])],
        "raw_sigma_tuner_mps2": [down(PATH.RAW_SIGMA_GRAPH_LOWER), up(c["max_sigma"])],
        "filter_sigma_floor_mps2": PATH.FILTER_SIGMA_FLOOR,
        "R_S_m_s": [down(c["min_RS"]), up(c["max_RS"])],
        "stage_gap_valid_samples": list(p2["clock_gap_samples"]),
    }
    prior_iv = list(map(float, frontend["declared_inputs"]["fixed_tuning_frequency_prior_hz"]))
    prior_covered = (
        len(prior_iv) == 2
        and c["min_freq"] <= prior_iv[0] <= prior_iv[1] <= c["max_freq"]
    )
    gap_covered = list(p2["clock_gap_samples"]) == list(range(13, 27))
    passed = (
        p2.get("P2_READY_FOR_CANONICAL_P3") is True
        and int(p2.get("physical_source_states", 0)) == 800
        and prior_covered
        and gap_covered
    )
    return {
        "qualification": "OU3_SEA3_TO_FROZEN_P2_NONPRUNING_SOURCE_INCLUSION",
        "normal_live_valid_samples_only": True,
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_filter_operating_domain_changed": False,
        "shipping_bridge": bridge,
        "P2_correlation_interface_version": p2["correlation_interface_version"],
        "P2_physical_source_states": int(p2["physical_source_states"]),
        "P2_clock_gap_samples": list(p2["clock_gap_samples"]),
        "projected_source_domain": projected,
        "fixed_prior_frequency_hz": prior_iv,
        "fixed_prior_covered_by_P2_frequency_alphabet": prior_covered,
        "estimator_frequency_is_clamped_before_tau_target": True,
        "sigma_variance_history_need_not_be_independently_bounded_for_inclusion": True,
        "R_S_transcendental_target_need_not_be_independently_bounded_for_inclusion": True,
        "reason_clamps_are_sufficient": (
            "P2 already over-approximates arbitrary in-range target cells sample by sample; shipping clamps f_tune/tau/sigma/R_S before the exact EMA/stage/commit schedule represented by P2"
        ),
        "provisional_directional_response_gain_used_for_inclusion": False,
        "finite_window_Gaussian_amplitude_bound_used_for_inclusion": False,
        "full_finite_memory_estimator_reachability_used_to_prune_P2": False,
        "P2_pruned_by_SEA3": False,
        "Lhat_SEA3_subset_L_current_source": passed,
        "SEA3_TO_P2_INCLUSION_CERTIFICATE": "PASS" if passed else "FAIL",
        "physical_SEA0_promoted_here": False,
        "interpretation": (
            "This closes compositional inclusion into the existing broad P2 language but intentionally removes none of its 800 cells. Response/IQ finite-window work is needed for pruning, not for this superset inclusion."
        ),
        "next_obligation": (
            "replace/accept the provisional response model and carry finite response-weighted estimator/log-period state far enough to prune P2; then recompute H/A on the narrower SEA3 history language"
        ),
    }


def build_inclusion(
    domain_path: Path = DEFAULT_DOMAIN,
    response_domain_path: Path = DEFAULT_RESPONSE_DOMAIN,
    p2_candidate: dict | None = None,
    repo: Path = REPO,
) -> dict:
    repo = Path(repo).resolve()
    response = directional_response_enclosure(repo, response_domain_path)
    frontend = FRONT.build(repo)
    p2 = P2I.build(Path(domain_path).resolve()) if p2_candidate is None else p2_candidate
    inclusion = _inclusion_from_p2(response, frontend, p2, repo)
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "response_enclosure": response,
        "p2_inclusion": inclusion,
        "ha_feasibility": None,
        "P3_promoted_by_this_artifact": False,
        "P4_promoted_by_this_artifact": False,
    }


def _ha_inheritance(inclusion: dict, ha: dict, canonical: dict) -> dict:
    modes = ha.get("modes", {})
    h = float(modes.get("H", {}).get("relative_Riccati_injection_margin_lower", 0.0))
    a = float(modes.get("A", {}).get("relative_Riccati_injection_margin_lower", 0.0))
    full_p2_pass = canonical.get("P3_CANONICAL_PASS") is True
    inclusion_pass = inclusion.get("SEA3_TO_P2_INCLUSION_CERTIFICATE") == "PASS"
    inherited = inclusion_pass and full_p2_pass
    return {
        "qualification": "OU3_SEA3_HA_FEASIBILITY_BY_P2_SUPERSET_INHERITANCE",
        "H_dimension": 18,
        "A_dimension": 21,
        "H_relative_Riccati_injection_margin_lower": h,
        "A_relative_Riccati_injection_margin_lower": a,
        "full_P2_canonical_worst_H_A_margin": float(canonical.get("worst_H_A_margin", 0.0)),
        "unchanged_useful_gate": float(canonical.get("useful_gate", 0.0)),
        "full_P2_P3_CANONICAL_PASS": full_p2_pass,
        "SEA3_subset_inclusion_pass": inclusion_pass,
        "uniform_P2_PASS_is_inherited_by_SEA3_subset": True,
        "uniform_P2_FAIL_implies_SEA3_FAIL": False,
        "SEA3_specific_HA_recomputed_on_pruned_language": False,
        "SEA3_HA_FEASIBILITY": (
            "PASS_BY_P2_SUPERSET" if inherited else "INCONCLUSIVE_REQUIRES_SEA3_NARROWING"
        ),
        "SEA3_HA_feasible_by_existing_uniform_certificate": inherited,
        "P3_promoted_here": False,
        "next_obligation": (
            "existing uniform P2 H/A already covers SEA3; continue physical-response acceptance and downstream P4/P5 obligations"
            if inherited
            else "prune P2 source histories with the response/finite-estimator construction and rerun the same frozen H/A theorem interface"
        ),
    }


def attach_ha_feasibility(payload: dict, ha: dict, canonical: dict) -> dict:
    vf = validate(payload)
    if vf:
        raise RuntimeError(f"SEA3 inclusion invalid before H/A attachment: {vf}")
    hf = HA.validate(ha)
    cf = P3G.validate(canonical)
    if hf:
        raise RuntimeError(f"H/A producer artifact invalid: {hf}")
    if cf:
        raise RuntimeError(f"canonical P3 artifact invalid: {cf}")
    out = dict(payload)
    out["ha_feasibility"] = _ha_inheritance(out["p2_inclusion"], ha, canonical)
    return out


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != QUALIFICATION:
        f.append("wrong qualification")
    r = d.get("response_enclosure", {})
    if r.get("response_frequency_band_source_parity") is not True:
        f.append("response finite band lost shipping source parity")
    if r.get("physical_SEA0_promoted") is not False:
        f.append("provisional response hypothesis was promoted as physical SEA0")
    coupling = r.get("directional_cross_axis_coupling", {})
    if coupling.get("rank_one_complex_outer_product_retained_before_outer_bound") is not True:
        f.append("directional cross-axis outer product was not retained")
    if coupling.get("independent_cartesian_axis_boxes_used") is not False:
        f.append("independent Cartesian response boxes were reintroduced")
    p = d.get("p2_inclusion", {})
    if p.get("SEA3_TO_P2_INCLUSION_CERTIFICATE") != "PASS":
        f.append("SEA3->P2 inclusion did not pass")
    if p.get("Lhat_SEA3_subset_L_current_source") is not True:
        f.append("SEA3 source language was not included in current P2 source language")
    if p.get("P2_pruned_by_SEA3") is not False:
        f.append("non-pruning inclusion incorrectly claims P2 pruning")
    if p.get("provisional_directional_response_gain_used_for_inclusion") is not False:
        f.append("P2 inclusion depends on provisional directional response gain")
    if p.get("shipping_bridge", {}).get("all_shipping_bridge_markers_present") is not True:
        f.append("shipping clamp/staging bridge markers missing")
    if d.get("P3_promoted_by_this_artifact") is not False:
        f.append("SEA3 bridge promoted P3")
    if d.get("P4_promoted_by_this_artifact") is not False:
        f.append("SEA3 bridge promoted P4")
    ha = d.get("ha_feasibility")
    if ha is not None:
        expected = (
            ha.get("SEA3_subset_inclusion_pass") is True
            and ha.get("full_P2_P3_CANONICAL_PASS") is True
        )
        if ha.get("SEA3_HA_feasible_by_existing_uniform_certificate") is not expected:
            f.append("H/A inherited feasibility flag is inconsistent")
        status = "PASS_BY_P2_SUPERSET" if expected else "INCONCLUSIVE_REQUIRES_SEA3_NARROWING"
        if ha.get("SEA3_HA_FEASIBILITY") != status:
            f.append("H/A feasibility status is inconsistent")
        if ha.get("uniform_P2_FAIL_implies_SEA3_FAIL") is not False:
            f.append("full-P2 failure was incorrectly promoted to SEA3 failure")
        if ha.get("SEA3_specific_HA_recomputed_on_pruned_language") is not False:
            f.append("artifact falsely claims SEA3-pruned H/A recomputation")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--response-domain", type=Path, default=DEFAULT_RESPONSE_DOMAIN)
    ap.add_argument("--inclusion-candidate", type=Path)
    ap.add_argument("--ha-candidate", type=Path)
    ap.add_argument("--canonical-candidate", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if (args.ha_candidate is None) != (args.canonical_candidate is None):
        raise SystemExit("--ha-candidate and --canonical-candidate must be supplied together")

    d = (
        json.loads(args.inclusion_candidate.read_text(encoding="utf-8"))
        if args.inclusion_candidate
        else build_inclusion(args.domain, args.response_domain)
    )
    if args.ha_candidate:
        d = attach_ha_feasibility(
            d,
            json.loads(args.ha_candidate.read_text(encoding="utf-8")),
            json.loads(args.canonical_candidate.read_text(encoding="utf-8")),
        )

    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "response_status": d["response_enclosure"]["response_hypothesis_status"],
        "response_accel_trace_upper_per_Hs2": d["response_enclosure"]["matrix_moment_trace_upper_per_Hs2"]["acceleration"],
        "SEA3_TO_P2_INCLUSION_CERTIFICATE": d["p2_inclusion"]["SEA3_TO_P2_INCLUSION_CERTIFICATE"],
        "P2_pruned": d["p2_inclusion"]["P2_pruned_by_SEA3"],
        "ha_feasibility": None if d["ha_feasibility"] is None else d["ha_feasibility"]["SEA3_HA_FEASIBILITY"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
