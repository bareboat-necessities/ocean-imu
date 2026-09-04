#!/usr/bin/env python3
"""SEA3 directional RAO-family -> P2 inclusion -> H/A feasibility bridge.

This producer certifies a continuum of vessel response operators, not a nominal
hull and not a sampled RAO catalogue.  The CoG translational projection h of an
arbitrary finite six-DOF linear RAO is admitted whenever

    ||h(f,theta)||_2 <= G min(1, (f_c/f)^p),  f > 0,

for any G, f_c and p in the declared compact/one-sided parameter family.  Phase,
heading dependence and cross-axis coupling are otherwise arbitrary.  The
moment proof is monotone in G and f_c and worst at p=2, so certifying the single
*envelope corner* (G_max, f_c,max, p_min) certifies every RAO below it.  That
corner is a set bound, not a representative RAO.

The p>=2 roll-off is important: it cancels the omega^4 acceleration weighting
at high frequency, so the acceleration-response moment is finite without
pretending the JONSWAP/PM tail is hard-truncated at 6 Hz.

The right SEA3 -> P2 inclusion remains independent of the response-envelope
numbers because shipping clamps the tuner targets before the exact EMA/staging
language represented by P2.  If canonical full-P2 H/A passes, that stronger
certificate transfers to every member of the SEA3 RAO family.  A full-P2 fail
is only inconclusive until SEA3-specific source pruning is constructed.
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
SCHEMA = 3
QUALIFICATION = "OU3_SEA3_RAO_ENVELOPE_FAMILY_P2_HA_FEASIBILITY"
RESPONSE_SCHEMA = "OU3_SEA3_DIRECTIONAL_RESPONSE_DOMAIN_V3"

# Decimal enclosure of mathematical pi.  Endpoints deliberately straddle pi by
# much more than one binary64 ulp; nextafter widens once more below.
PI_LO = 3.141592653589793
PI_HI = 3.141592653589794


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _source_const(text: str, name: str) -> float:
    return float(SOURCE.parse_const(text, name))


def _load_response_domain(path: Path = DEFAULT_RESPONSE_DOMAIN) -> dict:
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    if d.get("schema_version") != RESPONSE_SCHEMA:
        raise RuntimeError("unexpected SEA3 directional response-domain schema")
    if int(d.get("sea_modes_max", 0)) != 3:
        raise RuntimeError("SEA3 response domain must retain M_max=3")

    r = d.get("response_contract", {})
    gain = r.get("peak_translation_gain_range")
    corner = r.get("rolloff_corner_hz_range")
    if not (isinstance(gain, list) and len(gain) == 2):
        raise RuntimeError("RAO family lost peak-gain range")
    if not (isinstance(corner, list) and len(corner) == 2):
        raise RuntimeError("RAO family lost roll-off-corner range")
    g0, g1 = map(float, gain)
    f0, f1 = map(float, corner)
    pmin = float(r.get("high_frequency_rolloff_power_min", math.nan))
    if not (g0 == 0.0 and math.isfinite(g1) and g1 > 0.0):
        raise RuntimeError("invalid RAO peak-gain range")
    if not (0.0 < f0 <= f1 and math.isfinite(f1)):
        raise RuntimeError("invalid RAO roll-off-corner range")
    if not (math.isfinite(pmin) and pmin >= 2.0):
        raise RuntimeError("RAO acceleration-moment theorem requires p>=2")

    if r.get("worst_member_dominates_entire_parameter_box") is not True:
        raise RuntimeError("RAO family lost monotone worst-envelope contract")
    if r.get("single_nominal_RAO_used") is not False:
        raise RuntimeError("RAO theorem must not select one nominal hull")
    if r.get("finite_RAO_grid_used") is not False:
        raise RuntimeError("RAO theorem must not use a finite RAO proof grid")
    if r.get("phase_quantifier") != "arbitrary complex phase":
        raise RuntimeError("RAO phase quantifier changed")
    for key in (
        "six_dof_parent_RAO_allowed",
        "arbitrary_frequency_dependence_below_envelope",
        "arbitrary_directional_dependence",
        "arbitrary_cross_axis_coupling_subject_to_PSD",
        "unbanded_acceleration_moment_is_finite_from_response_rolloff",
    ):
        if r.get(key) is not True:
            raise RuntimeError(f"RAO family lost {key}")

    worst = r.get("worst_envelope_member", {})
    if not (
        math.isclose(float(worst.get("peak_translation_gain", math.nan)), g1)
        and math.isclose(float(worst.get("rolloff_corner_hz", math.nan)), f1)
        and math.isclose(float(worst.get("high_frequency_rolloff_power", math.nan)), pmin)
    ):
        raise RuntimeError("declared worst RAO envelope is not the parameter-box corner")
    return d


def _member_moment_coefficients(gain: float, corner_hz: float, power: float) -> dict[str, list[float]]:
    """Outward c_q in tr(M_q) <= c_q H_s^2 for one envelope member.

    With m0=H_s^2/16 and ||h||<=G min(1,(f_c/f)^p), p>=2,

      max ||h||^2                         <= G^2,
      max (2*pi*f)^2 ||h||^2             <= (2*pi*f_c)^2 G^2,
      max (2*pi*f)^4 ||h||^2             <= (2*pi*f_c)^4 G^2.

    The last inequality is why p>=2 is sufficient for an unbanded acceleration
    moment even though an unbanded JONSWAP surface-acceleration moment diverges.
    """
    g = float(gain)
    fc = float(corner_hz)
    p = float(power)
    if not (math.isfinite(g) and g >= 0.0):
        raise ValueError("RAO peak gain must be finite and nonnegative")
    if not (math.isfinite(fc) and fc > 0.0):
        raise ValueError("RAO roll-off corner must be finite and positive")
    if not (math.isfinite(p) and p >= 2.0):
        raise ValueError("RAO high-frequency roll-off power must satisfy p>=2")

    g2 = up(g * g)
    omega_c = up(2.0 * PI_HI * fc)
    disp = up(g2 / 16.0)
    vel = up(disp * omega_c * omega_c)
    acc = up(vel * omega_c * omega_c)
    return {
        "displacement": [down(0.0), disp],
        "velocity": [down(0.0), vel],
        "acceleration": [down(0.0), acc],
    }


def evaluate_rao_envelope_member(
    gain: float,
    corner_hz: float,
    power: float,
    response_enclosure: dict | None = None,
) -> dict[str, list[float]]:
    """Evaluate any member of the already-proved continuum parameter box."""
    if response_enclosure is not None:
        box = response_enclosure["rao_envelope_parameter_box"]
        gr = list(map(float, box["peak_translation_gain"]))
        fr = list(map(float, box["rolloff_corner_hz"]))
        pmin = float(box["high_frequency_rolloff_power_min"])
        if not (gr[0] <= float(gain) <= gr[1]):
            raise ValueError("RAO gain lies outside certified parameter box")
        if not (fr[0] <= float(corner_hz) <= fr[1]):
            raise ValueError("RAO corner lies outside certified parameter box")
        if float(power) < pmin:
            raise ValueError("RAO roll-off power lies outside certified parameter box")
    return _member_moment_coefficients(gain, corner_hz, power)


def directional_response_enclosure(
    repo: Path = REPO,
    response_domain_path: Path = DEFAULT_RESPONSE_DOMAIN,
) -> dict:
    """Uniform matrix-moment theorem for every RAO in the declared envelope."""
    repo = Path(repo).resolve()
    cfg = _load_response_domain(Path(response_domain_path).resolve())
    r = cfg["response_contract"]
    wrapper = (repo / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h").read_text(
        encoding="utf-8"
    )

    shipping_band = [
        _source_const(wrapper, "SIGMA_BAND_MIN_HZ_DEFAULT"),
        _source_const(wrapper, "SIGMA_BAND_MAX_HZ_DEFAULT"),
    ]
    declared_band = list(map(float, r["shipping_sigma_band_hz"]))
    parity = all(
        math.isclose(a, b, rel_tol=0.0, abs_tol=1.0e-12)
        for a, b in zip(shipping_band, declared_band)
    )
    if not parity:
        raise RuntimeError(
            f"declared shipping sigma band {declared_band} drifted from source {shipping_band}"
        )

    gr = list(map(float, r["peak_translation_gain_range"]))
    fr = list(map(float, r["rolloff_corner_hz_range"]))
    pmin = float(r["high_frequency_rolloff_power_min"])
    worst = _member_moment_coefficients(gr[1], fr[1], pmin)

    # Diagnostic only: compare with the discarded flat-gain-at-6-Hz outer
    # corner.  The theorem does not use the 6-Hz endpoint for acceleration once
    # the physical response roll-off is retained.
    omega_6 = up(2.0 * PI_HI * declared_band[1])
    flat_disp = up(up(gr[1] * gr[1]) / 16.0)
    flat_acc = up(flat_disp * omega_6**4)
    shaped_acc = float(worst["acceleration"][1])
    reduction = down(flat_acc / shaped_acc)

    return {
        "qualification": "OU3_SEA3_UNIFORM_COMPLEX_DIRECTIONAL_RAO_ENVELOPE",
        "sea_modes_max": 3,
        "trajectory_replay_used": False,
        "measured_hull_RAO_used": False,
        "single_nominal_RAO_used": False,
        "finite_RAO_grid_used": False,
        "physical_SEA0_promoted": False,
        "response_hypothesis_status": "ROBUST_CONTINUUM_ENVELOPE_FAMILY",
        "response_definition": r["definition"],
        "reasonable_RAO_range_definition": cfg["reasonable_RAO_range_definition"],
        "engineering_basis": cfg.get("engineering_basis", {}),
        "rao_envelope_parameter_box": {
            "peak_translation_gain": gr,
            "rolloff_corner_hz": fr,
            "high_frequency_rolloff_power_min": pmin,
            "complex_phase": "arbitrary",
            "heading_dependence": "arbitrary",
            "frequency_dependence_below_envelope": "arbitrary",
            "cross_axis_coupling": "arbitrary subject to PSD",
        },
        "worst_envelope_member": {
            "peak_translation_gain": gr[1],
            "rolloff_corner_hz": fr[1],
            "high_frequency_rolloff_power": pmin,
            "is_set_envelope_not_representative_RAO": True,
        },
        "single_worst_envelope_proves_entire_parameter_box_by_monotonicity": True,
        "six_dof_parent_RAO_allowed": r["six_dof_parent_RAO_allowed"],
        "zero_lever_arm_translation_projection_used": True,
        "rotational_parent_scope": (
            "rotational RAO remains arbitrary subject to the separately declared P1 Normal-Live body-rate/attitude source bounds; zero lever arm prevents rotational response from being injected as CoG translational acceleration"
        ),
        "shipping_sigma_band_hz": declared_band,
        "shipping_sigma_band_source_parity": parity,
        "directional_cross_axis_coupling": {
            "rank_one_complex_outer_product_retained_before_outer_bound": True,
            "response_spectral_matrix_PSD": True,
            "independent_cartesian_axis_boxes_used": False,
            "entrywise_cross_terms_constrained_by_PSD_Cauchy_Schwarz": True,
            "arbitrary_phase_is_covered": True,
        },
        "uniform_moment_theorem": {
            "envelope": "||h||_2 <= G min(1,(f_c/f)^p)",
            "parameter_quantifier": "for every G in [0,G_max], every f_c in [f_c,min,f_c,max], and every p>=2",
            "analytical_not_sampled": True,
            "unbanded_acceleration_moment_finite": True,
            "proof_corner": "G=G_max, f_c=f_c,max, p=2",
        },
        "worst_envelope_trace_upper_per_Hs2": worst,
        "flat_gain_6Hz_acceleration_trace_upper_per_Hs2_diagnostic": [down(0.0), flat_acc],
        "acceleration_moment_tightening_vs_flat_6Hz_corner_lower": reduction,
        "numerical_response_member_is_allowed_to_prune_P2": False,
        "left_physical_inclusion_obligation": (
            "SEA0 must still prove that the admitted physical sea/ship population lies inside this declared RAO envelope and that its finite-window realization satisfies the P1 Normal-Live hard source bounds. A Gaussian spectrum alone is not an infinite-time pointwise certificate."
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


def _inclusion_from_p2(response: dict, frontend: dict, p2: dict, repo: Path = REPO) -> dict:
    pf = P2I.validate(p2)
    ff = FRONT.validate(frontend)
    if pf:
        raise RuntimeError(f"P2 correlation interface invalid: {pf}")
    if ff:
        raise RuntimeError(f"wave-period front-end invalid: {ff}")
    if response.get("shipping_sigma_band_source_parity") is not True:
        raise RuntimeError("response bridge lost shipping source parity")
    if response.get("single_worst_envelope_proves_entire_parameter_box_by_monotonicity") is not True:
        raise RuntimeError("RAO continuum envelope theorem is missing")

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
        "qualification": "OU3_SEA3_TO_FROZEN_P2_RAO_FAMILY_NONPRUNING_INCLUSION",
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
        "RAO_parameter_box_consumed": response["rao_envelope_parameter_box"],
        "single_RAO_selected_for_inclusion": False,
        "estimator_frequency_is_clamped_before_tau_target": True,
        "sigma_variance_history_need_not_be_independently_bounded_for_right_inclusion": True,
        "R_S_transcendental_target_need_not_be_independently_bounded_for_right_inclusion": True,
        "finite_window_Gaussian_amplitude_bound_used_for_inclusion": False,
        "full_finite_memory_estimator_reachability_used_to_prune_P2": False,
        "P2_pruned_by_SEA3": False,
        "Lhat_SEA3_subset_L_current_source": passed,
        "SEA3_TO_P2_INCLUSION_CERTIFICATE": "PASS" if passed else "FAIL",
        "physical_SEA0_promoted_here": False,
        "interpretation": (
            "The right inclusion is uniform over the complete RAO envelope parameter box. It intentionally removes no P2 cells; the envelope is a robust family constraint, not one response trace."
        ),
        "next_obligation": (
            "close the physical left inclusion for the declared RAO/sea family; only if full-P2 H/A is inadequate, propagate response-weighted finite-memory estimator state to obtain SEA3-specific pruning"
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
        "qualification": "OU3_SEA3_HA_FEASIBILITY_BY_RAO_FAMILY_P2_SUPERSET_INHERITANCE",
        "H_dimension": 18,
        "A_dimension": 21,
        "H_relative_Riccati_injection_margin_lower": h,
        "A_relative_Riccati_injection_margin_lower": a,
        "full_P2_canonical_worst_H_A_margin": float(canonical.get("worst_H_A_margin", 0.0)),
        "unchanged_useful_gate": float(canonical.get("useful_gate", 0.0)),
        "full_P2_P3_CANONICAL_PASS": full_p2_pass,
        "SEA3_subset_inclusion_pass": inclusion_pass,
        "uniform_over_entire_RAO_parameter_box": True,
        "uniform_P2_PASS_is_inherited_by_SEA3_subset": True,
        "uniform_P2_FAIL_implies_SEA3_FAIL": False,
        "SEA3_specific_HA_recomputed_on_pruned_language": False,
        "SEA3_HA_FEASIBILITY": (
            "PASS_BY_P2_SUPERSET" if inherited else "INCONCLUSIVE_REQUIRES_SEA3_NARROWING"
        ),
        "SEA3_HA_feasible_by_existing_uniform_certificate": inherited,
        "P3_promoted_here": False,
        "next_obligation": (
            "existing uniform P2 H/A covers the complete SEA3 RAO envelope family; continue physical left inclusion and downstream P4/P5"
            if inherited
            else "derive response/finite-estimator source pruning and rerun the same frozen H/A theorem interface"
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
    if r.get("shipping_sigma_band_source_parity") is not True:
        f.append("shipping sigma-band source parity lost")
    if r.get("physical_SEA0_promoted") is not False:
        f.append("RAO response-family bridge was promoted as physical SEA0")
    if r.get("single_nominal_RAO_used") is not False:
        f.append("one nominal RAO was reintroduced")
    if r.get("finite_RAO_grid_used") is not False:
        f.append("finite RAO proof grid was reintroduced")
    if r.get("single_worst_envelope_proves_entire_parameter_box_by_monotonicity") is not True:
        f.append("RAO parameter-box monotonicity theorem missing")
    box = r.get("rao_envelope_parameter_box", {})
    if box.get("peak_translation_gain") != [0.0, 4.0]:
        f.append("certified RAO gain range changed")
    if box.get("rolloff_corner_hz") != [0.03, 1.2]:
        f.append("certified RAO roll-off-corner range changed")
    if float(box.get("high_frequency_rolloff_power_min", math.nan)) < 2.0:
        f.append("RAO roll-off is too weak for finite acceleration moment")
    coupling = r.get("directional_cross_axis_coupling", {})
    if coupling.get("rank_one_complex_outer_product_retained_before_outer_bound") is not True:
        f.append("directional cross-axis outer product was not retained")
    if coupling.get("independent_cartesian_axis_boxes_used") is not False:
        f.append("independent Cartesian response boxes were reintroduced")
    if coupling.get("arbitrary_phase_is_covered") is not True:
        f.append("arbitrary complex RAO phase is not covered")
    if r.get("uniform_moment_theorem", {}).get("unbanded_acceleration_moment_finite") is not True:
        f.append("response roll-off did not close the acceleration moment")

    p = d.get("p2_inclusion", {})
    if p.get("SEA3_TO_P2_INCLUSION_CERTIFICATE") != "PASS":
        f.append("SEA3->P2 inclusion did not pass")
    if p.get("Lhat_SEA3_subset_L_current_source") is not True:
        f.append("SEA3 source language was not included in current P2 source language")
    if p.get("P2_pruned_by_SEA3") is not False:
        f.append("non-pruning inclusion incorrectly claims P2 pruning")
    if p.get("single_RAO_selected_for_inclusion") is not False:
        f.append("P2 inclusion selected one RAO")
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
        if ha.get("uniform_over_entire_RAO_parameter_box") is not True:
            f.append("H/A inheritance is not uniform over the RAO parameter box")
        if ha.get("SEA3_HA_feasible_by_existing_uniform_certificate") is not expected:
            f.append("H/A inherited feasibility flag is inconsistent")
        status = "PASS_BY_P2_SUPERSET" if expected else "INCONCLUSIVE_REQUIRES_SEA3_NARROWING"
        if ha.get("SEA3_HA_FEASIBILITY") != status:
            f.append("H/A feasibility status is inconsistent")
        if float(ha.get("unchanged_useful_gate", math.nan)) != 1.0e-18:
            f.append("H/A useful gate changed from 1e-18")
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
        "RAO_parameter_box": d["response_enclosure"]["rao_envelope_parameter_box"],
        "worst_acceleration_trace_upper_per_Hs2": d["response_enclosure"]["worst_envelope_trace_upper_per_Hs2"]["acceleration"],
        "tightening_vs_flat_6Hz": d["response_enclosure"]["acceleration_moment_tightening_vs_flat_6Hz_corner_lower"],
        "SEA3_TO_P2_INCLUSION_CERTIFICATE": d["p2_inclusion"]["SEA3_TO_P2_INCLUSION_CERTIFICATE"],
        "P2_pruned": d["p2_inclusion"]["P2_pruned_by_SEA3"],
        "ha_feasibility": None if d["ha_feasibility"] is None else d["ha_feasibility"]["SEA3_HA_FEASIBILITY"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
