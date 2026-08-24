#!/usr/bin/env python3
"""Validate rigorous OU-III source-word enclosure and promotion margins.

Schema 4 binds P3 and P4 to one source information geometry.  P4 uses the exact
Cayley coordinate c(R)=2 tan(theta/2)u on theta<pi and the metric

    W_g = s_m [c(R);xi]^T Sigma_KF(g)^-1 [c(R);xi],

where one positive scalar s_m is shared by every source node of a fixed mode.
That scalar is a numerical normalization only: it changes no physical level set
or generalized word contraction.  Full attitude--linear cross terms are kept.
The retired block-diagonal group metric, common-Euclidean route, replay-fitted
promotion and repeated one-sample contraction are not accepted as fallbacks.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_numerical_certificate as BASE
import ou3_source_domain_contract as SOURCE_DOMAIN

SCHEMA = 4
SOURCE_NOISE_SCHEMA = 1
ANCHOR_REL_TOL = 5.0e-6
REQUIRED_HYBRID_ROWS = set(SOURCE_DOMAIN.HYBRID_OBLIGATIONS) - {"periodic_aw_covariance_sync"}


def finite_positive(value) -> bool:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(x) and x > 0.0


def finite_nonnegative(value) -> bool:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(x) and x >= 0.0


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _validate_cayley_information_metric(payload: dict, sigma_lo, sigma_hi) -> tuple[dict, list[str]]:
    failures: list[str] = []
    if payload.get("kind") != "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC":
        failures.append("path metric kind is not CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC")
    for key, msg in (
        ("source_covariance_inverse", "Cayley metric is not source covariance information geometry"),
        ("node_dependent", "Cayley metric is not source/node dependent"),
        ("full_attitude_linear_cross_terms_retained", "Cayley metric discards attitude-linear cross terms"),
        ("local_coordinate_matches_P3_delta_theta", "Cayley local coordinate does not match P3 delta-theta"),
        ("local_quadratic_is_positive_scalar_multiple_of_P3_information_metric", "P4 local quadratic is not the P3 information geometry"),
        ("endpoint_metric_must_match_endpoint_source_covariance", "endpoint metric/source covariance correlation is not required"),
        ("joint_source_reachability_required", "Cayley metric does not require joint source reachability"),
        ("same_scale_on_every_source_node_in_mode", "P4 permits arbitrary node-wise metric rescaling"),
    ):
        if payload.get(key) is not True:
            failures.append(msg)
    if payload.get("block_diagonal_metric_used") is not False:
        failures.append("retired block-diagonal metric is active")
    if payload.get("common_Euclidean_metric_used") is not False:
        failures.append("common Euclidean metric fallback is active")
    chart = str(payload.get("chart_domain", ""))
    if "theta" not in chart or "pi" not in chart:
        failures.append("Cayley chart theta<pi is not explicit")

    scale = payload.get("mode_global_positive_scale")
    mlo = payload.get("metric_lambda_min_lower")
    mhi = payload.get("metric_lambda_max_upper")
    if not finite_positive(scale):
        failures.append("mode_global_positive_scale is not finite positive")
    if not finite_positive(mlo):
        failures.append("metric_lambda_min_lower is not finite positive")
    if not finite_positive(mhi):
        failures.append("metric_lambda_max_upper is not finite positive")
    if finite_positive(mlo) and finite_positive(mhi) and float(mhi) < float(mlo):
        failures.append("metric eigenvalue upper bound is below lower bound")

    # M=s_m Sigma^-1.  Reported source-uniform bounds may be wider, never tighter.
    if finite_positive(scale) and finite_positive(sigma_hi) and finite_positive(mlo):
        exact_lower = float(scale) / float(sigma_hi)
        tol = ANCHOR_REL_TOL * max(abs(exact_lower), 1.0e-300)
        if float(mlo) > exact_lower + tol:
            failures.append("metric lower bound is optimistic relative to s_mode/Sigma upper")
    if finite_positive(scale) and finite_positive(sigma_lo) and finite_positive(mhi):
        exact_upper = float(scale) / float(sigma_lo)
        tol = ANCHOR_REL_TOL * max(abs(exact_upper), 1.0e-300)
        if float(mhi) + tol < exact_upper:
            failures.append("metric upper bound is optimistic relative to s_mode/Sigma lower")

    return {
        "kind": payload.get("kind"),
        "chart_coordinate": payload.get("chart_coordinate"),
        "chart_domain": payload.get("chart_domain"),
        "exact_group_metric": payload.get("exact_group_metric"),
        "source_covariance_inverse": payload.get("source_covariance_inverse"),
        "node_dependent": payload.get("node_dependent"),
        "mode_global_positive_scale": scale,
        "same_scale_on_every_source_node_in_mode": payload.get("same_scale_on_every_source_node_in_mode"),
        "full_attitude_linear_cross_terms_retained": payload.get("full_attitude_linear_cross_terms_retained"),
        "block_diagonal_metric_used": payload.get("block_diagonal_metric_used"),
        "common_Euclidean_metric_used": payload.get("common_Euclidean_metric_used"),
        "local_quadratic_is_positive_scalar_multiple_of_P3_information_metric": payload.get("local_quadratic_is_positive_scalar_multiple_of_P3_information_metric"),
        "endpoint_metric_must_match_endpoint_source_covariance": payload.get("endpoint_metric_must_match_endpoint_source_covariance"),
        "metric_lambda_min_lower": mlo,
        "metric_lambda_max_upper": mhi,
    }, failures


def validate_mode(mode: str, payload: dict, contract_mode: dict,
                  sampled_mode: dict | None = None) -> dict:
    failures: list[str] = []
    if payload.get("source_complete") is not True:
        failures.append("source family is not declared source-complete")
    if payload.get("outward_rounded") is not True:
        failures.append("mode bounds are not declared outward-rounded")
    if payload.get("joint_source_reachability") is not True:
        failures.append("mode enclosure is not jointly source-reachable")
    if payload.get("one_sample_decrease_used") is not False:
        failures.append("one-sample decrease/repeated-step shortcut is active")
    if payload.get("source_replay_used", False) is not False:
        failures.append("source replay is used to establish P4")

    horizon = payload.get("word_horizon_s")
    if not finite_positive(horizon):
        failures.append("word_horizon_s is not finite positive")
    eta = payload.get("word_endpoint_relative_Riccati_injection_margin_lower")
    sigma_lo = payload.get("Sigma_lambda_min_lower")
    sigma_hi = payload.get("Sigma_lambda_max_upper")
    prefix = payload.get("prefix_information_gain_upper")
    if not finite_positive(eta) or not float(eta) < 1.0:
        failures.append("word-endpoint Riccati injection lower bound is not in (0,1)")
    if not finite_positive(sigma_lo):
        failures.append("Sigma_lambda_min_lower is not strictly positive")
    if not finite_positive(sigma_hi):
        failures.append("Sigma_lambda_max_upper is not finite positive")
    if finite_positive(sigma_lo) and finite_positive(sigma_hi) and float(sigma_hi) < float(sigma_lo):
        failures.append("Sigma eigenvalue upper bound is below lower bound")
    if not finite_positive(prefix):
        failures.append("prefix_information_gain_upper is not finite positive")

    linear_failures = list(failures)
    metric, metric_failures = _validate_cayley_information_metric(
        payload.get("path_metric", {}), sigma_lo, sigma_hi
    )
    failures.extend(metric_failures)

    theta = payload.get("theta_star")
    level = payload.get("certified_level_W")
    decrease = payload.get("endpoint_relative_W_decrease_lower")
    mu = payload.get("mu_W_lower")
    if not finite_positive(theta) or not float(theta) < math.pi:
        failures.append("theta_star is not in (0,pi)")
    if not finite_positive(level):
        failures.append("certified_level_W is not finite positive")
    if not finite_positive(decrease) or not float(decrease) < 1.0:
        failures.append("endpoint_relative_W_decrease_lower is not in (0,1)")
    if not finite_positive(mu):
        failures.append("mu_W_lower is not finite positive")
    elif finite_positive(decrease) and finite_positive(metric.get("metric_lambda_min_lower")):
        derived = down(float(decrease) * float(metric["metric_lambda_min_lower"]))
        if float(mu) > derived and not math.isclose(float(mu), derived, rel_tol=5e-15, abs_tol=0.0):
            failures.append("mu_W_lower exceeds decrease*metric_lambda_min_lower")
    if payload.get("all_word_prefixes_safe") is not True:
        failures.append("nonlinear word prefixes are not validated safe")
    if payload.get("accepted_correction_uses_source_series_branch") is not True:
        failures.append("accepted correction branch is not the certified deployed series branch")
    q = payload.get("prefix_canonical_error_norm_upper")
    qlim = payload.get("cayley_norm_limit")
    if not finite_nonnegative(q) or not finite_positive(qlim) or not float(q) < float(qlim):
        failures.append("Cayley prefix bootstrap does not stay inside chart")
    dcorr = payload.get("accepted_correction_norm_prefix_upper")
    if not finite_nonnegative(dcorr) or not float(dcorr) < 1.0e-2:
        failures.append("accepted correction can leave the certified quaternion branch")
    if mode == "A":
        proj = payload.get("active_bias_projection", {})
        if proj.get("projection_surface_reached_in_certified_funnel") is not False:
            failures.append("A-mode inner funnel reaches accelerometer-bias projection surface")
        if proj.get("exact_projection_branch_in_certified_funnel") != "identity_interior_branch":
            failures.append("A-mode projection branch is not exact interior identity")
    if sampled_mode and finite_positive(level):
        first_fail = sampled_mode.get("first_fail_W")
        if finite_positive(first_fail) and float(level) >= float(first_fail):
            failures.append(f"certified nonlinear level reaches/exceeds a sampled failure ({level} >= {first_fail})")

    # Preserve tiny strict gaps directly.  1-gap is diagnostic only when representable.
    ratio = None
    if finite_positive(decrease) and float(decrease) > math.ulp(1.0):
        candidate = up(1.0 - float(decrease))
        if candidate < 1.0:
            ratio = candidate

    return {
        "mode": mode,
        "linear_pass": not linear_failures,
        "nonlinear_pass": not failures,
        "pass": not failures,
        "failures": failures,
        "linear_failures": linear_failures,
        "word_horizon_s": horizon,
        "word_endpoint_relative_Riccati_injection_margin_lower": eta,
        "Sigma_lambda_min_lower": sigma_lo,
        "Sigma_lambda_max_upper": sigma_hi,
        "P3_inverse_covariance_conditioning_lambda_min_lower": 1.0/float(sigma_hi) if finite_positive(sigma_hi) else None,
        "P3_inverse_covariance_conditioning_lambda_max_upper": 1.0/float(sigma_lo) if finite_positive(sigma_lo) else None,
        "prefix_information_gain_upper": prefix,
        "path_metric": metric,
        "theta_star": theta,
        "endpoint_relative_W_decrease_lower": decrease,
        "endpoint_W_ratio_upper": ratio,
        "mu_W_lower": mu,
        "certified_level_W": level,
        "all_word_prefixes_safe": payload.get("all_word_prefixes_safe"),
    }


def validate_hybrid(payload: list[dict], modes: dict) -> dict:
    required = set(REQUIRED_HYBRID_ROWS)
    seen: set[str] = set()
    failures: list[str] = []
    rows = []
    for i, jump in enumerate(payload):
        kind = str(jump.get("kind", ""))
        if kind not in required:
            failures.append(f"hybrid[{i}] unknown or non-row obligation name: {kind}")
            continue
        seen.add(kind)
        row_failures: list[str] = []
        if jump.get("source_complete") is not True:
            row_failures.append("source family is not source-complete")
        if jump.get("outward_rounded") is not True:
            row_failures.append("bounds are not outward-rounded")
        source, gain = jump.get("source_level_W_upper"), jump.get("jump_gain_upper")
        additive, dest = jump.get("additive_W_upper"), jump.get("destination_level_W")
        new_coord = jump.get("new_coordinate_W_upper", 0.0)
        for label, value in (("source_level_W_upper",source),("jump_gain_upper",gain),
                             ("additive_W_upper",additive),("new_coordinate_W_upper",new_coord)):
            if not finite_nonnegative(value):
                row_failures.append(f"{label} missing/nonfinite/negative")
        if not finite_positive(dest):
            row_failures.append("destination_level_W missing/nonfinite/nonpositive")
        dest_mode = str(jump.get("destination_mode", ""))
        if dest_mode not in modes:
            row_failures.append("destination_mode must be H or A")
        elif finite_positive(dest):
            mode_level = modes[dest_mode].get("certified_level_W")
            if not finite_positive(mode_level) or float(dest) > float(mode_level):
                row_failures.append("destination level is outside certified same-mode level")
        if kind == "held_to_active":
            if jump.get("dimension_change_handled_by_embedding") is not True:
                row_failures.append("held_to_active must use explicit dimension embedding")
            if jump.get("source_dimension") != 18 or jump.get("destination_dimension") != 21:
                row_failures.append("held_to_active dimensions must be 18 -> 21")
            if "new_coordinate_W_upper" not in jump:
                row_failures.append("held_to_active must bound new-coordinate energy")
        if kind == "tilt_reset":
            if jump.get("discarded_pre_reset_tilt_excluded_from_multiplicative_gain") is not True:
                row_failures.append("tilt reset charges discarded tilt in multiplicative gain")
            if jump.get("reset_to_funnel_exact_map") is not True:
                row_failures.append("tilt reset is not exact reset-to-funnel map")
        if kind == "cooldown_reentry":
            if jump.get("reachable_word_product_used") is not True:
                row_failures.append("cooldown does not use reachable word products")
            if jump.get("global_worst_word_power_used") is not False:
                row_failures.append("cooldown uses global worst-word power")
        post = margin = None
        if (finite_nonnegative(source) and finite_nonnegative(gain) and finite_nonnegative(additive)
                and finite_nonnegative(new_coord) and finite_positive(dest)):
            post = float(gain)*float(source)+float(additive)+float(new_coord)
            margin = float(dest)-post
            if not margin > 0.0:
                row_failures.append("recomputed inward margin is not strictly positive")
        failures.extend(f"hybrid[{i}] {kind}: {x}" for x in row_failures)
        rows.append({"kind":kind,"destination_mode":dest_mode,"source_level_W_upper":source,
                     "jump_gain_upper":gain,"additive_W_upper":additive,"new_coordinate_W_upper":new_coord,
                     "post_jump_W_upper":post,"destination_level_W":dest,
                     "inward_margin_lower":margin,"pass":not row_failures})
    missing = sorted(required-seen)
    if missing:
        failures.append(f"missing hybrid obligations: {missing}")
    return {"pass":not failures,"failures":failures,"seen":sorted(seen),
            "analytic_obligation_separate":"periodic_aw_covariance_sync","bounds":rows}


def source_noise_moments(payload: dict) -> tuple[dict, list[str]]:
    failures=[]
    if payload.get("schema") != SOURCE_NOISE_SCHEMA:
        failures.append("source noise certificate schema mismatch")
    if payload.get("source_generated_not_trajectory_fit") is not True:
        failures.append("noise certificate is not source-generated")
    std=payload.get("standardized_increment",{})
    d=std.get("dimension")
    if not isinstance(d,int) or d<=0:
        failures.append("standardized primitive dimension is not positive integer"); d=0
    if std.get("covariance_upper_identity") is not True:
        failures.append("standardized covariance is not certified <= identity")
    return {"dimension":d,"Sigma_bar_norm_upper":1.0 if d else None,
            "trace_Sigma_bar_upper":float(d) if d else None,
            "trace_Sigma_bar_squared_upper":float(d) if d else None,
            "s2_upper":float(d) if d else None,
            "s4_upper":float(d*d+2*d) if d else None}, failures


def _geometric_square_sum(L: float, n: int) -> float:
    if n<=0: return 0.0
    q=L*L
    if abs(q-1.0)<=1e-12: return float(n)
    try: ans=(q**n-1.0)/(q-1.0)
    except OverflowError: return math.inf
    return ans if math.isfinite(ans) and ans>=0.0 else math.inf


def _gaussian_t_star(w_star: float, m: float, v: float, b: float) -> float|None:
    w2=w_star*w_star
    if not (w2>m and b>0.0 and v>=0.0): return None
    root=math.sqrt(v+2.0*b*(w2-m)); y=(root-math.sqrt(v))/(2.0*b)
    return y*y if y>0.0 and math.isfinite(y) else None


def _stochastic_lambda_upper(mode: dict) -> float|None:
    ratio=mode.get("endpoint_W_ratio_upper")
    if finite_nonnegative(ratio) and float(ratio)<1.0: return float(ratio)
    gap=mode.get("endpoint_relative_W_decrease_lower")
    if finite_positive(gap) and float(gap)>math.ulp(1.0):
        r=up(1.0-float(gap)); return r if r<1.0 else None
    return None


def validate_stochastic(payload: dict, modes: dict, source_noise: dict) -> dict:
    failures=[]
    for key,msg,want in (
        ("source_complete","stochastic source family is not source-complete",True),
        ("outward_rounded","stochastic sensitivity bounds are not outward-rounded",True),
        ("localization_prefix_safe","localized word prefixes are not validated safe",True),
        ("gaussian_localization_used","Gaussian localization is not explicitly certified",True),
        ("freedman_excursion_used","Freedman excursion concentration is not explicitly certified",True),
        ("markov_union_fallback_used","retired Markov/union stochastic fallback is active",False)):
        if payload.get(key) is not want: failures.append(msg)
    moments,nf=source_noise_moments(source_noise); failures.extend(nf)
    vals={}
    for name in ("L_X_upper","G_bar_upper","c_zw_upper","r_star_upper","c_ww_upper","g_W_upper","h_W_upper"):
        x=payload.get(name)
        if not finite_nonnegative(x): failures.append(f"{name} missing/nonfinite/negative")
        else: vals[name]=float(x)
    w_star=payload.get("localization_radius_standardized")
    if not finite_positive(w_star): failures.append("localization_radius_standardized is not finite positive")
    else: w_star=float(w_star)
    ell=payload.get("word_samples_upper"); N=payload.get("finite_horizon_words")
    if not isinstance(ell,int) or ell<=0: failures.append("word_samples_upper must be a positive integer"); ell=0
    if not isinstance(N,int) or N<=0: failures.append("finite_horizon_words must be a positive integer"); N=0
    a=payload.get("funnel_level_a"); W0=payload.get("W0_upper")
    if not finite_positive(a): failures.append("funnel_level_a is not finite positive")
    else:
        a=float(a); levels=[m.get("certified_level_W") for m in modes.values()]
        if any(not finite_positive(x) for x in levels) or a>min(float(x) for x in levels):
            failures.append("stochastic funnel level exceeds a certified deterministic mode level")
    if not finite_nonnegative(W0): failures.append("W0_upper missing/nonfinite/negative")
    elif finite_positive(a) and not float(W0)<float(a): failures.append("W0_upper must be strictly inside funnel_level_a")
    else: W0=float(W0)
    ratios=[_stochastic_lambda_upper(m) for m in modes.values()]; lambda_W=None
    if any(r is None for r in ratios): failures.append("deterministic mode contraction factors are unavailable at stochastic arithmetic precision")
    else: lambda_W=max(float(r) for r in ratios)
    dt=source_noise.get("physical_scales",{}).get("imu_dt_s"); horizons=[m.get("word_horizon_s") for m in modes.values()]
    if finite_positive(dt) and ell>0 and all(finite_positive(h) for h in horizons):
        minimum=math.ceil(max(float(h) for h in horizons)/float(dt)-1e-12)
        if ell<minimum: failures.append(f"word_samples_upper {ell} is below source horizon minimum {minimum}")
    else: failures.append("cannot connect word_samples_upper to source sample period")
    nu1=nuW=lams=sigma2=bW=vW=VW=drift=xstar=tstar=ploc=pfree=ptot=None
    if not failures:
        s2,s4=moments["s2_upper"],moments["s4_upper"]
        g=vals["G_bar_upper"]+vals["c_zw_upper"]*vals["r_star_upper"]
        nu1=2*g*g*s2+2*vals["c_ww_upper"]**2*s4
        nuW=float(ell)*nu1*_geometric_square_sum(vals["L_X_upper"],ell)
        if not math.isfinite(nuW): failures.append("localized source-word noise moment bound overflowed")
        else:
            lams=0.5*(1.0+lambda_W)
            sigma2=(vals["h_W_upper"]+vals["g_W_upper"]**2*lambda_W/(2*(1-lambda_W)))*nuW
            fixed=sigma2/(1-lams); lamN=lams**N; WN=lamN*W0+(1-lamN)*fixed
            drift=max(W0,WN); xstar=a-drift
            if not xstar>0.0: failures.append("stochastic drift envelope reaches/exceeds funnel boundary")
            bW=a; vW=a*a/4; VW=vW/(1-lams*lams)
            tstar=_gaussian_t_star(w_star,moments["trace_Sigma_bar_upper"],moments["trace_Sigma_bar_squared_upper"],moments["Sigma_bar_norm_upper"])
            if tstar is None: failures.append("localization radius does not exceed source Gaussian RMS threshold")
            if xstar>0.0 and tstar is not None:
                exponent=-xstar*xstar/(2*(VW+bW*xstar/3)); pfree=float(N)*math.exp(exponent); ploc=float(N*ell)*math.exp(-tstar); ptot=pfree+ploc
                if not (math.isfinite(ptot) and 0<=ptot<1): failures.append("recomputed finite-horizon stochastic failure probability is not < 1")
    return {"pass":not failures,"failures":failures,"source_noise_moments":moments,
            "lambda_W_upper":lambda_W,"lambda_s_upper":lams,"nu1_upper":nu1,"nu_W_upper":nuW,
            "sigma_s2_upper":sigma2,"funnel_level_a":a,"W0_upper":W0,"b_W_upper":bW,
            "v_W_upper":vW,"V_W_upper":VW,"drift_envelope_max_upper":drift,
            "excursion_margin_x_star_lower":xstar,"gaussian_t_star_lower":tstar,
            "localization_failure_probability_upper":ploc,"freedman_failure_probability_upper":pfree,
            "finite_horizon_failure_probability_upper":ptot,"word_samples_upper":ell,"finite_horizon_words":N}


def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument("--certificate-dir",type=Path,default=BASE.DEFAULT_OUT); ap.add_argument("--enclosure",type=Path,required=True); ap.add_argument("--source-noise",type=Path,default=None); args=ap.parse_args()
    cert=args.certificate_dir.resolve(); inp=json.loads(args.enclosure.read_text())
    if inp.get("schema")!=SCHEMA: raise RuntimeError(f"unsupported validated enclosure schema: {inp.get('schema')}")
    contract=json.loads((cert/"information_enclosure_contract.json").read_text())
    if contract.get("schema")!=SCHEMA: raise RuntimeError("information enclosure contract schema mismatch")
    source_noise_path=args.source_noise or (cert/"source_noise_certificate.json")
    source_noise=json.loads(source_noise_path.read_text())
    sampled={}; radius=cert/"neighborhood_radius_search"/"neighborhood_radius_search.json"
    if radius.exists(): sampled=json.loads(radius.read_text()).get("modes",{})
    prov=inp.get("provenance",{})
    provenance_ok=all(prov.get(k) is True for k in ("validated_arithmetic","outward_rounding","source_generated_not_trajectory_fit","continuous_source_coverage","joint_source_reachability","exact_deployed_quaternion_injection"))
    modes={m:validate_mode(m,inp.get("modes",{}).get(m,{}),contract.get("modes",{}).get(m,{}),sampled.get(m)) for m in ("H","A")}
    hybrid=validate_hybrid(inp.get("hybrid",[]),modes); stochastic=validate_stochastic(inp.get("stochastic",{}),modes,source_noise)
    linear=provenance_ok and all(x["linear_pass"] for x in modes.values()); nonlinear=provenance_ok and all(x["nonlinear_pass"] for x in modes.values()); hybrid_pass=nonlinear and hybrid["pass"]; stochastic_pass=hybrid_pass and stochastic["pass"]; deployment=linear and nonlinear and hybrid_pass and stochastic_pass
    out={"schema":SCHEMA,"metric_policy":{"P3":"word-endpoint generalized information inequality","P4":"mode-global positive scalar times exact Cayley lift of matching source Sigma_KF^-1","attitude_linear_cross_terms":"RETAINED","fallbacks":"NONE"},"validated_enclosure_provenance_pass":provenance_ok,"modes":modes,"hybrid":hybrid,"stochastic":stochastic,"continuous_linear_information_certificate":"PASS" if linear else "FAIL","numerical_neighborhood_certificate":"PASS" if nonlinear else "FAIL","hybrid_funnel_certificate":"PASS" if hybrid_pass else "FAIL","stochastic_certificate":"PASS" if stochastic_pass else "FAIL","deployment_theorem_certificate":"PASS" if deployment else "FAIL","promotion_is_machine_verified":True}
    (cert/"validated_enclosure_check.json").write_text(json.dumps(out,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({"continuous_linear_information_certificate":out["continuous_linear_information_certificate"],"numerical_neighborhood_certificate":out["numerical_neighborhood_certificate"],"deployment_theorem_certificate":out["deployment_theorem_certificate"],"H":modes["H"],"A":modes["A"]},indent=2,sort_keys=True))
    return 0 if deployment else 2


if __name__=="__main__":
    raise SystemExit(main())
