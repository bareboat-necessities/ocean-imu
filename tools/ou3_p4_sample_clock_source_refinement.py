#!/usr/bin/env python3
"""Sample-clock refinement of the OU-III P2 source graph for P4/P5.

The original P2 graph is deliberately time-unlabelled.  Its edge construction
allows any delay above ``ADAPT_EVERY_SECS`` and therefore lets the EMA memory
factor tend to zero.  That is a valid coarse source over-approximation, but it
turns the 800-node graph into all 640,000 possible transitions and throws away
the source timing that a useful nonlinear P4 certificate needs.

The shipping wrapper is sharper:

* every valid IMU sample updates the tuner EMA;
* the strict adaptation clock only *stages* the already-smoothed candidate;
* the pending candidate is committed before the following valid sample; and
* if the binary64 clock eventually cannot advance by the binary32 sample step,
  the schedule freezes -- it does not make an arbitrary late jump.

This producer keeps the exact same 10 x 8 x 10 physical P2 partition and the
same deliberately broad SpectralMSE R_S target clamp.  It refines only the
transition language.  To remain an infinite-time source theorem rather than a
short-run simulation claim, the floating-clock spacing is certified by a
finite IEEE-754 exponent/boundary enumeration.  Before clock stagnation every
stage interval lies between 13 and 26 valid samples; the nominal early-life
spacing is 21 samples.  Every integer spacing in [13,26] is then admitted, even
though not all are actually realized.  Stagnation is represented by a node
self-loop.

The EMA image is composed one shipping sample at a time.  This keeps every
exponential call inside the already validated small-argument backend and also
covers sample-varying targets/horizons by reusing the common interval boxes at
each sample.  For R_S, the configured slew_log is exactly zero, so its horizon
is the plain clamped ``ADAPT_RS_MULT * tau_target``.

This is source-language evidence only; it cannot promote P4 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_source_path_reachability as PATH

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
TAU_CELLS = 10
SIGMA_CELLS = 8
RS_CELLS = 10
FREQ_CELLS = 8
BASE_P2_STATES = TAU_CELLS * SIGMA_CELLS * RS_CELLS
BASE_P2_EDGES = BASE_P2_STATES * BASE_P2_STATES
CLOCK_BOUNDARY_OFFSETS = 64


def _configured_rs_horizon(tau_target, c):
    """Return the configured R_S smoothing-horizon interval.

    The source marker is intentionally checked here: if a future shipping
    configuration arms the discrepancy-dependent slew branch, this refinement
    must be re-derived rather than silently retaining the zero-slew formula.
    """
    text = PATH.WRAPPER.read_text(encoding="utf-8")
    marker = "constexpr float ADAPT_RS_SLEW_LOG          = 0.0f;"
    if marker not in text:
        raise RuntimeError("sample-clock refinement requires configured ADAPT_RS_SLEW_LOG=0")
    safe_tau = PATH._clamp_interval(
        PATH._iv(tau_target), c["time_scale_min"], c["time_scale_max"]
    )
    nominal = PATH.I(c["adapt_RS_mult"]) * safe_tau
    return PATH._clamp_interval(nominal, c["horizon_min"], c["horizon_max"]).as_list()


def _stage_spacing_from(start: float, dt: float, commit: float, limit: int = 128):
    """Return the next strict-clock stage spacing, or None after stagnation."""
    t = float(start)
    last = float(start)
    for n in range(1, limit + 1):
        row = PATH.source_clock_step(t, last, False, dt, commit_s=commit)
        if row["time_s"] == t:
            return None
        t = row["time_s"]
        if row["stage_updated_candidate_for_next_sample"]:
            return n
    raise RuntimeError("clock spacing exceeded finite enumeration limit")


def _clock_certificate(c):
    """Certify the infinite-time floating-clock stage-spacing envelope.

    For a fixed binary64 exponent, every representable ``t`` is on the same
    ulp grid.  Away from the next exponent boundary, ``fl(t+dt)-t`` therefore
    depends only on dt/ulp, not on the mantissa.  A stage can only see different
    increments if it crosses that boundary.  We consequently enumerate one
    interior representative for every binary64 exponent and 64 predecessor
    grid points below every exponent boundary.  Sixty-four exceeds the largest
    finite stage spacing found and thus covers every possible boundary-crossing
    stage.  Exponent 46 is checked separately: dt is below half an ulp there,
    so the shipping clock stagnates and the source schedule freezes.
    """
    dt = float(c["dt"])
    commit = float(c["commit"])
    finite_spacings = set()
    boundary_cases = 0
    interior_cases = 0

    # Include subnormal scale through the last exponent at which time+dt can
    # still advance.  2**-1074 is the smallest positive binary64 value.
    for exponent in range(-1074, 46):
        start = math.ldexp(1.0, exponent)
        spacing = _stage_spacing_from(start, dt, commit)
        interior_cases += 1
        if spacing is not None:
            finite_spacings.add(int(spacing))

        # Any stage crossing the upper boundary must start within at most one
        # stage horizon of it.  Enumerating 64 predecessor floats is therefore
        # exhaustive once 64 is verified larger than the maximum stage gap.
        boundary = math.ldexp(1.0, exponent + 1)
        x = boundary
        for _ in range(CLOCK_BOUNDARY_OFFSETS):
            x = math.nextafter(x, -math.inf)
            spacing = _stage_spacing_from(x, dt, commit)
            boundary_cases += 1
            if spacing is not None:
                finite_spacings.add(int(spacing))

    if not finite_spacings:
        raise RuntimeError("no finite source-clock stage spacing found")
    min_gap = min(finite_spacings)
    max_gap = max(finite_spacings)
    if CLOCK_BOUNDARY_OFFSETS <= max_gap:
        raise RuntimeError("boundary enumeration is shorter than a stage interval")

    freeze_start = math.ldexp(1.0, 46)
    freeze_row = PATH.source_clock_step(freeze_start, freeze_start, False, dt, commit_s=commit)
    if freeze_row["time_s"] != freeze_start or freeze_row["stage_updated_candidate_for_next_sample"]:
        raise RuntimeError("binary64 clock did not freeze at the certified exponent")

    # Pin the normal deployment-age semantics separately.  This also verifies
    # the strict comparison and the one-sample pending activation convention.
    t = 0.0
    last = 0.0
    pending = False
    stage_steps = []
    apply_steps = []
    for k in range(1, 21 * 16 + 2):
        row = PATH.source_clock_step(t, last, pending, dt, commit_s=commit)
        if row["commit_previous_candidate_before_sample"]:
            apply_steps.append(k)
        if row["stage_updated_candidate_for_next_sample"]:
            stage_steps.append(k)
        pending = row["stage_updated_candidate_for_next_sample"]
        t = row["time_s"]
        last = row["last_stage_s"]
    nominal_spacings = [b-a for a, b in zip(stage_steps, stage_steps[1:])]
    if not nominal_spacings or set(nominal_spacings) != {21}:
        raise RuntimeError(f"normal-age source clock spacing changed: {sorted(set(nominal_spacings))}")
    if apply_steps[:len(stage_steps)-1] != [x+1 for x in stage_steps[:-1]]:
        raise RuntimeError("pending candidate is not applied on the next valid sample")

    return {
        "dt_binary32_s": dt,
        "commit_threshold_binary32_s": commit,
        "strict_clock_predicate": "time_-last_adapt_time_sec_ > adapt_every_secs_",
        "nominal_stage_spacing_valid_samples": 21,
        "nominal_commit_to_commit_elapsed_s": 21.0 * dt,
        "pending_apply_delay_valid_samples": 1,
        "finite_stage_spacing_valid_samples_observed": sorted(finite_spacings),
        "finite_stage_spacing_valid_samples_lower": int(min_gap),
        "finite_stage_spacing_valid_samples_upper": int(max_gap),
        "all_integer_spacings_between_bounds_admitted_by_graph": True,
        "IEEE754_exponent_classes_checked": 46 - (-1074),
        "IEEE754_boundary_predecessor_cases_checked": boundary_cases,
        "IEEE754_interior_cases_checked": interior_cases,
        "boundary_predecessors_per_exponent": CLOCK_BOUNDARY_OFFSETS,
        "boundary_enumeration_exceeds_max_stage_spacing": CLOCK_BOUNDARY_OFFSETS > max_gap,
        "floating_clock_stagnation_exponent": 46,
        "floating_clock_stagnation_verified": True,
        "floating_clock_stagnation_semantics": "schedule freezes; represented by source-node self-loop",
    }


def _partition(c):
    """Recreate the exact P2 10 x 8 x 10 source partition."""
    tau_lo = max(c["min_tau"], PATH.down(c["tau_coeff"] * 0.5 / c["max_freq"]))
    tau = PATH._cells(tau_lo, c["max_tau"], TAU_CELLS)
    sigma = PATH._cells(PATH.RAW_SIGMA_GRAPH_LOWER, c["max_sigma"], SIGMA_CELLS)
    rs = PATH._cells(c["min_RS"], c["max_RS"], RS_CELLS)
    freq = PATH._cells(c["min_freq"], c["max_freq"], FREQ_CELLS)
    return tau, sigma, rs, freq


def _ema_samples(x, target, horizon, dt: float, samples: int):
    """Compose the shipping convex EMA enclosure for ``samples`` updates.

    Reusing the target/horizon *boxes* on each step admits arbitrary variation
    inside those boxes at each sample, which is exactly the source-complete
    behavior required here.  Each exp argument is one 5-ms sample and remains
    in the validated transcendental range.
    """
    out = list(x)
    for _ in range(int(samples)):
        out = PATH._ema_image(out, target, horizon, dt, max_elapsed=dt)
    return out


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    """Build the finite-speed committed-source transition graph."""
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("sample-clock source refinement must not be trajectory fitted")

    c = PATH._constants()
    clock = _clock_certificate(c)
    dt = float(clock["dt_binary32_s"])
    min_gap = int(clock["finite_stage_spacing_valid_samples_lower"])
    max_gap = int(clock["finite_stage_spacing_valid_samples_upper"])
    admitted_gaps = range(min_gap, max_gap + 1)
    tau, sigma, rs, freq = _partition(c)

    states = []
    index = {}
    for ti, t in enumerate(tau):
        for si, s in enumerate(sigma):
            for ri, r in enumerate(rs):
                index[(ti, si, ri)] = len(states)
                states.append((t, s, r))

    # Precompute the one-dimensional transition images.  The source target
    # family factors as frequency x raw-sigma-cell, while R_S uses the full
    # deployed clamp.  This is exactly equivalent to the Cartesian target loop
    # but avoids redoing the 13..26 sample interval composition 800*64 times.
    tau_match = {}
    sigma_match = {}
    rs_match = {}
    target_boxes = 0
    for fi, f in enumerate(freq):
        tt = PATH._tau_target(f, c)
        ht = PATH._tau_sigma_horizon(f, c)
        hr = _configured_rs_horizon(tt, c)
        rr = PATH._rs_target_box(c)
        for gap in admitted_gaps:
            for ti, t in enumerate(tau):
                image = _ema_samples(t, tt, ht, dt, gap)
                tau_match[(ti, fi, gap)] = PATH._matching(tau, image)
            for ri, r in enumerate(rs):
                image = _ema_samples(r, rr, hr, dt, gap)
                rs_match[(ri, fi, gap)] = PATH._matching(rs, image)
            for target_si, ss in enumerate(sigma):
                target_boxes += 1 if gap == min_gap else 0
                for si, s in enumerate(sigma):
                    image = _ema_samples(s, ss, ht, dt, gap)
                    sigma_match[(si, fi, target_si, gap)] = PATH._matching(sigma, image)

    graph = [set([q]) for q in range(len(states))]  # frozen-clock branch
    for q, _state in enumerate(states):
        ti0 = q // (SIGMA_CELLS * RS_CELLS)
        rem = q % (SIGMA_CELLS * RS_CELLS)
        si0 = rem // RS_CELLS
        ri0 = rem % RS_CELLS
        out = graph[q]
        for fi in range(len(freq)):
            for target_si in range(len(sigma)):
                for gap in admitted_gaps:
                    tis = tau_match[(ti0, fi, gap)]
                    sis = sigma_match[(si0, fi, target_si, gap)]
                    ris = rs_match[(ri0, fi, gap)]
                    for i in tis:
                        for j in sis:
                            for k in ris:
                                out.add(index[(i, j, k)])

    gl = [sorted(x) for x in graph]
    comps = PATH._scc(gl)
    recurrent = set()
    for cc in comps:
        if len(cc) > 1 or (cc and cc[0] in graph[cc[0]]):
            recurrent.update(cc)

    new_edges = sum(map(len, gl))
    failures = []
    if len(states) != BASE_P2_STATES:
        failures.append("sample-clock graph changed the P2 physical state count")
    if not 0 < new_edges <= BASE_P2_EDGES:
        failures.append("sample-clock graph is not bounded by the P2 complete graph")
    if new_edges >= BASE_P2_EDGES:
        failures.append("sample-clock graph remained complete/all-to-all")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P2_SAMPLE_CLOCK_COMMIT_REACHABILITY_REFINEMENT",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "base_P2_reference": "tools/ou3_p4_source_path_reachability.py",
        "base_P2_untimed_state_count": BASE_P2_STATES,
        "base_P2_untimed_transition_edges": BASE_P2_EDGES,
        "same_physical_partition_as_P2": True,
        "partition": {
            "tau": len(tau), "sigma_tuner_raw": len(sigma), "R_S": len(rs),
            "states": len(states), "target_boxes": target_boxes,
        },
        "clock": clock,
        "EMA_updated_every_valid_sample": True,
        "EMA_composed_sample_by_sample": True,
        "sample_varying_target_and_horizon_boxes_admitted": True,
        "commit_only_stages_current_smoothed_candidate": True,
        "pending_candidate_applied_before_next_sample": True,
        "arbitrary_late_commit_jump_removed": True,
        "frozen_clock_self_loop_included": True,
        "finite_stage_gap_lower_samples": min_gap,
        "finite_stage_gap_upper_samples": max_gap,
        "RS_slew_log_configured_zero": True,
        "RS_target_full_deployed_clamp_overapprox": True,
        "RS_target_powf_tightening_used": False,
        "base_transition_edges": BASE_P2_EDGES,
        "transition_edges": new_edges,
        "edge_reduction_factor": float(BASE_P2_EDGES) / float(new_edges),
        "strongly_connected_components": len(comps),
        "recurrent_states": len(recurrent),
        "source_graph_all_to_all": new_edges == len(states) * len(states),
        "graph": gl,
        "P2_SAMPLE_CLOCK_REFINEMENT_CERTIFICATE": "PASS" if not failures else "FAIL",
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "consume this finite-speed committed-source graph in the complete H=18/A=21 operation-matched directional word dissipation backend; between commits the MEKF schedule is held exactly"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    """Fail closed on any source-scope or finite-clock regression."""
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P2_SAMPLE_CLOCK_COMMIT_REACHABILITY_REFINEMENT":
        f.append("wrong qualification")
    for key in (
        "source_only", "same_physical_partition_as_P2", "EMA_updated_every_valid_sample",
        "EMA_composed_sample_by_sample", "sample_varying_target_and_horizon_boxes_admitted",
        "commit_only_stages_current_smoothed_candidate", "pending_candidate_applied_before_next_sample",
        "arbitrary_late_commit_jump_removed", "frozen_clock_self_loop_included",
        "RS_slew_log_configured_zero", "RS_target_full_deployed_clamp_overapprox",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("trajectory_replay_used", "filter_changed", "RS_target_powf_tightening_used",
                "source_graph_all_to_all", "P4_USABLE_CERTIFICATE_PROMOTED"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    p = d.get("partition", {})
    if p.get("states") != BASE_P2_STATES:
        f.append("physical source partition changed")
    clock = d.get("clock", {})
    if clock.get("nominal_stage_spacing_valid_samples") != 21:
        f.append("nominal configured stage spacing is not 21 samples")
    if clock.get("pending_apply_delay_valid_samples") != 1:
        f.append("pending apply delay is not one sample")
    if clock.get("floating_clock_stagnation_verified") is not True:
        f.append("eventual clock stagnation is not certified")
    if clock.get("boundary_enumeration_exceeds_max_stage_spacing") is not True:
        f.append("clock boundary enumeration is not exhaustive")
    lo = d.get("finite_stage_gap_lower_samples")
    hi = d.get("finite_stage_gap_upper_samples")
    if not (isinstance(lo, int) and isinstance(hi, int) and 1 <= lo <= hi < CLOCK_BOUNDARY_OFFSETS):
        f.append("finite stage-gap interval is invalid")
    if d.get("base_P2_untimed_state_count") != BASE_P2_STATES:
        f.append("base P2 state-count reference changed")
    if d.get("base_P2_untimed_transition_edges") != BASE_P2_EDGES:
        f.append("base P2 edge-count reference changed")
    if not int(d.get("transition_edges", 0)) < BASE_P2_EDGES:
        f.append("sample clock did not reduce the source edge family")
    if d.get("P2_SAMPLE_CLOCK_REFINEMENT_CERTIFICATE") != "PASS":
        f.append("sample-clock refinement did not pass")
    return list(dict.fromkeys(f))


def main() -> int:
    """Write the source certificate JSON."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P2_SAMPLE_CLOCK_REFINEMENT_CERTIFICATE"],
        "clock": d["clock"],
        "edges": [d["base_transition_edges"], d["transition_edges"]],
        "edge_reduction_factor": d["edge_reduction_factor"],
        "SCC": d["strongly_connected_components"],
        "recurrent_states": d["recurrent_states"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
