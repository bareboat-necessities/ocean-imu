"""Declared-domain coordinate retention of the same attached complete word.

Route 1 of `ou3_p4_storage_routes.py` converts every reachable state into one
scalar information storage and then converts that scalar back into a 30-degree
attitude excursion through the worst direction of the metric. That second
conversion charges the whole bias-driven storage -- which lives in velocity and
displacement -- to attitude, and is where its sufficient chart bound 6.41 (H18)
and 2.27 (A21) is manufactured. This experiment deletes the scalarization
instead of tightening it: the reachable set is propagated in the declared
physical coordinates and each coordinate group is compared against its own
declared bound.

Four initial sets are reported against the same word. The declared product of
operating-domain balls in `tools/stability/ou3_proof_operating_domain.json`,
with attitude at the same 30-degree Cayley radius route 1 uses; the closed
bias ball alone, which is route 1's budget in these coordinates; the same
declared set without its integral-displacement ball, which says how much of
the failure that one ball carries; and the common forcing template on its own.
Each carries one amplitude on the same physical forcing template.

`ellipsoid_retention` then reports the correlated alternative: the covariance
of the word's initial point, which is a single convex set and needs no
subadditive step at all.

Nothing here is fitted to a replay, no domain is reduced, no filter
coefficient changes, and no route is promoted.

For each prefix the report gives both sides of the enclosure: a certified
upper bound (subadditive over the initial groups) and an attained lower bound
(a maximizing unit functional). A group is retained only when its upper bound
stays inside its own declared ball; it is definitely violated only when the
attained lower bound leaves it. Anything between the two is an enclosure gap,
not a result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import ou3_p4_motion_gain as G
import ou3_p4_storage_routes as R

REPO = Path(__file__).resolve().parents[2]
DOMAIN = REPO / "tools/stability/ou3_proof_operating_domain.json"

# One relative margin for both the enclosure check and the violation
# classification. A group whose attained bound sits within it of its declared
# radius is at the boundary in binary64, not demonstrably outside it.
MARGIN = 1e-9

# The 21 error coordinates in their declared physical groups. Offsets follow
# Kalman3D_Wave_OU_III: attitude, gyro bias, v, p, S, latent a_w, accel bias.
GROUPS = (
    ("attitude", 0, "chart"),
    ("gyro_bias", 3, "gyro_bias_error_norm_upper_rad_s"),
    ("velocity", 6, "velocity_error_norm_upper_mps"),
    ("position", 9, "position_error_norm_upper_m"),
    ("integral_displacement", 12, "integral_displacement_error_norm_upper_m_s"),
    ("latent_acceleration", 15, "latent_acceleration_error_norm_upper_mps2"),
    ("accelerometer_bias", 18, "accelerometer_bias_error_norm_upper_mps2"),
)


def chart_radius(degrees=30.):
    """Cayley radius of the declared attitude chart, as route 1 uses it."""
    return 2*np.tan(np.deg2rad(degrees)/2)


def declared_radii():
    """The declared operating-domain radius of each coordinate group.

    Every value comes from `ou3_proof_operating_domain.json` except attitude,
    which uses the same 30-degree Cayley chart radius route 1 uses. None is
    fitted to a replay.
    """
    bounds = json.loads(DOMAIN.read_text())["startup"]["physical_handoff_coordinate_bounds"]
    radii = {}
    for name, _, key in GROUPS:
        radii[name] = chart_radius() if key == "chart" else float(bounds[key])
        if not np.isfinite(radii[name]) or radii[name] <= 0:
            raise ValueError("declared domain radius must be positive and finite: "+name)
    return radii


def attained_lower_bound(blocks, radii, forcing, restarts=24, iterations=400):
    """max_{|u|=1} sum_H r_H |B_H^T u| + |u.f|, by normalized fixed point.

    The objective is convex, so every fixed point is an attained value of the
    exact supremum and never overstates it.
    """
    rng = np.random.default_rng(20260908)
    n = len(forcing)
    best = 0.
    starts = [forcing.copy()] if np.linalg.norm(forcing) > 0 else []
    starts += [np.linalg.svd(b)[0][:, 0] for b in blocks.values()]
    starts += [rng.normal(size=n) for _ in range(restarts)]
    for start in starts:
        u = np.asarray(start, dtype=float)
        if np.linalg.norm(u) == 0:
            continue
        u /= np.linalg.norm(u)
        for _ in range(iterations):
            gradient = np.zeros(n)
            for name, block in blocks.items():
                projected = block.T @ u
                norm = np.linalg.norm(projected)
                if norm > 0:
                    gradient += radii[name]*(block @ projected)/norm
            gradient += forcing*np.sign(u @ forcing or 1.)
            norm = np.linalg.norm(gradient)
            if norm == 0:
                break
            candidate = gradient/norm
            if np.linalg.norm(candidate-u) < 1e-14:
                u = candidate
                break
            u = candidate
        value = sum(radii[name]*np.linalg.norm(block.T @ u) for name, block in blocks.items())
        best = max(best, value + abs(u @ forcing))
    return float(best)


def stacked_blocks(transitions, responses, offset, active):
    """Every prefix's row block and forcing for one output group."""
    rows = slice(offset, offset+3)
    starts = [start for name, start, _ in GROUPS if active[name] > 0]
    weights = np.array([active[name] for name, _, _ in GROUPS if active[name] > 0])
    forcing = np.stack([response[rows] for response in responses])
    # An initial set with every ball at zero is legal: it is the forcing
    # template on its own, and it stacks to an empty group axis rather than
    # to nothing at all.
    blocks = np.empty((len(transitions), len(starts), 3, 3))
    for index, transition in enumerate(transitions):
        for position, start in enumerate(starts):
            blocks[index, position] = transition[rows, start:start+3]
    return blocks, weights, forcing


def attained_every_prefix(blocks, weights, forcing, restarts=6, iterations=120):
    """The same functional as `attained_lower_bound`, at every prefix at once.

    Selecting only the prefix with the largest certified bound can miss a
    prefix whose attained bound is larger, and it is the attained bound that
    demonstrates a violation. This runs the fixed point on every prefix so the
    reported maximum is over all of them.
    """
    rng = np.random.default_rng(20260908)
    count = len(forcing)
    fallback = np.tile(np.array([1., 0., 0.]), (count, 1))
    starts = [forcing] + [np.tile(axis, (count, 1)) for axis in np.eye(3)]
    starts += [rng.normal(size=(count, 3)) for _ in range(restarts)]
    best = np.zeros(count)
    for start in starts:
        u = np.array(start, dtype=float)
        norm = np.linalg.norm(u, axis=1, keepdims=True)
        u = np.where(norm > 0, u/np.where(norm > 0, norm, 1.), fallback)
        for _ in range(iterations):
            projected = np.einsum("pgji,pj->pgi", blocks, u)
            norms = np.linalg.norm(projected, axis=2, keepdims=True)
            unit = np.divide(projected, norms, out=np.zeros_like(projected),
                             where=norms > 0)
            gradient = np.einsum("g,pgij,pgj->pi", weights, blocks, unit)
            sign = np.sign(np.einsum("pi,pi->p", u, forcing))
            gradient = gradient + forcing*np.where(sign == 0, 1., sign)[:, None]
            norm = np.linalg.norm(gradient, axis=1, keepdims=True)
            u = np.where(norm > 0, gradient/np.where(norm > 0, norm, 1.), u)
        projected = np.einsum("pgji,pj->pgi", blocks, u)
        value = np.einsum("g,pg->p", weights, np.linalg.norm(projected, axis=2))
        best = np.maximum(best, value + np.abs(np.einsum("pi,pi->p", u, forcing)))
    return best


def source_contributions(blocks, radii, forcing):
    """Per-source share of the subadditive bound, in the output group's unit."""
    shares = {name: radii[name]*float(np.linalg.norm(block, 2))
              for name, block in blocks.items()}
    shares["forcing_template"] = float(np.linalg.norm(forcing))
    return shares


def certified_upper_bound(blocks, radii, forcing):
    """Subadditive over the independent declared balls; never understates."""
    return float(sum(source_contributions(blocks, radii, forcing).values()))


def retention(transitions, responses, radii, keep=None):
    """Per-group retention of the declared domain over every prefix.

    The certified bound answers retention and the attained bound answers
    violation, and the two are maximized over the prefixes independently: the
    prefix that reaches furthest need not be the prefix whose subadditive
    bound is loosest. Each is reported with the prefix that attains it.
    """
    active = {name: (radii[name] if keep is None or name in keep else 0.)
              for name, _, _ in GROUPS}
    result = {}
    for name, offset, _ in GROUPS:
        blocks, weights, forcing = stacked_blocks(transitions, responses, offset, active)
        certified = (np.einsum("g,pg->p", weights, np.linalg.svd(blocks, compute_uv=False)[..., 0])
                     + np.linalg.norm(forcing, axis=1))
        certified_at = int(np.argmax(certified))
        attained = attained_every_prefix(blocks, weights, forcing)
        # Refine with the multi-restart search at the two prefixes that can
        # carry the maximum, so the reported value is never below the coarse
        # sweep and the certified prefix is always evaluated.
        attained_at = int(np.argmax(attained))
        best = float(attained[attained_at])
        for index in {attained_at, certified_at}:
            named = {source: transitions[index][offset:offset+3, start:start+3]
                     for source, start, _ in GROUPS if active[source] > 0}
            refined = attained_lower_bound(named, active, responses[index][offset:offset+3])
            if refined > best:
                best, attained_at = refined, index
        upper = float(certified[certified_at])
        if best > upper*(1+MARGIN):
            raise ValueError("attained bound exceeded its certified enclosure: "+name)
        named = {source: transitions[certified_at][offset:offset+3, start:start+3]
                 for source, start, _ in GROUPS if active[source] > 0}
        result[name] = {
            "certified_prefix_index": certified_at,
            "attained_prefix_index": attained_at,
            "declared_radius": radii[name],
            "certified_upper_bound": upper,
            "attained_lower_bound": best,
            "certified_retention_ratio": upper/radii[name],
            "attained_retention_ratio": best/radii[name],
            "source_contributions": source_contributions(
                named, active, responses[certified_at][offset:offset+3])}
        result[name]["dominant_source"] = max(result[name]["source_contributions"],
                                              key=result[name]["source_contributions"].get)
    return result


def ellipsoid_retention(transitions, responses, radii, covariance):
    """Retention of the filter's own correlated initial set, exactly.

    The declared product box lets every coordinate sit at its own bound
    independently, which this filter's own kinematics cannot produce: the
    integral state is the running integral of the position state. The
    covariance of the word's initial point is the correlated alternative --
    it carries exactly the position/integral/attitude cross terms the box
    discards. It is also one convex set rather than a product, so its image
    needs no subadditive step and the bound below is exact.

    For `{x : x^T P^-1 x <= c}` with `P = L L^T`, the reachable excursion of
    group G at prefix j is `sqrt(c)*||Pi_G T_j L||_2 + |alpha|*|Pi_G r_j|`.
    The critical level is the largest `sqrt(c)` keeping every group inside its
    declared radius at every prefix, in units of the initial covariance, so 1
    is the filter's own one-sigma set.

    This exchanges one unproved premise for another: the box was never a
    reachable set, and this ellipsoid is the covariance the filter believes
    rather than a qualified bound on its actual error. The runtime audit
    records that actual covariances differ from the frozen P3 premises. The
    level below is therefore a requirement on covariance consistency, not a
    certificate.
    """
    root = np.linalg.cholesky(G.sym(covariance))
    report, levels = {}, {}
    for name, offset, _ in GROUPS:
        rows = slice(offset, offset+3)
        gains = np.array([float(np.linalg.norm(transition[rows] @ root, 2))
                          for transition in transitions])
        forcing = np.array([float(np.linalg.norm(response[rows]))
                            for response in responses])
        headroom = radii[name] - forcing
        # A zero gain places no bound on the level, so it cannot be limiting.
        allowed = np.divide(headroom, gains, out=np.full_like(gains, np.inf),
                            where=gains > 0)
        # A prefix whose forcing alone leaves the ball admits no level at all.
        level = 0. if np.any(headroom <= 0) else float(np.min(allowed))
        worst = int(np.argmin(allowed)) if level > 0 else int(np.argmax(forcing))
        marginal = float(np.max(np.linalg.eigvalsh(covariance[rows, rows])))
        report[name] = {
            "declared_radius": radii[name],
            "critical_initial_sigma_level": level,
            "limiting_prefix_index": worst,
            "gain_at_limiting_prefix": float(gains[worst]),
            "forcing_at_limiting_prefix": float(forcing[worst]),
            "maximum_forcing_excursion": float(np.max(forcing)),
            # How many initial standard deviations the declared radius is in
            # this coordinate, which is what makes the two sets comparable.
            "declared_radius_in_initial_sigma": (
                radii[name]/np.sqrt(marginal) if marginal > 0 else float("inf")),
            "one_sigma_prefix_index": int(np.argmax(gains+forcing)),
            "excursion_at_one_sigma": float(np.max(gains+forcing)),
            "retention_ratio_at_one_sigma": float(np.max(gains+forcing)/radii[name])}
        levels[name] = level
    limiting = min(levels, key=levels.get)
    return {"groups": report,
            "critical_initial_sigma_level": levels[limiting],
            "limiting_group": limiting,
            "exact_no_subadditive_step": True,
            "covariance_consistency_is_an_unproved_premise": True,
            "P4_PASS": False}


def audit_mode(root, rows, points, mode):
    """Run every initial set of this experiment against one attached word.

    Reports the four product-set subsets and the correlated covariance
    ellipsoid, together with the limiting group and the largest box of the
    declared shape the word retains. Promotion flags stay false throughout.
    """
    steps, _, _, defects, counts = G.build_word(root, rows, points, mode)
    if defects.failures:
        raise ValueError("finite factorization failed: "+repr(defects.failures[:1]))
    events = [r for r in rows if r["word"] == mode]
    transitions, responses, _, _ = R.compose(steps, events)
    point = next(r for r in points if r["word"] == mode)
    _, covariance = R.C.SOURCE.error_and_covariance(point, 21)
    radii = declared_radii()
    names = [name for name, _, _ in GROUPS]
    subsets = {
        # The declared product box exactly as written.
        "full_declared_initial_set": None,
        # Route 1's budget in these coordinates: bias ball plus the template.
        "bias_ball_and_template_only": {"accelerometer_bias"},
        # Diagnostic only. It reports how much of the failure the declared
        # 300 m*s integral ball carries; it does not reduce any domain.
        "declared_set_without_integral_displacement_ball":
            {n for n in names if n != "integral_displacement"},
        # The common forcing template alone, with every initial ball at zero.
        # This is the globally bounded particular solution the critic listed:
        # if the response were the obstruction it would show here.
        "forcing_template_only": set(),
    }
    result = {"counts": dict(counts),
              "declared_radii": radii,
              "declared_domain_sha256": hashlib.sha256(DOMAIN.read_bytes()).hexdigest(),
              "scalarized_storage_conversion_used": False,
              "operating_domain_reduced": False,
              "P4_PASS": False}
    for label, keep in subsets.items():
        result[label] = retention(transitions, responses, radii, keep=keep)
    full = result["full_declared_initial_set"]
    limiting = max(full, key=lambda k: full[k]["certified_retention_ratio"])
    ratio = full[limiting]["certified_retention_ratio"]
    result.update({
        "limiting_group": limiting,
        "limiting_certified_retention_ratio": ratio,
        # Every bound is positively homogeneous in the declared radii and the
        # template amplitude, so one similarity factor scales them all.
        "largest_retained_similar_box_fraction": 1./ratio,
        "declared_domain_retained": all(
            item["certified_retention_ratio"] <= 1. for item in full.values()),
        "violation_margin": MARGIN,
        "definitely_violated_groups": sorted(
            name for name, item in full.items()
            if item["attained_retention_ratio"] > 1.+MARGIN)})
    result["covariance_ellipsoid_initial_set"] = ellipsoid_retention(
        transitions, responses, radii, covariance)
    return result


def main():
    """Verify the capture hashes, audit both modes and write the evidence.

    The trace must recover its baseline bit for bit and every hash must match
    the recorded attachment, so a detached or regenerated capture fails here
    rather than producing a report.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--attachment", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    attachment = json.loads(args.attachment.read_text())
    paths = {suffix: Path(str(args.prefix)+suffix) for suffix in
             (".root.json", ".inputs.csv", ".prefixes.jsonl", ".events.jsonl")}
    hashes = {suffix: hashlib.sha256(path.read_bytes()).hexdigest() for suffix, path in paths.items()}
    if not attachment["read_only_trace_recovers_baseline_bit_for_bit"]:
        raise ValueError("passive source trace parity is required")
    if hashes[".events.jsonl"] != attachment["event_trace_sha256"]:
        raise ValueError("detached event trace")
    if any(hashes[k] != attachment["baseline_capture_sha256"][k] for k in
           (".root.json", ".inputs.csv", ".prefixes.jsonl")):
        raise ValueError("detached source/root/prefix")
    rows = [json.loads(line) for line in paths[".events.jsonl"].read_text().splitlines()]
    points = [json.loads(line) for line in paths[".prefixes.jsonl"].read_text().splitlines()]
    root = json.loads(paths[".root.json"].read_text())
    report = {"experiment": "DECLARED_DOMAIN_COORDINATE_RETENTION",
              "capture_sha256": hashes,
              "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "frozen_coefficients_only": True, "same_history_forcing_preserved": True,
              "operating_domain_reduced": False, "replay_fitted_radius": False,
              "physical_source_admission_pass": False,
              "P4_PASS": False, "P5_MAY_START": False, "modes": {}}
    for mode in ("H18", "A21"):
        if attachment["modes"][mode]["decision"] != "CONNECTED_POINT_ATTACHMENT_PASS":
            raise ValueError("unattached word: "+mode)
        report["modes"][mode] = audit_mode(root, rows, points, mode)
        print("DOMAIN_RETENTION", mode,
              report["modes"][mode]["limiting_group"],
              report["modes"][mode]["limiting_certified_retention_ratio"], flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
