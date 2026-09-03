#!/usr/bin/env python3
"""Full-matrix translation lower transfer for one frozen P2-correlation segment.

This is the numerical primitive used by the source-varying P3 consumer.  It is
bound to ``OU3_P2_CORRELATED_STAGE_TRANSFER_V1`` and deliberately uses a fixed
physical scaling

    z = [v/h, p/h^2, S/h^3, a_w],

rather than a sigma-dependent state scaling.  A source transition therefore
changes OU process intensity but never silently rescales the state or metric.

The first implementation propagated an *interval enclosure* of the Riccati
posterior with

    P+ = P - P e (e' P e + R)^-1 e' P.

That formula is mathematically SPD, but natural interval subtraction loses the
repeated dependency between the two occurrences of P.  Over only a few samples
the interval diagonal could cross zero even though no true covariance did.
That is an enclosure-formulation failure, not a physical loss of excitation.

This implementation propagates a different object: one deterministic matrix L
that is a common Loewner lower bound for every covariance represented by the
source cell.  For an interval enclosure A of a symmetric prediction candidate,
write C for its binary64 midpoint matrix and let epsilon be an outward-rounded
maximum absolute row-sum bound for A-C.  Every represented symmetric A_real
then satisfies

    ||A_real-C||_2 <= epsilon,
    A_real >= C - epsilon I.

We therefore collapse every prediction to the point lower
``C-epsilon I`` before applying measurements.  Measurement monotonicity gives

    P >= L, R >= R_min  =>  M(P,R) >= M(L,R_min),

where M is the scalar-row Kalman covariance map.  The point posterior M(L,Rmin)
is itself evaluated with outward interval arithmetic and collapsed again by the
same spectral-radius construction.  Thus interval widths do not recursively
feed the subtractive Riccati formula.

For process intensity, the lower uses the committed source cell's sigma lower
endpoint.  Since Q(x) is PSD and sigma^2 is a positive scalar, this is a valid
Loewner lower for every sigma in that same physical P2 cell.  Tau uncertainty
remains interval-enclosed through x=h/tau; tau, sigma and R_S still come from
one common source node.

The two extra accepted translation measurements (a_w and S=0) at every sample
are deliberately stronger than the shipping branch language.  The optimal
posterior with those measurements is a lower covariance for any implemented
Kalman gain/accepted-rejected pattern, so proving this lower is sufficient for
P3 process excitation.

This module is one-segment proof machinery only.  It cannot promote P3/P4/P5 by
itself.
"""
from __future__ import annotations

import argparse
import functools
import json
import math
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_p2_correlation_path_memory as CORR
import ou3_p3_frozen_full_matrix_translation as FROZEN
import ou3_p3_scaled_process as SCALED
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 3
X_SUBCELLS = 4
MAX_ADAPTIVE_X_DEPTH = 12


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _point(x: float) -> Interval:
    return Interval.point(float(x))


def _scale(A, s: Interval):
    return [[s * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


@functools.lru_cache(maxsize=1)
def _acc_measurement_variance() -> float:
    vc = BASE.VECTOR.build()["configured_measurement_bounds"]
    return BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)


@functools.lru_cache(maxsize=1)
def _strongest_S_axis_factor() -> float:
    return float(min(BASE.source_rs_axis_std_factors()))


def _physical_measurement_variances(node: dict, h: float) -> tuple[float, float]:
    R_aw = _acc_measurement_variance()
    rs = Interval(*map(float, node["R_S_filter_std"]))
    rS = BASE.down(rs.lo * _strongest_S_axis_factor())
    # z_S=S/h^3, hence H_z=[0,0,h^3,0] for the physical S measurement.
    # Equivalently use coordinate z_S with R_z=R_S/h^6.
    R_S_z = BASE.down(BASE.down(rS * rS) / BASE.up(h ** 6))
    if not (R_aw > 0.0 and R_S_z > 0.0):
        raise RuntimeError("physical-scaled measurement variance lost positivity")
    return R_aw, R_S_z


def _node_x_subcells(node: dict, h: float, count: int = X_SUBCELLS):
    tau = Interval(*map(float, node["tau_s"]))
    xlo = BASE.down(h / tau.hi)
    xhi = BASE.up(h / tau.lo)
    return FROZEN._split_x(xlo, xhi, count)


def _center_and_radius(A):
    """Symmetric binary64 center/radius decomposition of an interval matrix."""
    A = FROZEN._sym(A)
    n = len(A)
    if n == 0 or any(len(row) != n for row in A):
        raise ValueError("square nonempty interval matrix required")
    center = [[0.0] * n for _ in range(n)]
    radius = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            a = A[i][j]
            # Stable midpoint; clamp protects a one-ulp endpoint accident.
            c = float(a.lo + 0.5 * (a.hi - a.lo))
            c = min(max(c, a.lo), a.hi)
            center[i][j] = c
            radius[i][j] = BASE.up(max(c - a.lo, a.hi - c))
    # A is symmetrically hulled, but make C exactly symmetric in binary64 too.
    for i in range(n):
        for j in range(i + 1, n):
            c = 0.5 * (center[i][j] + center[j][i])
            center[i][j] = c
            center[j][i] = c
            r = BASE.up(max(radius[i][j], radius[j][i]))
            radius[i][j] = r
            radius[j][i] = r
    return center, radius, n


def _absolute_row_sum_epsilon(radius, n: int) -> float:
    eps = 0.0
    for i in range(n):
        row = 0.0
        for j in range(n):
            row = BASE.up(row + radius[i][j])
        eps = max(eps, row)
    return BASE.up(eps)


def _scaled_row_sum_epsilon(radius, scale, n: int) -> float:
    """max_i sum_j r_ij/(d_i d_j), rounded toward +infinity.

    The denominator is rounded toward zero so the quotient stays an upper
    bound.  A denominator that underflows leaves no usable relative scaling.
    """
    eps = 0.0
    for i in range(n):
        row = 0.0
        for j in range(n):
            denom = BASE.down(scale[i] * scale[j])
            if not denom > 0.0:
                raise ArithmeticError("relative Loewner scaling underflowed")
            row = BASE.up(row + BASE.up(radius[i][j] / denom))
        eps = max(eps, row)
    return BASE.up(eps)


def _common_point_lower(A):
    """Return point L with A_real >= L for every symmetric A_real in A.

    ``A`` is first symmetrically hulled.  For each entry choose a binary64
    center c_ij inside its interval and an outward radius r_ij, so every
    represented symmetric A_real satisfies A_real=C+E with |E_ij|<=r_ij.

    For any positive diagonal D=diag(d_i),

        ||D^-1 E D^-1||_2 <= ||D^-1 E D^-1||_inf
                          <= max_i sum_j r_ij/(d_i d_j) = eps_D,

    hence D^-1 E D^-1 >= -eps_D I and therefore

        A_real >= C - eps_D D^2.

    The translation states are v/h, p/h^2, S/h^3 and a_w, whose certified
    covariance magnitudes differ by many orders.  A single absolute shift
    (D=I) is dominated by the largest block and destroys strict positivity of
    the smallest one, so the shift is taken in the natural relative scaling
    d_i=sqrt(c_ii) whenever the center diagonal is strictly positive.  Then
    the correction is exactly the relative shave

        L_ii = c_ii (1-eps_D),   L_ij = c_ij  (i!=j),

    which respects the dynamic range instead of flattening it.  D=I is
    retained as the fallback whenever some center diagonal is nonpositive,
    where no relative scaling exists.  Both branches are non-recursive and
    independently auditable, and the diagonal subtraction is rounded toward
    -infinity.
    """
    center, radius, n = _center_and_radius(A)
    scale = []
    relative = True
    for i in range(n):
        if center[i][i] <= 0.0:
            relative = False
            break
        scale.append(BASE.down(math.sqrt(center[i][i])))
        if not scale[-1] > 0.0:
            relative = False
            break

    if relative:
        try:
            eps = _scaled_row_sum_epsilon(radius, scale, n)
            shift = [BASE.up(eps * BASE.up(scale[i] * scale[i])) for i in range(n)]
        except ArithmeticError:
            relative = False
    if not relative:
        eps = _absolute_row_sum_epsilon(radius, n)
        shift = [eps] * n

    L = [[_point(center[i][j]) for j in range(n)] for i in range(n)]
    for i in range(n):
        L[i][i] = _point(BASE.down(center[i][i] - shift[i]))
    return FROZEN._sym(L), eps


def _strict_spd(A) -> bool:
    return bool(symmetric_positive_definite_ldlt(A)[0])


def _point_measurement_lower(L, coordinate: int, R_lower: float):
    """Lower one scalar optimal covariance update without recursive intervals."""
    if not _strict_spd(L):
        raise RuntimeError("Loewner measurement input is not strict SPD")
    # L is a point lower.  The interval expression encloses the exact M(L,R),
    # then the spectral collapse returns a deterministic matrix below it.
    post_interval = FROZEN._measurement_update(L, int(coordinate), float(R_lower))
    post, eps = _common_point_lower(post_interval)
    if not _strict_spd(post):
        raise RuntimeError("Loewner measurement lower lost strict SPD")
    return post, eps


def one_step(P, node: dict, x: Interval, h: float):
    """One rigorous strongest-measurement Loewner-lower covariance step."""
    # Interpret the incoming interval matrix only as an enclosure of a known
    # lower matrix and collapse it once.  Thereafter every operation returns a
    # point lower, so interval dependency cannot recursively accumulate.
    L0, _ = _common_point_lower(P)

    sigma = Interval(*map(float, node["sigma_filter_committed_mps2"]))
    sigma2_lower = BASE.down(sigma.lo * sigma.lo)
    if not sigma2_lower > 0.0:
        raise RuntimeError("source-cell sigma lower lost positivity")

    Fm = FROZEN._transition(x)
    Ft = FROZEN._transpose(Fm)
    # Qbase(x) is PSD.  sigma^2 >= sigma_lo^2, hence using sigma_lo^2 is a
    # source-faithful Loewner lower without intervalizing the scalar intensity.
    Qz = _scale(FROZEN._scaled_Q(x), _point(sigma2_lower))
    prediction_interval = FROZEN._sym(
        FROZEN._add(FROZEN._mul(FROZEN._mul(Fm, L0), Ft), Qz)
    )
    pred, pred_eps = _common_point_lower(prediction_interval)
    if not _strict_spd(pred):
        raise RuntimeError("Loewner prediction lower lost strict SPD; split x cell")

    R_aw, R_S_z = _physical_measurement_variances(node, h)
    post_aw, aw_eps = _point_measurement_lower(pred, 3, R_aw)
    post_s, s_eps = _point_measurement_lower(post_aw, 2, R_S_z)
    return post_s


def propagate_subcell(P, node: dict, samples: int, x: Interval, h: float):
    out = [[P[i][j] for j in range(4)] for i in range(4)]
    for _ in range(int(samples)):
        out = one_step(out, node, x, h)
    return out


def _split_request(exc: Exception) -> bool:
    """Both scaled-process producers ask for a narrower x cell by message.

    ``one_step`` raises ``RuntimeError`` when the collapsed Loewner lower is
    not strict SPD, and ``ou3_p3_frozen_full_matrix_translation._scaled_Q``
    raises ``ValueError`` when the cell straddles a scaled-process series
    branch.  Both are requests for subdivision, not proof failures.
    """
    text = str(exc)
    if isinstance(exc, RuntimeError):
        return "split x cell" in text
    if isinstance(exc, ValueError):
        return "split before evaluation" in text
    return False


def _split_x(x: Interval):
    """Subdivide exactly as the retained scaled-process cell splitter does.

    A cell that straddles a series branch is cut at that branch first, so one
    subdivision resolves it; otherwise the geometric midpoint is used.
    """
    for cut in (SCALED.BRANCH_X, SCALED.NEAR_EXACT_SERIES_MAX_X):
        if x.lo < cut < x.hi:
            return (
                Interval(x.lo, math.nextafter(cut, -math.inf)),
                Interval(cut, x.hi),
            )
    mid = math.sqrt(x.lo * x.hi)
    if not (x.lo < mid < x.hi):
        return None
    return Interval.outward_bounds(x.lo, mid), Interval.outward_bounds(mid, x.hi)


def _adaptive_image(P, node: dict, samples: int, x: Interval, h: float, depth: int):
    try:
        return [(x, propagate_subcell(P, node, samples, x, h))]
    except (RuntimeError, ValueError) as exc:
        if not _split_request(exc) or depth >= MAX_ADAPTIVE_X_DEPTH:
            raise
        halves = _split_x(x)
        if halves is None:
            raise RuntimeError(f"cannot further split failing x cell {x.as_list()}") from exc
        left, right = halves
        return (
            _adaptive_image(P, node, samples, left, h, depth + 1)
            + _adaptive_image(P, node, samples, right, h, depth + 1)
        )


def segment_images(P, source_node: int, gap: int,
                   runtime: dict | None = None,
                   domain_path: Path = DEFAULT_DOMAIN,
                   x_subcells: int = X_SUBCELLS):
    """Return rigorous point-lower images for one V1 applied-source segment."""
    rt = CORR.runtime(Path(domain_path).resolve()) if runtime is None else runtime
    s = int(source_node)
    g = int(gap)
    if not 0 <= s < len(rt["nodes"]):
        raise IndexError("source node outside P2 correlation partition")
    if g not in rt["gaps"]:
        raise ValueError("gap outside certified P2 correlation alphabet")
    node = rt["nodes"][s]
    h = float(rt["clock"]["dt_binary32_s"])
    rows = []
    for x in _node_x_subcells(node, h, x_subcells):
        for xx, posterior in _adaptive_image(P, node, g, x, h, 0):
            rows.append({
                "x_h_over_tau": xx.as_list(),
                "posterior": posterior,
                "posterior_is_common_Loewner_lower": True,
            })
    return rows


def build(domain_path: Path = DEFAULT_DOMAIN,
          representative_nodes=(0, 137, 729, 799),
          representative_gaps=(13, 21, 26)) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("correlated segment transfer must not be trajectory fitted")
    rt = CORR.runtime(path)
    corr = CORR.build(path)
    cf = CORR.validate(corr)
    if cf:
        raise RuntimeError(f"P2 correlation interface failed: {cf}")
    if corr.get("interface_version") != CORR.INTERFACE_VERSION:
        raise RuntimeError("P2 correlation interface version mismatch")

    zero = FROZEN._mat_zero()
    rows = []
    for s in map(int, representative_nodes):
        node = rt["nodes"][s]
        for g in map(int, representative_gaps):
            images = segment_images(zero, s, g, rt, path)
            rows.append({
                "source_node": s,
                "gap_samples": g,
                "tau_s": node["tau_s"],
                "sigma_filter_committed_mps2": node["sigma_filter_committed_mps2"],
                "R_S_filter_std": node["R_S_filter_std"],
                "same_source_cell_for_tau_sigma_R_S": True,
                "x_subcells": len(images),
                "all_outputs_interval_spd": all(
                    symmetric_positive_definite_ldlt(r["posterior"])[0]
                    for r in images
                ),
                "all_outputs_common_Loewner_lowers": all(
                    r.get("posterior_is_common_Loewner_lower") is True for r in images
                ),
            })

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_FULL_MATRIX_TRANSLATION_P2_CORRELATED_SEGMENT_TRANSFER",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "P2_correlation_interface_consumed": True,
        "P2_correlation_interface_version": CORR.INTERFACE_VERSION,
        "fixed_physical_scaling": "z=[v/h,p/h^2,S/h^3,a_w]",
        "sigma_dependent_state_rescaling_used": False,
        "tau_sigma_R_S_same_source_cell_per_segment": True,
        "full_4x4_translation_matrix_retained": True,
        "strongest_translation_measurements_every_sample": True,
        "fixed_measurement_constants_cached": True,
        "recursive_natural_interval_Riccati_subtraction_used": False,
        "deterministic_Loewner_lower_propagated": True,
        "prediction_interval_collapsed_by_validated_spectral_radius": True,
        "measurement_monotonicity_used": True,
        "sigma_process_intensity_lower_from_same_source_cell": True,
        "adaptive_x_subdivision_fail_closed": True,
        "gap_alphabet_samples": list(rt["gaps"]),
        "representative_rows": rows,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "next_obligation": (
            "use this stable Loewner-lower segment primitive in the all-source stage/phase P3 producer; "
            "if the unchanged 1e-18 gate is still limiting, tighten the exact limiting source/history/phase class rather than reintroducing interval posterior subtraction"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_FULL_MATRIX_TRANSLATION_P2_CORRELATED_SEGMENT_TRANSFER":
        f.append("wrong qualification")
    for key in (
        "source_only", "P2_correlation_interface_consumed",
        "tau_sigma_R_S_same_source_cell_per_segment",
        "full_4x4_translation_matrix_retained",
        "strongest_translation_measurements_every_sample",
        "fixed_measurement_constants_cached", "deterministic_Loewner_lower_propagated",
        "prediction_interval_collapsed_by_validated_spectral_radius",
        "measurement_monotonicity_used", "sigma_process_intensity_lower_from_same_source_cell",
        "adaptive_x_subdivision_fail_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "sigma_dependent_state_rescaling_used", "recursive_natural_interval_Riccati_subtraction_used",
        "P3_PROMOTED", "P4_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("P2_correlation_interface_version") != CORR.INTERFACE_VERSION:
        f.append("segment transfer is not bound to frozen P2 correlation version")
    if d.get("gap_alphabet_samples") != list(range(13, 27)):
        f.append("segment transfer lost 13..26 gap alphabet")
    rows = d.get("representative_rows", [])
    if not rows:
        f.append("no representative segment transfers emitted")
    for row in rows:
        if row.get("same_source_cell_for_tau_sigma_R_S") is not True:
            f.append("representative row lost source-cell correlation")
        if int(row.get("x_subcells", 0)) <= 0:
            f.append("representative row has no x subdivisions")
        if row.get("all_outputs_interval_spd") is not True:
            f.append("representative zero-start segment lost strict SPD lower")
        if row.get("all_outputs_common_Loewner_lowers") is not True:
            f.append("representative row did not emit common Loewner lowers")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True, default=lambda x: x.as_list()), encoding="utf-8")
    print(json.dumps({
        "P2_correlation_interface_version": d["P2_correlation_interface_version"],
        "representative_segments": len(d["representative_rows"]),
        "Loewner_lower": d["deterministic_Loewner_lower_propagated"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
