#!/usr/bin/env python3
"""What one vector Joseph update can and cannot bound in the attitude block.

This module exists because a bound asserted earlier in this proof line is FALSE
for the deployed filter, and the refutation matters more than the bound.

The claim was: a scalar angle measurement with gain ``|f|`` and noise variance
``r`` gives ``P^+ = P^- r/(P^- |f|^2 + r) <= r/|f|^2`` for every prior, so the
accelerometer Joseph update caps the reachable attitude variance transverse to
the specific force at ``R_acc/|f|^2`` prior-independently.

That is true only when attitude is the ONLY state in the residual.  The deployed
accelerometer residual is not of that form.  ``measurement_update_acc_only``
builds ``J_att = -skew(f_cog_b)`` alongside a latent-acceleration block and an
accelerometer-bias block, so the same residual is produced by a transverse
attitude error and by an ``a_w`` error.  One update cannot separate them, and the
attitude MARGINAL of the full-state posterior is not capped by anything.
``refutation()`` exhibits that directly: at a latent-acceleration prior of 1e8
the transverse attitude marginal stays above 1e6 rad^2, against a claimed cap of
4.16e-4.

What survives is a CONDITIONAL cap.  Writing the posterior through its
variational form

    c^T P^+ c = min_k [ (c - H^T k)^T P (c - H^T k) + k^T R k ],

take ``c = (e,0,0)`` with ``e`` a unit vector orthogonal to ``f`` and choose

    k = -(f x e)/|f|^2,   so   H_theta^T k = e  and  |k| = 1/|f|.

That kills the attitude component exactly, leaving ``c - H^T k = (0,-R^T k,-k)``
with ``R`` the rotation carried by the latent block, so ``|R^T k| = |k|``, and

    e^T P^+_theta,theta e  <=  (sigma_a^2 + 2 lambda_max(P_(a_w,b_a))) / |f|^2 .

The bound is conditional on the joint latent/bias covariance block, and it is
tight rather than merely true.  ``sweep()`` searches the full deployed nine-state
residual structure ``H = [-skew(f), R^T, I]`` with random specific force, random
carried rotation and prior condition numbers spanning many decades; it runs on
every build, finds no violation, and attains about 91% of the bound.  Cases where
the variational identity fails to reproduce ``c^T P^+ c`` to a relative 1e-8 are
rejected as numerically invalid rather than counted either way -- at these
condition numbers the posterior subtraction can lose the answer outright, and an
ill-conditioned case is evidence of nothing.  A larger out-of-repo sweep of
188609 validated cases attains 99.75%.

Consequence, and the reason this is worth recording rather than quietly
dropping: with the deployed constants the cap is of order 1 rad^2, not 1e-3.  It
is a large improvement on the retained endpoint-referenced envelope but it does
NOT bring the same-cell Joseph correction ceiling inside the reset utility
domain, and it does not show that envelope to be loose.  The one-shot
information route is therefore a dead end for the correction/reset blocker; a
bound on the attitude covariance has to come from the uniform
observability/detectability machinery instead, where the latent and bias states
are separated over a window rather than at one event.
"""
from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
QUALIFICATION = "OU3_P4_ATTITUDE_MEASUREMENT_CAP_V1"

# The in-repo sweep must attain at least this fraction of the bound.  A bound
# nothing approaches is not evidence that the bound is right, so this is a floor
# on tightness, not a target: if a future edit loosens the cap the sweep stops
# reaching the floor and validate() fails instead of quietly passing.
CONDITIONAL_BOUND_TIGHTNESS_FLOOR = 0.85
SWEEP_CASES = 3000


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _skew(v):
    return [[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]]


def _matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
            for i in range(len(A))]


def _transpose(A):
    return [list(col) for col in zip(*A)]


def _solve(A, b):
    """Gaussian elimination with partial pivoting, for the 3x3 innovation solve."""
    n = len(A)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-300:
            raise ZeroDivisionError("singular innovation covariance")
        M[col], M[piv] = M[piv], M[col]
        for r in range(n):
            if r == col:
                continue
            factor = M[r][col] / M[col][col]
            for c in range(col, n + 1):
                M[r][c] -= factor * M[col][c]
    return [M[i][n] / M[i][i] for i in range(n)]


def refutation(latent_prior: float = 1.0e8, attitude_prior: float = 1.0e8,
               f_norm: float = 9.80665, sigma_acc: float = 0.2) -> dict:
    """The prior-independent transverse cap fails once ``a_w`` shares the residual.

    Six states ``(theta, a_w)``, the deployed ``H_theta = -skew(f)`` and an
    identity latent block, one exact Joseph update.  The transverse attitude
    marginal is reported against the claimed cap ``sigma_a^2/|f|^2``.
    """
    f = [0.0, 0.0, -float(f_norm)]
    Hth = [[-x for x in row] for row in _skew(f)]
    H = [Hth[i] + [1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    P = [[0.0] * 6 for _ in range(6)]
    for i in range(3):
        P[i][i] = float(attitude_prior)
        P[3 + i][3 + i] = float(latent_prior)
    r = float(sigma_acc) ** 2

    HP = _matmul(H, P)
    S = _matmul(HP, _transpose(H))
    for i in range(3):
        S[i][i] += r
    # P^+ = P - (HP)^T S^{-1} (HP)
    cols = [_solve(S, [HP[i][j] for i in range(3)]) for j in range(6)]
    Pplus = [[P[i][j] - sum(HP[k][i] * cols[j][k] for k in range(3)) for j in range(6)]
             for i in range(6)]
    # transverse directions: e1, e2 are orthogonal to f = -|f| e3
    transverse = max(Pplus[0][0], Pplus[1][1])
    claimed = r / (float(f_norm) ** 2)
    # The covariance subtraction P - (HP)^T S^{-1} (HP) cancels against the prior,
    # so a verdict has to be taken outside the cancellation noise floor: the
    # absolute error is O(ulp(largest prior)).  Without this the degenerate
    # latent_prior=0 case, where the scalar identity is exactly right, reads as a
    # refutation at a relative 2e-5 and the refutation would be worthless.
    noise = 64.0 * math.ulp(max(float(attitude_prior), float(latent_prior), 1.0))
    return {
        "model": "theta and a_w share one accelerometer residual, H = [-skew(f), I]",
        "attitude_prior": float(attitude_prior),
        "latent_acceleration_prior": float(latent_prior),
        "claimed_prior_independent_cap": claimed,
        "achieved_transverse_attitude_marginal": transverse,
        "cancellation_noise_floor": noise,
        "claim_holds": bool(transverse <= claimed + noise),
        "exceedance_factor": up(transverse / down(claimed)),
        "refutation_is_outside_the_noise_floor": bool(transverse > claimed + noise),
    }


def _lambda_max(M, iters: int = 400) -> float:
    """Power iteration for the top eigenvalue of a symmetric PSD block."""
    n = len(M)
    v = [1.0 / math.sqrt(n)] * n
    lam = 0.0
    for _ in range(iters):
        w = [sum(M[i][j] * v[j] for j in range(n)) for i in range(n)]
        nw = math.sqrt(sum(x * x for x in w))
        if nw <= 0.0:
            return 0.0
        v = [x / nw for x in w]
        lam = nw
    # Rayleigh quotient is the accurate readout once v has converged.
    Mv = [sum(M[i][j] * v[j] for j in range(n)) for i in range(n)]
    return max(lam, sum(v[i] * Mv[i] for i in range(n)))


def sweep(cases: int = 3000, seed: int = 20260910) -> dict:
    """Adversarial search for a violation of the conditional cap.

    Builds the FULL deployed residual structure ``H = [-skew(f), R^T, I]`` over
    nine states with random specific force, random carried rotation and random
    priors spanning many decades of conditioning, forms the exact posterior, and
    compares the transverse attitude marginal against ``conditional_cap``.

    A random prior at condition number 1e12 can make the posterior subtraction
    lose the answer entirely, and such a case must not be reported either as a
    violation or as a confirmation.  Each case is therefore gated on the
    variational identity reproducing ``c^T P^+ c`` to a relative 1e-8; cases that
    fail the gate are rejected and counted, not silently kept.
    """
    rng = random.Random(seed)
    kept = rejected = violations = 0
    worst = 0.0
    worst_case = None
    for _ in range(int(cases)):
        fmag = rng.uniform(8.0, 12.0)
        fd = [rng.gauss(0.0, 1.0) for _ in range(3)]
        nf = math.sqrt(sum(x * x for x in fd))
        fd = [x / nf for x in fd]
        f = [fmag * x for x in fd]
        Rrot = _random_rotation(rng)
        Hth = [[-x for x in row] for row in _skew(f)]
        RT = _transpose(Rrot)
        H = [Hth[i] + RT[i] + [1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
        sig = rng.uniform(0.05, 0.5)
        r = sig * sig

        scale = [10.0 ** rng.uniform(-3.0, 4.0) for _ in range(9)]
        A = [[rng.gauss(0.0, 1.0) for _ in range(9)] for _ in range(9)]
        C = _matmul(A, _transpose(A))
        P = [[scale[i] * C[i][j] * scale[j] / 9.0 for j in range(9)] for i in range(9)]
        tr = sum(P[i][i] for i in range(9))
        for i in range(9):
            P[i][i] += 1.0e-12 * tr / 9.0

        # unit e orthogonal to f
        w = [rng.gauss(0.0, 1.0) for _ in range(3)]
        d = sum(w[i] * fd[i] for i in range(3))
        w = [w[i] - d * fd[i] for i in range(3)]
        nw = math.sqrt(sum(x * x for x in w))
        if nw < 1e-6:
            rejected += 1
            continue
        e = [x / nw for x in w]

        HP = _matmul(H, P)
        S = _matmul(HP, _transpose(H))
        for i in range(3):
            S[i][i] += r
        try:
            cols = [_solve(S, [HP[i][j] for i in range(3)]) for j in range(9)]
        except ZeroDivisionError:
            rejected += 1
            continue
        Pp = [[P[i][j] - sum(HP[k][i] * cols[j][k] for k in range(3)) for j in range(9)]
              for i in range(9)]
        achieved = sum(e[i] * Pp[i][j] * e[j] for i in range(3) for j in range(3))

        # validity gate via the variational identity at its own minimiser
        c = e + [0.0] * 6
        rhs = [sum(HP[i][j] * c[j] for j in range(9)) for i in range(3)]
        try:
            kstar = _solve(S, rhs)
        except ZeroDivisionError:
            rejected += 1
            continue
        dv = [c[j] - sum(H[i][j] * kstar[i] for i in range(3)) for j in range(9)]
        obj = (sum(dv[i] * P[i][j] * dv[j] for i in range(9) for j in range(9))
               + r * sum(x * x for x in kstar))
        if not (abs(obj - achieved) <= 1e-8 * max(abs(obj), abs(achieved), 1e-300)):
            rejected += 1
            continue

        kept += 1
        lam = _lambda_max([row[3:] for row in P[3:]])
        bound = conditional_cap(sig, fmag, lam)
        if achieved > bound:
            violations += 1
        ratio = achieved / bound
        if ratio > worst:
            worst = ratio
            worst_case = {"specific_force_norm": fmag, "sigma_accelerometer": sig,
                          "joint_latent_bias_lambda_max": lam,
                          "achieved_transverse_marginal": achieved, "bound": bound}
    return {
        "cases": int(cases), "seed": int(seed),
        "numerically_valid_cases_kept": kept,
        "ill_conditioned_cases_rejected": rejected,
        "violations": violations,
        "max_achieved_over_bound": worst,
        "worst_case": worst_case,
        "residual_structure": "H = [-skew(f), R^T, I] over (theta, a_w, b_a)",
    }


def _random_rotation(rng):
    """Gram-Schmidt on a random 3x3, then fix the sign to keep det = +1."""
    cols = []
    while len(cols) < 3:
        v = [rng.gauss(0.0, 1.0) for _ in range(3)]
        for u in cols:
            d = sum(v[i] * u[i] for i in range(3))
            v = [v[i] - d * u[i] for i in range(3)]
        n = math.sqrt(sum(x * x for x in v))
        if n < 1e-6:
            continue
        cols.append([x / n for x in v])
    M = _transpose(cols)
    det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
           - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
           + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    if det < 0.0:
        for i in range(3):
            M[i][2] = -M[i][2]
    return M


def conditional_cap(sigma_acc: float, f_norm_lower: float, latent_bias_lambda_max: float) -> float:
    """Outward ``(sigma_a^2 + 2 lambda_max(P_(a_w,b_a)))/|f|^2``.

    Valid for every prior and every unit direction orthogonal to the specific
    force, by the variational argument in the module docstring.
    """
    if not (sigma_acc > 0.0 and f_norm_lower > 0.0 and latent_bias_lambda_max >= 0.0):
        raise ValueError("positive sigma/|f| and nonnegative latent-bias bound required")
    numerator = up(up(sigma_acc * sigma_acc) + up(2.0 * latent_bias_lambda_max))
    return up(numerator / down(f_norm_lower * f_norm_lower))


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    sigma_acc = min(map(float, domain["configured_runtime"]["measurement_noise_std"]["accelerometer_mps2"]))
    f_lower = float(domain["normal_live"]["specific_force_norm_lower_mps2"])

    # The joint (a_w, b_a) block bound: for a PSD block matrix,
    # lambda_max([[A,B],[B^T,D]]) <= lambda_max(A) + lambda_max(D), because
    # |u^T B v| <= sqrt(u^T A u) sqrt(v^T D v) gives
    # x^T M x <= (sqrt(a)|u| + sqrt(d)|v|)^2 <= (a+d)(|u|^2+|v|^2).
    import ou3_brmm_riccati_tube_factored as TUBE
    tube = TUBE.build()
    tube_failures = TUBE.validate_covariance_ceiling(tube)
    if tube_failures:
        raise RuntimeError("covariance ceiling invalid: " + repr(tube_failures))
    diag_a = list(map(float, tube["modes"]["A"]["Pbar_diagonal_variance_upper"]))
    latent_lambda = max(diag_a[15:18])
    bias_lambda = max(diag_a[18:21])
    joint_lambda = up(latent_lambda + bias_lambda)

    cap = conditional_cap(sigma_acc, f_lower, joint_lambda)
    ref = refutation()
    sw = sweep(SWEEP_CASES)

    return {
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "deployed_attitude_jacobian": "J_att = -skew(f_cog_b) in measurement_update_acc_only",
        "residual_also_carries_latent_and_bias_blocks": True,

        "prior_independent_transverse_cap_claimed": down(sigma_acc * sigma_acc / up(f_lower * f_lower)),
        "prior_independent_transverse_cap_holds_for_the_deployed_filter": False,
        "refutation": ref,

        "conditional_cap_statement": (
            "for every prior and every unit e orthogonal to f, "
            "e^T P^+_theta,theta e <= (sigma_a^2 + 2 lambda_max(P_(a_w,b_a)))/|f|^2"),
        "conditional_cap_argument": (
            "variational posterior c^T P^+ c = min_k [(c-H^T k)^T P (c-H^T k) + k^T R k] "
            "evaluated at k = -(f x e)/|f|^2, which gives H_theta^T k = e and |k| = 1/|f|"),
        "conditional_cap_is_prior_independent_in_attitude_only": True,
        "joint_latent_bias_lambda_max_argument": (
            "lambda_max([[A,B],[B^T,D]]) <= lambda_max(A) + lambda_max(D) for PSD blocks"),
        "latent_lambda_max": latent_lambda,
        "bias_lambda_max": bias_lambda,
        "joint_latent_bias_lambda_max": joint_lambda,
        "sigma_accelerometer": sigma_acc,
        "specific_force_norm_lower": f_lower,
        "conditional_transverse_cap_rad2": cap,
        "sweep": sw,
        "conditional_bound_attained_fraction": sw["max_achieved_over_bound"],
        "conditional_bound_tightness_floor": CONDITIONAL_BOUND_TIGHTNESS_FLOOR,

        "cap_closes_the_correction_reset_domain": False,
        "one_shot_information_route_is_a_dead_end_for_the_correction_domain": True,
        "attitude_covariance_bound_needs_uniform_observability_over_a_window": True,
        "envelope_shown_to_be_loose_by_this_argument": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "bound the attitude covariance through the uniform observability/detectability "
            "machinery, where the latent-acceleration and bias states are separated over a "
            "window instead of at a single event"),
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("qualification") != QUALIFICATION:
        f.append("qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("source changed")
    for k in ("residual_also_carries_latent_and_bias_blocks",
              "conditional_cap_is_prior_independent_in_attitude_only",
              "one_shot_information_route_is_a_dead_end_for_the_correction_domain",
              "attitude_covariance_bound_needs_uniform_observability_over_a_window"):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in ("prior_independent_transverse_cap_holds_for_the_deployed_filter",
              "cap_closes_the_correction_reset_domain",
              "envelope_shown_to_be_loose_by_this_argument", "P4_promoted_here"):
        if d.get(k) is not False:
            f.append(k + " not false")
    # The refutation must actually refute: the achieved marginal has to exceed
    # the claimed cap, or the claim is not disproved and must not be reported as
    # disproved.
    ref = d.get("refutation", {})
    if ref.get("claim_holds") is not False:
        f.append("refutation does not refute the claimed cap")
    if not (float(ref.get("exceedance_factor", 0.0)) > 1.0):
        f.append("refutation exceedance factor is not above one")
    if ref.get("refutation_is_outside_the_noise_floor") is not True:
        f.append("refutation is inside the covariance-cancellation noise floor")
    if float(ref.get("achieved_transverse_attitude_marginal", -1.0)) <= \
            float(ref.get("claimed_prior_independent_cap", 0.0)):
        f.append("refutation numbers do not exceed the claim")
    cap = float(d.get("conditional_transverse_cap_rad2", -1.0))
    if not (math.isfinite(cap) and cap > 0.0):
        f.append("conditional cap is not a positive finite bound")
    # Reporting the cap as a closure would be the same error again.
    if cap <= 1.0e-3:
        f.append("conditional cap is suspiciously small; recheck before consuming it")
    # The sweep is the evidence that the conditional bound is both true and not
    # vacuous. Both halves are enforced: no violation, and the bound is actually
    # approached.
    sw = d.get("sweep", {})
    if int(sw.get("violations", -1)) != 0:
        f.append("sweep found a violation of the conditional cap")
    kept = int(sw.get("numerically_valid_cases_kept", 0))
    if kept < int(sw.get("cases", 0)) // 2:
        f.append("sweep rejected too many cases to be evidence of anything")
    frac = float(d.get("conditional_bound_attained_fraction", 0.0))
    if not (frac <= 1.0):
        f.append("sweep attained more than the bound; the bound is violated")
    if frac < CONDITIONAL_BOUND_TIGHTNESS_FLOOR:
        f.append("conditional cap is no longer tight; the sweep stopped approaching it")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "prior_independent_claim_holds": d["prior_independent_transverse_cap_holds_for_the_deployed_filter"],
        "claimed_cap": d["prior_independent_transverse_cap_claimed"],
        "refutation_achieved": d["refutation"]["achieved_transverse_attitude_marginal"],
        "refutation_exceedance": d["refutation"]["exceedance_factor"],
        "conditional_cap_rad2": d["conditional_transverse_cap_rad2"],
        "closes_correction_domain": d["cap_closes_the_correction_reset_domain"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
