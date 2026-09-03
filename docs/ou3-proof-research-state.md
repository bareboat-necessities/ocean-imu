# OU-III proof research state

Read this before modifying `tools/ou3_*.py` or the OU-III proof documents.
The working protocol is in `AGENTS.md`, section "OU-III Proof Research
Protocol". This file is the research state, not a changelog: it records what is
currently believed, what has been ruled out, and what experiment comes next.

## Immutable goal

For both fixed-dimension modes (H=18, A=21), on the declared operating domain,
establish or falsify

    sup_w rho_w < 1,     rho_w = V_after(F_w(x)) / V_before(x)

over the source-complete word family. Every construction below is a hypothesis
about how to reach that, never part of the goal.

## Current limiter

Canonical P3 cannot emit an artifact. `ou3_p3_p2_v1_stage_phase_translation`
aborts inside `common_boundary_floor`, which calls
`ou3_p3_correlated_translation_segment.segment_images`. That producer
propagates a covariance *lower* bound sample by sample and requires strict SPD
after every 5 ms step; on the small-x (large-tau) branch the intermediate lower
is genuinely near-singular in the S direction, and the adaptive x splitter
exhausts its depth without recovering.

Because P3 emits nothing, no real H/A margin exists, and therefore no P4
verdict is meaningful yet.

## Evidence

* The enclosure is **not** the binding constraint. Instrumented on node 137 /
  gap 13 with the congruence-scaled collapse: per-step `eps_D` is 5.3e-16 to
  1.9e-5 against correlation `lambda_min` 1.33e-3 to 1.48e-3 -- two to thirteen
  orders of margin, strict SPD true at every recorded step.
* Subdivision does not converge. The final traceback recurses down the left
  (small-x) branch to `MAX_ADAPTIVE_X_DEPTH` and still fails. Subdivision
  reduces interval dependency width; the near-singularity is structural, so the
  quantity that fails does not scale with cell width.
* `ou3_p3_frozen_full_matrix_translation._measurement_update` is in covariance
  form, `P - P e (e'Pe + R)^-1 e'P`. It requires only `P_cc + R > 0`, never
  strict SPD. Randomised checks: the update is monotone with a PSD lower
  (0/694 violations) and with an indefinite lower (0/4484). `F L F' + Q` stays
  PSD for singular `L` by congruence.
* That sibling producer already checks SPD only at `step in (1, steps)`.
  The segment producer is stricter than its own sibling, and nothing in the P3
  theorem requires per-step strict SPD -- P3 consumes the segment *endpoint*
  floor.

## Verified theorem steps

Randomised verification against the repo's own `symmetric_positive_definite_ldlt`.
All hold; none of these are suspect.

| Step | Result |
| --- | --- |
| `\|\|(R-I)v\|\|^2 = 4/(4+q^2) \|\|[c]x v\|\|^2` | exact, max err 4.7e-13 |
| `\|\|eta\|\| = sin(theta/2) \|\|h\|\|` | exact, max err 6.3e-13 |
| `\|\|R-I-[c]x\|\| <= (3/4) q^2` | 0 violations |
| `J <= n blockdiag(J_ii)` for PSD J | 0 violations |
| `K_theta S K_theta' <= P_theta_theta` | 0/1200 |
| `Phi_s' Sigma_s^-1 Phi_s <= Sigma_0^-1` | 0/1200 |
| `(1-3d/8)^2 <= 1-d/2` | holds for all d <= 16/9 |

Gap in the stated hypotheses: prefix nonexpansiveness needs `Phi_s` invertible.
It is, via the Joseph form and invertible OU prediction, but the document does
not say so.

## Available relaxations (not yet applied)

None of these change `rho = 1 - delta/2`; they only enlarge the certified
funnel `W_*`. Apply only when a real `delta` exists to work with.

* `lambda_max(Sigma) <= sum_g U_g` instead of `n_g max_g U_g` (verified, 0
  violations). With `U_S ~ 9e4` dominating: 9.04e4 vs 6.3e5, i.e. 7x.
* Attitude corrections need only the attitude marginal:
  `|dtheta| <= sqrt(U_theta / lambda_min(R)) |y|`, so sqrt(0.25) replaces
  sqrt(9e4) -- 600x on the correction gain.
* Sharp Cayley remainder `q^2(q+2)/(4+q^2)` is uniformly 0.805x the `(3/4)q^2`
  bound.
* The `4` in `B_m = 4 N_op sqrt(m_+) C / m_-` comes from a `W_s <= 4 W_0`
  bootstrap; prefix nonexpansiveness gives `W_s <= (1+d/8)^2 W_0 <= 1.27 W_0`,
  so 2 suffices.

## Dead ends

A rejected route may not be resurrected without stating which new mathematical
fact invalidates its rejection.

* **REJECTED -- recursive absolute Loewner point lower.** One absolute shift
  `eps = max_i sum_j r_ij` subtracted from every diagonal. The translation
  states span many orders (`v/h, p/h^2, S/h^3, a_w`), so the shift is set by
  the largest block and drives the smallest negative. Subdivision cannot fix it
  because the shift does not scale with cell width.
* **REJECTED -- congruence-scaled collapse as the fix for the same mechanism.**
  `A >= C - eps_D D^2` with `D = diag(sqrt(c_ii))` is correct and measurably
  good (see Evidence), and it did fix one of three failing tests. It did not
  fix the mechanism: 2 of 4 segment tests still error. Second failure of the
  same mechanism, so per the two-strike rule per-step-strict-SPD forward
  propagation is frozen. The code is retained because it is strictly better
  than the absolute form, but it is no longer load-bearing.
* **REJECTED -- `P_k >= (G_c^-1 + G_o)^-1`** as a segment endpoint floor
  (reachability and observability Gramians). Disproved numerically: 223/352
  violations on integrator-chain systems started from `P_0 = 0`. The
  information-form derivation needs a finite `Y_0`, which `P_0 = 0` does not
  provide. Recorded because it is superficially plausible and was nearly
  implemented.
* **REJECTED -- increasing `MAX_ADAPTIVE_X_DEPTH`.** The failing quantity does
  not scale with cell width, so no depth suffices.

## Retained

* Same-history P2-V1 source correlation.
* The exact Joseph information identity.
* The exact Cayley lift, residual identity and remainder sector at 0.8 rad.
* The covariance-form scalar measurement update and its monotonicity.
* Prefix nonexpansiveness in matched source information metrics.

## Known contract inconsistency

`tools/ou3_p4_canonical_gate.py` requires the candidate's `outer_angle_rad` to
equal 0.8 exactly, bound to the Cayley and remainder artifacts. The operating
domain declares a P4 certificate search over
`p4_complete_word_full_attitude_candidate_deg = [30, 25, 20, 15]`. If the
0.8-rad formulation is abandoned in favour of a narrower cell, the gate rejects
every candidate on that list. Resolve before relying on the search list.

## Out of order

`tools/ou3_p4_complete_word_dissipation.py` is a rigorous complete-word
producer written before any high-precision feasibility diagnostic of `rho_w`
existed. Under the protocol the diagnostic comes first. It is gated behind
canonical P3 in CI so it cannot produce a false PASS, but it must not be
treated as justified work until a diagnostic shows `rho_w < 1`.

## Next falsifiable experiment

Three qualitatively different alternatives to per-step-strict-SPD forward
propagation, to be compared before any is implemented:

1. **Endpoint-only SPD.** Keep forward propagation but drop the intermediate
   strict-SPD gates, carrying a PSD (possibly singular) lower and requiring
   strict SPD only at the segment endpoint. Justified by the covariance-form
   monotonicity above. Prediction: the split tree collapses and the run gets
   *faster*. If it does not, the diagnosis is wrong and this is abandoned, not
   deepened.
2. **Verified generalized-eigenvalue formulation.** Certify
   `P_endpoint >= rho_z D_h^2` directly as a generalized eigenvalue problem
   against the target floor, rather than propagating a matrix lower bound.
   Never forms an intermediate lower, so per-step conditioning cannot arise.
3. **Different Lyapunov representation.** Certify the endpoint floor in the
   `z = [v/h, p/h^2, S/h^3, a_w]` scaled coordinate with the floor built into
   the metric, so the quantity being propagated is dimensionless and the
   dynamic range that breaks the entrywise collapse is removed by construction.

Before P4 resumes: a non-promoting, structure-exact, high-precision `rho_w`
diagnostic, interpreted by the thresholds in `AGENTS.md`.
