# OU-III P4 joint bias/source continuation

## Result

This continuation removes two relaxations that blocked the bounded-bias route:

1. accelerometer-bias error is no longer replaced by a fresh arbitrary vector at each sample; and
2. the physical bias driver is no longer rejected or separated from the source history.

The proof graph carries `(e_H,e_b,beta_true)` in 24 dimensions, while the
shipping filter itself remains the original 21-error estimator. The conditional
physical source used by the diagnostic is

    beta(t)=beta0 exp(-t/tau_true)+a sin(omega t+phase),

and every prediction obeys the exact serialized recurrence

    beta_i = phi_true beta_(i-1) + w_i,
    e_b,i = phi_hat e_b,(i-1)
            + (phi_true-phi_hat) beta_(i-1) + w_i.

The same `w_i` is therefore present once in physical truth and once in corrected
bias error, as required by `e_b=beta_true-b_hat`; it is not an independent ISS
input. Full shipping corrections, covariance cross-feedback, actual anisotropic
`R_S`, finite attitude reset, source forcing and projection remain attached to
the same history.

A source-indexed cyclic motion-storage sequence satisfies

    M_i = A_i^T M_(i+1) A_i + Q,
    M_N = M_0,

on the attached sample continuation. This is a compatible consecutive storage,
not one unrelated metric selected for every isolated word. Every completed
event also gets a direct physical-coordinate retention budget about the driven
center using one correlated covariance root.

## Driven point evidence

Dedicated CI run `34308574484` rebuilt both baseline and passive-trace shipping
observers with the same nonzero deterministic physical-bias driver. Passive
instrumentation remained bit-for-bit attached and the connected finite factors
passed in both H18 and A21.

The declared physical-bias model, its affine recurrence and the proof-only
24-D `beta_true` state reconstruct with zero reported defect. All 600
predictions have a nonzero driver increment.

| Quantity | H18 | A21 |
| --- | ---: | ---: |
| completed-event prefixes | 2591 | 2570 |
| one-sigma motion coordinates retained | yes | yes |
| limiting coordinate | velocity | latent acceleration |
| critical initial covariance level | 4.42138 sigma | 111.90580 sigma |
| worst one-sigma coordinate ratio | 0.31214 | 0.02872 |
| max sample compatible-storage ratio | 0.9999986022 | 0.9997983430 |
| max driver increment norm | 8.57e-7 | 7.92e-7 m/s^2 |

The canonical evidence is
`reports/results/rao_stability/joint-nonzero-bias-driver.json`.

These numbers remain a frozen source/coefficient point experiment. In
particular, the covariance ellipsoid is the filter's believed covariance, not
a qualified hard error set.

## Radial projection no longer needs a saturation replay

The driven center did not saturate: every recorded projection multiplier was
one. That does not force the nonlinear proof to enumerate or deliberately
manufacture a saturation trajectory.

Let `Pi_R` be Euclidean projection onto the estimate ball and define the error
projection map

    F_R(e,beta) = beta - Pi_R(beta-e).

For every two joint inputs,

    ||Delta F_R||^2 <= ||Delta e||^2 + ||Delta beta||^2.          (1)

Away from the projection boundary, let `J=D Pi_R`. Inside the ball `J=I`.
Outside, for `x=r u`,

    J=(R/r)(I-u u^T),

so `J` is symmetric and `0 <= J <= I`. Therefore

    dF = J de + (I-J) d beta,

and

    [J,I-J][J,I-J]^T
      = J^2 + (I-J)^2
      = I - 2 J(I-J) <= I.

At the boundary, each Clarke generalized Jacobian is a convex combination of
the inside and outside limits and still has spectrum in `[0,1]`; the same
inequality holds. Integrating the almost-everywhere derivative along a line
segment proves (1) globally. This covers unsaturated, saturated and boundary
crossings without freezing the point multiplier `s`.

For fixed physical truth, (1) reduces to ordinary nonexpansiveness in corrected
bias error. It also proves exact-real estimate-ball invariance immediately.
The executable contract and regression tests are in
`tools/stability/ou3_projection_sector.py` and
`tests/validation/test_ou3_projection_sector.py`; CI run `34308993332` passes.
Canonical output is `reports/results/rao_stability/projection-sector.json`.

This closes the **exact-real radial projection operator sector**, not the whole
P4 theorem. Floating-point rounding of the projection remains a separate
small enclosure.

## Remaining P4 closure

The dominant mathematical gap is now upstream of the radial clamp rather than
inside it. A complete P4 proof still needs a source-uniform enclosure that
turns the pointwise source-indexed storage into a consecutive inequality over
the admitted nonlinear coefficient family, together with a qualified word-entry
set and the declared physical/source premises. In particular:

- the Kalman/reset coefficient dependence on the retained state/source tube
  needs outward covering, including actual covariance and accepted-vector
  histories;
- the BIAS1 physical root/driver family needs a qualified hard bound rather
  than one deterministic point trajectory;
- the covariance entry ellipsoid needs a consistency or qualification theorem;
- canonical P3 execution/source admission remains separate and unchanged;
- projection floating-point rounding needs enclosure, although its exact-real
  nonlinear sector is now available as (1).

No deployed filter parameter, quality gate, physical domain or frozen P3
`delta=1e-18` is changed. `P4_MOTION_PASS=false`, `P4_PASS=false`, and neither
P5 path may start.
