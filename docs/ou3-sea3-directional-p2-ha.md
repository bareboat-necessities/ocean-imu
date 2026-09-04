# OU-III SEA3 directional response -> P2 -> H/A increment

This increment starts from `main` after PR #481 and implements the next proof
chain without redefining the frozen P2/P3 theorem interfaces.

## 1. Directional response enclosure

`tools/ou3_sea3_directional_response_domain.json` declares a provisional
three-mode directional response hypothesis.  For a complex three-axis response
vector `h(omega,theta)` it uses

```
a_hat(omega,theta) = omega^2 h(omega,theta) eta_hat(omega,theta)
```

and retains `h h*` as a positive-semidefinite matrix before taking the outer
norm bound.  The finite response band is tied by source parity to the deployed
sigma-band guard, 0.01--6 Hz.  With `||h||_2 <= G` and
`sum_r H_r^2 = H_s^2`, the enclosure records

```
tr M_a <= G^2 omega_hi^4 H_s^2 / 16.
```

The current `G=4` value is deliberately marked as a provisional feasibility
hypothesis, not a measured hull RAO.  It does not promote physical SEA0 and it
is forbidden from pruning P2.  A repository-owned vessel/IMU response model or
an explicitly accepted deployment response bound is still required for that
promotion.

## 2. SEA3-to-P2 inclusion

The inclusion result is independent of the provisional response gain.  The
shipping path already clamps the period-derived tuning frequency and the
resulting `tau`, `sigma_aw`, and `R_S` targets before the same samplewise EMA,
staging, one-sample pending apply, and finite 13--26-sample stage clock carried
by the existing P2 source graph.

Therefore every finite valid normal-Live SEA3 source realization projects into
the current broad P2 language.  The producer certifies

```
Lhat_SEA3 subset L_current_source
```

as a **non-pruning inclusion**.  It intentionally does not claim that any of
the 800 current P2 tuner cells are unreachable.  The finite response-weighted
WavePeriodEstimator/log-period/variance construction is still needed to obtain
that useful narrowing.

## 3. H/A feasibility

The second CI job consumes the inclusion artifact and builds the existing
source-complete P2-V1 translation candidate, the 18-state H / 21-state A
precision-block join, and the unique canonical P3 verdict.

The logic is one-way and fail-closed:

* if canonical P3 passes uniformly on the full P2 language, that result is
  inherited by the SEA3 subset;
* if canonical P3 fails on the broader P2 language, the SEA3 result is
  **inconclusive**, not failed, because a later SEA3-pruned source language may
  have a larger margin.

No P3 or P4 promotion authority is added by the SEA3 bridge.

## CI artifacts

`.github/workflows/ou3-sea3-directional.yml` produces:

* `ou3-sea3-directional-p2-inclusion` -- provisional matrix response enclosure
  plus the mechanical non-pruning SEA3-to-P2 certificate;
* `ou3-sea3-directional-ha-feasibility` -- the frozen full-P2 H/A candidate,
  canonical P3 verdict, and the SEA3 inheritance/inconclusive interpretation.

## Next proof step

Replace or explicitly accept the provisional response bound, propagate the
response-weighted finite-window moments through the actual WavePeriodEstimator,
log-period state, variance estimator, target libm/clamps, and staging state, and
use that construction to prune P2 source histories.  Then rerun the unchanged
H/A theorem interface on the narrower SEA3 language and report the limiting sea
parameters, direction/phase state, error direction, and distance to the
`rho=1` boundary.
