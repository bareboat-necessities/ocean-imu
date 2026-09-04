# SEA3 finite-horizon Gaussian concentration certificate

The hard Normal-Live P1 caps and an unbounded Gaussian sea model cannot be joined by a deterministic statement about **all** Gaussian realizations: any non-degenerate Gaussian sample has unbounded support. The existing operating-domain contract already declares a finite-horizon stochastic failure-probability budget of `0.05`, so the physically meaningful bridge is a finite-horizon good-event statement, not an infinite-horizon Gaussian hard-bound claim.

`tools/ou3_sea3_finite_horizon_concentration.py` adds a replay-free concentration kernel for that bridge.

## Bound

For a centered `d`-dimensional Gaussian response sample `X` with

`tr Cov[X] <= v`,

each component variance is at most `v`. If `||X||_2 > A`, at least one component exceeds `A/sqrt(d)` in magnitude. The one-dimensional Gaussian Chernoff bound plus a coordinate union bound gives

`P(||X||_2 > A) <= 2 d exp(-A^2/(2 d v))`.

A union bound across `N` sampled instants therefore gives

`P(max_k ||X_k||_2 > A) <= 2 d N exp(-A^2/(2 d v))`.

No independence between time samples is required, and arbitrary cross-axis covariance is allowed.

The producer chooses the smallest integer `t` for which a validated upper enclosure of

`2 d N exp(-t)`

fits inside the allocated event budget. It then emits the sufficient covariance-trace condition

`v <= A^2/(2 d t)`.

The proof decision does not call ordinary `libm` exponentials: `exp(-t)` is built from repeated outward interval products of the repository's validated `exp(-1/2)` primitive. The square-root RMS numbers in the artifact are diagnostic only; the PASS threshold is the covariance-trace inequality above.

## Budget allocation

The current kernel splits the declared `0.05` finite-horizon budget equally between the non-gravitational CoG acceleration event and the body-rate event. This is conservative and can later be optimized without changing the theorem.

For orientation, the integer construction gives approximately these diagnostic thresholds:

| Horizon | Samples at 200 Hz | integer `t` | acceleration trace-covariance cap | acceleration trace RMS |
| --- | ---: | ---: | ---: | ---: |
| 1 s | 200 | 11 | about `0.2424 (m/s^2)^2` | about `0.492 m/s^2` |
| 20 min | 240000 | 18 | about `0.1481 (m/s^2)^2` | about `0.385 m/s^2` |

The exact machine artifact uses outward-rounded values.

## What this closes and what it does not

This closes a reusable **finite-horizon concentration kernel**. It does not yet close the SEA0 left inclusion because no validated coupled JONSWAP/RAO producer is attached for the acceleration and body-rate covariance traces.

The next numerical proof task is therefore much sharper than before: derive uniform validated trace-covariance upper bounds for the coupled sea/RAO family and compare them with these thresholds. If the current robust RAO envelope is too broad, the coupling must be expressed physically rather than by reverting to an invalid independent sea x RAO Cartesian product.

Even after such a candidate passes, the theorem statement must remain finite-horizon/high-probability (or explicitly conditioned on the good event). This certificate deliberately does not claim deterministic or infinite-horizon containment of an unbounded Gaussian process inside finite hard P1 caps.
