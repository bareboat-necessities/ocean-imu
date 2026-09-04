# SEA3 acceleration covariance coupling

The continuum RAO family retained after #482 has

`||h(f,theta)|| <= G min(1,(fc/f)^p),  p >= 2`.

That structure is enough to obtain a useful acceleration covariance certificate without adding another floating JONSWAP quadrature.

For every frequency and heading,

`(2*pi*f)^4 ||h(f,theta)||^2 <= G^2 (2*pi*fc)^4`.

The retained SEA3 spectral contract uses `m0_r = H_r^2/16` and `H_s^2 = sum_r H_r^2`, so after summing up to three partitions,

`tr Cov[a] <= (H_s^2/16) G^2 (2*pi*fc)^4 = pi^4 H_s^2 G^2 fc^4`.

`tools/ou3_sea3_acceleration_covariance_coupling.py` evaluates this inequality with outward interval arithmetic and compares it directly, without a square root in the proof decision, against the finite-horizon Gaussian concentration threshold.

## What it tells us

The old independent Cartesian sea x RAO box is still correctly rejected: `H_s=8.5 m`, `G=4`, `fc=1.2 Hz` is vastly outside the finite-horizon covariance threshold.

But the coupled domain is not empty. For a 20-minute, 200 Hz horizon (`N=240000`) and the current equal split of the declared 5% stochastic budget, the acceleration trace-covariance threshold is about `0.1481 (m/s^2)^2`. The same `H_s=8.5 m`, `G=4` tuple at the admitted low corner `fc=0.03 Hz` has the conservative shape-independent bound about `0.09121 (m/s^2)^2`, so that tuple passes the acceleration covariance condition.

This is important because it converts the vague “couple sea and RAO somehow” requirement into an explicit monotone predicate on `(H_s,G,fc)`. A real vessel can be qualified against that predicate, or against a tighter response-weighted spectral certificate if the uniform envelope is too conservative.

## Promotion boundary

This artifact still does not promote the SEA0 left inclusion. It proves a sufficient covariance condition and exhibits a nonempty passing part of the declared continuum response box, but it does not establish that the actual vessel/sea population lies in that part. Body-rate qualification also remains a separate P1 obligation because the current SEA3 RAO contract treats rotational response separately.

Once an actual vessel-response qualification is supplied, the acceleration covariance result can compose with the finite-horizon concentration kernel to produce the acceleration part of the Normal-Live good-event probability. Deterministic infinite-horizon containment of a non-degenerate Gaussian sea under a finite hard cap remains intentionally unclaimed.
