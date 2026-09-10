# OU-III continuation note — end-to-end shipping stability

Canonical handoff. Read `docs/ou3-end-to-end-stability-theorem.md` first; this
file is the short list of what to do next and what not to redo.

## The target changed

The goal is no longer "close P4". It is the end-to-end theorem

```
declared startup uncertainty
  --[P5 capture, finite time, real Mahony/proxy path]-->  R_P4
  --[P4 invariance]-->  practical ISS forever
```

for every admitted BRMM motion, every admitted BIAS0/BIAS1/BIAS2 history and
every admitted sensor disturbance. `ou3_end_to_end_stability_gate.py` holds the
composition fail-closed. `END_TO_END_STABILITY_PASS = false`.

## The two rules that decide most design choices

1. **Maximise the basin, never shrink it.** Width certified in P4 is width the
   capture stage does not have to achieve. A tiny entrance set does not make the
   theorem easier, it makes P5 harder and can demand unrealistic settling times.
   Use `ou3_p4_basin_frontier.py`: the chart hypothesis is a half-space in
   scaled-radius space, so the basin is a correlated polytope and there is a
   whole frontier to choose from, not one number.
2. **Never assert unreachability, prove it.** Removing coordinate combinations
   the implementation cannot generate is legitimate and is NOT domain shrinking,
   but only with a reachability argument. `ou3_p4_live_entry_reachability.py`
   is the template: it reads the shipping sources and checks parity strings.

## What is already settled — do not redo it

* **The Live-entry state is pinned by the deployed code.** `updateFrontEnd`
  passes `drive_mekf=false` and every MEKF drive call is inside that guard, so
  the translational block, the latent state and both bias states are never
  propagated before `goLive` and stay at the constructor zero. `goLive` zeroes
  the whole attitude-to-linear cross covariance. Hence at entry every
  non-attitude coordinate is the negated physical truth of one BRMM history at
  one instant, `P_theta,lin = 0` and `K_theta,S = 0`.
* Therefore the 300 m*s integral ball is not the object of study at all. The
  reachable integral entry radius is the BRMM primitive `S_m`.
* **All three bias families are closed** on the same-history graph, and `mu_sep`
  is an optional sharpener, not a blocker.
* The conditional binary32 additive ISS enclosure is closed. Target-toolchain
  qualification is a separate deployment blocker and is not a mathematical one.
* The adaptive coefficients are pinned: `T_S = clamp(pseudo_ratio*tau)`,
  `R_S_target = SpectralMSE(tau,sigma)`, applied `R_S` is the ray
  `(.72r,.72r,r)`, and the clamped target box is forward invariant for the
  candidate and active schedules given initial membership, which the producer
  checks (`candidate_ema` does not clamp, so the premise is not free).
* No class A counterexample and no rigorous class B infeasibility exists. The
  4.5788 frozen-word Cayley value is not a lower bound on the true nonlinear
  trajectory and must not be quoted as a disproof.

## Next moves, in order

1. **Instantiate `V_m`, `P_m`, `S_m`** in the physical BRMM contract from the
   declared source class. They are the entry radii now, and they are `None`.
   The certificate currently admits `S_m <= 7.3537 m*s` from the same-cell
   correction ceiling and `4.566 m*s` from the chart frontier with the other
   coordinates at their declared radii, rising to `11.03 m*s` at the
   maximum-volume basin point. If the qualified `S_m` lands above those, the
   magnitude-only route has no headroom and the proof must retain the same-cell
   attitude/S cross-covariance direction instead of the Cauchy-Schwarz product.
2. **Repair the attitude covariance envelope.** The endpoint-referenced
   Lagrange/Vandermonde inversion gives 3.99983e8 rad^2 per axis, 6.4485e11
   above the prior-independent accelerometer posterior cap `r/|f|^2`. That one
   number is what blocks the correction/reset domain, not the entry set. Carry
   the per-sample accelerometer cap into the envelope and keep the same-cell
   `S^{-1}` instead of relaxing to `R^{-1}`.
   The third direction, rotation about the specific force, has a route rather
   than a wall: the same prior-independent argument applies to the magnetometer
   event in the two directions transverse to `m`, and the declared
   `vector_sine_separation_lower = .1` means `f` and `m` are never parallel, so
   the two transverse pairs span all three directions. What is missing is a
   maximum magnetometer inter-event gap. The domain declares the magnetometer an
   asynchronous family with no cadence upper bound; if the
   `vector_pe_recurrence_window_s = 1.0` assumption is read as guaranteeing an
   accepted magnetic packet in every window, the gap is bounded and the yaw
   covariance obeys `cap_at_last_mag_event + growth over the gap`. Settle that
   reading first, because the growth term needs a gyro-bias covariance bound and
   the retained one has the same conditioning defect being repaired.
3. **Finish the source cover.** `ou3_p4_coefficient_dependency_cover.py` has the
   pinnings, the forward-invariant coefficient box and the deployed EMA rate
   limits, which are the adjacency data a consecutive-word cell family needs.
   What is missing is the reachable `(f,sigma)` set over admitted BRMM
   continuations and the composition with the Riccati/Joseph/reset graph.
4. **Then the capture half.** The deployed machine is
   `Cold -> TunerWarm -> TunerReady -> (external bootstrap goLive) -> Live`.
   Only the `Cold -> TunerWarm` elapsed-time warmup is discharged. The open
   ones are tuner/usable-period attainment, bootstrap tilt and north
   acquisition, the certified proxy attitude radius, basin membership at entry,
   and the H18-to-A21 release.

## Dead ends

Do not spend time on these unless new structure materially changes them.

* isolated spectral radius of a captured word;
* one independently selected Lyapunov metric per word without consecutive-word
  compatibility;
* detached per-sample bias-energy ports;
* marginal `Pbar` plus hard residual boxes for reset-domain closure;
* rowwise independent `K` boxes for contraction or reset;
* scalar nonlinear eta budgeting against the tiny global Riccati/P3 margin;
* independent `(f,sigma)` tuner rectangles as a theorem family, and independent
  `(tau,sigma,T_S,R_S)` rectangles generally: the deployed map pins
  `T_S = clamp(pseudo_ratio*tau)` and `R_S_target = SpectralMSE(tau,sigma)`;
* covariance ellipsoids used as if they were hard physical entry sets;
* replay or a pinned RAO history used as a source-uniform proof;
* anchoring `e_S` through the `S=0` regulation's contraction. It is a
  non-expansion in the `R_S`-weighted norm, but the source-uniform `P_SS` floor
  gives `mu >= 1.3e-26`, above 1e25 events per e-fold. The entry value comes
  from the deployed pinning instead.

## Suggested first prompt

> Continue OU-III from current `main` using
> `docs/ou3-end-to-end-stability-theorem.md` as the canonical handoff. The
> target is the end-to-end theorem, not P4 alone. Maximise the certified basin
> rather than shrinking it, instantiate the BRMM `V_m`/`P_m`/`S_m` primitives,
> repair the attitude covariance envelope with the prior-independent
> accelerometer posterior cap, then finish the source cover and the capture
> half. Keep every promotion bit false unless the full certificate closes.
