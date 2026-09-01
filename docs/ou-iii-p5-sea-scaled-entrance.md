# OU-III P5 entrance and reduced P4 search domain

This note records the deployment entrance assumptions used to reduce the remaining OU-III P4/P5 certificate search without changing the shipping filter.

## Declared P5 entrance

The P5 capture theorem is presented with the following initial error set:

- on a gauged branch, the full attitude error is an **SO(3) geodesic error of at most 45 deg**;
- this is not an independent `+/-45 deg` roll/pitch/yaw component box;
- on an ungauged timeout branch, the 45 deg condition is used only for the gravity-direction/tilt coordinate until magnetic regauging; yaw remains a gauge;
- position is a truth-error condition, not an estimate-magnitude condition:
  `|delta p_x|, |delta p_y|, |delta p_z| <= 0.5 Hs`;
- `Hs > 0` is required.

For 45 deg the Cayley boundary is exact:

`q_45 = 2 (sqrt(2) - 1) ~= 0.82842712474619`.

The existing broad operation-matched geometry envelope is `0.80 rad = 45.84 deg`, with `q ~= 0.8456`, so the declared 45 deg entrance is immediately inside that geometry envelope. This does **not** require the complete 18/21-state P4 contraction proof to cover the whole 0.80 rad sector.

The componentwise position assumption implies

`||delta p||_2 <= (sqrt(3)/2) Hs ~= 0.8660254 Hs`.

Compared with the conservative P1 handoff bound `||delta p|| <= 20 m`, the new P5 entrance position norm is smaller whenever `Hs` is below about `23.1 m`. The old P1 box is deliberately retained: this PR does not claim that a preceding startup interval automatically preserves the `0.5 Hs` entrance condition.

## Sea-scaled translation coordinates

The complete-word P4/P5 search should use dimensionless translation coordinates

- `delta p / Hs`;
- `delta v Ts / Hs`;
- `delta S / (Hs Ts)`;
- `delta a_w Ts^2 / Hs`.

Only position receives a new hard entrance bound in this change. No hard sea-scaled limits on `v`, `S`, or `a_w` are invented. Their normalized forms are conditioning/partition coordinates; future proof work must intersect them with source-faithful physical bounds or derive additional bounds explicitly.

## Narrower P4 complete-word search

The broad `0.80 rad` finite-angle geometry certificate remains useful as an outer safety envelope. The expensive complete 18/21-state word does not need to start there. The search ladder is

`30 deg -> 25 deg -> 20 deg -> 15 deg`.

The rule is to attempt the widest candidate first and promote only the widest candidate for which source-complete, outward-rounded full-state dissipation closes. Decreasing the candidate angle improves the exact Cayley vector monotonicity and reduces the nonlinear eta/information ratio, so the ladder gives a controlled way to trade P4 search difficulty against the later P5 finite-capture distance.

This ladder is a certificate-search strategy, not an extra deployment theorem assumption.

### Measured limit of the descending ladder

`tools/ou3_p4_first_accel_sector_budget.py` evaluates the ladder against the
first deployed accelerometer operation.  Narrowing the candidate does widen the
correction budget and does lower the nuisance-over-budget ratio -- from `5.05`
at 30 deg to `2.12` at 15 deg -- but the producer's non-candidate limit probes
show the ratio saturating at `1.34` for a 1 deg candidate.  The floor is set by
the declared `0.3 g` latent-acceleration error over the lowest admitted specific
force, `2.941995 / 5.0 = 0.5884`, which no candidate angle changes.

The ladder therefore stays a useful conditioning strategy for the complete-word
search, but it is not the route that closes the first accelerometer operation.
That still needs the operation-matched information decrease and a directional
block margin.  The producer reports a distance, never a verdict.

## Proof obligations after this change

1. If the theorem starts before P5, propagate the declared 45 deg / `0.5 Hs` entrance set through the exact preceding startup/source interval; do not substitute it directly for the P1 handoff box.
2. Run complete H/A word propagation on the sea-scaled translation coordinates for the 30 deg P4 candidate first, then move inward only if a validated full-state cross-block/dissipation inequality fails.
3. Once one P4 candidate closes, derive a finite P5 word count from the 45 deg entrance to that certified P4 set.
4. Continue from that set to the existing inner stochastic localization level.

No filter gains, covariance equations, schedules, source branches, or estimator state dimensions are changed by this entrance/search reduction.
