# BRMM stability proof handover

This is the canonical handover for the continuation after PR #516. Read
`AGENTS.md`, this file, and `docs/ou3-proof-research-state.md` before changing
proof code. Start the next continuation from the latest `main` after #516 and
open a new PR.

## Shipping/runtime invariants

The shipping filter was not changed to accommodate the old constant-position
witness. Preserve the actual Mahony/proxy startup, one-time Live S-origin,
BIAS0/BIAS1/BIAS2 ancestry, same-signal WPE -> sigma -> tau -> T_S -> R_S
adaptation, staged commit/scheduler state, full P/H/R/K Joseph/reset/projection
path, and joint24 error/true-bias storage. P3 remains `delta = 1e-18`. Do not
resurrect the failed A21 18-state marginal-motion storage. Do not use replay,
finite seed sets, independent coefficient rectangles, packetwise S resets,
Gaussian good-event assumptions, or covariance consistency as hard state
membership.

## Corrected COMPLETE-BRMM physics

BRMM describes bounded oscillatory wave-induced vessel displacement relative to
a local equilibrium. Current, propulsion, leeway, secular/global translation
and arbitrary global-position origin offsets are outside the BRMM wave
coordinate. The primary physical pathwise property is a bounded centered
primitive,

`|| integral_u^t p_wave(s) ds || <= D_S`,

for every admitted continuation. Harmonic/spectral and bounded shaping-state
realizations are sufficient certificate constructions for this physical
condition; they do not define the physics.

The old exact history `p=d != 0, v=0, a=0` remains a regression for the old
finite-window-only source. Historical classification: obstruction B under the
old source surface, source-specification omission E under the intended physical
semantics, excluded after the corrected physical theorem. It is not a nominal
filter instability.

The current padded deterministic COMPLETE-BRMM engineering envelope is anchored
to the existing Hs=8.5 m reference family with 10% amplitude/kinematic padding
and 10% outward frequency-support padding:

- `Hs <= 9.35 m`
- `||p_wave|| <= 8.10 m`
- `||v_wave|| <= 5.50 m/s`
- `||a_wave|| <= 8.80 m/s^2`
- `||omega_body|| <= 35 deg/s`
- `f in [0.018, 0.88] Hz`
- `D_S <= 1100 m*s`

The finite-harmonic 28-ft reference gives `D_S <= 863.7794 m*s`; the padded
bound is `<1056 m*s` and is rounded outward to 1100. The legacy 300 m*s number
is not a physical source assumption. Keep physical source bounds distinct from
estimator error-entry/working radii.

The one-time Live transformation remains
`S_L(t)=S(t)-S(t_L)`, so `S_L(t_L)=0`. Do not wordwise re-zero S and do not
re-anchor position.

## Correlated innovation repair

The canonical measurement proof now carries the same `(P,H,R)` object through

`PHt = P H^T`,
`S = H P H^T + R`,
`S^-1`,
`K = P H^T S^-1`,
then Joseph/reset/projection.

`S` is no longer an independently selectable entrywise rectangular matrix in
this path. The structural Riccati facts `P >= 0` and `R > 0` exclude the false
singular innovation matrices that existed only in the rectangular hull.

## H18 information result closed at frozen gate

The former post-correlation H18 lower was
`7.092471820569811e-19`, only 0.7092471821 of the frozen 1e-18 gate. The limiting
term was the `eta6/a_w` cross-information enclosure.

Two proof-only tightenings close that numerical obstruction:

1. The four guaranteed S=0 witnesses are selected from
   `[0,g]`, `[4g,5g]`, `[8g,9g]`, `[12g,13g]`. The last selected window ends by
   `1.9499999564 s`, inside the same canonical 3 s word. Every other due S
   update remains in the literal shipping word, and the longer-horizon OU and
   process-noise penalties are retained.
2. The second PE accelerometer occurrence retains the same-history homogeneous
   OU attenuation of the initial `a_w` coordinate instead of assigning unit
   sensitivity twice.

Canonical CI now certifies

`lambda_min(H18 information) >= 4.253919518541475e-18`,

with gate ratio `4.253919518541474`. The coupled `eta6/a_w` lower is
`4.253919518541476e-18`; the cross-norm-squared upper is
`3767421.6507477993`. Directional translation information lowers are:

- `S`: `1.8434197928164813e-11`
- `g*p`: `2.3086404646316095e-11`
- `g^2*v`: `4.447787288026003e-10`
- `g^3*a_w`: `3.615365438779952e-08`

Thus the eta6/a_w H18 information bottleneck is closed without changing the
filter, source family, or 1e-18 gate.

## Fail-closed state and next work

`P4_PASS=false` and `P5_MAY_START=false`. Do not infer theorem closure from the
H18 information lemma or green CI.

Immediate open items for the next PR are:

1. Reconcile the current metric-memory/PE domain consistency failure:
   `declared PE does not refine vector certificate` in
   `ou3_brmm_riccati_tube.py::_declared_vector_alpha6`. Treat it as a proof
   representation/domain issue until evidence shows otherwise.
2. Requalify the actual private Mahony/proxy startup invariant for the padded
   `||a_wave|| <= 8.8 m/s^2` family and preserve both measured-period takeover
   and prior-frequency timeout paths.
3. Complete same-history physical source attachment through frontend/tuner and
   actual P/H/R/K lineage for every H18/A21 word.
4. Close downstream H18/A21 prior-free and finite-bias composition for
   BIAS0/BIAS1/BIAS2.
5. Close compatible consecutive joint24 storage with uniform coercivity,
   endpoint `W_next <= rho W + C` with `rho < 1`, and no hidden metric jump.
6. Close literal every-prefix augmented LDLT, including prediction, covariance
   floor, S Joseph, accelerometer, magnetometer, reset, projection,
   rejected/not-due and tuner/scheduler transitions.
7. Prove first-exit/domain retention over a physics-compatible working tube,
   determine the largest rigorously certifiable P4 basin, then prove finite H18
   capture and H18->A21 transport.
8. Close deployment finite-precision enclosure and compose the indefinite
   end-to-end theorem.

Classify future failure as A/B/C/D/E/F/G using the established taxonomy. The old
constant-position witness may not be reused as A or B because it is outside the
corrected physical COMPLETE-BRMM source.
