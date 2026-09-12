# BRMM stability proof handover

This is the canonical handover for the non-ALT (BRMM) route. Read `AGENTS.md`,
this file, and `docs/ou3-proof-research-state.md` before changing proof code.
Start the next continuation from the latest `main` and open a new PR.

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

## Named magnetometer theorem classes

The route carries two named magnetic classes with their own modules and no ALT
dependency. Keep the value class and the call-cadence class separate.

`MAG-BMM150-DET-v1` (`ou3_brmm_magnetic_source_envelope.py`) admits a
commissioned installation with `20 <= ||B_W|| <= 75 uT`, horizontal field
`>= 15 uT`, `||b_HI|| <= 5 uT` and `||n_m|| <= 2 uT` per theorem sample, for the
exact identity `m_body = R_true B_W + b_HI + n_m`. It is an engineering source
requirement, not a datasheet guarantee. `Rmag` is a model covariance and is
never read as a deterministic bound.

`MAG-CALL-SCHEDULE-v1` (`ou3_brmm_magnetic_call_schedule.py`) requires the first
post-Live `updateMag()` call within 0.04 s and every later gap at most 0.04 s.
That is 25 Hz, deliberately weaker than the sensor.

What they close, at unchanged gates and with no filter change:

- the declared magnetic PE band becomes a consequence. `||m_body||` lies in
  `[13, 82] uT`, so `normal_live.magnetic_vector_norm_lower_uT = 10` and upper
  `200` follow and the shipping `mag_init_min_mag_norm` guard clears
  unconditionally;
- the deterministic startup yaw capture. At a supplied tilt-frame error
  `E = ||b_HI|| + ||n_m|| + 2 Bmax sin(delta/2)` and
  `|sin(theta_yaw)| <= E/H_min`. At the declared 0.02 rad startup direction
  error this is `E <= 8.5 uT`, `|sin(theta_yaw)| <= 17/30`, `theta_yaw < 0.61`
  rad and a full attitude error `< 0.63` rad `= 36.0962 deg`, inside the
  declared 45 deg `initial_filter_entrance.attitude` set. No `1/sqrt(N)`
  statistical reduction is used anywhere on this path;
- the finite `H18 -> A21` release *timing*. The shipping counter increments on
  every post-delay `updateMag()` call regardless of innovation acceptance, so
  250 counts land within `10.0 s` of Live and the elapsed `9.96 s` strictly
  exceeds the 1 s guard.

What they do not close:

- universal per-sample admission at `max_sample_norm_ratio_from_mean = 0.35`
  needs combined hard iron + residual `<= 7 Bmin/47 = 2.9787 uT` against the
  declared `7 uT`;
- the vector sine separation stays a declared PE hypothesis; it is not
  derivable from this class;
- eventual A21 under an arbitrary external `acc_bias_hold_`, and the
  `H18 -> A21` joint24 covariance transport itself.

## Startup tilt supply is the limiting quantity

The shipping accumulation frame is the private observer's own tilt
(`tiltOnlyQuatFromBoatQuat_(attitudeReferenceQuat_())`) and the handoff seed is
`boatQuatWithAbsoluteYaw_(q_proxy, pending_yaw_abs_rad_)`, so the proxy tilt
drives both the gauge error and the seed tilt. Derived requirements, exact over
the rationals:

- `min_horizontal_fraction = 0.05` binds at `delta <= 0.0495289` rad
  `= 2.8378 deg`;
- non-vanishing north binds at `delta < 0.1067173` rad `= 6.1145 deg`;
- the certified proxy tilt is `86.2567 deg`: shortfall `30.3957x` and
  `14.1094x`, and `1.9168x` against the declared 45 deg entrance on its own.

The declared `world_averaged_gravity_direction_error = 0.02 rad` satisfies the
binding requirement but is the low-passed world-gravity direction error, not the
private observer's tilt error. Do not identify the two.

The private Mahony level-set route is frozen. At the padded `8.8 m/s^2`
envelope the seed angle is `asin(8.8/9.80665) = 1.1137 rad`, the level needed to
contain it is `1.8540056`, the largest level inside the 87 deg chart is
`1.4912551`, and the admissible level window is empty with factor `1.2432518`.
The boundary flow still closes with margin `0.0357940`, so the obstruction is the
level-set formulation, not the enclosure. A metric-free floor of
`m/s_min + 0.1*xi = 0.17603 rad = 10.0855 deg` — `3.5540x` the binding
requirement — rules out any further metric or subdivision work. The three
permitted alternatives are recorded in `docs/ou3-proof-research-state.md`.

## Fail-closed state and next work

`P4_PASS=false` and `P5_MAY_START=false`. Do not infer theorem closure from the
H18 information lemma, the magnetic certificates, or green CI.

Immediate open items for the next PR are:

1. Pick one of the three recorded startup-tilt alternatives and requalify the
   private Mahony/proxy startup for the padded `||a_wave|| <= 8.8 m/s^2` family,
   preserving both measured-period takeover and prior-frequency timeout paths.
   Do not reopen the metric route.
2. Discharge or replace the private-observer accumulation tilt-frame supply, and
   close the ungauged 150 s timeout branch: `ready_by_timeout` does not require
   `north_ready`, so that handoff still takes
   `proxy_handoff_yaw_sigma_free_rad` with no yaw gauge at all.
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
   capture and the H18->A21 joint24 covariance transport; only its release
   timing is closed.
8. Close deployment finite-precision enclosure and compose the indefinite
   end-to-end theorem.

Classify future failure as A/B/C/D/E/F/G using the established taxonomy. The old
constant-position witness may not be reused as A or B because it is outside the
corrected physical COMPLETE-BRMM source.
