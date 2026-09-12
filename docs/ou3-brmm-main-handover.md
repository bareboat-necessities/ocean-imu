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
- the startup yaw *algebra*. At a supplied accumulation-frame excursion
  `E = ||b_HI|| + ||n_m|| + 2 Bmax sin(delta/2)` and
  `|sin(theta_yaw)| <= E/H_min`. At `delta <= 0.02` rad this is `E <= 8.5 uT`,
  `|sin(theta_yaw)| <= 17/30`, `theta_yaw < 0.61` rad and a full attitude error
  `< 0.63` rad `= 36.0962 deg`, inside the declared 45 deg
  `initial_filter_entrance.attitude` set. No `1/sqrt(N)` statistical reduction
  is used anywhere on this path;
- the `H18 -> A21` release time **after north lock**: `<= 10.0 s`, by a case
  split on whether the 250-count or the strict 1 s guard binds last.

What they do not close:

- universal per-sample admission at `max_sample_norm_ratio_from_mean = 0.35`
  needs combined hard iron + residual `<= 7 Bmin/47 = 2.9787 uT` against the
  declared `7 uT`;
- the vector sine separation stays a declared PE hypothesis; it is not
  derivable from this class;
- eventual A21 under an arbitrary external `acc_bias_hold_`, and the
  `H18 -> A21` joint24 covariance transport itself;
- the release **from Live**. The counter-owning call sits behind
  `if (mag_ref_set_ && stage_ == Stage::Live)`, so the admitted ungauged timeout
  path never advances it. North lock is the shared unmet prerequisite of the
  yaw gauge and the A21 release;
- the *supply* for the yaw algebra. `tiltOnlyQuatFromBoatQuat_` strips the
  estimator's yaw, not the vessel's, so the accumulation frame turns with the
  boat and a legal heading excursion smears the mean. The parameter is the total
  excursion `delta_tilt + delta_heading`, and nothing declared bounds the
  heading part, so `DECLARED_SUPPLY_ENTRANCE_CLOSED` is false. Never charge a
  heading excursion as tilt, and never read
  `world_averaged_gravity_direction_error` as a total excursion.

## Two-phase private Mahony certificate

The single-level formulation was infeasible at the padded `8.8 m/s^2` envelope:
with the seed angle taken from the source's own `asin(8.8/g) = 1.1137 rad`, seed
containment needed `C >= 1.8540056` while the 87 deg chart allowed
`C < 1.4912551`. The boundary flow still closed with margin `0.0357940`, so the
formulation, not the enclosure, was at fault.

Because `q` and `sup` in

`Vdot = -sqrt(C) [ sqrt(C) q - 2 sup ]`

are level-*independent*, inward flow is a lower bound on the level and holds on
*every* level above `W_in = 2 sup_max/q_min`. The certificate therefore uses two
nested levels of one better-conditioned metric `p=3, c=11`
(`P=[[1,-3],[-3,130]]`, `det=121`):

- seed `1.6707881129` <= outer `1.7689` < chart ceiling `1.8726722863`
  (headroom `5.5414%`), outer margin `0.0262974`;
- inner level `sqrt(C_in)=1`, margin `0.0134738`, `W_in = 0.8578833`;
- capture rate `q_min/2 = 0.0193788 /s`, time constant `51.6029 s`.

Certified tilt: all-time `84.7161 deg`, at the deployed 150 s horizon
`58.2102 deg`, asymptotic `56.6779 deg`. The dependent chain closes again —
binary32/discrete charge, frontend transition, `ou3_startup_timeout_capture`
(branch margin `4.1379 deg`) and `ou3_p4_brmm_frontend_predecessor_invariant`.

Keep the metric read from the continuous certificate, never hardcoded, and keep
the first-order decrease at **one** factor of `sqrt(C)`.

## The private observer cannot supply the magnetic gates

The shipping accumulation frame is the private observer's own tilt
(`tiltOnlyQuatFromBoatQuat_(attitudeReferenceQuat_())`) and the handoff seed is
`boatQuatWithAbsoluteYaw_(q_proxy, pending_yaw_abs_rad_)`, so the proxy tilt
drives both the gauge error and the seed tilt. Derived requirements, exact over
the rationals:

- `min_horizontal_fraction = 0.05` binds at `delta <= 0.0495289` rad
  `= 2.8378 deg`;
- non-vanishing north binds at `delta < 0.1067173` rad `= 6.1145 deg`;
- certified all-time `84.7161 deg` (shortfall `29.8528x` / `13.8551x`) and
  asymptotic `56.6779 deg` (`19.9725x` / `9.2695x`); the certified tilt alone
  exceeds the declared 45 deg entrance by `1.8826x`, so no yaw gauge, however
  accurate, repairs the entrance from this supply.

The declared `world_averaged_gravity_direction_error = 0.02 rad` satisfies the
binding requirement but is the low-passed world-gravity direction error, not the
private observer's tilt error. Do not identify the two.

This is a dead end for the whole route, not just for a metric choice. The
certificate quantifies over the sector value `s` as an independent parameter, so
its conclusion must hold for `s=1`, where the deployed gains give

`theta/r = -(0.01 + 0.1 j w)/(0.01 s - w^2 + 0.1 s j w)`.

An admitted primitive-bounded sinusoid peaks *inside* the declared band, at
`0.0186 Hz` against `[0.018, 0.88] Hz`, contributing `0.14678877` rad; an
admitted DC mean chord superposes `0.05` rad. So no metric, level or subdivision
depth can certify a tilt below `0.19678877 rad = 11.2752 deg` — `3.9732x` the
binding requirement and `1.8440x` north capture. Do not retry a static quadratic
metric for the magnetic gates.

## Fail-closed state and next work

`P4_PASS=false` and `P5_MAY_START=false`. Do not infer theorem closure from the
H18 information lemma, the magnetic or startup certificates, or green CI.

Immediate open items for the next PR are:

1. Attempt the tilt certificate with a frequency-dependent multiplier/IQC on the
   primitive channel and measure the distance to the `11.2752 deg` floor. The
   existing `ou3_p4_affine_hard_tube_iqc.py` and
   `ou3_brmm_acceleration_moment_iqc.py` carry that machinery.
2. Decide the magnetic side: narrow `MAG-BMM150-DET-v1` to a commissioned band
   that forces the gates at the achieved tilt (`<= 55 uT` total needs
   `>= 21.4462 uT` horizontal at the floor; all four tabulated rows in the
   research state are feasible), or tighten the forcing qualification from the
   COMPLETE-BRMM spectral support. Then close the ungauged 150 s timeout branch:
   `ready_by_timeout` does not require `north_ready`, so that handoff still
   takes `proxy_handoff_yaw_sigma_free_rad` with no yaw gauge at all.
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
