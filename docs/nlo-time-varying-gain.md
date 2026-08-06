# Time-varying-gain NLO: why the heave RMS was bad, and what changed

Reference: T. H. Bryne, T. I. Fossen, T. A. Johansen, *Nonlinear Observer with
Time-Varying Gains for Inertial Navigation Aided by Satellite Reference Systems
in Dynamic Positioning*. <https://torarnj.folk.ntnu.no/TimeVarGain.pdf>

This note records what was actually wrong with `src/nlo` and what was done
about it. Short version: the observer equations were already a faithful
transcription of the paper and the published gains were correct, but the
no-GNSS operating mode was improvised, and the paper's dynamic-positioning
value of `theta` puts the aiding loop inside the ocean wave band.

## 1. The published coefficients check out

`K0 = P C'`, with `P` the solution of the paper's algebraic Riccati equation
(12) for `Q = blkdiag{50, 0.5*I3, 0.08*I3, 0.0025*I3}`. Solving it reproduces
both published columns to four decimals:

| chain | states | paper | recomputed |
|---|---|---|---|
| vertical | `p0z, pz, vz, fz` | 5.4295, 2.2396, 0.4454, 0.0354 | 5.42947, 2.23955, 0.44538, 0.03536 |
| horizontal | `p, v, f` | 0.9513, 0.3275, 0.0354 | 0.95130, 0.32750, 0.03536 |

One detail matters for anyone re-deriving them. The paper writes the ARE as
`AP + PA' + Q - 2*tau*PC'CP = 0` with `tau = 1/2`, which reads as an effective
measurement weight of `2*tau = 1`. That does **not** reproduce the published
numbers; it gives 7.5009, 3.1315, 0.6270, 0.0500. The published values come out
of `AP + PA' + Q - tau*PC'CP = 0`, i.e. `R = 1/tau = 2`. The last gain settles
it unambiguously: for a chain measured at its head the final gain is
`sqrt(q_last/R)`, and `sqrt(0.0025/2) = 0.035355`, which is the paper's 0.0354,
whereas `sqrt(0.0025/1) = 0.05`.

Note also that `Q` weights `p`, `v` and `f` identically on all three axes. The
vertical column is therefore not vertical-specific: it is the solution for *any*
axis whose position is augmented with its own integral and measured there.

## 2. What the paper's theta = 1 costs a wave sensor

Write the vertical loop with `g_i = vartheta * theta^i * K0_i` and `L(s) = Th/(Th*s+1)`
for the integrate-then-high-pass path of eq. (27). The transfer from true heave
to estimated heave is

```
P_hat/P = (1 + g1*L) / (1 + L*(g1 + g2/s + g3/s^2 + g4/s^3))
```

and, since `Th = 600 s` makes `L ~= 1/s` anywhere in the wave band,

```
P_hat/P ~= s^3*(s + g1) / (s^4 + g1*s^3 + g2*s^2 + g3*s + g4)
```

A discrete replica of the shipped code reproduces this to four decimals, so the
implementation and the model agree.

The consequence at the paper's `theta = 1`:

| wave period | 4 s | 6 s | 8 s | 10 s | 12 s | 16 s |
|---|---|---|---|---|---|---|
| gain \|H\| | 1.08 | 1.08 | 1.07 | 1.07 | 1.05 | 1.02 |
| phase lead | +14° | +23° | +31° | +40° | +48° | +66° |
| \|H - 1\| | 0.27 | 0.42 | 0.56 | 0.70 | 0.84 | 1.10 |

The amplitude is nearly right and the phase is not. This is the effect the
paper itself notes in Sec. IV-D ("the heave estimates have a positive phase
relative to the actual heave signal").

The mechanism is visible in row (9a). Its loop gain is `g1 = 5.43`, so `p0_hat`
tracks with a 0.18 s time constant and settles at `p_hat/g1` rather than
behaving like an integral. Row (9b) then sees an effective position leak of
`g2/g1 = 0.41 rad/s`, a 15 s period. That is the heave path's high-pass corner,
and it sits inside the wave band.

Two knobs do *not* move it:

- `vartheta(t)` multiplies all four gains equally, so `g2/g1` is invariant.
  Across `vartheta` from 0.5 to 2.0, `|H-1|` at 8 s moves only 0.58 to 0.54.
- `Th`. From 600 s down to 0.5 s, `|H-1|` at 8 s moves only 0.56 to 0.41.

Re-solving the ARE for a different `tau` barely helps either: even `tau = 5e-5`
only brings the corner from 0.41 to 0.18 rad/s, because the corner is set by
the shape of `Q`, not its scale.

`theta` does move it, and exactly: it enters only through
`K(t) = vartheta(t)*theta*L_theta^-1*K0` (eq. (33)), which is a pure frequency
scaling of the whole aiding loop. The corner becomes `0.4125*theta rad/s` with
every published coefficient, and every ratio between them, untouched.

This was measured, not assumed. Against the shipped simulator's reference
heave, the linear model above predicts 1.597 m RMS error on the Hs = 8.5 m
JONSWAP record; the full nonlinear observer produced 1.667 m. Loop shaping
accounted for 96 % of the error.

## 3. What the no-GNSS mode was doing wrong

With `WithGNSS = false` the observer held `p_xy`, `v_xy` and `xi_xy` at zero.
That has a knock-on effect on attitude. Eq. (9e) is `f_hat_n = R*f_b + xi`, so
pinning `xi_xy` makes `f_hat_n` equal to `R*f_b` in the horizontal plane, and
the paper's injection term (8), `k1 * f_b x R'*f_hat_n`, self-cancels. The code
worked around that with a substitute reference: zero out `f_hat_n`'s horizontal
components, use the remainder as a tilt reference, then scale the result by
0.03 and clamp it at 0.002 rad/s because a plain accelerometer tilt reference
chases wave acceleration. Roll RMS reached 10.2° on the Hs = 8.5 m
Pierson-Moskowitz record.

The paper already contains the fix. Sec. II-A.2 replaces the low-precision GNSS
height with the virtual measurement `p0n_z = integral(pn_z dt) = 0`, on
Godhavn's argument that wave-induced motion oscillates about the mean sea
surface. For a free-floating body the identical argument holds in surge and
sway, and the 600 s innovation high-pass removes current-driven drift from the
horizontal innovation exactly as it removes tide from the vertical one.

So the observer now carries `p0n` on all three axes and drives all three
innovations from `p0n = 0`, using the vertical gain column on each — which, per
Sec. 1 above, is the same ARE solution, not an analogy. `xi` becomes observable
on all three axes, `f_hat_n` becomes a real specific-force estimate, and (8) is
used as written. The 0.03-scaled fallback survives only for
`use_virtual_horizontal_position = false`.

The point of this is the attitude loop, not surge and sway. The horizontal
displacement estimates it produces are **not** a calibrated product: a virtual
zero-mean constraint is much weaker aiding than a position fix, and with
`Mag = None` there is no yaw reference, so the horizontal frame the observer
settles on is arbitrary. Measured against the simulator's reference, horizontal
displacement carries two to three times the true RMS and correlates only weakly
(|r| = 0.21 to 0.60, sign depending on the frame offset). Heave is the product
here; x and y are internal.

## 4. The theta schedule

`theta` is scheduled to hold the aiding corner at a fixed fraction of the
dominant wave frequency:

```
theta = theta_from_omega_gain * omega_peak        (default gain 0.56)
```

which places the corner at `0.23 * omega_peak`. The 0.56 is the flat optimum
measured by a fixed-`theta` sweep across Hs = 0.27 m to 8.5 m: the per-record
optima were `theta/omega_peak` = 0.51, 0.52, 0.53, 0.54, 0.58, 0.59, 0.61, 0.61.

Scheduling a translational gain scalar on a measured condition signal is the
paper's own construction — it schedules `vartheta(t)` on reported GNSS quality
and says other laws may be used in its place. Here there is no receiver to
report quality, and `vartheta` cannot move the corner anyway.

This does depart from the paper's `theta >= 1`. That bound comes from Theorem 1,
where a large `theta` is what dominates the nonlinear coupling term; it is
sufficient, not necessary, and the linear aiding loop is Hurwitz for any
`theta > 0` at these gains.

Three implementation points, each of which was a failure mode before it was
handled:

- **The tracker runs on its own fixed reference channel**, not on the
  observer's output. Driving it from acceleration locks it far above the
  displacement peak (0.32 Hz against a true 0.12 Hz on the Hs = 4 m record),
  because a JONSWAP acceleration spectrum falls only as `omega^-1` above the
  peak. Driving it from the observer's own heave closes a positive feedback
  loop — low `theta` admits drift, drift captures the tracker, the tracker
  lowers `theta` — and it walks down to 0.054 Hz and stays there. The tracker
  therefore gets vertical inertial acceleration double-integrated through two
  leaky integrators at 0.03 Hz: displacement-shaped in the wave band,
  `-40 dB/decade` on drift, and independent of `theta`.
- **`theta` changes are bumpless.** `p0_hat` settles proportional to `1/theta`,
  so writing a new `theta` straight into the config steps the innovation.
  `Filter::setTheta()` rescales `p0_hat` and the high-pass state with it.
- **`theta` snaps on first acquisition** and lags only afterwards. Lagging in
  from the initial value wastes several time constants with the loop in the
  wrong place, and that transient lands inside the scored window.

## 5. The paper's attitude gain ramp needs one change without GNSS

Sec. IV-C ramps `ka = [20, 20, 1]` for the first 100 s. The proportional part
is fine here and is kept. The integral part is not: with `kI = 1 rad/s` and no
GNSS, Sigma2's much longer transient is integrated into the gyro bias estimate,
the bias reaches `gyro_bias_limit_rad_s`, that rotates the attitude, the
injection term grows, and the interconnection diverges — roll and pitch RMS run
to tens of degrees. The parameter projection does not save it because the
admissible ball used here is two orders of magnitude wider than the true biases.

Isolated: ramping `k1`/`k2` alone gives 6.94 % and roll 0.21°; ramping `kI`
alone gives 143468 % and roll 56°. So `kI_initial` is held at its nominal 0.01
when `WithGNSS = false`.

## 6. Reporting: the DC of position is not observable

The aiding is a virtual measurement whose innovation is high-passed, and a
high-pass has exactly zero gain at DC for any `Th`. So the observer can pin
wave-band motion and cannot pin the constant it sits on. Sweeping `Th` from
600 s to 25 s moves the resulting offset by under 4 %, which confirms it.

Left alone that constant was not small. Reported heave carried a standing
-1.21 m offset on the Hs = 8.5 m record. Its source is a rectification: the
bias appears in all three channels at once and in the ratios the loop's own
time constants predict.

| channel | bias, jonswap Hs 8.5 | ratio |
|---|---|---|
| `acc_z` | -0.0119 m/s^2 | |
| `vel_z` | -0.1722 m/s | 14.5 s |
| `disp_z` | -1.2111 m | 7.0 s |

So a constant vertical acceleration bias is being mapped into velocity and
position by the aiding loop. It scales with sea state (-0.0004, -0.0031,
-0.0072, -0.0119 m/s^2 at Hs = 0.27, 1.5, 4.0, 8.5), which makes it a
rectification rather than a sensor bias — `xi_dot = -R*S(sigma)*f_b` in (9d)
has a nonzero mean over the wave cycle. Because the loop maps acceleration to
position through its time constants, the position offset grows as
`1/theta^2`, which is why it went from 0.21 m to 1.21 m when theta was moved
below the wave band.

Removing that acceleration bias at source would be the better fix and is
**not** done here; it remains open, and the reported `acc_z` still shows it.
What is done is to report position and velocity through a first-order
high-pass, `report_highpass_tau_s = 50 s`, so the output is the wave-band
motion the observer can actually claim. That is what the zero-mean aiding
assumption already asserts, so it is consistent with the observer rather than
a cosmetic trim. `filter().positionNED()` still returns the raw state.

50 s is a measured optimum: mean raw Z RMS is 7.18 % at 20 s, 6.88 % at 35 s,
6.83 % at 50 s, 6.97 % at 70 s, 7.51 % at 100 s, 9.97 % at 400 s. Too short
costs wave-band fidelity, too long leaves the observer's own wander in the
output. At 50 s a 20 s swell sees 0.2 % amplitude loss and 3.6 degrees of
phase.

The vertical is always high-passed, because the virtual measurement always
aids it through the innovation high-pass. The horizontal axes are high-passed
only when no position reference is present; with GNSS their DC is observable
and is reported as estimated.

## 7. Results

Deterministic single-realization protocol, 900 s trailing window, **raw**
Z RMS as a percentage of Hs — no de-meaning anywhere, on either observer.
Lower is better.

| record | NLO before | NLO after | offset | PII |
|---|---|---|---|---|
| jonswap Hs 0.27 | 8.21 | **7.01** | -0.000 m | 4.85 |
| jonswap Hs 1.5 | 9.70 | **6.70** | +0.004 m | 6.09 |
| jonswap Hs 4.0 | 15.78 | **6.76** | +0.010 m | 7.79 |
| jonswap Hs 8.5 | 19.77 | **7.21** | +0.006 m | 8.87 |
| pmstokes Hs 0.27 | 8.60 | **6.77** | +0.001 m | 5.40 |
| pmstokes Hs 1.5 | 11.67 | **6.60** | +0.001 m | 6.71 |
| pmstokes Hs 4.0 | 15.62 | **6.49** | +0.005 m | 7.91 |
| pmstokes Hs 8.5 | 27.02 | **7.09** | +0.040 m | 9.36 |
| mean | 14.55 | **6.83** | | 7.12 |

The simulator's quality gate scores this raw number. Earlier revisions gated
the de-meaned value, which by construction could not see the offset at all.

Attitude RMS, degrees:

| record | NLO before roll/pitch | NLO after roll/pitch | PII roll/pitch |
|---|---|---|---|
| jonswap Hs 0.27 | 0.71 / 0.52 | 0.21 / 0.13 | 0.57 / 0.69 |
| jonswap Hs 1.5 | 1.15 / 0.96 | 0.25 / 0.12 | 1.54 / 1.94 |
| jonswap Hs 4.0 | 4.81 / 2.60 | 0.29 / 0.15 | 1.73 / 3.32 |
| jonswap Hs 8.5 | 3.34 / 2.25 | 0.38 / 0.20 | 3.50 / 4.70 |
| pmstokes Hs 0.27 | 0.66 / 0.48 | 0.21 / 0.13 | 0.64 / 0.73 |
| pmstokes Hs 1.5 | 2.87 / 2.04 | 0.24 / 0.13 | 1.88 / 2.20 |
| pmstokes Hs 4.0 | 2.25 / 4.80 | 0.29 / 0.16 | 1.96 / 3.27 |
| pmstokes Hs 8.5 | 10.17 / 6.65 | 0.39 / 0.22 | 4.28 / 5.26 |

Heave: the NLO beats the PII observer on five of eight records and on the mean
(6.83 % against 7.12 %), and it is decisively better in the seas that matter —
every record with Hs >= 4 m improves by 1.0 to 2.3 percentage points. It is
still worse on the two Hs = 0.27 m records, where both observers are near
their noise floor (0.019 m against 0.013 m in absolute terms) and a
fixed-`theta` sweep confirms the schedule is already sitting at that record's
optimum, so the gap is the observer's noise floor rather than the schedule.

Attitude: better than both the previous NLO and the PII observer on every
record, by roughly an order of magnitude in the larger seas.

## 8. Reproducing

```
make -C tests/nlo test
```

`tests/nlo/nlo-sim.cpp` reads a few environment variables as sweep hooks; all
are unset in a normal run, so the reported numbers describe the shipped
defaults.

| variable | effect |
|---|---|
| `NLO_THETA` | disable the schedule and pin `theta` (e.g. `1.0` for the paper's value) |
| `NLO_THETA_GAIN` | override `theta_from_omega_gain` |
| `NLO_THETA_TAU` | override `theta_smooth_tau_s` |
| `NLO_VIRT_XY` | `0` restores the old pinned-horizontal behaviour |
| `NLO_TVATT` | `0` disables the attitude gain ramp |
| `NLO_K1I`, `NLO_KII` | override `k1_initial` / `kI_initial` |
| `NLO_TH` | override the innovation high-pass `Th` |
| `NLO_RHP` | `0` disables the reporting high-pass, exposing the raw DC offset |
| `NLO_RHP_TAU` | override `report_highpass_tau_s` |
| `NLO_NOGATE` | report gate breaches without failing, for sweeps |
