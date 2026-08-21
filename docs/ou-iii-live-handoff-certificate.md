# The Live-handoff certificate, and the gap it does not close

Phase 1 of this work gave the Live error dynamics a block-local ISS theorem.
Phase 2 made its linear half constructive. Phase 3, which is what this note is
about, tries to make the whole entrance condition a number a running filter can
evaluate — and reports what happens when you do.

The short version: the theorem's constants improved by about three orders of
magnitude, the handoff bounds are now constructed rather than assumed, the true
handoff error respects every one of them across 91 truth-based scenarios, and
the basin inequality still does not close — by a factor of about 50 at the
calmest reference sea and about 2100 at the roughest. That factor is reported
rather than tuned away.

## What changed in the mathematics

Phase 2 stated its uniform-exponential constants in a fixed dimensionless
coordinate system, so the covariance contraction the deployed schedule actually
exhibits had to be converted back into those coordinates. That conversion costs
`sqrt(pbar/punder)`, which the diagnostic measures as 556, and it is what made
the Phase-2 basin radius meaningless.

Phase 3 stays in the covariance metric. Two facts make that possible.

**The metric norm never increases.** For a Joseph-form update with
`A = (I-KC)F`,

    A P A^T = P+ - K R K^T - (I-KC) Q (I-KC)^T  <=  P+,

so `||L+^{-1} A L||_2 <= 1` at every sample, with no hypothesis about the sea.
The prefix constant of the horizon certificate is therefore exactly 1, and the
transition constant `M_H = rho_H^{-(H-1)}` is 1.020 where Phase 2 had 566. The
diagnostic checks `alpha_max <= 1` at all eight reference operating points and
fails the build if it is ever exceeded, because a violation would mean the
implementation's update is no longer the update the lemma is about.

**The metric is scale-invariant.** `||D e||_{(D P D')^{-1}} = ||e||_{P^{-1}}`,
so the certificate does not depend on the choice of physical scales at all.
The Phase-2 scaling survives only as a reporting convention.

On top of that, the nonlinear remainder is no longer priced at the slowest
mode's memory, which at the deployed schedule is 3.45e5 samples. Each remainder enters through a known three-dimensional channel
— the MEKF reset through attitude, the measurement-model remainders through
their own Kalman gains — so each gets its own `l1` injection gain, computed by
forward impulse propagation with nothing to store. The scalar alternative
`M_H/(1-rho_H)` is 3.45e5 samples; the measured directional gains are 1.35e7
per radian for the attitude channel, 2.0e3 per m/s^2 for the accelerometer
channel and 3.9e2 per uT for the magnetometer channel.

## What changed in the filter

Only what a constructed handoff bound needs.

- **The integral epoch restarts at the Live transition.** `S` is the running
  integral of `p` and its lower limit is free: `p`, `v` and the `S=0`
  pseudo-measurement are all invariant under `S -> S + const`. Declaring the
  epoch to be the handoff makes the handoff error in that coordinate exactly
  zero, which is the only bound in the whole certificate that is an identity
  rather than an assumption. The marginal is deliberately left at the
  constructor's value: a zero numerator contributes nothing whatever the
  denominator is, so collapsing it would change the pseudo-measurement's
  transient for nothing the theorem can use.

- **`a_w` is seeded from the proxy-levelled specific force.** The filter's own
  measurement model is `f_b = R_wb (a_w - g)`, so `a_w = R_bw f_b + g` in the
  proxy attitude is a pure function of the measurements — no MEKF state enters,
  so the startup path stays exogenous. Entering Live with `a_w = 0` declares
  the platform still, which in a seaway it is not. Measured against truth the
  seed leaves 0.04 to 0.51 m/s^2 against a constructed bound of 1.95 m/s^2.

- **The certificate itself**, `src/kalman_ou_iii/LiveEntranceCertificate.h`,
  with every term exposed: the constructed bound for each block, what that
  bound rests on, both sides of the basin inequality, and the failure reason.

- **`LiveCertified` vs `LiveUncertified`.** Evaluating the certificate never
  withholds a handoff in the deployed configuration; what it decides is whether
  the interval may be *reported* as covered by the theorem. A timeout-forced
  handoff is never certified. Certification is cleared at the interval
  boundaries the paper already recognises: the 70 degree tilt re-lock, the
  second magnetic re-gauging, and any return to Cold.

The scored accuracy of the filter is unchanged. Against `main` on the eight
reference records the roll, pitch and yaw RMS agree to the fourth decimal.

## Two things that were tried and rejected on measurement

**Seeding the gyro bias from the startup proxy.** The Mahony integral term does
estimate a constant bias at rest. In a seaway it also winds up against the
horizontal orbital acceleration, and the wind-up is the larger of the two by an
order of magnitude: seeding from it leaves 1.6e-3 to 3.3e-3 rad/s of true
gyro-bias error at handoff, against 4e-5 to 1e-4 rad/s when the state is left at
zero, for turn-on draws of order 1e-4. The estimate is worse than the thing it
estimates. The seed is implemented and disabled, and the bound stays the
declared turn-on envelope.

**Making the seed covariance agree with the declared envelope.** The
constructor seeds the residual accelerometer bias at 0.004 m/s^2 while the
declared calibration envelope only bounds the true error at 0.05, and
re-synchronises `a_w` to an OU stationary spread an order of magnitude below
its constructed bound. That is overconfidence in the ordinary estimator sense
and it is a real finding. Fixing it, however, moves the scored pitch RMS on the
pmstokes H4.0 record from 0.1811 to 0.2102 degrees against a committed limit of
0.1975 — the inflated bias marginal lets the filter chase that bias through the
wave band. Since the certificate does not close either way, taking a measured
accuracy regression to buy a factor the theorem still cannot use would be
paying for nothing. The switch exists (`seed_covariance_from_envelope`) and is
off.

**A translational bootstrap for `v` and `p`.** Not built, because it was
measured first. With both entering Live at zero the true errors at handoff over
the reference seas are 0.13 to 3.33 m/s and 0.04 to 4.12 m, inside a
wave-orbital envelope that has to cover the 8.5 m sea in any case, and the
kinematic block enters the basin inequality through `M_xi_ell` rather than
through the tube.

## The numbers, and the gap

Small-gain slope `c_eff` over the eight reference operating points, with
envelopes taken over a 600 s established interval:

| point | c_eff | r_cert | budget |
|---|---|---|---|
| J0.27 | 0.694 | 1.441 | 0.353 |
| J1.50 | 1.067 | 0.937 | 0.230 |
| J4.00 | 12.53 | 0.0798 | 0.0196 |
| J8.50 | 29.28 | 0.0341 | 0.0084 |
| P0.27 | 0.681 | 1.468 | 0.360 |
| P1.50 | 0.988 | 1.012 | 0.248 |
| P4.00 | 7.007 | 0.1427 | 0.0350 |
| P8.50 | 26.06 | 0.0384 | 0.0094 |

Against those, the constructed handoff error in metric units puts the left side
of the entrance inequality at 18.17 with the deployed seed, so the inequality
would close only for `c_eff < 0.0138`; with the coherent seed the same bounds
give 4.24 and `c_eff < 0.0589`. Both are printed by
`live_entrance_certificate-test` rather than typed here. The certificate is
short by a factor of about 50 at the calmest reference point and about 2100 at
the roughest with the deployed seed — about 12 and 500 with the coherent one —
and it certifies nothing in the truth-based sweep.

What is left is two structural terms rather than slack. `Gamma_theta`, the
roughly 100 s memory of the attitude channel, multiplied by the attitude
correction that the `S` pseudo-measurement and the magnetometer produce; and
`Gamma_a * (||dtheta|| ||da_w||)`, the accelerometer residual remainder. Both
grow with sea state, which is why the calm and rough reference points differ by
a factor of 43. Closing the rest needs either a structural change to the
`S`-to-attitude cross-gain — the same change the Phase-1 remark identifies as
what would move `S` out of the sensitive block — or a sharper treatment of the
quadratic remainder than one isotropic coefficient. Neither is attempted here.

## What the validation actually establishes

`tests/kalman_ou_iii/live_handoff_validation.cpp` runs 91 scenarios: eight seas
synthesised at the committed reference `(Hs, Tp)` pairs by three seeds by three
wave phases with bounded gyro- and accelerometer-bias draws, plus large initial
roll/pitch/yaw, attitudes 170 degrees from level, delayed magnetometer
availability, sea state changing during startup, magnetic reference errors at
the edge of the declared envelope, and runs with no magnetometer at all.

Because the certificate currently certifies nothing, checking bounds only on
certified handoffs would be vacuous. The test therefore checks the
unconditional form: on every handoff, each finite constructed bound must
actually bound the true error. It does, on all of them, with the worst
conservatism ratio at 1.07. The `S` identity is checked separately and is
exactly zero everywhere, and a handoff with no magnetic gauge is required never
to be certified.

The declared envelope in that test is wider than the library defaults on
purpose. An envelope is a declaration, and a deployment has to declare one that
covers what it will see; declaring one the scenarios then leave is not a
stricter test, it is applying the theorem outside its own hypotheses. The bias
draws are uniform inside a ball rather than Gaussian for the same reason — a
Gaussian leaves any finite envelope eventually, and a CI failure for that reason
would be the theorem correctly declining to cover a case, not a defect.
