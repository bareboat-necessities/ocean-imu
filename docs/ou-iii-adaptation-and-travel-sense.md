# OU-III adaptation loss and wave travel-sense reporting

Investigation record and remediation plan for the two reviewer findings on the
OU-III manuscript:

1. OU-III improves vertical RMS but degrades 3D RMS in moderate and large seas,
   and adapting `tau`/`sigma_aw` is worse than adapting `r_S` alone.
2. The reported travel-sense classes are defined against the estimator's own
   axis representative rather than against physical propagation.

Everything below is measured on the versioned simulation records
(`bareboat-necessities/oceanography-waves-lib`, `v1.1.3`) with
`tests/kalman_ou_iii/kalman_ou_iii-sim`, scoring the trailing 900 s.

## 1. Finding 1: where the adaptation actually loses

### 1.1 The leading hypothesis was posterior/process covariance inconsistency

It is not the cause. Two independent checks:

- The committed study already prices the policy. `Adaptive` versus
  `AdaptiveHeldCovariance` (identical tuning, periodic re-alignment of the
  posterior `a_w` marginal switched off) differ by 0.002 m in 3D RMS at
  `H_s = 8.5` m — 0.07% of a 3.10 m error.
- A direct probe (40 000 steps, `sigma_aw` swept over its adaptation range,
  marginal re-aligned every 0.1 s) shows the minimum eigenvalue of `P` staying
  at float round-off both before and after the re-alignment
  (`-6.8e-05` before, `-2.3e-05` after, against diagonal entries of order 1).
  The re-alignment *inflates* the marginal, typically by 40-100x, and adding a
  PSD term to a diagonal block cannot destroy PSD.

The re-alignment does discard a large amount of posterior information: the
inflation factor is not "bounded" in any useful sense.

Making the operation *consistent* makes it worse, which settles the question
from the other side. `W3D_AW_COV_SYNC=congruent` re-aligns the marginal through
a congruence `x_aw -> A x_aw`, `A = L_new L_old^-1`, which reaches the same
marginal, leaves the whitened cross-covariance untouched and stays PSD by
construction. On `jonswap H8.500`:

| a_w sync policy | disp Z | 3D |
| --- | --- | --- |
| `periodic` (deployed overwrite) | 0.645 | 2.767 |
| `congruent` | 0.678 | 3.618 |
| `congruent`, heading rotated 45 deg | 3.655 | 11.849 |

The congruence propagates the marginal inflation into the cross-covariances,
which is what consistency demands and which the filter cannot absorb. The
deployed overwrite is not better because it is more correct; it is better
because it silently damps an intervention that is too large. The conclusion is
that the periodic re-alignment should be retired or bounded, not made
self-consistent - and the `HeldCovariance` ablation already shows retiring it
costs nothing.

### 1.2 The frequency that drives `tau` does not track the sea state

`tau_target = tau_coeff * 0.5 / f`, with `f` from the frequency tracker running
on the low-passed body-Z acceleration proxy. Acceleration spectra are
`(2*pi*f)^4` weighted, so their apparent frequency sits far above the elevation
peak, and the offset grows with the sea state:

| record | `f_p` (elevation) | `f_mean` (acceleration) | ratio |
| --- | --- | --- | --- |
| jonswap H0.270 L14.047 | 0.335 Hz | 0.548 Hz | 1.63 |
| jonswap H1.500 L50.710 | 0.177 Hz | 0.474 Hz | 2.68 |
| jonswap H4.000 L112.766 | 0.118 Hz | 0.440 Hz | 3.72 |
| jonswap H8.500 L202.839 | 0.085 Hz | 0.417 Hz | 4.90 |

The acceleration-band frequency is nearly constant across a 4x change in sea
state. The channel the filter uses to learn "how long the waves are" carries
almost no sea-state information — only estimation noise. Across the ten seed
triplets of the committed study at `H_s = 8.5` m, `tau_applied` is
`1.660 +/- 0.431` s against the `2.330` s the same tuner produces on the
noise-free record, and `r_S = clip(1.2 * sigma_aw * tau^3)` cubes that error:
`11.5 +/- 4.3` against `34.7`.

### 1.3 The integral regularizer is therefore far too strong in developed seas

The `S = 0` pseudo-measurement is a high-pass on displacement whose corner is
set by `r_S`. With `r_S` too small the corner sits inside the wave band and the
filter removes the signal it is supposed to estimate. Scoring
`jonswap H8.500` (900 s window, fixed `tau = 2.33`, `sigma_a = 2.283`):

| `r_S` (m*s) | `R_S_xy` | disp X | disp Y | disp Z | 3D |
| --- | --- | --- | --- | --- | --- |
| 3 | 0.36 | 3.478 | 2.170 | 1.489 | 4.362 |
| 8 | 0.36 | 2.857 | 1.834 | 0.841 | 3.498 |
| 15 | 0.36 | 2.221 | 1.476 | 0.590 | 2.732 |
| 35 | 0.36 | 1.475 | 1.080 | 0.405 | 1.873 |
| 3 | 1.00 | 2.835 | 1.798 | 1.343 | 3.616 |
| 8 | 1.00 | 1.848 | 1.266 | 0.782 | 2.372 |
| 15 | 1.00 | 1.353 | 1.023 | 0.568 | 1.789 |
| 35 | 1.00 | 0.964 | 0.902 | 0.402 | 1.380 |

The response is monotone up to `MAX_R_S = 35`, which is where `FixedOracle`
sits (34.66) and is exactly why `FixedOracle` wins: the clamp binds at the
operating point the filter needs.

### 1.4 The horizontal axes are pushed further into the same failure

`R_S_xy_factor = 0.36` scales the horizontal pseudo-measurement standard
deviation down, i.e. it makes the horizontal high-pass 2.8x *stronger* than the
vertical one, in the channel where the corner is already too high. Adaptive
mode, sweeping only that factor:

| record | 3D at `R_S_xy = 0.36` | 3D at `R_S_xy = 1.00` |
| --- | --- | --- |
| jonswap H0.270 | 0.068 | 0.107 |
| jonswap H1.500 | 0.289 | 0.291 |
| jonswap H4.000 | 1.188 | 0.891 |
| jonswap H8.500 | 2.767 | 1.829 |

0.36 is a small-sea optimum applied to every sea state. It costs 25% of the 3D
RMS at `H_s = 4` m and 34% at `H_s = 8.5` m.

By contrast `S_factor = 1.87` (the horizontal stationary acceleration scale) is
mis-specified relative to physics — the records give per-axis
`sigma_ax/sigma_az = 0.81` and `sigma_ay/sigma_az = 0.55`, with
`sqrt(sigma_ax^2 + sigma_ay^2)/sigma_az = 0.99` as deep-water theory requires —
but it is *not* what drives the loss: sweeping it 1.87 -> 0.80 -> 0.60 moves
3D RMS by under 2%.

### 1.5 Mechanism, stated once

The tuner learns the sea state through a frequency estimate taken in the wrong
spectral band. `tau` therefore stays roughly constant while the true wave period
grows, `r_S ~ sigma * tau^3` stays far below what the sea requires, the
`S = 0` regularizer's corner ends up inside the wave band, and the horizontal
axes get an extra 2.8x of the same medicine. The vertical channel still beats
the fixed reference because `r_S` at least moves with `sigma`; the horizontal
channel does not, so 3D RMS degrades exactly where the paper reports it.
Freezing `tau`/`sigma_aw` at the nominal point (`AdaptiveRSOnly`) helps because
it stops the OU model from following a frequency signal that carries no
information, while `r_S` keeps tracking `sigma`.

## 1.6 The fix

`src/tuner/WavePeriodEstimator.h` estimates the zero-crossing period
`T_z = 2*pi*sqrt(m0/m2)` of the elevation from leveled vertical acceleration.
Two first-order high-pass stages and two leaky integrators, all at the same
corner `lambda`, give band-limited velocity and elevation proxies. They differ
by exactly one integrator, so

    sigma_v / sigma_eta = sqrt(omega^2 + lambda^2)

holds for a narrow band and inverts exactly; a broadband input returns the same
relation averaged over whatever the shared response passes. The relation is
invariant to any filtering the two proxies have in common, which is what makes
the high-pass stages free: double integration weights a spectrum by
`omega^-4`, so without them the sub-band energy of a real strapdown signal
dominates the elevation proxy. Measured against the record truth:

| record | `T_z` truth | `T_z` estimated | error |
| --- | --- | --- | --- |
| jonswap H0.270 | 2.51 s | 2.60 s | +3.8% |
| jonswap H1.500 | 4.48 s | 4.37 s | -2.4% |
| jonswap H4.000 | 6.60 s | 7.09 s | +7.4% |
| jonswap H8.500 | 8.55 s | 8.43 s | -1.5% |

Two implementation details mattered more than expected.

The input is the *leveled* vertical acceleration the direction stage already
forms, not the body-Z proxy. Fed the body-Z proxy the estimator reported
6.8-10.0 s across a family whose true `T_z` runs 2.5-8.6 s: the sub-band
gravity leakage of a tilting platform swamps a doubly integrated signal. The
leveled signal depends on attitude but not on the linear block, so it does not
close a loop through the quantity being tuned. Reading the filter's own
displacement would: the integral pseudo-measurement high-passes displacement,
which raises the apparent frequency, which shortens `tau`, which strengthens the
pseudo-measurement.

### Opening the loop

Levelling with the filter's own attitude is what puts the tuner inside a loop.
Two inputs open it, both pure functions of the measurements, and
`tools/ou_wave_period_input_study.py` replays all eight reference records under
each of them (`W3D_WAVE_PERIOD_INPUT`, `setWavePeriodInput()`) on the
deterministic single-seed protocol of `tools/ou_sim_table.py`. The second one
is now the default; ratios below are still taken against attitude levelling,
since the question is what changed relative to it.

**`body_z`** is the raw proxy the frequency tracker already runs on. It opens
the loop and is not usable:

| filter | vertical RMS | 3D RMS | worst record |
| --- | --- | --- | --- |
| OU-II | 1.88x | 1.63x | +350% vertical, PM-Stokes `H_s = 0.27` m |
| OU-III | 2.51x | 2.24x | +625% vertical, PM-Stokes `H_s = 0.27` m |

The first two columns are geometric means of the per-record ratio
body-Z / leveled. Degradation is monotone in sea state and worst in the calmest
records, which is the reported-period error read straight through: the body-Z
estimate is pinned near 6.8-10.0 s whatever the sea does, so it is 2.5-3.5x too
long at `H_s = 0.27` m and only 15-21% too long at `H_s = 8.5` m. `tau` is a
fixed multiple of `T_z`, so a period that never moves is the uninformative
operating point of section 1.2 with extra steps. Attitude RMS is essentially
unmoved either way, confirming the loss is the operating point and not the
attitude solution.

**`complementary`** is the measurement-only *leveled* signal the paragraph above
asks for, and is what the filter now ships.
`src/tuner/VerticalAccelComplementary.h` runs a private Mahony observer on the
raw gyro and accelerometer - never the calibrated values, never any filter
state - and reports `-((R f)_z + g)` from its own quaternion. It costs nothing:

| filter | vertical RMS | 3D RMS | largest per-record deviation |
| --- | --- | --- | --- |
| OU-II | 1.000x | 1.000x | 0.2% |
| OU-III | 1.000x | 1.000x | 0.1% |

Across all sixteen record-filter pairs the reported `T_z` matches the leveled
value to 0.01 s and the RMS to within 0.2%, which is replay noise. The loop is
gone for free.

The gain matters and the direction of the sensitivity confirms why. `two_kp`
sets the rate at which the accelerometer corrects the gyro, roughly a corner at
`two_kp/2` rad/s, and it must sit *below* the wave band: a fast observer chases
horizontal orbital acceleration into tilt, which is the body-Z leakage failure
by another route. The default 0.2 puts the corner near 0.016 Hz against a
0.11-0.42 Hz band. Sweeping it on the `H_s = 8.5` m JONSWAP record moves
reported `T_z` 8.42 -> 8.63 s and vertical RMS 0.3337 -> 0.3399 m as `two_kp`
goes 0.05 -> 2.0, i.e. monotonically worse as the corner climbs into the band,
and flat below 0.5. Bias is deliberately not estimated (`two_ki = 0`): at the
reference bias range a constant gyro bias leaves a static tilt error near
0.5 deg, which the estimator's two high-pass stages reject anyway, while an
integral term is another slow state that can wind up against a sustained
horizontal acceleration.

`tests/kalman_ou_iii/tuner_coupling-test.cpp` asserts the exogeneity
bit-for-bit rather than by tolerance, and does so on the default path with
nothing selected: displacing *every* estimator state - attitude and all linear
states - leaves the reported period and `tau_target` numerically identical.
The `Leveled` ablation, which must now be selected explicitly, moves them
8.05 -> 10.32 s and 1.28x under the same displacement.

This closes the item `app:iss-tuner` in the stability appendix was written
around. Both tuner channels are now exogenous, so the estimator and the tuner
do not form a feedback interconnection and no composite Lyapunov function or
small-gain condition is needed to describe the deployed schedule. The appendix
records the loop as it stood, and what removing it changed.

The estimator settles in about a minute, so the tracker frequency is used until
it is ready.

With `f_tune = 1/T_z` the existing law `tau = c_tau/(2 f_tune)` becomes
`tau = c_tau * T_z / 2`, and the same substitution fixes the tuner's variance
horizon, which is a fixed number of periods and had been shorter than one wave
period. The coefficients were re-fitted on the four stationary JONSWAP records
over a 180-point grid:

| parameter | was | now | why |
| --- | --- | --- | --- |
| `tau_coeff` | 1.38 | 1.0 | `tau = T_z/2`, the documented intent |
| `R_S_coeff` | 1.2 | 0.35 | re-fitted against the per-record optimum |
| `R_S_xy_factor` | 0.36 | 1.0 | the anisotropy was a small-sea optimum |
| `MAX_TAU_S` | 3.0 s | 12.0 s | 3.0 s bound at `H_s = 8.5` m |
| `MAX_R_S` | 35 m*s | 400 m*s | 35 clipped the operating point the sea needs |
| accel-bias RW | 1e-3 | 5e-4 | see below |

### The accelerometer-bias prior, and how not to justify it

The bias random walk came down from `1e-3` to `5e-4` per sqrt(s). It is worth
being precise about why, because the way it was found is not the way it is
justified.

It was found by sweeping until the historical accelerometer-bias gate went
green, which is fitting to the test. The justification is that
`5e-4` is exactly the bias random walk the reference simulation generates
(`acc_bias_rw` in `src/util/W3dSimCommon.h`, applied as `sigma*sqrt(dt)`), so
the filter's prior was previously twice the true process noise with nothing
justifying the excess. A correctly specified prior is the right default on its
own terms. On hardware it should come from the instrument's Allan variance
rather than from the simulator, so it is a default to override, not a constant.

The claim that has to survive scrutiny is that the displacement improvement
comes from the wave band and not from this. Decomposed on the four stationary
JONSWAP records (3D RMS in metres; "old band" is the acceleration-band
operating point restored through `W3D_TUNING_BAND=acceleration`, with the
raised clamps in both arms):

| record | old band, old prior | old band, new prior | wave band, old prior | wave band, new prior |
| --- | --- | --- | --- | --- |
| H0.270 | 0.0689 | 0.0734 | 0.0619 | 0.0638 |
| H1.500 | 0.2884 | 0.2957 | 0.2634 | 0.2711 |
| H4.000 | 1.1509 | 1.1436 | 0.7821 | 0.7847 |
| H8.500 | 2.3993 | 2.3731 | 1.4613 | 1.4880 |

The bias prior alone moves 3D RMS by under 1% and in both directions. The wave
band alone carries the entire improvement, and the new prior costs a little of
it back.

### What the accelerometer-bias gate actually measures

The gate the prior was swept against is not measuring bias observability
against the latent OU acceleration. At `H_s = 8.5` m:

| quantity | value |
| --- | --- |
| mean roll error | -1.114 deg |
| apparent specific force from that tilt | -0.1906 m/s^2 |
| mean accelerometer-bias error (x) | +0.1847 m/s^2 |

A 3% match, and 0.5% on a second run of the same record. The bias state is
absorbing a persistent tilt error, which is what a bias state is for, and the
same tilt error is present in the pre-change build (mean roll error -1.136 deg,
roll RMS 1.139 deg new against 1.165 deg old, so the error is almost entirely a
constant offset in both). The DC hypothesis it replaces -- that a longer `tau`
lets the OU process, whose spectrum peaks at DC, steal the bias -- is refuted by
the same measurement: the `a_w` DC offset is 0.028 m/s^2 against a bias error of
0.185, and the two add rather than trade.

So every choice of this prior only decides which state absorbs a tilt error
neither of them causes. Three configurations were measured and none dominates:

| configuration | gates | trade |
| --- | --- | --- |
| `Q_bacc` 5e-4 (shipped) | all 8 pass | yaw +0.4 to +0.5 deg in the small seas |
| `Q_bacc` 1e-3 (historical) | 4 fail | better small-sea yaw and 2% better 3D, but the yaw gate fails at `H_s` 4.0 and 8.5 and the vertical-bias gate fails at 8.5 |
| both bias priors matched to the simulator | 1 fails | best small-sea yaw, vertical-bias gate fails at 8.5 |

The shipped choice is the only one that clears every historical gate, and it is
now justified by the noise model rather than by the sweep that found it. The
1.1-degree static roll bias in steep seas is the real finding here: it is
pre-existing, it is what the bias gate has always been scoring, and it is out of
scope for this change.

One further attribution: the two raised clamps are worth about 13% of 3D RMS at
`H_s = 8.5` m on their own, because `MAX_TAU_S = 3` and `MAX_R_S = 35` were
binding under the old operating point too. They are raised as part of the
wave-band change, since a wave-band `tau` needs the headroom, but the deployed
baseline they are compared against had them binding.

## 1.7 The same defect, and the same fix, in OU-II

OU-II was left on the acceleration-band tracker when OU-III moved off it. That
made the family comparison unfair in a way that flattered neither filter: the
two differ by the extra integral-displacement state, but they were also being
tuned from two different spectral bands, so every OU-III-minus-OU-II difference
mixed the architecture change with the band change.

The defect is the same one described in section 1.2. OU-II's laws are
`tau = c_tau/(2 f_tune)`, `r_p0 = c_p0 sigma_aw tau^2` and
`r_v0 = c_v0 sigma_aw tau`; only the exponent differs from OU-III's `tau^3`, so
an uninformative `f_tune` is squared rather than cubed on the way into the
regularizer. It is the same failure with a smaller multiplier.

`WavePeriodEstimator` is now wired into `SeaStateFusionFilter_OU_II` exactly as
it is in OU-III: fed the leveled vertical acceleration the direction stage
already forms, consulted only once it reports ready, and ablatable back to the
old behaviour with `setWaveBandTuning(false)` / `W3D_TUNING_BAND=acceleration`.
The estimator itself is unchanged and shared, so the `T_z` accuracy table in
section 1.6 applies unaltered, as does the input study beside it -
`setWavePeriodInput()` exists on both filters and the eight-record numbers for
both are in section 1.6.

Coefficients were re-fitted the same way, one factor at a time over the four
stationary JONSWAP records with three seed triplets each, scored on the mean
normalized vertical error with 3D RMS as the secondary check:

| parameter | was | now | why |
| --- | --- | --- | --- |
| `tau_coeff` | 1.5 | 1.0 | `tau = T_z/2`, and the optimum of the scan |
| `R_p0_coeff` | 1.6 | 0.6 | the same law now sees a 2-3x longer `tau` |
| `R_v0_coeff` | 1.4 | 1.1 | shallow optimum; moved with `R_p0_coeff` |
| `R_p0_xy_factor` | 0.31 | 1.0 | the anisotropy was an acceleration-band optimum |
| `sigma_coeff` | 0.85 | 0.85 | scan optimum, unchanged |
| `P_factor` | 1.5 | 1.5 | moves the score by under 1%; left alone |
| accel-bias RW | 1e-3 | 5e-4 | see below |
| `MAX_TAU_S` | 3.0 s | 12.0 s | 3.0 s bound at `H_s = 8.5` m |
| `MAX_R_p0_std` | 18 m | 150 m | 18 bound exactly at `H_s = 8.5` m |
| `MAX_R_v0_std` | 6 m/s | 40 m/s | raised with `r_p0` for the same reason |

Two of these are worth stating separately.

`tau_coeff` and `R_p0_coeff` cannot be fitted one at a time. Along the good
ridge the conserved quantity is `R_p0_coeff * tau_coeff^2`, at about 0.57, which
is what the law says it should be: `r_p0 = c * sigma_aw * tau^2`. The scan was
therefore run as a small 2D grid, and `tau_coeff = 1` was chosen from the ridge
because it is the documented intent, `tau = T_z/2`, and OU-III's value.

The accelerometer-bias random walk had to move with the operating point, and
this was not anticipated. OU-II assumed `1e-3` per sqrt(s) where the reference
simulation generates `5e-4`; OU-III had already been corrected. Left alone it
was not merely inconsistent but actively harmful, because tying `tau` to the
wave band moves the OU corner down toward the bias band and the two states then
compete for the same low-frequency content. With the loose prior the bias state
won that competition and absorbed the persistent tilt-induced specific force.
Across the ten-seed ensemble that showed up as:

| metric, `H_s = 8.5` m | acceleration band | wave band, `1e-3` | wave band, `5e-4` |
| --- | --- | --- | --- |
| pitch RMS [deg] | 0.726 | 0.988 | see study |
| yaw RMS [deg] | 2.93 | 3.63 | see study |
| accel-bias 3D RMS [m/s^2] | 0.151 | 0.232 | see study |

and on the deterministic sentinel the accel-bias gate failed outright at
`H_s = 4.0` and `8.5` m (336% and 325% against a 242% limit, from 241% before).
That gate had 0.6 percentage points of margin in the old build, so it is exactly
the realization-specific sentinel this document argues should eventually be
replaced -- but the underlying regression was real and visible in the ensemble,
not an artifact of where the threshold sits, so the fix is the prior and not the
threshold. The tighter prior costs about 1.3 percentage points of mean
normalized vertical error across the four reference seas. It is chosen because
it is the true process noise and because OU-III uses it, so the family
comparison does not confound the translational state structure with two
different bias priors.

The ten-seed statement is in `reports/results/ou_validation/`.

## 2. Finding 2: the travel-sense classes are a gauge label

### 2.1 Controlled experiment

`tools/rotate_record_heading.py` rewrites a record so the vessel heading is
rotated by a fixed angle about the world vertical: the wave field, the world
specific force and the seeds are untouched, the attitude trajectory becomes
`R_new(t) = Rz(psi) R_old(t)`, body specific force is re-resolved in the new
body frame, and body angular rate is left alone (it is invariant under a
constant world-frame rotation of the attitude). Running the same record at
`psi = 0` and `psi = 180` degrees:

| quantity | heading 0 deg | heading 180 deg |
| --- | --- | --- |
| disp Z RMS | 0.645 m | 0.654 m |
| disp 3D RMS | 2.767 m | 2.706 m |
| axis error vs generator | 4.70 deg | 5.19 deg |
| FORWARD share | 95.5% | 1.8% |
| REVERSE share | 1.8% | 94.6% |

Nothing physical changed and the axis estimate is unaffected, but the reported
class inverts completely. `KalmanWaveDirection` seeds its representative from
boat `+X`, so `FORWARD` means "along whichever end of the axis currently points
forward". The class is a gauge label.

### 2.2 The directed output, however, is correct

The estimator's directed propagation vector, with the vessel heading removed,
scored against the record's physical propagation direction (`A + 180`) by the
new `dir_travel_*` metrics:

| heading | travel error | correct | wrong | unresolved | FORWARD | REVERSE |
| --- | --- | --- | --- | --- | --- | --- |
| 0 deg | +4.56 deg | 95.5% | 1.1% | 3.4% | 95.5% | 1.1% |
| 45 deg | +5.43 deg | 95.5% | 1.2% | 3.3% | 95.5% | 1.2% |
| 180 deg | +4.91 deg | 94.6% | 1.8% | 3.6% | 1.8% | 94.6% |

The correctness rate is flat across headings while the class shares invert. The
propagation sense the estimator resolves is physical and heading-invariant; only
the exported enumeration and the metric built from it were not. The paper
therefore understates the result: 95% is a correctness rate, not just a
commitment rate.

Scoring the 45-degree record also exposed a latent bug in the *axis* metric: it
compared a boat-frame angle against a world-frame azimuth, so it reported a
39.5-degree axis error for an axis that was correct to 5.5 degrees. Every
shipped record has heading 0, which is why it never showed. The vessel heading
is now removed before the comparison.

### 2.3 Record convention

The generator azimuth `A` in the record filename is the direction the waves come
**from**: the propagation-to vector recovered from the truth channels (principal
axis of horizontal displacement, sense from
`sign(<xi_parallel * d(zeta)/dt>) < 0`) lies at `A + 180` in every record. This
was never asserted anywhere, which is why the 180-degree question could not be
settled from the existing metrics.

## 3. Remediation plan

Stages are ordered so each one is separately reviewable and separately
falsifiable. Stages A and B change no filter behaviour.

### Stage A - make the failure visible (no filter behaviour change)

- A1. `tools/rotate_record_heading.py`: the heading-rotation transform used in
  section 2.1, with the angular-rate invariance documented. (done)
- A2. Truth-referenced travel-sense metrics in the simulator harness:
  `dir_travel_error_deg`, `dir_travel_rmse_deg`, `dir_travel_correct_pct`,
  `dir_travel_wrong_pct`, `dir_travel_unresolved_pct`, scored against `A + 180`
  after removing the reference heading. `dir_sense_dominant_pct` is kept so the
  existing series stays readable. (done)
- A3. Assert the record convention (section 2.3) and the invariants of the
  heading transform in `tests/validation/test_record_conventions.py`. (done)
- A4. Fix the axis metric's own heading dependence: `dir_deg_hist` compared a
  boat-frame angle against a world-frame azimuth, which is only valid because
  every shipped record has heading 0. Scoring the 45-degree record gave a
  39.5-degree "error" for an axis that was in fact correct to 5.5 degrees. The
  vessel heading is now removed before the comparison. (done)
- A5. Still to do: adaptation diagnostics behind `W3D_WRITE_ADAPT_LOG` -
  `tau`, `sigma_aw`, `r_S`, accelerometer NIS, `trace(P)`, `trace(P_aw)` and the
  implied `S`-loop corner frequency - plus an offline splitter that separates
  displacement error into drift, wave-band and high-band parts.

### Stage B - price the covariance-synchronisation policy (done)

- B1. `synchronize_aw_covariance_to_stationary_congruent()` implements the
  consistent congruence re-alignment alongside the deployed overwrite, selected
  by `W3D_AW_COV_SYNC=congruent`. The deployed default is unchanged.
- B2. Unit test in `tests/kalman_ou_iii/kalman_ou_common-test.cpp`: the target
  marginal is reached exactly, the whitened cross-covariance is unchanged, the
  non-`a_w` marginals are untouched, and the joint covariance stays PSD, for
  marginal changes of both 16x and 1/16.
- B3. Measured (section 1.1): the consistent version is worse, so the remaining
  work is to retire or bound the periodic re-alignment rather than to fix its
  algebra. That decision belongs with the Stage E re-run, because
  `AdaptiveHeldCovariance` is already the mode that measures it.

### Stage C - estimate the wave period in the wave band (done)

Implemented as described in section 1.6: `WavePeriodEstimator`, wired into the
tuner path, with the clamps and coefficients re-fitted.  `setWaveBandTuning(false)`
restores the acceleration-band behaviour so the change can be ablated.
`tests/kalman_ou_iii/wave_period-test.cpp` pins the estimator against
monochromatic and broadband signals with known `T_z`, against a constant
acceleration offset, and against sample-rate changes.

### Stage D - per-axis regularisation instead of scalar anisotropy constants

Superseded in part.  `R_S_xy_factor` is now 1, so the integral regularisation is
isotropic and no constant imposes an anisotropy; the remaining
`S_factor = 1.87` on the stationary acceleration is still mis-specified relative
to the records (per-axis 0.81 and 0.55 of vertical) but moves 3D RMS by under
2%, so it is left alone rather than changed without an effect to justify it.
Estimating the horizontal stationary covariance per axis, and feeding the full
3x3 through the existing `set_aw_stationary_cov_full`, remains the natural next
step and would make the OU process direction-aware.

### Stage E - re-validation and manuscript (done)

Both committed evidence bundles pin a SHA-256 of the simulator sources they
were produced from, so *any* change under `src/util/` or `tests/*-sim.cpp` -
including the Stage A metric additions - invalidates them by design and
`tests/validation` fails until they are regenerated:

```
python3 tools/ou_robustness.py --mode full --bootstrap-resamples 10000 \
    --output-dir reports/results/ou_robustness
python3 tools/ou_validation.py --mode full
```

That is roughly 310 and 840 twenty-minute simulator replays respectively, so it
is a CI job, not a workstation job. Both live in the `regenerate` job of
`.github/workflows/ou-validation.yml`, which is dispatch-only and runs one job
per bundle: together they do not fit in a single runner lifetime with any
margin, and splitting them means a failure in one does not throw away the
other's hours. The `validate` job that runs on every pull request stays in
smoke mode and never produces numbers anyone cites.

- E1. Re-fit coefficients on the stationary records, then re-run both bundles
  above and regenerate the generated `.tex` parts. (done for OU-III; repeated
  for OU-II with the section 1.7 change)
- E2. Update `doc/kalman_ou_iii/w3d-direction-methods.tex-part` and the results
  caption: report the travel-sense correctness rate and the heading-invariance
  experiment instead of the "not a correctness rate" caveat. (done)
- E3. Update `docs/ou-validation.md` and the covariance-policy discussion to
  match the measured picture in section 1.1. (done)

### Stage F - port the wave-band operating point to OU-II (done)

Section 1.7. Without it the family comparison confounds the extra
integral-displacement state with the spectral band the two filters are tuned
from. `W3D_TUNING_BAND=acceleration` restores the old OU-II behaviour so the
change stays ablatable on both sides of the comparison.

### Still open

- The periodic `P_aw_aw` re-alignment is still the deployed default. Section 1.1
  establishes that it did not cause the loss and that the congruent alternative
  is worse; retiring or bounding it is a separate change with its own evidence.
- The Stage A5 adaptation telemetry (`W3D_WRITE_ADAPT_LOG`: NIS, covariance
  traces, operating-point history, regularizer corner, band-separated
  displacement error) is unimplemented.
- Stage D's per-axis stationary acceleration covariance is unimplemented.
- The ~1.1 degree static roll offset absorbed by the accelerometer-bias state
  predates all of this and is unfixed.
- The executable gates remain calibrated to one realization and
  in at least one case had under one percentage point of margin. They should be
  replaced by separately calibrated ensemble acceptance criteria and kept only
  as deterministic sentinels.

## 4. What this changes about the manuscript's claims

- "Adapting `tau` and `sigma_aw` makes things worse" is reproduced, but the
  cause is not adaptation. It is that the adaptation input is measured in the
  wrong spectral band, so the OU time constant is uninformative and `r_S`,
  which depends on its cube, lands far from what the sea requires.
- The horizontal-versus-vertical asymmetry is a tuning constant
  (`R_S_xy_factor`), not a property of the OU model or of the MEKF.
- The travel-sense result is stronger than reported, not weaker: the directed
  estimate is heading-invariant and correct to about 5 degrees. Only the
  exported class label was gauge-dependent.
