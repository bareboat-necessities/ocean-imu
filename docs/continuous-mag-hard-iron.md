# Continuous magnetic hard-iron correction

`Config::mag_continuous_hard_iron` runs a body-fixed magnetometer offset
estimator for the whole life of the filter and feeds its answer back into the
magnetic measurement model. It is on by default in **both OU families** --
`SeaStateFusion_OU_III` and `SeaStateFusion_OU_II` -- on identical settings,
and each simulator exposes it as `SF_MAG_CONT_HI=0|1`, which is also the
matched ablation.

Everything from here to the OU-II section was measured on OU-III, which is
where the mechanism was worked out. None of it is specific to the OU-III
translational model: what the estimator reads is the private Mahony observer's
tilt and the raw magnetometer, and what it writes is the measurement model.
The OU-II numbers are at the end.

Scored on the eight reference records, the deterministic 900 s protocol, it
takes yaw RMS from 1.835 deg mean / 2.162 deg worst to **0.822 deg mean /
1.063 deg worst**. Every record improves. Pitch improves, roll and the
displacement channels do not move.

## What the standing yaw error actually was

The OU-III yaw error is not a tracking error. On `pmstokes H1.5` the scored
window has mean 2.146 deg and standard deviation 0.265 deg: the filter tracks
heading well and points 2.1 deg away from magnetic north the whole time.

That offset is a **gauge**. The startup acquisition averages the magnetometer
in a yaw-stripped tilt frame and declares the horizontal direction of that
average to be north:

    ref = mean(R_tilt m) = B_hull + mean(R_tilt) b + (distortion),

for a hull at fixed heading, where `B_hull` is the true field in the hull's own
level frame and `b` is the body-fixed offset. The filter's north is the
direction of `ref`; true north is the direction of `B_hull`; the yaw error is
the angle between them, which for the simulator's offset of about 1.15 uT
against a 20.7 uT horizontal field is a couple of degrees. Nothing downstream
can recover it, because the world frame *is* that direction.

Measured directly, by removing each part of the simulator's magnetometer error
with simulator truth and rescoring (branch-only oracle, not deployable):

| oracle correction | yaw mean | yaw worst |
| --- | --- | --- |
| none (baseline) | 1.835 | 2.162 |
| soft iron and misalignment only | 1.811 | 2.121 |
| hard iron only | 0.496 | 1.084 |
| both | 0.479 | 1.059 |

So on the shipped seed essentially all of the standing error is the hard-iron
offset, and removing it perfectly would leave 0.5 deg. That is the headroom the
estimator is chasing.

## Why a startup window cannot get it

`MagAutoTuner` already solves this identification problem, and its
`estimate_hard_iron` flag has always been off by default. Its own tests say
why: at one attitude an offset and a rotated world field are the same reading,
and only tilt breaks the tie. The normal matrix that decides the tie is, to
second order in roll and pitch standard deviations,

    I - mean(R)^T mean(R)  ~  diag(sigma_p^2, sigma_r^2, sigma_r^2 + sigma_p^2),

which is of order 1e-3 for a hull working through a few degrees. Fifteen
seconds of that is not enough, and the fix is not a longer startup window --
it is not closing the window at all.

## The estimator

`src/tuner/ContinuousMagHardIronEstimator.h` keeps exponentially weighted
sufficient statistics of `(R_i, m_i)` and solves

    (I - Abar^T Abar) b = mbar - Abar^T wbar,      Abar = mean(R_i),

periodically, eliminating the world field analytically so it never becomes a
filter state. Inputs are the raw magnetometer and the yaw-stripped tilt
quaternion of the private Mahony observer. Both are exogenous: no MEKF state is
read, no loop is closed through the filter, and the ISS argument, the state
dimension, `F`, `Q`, `P` and the gain computation are untouched.

Stripping yaw from the accumulation frame matters twice. It makes the frame a
pure function of gravity, so the observer's unobservable, drifting heading
cannot leak in. And it ties the frame to the hull's own heading, so a hull that
turns makes the world field move inside the frame and the model residual rises
immediately -- which is what the residual gate is for.

### Why the solve is regularised, and why the ridge scales

Inverting a matrix of order 1e-3 multiplies everything the model does *not*
carry by up to a thousand. What it does not carry is soft iron and sensor
misalignment, and those are modulated by the very tilt the offset is read from.
The consequence is not noise and does not average away: fitted over the full
1200 s against truth attitude, the unregularised solve returns about 1.8x the
true offset on every well-excited record, and applying that overshoots the
correction into a yaw error of the opposite sign.

The size of that aliased error is roughly `eps*|B|` — the distortion times the
field — *independently of the sea state*, because it enters the right-hand side
through the same tilt that carries the signal. A hull working through thirty
degrees is not fitting a cleaner offset than one working through five; it is
fitting the same wrong one with more confidence. A fixed ridge therefore
shrinks hardest exactly where the fit is most trustworthy, which is backwards,
and the sweep shows it: at a fixed ridge the big seas over-apply and regress
while the calm seas are barely corrected.

`Config::model_ridge_relative` scales the ridge with the mean eigenvalue of the
normal matrix, so the fraction returned is a property of the sensor rather than
of the weather. `Config::model_ridge` remains as an absolute floor for a window
with almost no excitation at all. `continuous_mag_hard_iron-test` pins the
difference: replaying one offset and one distortion through a small and a large
sea, a fixed ridge alone returns 0.27 and 0.88 of the offset, and the relative
ridge returns 0.24 and 0.66.

## Why the reference cannot simply follow the offset

The obvious way to keep the correction self-consistent is to subtract `b` from
every sample and subtract the matching `mean(R) b` from the magnetic reference.
That is exactly wrong, and it took a measurement to see it: it is a **no-op**.
The innovation

    m - b - R^T (ref - mean(R) b)

is identical to the uncorrected one at the attitude the filter already holds,
because `R^T mean(R) ~ I`. Forcing a 10 uT offset through that path moves the
scored yaw by 0.3 deg. The error being removed is a gauge, and a change that
preserves the innovation preserves the gauge.

What the filter does instead:

- the offset is subtracted from the magnetometer stream, and
- the reference stays in `MagAutoTuner`'s canonical form -- horizontal
  magnitude on `+X`, vertical below it -- with only its magnitude and vertical
  component moved, and only by the amount the offset changes them.

The corrected field then implies a different heading than the reference asserts,
and the magnetometer update walks the yaw onto it over its own time constant.
No attitude state is written: the correction stays a change of measurement-model
parameters, which is what keeps it out of the stability argument.

The magnitude and dip are moved by a *delta* against the value the startup
acquisition gated, not adopted wholesale from the estimator's own window. That
window is longer and less selective, and simply adopting its magnitude and dip
costs a fifth of the roll accuracy on the moderate seas while the offset
correction itself costs none of it. With the delta form, an estimator that
never validates leaves all 56 scored metrics bit-identical.

Nothing starts until the two-stage startup acquisition has finished, so first
heading, handoff and the refinement are exactly as they were. Accumulation,
however, begins with the first magnetometer sample, so by the time the
correction is armed there is already a window to solve.

## Results

Deterministic 900 s protocol, one realization per record, default seeds. RMS in
degrees; 3D is displacement RMS in metres.

| Record | Hs (m) | yaw off | yaw on | roll off | roll on | pitch off | pitch on | 3D off | 3D on |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JONSWAP | 0.27 | 2.075 | 1.063 | 0.279 | 0.281 | 0.174 | 0.175 | 0.0586 | 0.0588 |
| JONSWAP | 1.50 | 1.825 | 1.050 | 0.256 | 0.257 | 0.169 | 0.152 | 0.2689 | 0.2688 |
| JONSWAP | 4.00 | 1.298 | 1.041 | 0.481 | 0.491 | 0.358 | 0.330 | 0.7902 | 0.7941 |
| JONSWAP | 8.50 | 1.104 | 0.709 | 0.439 | 0.447 | 0.469 | 0.445 | 1.5260 | 1.5169 |
| PM--Stokes | 0.27 | 2.155 | 0.741 | 0.265 | 0.266 | 0.200 | 0.200 | 0.0569 | 0.0573 |
| PM--Stokes | 1.50 | 2.162 | 1.026 | 0.260 | 0.259 | 0.223 | 0.208 | 0.2549 | 0.2550 |
| PM--Stokes | 4.00 | 2.053 | 0.473 | 0.405 | 0.401 | 0.278 | 0.235 | 0.7095 | 0.7199 |
| PM--Stokes | 8.50 | 2.007 | 0.476 | 0.328 | 0.317 | 0.259 | 0.196 | 1.3043 | 1.2887 |
| **mean** | | **1.835** | **0.822** | 0.339 | 0.340 | 0.266 | 0.243 | 0.6212 | 0.6199 |
| **worst** | | **2.162** | **1.063** | 0.481 | 0.491 | 0.469 | 0.445 | 1.5260 | 1.5169 |

Across all 56 scored metrics the geometric mean of the ratio is 0.877 and the
single worst regression is +4.8%, on the Y displacement of `pmstokes H4.0`.

## The limit, stated plainly

The estimator corrects one of the two things that shift the gauge. It cannot
correct the other, and there is no local signal that says which one it is
facing.

Repeating the study over five draws of the magnetometer calibration
(`W3D_INIT_SEED`, which redraws both the offset and the distortion), 40
record/seed pairs:

| | yaw mean | yaw worst | pairs with worse yaw |
| --- | --- | --- | --- |
| off | 2.430 | 3.827 | -- |
| on | 1.726 | 4.510 | 5 / 40 |

All five regressions are one draw, `seed=23`, and the oracle says why: on that
draw the standing error is misalignment-dominated, not offset-dominated —
removing its hard iron perfectly takes yaw from 2.98 to 2.74, while removing
its soft iron and misalignment perfectly takes it to 1.02. A hard-iron
estimator has nothing useful to do there, and the fraction of the distortion it
aliases into its answer moves the gauge the wrong way, by up to 28%.

That case is not detectable from inside. The distortion's aliased part is, by
construction, the part the model *explains*, so it does not raise the residual:
the estimator reports 1.683 uT on the shipped draw and 1.685 uT on `seed=23`,
against a white-noise floor of 1.386 uT. The residual gate catches a turning
hull, not a misaligned sensor.

Two things follow. The shrinkage is not a tuning convenience, it is the bound
on that failure mode, and it is why the ridge is left large enough to give up
real accuracy on the good draws. And the way to actually fix `seed=23` is to
identify the distortion, which needs heading excitation the records do not
contain — a hull that turns. The estimator stands itself down on a turn today
rather than exploiting it; using the turn is the obvious next piece of work and
is deliberately not attempted here.

## Quality gates

`FAIL_LIMITS` in `kalman_ou_iii-sim.cpp` is re-derived, by the documented rule
of worst observed plus about half a percent rounded up to the next tenth:

| gate | was | now | worst observed |
| --- | --- | --- | --- |
| yaw deg | 2.2 | **1.068** | 1.0627 (jonswap H0.27) |
| 3D % PM--Stokes | 20.6 | **20.9** | 20.72 (pmstokes H4.0) |
| acc Z bias % | 5.6 | **5.0** | 4.91 (jonswap H8.5) |
| acc 3D bias % | 95.8 | **98.4** | 97.89 (jonswap H4.0) |
| Z %Hs JONSWAP | 4.8 | 4.8 | 4.70 (jonswap H0.27) |
| Z %Hs PM--Stokes | 4.7 | 4.7 | 4.66 (pmstokes H0.27) |
| 3D % JONSWAP | 21.1 | 21.1 | 20.94 (jonswap H1.5) |

Two go up. The correction walks the heading onto the corrected field during the
run, and the horizontal accelerometer bias absorbs some of that motion. That
quantity is the least observable one scored here — its error exceeds the true
bias under every configuration this filter has shipped, which is what a figure
near 100% means — so a 2.6 point move on it is not a meaningful loss of
accelerometer-bias accuracy, but the sentinel has to admit it rather than hide
it.

Running the `SF_MAG_CONT_HI=0` ablation now exceeds the yaw gate, as it should:
that limit is fitted to the filter that ships. Score the ablation with
`W3D_COLLECT_ALL_GATES=1`.

## Configuration

| field | default | what it does |
| --- | --- | --- |
| `mag_continuous_hard_iron` | `true` | master switch (`SF_MAG_CONT_HI`) |
| `mag_hi_memory_sec` | 600 | exponential memory of the statistics |
| `mag_hi_model_ridge` | 4e-3 | absolute ridge floor |
| `mag_hi_model_ridge_relative` | 0.5 | ridge as a multiple of the mean excitation |
| `mag_hi_min_information` | 2.0 | `weight * lambda_min` the window must reach |
| `mag_hi_min_effective_weight` | 500 | effective samples before any answer |
| `mag_hi_max_residual_rms_uT` | 3.0 | model-fit gate; this is what a turn trips |
| `mag_hi_max_bias_fraction` | 0.35 | plausibility bound against the field norm |
| `mag_hi_apply_fraction` | 1.0 | blunt shrinkage on top of the ridge |
| `mag_hi_slew_tau_sec` | 45 | how fast the applied offset moves |

The simulator reads each as `SF_MAG_HI_*`; `W3D_MAG_HI_TRACE=1` prints the fit,
the applied offset, the information and the residual once a minute to stderr.

The table is identical for both families: OU-II takes the same defaults, and a
sweep that moves one and not the other would be comparing two calibrations
rather than two filters.

## OU-II

OU-II carries the correction on the same settings, through the same
`ContinuousMagHardIronEstimator`, wired at the same three points of
`updateMag()`: accumulate ahead of every startup gate, apply after the
second-stage refinement has landed, subtract from the stream the MEKF sees.
The port is mechanical because nothing in the mechanism touches the
translational model — the estimator reads the private Mahony observer's
yaw-stripped tilt and the raw magnetometer, and writes measurement-model
parameters.

Deterministic 900 s protocol, one realization per record, default seeds:

| Record | Hs (m) | yaw off | yaw on | pitch off | pitch on | Z%Hs off | Z%Hs on |
| --- | --- | --- | --- | --- | --- | --- | --- |
| JONSWAP | 0.27 | 2.069 | 1.062 | 0.193 | 0.195 | 6.860 | 6.864 |
| JONSWAP | 1.50 | 1.858 | 1.089 | 0.200 | 0.180 | 6.563 | 6.583 |
| JONSWAP | 4.00 | 1.401 | 0.937 | 0.370 | 0.319 | 6.351 | 6.408 |
| JONSWAP | 8.50 | 1.327 | 0.583 | 0.418 | 0.362 | 6.263 | 6.270 |
| PM--Stokes | 0.27 | 2.134 | 0.726 | 0.208 | 0.209 | 6.806 | 6.806 |
| PM--Stokes | 1.50 | 2.161 | 1.033 | 0.248 | 0.230 | 6.322 | 6.318 |
| PM--Stokes | 4.00 | 2.086 | 0.504 | 0.341 | 0.288 | 6.190 | 6.186 |
| PM--Stokes | 8.50 | 2.060 | 0.569 | 0.333 | 0.253 | 6.277 | 6.256 |
| **mean** | | **1.887** | **0.813** | 0.289 | 0.255 | 6.454 | 6.461 |
| **worst** | | **2.161** | **1.089** | 0.418 | 0.362 | 6.860 | 6.864 |

Every record improves, mean yaw falls 57%, and the worst record lands at
1.089 deg against OU-III's 1.063 — the two families end up within three
hundredths of a degree of each other, which is what should happen when the
error being removed belongs to the magnetometer rather than to the filter.

The price is the same one OU-III paid, in the same place:

| | off | on |
| --- | --- | --- |
| acc Z bias, worst % of true | 5.33 | 5.41 |
| acc 3D bias, worst % of true | 91.66 | 93.90 |
| 3D displacement, worst % | 21.02 | 21.19 |
| vertical, mean %Hs | 6.454 | 6.461 |
| roll, mean deg | 0.333 | 0.332 |

The horizontal accelerometer bias absorbs part of the heading motion, and it is
the least observable quantity in the set — an error above 90% of the true bias
means the error is larger than the thing being estimated, which was true before
this change and is true after it. Displacement does not move. Four of OU-II's
seven gates were re-derived accordingly; see `docs/quality-gate-regauge.md`.

With `SF_MAG_CONT_HI=0` the filter reproduces its pre-correction self to within
2.6e-4 relative, which is the same order as rebuilding it at a different
`-march`. The ablation is not bit-identical because adding the estimator's
members to the wrapper changes what the optimizer does with the rest; it is
identical in behaviour.

The limit stated for OU-III applies unchanged here, because it is a property of
the identification problem and not of either filter: a calibration draw whose
standing error is misalignment-dominated rather than offset-dominated has
nothing for this estimator to win, and the fraction of the distortion it aliases
moves the gauge the wrong way. The residual gate cannot see that case.
