# Continuous magnetic hard-iron correction

`Config::mag_continuous_hard_iron` runs a body-fixed magnetometer offset
estimator for the whole life of the filter and feeds its answer back into the
magnetic measurement model. It is on by default in **both OU families** --
`SeaStateFusion_OU_III` and `SeaStateFusion_OU_II` -- and in
`SeaStateFusion_TFG`, on identical settings. The two OU simulators expose it as
`SF_MAG_CONT_HI=0|1`, which is also the matched ablation; the TFG simulator
does not expose it, so scoring TFG against a setting means rebuilding it.

Everything from here to the OU-II section was measured on OU-III, which is
where the mechanism was worked out. None of it is specific to the OU-III
translational model: what the estimator reads is the private Mahony observer's
tilt and the raw magnetometer, and what it writes is the measurement model.
The OU-II numbers are at the end.

Scored on the eight reference records, the deterministic 900 s protocol, it
takes yaw RMS from 2.065 deg mean / 2.133 deg worst to **0.617 deg mean /
0.878 deg worst** on the filter as it now stands. Every record improves. The
tilt channels rise by about two percent as it does -- roll 0.264 -> 0.267 deg,
pitch 0.161 -> 0.164 deg -- and displacement does not move, 0.4108 -> 0.4123 m
of 3-D RMS. `doc/kalman_ou_iii/ins-startup.tex` quotes this table.

The correction landed on a different baseline: 1.835 deg mean / 2.162 deg worst
uncorrected, 0.822 / 1.063 corrected. Its calibration was re-measured against
the filter as it now stands, after the changes that landed under it; the shipped
defaults, the one that moved, and the current results are in
[the re-tune](#the-re-tune). The mechanism sections between here and there are
written against the numbers it landed with and still describe what ships.

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
sea, a fixed ridge of 4e-3 returns 0.27 and 0.88 of the offset, and the
relative ridge returns 0.58 and 0.71.

The floor only does its job while it stays *under* the excitation a working
hull produces. It did not, and [the re-tune](#the-re-tune) below is that
correction.

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

## The re-tune

The calibration above was fitted when the correction was introduced, and the
filter has changed underneath it since -- the reduced physical-MSE integral
regularizer, the `S_factor = 1` anisotropy constant, and the MEKF
sensor-variance sweep, which doubled `sigma_m` and so changed how hard the
magnetometer pulls the heading onto the corrected field. The results table
above no longer describes what ships: measured on this tree with the shipped
`mag_hi_*`, yaw is 0.768 deg mean and 1.265 worst, against the 0.822 / 1.063
recorded when the correction landed. The worst record got worse while the mean
got better.

So the knobs were re-measured against the filter as it now stands. The dominant
finding is **not** a consequence of any of those changes, which is worth saying
plainly: one of the two ridge terms had been mis-scaled from the day it was
written, and what the sweep did was make that visible.

`tools/mag_hard_iron_retune.py` is the sweep. It scores yaw within a draw and
pools only after that normalization, over two seed axes that are not
interchangeable: `--axis init` redraws the magnetometer calibration itself
(`W3D_INIT_SEED`, the thing being estimated), `--axis imu` redraws the sensor
noise with the calibration held (`W3D_IMU_SEED`).

### What was wrong

**The absolute floor had become the whole regularisation.** The ridge is
`model_ridge + model_ridge_relative * mean_eigenvalue(M)`, and the design
intent -- stated in the section above and in the header -- is that the second
term governs and the first only catches a window with no excitation at all.
Measured on the eight records, the excitation is:

| record | `lambda_min` | reported information |
| --- | --- | --- |
| jonswap H0.27 | 3.4e-4 | 4.4 |
| jonswap H1.5 | 1.4e-3 | 18.1 |
| jonswap H4.0 | 2.5e-3 | 32.9 |
| jonswap H8.5 | 4.6e-3 | 58.9 |

(`W3D_MAG_HI_TRACE=1`, end of the replay, effective weight 12935.)

Against that, a floor of 4e-3 is ten times the calmest record's excitation and
comparable to the roughest one's. It was not a floor; it was a fixed ridge
sitting on top of the relative one, and it shrank hardest exactly where the
standing yaw error is largest -- which is the failure mode the whole "why the
ridge scales" argument above exists to avoid.

The consequence is visible in the fitted offsets. Unregularised, over the full
record against a true offset of 1.145 uT, the solve returns 2.73, 2.44, 3.36
and 3.18 uT on the four JONSWAP records -- inflated by a factor of 2.1 to 2.9,
in every sea, roughly constant across them, exactly as the aliasing argument
predicts. The shipped ridge then returned 0.34 of that on the calmest record
and 0.53 on the roughest, when a constant fraction was what the physics asked
for.

### What changed

One number:

    model_ridge:  4.0e-3  ->  5.0e-4

`model_ridge_relative` stays at 0.5, and so does everything else. The new value
is set by the information gate rather than by a sweep: at the saturated weight
of a 600 s exponential window at 25 Hz, `min_information = 2` admits a
direction of eigenvalue 1.3e-4, so a floor of 5e-4 binds only on windows at or
below what the gate already declines. The records barely notice where in that
range it lands: the paired yaw ratio is 0.871 at a floor of zero, 0.889 at
5e-4 and 0.905 at 1e-3, against 1.000 at the shipped 4e-3. A three-point
spread across two decades, under a fifteen-point step out of the old value, is
what a floor that has stopped competing with the relative term looks like.

### Why nothing else moved

The one-factor sweep, five calibration draws by eight records, pooled yaw
change against the shipped configuration:

| knob | values, pooled yaw change |
| --- | --- |
| `model_ridge` | 0: **-10.3%**, 1e-3: -8.4%, 4e-3 shipped, 8e-3: +6.7%, 2e-2: +17.2% |
| `model_ridge_relative` | 0.15: -4.2%, 0.25: -3.3%, 0.5 shipped, 1.0: +5.8%, 1.5: +10.7% |
| `memory_sec` | 150: +13.5%, 300: +1.9%, 1200: -0.5%, cumulative: -0.8% |
| `apply_fraction` | 0.5: +16.8%, 0.75: +6.3%, 1.25: -1.7%, 1.5: -0.1% |
| `slew_tau_sec` | 10: -0.4%, 20: -0.2%, 90: +0.3%, 180: +2.2% |
| `min_information` | 0.5: -0.8%, 8: +12.6% |
| `min_effective_weight`, `max_residual_rms_uT` | inert |

The floor dominates the list. `apply_fraction` and `model_ridge_relative` are
the same lever seen twice -- both scale how much of the fit is applied -- and
both point the same way the floor does, which is what a single over-shrinkage
looks like from three directions. The two gates are inert because these records
carry no heading change and the window saturates long before scoring opens;
they are safety gates, not tuning knobs, and are left alone.

A 4x4 joint grid over the floor and the relative term found the pooled optimum
at `ridge_rel = 0.3`, better than 0.5 by three points of pooled yaw on the
calibration axis. **It did not replicate.** On the IMU-noise axis the same
point is worth -1.4% against the floor cut's -11.8%: it was fitting the five
calibration draws, not the filter. The floor cut holds on both axes, which is
why it is the only thing that ships.

### What it buys

Deterministic protocol, default seeds, the eight scored records, yaw RMS in
degrees:

| Record | Hs | OU-III was | now | OU-II was | now | TFG was | now |
| --- | --- | --- | --- | --- | --- | --- | --- |
| JONSWAP | 0.27 | 1.135 | **0.631** | 1.067 | **0.644** | 1.475 | **0.600** |
| JONSWAP | 1.50 | 1.265 | **0.878** | 1.040 | **0.653** | 1.316 | **0.950** |
| JONSWAP | 4.00 | 0.303 | 0.504 | 0.832 | 1.078 | 0.380 | 0.506 |
| JONSWAP | 8.50 | 0.755 | **0.531** | 0.391 | 0.491 | 1.582 | **1.345** |
| PM--Stokes | 0.27 | 0.705 | **0.695** | 0.690 | 0.712 | 0.926 | **0.427** |
| PM--Stokes | 1.50 | 0.975 | **0.562** | 0.977 | **0.565** | 0.765 | **0.405** |
| PM--Stokes | 4.00 | 0.493 | 0.734 | 0.487 | 0.719 | 0.978 | **0.739** |
| PM--Stokes | 8.50 | 0.515 | **0.400** | 0.528 | **0.429** | 1.472 | **1.247** |
| **mean** | | 0.768 | **0.617** | 0.751 | **0.661** | 1.112 | **0.777** |
| **worst** | | 1.265 | **0.878** | 1.067 | 1.078 | 1.582 | **1.345** |

Mean yaw falls 20% on OU-III, 12% on OU-II and 30% on TFG. No other channel
moves outside the third digit on either OU family: OU-III's vertical goes
4.5214 -> 4.5218 %Hs, 3D 13.6405 -> 13.6518 %, roll 0.3608 -> 0.3616 deg,
pitch 0.1949 -> 0.1957 deg.

Over five magnetometer-calibration draws and all eight records -- 40 paired
replays -- the paired geometric-mean yaw ratio is **0.889**, the worst record
of each draw improves by 10.1% on average, and the largest yaw across all 40
pairs falls from 4.155 to 3.708 deg. Nine of the 40 pairs get worse, the worst
of them by 0.38 deg, on the one draw whose standing error is
misalignment-dominated rather than offset-dominated -- the failure mode "the
limit, stated plainly" describes, which this change does not fix and does not
make qualitatively worse. On the IMU-noise axis, which holds the calibration
and redraws the sensor noise, pooled yaw falls 11.8% and the worst record
16.5%.

### Two records get worse, and it is the same two

`jonswap H4.0` and `pmstokes H4.0` regress on both OU families. They are the
mid-sea records, where the old floor was shrinking least and the fit was
therefore already close to fully applied; handing back the rest of it there
adds more aliased distortion than offset. This is the sea-state dependence the
relative ridge is supposed to remove and does not remove completely, and it is
why `model_ridge_relative` is left where it is rather than cut further.

OU-II's worst record changes identity because of it -- `jonswap H0.27` at 1.067
was binding and is now 0.644, while `jonswap H4.0` becomes the worst at 1.078 --
so OU-II's yaw gate goes *up* by one percent while its mean falls 12%. That is
reported rather than hidden: see the `FAIL_LIMITS` comment in
`kalman_ou_ii-sim.cpp`.

### Measuring the offset directly

The simulators used to report `Bias error RMS (mag, uT)` from the MEKF's
magnetometer-bias state, which does not exist, so the number was the injected
offset itself no matter what the correction did. The three OU simulators now
report the wrapper's applied hard iron in that slot, which makes the offset
error readable in uT: 1.403 uT with the correction off, 1.052 with it on at the
old ridge, 0.784 at the new one on the IMU axis.

Note that this quantity does **not** track yaw, and the reason is worth
stating: the applied offset is a heading-gauge correction, not an offset
estimate. Because the fit absorbs part of the distortion, applying more of it
can move the gauge the right way while moving the offset vector further from
truth. Yaw is the objective; this number is a diagnostic.

## Quality gates

`FAIL_LIMITS` in `kalman_ou_iii-sim.cpp` was re-derived when the correction
landed, by the documented rule of worst observed plus about half a percent
rounded up to the next tenth. **This table is that cut, not the bars in force
now** -- the filter has been re-gauged twice since, and the current values are
in the `FAIL_LIMITS` comment in the simulator (yaw is 0.8827 there, not the
1.068 below):

| gate | before that cut | at that cut | worst observed then |
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

That table is the cut the correction landed with. All three families were
re-gauged again for the re-tune, and the current bars are in each simulator's
`FAIL_LIMITS` comment. OU-III's was re-derived once more when this change met
the shortened-transition `r_S` horizon refit in a merge, which moved its bars
in the fourth digit (yaw 0.8824 to 0.8827); the comment there carries that
pass too. OU-III's yaw bar comes down from 1.27 to 0.8824 and its
two displacement bars come down as well -- though the displacement move is
drift already in the tree from the integral-regularizer schedule rather than
anything this change did, and the comment there separates the two. OU-II's yaw
bar goes up one percent for the reason given in "two records get worse"; TFG's
accelerometer Z-bias bar goes up two percent for the reason this section
already gives about the horizontal accelerometer bias.

## Configuration

| field | default | what it does |
| --- | --- | --- |
| `mag_continuous_hard_iron` | `true` | master switch (`SF_MAG_CONT_HI`) |
| `mag_hi_memory_sec` | 600 | exponential memory of the statistics |
| `mag_hi_model_ridge` | 5e-4 | absolute ridge floor; see [the re-tune](#the-re-tune) |
| `mag_hi_model_ridge_relative` | 0.5 | ridge as a multiple of the mean excitation |
| `mag_hi_min_information` | 2.0 | `weight * lambda_min` the window must reach |
| `mag_hi_min_effective_weight` | 500 | effective samples before any answer |
| `mag_hi_max_residual_rms_uT` | 3.0 | model-fit gate; this is what a turn trips |
| `mag_hi_max_bias_fraction` | 0.35 | plausibility bound against the field norm |
| `mag_hi_apply_fraction` | 1.0 | blunt shrinkage on top of the ridge |
| `mag_hi_slew_tau_sec` | 45 | how fast the applied offset moves |

The simulator reads each as `SF_MAG_HI_*`; `W3D_MAG_HI_TRACE=1` prints the fit,
the applied offset, the information and the residual once a minute to stderr.

The table is identical for all three families that carry the correction: OU-II
and TFG take the same defaults, and a sweep that moves one and not the others
would be comparing calibrations rather than filters. `tools/mag_hard_iron_retune.py`
is the sweep that produced them.

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

(Those are the numbers OU-II landed with; for what it produces on the re-tuned
calibration, see [the re-tune](#the-re-tune).)

Every record improves, mean yaw falls 57%, and the worst record lands at
1.089 deg against OU-III's 1.063 — the two families end up within three
hundredths of a degree of each other, which is what should happen when the
error being removed belongs to the magnetometer rather than to the filter.

On the re-tuned calibration the mechanism is the same and the agreement is
looser: OU-II goes 1.869 -> 0.662 deg mean and 2.134 -> 1.078 worst, against
OU-III's 2.065 -> 0.617 and 2.133 -> 0.878. Both families still lose about two
thirds of their mean yaw to the correction from starting points two tenths of a
degree apart, but the worst records no longer coincide; what separates them is
the residue each family's own attitude solution leaves for the correction to
work against.

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
