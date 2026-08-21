# The startup gravity gate, and why it never closed under way

`SeaStateFusion_OU_II` and `SeaStateFusion_OU_III` do not start their MEKF
until a gate certifies that the startup observer is levelled. That gate is
what decides two things: when the magnetometer may begin averaging the world
reference, and when the finished attitude may be handed to the MEKF, which is
the moment `isLive()` becomes true and the device starts reporting heading,
heave and wave direction.

Measured on the eight reference records, time to `isLive()` was:

| record | before | after |
|---|---|---|
| jonswap H0.27 | 22.0 | 22.0 |
| jonswap H1.5 | 52.1 | 22.0 |
| jonswap H4.0 | 109.1 | 22.0 |
| jonswap H8.5 | 107.1 | 23.8 |
| pmstokes H0.27 | 22.0 | 22.0 |
| pmstokes H1.5 | 72.2 | 22.0 |
| pmstokes H4.0 | 76.4 | 22.5 |
| pmstokes H8.5 | 150.0 | 32.5 |

150 s is the handoff timeout. On that record the gate never closed at all and
startup ended by giving up on quality rather than by reaching it.

Both families are affected identically; the numbers above are OU-III and OU-II
agrees to the sample.

## What the gate was measuring

The residual was `||u_a x u_g||`, the sine of the angle between the low-passed
body-frame specific force and the direction the attitude predicts gravity
should lie in, with a 1 s low pass on the accelerometer and a threshold of
0.075 (4.3 deg) that had to hold for 2 s.

At anchor that is exactly right. Under way it is not a levelling measurement at
all. The measured specific force is gravity plus the orbital acceleration of
the wave, and the low pass that is supposed to remove the second term runs in
the *body* frame -- the frame that is rolling and pitching through the wave
period. The orbital term is zero mean in the world frame, not in a frame that
is itself oscillating with it, so what the average leaves behind is a
wave-correlated angle rather than nothing.

Measured on jonswap H8.5, the body-frame residual runs between 0.03 and 0.45
for the entire record and spends most of it above the 0.075 threshold. The
gate could only close when a couple of quiet seconds happened to line up, which
is why the times above look random rather than ordered by sea state: pmstokes
H4.0 (76 s) beat jonswap H4.0 (109 s) not because it was easier but because it
got luckier.

Meanwhile the observer's actual tilt error on that record is about 0.85 deg by
40 s (`docs/ou-iii-startup-init.md` has the measurement against truth). The
gate was rejecting an attitude that was already five times better than its own
threshold.

## What it measures now

The specific force is rotated into the estimator's own world frame first, and
the average is taken there:

```
a_w   = q_bw * a_body                       accWorldFromBody()
a_lp  = lowpass(a_w, tau)                   tau = 12 s
sin   = |horizontal(a_lp)| / |a_lp|         gravityAlignResidualSinWorld()
```

Orbital acceleration *is* zero mean in that frame, so a long enough average
leaves gravity, and the angle between it and world down is the estimator's tilt
error and nothing else. On jonswap H8.5 that residual settles below 0.05 within
about twenty seconds and stays there, so the gate closes on tilt quality, when
tilt is actually good, in every sea state.

Only the tilt part of `q_bw` matters -- both the residual and the branch test
read the world z axis and the horizontal magnitude, and both are invariant
under `q_bw -> Rz(psi) q_bw` -- so the observer's arbitrary, drifting heading
cannot leak into the gate.

The branch test comes with it. A sine reads the same at an angle and at its
supplement, so on its own it accepts an attitude flipped through 180 degrees;
in the world frame the branch is simply the sign of the down component, since
a specific force at rest points up. As before, the branch is a proof condition
rather than a threshold and has no knob, and the handoff timeout is held to it
so a forced handoff is delayed rather than seeding the filter upside down.

### The two knobs

| setting | default | what it does |
|---|---|---|
| `mag_gravity_align_world_tau_sec` | 12 | horizon of the world-frame average |
| `mag_gravity_align_world_warmup_sec` | 5 | how long that average must run before its verdict counts |

The horizon has to span whole wave periods for the orbital term to cancel, so
it is set against the longest swell the device is expected to start up in
rather than against the sea it is in. The frequency tracker is not converged
this early and being conservative costs settling time rather than accuracy, so
it is a constant: 12 s covers the band these filters work in. It barely
matters. Swept over 8 s and 16 s against warmups of 6 s and 12 s, every scored
metric on every record moves in the fourth digit -- the change is structural,
not fitted.

The warmup exists because the average and the observer are seeded from the
*same* accelerometer sample. Until the average has moved off that seed, a small
residual only says the two agree about the instant they both started from,
which they do by construction even when the boat was mid-wave and both are
wrong. It is set so the gate's earliest possible verdict lands with the
magnetometer's first eligible sample -- `mag_delay_sec` less the 2 s hold the
gate serves anyway -- which is the last moment at which it is free. Past that
it delays a calm start one second for one: the magnetometer cannot begin
averaging before `mag_delay_sec` however early the gate closes. At 5 s the two
calm records come out bit-identical to the old gate.

## What it costs

Five of the eight records are unchanged to four significant figures and the two
calm ones are bit-identical, because there the old gate already closed on the
first quiet stretch. The big-sea records re-mix, and not all in the same
direction.

OU-III, default seed:

| record | roll | pitch | acc 3D bias | 3D % |
|---|---|---|---|---|
| jonswap H4.0 | 0.3607 -> 0.3148 | 0.1327 -> 0.1252 | 74.64 -> 65.13 | 12.171 -> 12.156 |
| jonswap H8.5 | 0.3102 -> 0.1678 | 0.1517 -> 0.1122 | 65.36 -> 34.73 | 13.626 -> **13.919** |
| pmstokes H8.5 | 0.2719 -> 0.1619 | 0.1835 -> 0.1129 | 62.95 -> 31.37 | 14.599 -> 14.431 |

OU-II, default seed:

| record | roll | pitch | acc 3D bias | 3D % |
|---|---|---|---|---|
| jonswap H4.0 | 0.4345 -> 0.3876 | 0.2122 -> 0.2193 | 91.83 -> 83.77 | 15.051 -> 14.994 |
| jonswap H8.5 | 0.3665 -> 0.1933 | 0.2804 -> **0.3280** | 81.71 -> 64.19 | 16.711 -> **17.092** |
| pmstokes H8.5 | 0.2774 -> 0.1667 | 0.1772 -> 0.1276 | 62.42 -> 34.49 | 18.218 -> 17.850 |

jonswap H8.5 is the record that pays. In OU-II its pitch is 12 percent worse
and that is systematic rather than a re-draw: paired over six IMU seeds the
ratio is between 1.087 and 1.155 on every one of them. In OU-III its 3D
displacement is 1.0 percent worse pooled over the same six seeds, while its
sister record pmstokes H8.5 is 0.8 percent better.

Pooled over all eight records and six IMU seeds, geometric mean of new/base
with a 95% interval:

| channel | OU-III | OU-II |
|---|---|---|
| Z %Hs | 0.9977 [0.9946, 1.0008] | **0.9966 [0.9940, 0.9992]** |
| 3D % | 0.9998 [0.9972, 1.0023] | 0.9994 [0.9958, 1.0030] |
| roll | 1.0463 [0.9430, 1.1611] | 1.0328 [0.9468, 1.1266] |
| pitch | 0.9789 [0.9503, 1.0083] | 1.0112 [0.9900, 1.0327] |
| yaw | 0.9980 [0.9786, 1.0178] | 0.9853 [0.9696, 1.0012] |
| acc Z bias | 0.9973 [0.9914, 1.0032] | 0.9989 [0.9925, 1.0054] |
| acc 3D bias | 1.0375 [0.9427, 1.1418] | 1.0184 [0.9582, 1.0825] |
| gyro 3D bias | **0.9650 [0.9410, 0.9896]** | **0.9628 [0.9374, 0.9889]** |

Nothing degrades at 95 percent in either family. Gyro bias improves in both,
vertical displacement improves in OU-II, and everything else is flat. The
per-record spread the ratio columns show is dominated by seed scatter: roll on
these records spans 0.04 to 0.72 deg across seeds under *either* filter, so a
ratio taken against a small draw says more about the draw than about the
filter. Pooled RMS across seeds is the honest number, and by it roll moves
+9% on jonswap H8.5, +4% on pmstokes H8.5, -9% on jonswap H4.0 and 0% on the
rest.

So this is a trade, not a free win: one record loses a channel, three gain
more, and the aggregate is neutral-to-better. It is taken because the thing
being bought is not a fourth digit. A device that reports nothing for two and
a half minutes after power-on -- and, on the record where that happened, was
never going to report on quality at all -- is a worse instrument than one whose
worst-case attitude on a single 8.5 m sea is a tenth of a degree softer.

## TFG has the same construction and is not changed here

`SeaStateFusionFilter_TFG::updateProxyGravityQuality_()` is the same gate --
its own comment says so, "that is OU-III's construction" -- against the same
body-frame low pass, so the defect is there too. It is left alone: the change
would be mechanical but its regression bars are a separate deterministic
protocol that has to be re-derived and argued on its own numbers, and the
OU-II/OU-III re-gauge above is what this change was measured against. That is
a scope decision rather than a judgement that TFG is fine.

## Gates

Both families' `FAIL_LIMITS` blocks were re-derived with
`tools/ou_regauge_gates.py` under the usual deterministic protocol. Seven of
ten bars come down in OU-III (the gyro-bias aggregate by 16 percent) and seven
of nine in OU-II (the accelerometer-bias aggregate by 9 percent). The bars that
go up are jonswap H8.5's, named above, and the per-record and per-seed
measurements behind each are recorded in the comment blocks next to the limits.
