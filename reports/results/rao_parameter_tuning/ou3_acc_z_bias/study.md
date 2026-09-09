# OU-III accelerometer Z-bias gate

`make all` exits 2 at two OU-III records: JONSWAP Hs=8.5 accelerometer Z-bias
error 4.39209% and PM-Stokes Hs=8.5 4.46845%, both against the unchanged
4.3% bar. This study screens the estimator coefficients that move that metric
and records why no screened profile closes the gates without paying for it
elsewhere. **No deployed coefficient, executable limit, sensor injection,
scoring window or quality bar is changed by this study.** The two gates remain
open.

## Protocol

All runs use the pinned v1.2.1 vessel records, full 1200 s scored over the
final 900 s, all eight JONSWAP/PM-Stokes records. Screening used the default
draw. Validation uses the eight unseen paired seeds 26317, 27509, 28607,
29717, 17011, 18121, 19333 and 20507, giving 64 paired records and 512 paired
gate cells per profile. Sensor and initialization draws are paired between
profiles. Every screened profile and every failure is retained in
`configs.json` and `paired_gates.csv.gz`.

Per-axis accelerometer-bias driving noise is not an override the committed
simulator exposes, so the screen applied `per-axis-bias-source.patch`, which
adds `OU_III_ACC_BIAS_RW_X/Y/Z` in the same shape the OU-II vessel adapter
already uses. That patch is study-only and is not committed to the simulator.

## What moves the metric

The scored quantity is the RMS error of the Z bias estimate over the window,
divided by the largest true Z bias in that window. It is insensitive to the
integral regularizer, the horizontal anisotropy factor, the OU sigma and tau
coefficients, the accelerometer noise floor and the moment-averaging horizon:
each of those moves it by less than .3%. Only three coefficients move it.

| Coefficient | Direction that helps | Default-seed effect |
|---|---|---:|
| Initial accelerometer-bias std | .004 -> .0453, the std of the injected uniform draw | -0.45% |
| Bias correlation time | 5000 s -> 2e5 s, toward the injected pure random walk | -1.0% |
| Horizontal bias driving noise | 5e-4 -> 5e-5, pricing the tilt/bias confound | -1.7% |
| Vertical bias driving noise | 5e-4 -> 4e-4, below the injected 5e-4 | -0.8% |

The first two are model corrections: the harness draws the initial bias
uniformly on +/-8e-3 g, whose standard deviation is .0453 m/s^2 rather than
the .004 the filter assumes, and it drives the bias as a pure random walk with
no mean reversion. Both help, and together they reach 4.330 / 4.393 — still
above the bar. Closing the remaining 2.2% requires the last two rows, which
are not model corrections: they deliberately drive the filter's bias noise
below the value the harness injects.

## Paired validation

Mean paired ratio against the committed profile over 64 fresh records; below
1 is better. `viol` counts individual gate violations out of 512 paired cells.

| Profile | viol | z | 3D | yaw | roll | pitch | accz | acc3d | gyro3d |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| committed | 158 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| matched prior | 156 | 1.0000 | .9997 | 1.0014 | 1.0043 | .9217 | .9986 | .9571 | .9941 |
| prior + long tau | 158 | .9956 | .9984 | 1.0048 | 1.0360 | .9237 | .9889 | .9557 | .9945 |
| horizontal suppression | 168 | 1.0009 | 1.0213 | 1.0038 | .9580 | 1.0142 | .9969 | 1.0341 | 1.0024 |
| suppression + long tau | 173 | 1.0109 | 1.0244 | 1.0022 | .9523 | 1.0138 | .9964 | 1.0342 | 1.0023 |
| both + matched prior | 172 | 1.0109 | 1.0241 | 1.0017 | .9215 | .9284 | .9953 | .9769 | .9966 |
| the same, true vertical noise | 168 | .9964 | 1.0190 | 1.0011 | .9221 | .9320 | .9832 | .9803 | .9967 |
| clears every default gate | 177 | 1.0108 | 1.0249 | .9948 | .9290 | .9721 | .9957 | 1.0093 | .9987 |

Only the two model corrections are paired improvements. The matched prior is
the single clean result: 156 violations against 158, pitch 7.8% better,
accelerometer 3D bias 4.3% better, and nothing worse by more than half a
percent. Adding the long correlation time keeps that and improves the scored
metric by 1.1%, at 3.6% worse roll.

## Why the gates stay open

Exactly one screened profile clears all eight default-seed gates with every
bar unchanged: horizontal driving noise 5e-5, vertical 4e-4, correlation time
2e5 s, initial std .025. It reaches 4.29673 and 4.29670 against the 4.3 bar,
with pitch .9715 -> .7265, roll .8969 -> .5815 and accelerometer 3D bias
.9175 -> .5953 of their bars.

It is not committed, because the paired evidence says the pass is
realization-specific rather than an improvement in the scored quantity. Over
the eight fresh seeds that profile changes the Z-bias metric by -0.4%, which
is inside its own scatter, and leaves its violation count at 22 of 64 exactly
as the committed profile does. What it does change is low-wave vertical
accuracy, 1.1% worse with 7 new violations where the committed profile has
none, and 3D displacement, 2.5% worse on 62 of 64 records. Both regressions
track the vertical driving-noise reduction, which is also what buys the last
half percent on the gate.

The metric's own scatter is the reason a default-seed pass proves little: over
the eight fresh seeds the same committed profile scores between 2.51% and
14.58% on these two records, so the 4.3 bar, cut from one realization, is
violated on 22 of 64 fresh cells before and after every profile screened here.

## Conclusion

The two gates are open, and closing them at the default seed requires driving
the filter's accelerometer-bias noise below the value the harness injects,
which costs low-wave vertical and 3D accuracy on unseen seeds. Two model
corrections — an initial bias prior matched to the injected draw and a
correlation time matched to its pure random walk — are real paired
improvements and remain available, but they leave the worst records at 4.393
and 4.330. Neither the bar nor the profile is moved by this study.
