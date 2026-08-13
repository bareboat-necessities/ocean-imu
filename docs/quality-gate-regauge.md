# Re-gauging the eight-record quality gates

Every simulator that scores the eight reference records — OU-II, OU-III, NLO,
PII observer, TFG — carries regression sentinels fitted by one rule:

> the worst value the filter currently produces across the scored records, plus
> about half a percent, rounded up in the last digit the channel is quoted in.

The rule only works if it is re-applied after the filter changes. A sentinel
that keeps ten points of slack is not catching anything. This is the record of
re-applying it to all five families against the tree that produced build run
3558, and of what had to move.

The rounding used to be "up to the next tenth" for every channel of every
family. One quantum for values spanning 1 to 400 cannot deliver one margin: a
tenth is 3.5% of OU-III's yaw gate, 2% of its vertical gate and 0.03% of its
bias gate. The small-valued channels were therefore carrying two to seven times
the margin the rule asks for, and the large ones none of it.

Every gate is now written to whatever precision delivers about half a percent —
a thousandth where the value is near 1, a hundredth for single digits, a tenth
where a tenth is already fine enough. **All 23 gates across the five families
now sit between 0.50% and 0.65% above what the filter produces**, and none of
that rounding change involved a filter change.

## What was measured

The deterministic protocol each simulator gates on: default seeds, the final
900 s of each 1200 s replay, the eight JONSWAP and PM-Stokes records at
`H_s` 0.27, 1.50, 4.00 and 8.50 m. The three short record families (cnoidal,
fenton, gerstner) are shorter than the scoring window and are reported
`QUALITY_GATE: SKIPPED` by every simulator, so they are not part of this.

Runs are `W3D_COLLECT_ALL_GATES=1` so a breach scores the remaining records
instead of exiting at the first one.

    ./kalman_ou_iii-sim --input wave_data_jonswap_H0.270_....csv     # per record
    ./nlo-sim                                                        # whole directory

## Result

| family | gates | verdict |
| --- | --- | --- |
| OU-II | 7 | four re-derived for the continuous hard-iron correction, three cut finer |
| OU-III | 7 | six cut finer; bias 3D already at the rule |
| NLO | 2 gated (Z only) | both cut finer |
| PII observer | 3 gated (Z, yaw) | all three cut finer |
| TFG | 7 | all seven re-derived, then five cut finer |

OU-II's four moves are a filter change, not a rounding change: it now carries
the continuous hard-iron correction, on OU-III's settings, and pays OU-III's
price for it. The rest of the moves in the table are the change of quantum
described above, with no filter behind them.

OU-III's gates were re-derived when it took the correction, NLO's and PII's
have not been disturbed since the 900 s window landed, and TFG's were fitted
before the adaptation-policy change that brought its orchestrator up to
OU-III's and were not revisited afterwards.

## Measured worst values, per family

Vertical is `Z RMS` as a percentage of `H_s`; 3D is displacement RMS as a
percentage of max `|disp_ref|_3D`; bias figures are RMS error as a percentage
of the maximum true bias in the window. Bracketed record is where the worst
value occurs.

### OU-II — four re-derived for the continuous hard-iron correction, three cut finer

Limits below are for the shipped filter, which now runs the correction. The
"was" column is the same filter without it, which is the `SF_MAG_CONT_HI=0`
ablation and exceeds the yaw gate by a factor of two, as it should.

| gate | was | now | worst observed | margin |
| --- | --- | --- | --- | --- |
| Z %Hs JONSWAP | 6.9 | 6.9 | 6.8638 (H0.27) | 0.53% |
| Z %Hs PM-Stokes | 6.9 | **6.85** | 6.8061 (H0.27) | 0.65% |
| yaw deg | 2.18 | **1.095** | 1.0895 (jonswap H1.5) | 0.50% |
| 3D % JONSWAP | 21.1 | 21.1 | 20.9870 (H1.5) | 0.54% |
| 3D % PM-Stokes | 21.2 | **21.3** | 21.1940 (H8.5) | 0.50% |
| acc Z bias % | 5.4 | **5.44** | 5.4059 (jonswap H8.5) | 0.63% |
| bias 3D % | 92.2 | **94.4** | 93.8996 (jonswap H4.0, accel) | 0.53% |

Two of those moves are the filter and two are the quantum. Yaw, 3D PM-Stokes,
acc Z and bias 3D were re-derived because the correction moved what the filter
produces. Z PM-Stokes and acc Z were then also cut finer, from a tenth to a
hundredth, because a tenth is 1.5% of a 6.8 and 1.9% of a 5.4 — rounding a
half-percent margin up to a tenth was handing back three times the rule. Every
OU-II gate now sits between 0.50% and 0.65% above what the filter produces, and
the whole set passes on an `-march=x86-64` rebuild as well as a native one, with
0.49% of headroom at the tightest point.

Yaw halves and three go up, which is the trade OU-III recorded when it took the
same change: the correction walks the heading onto the corrected field during
the run and the horizontal accelerometer bias — the least observable quantity
scored, its error already above 90% of the true bias — absorbs part of that
motion. Displacement does not move (vertical mean 6.454 → 6.461 %Hs, 3D mean
18.90 → 19.03 %) and pitch improves, 0.289 → 0.255 deg.

### OU-III — six cut finer, then all seven re-derived for `S_factor = 1`

The first pass, at the then-shipped `S_factor = 1.87`:

| gate | was | then | worst observed | margin |
| --- | --- | --- | --- | --- |
| Z %Hs JONSWAP | 4.8 | **4.72** | 4.6952 (H0.27) | 0.53% |
| Z %Hs PM-Stokes | 4.7 | **4.69** | 4.6600 (H0.27) | 0.64% |
| yaw deg | 1.1 | **1.068** | 1.0627 (jonswap H0.27) | 0.50% |
| 3D % JONSWAP | 21.1 | **21.05** | 20.9361 (H1.5) | 0.55% |
| 3D % PM-Stokes | 20.9 | **20.83** | 20.7197 (H4.0) | 0.53% |
| acc Z bias % | 5.0 | **4.93** | 4.9054 (jonswap H8.5) | 0.50% |
| bias 3D % | 98.4 | 98.4 | 97.8908 (jonswap H4.0, accel) | 0.52% |

Then re-derived once more when the horizontal stationary acceleration scale
went from 1.87 to the records' own value of 1
([`ou-iii-anisotropy-consistency.md`](ou-iii-anisotropy-consistency.md)), which
moved every gated quantity:

| gate | then | now | worst observed | margin |
| --- | --- | --- | --- | --- |
| Z %Hs JONSWAP | 4.72 | **4.74** | 4.7106 (H0.27) | 0.63% |
| Z %Hs PM-Stokes | 4.69 | 4.69 | 4.6580 (H0.27) | 0.69% |
| yaw deg | 1.068 | **1.297** | 1.2896 (jonswap H1.5) | 0.57% |
| 3D % JONSWAP | 21.05 | **20.95** | 20.8367 (H1.5) | 0.54% |
| 3D % PM-Stokes | 20.83 | **20.86** | 20.7483 (H4.0) | 0.54% |
| acc Z bias % | 4.93 | **4.63** | 4.6004 (jonswap H8.5) | 0.64% |
| bias 3D % | 98.4 | **81.84** | 81.4268 (pmstokes H4.0, accel) | 0.51% |

Five come down, three of them materially: JONSWAP 3D, accelerometer Z bias, and
the 3D bias limit, which falls by 17 points because the isotropic prior stops
the horizontal bias absorbing so much of the sea. Two go up with small-sea
losses that the sweep priced at 0.1 to 0.2% of 3D RMS.

The yaw sentinel moving up 21% needs its own note, since a loosened sentinel
that hides a regression is worse than no sentinel. It is not hiding one. Yaw on
the binding record spans 1.05 to 6.57 deg over five IMU seeds under the *old*
constant, so the default-seed value this protocol scores is one draw from a wide
distribution, not a measure of yaw quality. Paired over those seeds and all
eight records, the new constant lowers yaw RMS by 3.2% pooled and improves four
of five seeds on the binding record itself; the deployed draw is one of the few
that moves the other way. The rule is applied to the protocol as written, and
the quality claim rests on the seeds — `reports/results/ou_anisotropy` carries
both.

### NLO — both cut finer

Yaw is free and ungated here, and the 3D and bias limits are open by design.

| gate | was | now | worst observed | margin |
| --- | --- | --- | --- | --- |
| raw Z %Hs JONSWAP | 7.3 | **7.26** | 7.2143 (H8.5) | 0.63% |
| raw Z %Hs PM-Stokes | 7.2 | **7.13** | 7.0865 (H8.5) | 0.61% |

### PII observer — all three cut finer

| gate | was | now | worst observed | margin |
| --- | --- | --- | --- | --- |
| Z %Hs JONSWAP | 9.0 | **8.91** | 8.8651 (H8.5) | 0.51% |
| Z %Hs PM-Stokes | 9.5 | **9.41** | 9.3622 (H8.5) | 0.51% |
| yaw deg | 10.9 | **10.84** | 10.7801 (pmstokes H8.5) | 0.56% |

### TFG — all seven re-derived, then five cut finer

| gate | was | now | worst observed | margin |
| --- | --- | --- | --- | --- |
| Z %Hs JONSWAP | 5.5 | **5.24** | 5.2090 (H0.27) | 0.60% |
| Z %Hs PM-Stokes | 5.4 | **5.13** | 5.1004 (H0.27) | 0.58% |
| yaw deg | 3.3 | **2.938** | 2.9230 (pmstokes H4.0) | 0.51% |
| 3D % JONSWAP | 30.6 | **21.1** | 20.9914 (H1.5) | 0.52% |
| 3D % PM-Stokes | 68.0 | **25.91** | 25.7764 (H4.0) | 0.52% |
| acc Z bias % | 9.5 | **8.89** | 8.8360 (jonswap H4.0) | 0.61% |
| bias 3D % | 415.0 | **400.3** | 398.2190 (jonswap H4.0, accel) | 0.52% |

Nothing in the simulator moved these; the filter did, in the adaptation-policy
commit that fixed the tuner's commit timing, moved the schedule onto the 0.1 s
tick, lowered the `r_S` floor from 0.4 to 0.15 and made the `S=0` cadence
self-similar in `tau`. The generated results table in
`doc/kalman_tfg/tfg-sim-results-generated.tex-part` already carried the new
numbers; only the gates were left behind.

The horizontal channel is where that shows. 3D error on PM-Stokes was 67.63%
when the bars were last fitted and is 25.78% now, and on JONSWAP it has landed
exactly on OU-III's own bar of 21.1. The "3D error degrades sharply on the
large-wave records" note in the simulator was true of the filter it was written
against and is not true of this one.

Two things did not improve and are still recorded rather than endorsed:

- **Accelerometer bias.** 398% of the true bias at worst, against OU-III's 98%,
  and above 100% on six of the eight records. An error larger than the quantity
  being estimated is not an estimate. Cause still unestablished.
- **Yaw.** 2.92 deg against OU-III's 1.06 and OU-II's 1.09. That gap is the
  continuous magnetic hard-iron correction, which both OU families carry and
  this one does not; see `docs/continuous-mag-hard-iron.md` for what it removes
  and why the error it removes is a gauge rather than a tracking error.

`bias_3d_percent` gates the gyro channel as well as the accelerometer, and the
accelerometer sets it on every family. TFG's worst gyro value is 125.60%
(pmstokes H4.0) against a bar of 400.3, so a gyro-bias regression has to more
than triple before that bar sees it. Splitting the two means changing
`W3dFailureLimits` and re-deriving for every family that uses it; not done here.

## Reproducing

```
make -C tests/<dir> build
cd tests/<dir> && W3D_COLLECT_ALL_GATES=1 W3D_WRITE_TIMESERIES=0 ./<sim>
```

with the eight `wave_data_{jonswap,pmstokes}_*.csv` records in the working
directory.

For OU-III the rule itself is mechanised, which removes the arithmetic from the
hand-application above:

```
python3 tools/ou_iii_regauge_gates.py
```

It runs the eight records under this protocol, reports the worst value and its
record per gate, applies the half-percent rule at the quantum the channel is
quoted in, flags any shipped gate the filter no longer clears, and prints a
`FAIL_LIMITS` body to paste. It reproduces all seven shipped limits from the
filter that produced them, which is the check that it implements the rule
rather than a rule. `--env OU_III_S_FACTOR=1.87` and friends re-gauge an
ablation without rebuilding.

## How tight is safe: the determinism budget, measured

The half-percent margin only works because the metrics are reproducible, and
the simulator comments used to put that reproducibility at 6e-6 relative across
`-march` levels. That number is now stale for the two filters that carry matrix
solves. Rebuilding each simulator at `-march=x86-64` instead of the host's
native `cascadelake` and rescoring all eight records:

| family | worst yaw drift | worst drift, any gated channel |
| --- | --- | --- |
| OU-III | **8.3e-4** (jonswap H8.5) | 8.3e-4 — yaw is the worst channel |
| TFG | 2.6e-4 (pmstokes H8.5) | 8.8e-4 (gyro bias 3D, pmstokes H4.0) |
| OU-II | 6.0e-4 (jonswap H8.5) | 7.3e-4 (acc Z bias, pmstokes H8.5) |
| PII observer | 0 — bit-identical | 1.5e-6 (Z, pmstokes H1.5) |
| NLO | not gated | 1.5e-6 (raw Z, pmstokes H4.0) |

OU-II's row is measured with the continuous hard-iron correction on, and it
roughly tripled when the correction landed — 2.4e-4 yaw and 4.4e-4 worst-channel
before it. That is the same mechanism as OU-III's row and is the strongest
evidence for the cause given below: the drift follows the solve, not the
family.

So 6e-6 still describes NLO and the PII observer, and understates the OU and
TFG families by two orders of magnitude. The likely cause is visible in the
code rather than inferred: OU-III's continuous hard-iron solve inverts a normal
matrix of order 1e-3, which multiplies the last bits of a `double` accumulation
by up to a thousand on the way into the applied offset, and the applied offset
walks the heading — so yaw is exactly where it should show up, and does.

That is the number every gate in this document is checked against, one by one,
rather than against the half-percent rule alone. The thinnest margin-to-drift
ratios in the whole set are OU-III's bias 3D at 15x, TFG's 3D PM-Stokes at 14x
and TFG's yaw at 25x; everything else is between 40x and several thousand. The
check that matters is not arithmetic: every family was rescored with an
`-march=x86-64` build against the tightened gates, and all five pass with at
least 0.48% of headroom at their tightest point. They are not
safe by an enormous factor any more, which is the reason this table exists:
cutting any of these further, or below half a percent on any channel, needs this
measurement redone first rather than the old 6e-6 quoted at it. The three
ratios named above are where that bites first.

Reproduce with:

```
make -C tests/<dir> clean
make -C tests/<dir> build CXXFLAGS="-O3 -std=c++20 -Wall -Wextra -Wshadow \
  -Wconversion -funroll-loops -fno-finite-math-only -I$PWD/src \
  -isystem /usr/include/eigen3 -march=x86-64"
```

and compare the `Angles RMS` and `XYZ RMS` lines against a native build.
