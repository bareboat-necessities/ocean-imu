# Re-gauging the eight-record quality gates

Every simulator that scores the eight reference records — OU-II, OU-III, NLO,
PII observer, TFG — carries regression sentinels fitted by one rule:

> the worst value the filter currently produces across the scored records, plus
> about half a percent, rounded up in the last digit the channel is quoted in.

The rule only works if it is re-applied after the filter changes. A sentinel
that keeps ten points of slack is not catching anything. This is the record of
re-applying it to all five families against the tree that produced build run
3558, and of what had to move.

The rounding used to be "up to the next tenth" for every channel. That is the
right quantum for the percentage channels, whose values sit between 4 and 400,
and the wrong one for yaw in degrees, whose values sit between 1 and 11: a
tenth is 3.5% of OU-III's yaw gate, seven times the margin the rule asks for.
Yaw is now cut to hundredths in all four families that gate it. Nothing about
any filter moved with it.

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
| OU-II | 7 | **four re-derived** for the continuous hard-iron correction, **two cut finer**; all seven now within 0.65% |
| OU-III | 7 | six at the rule and unchanged; **yaw 1.1 → 1.07** |
| NLO | 2 gated (Z only) | both already at the rule; unchanged |
| PII observer | 3 gated (Z, yaw) | Z unchanged; **yaw 10.9 → 10.84** |
| TFG | 7 | **all seven re-derived**; every one came down |

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

### OU-II — four re-derived for the continuous hard-iron correction

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

### OU-III — yaw re-cut, rest unchanged

| gate | limit | worst observed |
| --- | --- | --- |
| Z %Hs JONSWAP | 4.8 | 4.70 (H0.27) |
| Z %Hs PM-Stokes | 4.7 | 4.66 (H0.27) |
| yaw deg | **1.07** | 1.0627 (jonswap H0.27) |
| 3D % JONSWAP | 21.1 | 20.94 (H1.5) |
| 3D % PM-Stokes | 20.9 | 20.72 (H4.0) |
| acc Z bias % | 5.0 | 4.91 (jonswap H8.5) |
| bias 3D % | 98.4 | 97.89 (jonswap H4.0, accel) |

### NLO — unchanged

Yaw is free and ungated here, and the 3D and bias limits are open by design.

| gate | limit | worst observed |
| --- | --- | --- |
| raw Z %Hs JONSWAP | 7.3 | 7.21 (H8.5) |
| raw Z %Hs PM-Stokes | 7.2 | 7.09 (H8.5) |

### PII observer — yaw re-cut, rest unchanged

| gate | limit | worst observed |
| --- | --- | --- |
| Z %Hs JONSWAP | 9.0 | 8.87 (H8.5) |
| Z %Hs PM-Stokes | 9.5 | 9.36 (H8.5) |
| yaw deg | **10.84** | 10.7801 (pmstokes H8.5) |

### TFG — all seven re-derived

| gate | was | now | worst observed |
| --- | --- | --- | --- |
| Z %Hs JONSWAP | 5.5 | **5.3** | 5.21 (H0.27) |
| Z %Hs PM-Stokes | 5.4 | **5.2** | 5.10 (H0.27) |
| yaw deg | 3.3 | **2.94** | 2.9230 (pmstokes H4.0) |
| 3D % JONSWAP | 30.6 | **21.1** | 20.99 (H1.5) |
| 3D % PM-Stokes | 68.0 | **26.0** | 25.78 (H4.0) |
| acc Z bias % | 9.5 | **8.9** | 8.84 (jonswap H4.0) |
| bias 3D % | 415.0 | **400.3** | 398.22 (jonswap H4.0, accel) |

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

Half a percent still leaves 6x on the tightest gate in the set (OU-III yaw) and
better than 10x everywhere else — OU-II's yaw gate, the other one cut against a
filter carrying the solve, sits 0.96% above its binding record whose own drift
is 1.7e-4, a factor of 58 — so the re-cut values are safe. They are not
safe by an enormous factor any more, which is the reason this table exists:
cutting a gate finer than a hundredth of a degree, or below half a percent on
any channel, needs this measurement redone first rather than the old 6e-6
quoted at it.

Reproduce with:

```
make -C tests/<dir> clean
make -C tests/<dir> build CXXFLAGS="-O3 -std=c++20 -Wall -Wextra -Wshadow \
  -Wconversion -funroll-loops -fno-finite-math-only -I$PWD/src \
  -isystem /usr/include/eigen3 -march=x86-64"
```

and compare the `Angles RMS` and `XYZ RMS` lines against a native build.
