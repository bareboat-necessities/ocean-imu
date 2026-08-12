# Re-gauging the eight-record quality gates

Every simulator that scores the eight reference records — OU-II, OU-III, NLO,
PII observer, TFG — carries regression sentinels fitted by one rule:

> the worst value the filter currently produces across the scored records, plus
> about half a percent, rounded up to the next tenth.

The rule only works if it is re-applied after the filter changes. A sentinel
that keeps ten points of slack is not catching anything. This is the record of
re-applying it to all five families against the tree that produced build run
3558, and of what had to move.

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

## Result: four of five were already at the rule

| family | gates | verdict |
| --- | --- | --- |
| OU-II | 7 | all seven already at the rule; unchanged |
| OU-III | 7 | all seven already at the rule; unchanged |
| NLO | 2 gated (Z only) | both already at the rule; unchanged |
| PII observer | 3 gated (Z, yaw) | all three already at the rule; unchanged |
| TFG | 7 | **all seven re-derived**; every one came down |

OU-III's were re-derived with the continuous hard-iron correction and OU-II's
with the parity work, both recent; NLO's and PII's have not been disturbed
since the 900 s window landed. TFG's were fitted before the adaptation-policy
change that brought its orchestrator up to OU-III's, and were not revisited
afterwards.

## Measured worst values, per family

Vertical is `Z RMS` as a percentage of `H_s`; 3D is displacement RMS as a
percentage of max `|disp_ref|_3D`; bias figures are RMS error as a percentage
of the maximum true bias in the window. Bracketed record is where the worst
value occurs.

### OU-II — unchanged

| gate | limit | worst observed |
| --- | --- | --- |
| Z %Hs JONSWAP | 6.9 | 6.86 (H0.27) |
| Z %Hs PM-Stokes | 6.9 | 6.81 (H0.27) |
| yaw deg | 2.2 | 2.16 (pmstokes H1.5) |
| 3D % JONSWAP | 21.1 | 20.91 (H1.5) |
| 3D % PM-Stokes | 21.2 | 21.02 (H8.5) |
| acc Z bias % | 5.4 | 5.33 (pmstokes H8.5) |
| bias 3D % | 92.2 | 91.65 (jonswap H4.0, accel) |

### OU-III — unchanged

| gate | limit | worst observed |
| --- | --- | --- |
| Z %Hs JONSWAP | 4.8 | 4.70 (H0.27) |
| Z %Hs PM-Stokes | 4.7 | 4.66 (H0.27) |
| yaw deg | 1.1 | 1.06 (jonswap H0.27) |
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

### PII observer — unchanged

| gate | limit | worst observed |
| --- | --- | --- |
| Z %Hs JONSWAP | 9.0 | 8.87 (H8.5) |
| Z %Hs PM-Stokes | 9.5 | 9.36 (H8.5) |
| yaw deg | 10.9 | 10.78 (pmstokes H8.5) |

### TFG — all seven re-derived

| gate | was | now | worst observed |
| --- | --- | --- | --- |
| Z %Hs JONSWAP | 5.5 | **5.3** | 5.21 (H0.27) |
| Z %Hs PM-Stokes | 5.4 | **5.2** | 5.10 (H0.27) |
| yaw deg | 3.3 | **3.0** | 2.92 (pmstokes H4.0) |
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
- **Yaw.** 2.92 deg against OU-III's 1.06. That gap is the continuous magnetic
  hard-iron correction, which is an OU-III component TFG does not carry; see
  `docs/ou-iii-continuous-hard-iron.md` for what it removes and why the error
  it removes is a gauge rather than a tracking error.

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
directory. The metrics are deterministic: `-march=native`, `x86-64` and
`x86-64-v2` agree to within 6e-6 relative, which is why the sentinels can sit
half a percent above what the filter produces without failing every run.
