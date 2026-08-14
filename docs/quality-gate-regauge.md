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
now sit between 0.50% and 0.69% above what the filter produces**, and none of
that rounding change involved a filter change. (The band was 0.50–0.65% when
this was written; OU-III's re-derivation for `S_factor = 1`, below, put its
PM-Stokes vertical gate at 0.69% because a hundredth is the finest quantum that
channel is quoted in.)

That 0.69% is the residue of quoting the quantum in absolute terms at all. The
two OU families have since gone to a **relative** quantum — four significant
figures, so the last digit is always a thousandth of the value — which removes
the choice of decimal from the margin and lands all of their gates between
0.50% and 0.57%; see "The OU families" below, which also adds the first bars
this repository has ever put on roll and pitch. The other three families still
carry the absolute quantum and the 0.50–0.69% band above.

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
| OU-II | 7 → **9** | four re-derived for the continuous hard-iron correction, three cut finer; then four cut finer again on the relative quantum, and roll and pitch gated for the first time |
| OU-III | 7 → **9** | six cut finer; then all seven re-derived for `S_factor = 1`; then three cut finer again on the relative quantum, and roll and pitch gated for the first time |
| NLO | 2 gated (Z only) | both cut finer |
| PII observer | 3 gated (Z, yaw) | all three cut finer |
| TFG | 7 | all seven re-derived, then five cut finer; then re-derived again after the OU parity work |

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

Four of these were cut again, and two bars added, in "The OU families" below;
this section is the state before that pass.

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

Three of these were cut again, and two bars added, in "The OU families" below;
this section is the state before that pass.

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

### The OU families — a relative quantum, and two new attitude bars

Neither family's filter moved for this pass. Both sets of gates got stricter
anyway, in two independent ways.

#### The quantum goes relative

"Rounded up in the last digit the channel is quoted in" only delivers one
margin if that digit is a fixed *fraction* of the value. It was a fixed
absolute step — a hundredth for both families — which is 0.2% of OU-III's 4.7
vertical gate and 0.01% of OU-II's 94 bias gate, so the decimal point was
setting the margin instead of the rule. Quoting every channel to four
significant figures makes the quantum a thousandth of the value everywhere.
Seven gates move; the other seven were already on a four-figure value and are
unchanged:

| family | gate | was | now | worst observed | margin |
| --- | --- | --- | --- | --- | --- |
| OU-II | Z %Hs JONSWAP | 6.9 | **6.899** | 6.8644 (H0.27) | 0.52% → 0.50% |
| OU-II | Z %Hs PM-Stokes | 6.85 | **6.841** | 6.8062 (H0.27) | 0.64% → 0.51% |
| OU-II | acc Z bias % | 5.44 | **5.435** | 5.4073 (jonswap H8.5) | 0.61% → 0.51% |
| OU-II | bias 3D % | 94.4 | **94.37** | 93.8979 (jonswap H4.0, accel) | 0.53% → 0.50% |
| OU-III | Z %Hs JONSWAP | 4.74 | **4.735** | 4.7106 (H0.27) | 0.63% → 0.52% |
| OU-III | Z %Hs PM-Stokes | 4.69 | **4.682** | 4.6580 (H0.27) | 0.69% → 0.52% |
| OU-III | acc Z bias % | 4.63 | **4.624** | 4.6004 (jonswap H8.5) | 0.64% → 0.51% |

The 0.69% that this document had to explain away when it was written is gone,
and the widest margin left in either family is OU-III's yaw at 0.57%.

#### Roll and pitch are gated

Both simulators have printed roll and pitch RMS on every record since they were
written, and neither had a bar on either. That is the wrong shape for what
these filters have actually been doing: the Mahony-proxy startup policy, the
two-stage magnetic reference, the continuous hard-iron correction and OU-III's
`S_factor = 1` are attitude work, and yaw was the only attitude channel
carrying a sentinel. Mean pitch improved under both of the changes whose
ablations are still wired up — OU-II 0.2987 → 0.2545 deg for the startup
policy and 0.2890 → 0.2545 for the hard-iron correction, OU-III 0.2230 and
0.2200 → 0.1891 — with nothing watching it.

| family | gate | new bar | worst observed | margin |
| --- | --- | --- | --- | --- |
| OU-II | roll deg | **0.4778** | 0.4753 (jonswap H4.0) | 0.52% |
| OU-II | pitch deg | **0.3639** | 0.3620 (jonswap H8.5) | 0.52% |
| OU-III | roll deg | **0.42** | 0.4179 (jonswap H4.0) | 0.50% |
| OU-III | pitch deg | **0.2211** | 0.2200 (pmstokes H4.0) | 0.50% |

Fitted by the same rule on the same protocol, and they discriminate. Both
matched ablations — `SF_MAG_CONT_HI=0` and `W3D_STARTUP_INIT=staged_mekf`,
which are the filters these defaults replaced — breach the new pitch bar:
OU-III on four of the eight records under either ablation, OU-II on two. The
roll bars are the weaker of the pair, and only OU-III's discriminates at all:
`staged_mekf` breaches it on two records, `SF_MAG_CONT_HI=0` on one, and
neither ablation moves OU-II's roll enough to reach its bar.

| family | arm | worst roll | worst pitch | mean pitch |
| --- | --- | --- | --- | --- |
| OU-II | deployed | 0.4753 | 0.3620 | 0.2545 |
| OU-II | `SF_MAG_CONT_HI=0` | 0.4711 | 0.4176 | 0.2890 |
| OU-II | `staged_mekf` | 0.3998 | 0.4756 | 0.2987 |
| OU-III | deployed | 0.4179 | 0.2200 | 0.1891 |
| OU-III | `SF_MAG_CONT_HI=0` | 0.4215 | 0.2656 | 0.2200 |
| OU-III | `staged_mekf` | 0.5398 | 0.2626 | 0.2230 |

Those breaches are the point, in the same way the yaw breach under
`SF_MAG_CONT_HI=0` is: a bar fitted to the filter that ships should fail the
filter it replaced. Score the ablations with `W3D_COLLECT_ALL_GATES=1`.

Both bars are gated on the magnetometer-on protocol only, like yaw, because
`--nomag` is a different filter and these bars were not fitted to it. OU-II
loses the most: worst-case roll 0.4753 → 0.5382 and worst-case pitch
0.3620 → 0.7113 deg, both past its bars. OU-III's worst-case pitch goes
0.2200 → 0.2595, also past its bar, while its worst-case roll *improves* by
19% — the magnetic reference is not uniformly good for tilt, which is worth
knowing and is not something a gate should be deciding.

#### OU-II, re-derived again after the `(r_p0, r_v0)` re-fit

The two passes above were cuts against a filter that stood still. This one is
the other case: OU-II's two pseudo-measurement variance coefficients were
re-fitted — `R_p0_coeff` 0.6 → 0.65, `R_v0_coeff` 1.1 → 1.3, see
`docs/ou-ii-pseudo-variance-tuning.md` — and all nine bars followed the filter.

Every one of the nine still passed on its previous value, so nothing here was
forced. What forced the pass is that the previous bars had been cut to the rule
against the previous filter, so once the filter moved, five sat too loose and
four too tight — pitch at 0.0001 deg of margin, which is under this family's own
rebuild drift on that channel and therefore a bar a rebuild decides.

| gate | was | now | worst observed, before → after | margin |
| --- | --- | --- | --- | --- |
| Z %Hs JONSWAP | 6.899 | **6.865** | 6.8644 → 6.8300 (H0.27) | 0.51% |
| Z %Hs PM-Stokes | 6.841 | **6.848** | 6.8062 → 6.8139 (H0.27) | 0.50% |
| yaw deg | 1.095 | **1.089** | 1.0895 → 1.0833 (jonswap H1.5) | 0.52% |
| roll deg | 0.4778 | **0.4792** | 0.4753 → 0.4768 (jonswap H4.0) | 0.50% |
| pitch deg | 0.3639 | **0.3657** | 0.3620 → 0.3638 (jonswap H8.5) | 0.51% |
| 3D % JONSWAP | 21.1 | **20.92** | 20.9867 → 20.8140 (H1.5) | 0.51% |
| 3D % PM-Stokes | 21.3 | **21.03** | 21.1935 → 20.9203 (H8.5) | 0.52% |
| acc Z bias % | 5.435 | **5.324** | 5.4073 → 5.2969 (jonswap H8.5) | 0.51% |
| bias 3D % | 94.37 | **94.47** | 93.8979 → 93.9911 (jonswap H4.0, accel) | 0.51% |

Reverting the two coefficients through `OU_R_P0_COEFF=0.6 OU_R_V0_COEFF=1.1`
puts all nine back at the rule to the digit, which is the control that says
this set moved for the re-fit and for nothing else.

The margins land in the same 0.50–0.52% band as before, so the
margin-to-drift ratios of the determinism budget below are unchanged to within
their own precision — pitch stays the thinnest in the document, at 9.1x.

#### TFG, re-derived again for `S_factor` 1.20 → 1.00

The adaptation coefficients were re-swept on a multi-seed instrument
(`docs/tfg-adaptation-refit.md`); four of the five held and the horizontal
stationary acceleration prior went to the value the records measure. Six bars
move by under a percent. The seventh is the first outright breach this document
has had to record:

| gate | was | now | worst observed, before → after | margin |
| --- | --- | --- | --- | --- |
| Z %Hs JONSWAP | 4.803 | **4.807** | 4.7784 → 4.7830 (H0.27) | 0.50% |
| Z %Hs PM-Stokes | 4.709 | **4.707** | 4.6846 → 4.6833 (H0.27) | 0.51% |
| yaw deg | 1.536 | **1.59** | 1.5278 (pmstokes H8.5) → 1.5818 (jonswap H8.5) | 0.52% |
| 3D % JONSWAP | 21.14 | **21.13** | 21.0299 → 21.0169 (H1.5) | 0.54% |
| 3D % PM-Stokes | 20.71 | **20.74** | 20.6045 → 20.6322 (H4.0) | 0.52% |
| acc Z bias % | 5.026 | **5.022** | 5.0002 → 4.9961 (pmstokes H8.5) | 0.52% |
| bias 3D % | 167.6 | **164.5** | 166.688 → 163.607 (pmstokes H4.0, accel) | 0.55% |

The yaw bar goes up 3.5% and the filter no longer clears the old one, so it
gets the treatment OU-III's yaw sentinel got for the same reason. The binding
record moves to jonswap H8.5, whose yaw spans 1.03 to 4.56 deg across six IMU
seeds; four of the six improve under the new constant, the six-seed mean on
that record falls 2.685 → 2.476 deg, and pooled over all eight records yaw is
0.9578 of the shipped filter. The deterministic seed draws the smallest of the
six under the old constant, which is the whole of the 54% jump the sentinel
sees. The bar follows the protocol it is written against; the quality claim
rests on the seeds.

That is now twice that a wide-distribution yaw record has moved a deterministic
sentinel against the direction of the ensemble — OU-III's on jonswap H1.5, this
one on jonswap H8.5. Both were caught by pairing across seeds after the fact.
A yaw sentinel scored on one realization is the weakest bar in this document,
and the two OU families' figures should be read with that in mind.

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

#### TFG, re-derived again after the OU parity work

The two items recorded above as "did not improve" were the accelerometer bias
and yaw, and yaw was attributed to a specific missing feature. That feature —
along with the rest of the OU families' magnetic acquisition and the band-passed
sigma channel — is now carried by TFG, and the tuner coefficients were re-fitted
afterwards. See `docs/tfg-design.md` sections 8 and 9. All seven bars were
re-derived by the same rule:

| gate | was | now | worst observed | margin |
| --- | --- | --- | --- | --- |
| Z %Hs JONSWAP | 5.24 | **4.803** | 4.7784 (H0.27) | 0.51% |
| Z %Hs PM-Stokes | 5.13 | **4.709** | 4.6846 (H0.27) | 0.52% |
| yaw deg | 2.938 | **1.536** | 1.5278 (pmstokes H8.5) | 0.54% |
| 3D % JONSWAP | 21.1 | **21.14** | 21.0298 (H1.5) | 0.52% |
| 3D % PM-Stokes | 25.91 | **20.71** | 20.6045 (H4.0) | 0.51% |
| acc Z bias % | 8.89 | **5.026** | 5.0002 (pmstokes H8.5) | 0.52% |
| bias 3D % | 400.3 | **167.6** | 166.688 (pmstokes H4.0, accel) | 0.55% |

The JONSWAP 3D bar is the one that moves *up*, by 0.04, and it is worth saying
why rather than hiding it in the table: worst-case 3D fell from 25.78 to 21.03,
but the record that now sets the bar (jonswap H1.5) sits 0.19 above where the
old binding record happened to land, and the old 21.1 was already only 0.5%
above it. Every other bar comes down, three of them by more than a third.

Where that leaves the two standing findings:

- **Yaw** is no longer attributable to a missing feature. 1.53 deg against
  OU-III's 1.068 and OU-II's 1.09 — still behind, but by a factor of 1.4 rather
  than 2.8, and with the same correction now on both sides of the comparison.
  Ablating the continuous hard-iron correction alone puts it back to 3.30.
- **Accelerometer bias** is still the outlier and is still recorded rather than
  endorsed: 167% against OU-III's 98.4, above 100% of the true bias on two of
  the eight records rather than six. Most of what came out was the standing
  attitude error the old one-sample magnetic reference left behind, which the
  bias state was absorbing because the two are only weakly separable in waves.
  What remains has not been attributed.

`bias_3d_percent` still gates the gyro channel as well. TFG's worst gyro value
is now 46.58% (jonswap H8.5) against a bar of 167.6, so the accelerometer still
sets it and a gyro-bias regression still has to more than triple before that bar
sees it.

Both an `-march=native` and an `-march=x86-64` build pass all seven; the binding
records move by up to 3.5e-4 relative between them, so the thinnest margin is 14
times the spread it has to survive.

## Reproducing

```
make -C tests/<dir> build
cd tests/<dir> && W3D_COLLECT_ALL_GATES=1 W3D_WRITE_TIMESERIES=0 ./<sim>
```

with the eight `wave_data_{jonswap,pmstokes}_*.csv` records in the working
directory.

For the two OU families and TFG the rule itself is mechanised, which removes the
arithmetic from the hand-application above:

```
python3 tools/ou_regauge_gates.py --family ou_iii
python3 tools/ou_regauge_gates.py --family ou_ii
python3 tools/ou_regauge_gates.py --family tfg
```

TFG carries seven bars rather than nine — roll and pitch are gated for the OU
families and not for it — and the script skips whatever a family does not
carry.

It runs the eight records under this protocol, reports the worst value and its
record per gate, applies the half-percent rule at the quantum the channel is
quoted in, flags any shipped gate the filter no longer clears, and prints a
`FAIL_LIMITS` body to paste. It reproduces every shipped limit of all three
families from the filter that produced them, which is the check that it
implements the rule rather than a rule. `--env OU_III_S_FACTOR=1.87` and
friends re-gauge an ablation without rebuilding, and `--json` dumps every
channel of every record for the build-drift comparison below.

## How tight is safe: the determinism budget, measured

The half-percent margin only works because the metrics are reproducible, and
the simulator comments used to put that reproducibility at 6e-6 relative across
`-march` levels. That number is now stale for the two filters that carry matrix
solves. Rebuilding each simulator at `-march=x86-64` instead of the host's
native `cascadelake` and rescoring all eight records:

| family | worst yaw drift | worst drift, any gated channel |
| --- | --- | --- |
| OU-III | **8.3e-4** (jonswap H8.5) | 8.3e-4 — yaw is the worst channel |
| TFG | 3.5e-4 (jonswap H4.0) | 4.6e-4 (gyro bias 3D, jonswap H8.5) |
| OU-II | 6.0e-4 (jonswap H8.5) | 7.3e-4 (acc Z bias, pmstokes H8.5) |
| PII observer | 0 — bit-identical | 1.5e-6 (Z, pmstokes H1.5) |
| NLO | not gated | 1.5e-6 (raw Z, pmstokes H4.0) |

OU-II's row is measured with the continuous hard-iron correction on, and it
roughly tripled when the correction landed — 2.4e-4 yaw and 4.4e-4 worst-channel
before it. That is the same mechanism as OU-III's row and is the strongest
evidence for the cause given below: the drift follows the solve, not the
family. TFG's row was re-measured after it took the same correction; its yaw
drift is 3.5e-4 against 2.6e-4 before, which is the same mechanism again and a
smaller step because TFG already carried an `LDLT` per measurement update.

So 6e-6 still describes NLO and the PII observer, and understates the OU and
TFG families by two orders of magnitude. The likely cause is visible in the
code rather than inferred: OU-III's continuous hard-iron solve inverts a normal
matrix of order 1e-3, which multiplies the last bits of a `double` accumulation
by up to a thousand on the way into the applied offset, and the applied offset
walks the heading — so yaw is exactly where it should show up, and does.

That is the number every gate in this document is checked against, one by one,
rather than against the half-percent rule alone. The thinnest margin-to-drift
ratios in the whole set are OU-III's bias 3D at 15x and TFG's yaw at 15x;
everything else is between 40x and several thousand. TFG's 3D PM-Stokes, at 14x
when its gates were last cut, is now 90x — the channel improved by 20% while its
`-march` drift did not move. The
check that matters is not arithmetic: every family was rescored with an
`-march=x86-64` build against the tightened gates, and all five pass with at
least 0.48% of headroom at their tightest point. They are not
safe by an enormous factor any more, which is the reason this table exists:
cutting any of these further, or below half a percent on any channel, needs this
measurement redone first rather than the old 6e-6 quoted at it. The three
ratios named above are where that bites first.

### Re-measured for the OU families' nine gates

The relative quantum and the two new attitude bars, above, are a cut, so the
measurement was redone first rather than the paragraph above quoted at them.
Both simulators were rebuilt at `-march=x86-64` and all eight records rescored,
and the drift is taken **per binding record** — the record that actually sets
each bar — rather than as a worst-channel figure, since that is the number the
bar has to survive:

| family | gate | margin | drift on its binding record | ratio |
| --- | --- | --- | --- | --- |
| OU-II | Z %Hs JONSWAP | 0.50% | 2.2e-5 | 227x |
| OU-II | Z %Hs PM-Stokes | 0.51% | 1.0e-4 | 51x |
| OU-II | yaw deg | 0.51% | 1.7e-4 | 31x |
| OU-II | roll deg | 0.52% | 1.2e-4 | 44x |
| OU-II | pitch deg | 0.52% | 5.6e-4 | **9.3x** |
| OU-II | 3D % JONSWAP | 0.54% | 8.0e-6 | 675x |
| OU-II | 3D % PM-Stokes | 0.50% | 1.4e-4 | 37x |
| OU-II | acc Z bias % | 0.51% | 4.6e-5 | 111x |
| OU-II | bias 3D % | 0.50% | 2.8e-4 | 18x |
| OU-III | every gate | 0.51–0.57% | 2.2e-6 – 3.0e-5 | 169x – 2400x |

All eighteen pass on both builds. OU-III is comfortable everywhere: its worst
binding-record drift is 3.0e-5, on pitch, which is 169x inside that bar.

OU-II's pitch, at 9.3x, is now the thinnest margin-to-drift ratio in the whole
document — thinner than the two 15x cases named above — and it is thin because
the channel is drifty rather than because the bar is tight: 5.6e-4 relative
between builds against 2.2e-5 for the same family's vertical channel. It still
clears on both builds (0.48% of headroom on the `x86-64` one), and it is the
first bar to re-measure rather than re-cut if a rebuild ever breaches it. No OU
gate should be cut below half a percent, and this one should not be cut at all
without a fresh measurement.

Reproduce with:

```
make -C tests/<dir> clean
make -C tests/<dir> build CXXFLAGS="-O3 -std=c++20 -Wall -Wextra -Wshadow \
  -Wconversion -funroll-loops -fno-finite-math-only -I$PWD/src \
  -isystem /usr/include/eigen3 -march=x86-64"
```

and compare the `Angles RMS` and `XYZ RMS` lines against a native build. Note
that `make clean` in a test directory deletes `*.csv`, which includes the
unpacked wave records — re-fetch them with `make ensure-sim-data` afterwards, or
delete only `*.o *.d` and the binary. For the two OU families,
`tools/ou_regauge_gates.py --family ou_ii --json before.json` dumps every
channel of every record, so the two builds can be differenced directly instead
of by eye.
