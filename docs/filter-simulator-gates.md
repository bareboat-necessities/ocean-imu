# Filter simulator regression gates

External-data regression now uses oceanography-waves-lib v1.3.1,
`sim-data-files-vessel-rao-28ft-5kn.zip`: a 28 ft estimated hull response at
5 kn with heading oscillating ±3° about north every 120 s. The dataset change
has its own measured baseline. The accepted encounter fixes preserve wave
variance through stillness detection, normalize the shared CARMA vertical
band, and use the wave-matched TFG drift observation. Gate tightening locks
in their measured improvements; it does not itself improve estimator accuracy.

See [current measurements](moving-vessel-results.md) and the
[track/scoring convention](moving-vessel-replays.md). All physical sensor
inputs and turning acceleration remain intact. Only reference displacement
and velocity are expressed about the declared cruise track. In particular,
3D error percentages are normalized by wave displacement, not total travel.

`tests/common/CurrentReplayLimits.h` contains separately qualified production
and interval-rate limits. Each channel uses the worst of eight full records,
plus 0.5%, rounded upward to 0.00001, capped by the previous gate. No gate
may widen during this tightening. For the updated TFG/shared profiles, use
the maximum of retained native and available CI controls; shared TFG bars
also cover the tangent consumer. Both production and historical interval-rate
profiles are checked. `gate_tightening.json` retains those controls and every
before/after limit. `limits.json` records the exact compiled values.

The old drifting-vessel limits remain accuracy-recovery targets, not an
achievement implied by passing a higher moving-vessel gate. A reduction in
mean RMS does not establish recovery of the worst record or of attitude/bias
accuracy. The following historical comparison makes the unresolved gap
explicit; the old limits used v1.2.3, whereas current replays use v1.3.1's
corrected clocks. Same-executable speed comparisons are documented separately
in [moving-vessel-comparison.md](moving-vessel-comparison.md).


| Filter | Old drifting gate JON/PM (% Hs) | Previous 5 kn gate | Tightened 5 kn gate |
| --- | ---: | ---: | ---: |
| OU-II | 6.662/6.462 | 8.525/8.068 | 8.525/8.068 |
| OU-III | 4.271/4.234 | 8.277/7.878 | 8.277/7.878 |
| TFG | 4.405/4.356 | 7.858/7.513 | 5.745/5.486 |
| CARMA | 3.072/3.231 | 3.590/3.745 | 3.590/3.745 |
| Shared QMEKF | 2.467/2.415 | 3.568/3.720 | 3.473/3.549 |
| TFG-CARMA | 3.072/3.231 | 3.590/3.745 | 3.590/3.745 |
| Shared TFG / tangent | 2.467/2.415 | 3.568/3.721 | 3.474/3.551 |
| TFG legacy | 4.344/4.317 | 8.288/8.056 | 4.713/4.857 |

OU-II, OU-III and independent CARMA limits are unchanged here: their current
worst scores still need those ceilings. PII/NLO and heave-only limits were
already tightened with their accepted changes and remain separate; the
process-prior and lever-arm simulators have their own diagnostic contracts.
The present 61 reductions cover heave, 3D displacement, attitude and bias
channels where current evidence permits tightening. They do not establish
complete free-drifting accuracy recovery.

W3d scores use the trailing 900 s of each 1200 s, 200 Hz record. Heave is
RMS/incident Hs; 3D displacement is RMS/max reference wave-displacement norm.
Accelerometer Z bias uses max absolute true Z bias; each vector-bias score
uses max true vector-bias norm. Absolute displacement and bias RMS accompany
these different percentage normalizations. AHRS uses the trailing 60 s.
NLO gates raw heave, including its mean error. Adaptive PII has separately
measured vertical/yaw bars with the same 0.5% margin, rounded upward to
0.01 in its quoted units; it rejects nonfinite heave and attitude scores. Heave-only gates remain
separate for Mahony, proxy and ideal-attitude frontends and each method.

The inventory covers all native W3d filters, the legacy cubic TFG diagnostic,
AHRS, heave-only methods, NLO, the independently checked CARMA process prior,
and the lever-arm study instrument. Detrender inputs also use wave-relative
horizontal displacement from the moving records; cruise travel is not added
to their wave-reference target.

Nonfinite scored channels, missing/duplicate records, incomplete histories
and wrong input hashes remain failures. Full-rate samples determine RMS;
plot decimation never changes scoring. The conditional eight-record gates
are deterministic sentinels, not population confidence intervals, general
accuracy targets or guarantees for engine noise and alternate priors.

After building, reproduce the passing evidence with:

```sh
make all
python3 tools/filter_sim_gates.py --run --jobs 3 \
  --log-dir /tmp/filter-sim-gates \
  --output reports/results/filter_sim_gates/measurements.json
python3 tools/moving_vessel_publication.py
```

`--measure` permits collecting a changed experiment before gate qualification.
It retains failure exit statuses and still requires complete finite scores;
it neither changes limits nor produces a passing publication report. The
ordinary command requires every gate to pass. Both record input and binary
hashes and the source revision. Dense outputs and runner logs are omitted.
