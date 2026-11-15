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

Every default replay family now applies **per-sea limits**. A quiet or moderate
record cannot consume the error budget of the worst sea. The table in
`tests/common/RecordReplayLimits.inc` covers 256 family/profile/sea combinations
and 1,120 active scalar limits; 991 are lower than their previous global caps.
Coverage includes all OU and CARMA configurations with qualified default bars,
both TFG regularizers, shared QMEKF/TFG and the tangent consumer, AHRS, adaptive
PII, NLO, and all five heave-only methods on their supported frontends.

Each limit is the maximum retained native/CI score **for that record and
profile**, plus 0.5%, rounded upward to 0.00001, capped by its existing global
limit. No limit widens. `per_record_limits.json` retains controls, measured
maxima, previous caps, final limits and known pre-existing failures.
`RecordW3dLimits.h` adapts the common table to the existing W3d summary.
Unknown families, profiles, records and unmeasured channels fail closed within
this fixed-draw protocol. Explicit sensor-seed sweeps retain their original
global diagnostic bars; they are not used to qualify these deterministic caps.
Heave-only `HO_*` overrides that alter the active method or frontend likewise
retain the old diagnostic bars, including Kuchler body-Z input. With default
seeds, `NLO_NOGATE` suppresses numerical-limit failures only: unknown records
remain protocol errors and abort before a quality summary.

`CurrentReplayLimits.h` retains the global production and interval-rate ceilings
as additional upper bounds. Its 61 post-encounter reductions and their controls
are recorded in `gate_tightening.json`; `limits.json` records those exact global
values. Shared TFG and tangent use the same per-record envelope, containing
controls from both consumers. The generated table can be reproduced with
`python3 tools/generate_record_replay_limits.py`; `--check` verifies its bytes.

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

These are the maximum global ceilings, not the limit applied to every sea.
The following counts describe the additional per-record tightening:

| Filter | Lower per-record limits / active limits |
| --- | ---: |
| OU-II | 116 / 128 |
| OU-III | 108 / 128 |
| TFG | 108 / 128 |
| CARMA | 113 / 128 |
| Shared QMEKF | 106 / 128 |
| TFG-CARMA | 120 / 128 |
| Shared TFG / tangent | 106 / 128 |
| TFG legacy | 54 / 64 |
| AHRS | 24 / 24 |
| Adaptive PII | 16 / 16 |
| NLO | 8 / 8 |
| Heave PII | 16 / 16 |
| Godhavn | 24 / 24 |
| Extended Godhavn | 24 / 24 |
| Kuchler | 24 / 24 |
| Richter | 24 / 24 |

OU-III and TFG host fixtures use the fixed Linux x86-64-v3 (AVX2/FMA) CPU
target in `tests/common/HostReplayFlags.mk`. These host tests require a CPU
supporting x86-64-v3; other platforms retain their existing build settings.
Firmware flags and estimator source are unchanged. A different host arithmetic
profile can be built diagnostically with `make -B HOST_REPLAY_ARCH_FLAGS=...`.

Run 43734792425 exposed why the CPU target is part of the experiment:
`-march=native` selected AVX-512 on one runner and changed Eigen's arithmetic.
On the same local machine a clean native OU-III build gives Hs=4 m JONSWAP
yaw 0.886416°, while x86-64-v3 gives 0.871181° and restores the saved AVX2
result exactly. The interval result is 0.863419°. Both pass the unchanged
0.87554° / 0.86774° limits. The TFG comparison's Hs=8.5 m JONSWAP interval
pitch is 0.093262° versus the unchanged 0.09461° bar (CI native: 0.094639°).
All 40 OU-III, TFG and legacy-TFG sea/profile cases pass; source, binary and
input identities are in [the CPU-profile receipt](../reports/results/filter_sim_gates/host_cpu_profile.json).
The included build settings participate in OU and TFG evidence provenance,
so previous publication bundles require real regeneration, not restamping.

This repairs fixture reproducibility. It does not recover the old drifting
heave accuracy; that estimator limitation remains documented separately.

The [verification record](../reports/results/filter_sim_gates/per_record_verification.json)
retains all 32 passing native OU-II/OU-III scores, binary/input identities, and
completed hosted checks for the remaining families. The regenerated evidence
publication passed 1,210 tests with one skipped. The final heave-only guard passes all 112 default cases plus its body-Z
ablation control. Twelve focused protocol tests pass. Full local `make all`
attempts were interrupted by execution-service failures; none is reported as complete.
Hosted OU-III on source 8b5aebe passed 14/16 cases and reproduced the two
Hs=4 m JONSWAP yaw failures above. Its pre-regeneration checkout also failed
a publication-freshness check; the source-bound bundles are now refreshed.

W3d scores use the trailing 900 s of each 1200 s, 200 Hz record. Heave is
RMS/incident Hs; 3D displacement is RMS/max reference wave-displacement norm.
Accelerometer Z bias uses max absolute true Z bias; each vector-bias score
uses max true vector-bias norm. Absolute displacement and bias RMS accompany
these different percentage normalizations. AHRS uses the trailing 60 s.
NLO gates raw heave, including its mean error. Adaptive PII has separately
measured vertical/yaw bars with the same 0.5% margin; it rejects nonfinite heave and attitude scores. Heave-only gates remain
separate for Mahony, proxy and ideal-attitude frontends and each method.

The inventory covers all native W3d filters, the legacy cubic TFG diagnostic,
AHRS, heave-only methods, NLO, the independently checked CARMA process prior,
and the lever-arm study instrument. The latter two are not eight-sea filter
replays: process-prior spectral/derivative/variance checks and lever-arm numerical
contracts remain in force. Ungated CARMA research priors and shared-mode
ablations retain their diagnostic status; A41/T41/CSM caps are not fitted to
those different models. Detrender inputs also use wave-relative
horizontal displacement from the moving records; cruise travel is not added
to their wave-reference target.

Nonfinite scored channels, missing/duplicate records, incomplete histories
and wrong input hashes remain failures. Full-rate samples determine RMS;
plot decimation never changes scoring. The conditional eight-record gates
are deterministic sentinels, not population confidence intervals, general
accuracy targets or guarantees for engine noise and alternate priors.

After building, reproduce the evidence and gate outcomes with:

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
