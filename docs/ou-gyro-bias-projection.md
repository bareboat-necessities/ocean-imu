# Residual gyro-bias estimate protection in OU-II and OU-III

The estimator's gyro-bias state represents residual sensor bias, so it is
constrained to a deliberately generous physically plausible region. This also
prevents discrete-time attitude propagation from entering rotation aliases
caused by absurd internal bias estimates. The common fixed Euclidean radius
is **0.5 rad/s**, about 28.65 degrees/s. It is an estimator invariant, not a new
physical motion or sensor-bias assumption.

## Evidence and radius selection

The paired audit starts from `b0beb9c88737ffeb2c6f6eeefcbeb4bd22121524` on main,
which includes the current OU-III stability work. Read-only source taps record
every pre-projection corrected bias and every literal prediction argument
`dt * norm(last_gyr_bias_corrected)`. No attitude increment or measured rate is
clipped. The release dataset is v1.2.1, SHA-256
`6d6eb92db78e97f1d2456c6b92387993bf23be14a9a6cf26943f0a22fc6fee60`.

| Pre-change run | Max estimated bias, rad/s | Max prediction argument, rad |
|---|---:|---:|
| OU-II complete deterministic simulation suite | 0.00453739418 | 0.00532237681 |
| OU-III complete deterministic simulation suite | 0.00625276338 | 0.00532258244 |
| OU-II startup regression | 0.00232140642 | 0.00145344708 |
| OU-III startup regression | 0.00202080640 | 0.00145132385 |
| OU-II stationary-device stress | 0.02499966138 | 0.00049173543 |
| OU-III stationary-device stress | 0.02478606482 | 0.00049177997 |
| OU-II paired smoke validation (12 replays) | 0.00305672180 | 0.00207509899 |
| OU-III paired smoke validation (16 replays) | 0.00416276962 | 0.00207519865 |
| OU-II full paired validation | 0.01652885530 | 0.00545292124 |
| OU-III full paired validation | 0.03682743654 | 0.00545307943 |
| OU-III full robustness study | 0.01538185207 | 0.00283334343 |

The standalone suites and paired smoke replays have zero projection activations
and zero nonfinite observations after the change. The deterministic suites
change displacement RMS by at most 7.42 micrometers and roll/pitch/yaw RMS by
at most .000101 degrees. Paired smoke changes displacement RMS by at most
.000063 m and attitude RMS by .000047 degrees. These are compilation-scale
floating-point differences; the inactive helper itself leaves the stored
bias bit-for-bit unchanged. Existing simulator gates pass.

The full validation compares 840 scored rows plus 18 calibration replays.
Every generated input has its row count checked and its SHA-256 verified
before and after execution; paired input records and reference-motion metrics
match exactly. Prediction and correction counts also match. The largest
displacement RMS change is .000044 m and the largest attitude RMS change is
.000431 degrees, with zero projection activations or nonfinite observations.
The full OU-III robustness protocol adds 310 scored cases and one calibration
replay, including 72,631,734 prediction steps in each build. Its paired inputs,
reference metrics and operation counts also match exactly. Projection remains
inactive, with maximum changes of .000002861 m in 3D displacement RMS and
.0000716 degrees in attitude RMS. The unchanged study protocol retains
individual stress-case quality outcomes; these comparisons do not relax gates.

The stationary-device tests include a 0.025 rad/s imposed sensor offset,
rest/wave/rest transitions and fixed/adaptive pseudo cadence. Startup tests
include large heading transitions, magnetic acquisition and proxy handoff.
All four investigated radii, 0.1, 0.2, 0.5 and 1 rad/s, were inactive in these
baseline runs. The selected 0.5 radius leaves 25 times the qualified physical
residual, 20 times the largest stationary-device stress estimate, and about 80 times the largest
ordinary OU-III estimate in the deterministic suite. The broader full validation
still leaves over 13 times headroom. Smaller radii offer no needed alias benefit; 1 rad/s
offers extra headroom with no observed engineering need.

The physical **total post-calibration residual** qualification remains
`B_g = 0.02 rad/s`. It covers residual offset, temperature and calibration/scale
error on the assembled sensor over its qualified operating range. It is not a
claim about every uncalibrated device. In `GyroCalibration::apply`, the input is
`S * (raw - biasT.bias(temp))`. Both AtomS3R OU sketches pass `applyGyro` output
to fusion. Their extra stillness EMA serves the displayed rate of turn, not
the MEKF. The runtime calibration's slope limit is 0.001 rad/s/degree C per
axis; the generous estimator radius allows substantial temporary residual
error beyond the physical qualification. Calibration is unchanged.

## Implemented operation order

`KalmanOUCoreMath.h` supplies one radial projection helper. Construction and
truth initialization set gyro bias to zero. Both filters guard retained means
at attitude/accelerometer initialization, attitude handoff, frame retargeting,
every prediction entry, and the common correction-injection entry. All accepted
accelerometer, magnetic, position, velocity, vertical-velocity and OU-III
integral corrections reach that entry. Projection precedes the invalid-attitude
early return. Magnetic reference refinement and proxy handoff cannot bypass it.

Finite interior estimates are untouched. Outside estimates keep their direction
and are scaled to the sphere with a tiny inward floating-point allowance
(8 scalar epsilons in the target radius). Scaling before taking the norm avoids
overflow even for finite FLT_MAX/DBL_MAX components. Exact axial boundary values
are unchanged; an outward norm check protects the stored Euclidean invariant
near the boundary. A non-finite bias is reset to zero, consistently with the
existing accelerometer-bias recovery convention. Invalid accelerometer packets
or temperatures are rejected before a correction; magnetic packets already
have a finite-input guard.

This changes **only the gyro-bias mean**. It neither zeros nor projects nor
rescales covariance. Joseph updates and existing covariance/reset operations
remain in their original order. No noise, OU dynamics, cadence, quality gate,
heading acquisition, calibration or magnetic startup parameter changes.
TFG, NLO and PII do not call the new helper and are unchanged.

## Qualified prediction-angle bound and its limit

For the existing OU-III qualified source family, `constants.json` admits
`0.004 <= h <= 0.006 s`, physical `Omega_max = 0.6108652381980153 rad/s`,
total calibrated physical residual `B_g = 0.02 rad/s`, and fast gyro measurement
residual `N_g = 0.02 rad/s`. Deheeling is a rotation, so in real arithmetic

```
|omega_measured - b_hat_g| <= Omega_max + B_g + N_g + R_g
                         = 1.1508652381980153 rad/s
h |omega_measured - b_hat_g| <= 0.0069051914291880918 rad.
```

Both residual terms are charged; calibration residual is included in `B_g`,
not silently dropped or counted twice. The normalized small-angle quaternion
polynomial adds at most `4(phi^6/46080 + phi^7/645120)` to the real-source
rotation angle. Thus the conservative cap is **0.007 rad** (0.402 degrees),
with margins exceeding **3.13459 rad from pi** and **6.27618 rad from 2*pi**.
The exact `h=.005, b_hat_g=-400*pi*e_z` bias witness and its near-alias family
are outside the implemented estimator-state ball.

The polynomial allowance follows directly from the shipping fourth-order
cosine and fifth-order sine terms. Their Taylor remainders have combined
Euclidean norm at most `epsilon = phi^6/46080 + phi^7/645120`. Normalization
changes direction by at most `asin(epsilon/(1-epsilon))`; the corresponding
rotation-angle error is at most twice that value, hence below `4 epsilon`
for this domain (`epsilon < 1/4`). This is a real-source bound; finite-precision
word totality remains a separate proof obligation.

**The shipping device API does not enforce a maximum positive timestep.**
`SampleDtTracker` is constructed without `max_dt_s` in these sketches; the
fusion wrapper accepts any finite positive `dt`. Therefore no finite global
prediction-angle cap exists over every API-admitted device stall. The above
is a source-uniform bound over the already qualified 4--6 ms family, not a new
runtime timing guarantee. Nor does the marine angular-rate qualification cover
arbitrary measured gyro spikes. A global API guarantee needs a separate timing
and input qualification policy; this change does not alter those policies.
Some deterministic simulation motions also exceed the proof's marine rate
ceiling; their angles are reported as observations, not theorem evidence.

## OU-III proof consequence

This operation enters the existing historical AG readout prerequisite to the
finite-error service-word inequality; it does not create a second architecture.
For the real Rodrigues branch, with `phi=h|omega|`, the axial singular value
of the bias transport is `h`, and the two transverse values are

```
2 |sin(phi/2)| / |omega| = h sinc(phi/2)
 >= h (1 - 0.007^2/24)
 >= 0.003999991833333333... s > 0.
```

The separate shipping `|omega|<1e-7` polynomial branch has symmetric transverse
part `h(1-phi^2/6) I`, which gives a larger lower bound. The zero-rate value
is exactly `h I`. The rational certificate and its branch checks are in
`gyro_bias_projection.py` and `gyro-bias-projection.json`.

This closes the **one-prediction real-source gyro-transport separation**.
It does **not** close `inf_W Delta_gyr(W)>0` for the full signed temporal word:
chronological observation forcing, actual reset transport and projection
defects still need their joint uniform enclosure. Force/field collinearity,
historical action, full covariance upper bound, strict contraction, finite
capture/release, retention and whole-word float32 totality remain open.

For ideal Euclidean projection let `d_g=b_corr-Proj(b_corr)`. Since physical
`|b_g|<=B_g<R_g`, `e_g^+ . d_g <= -(R_g-B_g)|d_g|`, with `R_g-B_g=.48`.
This component sector is not contraction in the full covariance metric:
carry `E_g d_g` as a separate signed mean defect and retain all cross terms.
The implementation's inward rounding allowance is an additional arithmetic
defect, not a covariance Jacobian/reset. The construction mean-action observer
explicitly verifies inactivity before using its unprojected affine recursion.

## Reproduction

`tools/ou_gyro_bias_audit.py --family ii --output-dir /tmp/ou2-audit` (and `iii`)
builds source taps and records counts, maxima, candidate exceedances and source
hashes. Add `--source-ref b0beb9c88737ffeb2c6f6eeefcbeb4bd22121524` for the
baseline. To include the existing paired validation protocol, first build both
families, then pass `--validation-mode smoke --validation-peer-dir /tmp/ou2-audit`
on the OU-III command. The official aggregate requires both families.
For complete studies use `--validation-mode full`; add
`--validation-study robustness` for the existing robustness protocol.
`--no-build --replay-only` reuses the frozen binaries. Scalar replay checkpoints
are keyed by binary bytes, input bytes/name, arguments and all W3D/SF settings;
time-series diagnostics always execute again. A changed or incomplete generated
input fails the audit. `TMPDIR` can select storage for the large temporary inputs.
The committed paired report records numerical metric differences
and projection engagement. New regression targets run from both existing
`run_tests.sh` scripts; the ordinary suites and all publication gates remain
required. Finite replays do not establish an all-history stability theorem.

Reproduce the comparisons from retained audit directories with
`python3 tools/ou_gyro_bias_audit.py --compare BASELINE_DIR CANDIDATE_DIR --output-dir COMPARISON_DIR`.
This writes `COMPARISON_DIR/paired-audit.json`, including metric deltas, observation
maxima, projection activity, paired input hashes and measured frozen tuning
points. Run it for each standalone family and each full study; use
`--compare-part standalone` to compare only the native suite in a directory
that also contains replays. The command verifies retained stdout and rejects
incomplete/unverified replays, differing physical inputs, reference metrics or
protocols. Projection activity is reported even when nonzero, so an unexpected
activation cannot disappear from the diagnostic.
