# ALT commissioned Live input contract

The ALT theorem quantifies over post-Live raw sensor histories satisfying the
following bounds at every actual `SeaStateFusion_OU_III::update` call:

| Filter argument, body frame | Per-axis absolute bound |
| --- | ---: |
| Gyroscope, radians per second | 35 |
| Accelerometer specific force, metres per second squared | 160 |

The arguments must be finite, fresh, correctly mapped and converted, and use
commissioned calibration. Unsaturated acquisition is required when connecting
the sample to the physical source/error model. The bounds apply **after** any
application calibration and **before** the filter's vibration guard. They are
theorem admission conditions; shipping does not acquire a new clipping step.

Formally, the existing COMPLETE-BRMM and BIAS0/1/2 histories, startup contract,
and bounded ISS forcing condition are intersected with

`forall k after Live, max_i |gyro_body[k,i]| <= 35`

and

`forall k after Live, max_i |acc_body[k,i]| <= 160`.

The source and disturbance identities persist across the intersection. The
smaller commissioned startup residual caps and temporal direction conditions
remain in force during startup. This raw-input envelope does not enlarge them
or assert startup capture. The scalar finite-word proof still uses the same
canonical 5-ms source period and shared startup-plus-600-sample budget.

## Device and application basis

The [Bosch BMI270 datasheet, revision 1.6](https://www.bosch-sensortec.com/media/boschsensortec/downloads/datasheets/bst-bmi270-ds000.pdf)
specifies selectable maximum ranges of ±2000 degrees/s and ±16g. The
[MPU6886 datasheet, revision 1.1, hosted by M5Stack](https://m5stack.oss-cn-shenzhen.aliyuncs.com/resource/docs/datasheet/core/MPU-6886-000193%2Bv1.1_GHIC_en.pdf)
gives the same maxima. The selected bounds enclose those ranges: 2000 degrees/s
is approximately 34.9066 rad/s, and 16g is 156.9064 m/s².

The actual app mapping is `(gy,gx,-gz)*DEG2RAD` and
`(ay,ax,-az)*9.80665` in `src/AtomS3R/AtomS3R_ImuCal.h`.
The OU3 sketch passes `runtime_.applyGyro` and `runtime_.applyAccel` results
into the wrapper. Its source contains no range-selection setter; CI pins
M5Unified 0.2.13. The selected maximum-range envelope does not claim that a
particular connected board has a particular range register configured.
Calibration at the input boundary must satisfy the selected caps; arbitrary
calibration matrices are not inferred to preserve them.

The exact literal binary32 conversion proof leaves over 0.09 rad/s gyro margin
and over 3 m/s² accelerometer margin. Datasheet RMS noise or sensitivity
statistics are not being turned into deterministic error guarantees.

## Executable admission and proof use

`tools/stability/ou3_alt_live_input_domain.json` contains the fixed numerical
profile. `finite_live_input_contract.check_packet` checks the source-owned
packet and its exact binary32 API conversion. Both theorem-facing admitted IMU
entry points call it before prediction coefficients, state updates or forcing
results are constructed. The packet object and existing physical/bias/source
ordinals are preserved. The bounded-forcing history carries the same fixed
input profile. Freshness and physical restrictions continue through the
existing source-ordinal checks; magnitude checks do not prove those conditions.
The caps constrain the actual rounded API arguments. The exact comparison
shadow retains their full RNE cells, including the half-ulp interval above a
cap that still rounds to that cap; no extra real-shadow restriction is added.

The contract is used in a scalar private-Mahony induction on **every finite
initialized prefix**. It starts from the source-bound seed and retains the
integral across Live. The per-call integral increment is below
0.000130000010. At the binary32 value 4096, the upward half-ulp is
0.000244140625, so `RN(4096 + increment) = 4096`. Monotonicity of rounding and
the symmetric negative endpoint make `[-4096,4096]` forward invariant. This is
an invariant of the actual float recurrence, independent of elapsed duration.

The previous all-finite-word inverse-square-root lemma returns every defined
normalized successor to `||q||² < 1.112`. The raw acceleration cap keeps its
three-square sum below 76,801; the integral invariant keeps feedback-corrected
gyro components below 4131.261 rad/s and the Euler quaternion norm sum below
4597.799. The observer's elapsed float also has an exact rounding barrier at
131,072 seconds. Every state-changing scalar operation remains finite under
the named RNE, gradual-underflow, no-FMA profile. Clock saturation is not a
guarantee of indefinite time progress or correctness of all wrapper services.

The sharper 30,602-sample consequence remains available when its horizon
premise holds: integral below 3.986, corrected gyro below 39.246 rad/s, and
Euler norm sum below 7.530. Those sharper values are not used to prove the
arbitrary-length scalar totality result.

For the same quaternion and acceleration packet, the exact down-row norm is
`||d(q)|| = ||q||²`. Retaining that identity, its polynomial roundoff and the
same input gives a vertical output bound below 321.169 m/s². The WPE proof may
therefore use the common ceiling 512 without assuming tilt accuracy. This
argument is conditional on the inherited seed and dormant guard premises;
it does not replace startup by an installed Live state.

The former `2^80` gyro pulse is rejected before execution and remains an
outside-contract regression. The raw-input overflow issue is resolved for the
specified private-observer arithmetic. Full MEKF/covariance/LDLT, magnetic and
tilt-event deployment bounds, target correspondence and the complete carried
machine word remain separate obligations. In particular, an input gyro cap
does not itself bound the MEKF's estimated gyro bias.

The shipped AtomS3R OU3 application calls the three-argument `update`, whose
temperature default is exactly 35°C. Its internal MEKF thermal term is
therefore exactly zero; the source-checked certificate records that call-site
consequence. The general four-argument API and existing explicit thermal
history graph remain unrestricted by this observation. Their bound would
require the temperature/model term separately: the accelerometer bias clamp
does not constrain the added deterministic `k_a*(T-35)` term.
