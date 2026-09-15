# ALT target arithmetic qualification

The deployment proof uses the MCU build pinned in `.github/workflows/build.yml`:
Arduino ESP32 3.3.7, AtomS3/ESP32-S3, and Arduino Eigen 0.3.2. The official
[Espressif package index](https://raw.githubusercontent.com/espressif/arduino-esp32/gh-pages/package_esp32_index.json)
resolves this core to esp-x32 2511, GCC 14.2.0 build `esp-14.2.0_20251107`,
and esp32s3-libs 3.3.7. The verified Arduino Eigen archive contains Eigen 3.4.0.
The byte hashes and observed compiler options are recorded in
`tools/stability/ou3_alt_target_arithmetic_audit.json`. Package version labels
alone do not establish operation correspondence.

## What the target audit establishes

`target_arithmetic_audit.py` compiles the unchanged numerical Mahony, WPE and
Qaxis headers with the actual cross compiler and SDK C++ flags, including the
shipping `-Os -funroll-loops -fno-finite-math-only` settings. This is a standalone
numerical-header probe using `EIGEN_NON_ARDUINO` and `EIGEN_MPL2_ONLY`; it does
not compile the Arduino wrapper, execute the MCU or establish final firmware
link resolution. Its object contains ordinary add/multiply instructions and
`madd.s`/`msub.s`; GCC reports `-ffp-contract=fast`, with reassociation and
finite-math-only disabled. The no-contraction scalar observer lemma therefore
cannot simply be declared the deployed profile.

The tool records the multilib archive selected by the target driver and its
undefined math symbols. The library identifies itself as newlib 4.3.0. Neither
that identity nor a successful compile proves correct rounding of exp/log/sqrt.
The actual firmware link map, instruction semantics (including subnormal and
exception behavior), expression-specific contraction and Eigen reductions,
and attachment of all-input library bounds are still required. The exp/log
bounds below supply the library approximation step. WPE's two conditional
compiler tracks do not by themselves prove that either track matches this object.

Reproduce against the verified extracted packages:

```sh
python3 tools/stability/ou3_alt_contraction/target_arithmetic_audit.py \
  --toolchain /path/to/xtensa-esp-elf \
  --sdk /path/to/esp32s3-libs \
  --eigen /path/to/Eigen-0.3.2/ArduinoEigen \
  --output /tmp/ou3-target-audit.json
```

## Quantitative Qaxis library obligation

For the general Qaxis branch let `t=RN32(.01)` and `t<=x<=1`.
The existing same-argument real envelope is
`1-x <= exp(-x) <= 1-x+x²/2`. Alternating Taylor bounds give

`exp(-x)-(1-x) >= x²/2-x³/6 >= t²/3`,

`(1-x+x²/2)-exp(-x) >= x³/6-x⁴/24 >= t³/8`.

Both lower bounds exceed `2^-24`. Consequently, a separately proved absolute
error at most `2^-24` for each target `expf(-x)` call is sufficient to retain
the existing Qaxis enclosure uniformly, including the compiled threshold.
`general_branch_exp_error_budget()` checks the exact rational margins.
It neither asserts this target bound nor identifies the repeated calls with
the separate transition exp. Canonical arguments are within this domain;
startup WPE, BA decay and other transcendental calls have their own domains.

WPE retains its original correctly-rounded relation and also admits a named
tolerant profile that charges exp/log/sqrt errors throughout the actual
recurrence. The pinned-library bounds below supply the exp/log part of that
profile. Scalar-operation semantics and final-link attachment remain explicit
common target obligations; compilation alone cannot discharge them.

## Pinned expf and logf all-input certificates

The supplying algorithm is newlib commit
[`9a0d39153510ec5cbb51eb8c70cecbfeffdbb6ba`](https://github.com/espressif/newlib-esp32/tree/9a0d39153510ec5cbb51eb8c70cecbfeffdbb6ba),
tag `esp-4.3.0_20251107`. Its `ef_exp.c` and `ef_log.c` coefficient
words and branch graphs match the inspected objects in the pinned ESP32-S3
`libm.a`. The audit records the archive and member hashes; the source package
identity alone is not used as evidence of numerical accuracy.

| Target call and reached domain | Proved error, including operation rounding | Required WPE/Qaxis budget |
| --- | --- | --- |
| Qaxis `expf(-x)`, `RN32(.01) <= x <= .25` | absolute `< 4.736e-8` | absolute `2^-24` |
| WPE `expf(a)`, `-60 <= a <= 60` | relative `< 3.339e-7` | relative `2^-20` |
| WPE `logf(a)`, `2^-62 <= a <= 2^16` | absolute `< 5.249e-6` | absolute `2^-14` |

`target_qaxis_exp.py` compares the actual five-coefficient rational
approximant with a degree-12 Taylor polynomial using exact rational
coefficient arithmetic. It then charges every multiply, add, subtraction
and division. The resulting error fits both original Qaxis real-envelope
boundaries. Both repeated covariance calls satisfy that result independently;
they remain separate from the OU transition exp.

`target_wpe_libm.py` also carries exp's integer range reduction, the split
`ln2HI/ln2LO` constants, exact binary32 `k*ln2HI` products, both reconstruction
branches, and exact normal exponent scaling. For log, exact mantissa/exponent
normalization reduces the argument; both main source branches have the same
real expression `2*s+s*R`. The coefficient residual is bounded against the
convergent `2*atanh(s)` series. The small-input log and exp branches are
included. No sampled histories or host-libm accuracy claims enter these bounds.

These are conditional arithmetic theorems: RNE add/subtract/multiply,
either fused or separate RNE multiply/add, a division error of at most one
ulp, and selection of the identified library objects. The operation bound
already includes final binary32 rounding and does not require correctly
rounded exp or log. WPE's tolerant runtime profile can therefore consume it
once the common target scalar and firmware-link qualifications are attached.
Sqrt retains its separate target qualification. The library certificates do
not themselves qualify startup capture, Eigen reductions, or the full ALT
master.

Reproduce the exact arithmetic and pinned-object audits:

```sh
python3 -m tools.stability.ou3_alt_contraction.target_qaxis_exp \
  --libm /path/to/esp32s3/libm.a --ar /path/to/xtensa-esp32s3-elf-ar \
  --objdump /path/to/xtensa-esp32s3-elf-objdump \
  --output /tmp/ou3-target-qaxis-exp.json
python3 -m tools.stability.ou3_alt_contraction.target_wpe_libm \
  --libm /path/to/esp32s3/libm.a --ar /path/to/xtensa-esp32s3-elf-ar \
  --objdump /path/to/xtensa-esp32s3-elf-objdump \
  --log-source /path/to/pinned/newlib/libm/math/ef_log.c \
  --output /tmp/ou3-target-wpe-libm.json
```
