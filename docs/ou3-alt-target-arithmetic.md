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
compiler tracks are attached to this object by the expression-specific mapping
below. Its use by the final firmware remains an independent premise.

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

## Scalar profile and math namespace attachment

`target_scalar_profile.py` identifies the actual SDK ROM division entry
`0x40002274`, follows its thunk to `0x40056124`, and verifies the complete
instruction sequence against the pinned libgcc implementation. The
[Cadence ISA manual](https://www.cadence.com/content/dam/cadence-www/global/en_US/documents/tools/silicon-solutions/compute-ip/isa-summary.pdf)
specifies the IEEE divide and square-root sequences and final correction.
The checked profile explicitly requires `FCR.RM=0` throughout execution.
Initialization and preservation of that mode remain part of the complete
execution qualification.

`target_math_link.py` supplies a separate static namespace check. It compiles
actual C++ `std::exp`, `std::log`, and `std::sqrt`, links the pinned libraries
with the SDK ROM script, and verifies selection of the audited exp/log objects,
the IEEE square-root body, and the ROM divide entry. All text and literal
sections are placed together at the SDK instruction-segment origin
`0x42000020`, keeping Xtensa windowed longcalls in the same 1 GB segment
as the ROM entry.

The namespace probe is not a complete firmware. Bare libc reports the missing
`__getreent` runtime service; that diagnostic is retained, and runtime/full
firmware readiness remains false. It is not executed and does not supply
replacement math functions. Within the certified finite argument domains,
exp/log/sqrt wrappers return their kernel values without entering errno paths.

`target_wpe_libm.profile_correspondence()` composes the named scalar profile,
verified namespace, and all-input exp/log error certificates. It also attaches
the correctly-rounded IEEE sqrt result under the same profile. The narrow
library obligations are now closed under that explicit execution profile:
WPE exp/log approximation (qualification 8) and Qaxis exp error fitting the
original enclosure (qualification 9). Complete firmware compiler/link
correspondence is still qualification 5; this split does not promote it.

```sh
python3 -m tools.stability.ou3_alt_contraction.target_math_link \
  --toolchain /path/to/xtensa-esp-elf --sdk /path/to/esp32s3-libs \
  --output /tmp/ou3-target-math-link.json
```

## WPE compiler selection

`target_wpe_compiler.py` compiles the pinned WPE header with the actual MCU
flags and maps every contraction site to the existing binary32 relation.
The mapping includes the complete method's eight MADD and three MSUB sites;
source hashes, horizon-limit hashes, object identity, source-line register
operands and all exp/log/sqrt/divide call sites are recorded.

| WPE expression | Actual target evaluation |
| --- | --- |
| Two high-pass stages | Separate ADD, SUB, MUL |
| Velocity/elevation integration | Round `gain*input`, then FMA `decay*old` into it |
| Weight | FMA `(1-alpha)*old + alpha` |
| First moments | Round `alpha*value`, then FMA `(1-alpha)*old` into it |
| Second moments | Round `alpha*value`, round that times `value`, then FMA `(1-alpha)*old` into it |
| Variances and omega squared | Fused subtraction of the same already-computed operand square |
| Log-period EMA | Round `log_raw-log_previous`, then FMA with alpha and the previous log |

The deterministic target-step constructor derives the selected successors,
all early-return decisions, raw period, log update and usable latch from one
persistent state and the same-argument library outputs. It checks each chosen
moment transition through `finite_wpe_moment_binary32.step`, constructs the
existing `ModeWitnesses`, and `attach_to_product` checks that the successor is
exactly the product's persistent `.fma` projection. Reset initializes that
projection, so induction gives the same target history at every prefix.
The separate compiler track remains independently carried. The diagnostic
reset-horizon marker has no effect on any transition and is explicitly
projected out until the first horizon calculation.

This closes the named WPE compiler-selection obligation (qualification 6)
under the canonical 5 ms/default WPE configuration, the physical MEMS input
contract and the qualified scalar profile. It does not claim that the full
firmware uses this compiled object; that remains qualification 5. It also
does not replace the frontend's input-supply proof or establish startup
capture.

```sh
python3 -m tools.stability.ou3_alt_contraction.target_wpe_compiler \
  --toolchain /path/to/xtensa-esp-elf --sdk /path/to/esp32s3-libs \
  --output /tmp/ou3-target-wpe-compiler.json
```

## Tuner `powf` finite-range supply

The candidate tuner has two compiled exponents, `6/7` and `1/14`. The pinned
newlib `ef_pow.c` object is audited separately by
`target_powf_range.py`. For the complete source-owned base interval
`[2^-44, 2^17]`, its log2/reduction/exp2 graph has only finite normal
intermediates and returns a strictly positive result in the conservative
range `[2^-45, 2^18]`; both exponent constants are the actual binary32
values. This is a range/totality result, not an approximation-accuracy or
whole-firmware call-graph result. `finite_candidate_uniform_bounds.py`
consumes it to keep tau, sigma and `R_S` candidates finite at every prefix.
The target compiler/FCR and startup-root qualifications remain separate.
