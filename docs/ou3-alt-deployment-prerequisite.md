# ALT finite-word deployment prerequisite

## Saturating magnetic counters

The shipping inner magnetic-attempt counter and MagAutoTuner's accepted and
rejected sample counters saturate at `std::numeric_limits<int>::max()`. The
named deployment proof uses signed int32, M=2^31-1. Resets set each count to zero;
only a count strictly below M is incremented. No evaluated increment can exceed
M. All acquisition/refinement rejection paths use the same guarded helper.

For entry count c in [0,M] and n finite attempts on that counter,

`c_n = min(M, c+n)`.

Induction is exact: n=0 gives c. If c_n<M, the guarded successor is c_n+1<=M;
if c_n=M, the successor is M. These cases give the formula at n+1 and preserve
0<=c_n<=M for every finite n. No positive minimum call gap, uniform count cap,
or replay of billions of events is needed. A reset simply starts another such
prefix at zero. This closes integer safety, not floating-point totality.

For every representable unlock threshold u in [0,M],

`min(M,c+n) >= u  iff  c+n >= u`.

If c+n<M the sides coincide; otherwise both hold since u<=M. The implication
also holds after changing u through the public setter. Saturating at the
current configured threshold would not preserve later threshold increases.
MagAutoTuner uses the same fact with `max(1,min_samples)` and its positive-count
checks. Rejected count is diagnostic and has no sample-acceptance gate.

Only the counters stop increasing. Inner MEKF measurements, first-attempt time,
strict elapsed-time comparison, H18/A21 release and external-hold handling still
execute. Accepted magnetic samples continue updating sums, weight and window;
rejected samples still follow the same rejection return path. Floating sums and
windows need their separate arithmetic qualification.

## Source binding and composition

`finite_mag_counter_saturation.py` binds the invariant to the complete wrapper
and MagAutoTuner source hashes, including resets and configuration, and the
normalized inner update body. A source change fails the audit until re-reviewed.
The native regression checks the actual headers under signed-overflow UBSan:
it reaches gauged Live through public calls, installs M-1 only for the boundary
check, and verifies measurements and delayed release continue across saturation.
It separately crosses accepted/rejected tuner boundaries and checks statistics
continue. This is a boundary regression; the induction supplies arbitrary n.

The finite control, acquisition and refinement relations use the same saturated
successor. Proof bookkeeping retains an unbounded mathematical event ordinal;
its relation to the shipping counter is `min(M, attempted_calls)`, not equality.
The source-bound magnetic adapter consumes the shared event composer, retaining
this projection across a subsequent IMU transition.

A certified empty pre-Live magnetic prefix can enter the admitted source word
at an ungauged timeout. Waiting calls have physical packets but no inner
measurement forcing. Actual later north acquisition initializes the service
clock at that event, preserves the Live/S origin and BIAS history, and joins the
following source-owned IMU edge. Initial north uses MEKF tilt after Live;
continuous calibration/refinement use private-observer tilt.

## Actual AtomS3R caller

The committed sketch is
`sensors/full_marine_ins/atomS3R_ins_kalman_ou3/atomS3R_ins_kalman_ou3.ino`.
Its sequential `updateFilter_` contains one `fusion_.update` followed by at most
one `fusion_.updateMag`. The wrapper forwards at most once; the counter resets
to zero and has no other increasing write. Thus, on this caller, after K IMU
invocations the counter is at most K, as long as execution is defined up to
that point. For K<=INT_MAX this inductively rules out counter overflow itself.
It does not rule out unrelated arithmetic failures.

If startup ended by invocation 30,002 and the following word had 600 such
invocations, the total count would be at most 30,602, safely below INT_MAX.
This is a useful caller theorem with explicit premises. It does not prove that
startup ends then or attach this caller to the generic master. Counter safety
itself follows for every finite prefix from saturation, independent of K.

| Committed caller/build fact | Consequence for attachment |
| --- | --- |
| Nominal loop rate 200 Hz; dt from unsigned sample-timestamp differences | Nominal rate does not imply the proof's exact 5 ms grid. |
| dt falls back only when nonpositive or nonfinite | Positive jitter and long sample gaps remain possible. |
| Freshness gate uses 35 ms since last valid sample | An invalid sample resets the gate; 35 ms is not an unconditional minimum gap. |
| At most one magnetic call per IMU invocation | Gives the count bound above even when the freshness gate resets. |
| Caller sets magnetic delay 0 s and acquisition minimum norm 5 uT | These differ from the default-header profile used by the current master. |
| CI selects `esp32:esp32@3.3.7`, M5Unified 0.2.13, M5GFX 0.2.19, Eigen package 0.3.2 | The repository does identify platform/library versions. |
| CI board `esp32:esp32:m5stack_atoms3:CDCOnBoot=cdc,USBMode=hwcdc`; added flags `-funroll-loops -fno-finite-math-only` | These declarations alone do not prove compiler operation order or libm accuracy. |

The build facts come from `.github/workflows/build.yml`; Eigen's Arduino package
version is not asserted to be its upstream Eigen version. Host sanitizer tests
are not ESP32 compiler/libm qualification. Completing the actual-caller route
requires its timestamp/source, configuration and target arithmetic contracts to
be attached together. A one-call-per-step assertion alone cannot replace the
broader asynchronous master or silently import exact 5 ms timing.

## Remaining prerequisite boundary

`finite_master_guard` consumes the counter certificate as a closed prerequisite;
`falsified_prerequisites` is empty. MAG-CALL-SCHEDULE-v1 retains its existing
minimum-service semantics. Counter safety does not require replacing it with
the actual caller or adding an artificial maximum rate.

Universal startup/alignment/capture, Eigen and target/compiler/libm
correspondence, source-uniform arithmetic supplies and the complete same-history
600-step machine word remain unqualified. The master lists each dependency and
keeps storage/rho search and all ALT PASS flags false. In particular a native
host check cannot qualify ESP32 arithmetic, and source provenance alone cannot
prove all admitted histories reach Live. The binary32 wrapper clock's late-time
stall also remains a separate indefinite-continuation obligation.
