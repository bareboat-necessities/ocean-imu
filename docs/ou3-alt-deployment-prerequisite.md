# ALT finite-word deployment prerequisite

## Result

The current MAG-CALL-SCHEDULE-v1 language cannot satisfy the finite master's
requirement of defined shipping execution for every admitted event word. This
is a falsified prerequisite, not an outstanding interval estimate. Consequently
the current master cannot authorize rho estimation. The physical model, exact
conditional finite identities, attitude atlas and independent P2/P3/P4/P5 track
are unaffected. No dynamical instability claim follows.

## Exact finite schedule and counter induction

Consider an eligible gauged-Live endpoint: magnetometry is enabled, both delay
gates have passed, north is set, the inner MEKF exists and the signed counter is
c in [0, 2^31-1]. This is a conditional statement about any such entry, not a
new basin restriction. A quiet public-API startup is a native reachability
regression for this nonempty branch; it does not qualify all startup histories.

Let M=2^31-1 and N=M+1-c. At the existing physical endpoint, before IMU transition
1, make N magnetic calls. Then make calls every 40 ms forever. This satisfies
the literal schedule checker: first gap zero, burst gaps zero, tail gaps 40 ms.
The endpoint composer explicitly admits magnetic calls at the fresh Live
origin. No between-grid endpoint is introduced. All N calls use the same
physical sample and disturbance value; no new independent source is inserted.
A quiet zero-wave, zero-bias, zero-disturbance history can continue unchanged.

The schedule is locally finite. For t>=0 it has exactly
N+floor(t/0.04) calls through t relative to this origin; for t<0 it has zero.
Every compact interval therefore has finitely many calls, and the tail reaches
unbounded physical time. Local finiteness does not give a *uniform* count bound.
The burst lies in every positive finite horizon, including the 600-IMU word.

In the audited outer `updateMag`, continuous/refinement operations precede
`impl_.updateMag`, but no duplicate-time suppression or counter reset blocks
that call once gauged Live is reached. In the inner method the only early
returns are the enable/MEKF and delay gates. A returning magnetic measurement
is followed by the unchecked signed increment. Innovation rejection does not
prevent it. The later H18/A21 unlock predicates do not guard the increment.

Assume all N calls have defined execution. Induction gives counter c+j after
each of the first j returning calls. At j=N-1 it equals M. The last call then
requires M+1, which has no signed-int32 successor. This contradicts defined
execution. If a continuous-calibration, measurement or other operation fails
earlier, total execution already fails; the proof does not assume an arbitrary
number of safe covariance updates. A Python exception rejecting this edge
correctly models a partial relation but does not prove avoidance of the edge.

For c=0 the compressed witness has 2,147,483,648 calls. This is an event-language
counterexample, not a claim that the AtomS3R can physically execute that many
callbacks in three seconds. The finite-master guards may not conflate those
two deployment languages.

`finite_mag_counter_obstruction.py` retains this family without allocating or
executing billions of events. Comment/whitespace-independent fingerprints bind
the manual source-order argument to the inspected inner and outer methods.
The native sanitizer regression first reaches gauged Live through public API
calls on a constant physical source, then installs ONLY M-1 in the counter to
test the last defined increment and the following undefined increment. This
boundary installation is explicitly not a reset-to-overflow trajectory proof;
the induction above supplies the arbitrary-length argument.

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
startup ends then, attach this caller to the generic master, or cover indefinite
execution.

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

## Research decision

Freeze attempts to prove totality of unchanged MAG-CALL-SCHEDULE-v1. Interval
refinement, more trajectories and a different storage cannot change M+1>M.
The constructive alternative is a source-derived actual-caller contract with
measured-dt semantics, followed by the remaining startup/arithmetic composition.
A theorem conditional on defined execution is a different, weaker theorem and
does not satisfy the current guard. Changing the deployed counter is outside
this proof PR's immutable scope. No arbitrary call-rate cap is adopted.

`finite_master_guard` reports counter safety under `falsified_prerequisites`,
keeps the other unfinished qualifications under `open_obligations`, and leaves
all ALT PASS flags and `storage_search_allowed` false.
