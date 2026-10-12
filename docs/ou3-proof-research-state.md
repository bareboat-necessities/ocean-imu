Warning: truncated output (original token count: 31631)
Total output lines: 1892

# OU-III proof research state

## Current hypothesis

**GOVERNING STRATEGY:** symmetry/kernels → justified quotient/tube → exact
linked information dissipation → actual forcing → only necessary activation
bounds. [The persistent strategy and derivation](ou3-information-shear-strategy.md)
must accompany future handoffs. This supersedes “enclose the entire nominal
state first” and “wait for planar admission before general dissipativity.”
The general theorem assumes existing MAGNETIC SERVICE on admitted histories;
all-time service for the special planar witness remains a separate OPEN lemma.
No general theorem substitutes for that witness's admission proof.

Verified merged main: `a47da994b85f2ec844b199d27b619a7fcc553adb`.
PR #653 is merged. This continuation uses the NEW branch
`proof/ou3-receipt-schur-continuation` from that exact main.
Concurrent main `59171c1a4aef654bb9808629e8856e4a56edfcc8` is incorporated,
including PR #663, version 2.4.0 and the automatic evidence bundle.
Its separate mechanical repair is `9b951b6ac`; current publication and CI
are PR metadata. Upstream calibration, BMM150 compensation and saved magnetic
source selection are preserved. No merge is authorized for the follow-up.
The concise current handoff is [ou3-ai-handoff.md](ou3-ai-handoff.md).

Recomputed manifest: **25 PROVED, 4 CONDITIONAL, 1 OPEN**;
`theorem_closed=false`, `regional_practical_stability_claimed=false`,
`certified_capture_time=null`. Existing proved algebra is retained without
promoting uniform coercivity, entry or admission.

The stability proof uses exact algebra/rationals, rigorous compact inequalities
after derivation, and identity/source regressions. No new replay, secant,
sampled grid or empirical factor is promoted to a proof result. Runtime
engineering regressions and ordinary-aiding source reproductions below are
separate from this proof path. Before substantial calculations state the
lemma/quantifiers, finite-error role, literal dissipation structure, required
proved/open bounds and avoidance of recorded dead ends. Required CI remains.

Shipping magnetic configuration: the OU-III AtomS3R sketch sets
`mag_continuous_hard_iron = false`. After refinement the applied hard-iron offset
and the magnetic reference stay at their refined values, so the MR1--MR2
continuous-reference transport and its MR9 proxy-tilt premise no longer arise on
the shipping path; the facade default and simulators keep the continuous path.
A refined-reference error from calibration residuals (axis gain, soft iron,
latency) must be charged to the stated magnetic residual envelopes; that fit is
not verified here and no constant is changed. The
provenance pin is updated (classification E: stale binding, not mathematics).

## Rapid-motion runtime repair scope

Runtime qualification failure **E**: the shared protection's pending-capture
latch could suppress aiding indefinitely during continued motion. The pinned
OU-III startup/device-order regression at large physical heel holds aiding for
411.58 seconds, applies no magnetic corrections in its final 30 seconds, and
reaches 93.38 degrees of attitude error and 139.51 m raw heave with the original
policy. The corrected execution, with magnetic heading feedback retained during holds,
ends at 3.99 degrees and bounds raw heave to 1.52 m. This invalidates the claim that the previous rapid-motion core
tests qualified continued Live operation with a pending capture. It does not
refute a regional stability theorem: the handheld large-heel execution is
outside its retained domain. The capture-only hold is now bounded to two seconds;
the request remains available for later qualified capture. Existing core,
covariance and local proof identities are retained. Current limiter: AtomS3R
bench confirmation. Next falsifiable experiment: repeat the same physical
tilt/turn sequence while recording `[ATT]`, magnetic service, heading and raw
heave; verify no ordinary-bias non-dynamic pending hold starves accepted corrections.

The first timeout trial also released a grossly corrupted qMEKF gyro bias before
its noise-qualified capture closed, leaving a 20.49-degree inherited-state tail
(**E**, invalid timeout scope). The existing 0.08 rad/s gross-bias threshold now
retains that protected repair window; the two-second bound applies to ordinary
bias estimates. No measurement/covariance or capture-quality threshold is relaxed.
PR #671 also conflated force-derived tilt/bias protection with magnetic
validity in all seven backends. Healthy time-aligned magnetic heading feedback
now remains active through a force hold, while tilt/bias protection remains.
The first unrestricted 3D-magnetic trial exposed 34.13 degrees of excursion in
a 720 deg/s case with an unmodelled 40 ms field delay (**E**, invalid timing
qualification). This invalidates using a turn-long magnetic veto as evidence
that an incorrectly timestamped field is harmless. Protected latency tests
use the existing calibrated delay and transport; ordinary-aiding unmodelled
latency controls remain visible. A subsequent OU-III inherited-state trial ended at 16.52 degrees because
magnetic navigation gain rows carried corrupt turn cross-covariance into
retained navigation means (**E**, heading-only repair was incomplete). Held
magnetic corrections now retain physical navigation means as well as tilt
and bias. For TFG, the invariant rho=[x]phi correction cancels the finite
rotation of the retained physical navigation vector; the modified gain is
used in Joseph covariance and reset. A new held-field/absent-field comparison
checks informative magnetic service, tilt preservation and bias preservation
across the production backends. Its initial arbitrary 20% two-second response
criterion rejected OU-II despite a measured 17.19-to-14.45-degree improvement
(**E**, observer-independent rate assumption). The service check instead
requires a resolved 0.01 rad improvement without prescribing a common rate.
The NLO scheduling-gap case initially peaked at 17.80 degrees because an
end-of-gap field was reused as feedback during every intervening substep
(**E**, magnetic-epoch mismatch). Held NLO propagation now applies the current
field only on the final bounded substep rather than inventing magnetic
observations through the unobserved gap. Healthy NLO heading feedback also
resolved a 90-degree error before capture qualification; the old exactly-one
capture assertion rejected an already accurate 0.01-degree heading (**E**,
obsolete control-flow expectation). Large-angle tails still must be below
three degrees with at most one capture, and the zero-innovation antipode still
requires one independent capture. The PII sketch's actual update-body replay retains healthy magnetic service
through a 137-degree start and a full protected turn; a genuine field outage
still expires service and fresh fields restore it. Its graphical timeout
previously drew `0 / N` with no visible invalid-status strip. The shared
renderer now suppresses the unavailable needle and shows `--- / WAIT`; actual
zero-degree headings remain valid. These are runtime/output repairs, not proof
results. Hardware qualification remains the limiter. Local full-suite validation also
exhausted the temporary filesystem while writing reproducible replay CSV
outputs (**E**, infrastructure); completed outputs can be discarded between
backend runs without changing any estimator or quality gate. Required CI
correctly rejects stale proof/publication source bindings after these source
changes. Those records need exact regeneration before this draft can be merged;
their older ordinary-aiding results do not certify the new held-magnetic path.

Current recovery hypothesis: qualified two-vector/gyro bias repair,
physical-coordinate TFG covariance reset, and explicit treatment of corrupted
navigation feedback close the finite inherited-state recovery cases. qMEKF's
quiet seven-second tails are below 0.00003 degrees; TFG yaw/combined tails are
0.34495/0.92293 degrees. Ten qMEKF and ten TFG continuations include sensor
noise, resumed vessel motion, sixty-second quiet tails and calibrated 40 ms
magnetic delay; all remain below
five degrees in the fresh optimized build (maxima 0.48445/2.92141 degrees).
All fourteen seven-backend late-protection cases pass. Aiding resumes and
healthy biases survive. The 336 ordinary causal controls change by at most
0.00079 degrees across the combined-source rebuild; their mechanism conclusions
are unchanged. The native OU-III rerun reproduces every numerical value and both
native byte hashes exactly.
Fresh full-matrix evidence and full-suite status are recorded in PR metadata.

Failed recovery approaches / DEAD_ENDS (all class E, engineering qualification):
- Attitude-only capture retained 16.54773/17.24282-degree qMEKF tails and a
  13.96484-degree TFG yaw tail. Artificial-angle injection did not qualify
  carried bias/covariance/navigation recovery.
- Gyro-bias repair alone left TFG yaw at 13.77418 degrees; physical covariance
  mapping alone left 16.45352 degrees. Operation tracing isolated S=0 feedback.
  TFG now masks its attitude/bias gain rows after gross repair, retaining
  navigation updates, until explicit estimator initialization.
- Endpoint differentiation rejected a noisy stationary bias candidate; a
  short-arc least-squares fit had 0.04747 rad/s error. Regression residuals now
  qualify slope noise before replacement; slow observed motion is subtracted.
  A fixed two-second delay gave a 6.99196-degree TFG yaw tail; noise-qualified
  delay alone gave 8.16057 degrees. The retained acceleration confidence mattered.
- Adding `(2g)^2 I` world-acceleration uncertainty passed seven-second/noisy
  recovery but left 12.09322 degrees at sixty seconds. Adding only acceleration-
  bias uncertainty still left 5.93086 degrees. Quiet force cannot identify tilt
  and bias independently; further covariance tuning is rejected. The recovery
  now explicitly returns only a gross learned acceleration bias outside 4% of g
  to the calibrated zero-residual prior, with matching uncertainty. This prior
  is the same calibration assumption required by gravity capture, not a new
  observability theorem. Navigation means are retained and the actual runtime
  reset, rather than an artificial replay-boundary reseed, is tested.
- The independent finite-difference covariance test's absolute 1e-8 cross-term
  tolerance failed after the large PSD uncertainty charge. Its tolerance now
  scales with covariance norm and finite-difference roundoff; mean preservation,
  physical marginals, selective S rows, initialization and PSD remain checked.

Retained facts: qMEKF's independent SO(3) geometry fixes, protected replay bounds
and OU-III ordinary-aiding identities remain valid. No admitted shipping theorem
counterexample was found and no theorem flag is promoted. Current limiter is
hardware/calibration qualification, not the reproduced finite quiet tails.
Next falsifiable qualification: timestamped turns with calibrated force residuals,
raw gyro, magnetic epochs and recovery state; test the explicit calibration prior
and retained-state behavior on device. Main `c9c41de` (effective magnetic-delay
calibration) and intervening PR documentation `8e3e5d8` are incorporated. AUX
status-bit masking and the required latest calibration format are preserved.
The source-wide native diagnostic is rerun on this combined source tree.
The final local `make all EIGEN_DIR=... SIM_DATA_ZIP=...` stopped compiling
`tests/imu_calibrate/accel_cal-test.o` with `error writing to /tmp/...s: No space
left on device`; the qMEKF Arduino compile and subsequent git fetch hit the
same storage exhaustion (E: infrastructure, not a code regression). Dispose
only of ignored generated CSV outputs, retain pinned inputs/source/evidence,
then retry unchanged code while cleaning completed directories' transient output.
Both AtomS3R retries pass: qMEKF uses 823747 flash / 116276 RAM bytes; TFG
uses 1069043 flash / 45028 RAM bytes. Full native-suite retry is in progress;
its outcome and the required CI gates remain PR metadata.
The combined-tree `make -C tests/validation test` stops at
`python3 ../../tools/ou_evidence_contract.py --auto`: validation replay provenance
mismatches `src/ahrs/Mahony_AHRS.h`,
`src/kalman_ou_ii/Kalman3D_Wave_OU_II.h`,
`src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h`,
`src/tuner/VerticalAccelComplementary.h` and `src/util/MagneticRotation.h`;
robustness provenance mismatches the Mahony, vertical-complementary and magnetic-
rotation paths (E: stale full-replay evidence). The new rapid-motion controls and
ordinary-aiding native proof diagnostics pass their own source bindings, but do
not refresh these separate statistical replays. The current merge limiter is
genuine full OU/TFG evidence regeneration and validation of its published tree;
the branch's automatic full-evidence CI is queued. Do not replace old replay
hashes without rerunning their simulations or weaken the publication gate.
The direct Python suite runs 1137 tests: four failures and two errors concern
those stale replay bundles; one CMake check is skipped because CMake is absent.
A fifth failure expected one synthetic stale-source error but saw two because
this ledger changed during the test (E: overlapping edit/test execution).
All three provenance tests pass when rerun on a stable tree. The interrupted
parallel-build attempt could schedule top-level tests before build completion;
its output is discarded. The retry keeps top-level build/test ordering and
parallelizes only recursive compilation, with no source or threshold change.
The unchanged evidence builder rejects stale source/record hashes after the
merge (E); refresh only the bindings backed by this new replay and the unchanged
legacy-operation audit. No numerical record or theorem gate is weakened.
Integration tests caught duplicate AUX cache storage from the automatic merge
and fixture omissions for the new timing configuration (E). The shared storage
is deduplicated; tests retain new main's status-bit masking and verify calibrated
delay delivery and configuration mismatch, alongside gyro/acceleration freshness.

Documentation validation exposed missing sandbox `mpmath` and `numpy`
dependencies (classification E). The workflow dependencies are installed and
all three provenance tests pass. The source-tree gate also rejected stale byte
bindings after comment/docstring edits (E); refreshed source and linked-record
hashes pass the unchanged gate. Executable C++ tokens and Python ASTs excluding
docstrings match the prior revision, so the existing numerical evidence and
open theorem flags are retained. This rebind reuses the existing native run;
it does not claim a new replay or alter a mathematical hypothesis. The next
physical qualification remains the timestamped hardware turn below.
The full native build passes; dataset verification initially failed because
Python 3.9 lacks `hashlib.file_digest` (E: sandbox runtime mismatch).
Python 3.11 verifies the pinned archive. The full test run then stops in
`./heave_only-sim`: Küchler Mahony-frontend Z RMS is 17.2632% > 16.99% for
JONSWAP Hs=8.5 and 20.2382% > 19.98% for PMSTOKES Hs=8.5 (E: numerical
regression-gate failure). Clean PR code reproduces these values under GCC
11.5/Eigen 3.3.9; documentation did not introduce them. Existing source/token
identities and theorem flags remain valid. Current limiter: heave regression
qualification in this toolchain; next experiment is comparison under the CI
toolchain with the same pinned inputs. No threshold or estimator is changed.

The mathematical audit is recorded in `docs/rapid-motion-math.md`. Two qMEKF
implementation errors are independently reproduced (classification E): the
quaternion half-angle propagator block was used as an SO(3) attitude-error
transition, and covariance was not transported after finite error injection.
The no-noise/no-observation spin test incorrectly contracted unit axial
variance to 0.3726473193. Central differences of the physical right-error flow
and reset reject the original header and pass the corrected implementation.
These qMEKF changes do not alter OU-III's ordinary-aiding equations.
The offset-turn flip remains without protection after those repairs, so the
hypothesis that they alone explain the common failure is rejected. Matched
known-force controls isolate the invalid gravity assumption, while exact
antipodal controls isolate the local zero-tangent correction equilibrium.
The physical hardware cause and conversion epochs remain unqualified.
The source-wide secant gate requires a rerun after the qMEKF change; both
native byte hashes and all numerical results reproduce exactly on that rerun.
`make all` then stopped in `python3 -m unittest discover -v -p 'test_*.py'`
with eight provenance-test failures caused by the two newly regenerated
record hashes (E). The actual rerun records are rebound, preserving all numerical
payloads and open theorem flags; the unchanged evidence builder passes.
An interrupted Arduino compile left a zero-byte NLO object, causing
`undefined reference to setup` and `loop` at link (E). Removing that generated
object and recompiling the unchanged source passes. A diagnosis run also found
an old executable after an overlapping header edit; its incomplete controls
are discarded. The diagnosis driver now freshly compiles its targets. AHRS
simulation dependencies explicitly include the qMEKF header; the current
geometry and wave suites pass after rebuilding. An extra compile with
`-Wall -Wextra -Wshadow -Wconversion -Werror` failed on pre-existing `R`/`M`
shadowing in qMEKF (E: stricter diagnostic configuration). The repository's
`-Wall -Wextra -Wpedantic -Werror` flags pass without unrelated source changes.
All 336 matched causal controls and the seven-backend protection/recovery matrix
complete from fresh builds. Every exact-antipode local control stays at 180
degrees; each finite-angle recovery escapes with one capture. Full validation
status belongs in PR metadata. The next physical qualification is a timestamped
hardware turn. No admitted theorem counterexample is claimed; these runtime
cases remain finite engineering evidence.

The device motion guard, selective gain-row holds and settled attitude
reacquisition are new runtime operations. The existing local identities remain
applicable to the ordinary aiding branch only; no capture or stability theorem
is transferred to the guarded or recovery branch. Skipped corrections do not
establish MAGNETIC SERVICE. Covariance recovery clears attitude cross terms
while preserving navigation means and the established world reference. Bias
means are preserved except for the explicitly qualified qMEKF/TFG repair above.
Structures preserved: literal gyro propagation, measured specific force for
translation, OU/S and tuner chronology, magnetic calibration and frame.
Relaxations introduced: rapid-motion host tests use synthetic pivot-motion
histories, with ideal and seeded noisy sensor cases; they are engineering
regressions, not admitted theorem histories. There is no hardware qualification.
Classification E: the initial new replay build omitted the acquisition header
and referred to a nonexistent TFG proxy-initialization accessor. This invalidated
only that test build, retained all existing mathematical results, and was repaired
by including the production acquisition type and using the facade's tilt accessor.
The first protected qMEKF ramp reached 14.2504 degrees at 180 deg/s and a
15 cm pivot: classification E, incomplete motion detection. The rate-only
threshold fired after tangential acceleration had already entered the filter.
A 30 ms filtered angular-acceleration gate closes that interval; the same
unchanged 12-degree excursion bound is rerun. The next falsifiable check is the
common seven-backend motion replay followed by `make all`.
A second, distinct engineering failure was a 17.4363-degree OU-III excursion
2.886 s after the 720 deg/s turn ended. Literal-operation tracing located the
rotation in S=0 corrections after ordinary aiding resumed, while physical
acceleration was already quiet. A trial S-residual threshold still gave 17.1241 degrees (E); residual magnitude
alone did not bound the attitude cross-covariance effect. That trial is removed.
The architecture review instead keeps the selective hold until the independent
settled reacquisition window closes after **every** dynamic interval, even with
small heading error. Its attitude covariance reset removes the stale cross
terms before S can affect attitude again; navigation means remain intact. The
unchanged seven-backend excursion and recovery bounds passed for OU-III.
A vessel need not stop completely: extra settling hold is bounded at two seconds
unless a large error is pending. At re-entry all three Kalman cores remove
cross-correlations between held attitude/bias and driven navigation blocks,
preserving means, marginal covariances and PSD. The moving-after-turn case checks
that normal vessel motion regains aiding without a stationary reset.
Classification E also covers the initial test linkage omissions (`g_std`, C++
link driver), the parallel compiler OOM, and a toolchain archive checksum failure;
none is evidence about filter stability. Rebuild with less concurrency and
verify the archive against its pinned checksum before retrying. The new OU-II
covariance test initially called an OU-III-only integral-displacement accessor;
the build failure invalidated only that test. OU-II's navigation block has nine
states, not OU-III's twelve; its projection bounds and preservation checks use
that layout. A mixed-version M5GFX installation failed compilation and was
replaced by a clean installation of the CI-pinned dependency. These are E.
Compiler diagnostics also exposed a dangling Eigen expression in the new
noise-generator lambda. Its return now materializes a vector; earlier noisy
case results are discarded and all seven are rerun with the corrected generator.
The AW source diagnostic stopped with `shipping accelerometer anchor changed`
after adding the freshness guard (E: stale instrumentation anchor). Its tap now
retains that guard and records only a fresh correction attempt; unique-anchor
and untapped-control parity checks remain mandatory. The full evidence contract
also rejects the old replay/source hashes after these implementation edits;
full validation, robustness and TFG replays regenerate their source provenance.
The 1200-second compatibility probe's `g++ -O1` compilation exceeded its
90-second subprocess timeout under concurrent validation load (E). No replay
or mathematical comparison ran in that attempt; repeat the unchanged probe
after reducing compiler contention, retaining its existing parity checks.
The first full replay regeneration stalled in the rolling transition-plot
convolution with nine BLAS threads competing for eight CPU cores (E: execution
configuration). The same 240001-by-12000 convolution took 0.94 seconds with
one BLAS thread. Regeneration uses one numerical-library thread and six
independent simulator jobs; its protocol, seeds and acceptance criteria are
unchanged. Incomplete first-attempt outputs are not used as replay evidence.
The complete validation, robustness and TFG studies now reproduce locally.
`python3 tools/tfg_comparison.py --check --require-ci --publish` initially failed
with `ValueError: publication requires CI replay provenance` (E: missing CI
execution). The same bundle passes `--check`; no workflow identity is invented,
and no publication gate is weakened. PR #671's initial `evidence-contract` and
`ou-validation` jobs fail at this same command. Full branch replay
`38176210485` supplies the genuine TFG artifact, with source and byte hashes
verified; that artifact passes the unchanged `--check --require-ci` gate.
Full OU CI regeneration and final-head validation remain required.
The earlier `make all` passes the preceding
C++/simulation suites and stops at the same gate, invoked as
`python3 ../../tools/tfg_comparison.py --check --require-ci` in `tests/validation`.
The complete Python suite and the two downstream wave suites pass separately
(one CMake configuration check is skipped because CMake is unavailable).
The diagnostic-sketch host test first
failed extraction with `Missing definition: onSample` because the helper allowed
`static inline` but not `static` alone (E). The helper now accepts both valid
declarations; the actual extracted integration/freshness code passes its gap,
NaN and repeated-frame cases. No estimator inference depended on the failed
test invocation.
Source-bound native diagnostics were rerun with their untapped control and
carried-state checks. The signed-balance, construction-mean and planar-service
exports reproduced the committed trace bytes exactly. Their unchanged
deterministic high-precision reductions are retained on those identical inputs;
the native source fingerprints come from the new executions. This is not a
relabelled historical replay or qualification of the new guarded branch.
The broader Python suite exposed an obsolete unconditional legacy-operation
projection and a publication table still containing 18.371 instead of the
fresh replay's 18.372 (both E). The source audit now explicitly specializes
the always-enabled aiding/fresh-acceleration, no-external-reacquisition domain;
all three original token hashes remain unchanged. Mutation tests reject altered
defaults or ordinary operations. Its receipt binds the freshly regenerated
native records, without claiming that every numeric payload stayed unchanged.
The mutation test initially expected a token-hash mismatch for an edit that the
stricter unique-anchor check already rejected (E: test expectation). It now
asserts both rejection paths; neither audit rule was loosened.
After the table matched, the same publication check exposed stale generated
macros (including 12.523 versus 12.525, E: incomplete publication mirroring).
Both the table and macro definitions are copied from the new evidence. These repairs retain the original
ordinary-branch identities and do not qualify guarded/recovery execution.
PR review exposed two further E implementation failures: capture used held
acceleration without an age bound, and standalone backend scripts could hide a
failed first rapid-motion test because explicit `set -e` came later or was absent.
The qMEKF regression reproduced `missing acceleration caused global recovery`.
The shared policy now receives actual accelerometer freshness, holds aiding after
20 ms without a sample, and clears capture across any such gap. A valid 100 Hz
acceleration stream may be out of phase with 25 Hz magnetic observations. The
same three dropout/resumption cases run on all seven backends, including clock
wrap and a last contaminated force with gravity's norm. The initial intermittent
test allowed capture before the dropout began (E: test chronology); it now starts
the gap before the next magnetic frame. All six scripts enable errexit before
their first executable test; injected exit-23 failures verify propagation.
All seven backends pass their 22 rapid-motion, 14 large-angle and three dropout
cases after the repair. The dropout recovery tail is below 0.23 degrees in every
case, against the unchanged 3-degree bound. Full `make all` and firmware builds
are rerun on the repaired source; their completion remains a validation gate.
These repairs retain navigation/bias state and gyro propagation. No proof
identity or bound changes; only guarded-path engineering claims are involved.
The source-tree gate then reports `finite quotient/secant evidence verification
failed: shipping source tree changed` (E: stale diagnostic binding). Both native
secant runs were repeated on the new source tree with all fork/forcing controls;
their native byte hashes and numerical results reproduce exactly. The new
record and its receipt retain the finite-diagnostic scope and all open premises.
No shipping theorem counterexample is asserted. Existing algebraic identities
retain their original hypotheses; guarded intervals, covariance decorrelation
and global reacquisition need separate analysis. The next falsifiable runtime
check is the same replay on each final binary and an instrumented hardware turn
with fresh-sample, elapsed-time, magnetic-age and `[ATT]` diagnostics.

## Magnetic-input timing scope

The calibration-only follow-up exposed two binding mismatches in
`python3 -m tools.stability.ou3_theorem.build_evidence`: the edited ledger blob
and the native-secant whole-source digest (**E: implementation/provenance**).
The recorded digest matches immutable main `297a211c4105d0ae95c5f31f66d9d2c488268df4`.
All 36 shipping files in the actual native-probe compiler include closure are
byte-identical to that main; only the wizard and timing fitter changed outside
that closure. The existing binding receipt records the closure, old/new digests
and the baseline record hash with unchanged numerical payload. No replay, numerical result, acceptance
gate or qualification flag is changed by this binding refresh. The invalidated
hypothesis is that whole-tree identity survives calibration-only edits; retained
facts are the unchanged probe/dependencies and finite payload. Hardware timing
qualification remains the limiter. Next falsifiable check: rerun the unchanged
source/evidence validators and device capture diagnostics.

Calibration usability regressions (**E: implementation/diagnostic failure**)
reject the old hypothesis that scalar delay identification needs rank-two gyro
covariance: variable-rate single-axis rocking identifies an injected 40 ms lag.
A separate 180 ms pair needs 260 ms history at the largest candidate delay,
exceeding the former 250 ms query budget. Native noisy 100 Hz rocking/pauses and
180 ms pair regressions now pass; stillness, constant rate, field-axis rotation
and held-out disagreement remain rejected. Structures preserved: runtime
estimator/tuner/covariance, transport, v5 binding and all residual/profile/
held-out limits. Relaxations introduced: none to the stability argument;
synthetic tests remain finite diagnostics. No stability counterexample or
promoted theorem follows. Current limiter: no device trace identifies the
reported generic failure's exact gate. Next falsifiable experiment: capture
on-device reason/pair/gap diagnostics and verify a saved fit across reboot/load.

Calibration now identifies an effective magnetic delay relative to the delivered
gyro clock and stores it in the mandatory v5 record. Exact physical conversion
timing remains OPEN. CI rejected stale sketch/whole-source fingerprints after
this change (**E: implementation/diagnostic failure**); the hypothesis that the
old bindings described the current tree was invalid. Both complete native
secant replays reproduce identical raw bytes; original finite numerical payloads
and all non-promotion flags are retained. No mathematical counterexample or new
stability result follows. Structures preserved: the legacy core, covariance,
tuner, OU/S, resets and carried native histories. Relaxations introduced: none
beyond the existing finite off-orbit diagnostic. The limiter remains the timed
path's physical timestamp/service qualification; the next falsifiable check is
device timing under reboot and display/USB load, with applied-service records.

The AtomS3R OU-III frontend now transports known host-cache delay and suppresses
repeated host observations. The BMM150 conversion epoch/aperture is not known;
physical timestamp qualification remains OPEN. The 21-state core, Jacobian,
OU/S/SpectralMSE/tuner/BA dynamics and injection/reset are retained. Reference
and hard-iron accumulation use saved exogenous proxy snapshots. Timed MEKF
input and covariance depend on carried gyro-bias estimates/covariances: LR1's
zero corrected-input reverse port applies to the legacy API only. The timed
ports, physical latency envelopes and source-uniform actually applied ser…15631 tokens truncated…Hbar*y=L_pmu'*Dbar*mu`, `mu=L_pmu*y`, and
  `y'*Hbar*y=mu'*Dbar*mu>0`. Its first process loss is positive;
  no zero-loss condition on the measurement columns follows. These losses
  are already included in Hbar. Treating them as an additional reserve
  would double-charge the original gap, and qualitative service/kernel
  facts alone did not prove the required quantitative comparison.
  No actual feasible threshold mode was exhibited or excluded, and no
  shipping sign or impossibility is asserted. Retained: CR39's exact
  signed matrix and CR41--CR43's full cross cost. The ONE blocker remains
  a common positive-reserve sign for feasible `R_W(c)`; its receipt corner
  is the first unpaid part. No new completion, bound, storage or flag.
  Next falsifiable calculation: prove that the SAME shipping covariance
  action excludes every reconstructed feasible CR30 mode with `chi>=1`,
  then pay the actual CR43 cross subtraction. Another zero-action argument
  or conditional Schur identity cannot settle that comparison.
* **D, signed near-kernel implies zero-action shortcut:** exact failed
  inference: `y'(H-L'DL)y -> 0` implies every positive process/service
  loss tends to zero. The actual limiting equation is CR34,
  `y'H y=mu'D mu>0`; the stationary balance mode is CR35.
  Retained: all exact CR25--CR28 accounting, constrained equations,
  integrated same-tangent and directional service losses. Limiter: exclude
  or pay the feasible receipt-driven eigenmode at threshold one. Next:
  insert actual image/cone and linked service restrictions in CR30.
  No independent norms, tiny kernel reserve, toy witness or new storage.
* **E, local tooling:** corrected guessed appendix/optional filenames and
  a JavaScript call delimiter before edits; an unsupported encoded branch
  read did not replace the successful PR/head fetch. Missing local mpmath
  and ruff were installed in scratch for the unchanged required checks.
  The scratch appendix wrapper first lacked the article's norm macro;
  importing the existing macro let the full appendix compile unchanged.
  Native git push lacked a GitHub credential; publish the identical
  validated tree through the connected API with an expected-head check.
  Main advanced through PR #655. Two thin fetches reported unresolved
  deltas; missing repository blobs were restored through hash-verified API
  reads. Ledger/provenance conflicts were reconciled, preserving all proof
  work and the upstream device configuration. The regular planar equations
  and margins are unchanged; no shipping configuration was authored here.
  Full Schur continuation: guessed optional module paths were corrected
  to tools/stability/ou3_theorem. A render redirect found no scratch
  directory; creating that directory and rerunning rendered the appendix.
  Publication verification briefly read the remote ref before its fetch
  finished. After waiting for that fetch, the fetched tree matched the
  already validated/API-verified tree exactly; no branch update was forced.
  No mathematical premise or validation gate changed.
  Receipt-corner continuation: a guessed KalmanOU*.h filename matched no
  OU-III file; rg --files identified Kalman3D_Wave_OU_III.h and the actual
  floor, integrated prediction, S and accelerometer operations were read.
  No mathematical conclusion was taken from the failed read.
  An optional provenance.py path was absent; the actual build_evidence.py
  validator was used. The first `make -C tests/validation evidence-test`
  ran 1122 tests and failed eight provenance assertions because the edited
  ledger's bound blob hash had not yet been refreshed. The unchanged
  fail-closed gate detected that stale binding; refresh it and rerun the
  affected provenance tests. This is E, not a mathematical sign failure.
  The first refresh command used the test-directory cwd for repository-relative
  paths; it failed before writing. Its consequent test rerun again rejected
  the stale hash. Use the absolute repository root for the refresh/validator, and
  launch the focused test only after successful validation.
  Merge-CI continuation: the four failed merge jobs rejected the old finite
  secant whole-src binding after main's two AtomS3R header additions.
  Compiler dependency inspection proves those headers are not read by that
  probe; all 34 shipping dependencies are byte-identical. Rebind with an
  explicit retained original digest and source-identity audit, not another
  secant/replay run or a weakened validator. A guessed
  planar_native_secants.cpp read was corrected to the driver's literal
  planar_native_secant_probe.cpp path. No mathematical premise was inferred
  from the absent path. A later guessed docs/ou3-theorem-status.md read was
  corrected with rg --files to reports/results/ou3_stability/theorem-status.json.
  The optional repository dataset archive was absent; the unchanged required
  make all target remains responsible for fetching and verifying that input.
  Main's subsequent PR #658 is incorporated without discarding either
  history; its calibration/runtime edits do not enter the native probe's
  shipping include closure. The previous source-refresh record is retained.
  The calibration rebuild initially repeated a relative -C path inside that
  same test-directory cwd; use make build there instead. The failed invocation
  performed no build or test, and no validation gate was changed.
  The first make all session later became unavailable with no final exit
  status; its log stops during OU-III local-gravity compilation after the
  primary build completed. This is incomplete local validation, not a test
  pass or a mathematical failure. Rerun the unchanged target using its
  existing build artifacts and require a recorded exit status. That repeated
  make all completed with recorded exit status 0, including 1122 evidence
  tests (one existing skip). No test or validation gate was changed.
  A publication preparation found its cached base-file contents unavailable
  and stopped before creating a tree or updating the remote branch. Re-read
  immutable API file contents and require native/API tree equality before
  publishing the validation record.

* **D, process-short-only sign inference:** CR26's positive effective gap
  does not prove the complete nonzero receipt remainder nonnegative.
  Exact unpaid quantity: CR28, equivalently CR17a, on the actual inherited
  image/cones. Invalid inference: nonsingular receipt weights or positive
  conditional process action imply affordable regression cross work.
  Retained: CR25's source-uniform short, CR26's whole signed inequality,
  full prefix polarization and nonzero scalar zero-gap branches.
  Limiter: actual same-root residual process/acc/S/magnetic covariance-column
  payment and inherited feasibility. Next falsifiable calculation: evaluate
  the linked conditional sign and remaining reader on the literal generator,
  then its full endogenous Schur cross cost. No scalar tightening or replay.
* **E, nonzero-receipt tooling:** a read used a nonexistent guessed lin_source
  path; existing lin_path_certificate was read instead. Two focused runs
  exposed an existing strict-active error-message regression; the original
  wording was restored and all 43 focused tests pass. The first appendix
  render failed at CR25 with "Display math should end with $$"; its missing
  closing display delimiter was restored before final validation. The retained-square
  argument and validation thresholds were unchanged. Guessed optional render
  paths (doc/kalman_ou_iii/Makefile, tools/render* and tools/*latex*) were
  absent; the existing render workflow was read directly. No source premise
  was inferred from an absent path.

* **D, direct AW observation as automatic uniform floor absorption:** CR20
  pays eta>Cplus/s, but the full remaining reader is
  b+(eta/omega)t, with t generated by the SAME regression tangent.
  Invalidated implication: a positive local paid fraction proves the
  uniform complete-word covariance margin. Retained: CR18--CR21 exact
  shipping correction/floor payment and CR23 qualified signed reserve.
  Limiter: CR22's actual covariance service/readers and full CR17 cross
  remainder, plus inherited origin/rank. No admitted negative word or
  instability is established. Next: use CR24's full integrated generator
  and actual S/magnetic Dh columns to cross CR22's threshold, without
  independent maximization of eta, t or b.
* **E, continuation tooling:** a read requested nonexistent guessed status/
  manifest module paths; exact current paths were then resolved. TeX quote
  delimiters and mismatched document patch contexts failed before those edits;
  exact contexts and ordinary emphasis fixed the literals. A provenance sync
  initially replaced a summary audit; its existing fields were restored before
  validation, with only the new scoped result updated. Missing mpmath/ruff were
  installed in a scratch dependency directory; existing runtime Eigen is
  supplied through EIGEN_INCLUDE_DIR. No mathematical result or gate changed.
  The exact identities and fail-closed validation remain controlling.

* **D, sharp-reader payment with released coordinates:** (CR9)'s attempted
  Schur remainder has exact relative charge greater than one in (CR10)
  for nonzero AW cross regression. The invalidated implication is that
  paying each actual information-decrement reader independently suffices
  to pay the face. Retained: sharp same-correction Fisher reader, literal
  AA receipt, actual S/magnetic conditional equations and all integrated
  process action. Limiter: their shared covariance tangent must satisfy
  the uniform full linked Schur condition on the actual inherited image.
  Next: use (CR5),(CR7),(CR8) together, without releasing dC/dB/dT or
  classifying the endogenous reader as external supply.
* **E, read/edit tooling:** guessed proof/source paths and rejected patch
  invocations failed before modification (context mismatch and an
  unescaped Markdown delimiter in JavaScript). Exact paths/contexts were
  inspected and corrected; a trailing TeX space caught by diff validation
  was removed. The full evidence suite's architecture guard caught a
  covariance-test variable name colliding with a retired architecture
  token; it was renamed to describe prediction, with the guard unchanged.
  The publication payload gate compared Python Unicode scalar counts with
  JavaScript UTF-16 lengths; it was corrected to compare scalar counts.
  The final Git tree equality check remains mandatory and unchanged.
  No mathematical premise or test was changed to bypass a failure.
  Current shared-tangent derivation: guessed module paths and a TeX
  dollar/brace sequence interpreted as JavaScript interpolation failed
  before edits; exact paths and literal TeX were corrected. A proposed
  symbolic simplification could not import sympy in either Python runtime;
  it was not used as evidence or replaced by a numerical search. The
  exact-rational CR15--CR17a calculations and unchanged gates remain.
* **E, receipt-kernel certificate scope:** the initial CR14/CR16 statement
  did not explicitly require zero target-relative tangent at zero-gap
  faces, although the prefix helper rejects the unsupported crossing.
  Invalidated implication: zero strictly-active receipts automatically
  qualify a zero-gap event. Repaired using the literal directional map
  and the SAME zero-gap constraints in all partial matrices. The reserve
  is unchanged on that qualified intersection; nonzero crossing directions
  and its actual rank remain open. No physical or shipping assumption is
  strengthened. Next: qualify the actual directional history image before
  testing its receipt/full Schur sign; no cross direction is deleted.

* **D, complete causal-work domination still open:** the source audit removes
  homogeneous reference feedback and proves a non-circular reference amplitude
  bound. It does not bound the actual MEKF mean/covariance feedback, inherited
  comparison energy, or terminal source/gauge quadratic in (LR8). The target
  remains a positive uniform gap after those terms. Retain the exact zero
  port and paired-source cancellation; do not promote the amplitude bound to
  reference alignment, low NIS, or complete storage contraction.
* **E, local verification environment:** missing ruff/mpmath packages prevented
  the first lint/evidence invocations; installed the required dependencies and
  reran unchanged gates. Two rejected patch contexts changed no files and were
  corrected. No mathematical premise or test was weakened.
* **E, publication credentials:** `git push origin HEAD:refs/heads/proof/ou3-planar-service-phase`
  failed with `could not read Username for 'https://github.com': No such device or address`.
  Published the identical validated tree through the connected repository API
  using a non-forced expected-head update. Preserved the intervening evidence
  commit; no proof premise, shipping equation or check changed. Full local
  evidence contracts pass: 1,107 tests, two existing skips. Current-head CI
  must be inspected separately; neither a parent pass nor local validation is
  promoted to current-head CI success.

* **D, bordered-storage uniform-coercivity failure:** the new Fisher lift pays
  covered fixed-operand score work but accumulates `B_W` in its slack and can
  lose all uniform mean weight. A per-word budget cap does not bound its
  cumulative sum. The precise unpaid alternatives are the signed restoration
  charge (BC8) and actual causal work (BC4),(BC7), or a different uniformly
  coercive storage. No constant-source/covariance assumption is imposed on
  shipping. Retain the exact lift and two nonnegative correction losses;
  reject nonexpansion alone as physical comparison retention.
  The one permitted refinement now retains the exact restoration and restricts
  tangents to a causal lift. (BR1) recovers the original unpaid gap, and the
  propagated-root scope test refutes automatic coercivity from causality.
  Stop this route as a shortcut: another scalar reset, accumulated-budget cap,
  or isolated coefficient tightening supplies no new mechanism. Retain the
  identities; the next falsifiable obligation is a quantitative shipping
  bound on `A_B-F_B-R_W` or the actual-lift Gram inequality (BR2), including
  reference/source coordinates and endpoint gauge energy. Both remain OPEN.
* **E, storage-test lint:** an unused innovation binding was removed after
  ruff flagged it; algebra tests passed and no check or condition was weakened.
* **E, storage manuscript:** the standalone compile rejected an unavailable
  script-font command; standard mathcal notation replaces it, without adding
  a rendering dependency or changing the proof.

* **D, recurring comparison/work closure still not obtained:** composing the
  actual supported process balance with the applied acc/S/magnetic comparison
  identities gives exactly (HB7), not a new uniform bound. The unpaid quantity
  is `V_X,0-V_X,N+process_source_action+correction_defect_action-NIS+
  uncovered_signed_change`. The negative innovation and terminal terms must
  stay paired with their source terms; dropping them does not solve the root
  comparison problem. At S updates the defect remains `-S_physical`.
  Separately, (CP7) leaves the complete causal remainder in one signed
  cross-plus-square expression. Positive (CP5) does not bound its coefficient
  maps or the inherited comparison. This attempt discharges neither obligation
  and adds no theorem or count. First-Live, held LIN BIBO and magnetic physical
  probe service do not alone supply the missing recurring comparison-to-action
  estimate. Retained: the exact balances and all source/covariance cancellation.
  Architecture review is required before another refinement of this mechanism:
  exhibit a non-circular finite-error comparison retention estimate on the
  required loss projection, or a storage formulation that controls the actual
  endogenous ports without presupposing that estimate. Another boundary
  telescope, smaller auxiliary coefficient or independent-maxima enclosure is
  not a discharge. No radius, new numerical exploration or admitted
  counterexample results from this attempt.

* **E, manuscript tooling:** the new section had a missing EP4 closing delimiter
  and tags inside EP6/EP7 split environments. Standalone compilation caught
  these; delimiter/tag placement was repaired without changing the argument
  or any check. The inherited E_L/Schur-work limiter is unchanged;
  validation is a compile and source-evidence regeneration, not a new bound.

* **E, source-calculus shortcut rejected before promotion:** planar commuting
  rotations do not make the literal nominal pitch transition exactly equal
  to its covariance transition. Normalizing the implemented small-angle
  quaternion polynomial leaves the explicit (EP4)--(EP5) derivative defect.
  Retain it; the odd process block and causal source variations are not zero.

* **D, uniform budget closure not established:** physical process action and
  the first-Live bound do not supply a recurring bound on E_L in the exact
  B_W decomposition. Retained: signed cross cancellation and uniform held
  boundary-score bound. Limiter: inherited active comparison-loss projection
  and actual auxiliary/source work. Next is a non-circular bound on those
  quantities, not a radius or independent port maxima.
* **E, regression call:** the new effective-noise test initially omitted the
  gain differential's inverse-innovation operand. Corrected to the literal
  seven-operand identity; no proof condition was weakened.
* **E, regression path:** the full suite runs from tests/validation, exposing
  a working-directory-relative source path in the new test. It now resolves
  from the test file's repository root; no source-faithfulness check removed.

* **E, formatting gate:** `git diff --check` caught trailing whitespace in the
  new endpoint proof. Removed; no mathematical premise or check was weakened.

* **E, validation synchronization:** a certificate scope clarification was
  written while the full test process still held the older module in memory.
  Its fail-closed provenance tests rejected the mismatch. Freeze the source
  and certificate together before rerunning; no mathematical condition changes.

* **E, local test dependency:** the full evidence test initially could not find
  Eigen. The installed runtime Eigen 3.4.0 is supplied through the existing
  `EIGEN_INCLUDE_DIR` option; no source, test or proof condition is weakened.

* **E, corrected proof bookkeeping:** draft reader windows ending at the current
  correction could count that upcoming operation. They now start at prior IMU
  roots with an explicit `hmax` offset; separation, age and process charges are
  recomputed exactly. The strict covariance bounds survive. This is neither
  an admitted counterexample nor a weakened service requirement.

* **E, regression tooling:** the new word module had one unused import; removed
  after the lint gate identified it. No mathematical condition changed.

* **E, manuscript tooling:** the proof appendix exceeded IEEEtran's 26-letter
  subsection counter (`Counter too large` at the final implication). Scoped
  numeric appendix reference and display numbering removes that rendering
  limit (IEEEtran uses separate `thesubsection` and `thesubsectiondis` hooks); no
  mathematical condition is changed. The remaining proof limiter is the
  inherited prefix domain and positive complete-word gap.
* **E, corrected algebra:** the partial joint bound must retain the square
  gain-residual charge `NIS*L_P/2`. Its covariance loss coefficient is
  `(lambda-NIS)/2`, not `lambda/2`. An exact regression rejects the omitted
  square. This correction supplies no uniform bound on shipping NIS.
* **D, current incomplete sufficient inequality:** the exact linked storage
  identity and residual-loss bound do not provide a positive uniform epsilon
  without inherited nominal/innovation/precision bounds. The invalid shortcut
  would be to substitute the private frontend decrement or a finite factor gain.
  Retained: source-derived qualified frontend cell, real reference cone, exact
  residual-loss and complete-word identities. Next: bound the remaining linked
  MEKF domain and word terms, before attempting any radius or floor.

* **E:** physical acceleration substituted for nominal aw, point innovation floor
  used as a local proof constant, omitted dR/rotation-block terms. Corrected
  cell-dependent inequalities and complete linked identities replace them.
* **E tooling:** the new normalization regression initially matched the wrong
  source variable spelling; the test now binds the literal expression. No
  mathematical condition was changed.
* **D:** the fixed-root covariance ball requires radius >60 while its SPD domain
  requires radius <1. Do not retry it or subdivide that failed representation.
* **D:** per-operation/grouped AW norm products lost complete-word cancellation.
  Do not revive that route.
* **E:** F/(I-KH)/G is not the complete mean derivative; private Mahony gain
  .582705764 is not rho_mean. Native binary32 step halving did not stabilize
  the small mean-to-P block. No more secant sweeps.
* **E:** wrong-sign gauge, unscaled projector checks, AW Loewner monotonicity,
  scalar-radius promotion and orientation-free matrix ceiling were corrected.
* **D, analytical relaxation only:** if a cell permits all magnetic corrections
  to be rejected, its information lower bound is exactly zero. This does not
  prove rejection or exclusion for the reached planar family.
* Two-/three-epoch physical BA separation is refuted by exact compatibility.
* Full local TeX rendering lacks `newtxtext.sty`. A minimal 27-subsection
  regression with the CTAN IEEEtran class passes both numbering hooks; use the
  full manuscript-render CI for the actual document.

## Retained facts

The H18-to-A21 magnetic release gate counts only corrections the core
actually applied (`seastate::common::MagServiceCount`); rejected calls advance
neither the count nor its 1-s clock. This is the accepted-update count already
charged in [the release bound](ou3-h18-release-compactness.md); executions
whose forwarded magnetic calls are all applied are unchanged. For
x=h/tau<1 the integrated-OU Q and phi_pa/phi_Sa are Taylor series truncated
below 2^-60 relative (in float the old closed forms were wrong by orders
of magnitude in the S rows for x in [.01,1)); `series_branch_defect` proves
the literal defect lies inside the charged `small_x_source_defect` allowance
for every x<1, so no certificate value changes. Invalid initializer inputs return false instead of throwing.
All other shipping source, behavior, constants, gates and physical
assumptions are unchanged. Retain complete state/P/K, OU, BA/BG forcing/projection, actual
acc/mag/S correction, Joseph/reset/AW, Mahony/guard/frequency/variance, coupled
adaptation, reference acquisition/refinement and both same-history clocks.
The local qualified sqrt(V)<=.15 theorem and all existing proved components
remain. No genuine admitted shipping counterexample has been established.

## Alternatives

All-time admission of this family would rule out universal eventual retention
in the specified absolute V<.0225 ball on that profile and require a
compatibility quotient/tube for that target. It would not prove divergence or
exclude every larger absolute practical bound. Analytical service
failure would exclude the pair. Neither branch is selected from replay.
STILL retains compatible-tube boundedness with the existing centered-window
persistent SLOW+FAST source charge, without a restart at regime transitions.

## Next falsifiable experiment

Decide the actual feasible receipt corner S_c in CR41 with the remaining
linked action in the SAME Hbar. If positive, pay the full endogenous
T_c-B_c'*S_c^-1*B_c immediately. Equivalently decide the one combined
R_W(c)>=0 on the justified full inherited image for one common c>0.
CR42 retains the actual cone constraints when eliminating the kernel;
its unconstrained residual cannot be silently treated as reached.
No actual positive or negative receipt eigenmode has been established.
Another completion, positive-action statement or zero-kernel argument
without a sign for this remainder is not a next step. No new storage,
reset, proof frame, compact retention premise or empirical search.

Attack CR28/CR17a on feasible NONZERO receipts, using CR25--CR27's
joint process/face inequality and actual scalar directional cones.
Keep the completed receipt square, full residual process and directional
service loss on one inherited tangent. The effective gap is finite; the
signed receipt cross cost still has no positive uniform certificate.
Use CR22's actual weighted reader and CR24's integrated process
generator and the directional S/magnetic covariance losses. CR18's
regression alignment is regenerated by g and Q_ao, not a free independent
parameter. Pay the surviving T dB T' and linked dT B(h_o+T)' on ONE
inherited tangent, retaining the full cross remainder CR17. A uniform
innovation ratio alone would not pay that cross cost. No source completion
or activation promotion follows from CR20's positive fraction.

Use LR8/LR9's correlated homogeneous block (CR2), keeping E_L,
S_q and C_W together in (OF1a). Establish only the needed
construction-to-root image/rank facts in (CR12); propagate its literal
prefix generator into (CR13a). The actual receipt rows and service Grams
must satisfy (CR17a)'s signed threshold, where its inverse is justified,
and the full cross remainder (CR17). The proved nuisance receipt-kernel
reserve is already substituted; do not count its process loss again.
Do not repeat the released-coordinate bound (CR9): its sharp reader
already fails structurally. Do not infer a margin from precision
monotonicity, a fixed AA tangent or the formal failed word. The inherited
endogenous image and admissible physical source-origin map require proof;
neither can be replaced by arbitrary independent auxiliary noise.

Use the held-score boundary cancellation only before release. The independent
E_L cap remains an OPEN sufficient-budget route, not a mandatory first step
for the full coupled gap. Neither the fixed half-charge test nor base c_aw
is a net-margin milestone. Do not independently maximize ports or bound
source tangents by physical time-rate limits. No compact enclosure without
a derived domain and affordable threshold; no radius before the full gap.
Use the now-proved locked-Live causal split: no homogeneous reverse reference
charge on that lift, and the paired source reader (LR6)--(LR7) in the full
remaining work. Reacquisition and source variations do not inherit a zero port.

## Verification scope

Exact rational regressions cover variational NIS, coupled correction/prediction
identities, held-BA Joseph structure, first-boundary bounds and fail-closed
domain checks. The source-bound native diagnostics are reproduced against
the ordinary-aiding branch after the runtime changes above. Their local
identities do not cover the new motion guard or recovery transitions. Rapid
motion replays are engineering evidence only. Current test/CI results are
recorded in the PR metadata.

## Shipping-faithfulness handoff

1. Preserved: all 21 states, actual P/gain/Joseph/masks, OU/BG/BA and
   AW--v--p--S, applied acc/S/magnetic rows, inherited causal auxiliary/source
   history, staged-tuner/queued-target lag, reset/projection and endpoint gauge.
2. Relaxations: existing regular real-operation qualified planar covariance
   fibre; new nonzero scalar crossings use the actual directional cones.
   Previously retained kernel results keep their nonempty-rank qualification.
   (CR9) is a rejected relaxation that
   explicitly releases shared tangent coordinates.
   No source-origin graph, reachability, service coverage or storage is invented.
3. Failures: D for (CR9)--(CR10), automatic acc/floor uniform payment and
   process-short-only inference of complete receipt positivity and
   signed-near-zero-to-zero-action inference and the unsupported inference
   that a lower joint receipt weight alone gives a gap,
   not the actual
   complete shipping gap. E for corrected path/patch/quoting/whitespace
   tooling and the zero-gap certificate qualification; prior classified
   failures remain in DEAD_ENDS. No check weakened.
4. Genuine admitted shipping counterexample: none; no lossless impossibility
   or actual shipping negative eigen-direction has been established.
5. Retained: all 25 proved entries, CR25--CR27's actual-process/nonzero-
   receipt joint inequality and scalar zero-gap branches, locked-reference zero port/amplitude,
   suffix/generated-score, complementary Schur, coupled process/source,
   conditional AW/reader, frame, qualified magnetic, held and even-process
   results under their original hypotheses. No theorem/count promotion.
6. Next: decide CR41's feasible receipt corner using the linked residual
   process and transported acc/S/magnetic action, and its full Schur cost
   CR43, equivalently the ONE full R_W(c) comparison in CR42.
   The conditional algebra is not a sign or uniform margin; all inherited
   image and actual guard qualifications remain. Endogenous auxiliary
   columns cannot be relabeled external supply. General dissipativity and
   special planar all-time admission remain separate.

## Fresh-force INS runtime repair scope

**Current hypothesis.** Gyro-only frames changed the full INS frontend/heave
sampling history, and the compass force-norm latch could indefinitely suppress
INS aiding during legitimate waves. The repair retains every raw gyro increment
but clocks the full INS pipeline from fresh force; INS norm-latching is removed.

**Structures preserved.** Core estimator, covariance, bias, gain, OU/S and tuner
formulas, calibration, magnetic health gates and settled capture conditions
remain literal shipping operations. All five device pipelines retain their
states; raw gyro observations feed the guard and magnetic history. The private
frontend receives raw gyro, independent of learned main-filter bias.

**Relaxations introduced.** These are finite engineering replays with synthetic
bounded motion and sensor noise, not all-time MARINE/IMU/MAGNETIC qualification.
The force clock changes runtime chronology. No prior 5 ms frontend/release
certificate is promoted to the new 10 ms force clock; within-window bias/coning
modeling error and full device dropout paths remain unqualified.

**Failed approaches / DEAD_ENDS (E).** A separate NLO force-compensation equation
trial reduced a local identity residual but failed the existing recovery limit
(5.94893 degrees against 5 degrees). It was withdrawn; the equations are
unchanged. A reconstructed test initially used the wrong acceleration sign
for its displacement truth; its metrics were invalid until that sensor-generator
sign was corrected. The workspace reset interrupted validation and removed
uncommitted changes; the patch was restored and validation restarted. A parallel
OU-II test compile was killed by `cc1plus` memory exhaustion; its serial build
is the next check. The proof evidence command also stopped on missing `mpmath`
after the environment reset; dependencies are being restored.

**Retained facts / current limiter.** No genuine shipping counterexample or
end-to-end stability theorem is asserted. Established operation identities
remain controlling where their literal source binding applies; old device
source receipts require fresh binding, and prior statistical bundles are not
relabeled as new evidence. Hardware validation and CI remain necessary.

**Next falsifiable experiment.** Compare split and combined force-clock heave
on all five production facades, preserve existing heading/recovery limits,
then replay recorded device data and check motion/outage resumption.
