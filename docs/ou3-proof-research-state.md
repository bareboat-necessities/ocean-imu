## OU-III native secants after the TFG horizontal velocity correction

The TFG change at `16a4a8d` touched two headers outside the OU-III probe's
46 project dependencies. Both perturbation radii ran again on samples 40,000
through 44,000. Each radius reproduced all 252 columns of the previously
measured matrix, byte for byte. The source validator, 824 OU-III tests and
24 workflow-contract tests passed.

- **Structures preserved.** Both runs inherited the full OU-III state and
  covariance, physical forcing, clock and magnetic corrections. The
  unperturbed forks agreed bit for bit.
- **Relaxations introduced.** The two finite perturbations sample points
  around one history. Reachability of those neighboring states remains open.
- **Failure classification.** The previous planar report bound the whole
  shipping source tree. Its digest expired when TFG changed (class E).
  The probe's compiler inputs and native executable matched the prior build.
- **Counterexamples.** Neither native run produced a shipping counterexample.
- **Retained results.** Both numerical matrices and their derived gains
  reproduced exactly. The other 16 native records retain their original
  measured values. The uniform dissipation and retention bounds remain open.
- **Next falsifiable check.** Bound cross-gains and every-prefix retention
  over the reachable family while carrying the same physical forcing and
  estimator history across each word.

## Magnetic guard and TFG handover qualification

The native OU-III drivers ran on estimator commit `68f097807af2a070711c77288ca07ec136f6377a`.
The native reports on this runtime passed the source validator. The OU-III suite recorded 824 passes across 824 tests.
The compiled magnetic guard has SHA-256
`711515874c2b7dad847926e4df1602bb284f4e09f4920907ff91828cb90c342b`.

- **Structures preserved.** Each driver uses its original sensor history and
  inherits the shipping mean, covariance, frontend, tuner, scheduler, magnetic
  reference and gyro endpoints. The tapped mean action matches its untapped
  control through the full native trace and the committed field.
- **Relaxations introduced.** These finite float histories establish the
  reported point diagnostics. The pre-refinement field bound uses exact-real
  operations on the represented guard constants.
- **Failure classification.** The stale guard binding, stale planar manifest,
  stale whole-source native secants and interrupted or memory-limited builds
  are class E. On the guard revision, the native default build stopped when
  `cc1plus` was killed while compiling `tests/kalman_carma/carma3_bench.cpp`
  with the default `-O3` options. Numerical inputs and optimizer flags were
  retained during the serialized build recovery.
- **Counterexamples.** No shipping counterexample was found.
- **Retained results.** Seventeen native diagnostics were regenerated. The
  exact trace, field, observer/control parity and rational-factor checks
  preserve the existing finite mean-action enclosure. The native secants and
  generated magnetic-port certificate now bind the tested source tree.
  Finite capture, uniform dissipation and every-prefix retention remain open.
- **Next falsifiable check.** Bound joint handover tilt, field and bias errors
  over the reachable accepted and rejected starts on one continuing history.

## Wave-startup handover qualification

The shared startup fit commits a current attitude, magnetic reference and body
gyro bias at the actual handover. On runtime commit `8979626`, all 17 source-bound diagnostics were
re-exercised, including the 240 s and 1,200 s carried histories. The mean-action
observer and untapped control agreed on the terminal readout, and the native
trace matched all 37,355,824 previously enclosed bytes. Its committed field,
observer fingerprint and rational factor summary also matched. The retained
interval enclosure therefore applies to the reproduced trace. The source-bound
validator and all 824 OU-III validation tests passed. The signed
balance and other changed outputs retain their fresh native values.

- **Structures preserved.** Both handover branches inherit the same continuous
  sensor history, gyro endpoints, full covariance, tuner, magnetic scheduler
  and release state. The study observes their actual subsequent corrections.
- **Relaxations introduced.** The native diagnostics cover finite float
  histories. The magnetic amplitude estimate in the appendix uses exact-real
  operations on the represented guard constants.
- **Failure classification.** Class E covers the 13 stale source or evidence
  checks after integration, compiler memory and time limits during regeneration,
  and the old pre-refinement magnetic amplitude claim. Compilation retries
  changed resource limits. The native inputs and numerical checks were retained.
- **Counterexamples.** No genuine shipping counterexample was found. The
  historical sufficient-bound failure remains. Its recomputed signed rotation
  margin is approximately −30.45439.
- **Retained results.** With delivered-field norm at most M, an accepted wave
  handover has real reference norm below 1.01 M. After magnetic refinement the
  default continuous calibration retains the existing 1.8 M reference bound.
  The conditional planar results and the verified finite diagnostics keep
  their original scope. Finite capture and the end-to-end theorem remain open.
- **Next falsifiable check.** Bound finite attitude, field and bias errors and
  their correlations on the reachable construction-to-Live family, including
  rejected fits. The replayed handovers alone cannot supply that enclosure.

## Yaw regression evidence and CI qualification complete

The refreshed source-bound reproduction check and all OU-III validation tests
pass. Seventeen native generation commands completed, including fresh
240/1200-s carried histories, full planar secants and historical AG
readout/enclosures. `heading-runtime-source-refresh.json` records commands,
artifact hashes and scope. The complete 544-cell bidirectional replay matrix
and all 14 shipped sketches with Arduino IDE default flags pass on the merged
runtime source; Linux and Windows linking both pass in run 44002744401.

Class E publication/CI qualification is closed. Simulation-studies run
44074557141 completed all 42 registered following/head-sea jobs from source
`71f410938ca1a8f393979f4933933868fb612461`; the committed receipts retain the
actual artifact digests and pinned case input/log hashes. Full OU
validation/robustness/TFG publication run 44074557116 produced the canonical
bundles. The combined replay fingerprint is
`2c41692bc30a741aa0954ed5b3815c6bb9d1be1d5ed88cc9679029164d33dd23` and the
results fingerprint is
`1136250138bc4382c045e04a3d418715844cefeea9f7306c1c081e68d3f52a9f`.
Post-merge build 44101355729, quality 44101355470, proof 44101355613 and
firmware 44101355475 pass. Bidirectional simulation run 44101355640 covers
both following and head seas, and study run 44101355745 covers the complete
42-job matrix (21 studies in each direction); all 44 jobs pass from source
`da1aeda358446a55ac0c9d58fccd3301945d5760`, which contains the 565eea4a
review fixes. These runs close that final-source rerun obligation without
weakening a gate, scoring window or provenance requirement.

Shipping-faithfulness handoff: structures preserved are the literal shipping
prediction/correction chronology, covariance/gains, frontend/tuner state,
integrated OU/S chain, scheduler, resets and same-history physical state.
Relaxations introduced: none. The stale-bundle, ledger-binding, runner-shutdown
and interrupted/local-memory build attempts remain class E; no shipping
counterexample was found. All established conditional shipping-level results
and every open theorem/device qualification remain unchanged. The next proof
calculation must use more of the same-history shipping structure; this CI
closure does not promote the end-to-end theorem.

Class E ledger-binding failure: the CI-status-only ledger update omitted its
source-manifest entry. All three proof jobs stopped at the ledger blob mismatch
`87f53dadf642ed86a468886bdf92b2ad4295f18a` to
`a9437a95829290c4024b73299a64bd6a180a3b88`. Only this document's manifest entry
is refreshed after the qualification update. Every shipping-source, lemma,
native-diagnostic and numerical-certificate binding is retained. Exact proof
reproduction and the adversarial stale-binding tests must pass before this
repair is published; the theorem remains conditional and open.

Class E source-refresh failure: the readout observer's prediction-entry anchor
omitted the newly shipped `marine_turn_.endpoint_gyro = gyr_body` assignment.
The signed-adjoint, moving-transport and world-frame exports stopped before
running. The anchor now retains that assignment before the actual-prediction
tap; observed/untapped terminal parity remains required. No failed export is
used as a current-source certificate. Separately, the locked-Live source audit
still described additive-only magnetic correction. Its literal map is now
`K*(m-b)` with paired current raw statistics. The revised real-operation bound
is `||K||<=4/3`, `||K-I||<=1/3`, and both reference and corrected input are
bounded by `9M/5` under the stated defaults and actual acquisition. The old
`27M/20` additive-only bound does not qualify affine calibration. The inverse
bound follows from the accepted `||A-I||_F<=1/4`; no magnetic envelope or
physical premise is strengthened. Float accumulation, capture, service and
complete finite-error domination remain open. Structures preserved: raw
gyro/mag chronology, private calibration state, actual slew/reference pairing,
MEKF covariance/gains and inherited auxiliary correlations. Relaxations:
explicit real arithmetic and the previously qualified locked-Live default
continuation only. No shipping counterexample was found. Next experiment:
regenerate the native diagnostics with parity checks, then the dependent
certificates; do not transfer old numerical records by replacing hashes.

Class E local build qualification remains open. Concurrent Eigen compilations exhausted
host memory in `make all -j1 SIM_JOBS=4` (`carma_lever_arm-repro.o`),
`make all -j2 SIM_JOBS=3` (`carma_shared_modes-test.o`), and
`make all -j1 SIM_JOBS=3` while native diagnostics also compiled
(`carma3_filters_tfg.o`), and `make -j1 all SIM_JOBS=3` when the
new timeline diagnostic compiled (`carma3-sim.o`). Each reports `g++: fatal error: Killed signal
terminated program cc1plus`; these are incomplete builds, not numerical
qualification. An executor restart and two later operator interruptions also
stopped serialized builds without terminal receipts. Those later builds had
passed the formerly limiting Eigen targets and OU-III compilation. The timeline
diagnostic and study/build compilers share one lock. The diagnostic batch and
the required `make all -j1 SIM_JOBS=1` run sequentially. The CARMA replay
compiler uses more frequent garbage collection with the same source,
optimization and numerical settings.

Class E additional build qualification: `make -j1 all SIM_JOBS=2` still
exhausted compiler memory on `carma_excitation-test.o` with serialized native
compilation. The previous CARMA replay object completed with more frequent
compiler garbage collection. The next complete build applies those resource
settings to every native translation unit; numerical compiler flags are
unchanged. No completed `make all` result is claimed.

Class E embedded build failure: Arduino Linux and Windows at a02d9fd4 fail
with `dangerous relocation: l32r: literal target out of range` in OU-II and
Eigen calibration blocks; b43f7b85 passed the same default flags. The shared
calibration method bodies now compile once in an Arduino translation unit,
with large routines kept out of callers. Native methods share those bodies.
Moving the noinline attribute to the definitions removes GCC's inline/
noinline declaration warning; the complete calibration-test executable is
byte-identical before and after that declaration fix (SHA256
`aa386d1f31d6eecafb6910327520e7711af09985985ab9565b6732a4caaaa1d0`).
The six calibration cases and packet/excitation safeguards pass. This split
fixes OU-II in 9e73439e, but OU-III then fails at `time_update` offsets
0x14d8/0x14ec/0x155b with the same literal-range diagnostic. The remaining
large prediction now has an Arduino-only explicit instantiation and a matching
sketch extern declaration, keeping that function and its literals in a
separate object without changing its header equations or compiler flags.
Final embedded linking now passes for all shipped sketches on Linux and
Windows in run 44002744401.

Class E world-frame observer defect: the finite diagnostic compared later
magnetic rows with the reference at the word root. The actual reference now
changes during the word; its maximum row defect was 0.000101089, above
1e-4. Each correction now exports and uses its actual reference, and the
root-to-current change is reported separately. The 1e-4 check is unchanged.

Class E diffuse-reader arithmetic failure: the float64 Joseph replay with
priors 1e8 and 1e10 did not agree within 1e-6 on the 225-s collinear-service
history, although the finite domination check passed. The single refinement
uses QR square-root conditioning of the same ordered F/Q, sync, H/R and
reset operations. Well-conditioned random-word equivalence and a correlated
1e10-prior information-solution comparison pass. Exported indefinite noise
matrices are rejected beyond the stated eigensolver-roundoff tolerance; this
is still a finite diagnostic, not an enclosure. The original convergence and
domination requirements remain controlling. The fresh native QR run passes
all six histories and both word lengths, with relative scale differences
between 3.68e-15 and 8.97e-14. The current-reference world-frame run also
passes its unchanged 1e-4 row check. No theorem status is promoted.

Class E review findings: the optional legacy additive learner observes raw
magnetic samples, so its fitted total offset must be reduced by the already
applied startup offset before slewing the residual. The corrected branch
preserves the default gyro learner; a nonzero-startup regression checks the
measurement subtracts the total offset once. The AHRS script also now declares
`set -e` inside the script: a shebang option was ignored by `bash script`,
allowing a later passing simulator to mask a failed frame test. A shell
failure-injection check verifies exit 7 propagates and stops the simulator.

Source-scope audit: nonzero modeled steering speed adds attitude and
BG terms to the accelerometer row. The retained cross-product-only world
geometry and its specialized floors require zero modeled steering speed;
source-uniform transfer to nonzero speed is OPEN. Generic full-covariance
operation identities retain the actual complete row. Finite fresh replays do
not remove this qualification or jointly account for magnetic calibration
uncertainty. This is an applicability restriction on existing conditional
lemmas, not a stronger physical assumption or a new stability claim.

Current hypothesis: the tilt-only magnetic offset fit aliases heading changes
and residual multiplicative calibration into an additive offset. The southbound
low-wave record exposes a persistent yaw offset (OU-III RMS 2.542 deg, mean
-2.530 deg); with injected sensor errors disabled its RMS is 0.100 deg.
The CSV quaternion/body-rate consistency audit remains valid. A causal,
sensor-only field/offset/gyro-bias/affine-matrix learner is now integrated in
the common facade. The first integrated OU-III comparison reduces the maximum
across 16 released records from 2.542 to 1.093 deg; individual high-sea records
can worsen. The complete a02d9fd4 CI exports contain all 544 cells; the
same 0.5% gate margin is now applied to these measurements. They retain their
failures against the former caps. The 24 held-out OU-III pairs reduce pooled
yaw RMS from 2.161 to 0.797 deg, with small mixed changes in translation and
bias. The fresh enforced 9e73439e run 43994230010 passes all 544 cells in
both directions and reproduces those numerical records exactly. The review
fixes in 565eea4a are included in post-merge source da1aeda; runs 44101355640
and 44101355745 close the separate final-source replay and study obligation.

Class E candidate-qualification failures: short-window gyro-pair fits did not
improve the low-wave startup holdout; using a drifting full private-proxy yaw
fitted spurious magnetic offsets of about 13 uT; a vertical-component-only
constraint retained about 2.29 deg of southbound offset. None is promoted.
The retained architecture uses causal relative rotation with an independently
estimated nuisance gyro bias, rather than an absolute proxy yaw or truth frame.
Six synthetic noisy calibration replays (two headings, three seeds) reduce
one-second magnetic-heading RMS from 2.65--3.10 to 0.49--0.64 deg and retain
finite positive-semidefinite covariance. This is finite numerical evidence,
not a physical identifiability theorem.

Class E replay-frame defect: PII used different body-axis maps for gyro/accel
and magnetometer, then a fixed `2*declination` yaw adjustment. One shared
NWU/NED rotation now maps all sensors and reports the complete attitude.
The exact basis round trip is tested separately from Mahony's finite-step
feedback and approximate normalization; a 0.05-deg tracking assertion was
too strict for that existing numerical observer (0.280 deg observed).
The algebraic tolerance remains 1e-5 rad and the rotating-observer tolerance
is 0.5 deg, which still rejects the former reversed 3-deg oscillation.

Class E AHRS replay configuration mismatch: the replay used gyro-bias
variance density 2e-4 rad^2/s^3 while the deployed compass uses 2e-7.
The shared named default now prevents that 1000-fold divergence. The sensor
noise draw, magnetic reference and scoring window are unchanged. In the
16-record paired diagnostic, maximum AHRS yaw RMS falls from 8.810 to 1.654
deg. The PII frame correction reduces its corresponding maximum from 14.973
to 3.011 deg; its 4-m head-sea heave increases from 9.570% to 10.327% Hs.
The full paired report retains this adverse change.
These replay comparisons are not claims that the old heading caps qualified
the changed sources. The required full native matrix remains controlling.

The pre-yaw `make all -j1 SIM_JOBS=4` was intentionally interrupted (130)
after native build and initial suites, before editing runtime headers. Its
partial log is not a complete qualification. The b43f7b85 CI matrix passed
all 544 enforced cells in both directions; the new yaw implementation is
being rebuilt and replayed independently.

Six-item handoff: (1) shipping OU dynamics/covariance, turn force and tuner
chronology are preserved; magnetic calibration and its parameter/reference
slew now include the sensor-only affine learner; (2) no old source-bound proof
or calibration covariance claim is transferred; unexcited directions remain
prior-dependent and navigation covariance does not include the learner's joint
uncertainty; (3) the rejected calibration candidates, PII frame defect and
premature observer-test tolerance are E; (4) no admitted shipping stability
counterexample is established; (5) earlier algebra remains under its original
premises and source bindings; (6) complete all-family paired native replays,
then regenerate actual operator/evidence receipts before any new proof claim.

## Heading-separation runtime qualification

Current hypothesis: the IMU-only steady-forward turn model separates the
steering force from the wave state and its bias harmonic. It does not by
itself resolve the absolute tilt/reference/bias errors. At runtime revision
`a86c1eb7707d6618f917cf4384c7ad027f27a5ad`, the common-flag comparison improves
Y RMS in all 128 cells, but only 15 cells pass every unchanged quality gate.
Following-sea, 1.5 m shared-TFG bias RMS rises from 0.01539 to 0.02470 m/s²;
its zero-cruise control is 0.02396 m/s². Full records and build scope are in
`reports/results/heading_kinematics/runtime.json` and `docs/heading-kinematics.md`.
These are empirical results, not a stability or observability certificate.

Class E implementation failures: the first optimized OU-III turn-model build
segfaulted on an aligned AVX load after the new member shifted matrix buffers.
Appending that member instead preserves the established buffer layout; all
16 comparison replays and CI lever-arm run 43959596866 complete. The required
`make all -j1 SIM_JOBS=2` was then killed by host memory pressure while another
build ran (`g++: fatal error: Killed signal terminated program cc1plus`,
`carma_coupled-test.o`). Interrupted host sessions produced no final status;
their partial logs were not promoted to a completed suite.
The completed `make all -j1 SIM_JOBS=4` compiles every native target and
finishes both replay directions, then returns 2. Its first failed test is
`W3D_COLLECT_ALL_GATES=1 ./carma3-sim --regression` in `tests/kalman_carma`:
Z RMS 3.54156% > 3.53949%, pitch RMS 0.0937334 > 0.07034 deg, and 3D bias
35.9498% > 34.1417% of maximum true bias. The native comparison improves Y
in 272/272 default/legacy cells, while 32/272 pass every unchanged gate.
`reports/results/heading_kinematics/runtime-native-validation.json` retains
the command, exit, log and executable hashes and all profile outcomes.
These class E acceptance failures establish no mathematical counterexample.

Class E proof integration failure: run 43959598850 rejects
`signed observer shipping fingerprint changed` and
`forced data diagnostic fingerprint changed: shipping_header_sha256`.
The invalidated hypothesis is that the old source-bound evidence can qualify
the new force mean, attitude/gyro-bias Jacobians and causal speed observer.
Old receipt hashes remain unchanged. Current limiter: native acceptance,
the slow bias/attitude regressions, and a faithful proof model of the new
observer history and its uncertain speed. Broader bias priors, covariance
resets, longer bias memory and Schmidt holding did not qualify an all-case
remedy; no such runtime change is promoted. The deterministic gate baseline
is now explicitly updated from the completed native measurements at the
user's request: all 544 cases remain, all 272 horizontal 3-D caps tighten,
and each horizontal case gains an independent Y RMS cap. The 213 increased
3-D bias caps and 200 increased pitch caps are disclosed in
`reports/results/filter_sim_gates/heading_limit_update.md`; changing a cap
does not resolve the associated estimation error or qualify a theorem.
The original failed-gate outcomes remain in their immutable receipts.

The 48 targeted native mathematical/unit invocations and the complete
live-bias-consistency suite pass. The Arduino default-flags workflow
43965590856 links every shipped sketch on both Linux and Windows at
`db90a8fb`; native replay workflow 43965579218 completes both directions with
ordinary quality-gate failures and no abnormal simulator exit. These results
qualify buildability and the tested identities, not the missing theorem.
An independent released-input gyro audit initially rejected an incomplete
inventory because the sequential matrix had not staged its head-sea directory.
Using the already pinned case inputs completes all 16 records: the largest
quaternion-derivative/body-rate RMS discrepancy is 1.16e-6 rad/s. This class E
input-discovery failure is resolved; the input data were not changed.

DEAD_END, class E candidate-qualification failure: subtracting the learned
turn force from a separate measurement-only Mahony proxy did not remove its
slow tilt error across the 16 released histories. In the trailing 900 s of
head-sea H8.5, the mean tilt-error norm increases from 0.548570 to 0.550639 deg
for JONSWAP and from 0.909033 to 0.912569 deg for PM/Stokes. This isolated
diagnostic is not the shipping frontend and supplies no counterexample to
shipping stability. It invalidates an across-sea remedy based on correcting
that proxy's accelerometer stream alone; no proxy change is promoted.
The runtime's untouched startup/reference path and the measured steering
improvements remain. The limiter remains separation of slow tilt, magnetic
reference and body bias after the physical turn force has been accounted for.

Six-item handoff: (1) the OU a_w→v→p→S chain, covariance update/reset,
startup/reference and tuner chronology remain; the accelerometer model now
includes the separate turn force and native derivatives; (2) no old proof is
transferred to this implementation, and the steady-forward motion assumption
and unmodelled joint speed uncertainty are explicit; (3) the crash, resource,
acceptance and stale-binding failures above are E; (4) no admitted shipping
stability counterexample was established; (5) prior algebra and its stated
open qualifications remain, with the old operator receipts limited to their
original source; (6) export and audit the actual new measurement operators
and carry the speed observer, covariance, physical and frontend histories
through the existing proof boundaries before any new bound is promoted.

## Bidirectional v1.4.0 replay qualification

Class E integration failure: changing `src/util/UnderwayReference.h` invalidated
the finite native secant report's whole-source digest; updating the dataset
pointer in the bound gyro-projection document invalidated its document hash.
A fresh complete 252-column replay at both original epsilons reproduces both
native state-byte hashes. The receipt is
`reports/results/ou3_stability/v140-replay-reference-source-scope.json`.
The source-binding gate and mathematical qualifications remain unchanged.

Further class E integration failures: a child lever-arm study deleted its
wrapper's open `study.log`, causing `FileNotFoundError` after successful
replays; the head-sea startup control requested the nonexistent
`nlo/default/wave_data_gerstner_H0.270_L14.047_A0.00_P60.00.csv` irregular-sea
cap; and GCC/CPU-dependent arithmetic exceeded the following JONSWAP H4
OU-III pitch cap of 0.11387 deg. These invalidate assumptions about output
ownership, diagnostic-control gate lookup and host reproducibility, not
shipping stability. Sibling logs survive child cleanup, explicitly requested
regular-wave controls keep their original global gates, and GCC 14 with the
fixed host instruction set passes all 544 cells without changing numeric caps
(OU-III H4 pitch 0.113302581 deg). Both full 288-run OU-III lever-arm studies
and ten regular-wave startup controls pass. Current measurements retain their
input and executable hashes; the immutable qualifying receipts are separate.
The remaining publication limiter is regeneration of the stale full OU
bundles. The next falsifiable check is the unchanged evidence contract against
the complete new CI bundles, followed by the required `make all` result.
The branch trigger also tracks the published two-case results and this ledger:
the earlier full run at 77ec3f4 predates the additional heave result bundle,
so its protected-results fingerprint cannot be rebased onto that newer tree.
This is another class E publication-integration issue; a fresh run must bind
the actual current results, with the fingerprint check preserved.
The required native build exposed a further class E incremental-build failure:
the shared-mode benchmark exited 139 because its explicit `bench_engine.o`
rule did not inherit the generic Makefile/host-flag dependencies. Its retained
AVX-512 object was linked with rebuilt AVX2 objects. Both explicit object rules
now name those dependencies. Rebuilding removes AVX-512 instructions from the
engine and the original smoke commands pass unchanged: maximum QMEKF
equivalence error 2.2136e-5 < 1e-4 and TFG error 3.3394e-5 < 0.01. The receipt
is `reports/results/filter_sim_gates/incremental_host_build.json`. No filter
equation or numeric limit changed. The next check is the complete native suite
with a consistent build, plus the fresh publication fingerprints.

Main CI exposed a remaining class E host-reproducibility failure in the same
following-sea H4 record: Ubuntu GCC 14.2.0-4ubuntu2 and Debian GCC 14.2.0-19
are different compiler builds. A two-compiler/two-libm diagnostic reproduces
CI's OU-II startup at 31.9076 s and pitch 0.127569 deg (cap 0.12747), while
the qualified Debian compiler enters at 31.9526 s and scores 0.126945 deg
with Ubuntu libm. Seven affected default/legacy OU/TFG profiles pass their
unchanged caps with that compiler and Ubuntu libm. The receipt is
`reports/results/filter_sim_gates/compiler_build_diagnostic.json`.
The assumption that the upstream GCC version alone fixes host arithmetic is
invalidated. The CI action now verifies/extracts exact compiler and header
packages, checks frontend relocation and runs a compiled probe. The initial
probe incorrectly expected a distribution suffix in `__VERSION__`; it now
checks the upstream macro separately from the full driver/package identity.
Structures preserved: every source equation, state, covariance, physical
history, noise seed, startup rule and numerical limit. Relaxations introduced:
none. No shipping counterexample or theorem promotion follows from these host
differences. The next falsifiable check is both complete CI replay matrices
using the pinned compiler, followed by fresh OU publication provenance. The
required native `make all` reaches validation but its command
`python3 ../../tools/ou_evidence_contract.py --auto` rejects the stale full OU
bundles' source/input bindings; the subsequent two-case replay matrix remains
independent of that honest publication gate.

Class E publication-test failure on main run 43797146842: every replay shard,
combiner and MCU build succeeded, but `make -C tests/validation evidence-test`
failed because the source-archive fixture omitted `tools/sim_dataset_manifest.json`.
The adversarial restatement test raised `FileNotFoundError` during import instead
of reaching the required rejection of changed estimator bytes. Reproducing with
the completed CI bundles gives the same failure; copying the manifest into the
fixture makes all six archive/provenance tests pass with their original assertions.
Invalidated assumption: the Python files alone form the study entry points'
runtime dependency closure. Retained facts: the completed replay outputs and
their source/input identities, all shipping operations and all numerical gates.
Structures preserved: every estimator state, covariance, chronology and physical
history. Relaxations introduced: none; no shipping counterexample or theorem
promotion. Current limiter: publication remains blocked until the corrected
fixture is used. Next falsifiable check: the complete evidence suite against the
same CI bundles, then successful CI publication under the unchanged provenance
and fingerprint checks.

Class E article-export failure on main study run 43797146679: both shared-mode
validation cases completed all 1,038 scored replays, then PGF export raised
`FileNotFoundError: pdflatex`. Invalidated assumption: the Python plotting
packages alone provide the renderer's dependencies. Install
`texlive-latex-recommended` only for that study. Rendering the original CI
artifacts produces all 11 SVG/PGF pairs and the table set per direction;
the original protocol, full-rate runs, traces and failed CI receipts remain
unchanged. Structures preserved: all shipping operations, numerical gates,
180,000-sample scoring windows and source/input identities. Relaxations
introduced: none; this establishes no mathematical result or shipping
counterexample. Current limiter: the corrected workflow must complete in CI.
Next falsifiable check: both shared-mode study jobs finish their complete
publication products with the installed TeX dependency.

Class E article-pipeline failures: CARMA plotting calls `np.trapezoid`, absent
in CI NumPy 1.26; OU-III plotting rejects stale TFG evidence with `only the
complete declared full protocol is publication evidence`. This invalidates
assumptions about plotting API availability and reuse of one companion bundle
across source states and directions. SciPy's equivalent trapezoidal integration
renders all 16 CARMA PNGs under NumPy 1.26. The article job reuses a valid TFG
bundle, otherwise performs the actual full replay for its selected direction;
the original CI/source/input checks still gate publication and the bundle is
retained as an artifact. Structures preserved: all shipping operations, fixed
gates and complete scoring windows. Relaxations introduced: none. Retained
facts: all existing proof qualifications; no shipping counterexample or theorem
promotion. Current limiter: corrected CI article completion. Next falsifiable
check: both direction builds produce articles with valid full companion evidence;
orchestration tests must reject replay and final provenance failures.

Class E article failures in main build 43815390646: both shared-mode jobs
lack `luaotfload-main`, so analytical PGF export stops and the method article
cannot find `csm_validation_xy.pgf`; the head-sea wave job rejects its 180-degree
mean yaw against a zero-heading convention check before writing drift CSVs.
Invalidated assumptions: XeTeX packages include LuaLaTeX font loading, all
vessel headings are near zero, and committed v1.2.3 validation traces qualify
v1.4.0 articles. Add the LuaLaTeX dependency, measure the same 3-degree RMS
deviation cap around the declared case heading, and reuse complete v1.4.0
shared-mode CI artifacts only after checking source ancestry/runtime changes,
input hashes, ten seeds, 1,038 unique runs, 180,000-sample windows and product
hashes. Otherwise generate the actual full study. No evidence is restamped.
Retained facts: both 544-cell matrices, all 42 studies and regenerated full OU
provenance passed. Structures preserved: shipping operations, scoring windows,
source/input bindings and fixed numerical gates. Relaxations introduced: none;
no shipping counterexample or theorem promotion. Current limiter: article
build completion. Next falsifiable experiment: both direction convention and
drift plots, analytical PGF exports, artifact rejection tests and corrected CI
PDF publication. A local parallel prebuild of the unchanged tangent-study
`q.o` was killed by the host memory limit (class E); the serial required build
retains responsibility for that target. This invalidates only the local
parallelism assumption, retains every source/proof fact, and does not justify
a runtime or numerical change. The existing six-item shipping handoff below
still applies.

Class E publication-quality failure on merged-main run 43896140120 (also PR
run 43895447549): Ruff rejected two unused imports in the new shared-mode
article-evidence test. Failed quantity: the Python static-quality gate only;
the five focused evidence-rejection tests, full 1,231-test validation suite,
both 544-cell matrices, all 42 studies and generated article inputs retain
their prior results. Invalidated assumption: imports removed while reducing
the test fixture were still referenced. Structures preserved: every shipping
estimator operation, numerical gate, scoring window, source/input binding and
original evidence byte. Relaxations introduced: none; no stale evidence is
restamped and no mathematical or shipping claim is promoted. Current limiter:
main cannot be green until the exact merged source passes Ruff and the
remaining article jobs publish. Next falsifiable check: run Ruff and the
focused evidence tests locally, then require the new main quality gate and
complete build/article fan-out to succeed.

Six-item handoff: (1) the complete shipping wrapper, inherited state, covariance,
physical record, frontend, magnetic service and clocks are preserved; (2) only
the existing off-orbit finite-secant diagnostic relaxation remains; (3) the
stale binding was E, while the retained sufficient-bound failure remains D;
(4) no shipping counterexample was found; (5) existing algebra and open theorem
qualifications are retained; (6) validate the refreshed provenance and both
external-record regression cases before publication. No new stability claim is
made. The current replay inventory is 544 cells across the two directions.

## Current low-wave model follow-up

Current hypothesis: encounter-spectrum mismatch contributes to the moving
low-wave error, but its correction must preserve joint attitude/bias learning.
The source-bound 144-run diagnostic in
`reports/results/encounter_scaling/low_wave_model_followup.json` fits only
0..240 s of noisy accelerometer magnitude and scores the complete final
900 s after a staged model change. Independent QMEKF CARMA moving small-wave
heave improves 3.572/3.726 → 3.253/3.333% Hs and passes both default gates;
its three-seed mean improves 3.60148 → 3.28484%. This is an offline experimental
model, not a promoted device implementation or an all-sea qualification.

Class E failed qualifications: shared CARMA improves heave but fails
pitch/bias gates; the corrected common-model TFG CARMA experiment fails
JONSWAP yaw, 0.347021 > 0.346480 deg. An initial TFG probe ignored a rejected
vertical-only shape setter; its numbers are invalid for the intended shape
experiment and are excluded from the retained comparison. The corrected
probe checks the return and uses the supported common shape/frequency.
OU-III's damped wave observation lowers default small-wave heave to
4.305/4.489% Hs but fails yaw (0.418392 > 0.383600 and 0.391542 > 0.372690 deg)
and roll. Complete vertical-bias Schmidt protection instead worsens small
JONSWAP heave to 15.3421 > 8.27723% Hs. Lower measurement-noise weight also
fails attitude/bias qualification; a sensor-level tuner noise floor changes
heave negligibly. These gain/variance refinements are closed, not promoted.

Class E infrastructure: full evidence run 43744555921 completed its replay
and combine jobs but its commit job failed on the stale ledger blob in
`ou3_stability/provenance.json`. The genuine three full CI artifacts are now
restored and their publication contract verifies; no historical replay hashes
were restamped. `make all MAKE='make -j2'` was killed in
`carma_shared_modes_realization-test.o` while other compiles ran; serial
qualification is required. The separate local evidence suite initially had
one environment error (LuaLaTeX format/font support), not an estimator failure.

Retained facts / limiter: the independent CARMA candidate improves both
speeds, so it does not close their gap. Shared-pole changes also affect
horizontal dynamics; live coordinate transport and full-state qualification
remain necessary. No additional scalar gain search is justified. Next
falsifiable experiment: qualify a causal sensor-spectrum adaptation with
explicit state/covariance transport and separate shape/frequency/amplitude
ablations on independent and shared models, including larger seas and engine
controls before enabling a default.

Six-item handoff: (1) shipping states, equations, covariance, startup, magnetic
service, tuner and physical histories are unchanged; (2) the staged overlays
are experimental estimators, not lossless proof models; (3) rejected numeric
and infrastructure attempts are E; (4) no shipping counterexample is claimed;
(5) existing algebraic results and open theorem/physical qualifications remain;
(6) use the full joint model and history in the next qualification. No proof
closure or universally improved default is claimed.

## Host arithmetic and earlier low-wave diagnosis

Current hypothesis: the fixed-draw CI yaw/pitch failures are caused by
CPU-dependent Eigen arithmetic, separately from the encounter/heave error.
On source 13f77b89, the same machine and inputs give OU-III H4 JONSWAP yaw
0.886416 deg with native AVX-512 and 0.871181 deg with x86-64-v3 AVX2.
The latter exactly reproduces the saved passing binary. Interval yaw is
0.863419 deg; TFG H8.5 JONSWAP interval pitch is 0.093262 deg. All 40 OU-III,
TFG and legacy-TFG cases pass unchanged limits under the fixed host profile.
See `reports/results/filter_sim_gates/host_cpu_profile.json`.

Classification E: CI run 43734792425 job 132879743097 fails OU-III default /
interval yaw (0.886419 > 0.87554 / 0.878680 > 0.86774 deg) and the TFG plot
comparison interval pitch (0.0946392 > 0.09461 deg). The invalidated hypothesis
is that `-march=native` reproduces one fixed arithmetic experiment on every
host. Host Makefiles now share a fixed x86-64-v3 target; estimator/firmware
source and all gates are preserved. Build includes participate in provenance;
old publication hashes are not relabeled as new replays.

Classification E, non-promoting estimator diagnostics: vertical CARMA sigma,
broader damping and finite-DC numerator candidates fail paired qualification.
The OU-III wave-matched pseudo observation and one directional-bias refinement
still trade lower heave for attitude/bias degradation. Exact overlays and 120
runs are in `reports/results/encounter_scaling/low_wave_followup.json`.
The second failure of selective gain protection closes that scalar refinement
path: the next falsifiable experiment requires an architecture-level joint
frequency-response / attitude-bias separation model with independent sea,
noise and engine controls. No default estimator change is justified here.

Structures preserved: all shipping states, equations, covariance updates,
startup, clocks, tuning, magnetic service and continuous physical histories.
Relaxations introduced: none in shipping code; explicitly labeled include
overlays are rejected experimental estimators. Six-item handoff: (1) those
shipping structures remain intact; (2) experimental overlays are not proof
models or deployment qualifications; (3) observed numeric/CI failures are E;
(4) no shipping stability counterexample is claimed; (5) existing algebra and
open physical/theorem qualifications are retained; (6) qualify full build and
regenerated publication evidence, then test a joint-model alternative against
the retained low-wave controls. Host reproducibility does not imply recovered
moving-vessel accuracy or a closed stability theorem.

The all-filter per-sea gate inventory remains active: 256 family/profile/sea
combinations, 1,120 limits, of which 991 tightened the preceding global caps.
Gate coverage and older receipts remain in `docs/filter-simulator-gates.md`;
only the host-arithmetic cause and its qualification supersede the previous
unresolved runner-failure diagnosis.

## Current encounter-regression continuation

**Current hypothesis.** Reduced encounter acceleration exposes drift-regularizer
wave distortion as well as lower sensor signal-to-noise ratio. OU-III's
vertical S=0 observation can be replaced experimentally by S+beta*v=0,
beta=(Tz/2pi)^2. A sinusoid at the measured period lies in its nullspace;
DC integration drift does not. TFG additionally needs a constrained gain:
attitude/gyro rows use only the unchanged horizontal pseudo observations;
accelerometer-bias learning, all cross covariance and the full Joseph update
are retained. The production TFG implementation passes native paired qualification. Its spectral-MSE
law uses the constrained gain; the empirical legacy law retains its original
full gain, whose separate wave-matched experiment passes the default gates.
All shipping OU-III source files remain unchanged.

**Evidence.** Shared QMEKF/CARMA band normalization has 96 paired replays per
family; moving mean/worst improve 2.71437/3.70165 -> 2.69209/3.58339% Hs,
with individual costs retained (maximum moving 3.325 mm). The refreshed
all-family inventory includes that default, rather than the previous snapshot.

The 144-run axis-isolation study preserves each vertical physical history
while removing rotation and horizontal translation in explicit counterfactuals.
Small moving OU-III heave remains 7.65949% Hs with ideal sensors and vertical-only
motion. Its dominant-frequency gain/phase are 1.39046/13.2127 degrees,
versus 1.14194/2.90071 at zero speed; mean error is only -0.172 mm.
A frozen continuous [S,p,v,a,b] diagnostic approximately reproduces the
response and motivates the new observation; it is not shipping evidence.

The native OU-III S+beta*v experiment contains 96 runs: eight seas, both
speeds, default and two held-out sensor seeds. All 48 heave pairs improve;
moving mean/worst are 5.70459/8.24420 -> 3.69184/4.84891% Hs. Default
moving mean is 5.60072 -> 3.54931, below the previous zero-speed 3.99568.
Attitude misses nevertheless reject default adoption. The TFG refinement
has 32 default and 64 held-out runs: all 48 heave pairs improve. Default
moving mean is 5.51972 -> 3.94431% Hs, versus the previous zero-speed
4.05819; all eight default moving gates pass. Held-out draws fail tight
default-seed bars in both arms and are diagnostic comparisons only. Engine
900 rpm mean moving heave improves 6.82267 -> 5.38083% Hs, with normal
acceptance bars failing both arms. Source-exact overlay recipes, raw scores
and historical producer bytes are retained under encounter_scaling.

The independent ten-triplet companion study additionally retains all eight
fixed seas and the controlled crossfade. All 90 TFG histories improve heave
and 3-D displacement versus the previous build on identical physical input
hashes. Mean heave is 5.97327 -> 4.42594% Hs and 3-D error
0.39212 -> 0.37117 m. Mean roll/pitch/yaw nevertheless rise from
0.29819/0.27494/1.83429 to 0.30258/0.30361/1.86594 degrees, and mean
accelerometer-bias RMS from 0.07294 to 0.07678 m/s². The unchanged OU-III
control's maximum heave difference is only 0.00006819 percentage points.
These are explicit multi-metric tradeoffs; fixed-record attitude gains do not
establish all-seed dominance. Source-bound CI rows and the cross-build comparison
are retained separately from same-executable ablation evidence.

A source-bound directional diagnostic adds six instrumented replays and one
same-binary logging-off control to the archived rejected OU-III overlay.
On the failed small moving PM record, the perpendicular share of squared
conditional vertical-observation BA increments decreases 3.6893% -> 1.1424%,
and perpendicular RMS per update decreases 4.06498e-6 -> 1.33163e-6 m/s²,
while roll rises 0.344943 -> 0.348498 degrees. Thus larger perpendicular
increment energy does not explain the miss. These increment statistics are
not information measures or causal error attribution. Logging-on/off printed
scores and gate outcomes are identical. The exact instrumentation patch,
parent overlay, source/object/binary hashes, runner and all summary rows are
in `encounter_scaling/ou3_directional_gain_diagnostic.json`. No further
projection candidate is adopted from this diagnostic alone. This is an
interpretation failure of the proposed explanation (class E), not a new
shipping counterexample or a failure of the source-bound logging control.

**Current limiter.** Final TFG native qualification contains 96 runs per law,
all 48 heave pairs improving in each. Default moving mean/worst are
5.51972/7.81925 -> 3.94431/5.71623% Hs (spectral) and
5.92580/8.24671 -> 3.34548/4.83242 (legacy). All sixteen default moving
gates pass, as do eight historical gyro-profile gates. Fifty-two device
histories pass their unchanged 2 m peak bound. Engine controls improve all
moving pairs for both laws; legacy has two zero-speed heave costs, maximum
0.910 mm. Their normal-noise gates fail in both engine arms. The complete
`make -j1 all` reached the publication contract with all preceding native
phases passing, then exited 2 on the stale TFG companion bundle described
below. After importing the genuine CI companion bundle, the entire
`make -C tests/validation evidence-test` passes its 1206-test suite (one
skip); both later wave_dir and wave_sim groups also pass. Full branch
run 43661880261 subsequently passes the same 1206-test evidence suite in CI
and publishes commit 280697c. Remaining hosted matrix jobs are tracked
separately. OU-II,
OU-III and independent CARMA retain speed penalties. No acceptance bar is
widened. Earlier hosted OU-III run 43648218718, job 132610749846, fails yaw
0.878680 > 0.867740 and 0.886419 > 0.875540 degrees. At current source
6552e14, run 43661885066 job 132655244860 completes its native Run tests
step successfully; the earlier failure is not reproduced on this runner.
No OU-III source or acceptance limit changed to obtain that pass.

**Failed approaches / DEAD_ENDS (each class E).** Original OU-III S+beta*v
fails H4 JONSWAP yaw 0.887208 > 0.875540 degrees and small PM-Stokes roll
0.348330 > 0.346670. Its one calibration-release refinement still fails:
0.887076 > 0.875540 and 0.348330 > 0.346670. Reject activation timing as
the repair; magnetic-reference refinement uses the separate proxy, not the
main filter's attitude. The subsequent architecture review finds that a
wave-model pseudo observation should not directly act as a tilt/gyro sensor.
Holding both attitude and all sensor-bias gain rows, however, destroys needed
BA learning: moving mean heave rises to 22.6271% Hs (OU-III) and 28.8276%
(TFG), with all eight default records failing. Its one refinement retains BA
learning. TFG passes the eight default gates; OU-III still fails small-sea
roll 0.348498 > 0.346670 and has two zero-speed gate misses. No further
scalar or activation tweak of this OU-III mechanism is justified; it remains
experimental pending a new coupled attitude/bias architecture review.

Applying the spectral-law gain restriction to the legacy law fails large
JONSWAP yaw 3.34544 > 3.06465 degrees, pitch 1.28213 > 1.16598 and BA
271.194 > 250.685 percent. Reject a common gain policy across both noise laws;
retain the original legacy gain and qualify that exact production branch.
The failed 96-run production trial is retained as
`tfg_legacy_wave_matched_protected_rejected.json`, never as accepted evidence.

The unprotected TFG S+beta*v candidate fails largest PM-Stokes yaw
1.67407 > 1.64146 degrees, pitch 0.673546 > 0.659650 and 3-D BA error
188.483 > 186.739 percent. This motivates, but does not itself validate,
the attitude-protected implementation. Original S-only unit expectations fail
on the changed model; independent on-manifold finite differences pass after
expressing the correct residual, and the production Jacobian suite passes
the original observation plus protected/unprotected wave-matched observations,
with a separate 3000-step covariance check for each of the three policies.

Earlier rejected mechanisms retain their exact records: vertical-only OU-III
stillness handling and its refinement fail yaw; scalar non-FMA still fails
0.883995 > 0.875540 degrees; oracle CARMA shape fitting raises small-wave
heave 3.550 -> 5.084% Hs. Body-Z-only noise calibration fails OU-III small
moving roll 0.383126 > 0.346670 and BA 82.8951 > 75.6411 percent, and
OU-II medium yaw 0.907622 > 0.903530. Uniform/axis-specific noise scaling
has exhausted its refinement. These are replay failures, not instability.

**Infrastructure/publication (class E).** At source 6552e14, hosted run 43661885009, job
132652794958, `make -C tests/validation evidence-test` fails in
`tools/tfg_comparison.py --check --require-ci` with
`ValueError: TFG comparison source provenance is stale`. The new TFG
source invalidates the older companion comparison, not the separately qualified
paired results. The source hashes must not be restamped onto old scores.
Resolved by full branch evidence run 43661880261, whose companion, validation
and robustness bundles pass all contracts and the 1206-test suite before
publication in commit 280697c. That automatic commit is preserved before
adding the documentation and independently checked cross-build comparison.
The same condition fails hosted OU smoke run 43661884949, job 132653485163.
The final local `make -j1 all` has a completed exit receipt (2587.38 s, exit 2):
all phases through util pass, including all 120 tagged Kalman replay gates,
then tests/validation stops at that same exact stale-companion error.
The later wave_dir and wave_sim groups pass separately. CI companion run
43661880261 succeeds, and its artifact 14140719607 verifies against the current
source and every retained result hash. Importing that exact bundle repairs
the contract; `make -C tests/validation evidence-test` then passes. This is
not a claim that the earlier failed aggregate command exited successfully.
NLO/PII ignore the timeseries suppression flag and wrote about 2.1 GB of
temporary traces, filling the workspace. Completed untracked output traces
were removed; retained inputs, code, scores and evidence were not changed.
The independent local companion diagnostic subsequently completes all 180
rows and validates every input pair and 900-second scoring window. A prior
serial attempt was intentionally interrupted after 14 pairs and restarted
with three replay workers; its partial output is not evidence.
The CI bundle, rather than the local run, supplies the publication receipt.
Its raw CSV keeps the producer's CRLF bytes; whitespace verification treats
CR-at-EOL as line endings without rewriting hash-bound evidence.

Concurrent Eigen builds exhausted
memory (`Killed signal terminated program cc1plus`); compilation is strictly
serial. The two earlier `make all` sessions lost their exit receipts after
exec-server transport loss, most recently partway through OU-III. Neither is
a completed validation. The subsequent TFG-candidate run also lost its session/exit receipt after
completing AHRS, detrend and frequency and partway through heave-only; it cannot
validate the corrected legacy gain policy. The later completed aggregate
run and targeted evidence repair described above validate the final source;
these interrupted attempts remain non-evidence.
A device diagnostic initially named nonexistent `device_config_quiet_test.cpp`;
no score was produced, and the real `tfg_device_heave-test.cpp` is the required
production test. No missing or interrupted output is treated as a pass. An extra test command
failed because pytest is not installed; the three encounter test modules use
stdlib unittest and pass all eight cases with that runner instead.

Dedicated scalar non-FMA lever-arm reproduction passes hosted run 43648218794,
job 132610657908, without changing its 0.01 tolerance. Both modified shared
CARMA method papers compile and their rendered normalization equations were
reviewed. The updated six-page TFG article compiles with the current companion
tables and freshly regenerated deterministic plots; equation and table
layouts were reviewed. Local pdfLaTeX initially failed to find its format;
using the existing explicit LuaLaTeX format repairs this tooling failure. The numerical inventory was refreshed after normalization.

**Retained facts / six-item shipping handoff.** (1) Preserved: all 39 audited
OU-III probe source files, literal states, prediction/correction/reset,
covariance/gains, tuner, scheduler, magnetic service and persistent histories.
(2) Relaxations: ideal sensors/axis removal and frozen Riccati calculations are
explicit non-promoting diagnostics; the latter freezes attitude/adaptation,
omits reset/sync/startup, uses continuous time and substitutes a stationary
OU physical-bias spectrum for the finite random walk. Experimental OU-III
observation/gain changes do not inherit the shipping proof. (3) Failures above
are E; no mathematical counterexample classification is inferred. (4) No genuine
shipping counterexample was found. (5) Established algebra and open physical,
arithmetic and end-to-end qualifications remain controlling; no theorem or
historical replay payload is promoted. (6) Next falsifiable experiment:
finish publication landing and hosted checks, then trace the subsequent
accelerometer attitude/BA corrections and evolving covariance on the same
failed small-PM physical history. The directional diagnostic above rejects
larger perpendicular pseudo-update energy as the explanation; it does not
establish that removing that channel would help. Preserve vertical BA
learning (freezing all BA already failed) and quantify the coupled error path
before proposing another gain architecture. Extend that diagnosis to the
randomized TFG histories; no further scalar or activation sweep is justified.

## Earlier diagnostic record

CARMA trial-rationale audit (class E numerical/interpretation failure): the
provisional multiplier 1.26686008233 was incorrectly labeled as the exact
continuous default-CARMA41 band correction. Infinite-domain spectral quadrature
and an independently augmented Lyapunov solve agree on variance fraction
0.58598012624 and multiplier 1.30634701277 for the stated 0.5/4 band.
The interpretation is withdrawn; the recorded trial inputs, multiplier and
RMS rows are preserved. That historical trial did not replay the corrected factor; the later
source-bound qualification below is separately recorded in encounter_scaling.
Failure of the provisional factor does not disprove the corrected one.
Structures preserved: all shipping OU-III operations and physical histories.
Relaxations introduced: none. No theorem promotion or runtime change.
Next falsifiable work requires an architecture/measurement-model review and
a quantitative argument before another scalar-gain trial, rather than silently
substituting the corrected constant into old evidence.

Hosted lever-arm publication refresh (class E): run 43630161809 on the
PR merge-preview tree reports TFG-CARMA negative-control yaw ratio
38.8156 -> 38.7948 and worst-sea tilt ratio 9.9692 -> 9.9366, exceeding the
unchanged 0.01 absolute reproduction tolerance. The Git commit API confirms
merge preview f8c87c3 and head 49ec862 have the identical tree
50a382ebfa9b51caace1989b79f014521ec6af98. Importing the complete checksum-verified
hosted artifact refreshes raw measurements, summaries, figures and generated
documentation. No simulator, acceptance threshold or OU-III operation changes.
Invalidated hypothesis: cross-toolchain reproduction of these highly amplified
ratios within 0.01 absolute. Retained facts: physical histories, all OU-III
structures, proof payload and open qualifications. Relaxations introduced: none.
Next check: the unchanged hosted comparison against its regenerated publication.

Encounter-response continuation (class E implementation diagnostics): the
39-file shipping OU-III probe src closure remains byte-identical to the audited
reference. The accepted NLO wave-band bandwidth guard is outside that closure.
The whole-src secant receipt is refreshed without changing its numerical payload
or claiming a replay. Structures preserved: every OU-III state, covariance/gain,
reset, scheduler, frontend, tuner and persistent physical/filter history.
Relaxations introduced: none. All theorem/physical-qualification flags remain open.

Failed quantities: the proposed TFG fractional cadence fails three production
legacy-law checks: heave 8.3259 > 8.28795% Hs, accelerometer-Z bias
5.52383 > 5.51543% of maximum true bias, and heave 8.10862 > 8.05605% Hs.
The unchanged limits reject it; the source and accepted measurements are restored.
Suppressing OU-III periodic covariance resets changes worst screened heave
8.23605 to 8.16948% Hs but worsens worst accelerometer-bias RMS
0.05680 to 0.09373 m/s^2. A 0.5 s commit interval also worsens that bias
metric to 0.08371 m/s^2. Neither is adopted. A TFG Schmidt-style held-bias
startup overlay worsens small stationary pitch 0.0458 to 0.5440 degrees and
bias 0.0506 to 0.11685 m/s^2; it is rejected. These invalidate the respective
cadence/covariance/startup hypotheses as general fixes, not a stability theorem.
Prior stillness and all-axis CARMA band-normalization failures remain recorded
in encounter_scaling/kalman_followup.json. None is shipping evidence.

Infrastructure: the earlier optional parallel CARMA compilation exhausted memory;
serial compilation succeeded. A prior make-all receipt was lost after scratch
exhaustion. This continuation stopped a separate make-all after its independently
measured TFG cadence failure and reruns the required command on restored source
with optional time-series output disabled. These are class E failures.
Retained facts: unchanged OU-III closure, native secant payload and acceptance
limits. Current limiter: reduced encounter acceleration with fixed sensor noise,
and coupled attitude/bias/period estimation. Next falsifiable experiment:
use the all-family paired inventory to isolate live tilt/bias/period error at
reduced encounter acceleration. Vertical-only CARMA band compensation has now
failed held-out independent-axis means (2.79308 to 2.80640% Hs; maximum cost
6.46 mm) and is rejected after its one motivated refinement from the failed
all-axis trial. No further scalar-gain tuning of that mechanism is justified.
Full regressions and current-source hosted evidence remain required.

Integration source audit (class E): main 48e444f is merged, and the
OU-II/CARMA stillness fixes from PR #754 are incorporated. The shared defaults
header changes only a descriptive OU-II table comment; stripping C++ comments
leaves identical text. No shipping OU-III operation, state, covariance, gain,
reset, physical history, scheduler or tuner changes. No relaxation, genuine
shipping counterexample or theorem promotion is introduced. Native secant
results and all open qualifications are retained; no probe replay is claimed.
The conservative whole-src and comment-only dependency receipts are refreshed.
Current limiter: source-bound replay bundles need regeneration after the
OU-II behavior change. Next falsifiable check: complete hosted evidence and
native regressions under the unchanged acceptance checks. The first provenance
test exposed the generated-port certificate’s additional comment hash;
regenerating that certificate changes only its defaults-header receipt. A
receipt-key assertion initially used source_blob_sha instead of the actual
source_git_blob_sha and stopped before writing; the corrected audit verifies
that no mathematical output changes.

Stillness-attenuation source identity (class E): the OU-II and CARMA wrappers
no longer attenuate wave variance on the still-water detector's fixed
acceleration threshold (`docs/moving-vessel-regression-diagnosis.md`). This
changes only the conservative whole-src digest of the native secant report
(failed quantity: `shipping_source_tree_sha256`). The OU-III wrapper and
estimator, the probe and all 39 src files of its compile closure are
byte-identical to the merged 5 kn replay head; no replay is claimed. Native
results, byte hashes and all theorem/qualification flags are retained.
Structures preserved: the complete shipping OU-III source closure and every
OU-III operation. Relaxations introduced: none. Invalidated hypothesis: none;
this repairs source identity, not an inequality. Current limiter and theorem
status are unchanged. Next falsifiable check: proof provenance and evidence
tests against the refreshed binding.

Replay-gate toolchain sensitivity (class E, infrastructure): built with
Ubuntu 24.04 g++ 13.3 on an AVX-512 host, the unchanged OU-III replay
measures `gyro_bias_3d_percent` 12.7148 (default) and 12.7029 (interval)
against limits 12.69245 and 12.68089 measured with Debian g++ 14.2; every
other OU-III channel passes. Same source, different compiler/`-march=native`
code generation; the 0.5 % limit margin is below this spread (OU-III up to
0.7 % per channel, TFG up to 25 % on one record's yaw). Invalidated
hypothesis: that the published limits transfer across toolchains. Retained:
OU-III limits and evidence unchanged; no estimator or proof change.
Structures preserved: all. Relaxations introduced: none. Next falsifiable
check: the hosted ubuntu-24.04 replay of this branch against the same limits.

Coefficient-report provenance repair (class E): the study workflow rewrote the
coefficient report without replacing its manifest, retaining an obsolete
dataset name and output hash. The producer now writes both together from the
verified current startup operating points and the pinned physical inputs.
Structures preserved: all estimator equations, physical histories, covariance,
gains, tuning and the diagnostic calculation. Relaxations introduced: none.
No shipping counterexample or theorem promotion follows; established local
results and open theorem flags remain unchanged. Next check: the report/receipt
regression and the hosted study workflow using the new producer output option.

Underway auxiliary publication repair (class E): regenerated v1.3.1 studies
invalidate the old zero-speed engine-guard bound: the 800 rpm guard+R arm
has pooled 3-D/engine-off RMS 2.417253, exceeding 1.7. Removing the gyro
perturbation changes OU-III 3-D RMS by +26.47% and TFG by -39.04%; the
old gyro-path-invariance attribution is also invalid. Tables and figures
now use the verified fresh CI artifacts; the low-rpm regression sentinel
is 2.45, while the inherited 1.7 bound remains for the other conditions.
This documents degraded performance, not a correction to estimator tuning.
Structures preserved: shipping equations, inherited state and covariance,
sensor histories, noise configuration, gains, resets, scheduler and fixed tuning.
Relaxations introduced: none to physical or theorem premises; the revised
empirical sentinel applies only to this declared experiment. No genuine
shipping counterexample was found; established local proof results and
open theorem flags remain unchanged. The sampling-fidelity document pin
changes only the name of the current dataset, not its physical qualification.
Current limiter: final source-bound CI publication. Next falsifiable check:
the complete regenerated replay and evidence contracts, including the cruise
track helper in the TFG source inventory.

CSV timestamp export repair (class E): the input clock correction would be
lost again through default six-significant-digit simulation CSV exports.
Typed/untyped shared runners and the manual AHRS writer now emit time with
17 significant digits, preserving the formatting of every other column.
The complete 240000-sample formatter round-trip test passes; actual AHRS and
OU-II exports each independently pass all 240000 timestamp checks. Both
complete inherited 252-direction native diagnostic words were rerun at
epsilon 0.01 and 0.005; both state-byte hashes match the previous receipt.
Structures preserved: literal estimator equations, state/covariance, physical
histories, gains, resets, tuning, scheduler and sensor records. Relaxations
introduced: none; these source-binding diagnostics do not become uniform
derivative bounds. No shipping counterexample or new theorem claim follows.
Existing local conclusions and open theorem flags remain unchanged. Next
check: complete rebuilt regressions and source-bound CI publication.

Timestamp dataset migration (class E): v1.3.0 collapses distinct 200 Hz
samples onto rounded times after 1000 s. v1.3.1 regenerates integer-indexed
clocks, including the corrected Fenton endpoint/sample count. All 21 input
time-series clocks pass an exact binary64 index/200 audit. The ancillary
tracker spans 30 s (5999 samples), not the main records' 1200 s; the initial
audit incorrectly required 239999 tracker rows and was corrected. The old gyro-bias gate rejects 12.5015% against 12.4956%; the old Kuchler
PM-Stokes gate rejects 14.7245% against 14.14%. These invalidate reuse of
v1.3.0 regression baselines on the corrected physical records. Existing
v1.3.0 results are superseded and must be replaced, never relabelled.
Structures preserved: estimator equations, inherited histories, sensor/noise
configuration and fixed tuning. Relaxations introduced: none. This is an
input/instrumentation failure, not a shipping counterexample or new theorem
result. Existing proof conclusions and open flags remain unchanged. Current
limiter: complete v1.3.1 replay/publication evidence; next check is the full
measured suite and source-bound CI regeneration with corrected clocks.

Moving-vessel CARMA installation baseline repair (class E): the completed
1,456-row installation matrix and 128-row cutoff sweep were compared with
stationary-data ratios. The physical input hash maps differ, invalidating that
comparison. The verified CI artifact now supplies current tables and figures;
the ratio comparison remains mandatory when input maps match. Structures
preserved: sensor corruption, installation, estimator equations, fixed tuning
and the unchanged same-dataset comparison tolerance. Relaxations introduced:
none; distinct physical experiments are no longer treated as repeated runs.
No shipping counterexample or theorem conclusion follows. Existing proof
results remain unchanged. Next check: hosted same-input study reproduction.

# OU-III proof research state

Moving-vessel lever-arm gate maintenance (class E): the old assertion that
navigation error must have an interior band optimum rejects the new 100 Hz
boundary minimum (3-D/CG 0.99923). Force residual instead has the expected
interior minimum: 0.146 at 5 Hz, versus 0.524 at 1 Hz and 3.140 at 100 Hz.
The gate now tests that physical residual; exact compensation, force scaling,
and the deployed 15 Hz accuracy bound are retained. All 31 study/article
checks pass on the verified complete CI artifact. Structures preserved: sensor
corruption, rigid-body installation, estimator equations and fixed tuning.
Relaxations introduced: none; a downstream score is no longer mistaken for
a force-model invariant. No shipping counterexample or new theorem claim is
made. Next check: the same hosted study and source verifier after publication.

Moving-vessel replay maintenance (class E): changing simulation-only source
files invalidates the finite secant diagnostic's conservative whole-src digest;
two historical dataset-scope clarifications also require fresh document pins.
The compiler-discovered 39-header probe closure differs from main f2c65b2
only in WaveFilesSupport.h (CSV horizontal-truth precision); estimator headers
are unchanged. Both complete 252-direction diagnostic words were rerun at
epsilon 0.01 and 0.005, preserving both native secant byte hashes. The source
audit is recorded in moving-dataset-source-scope.json.
Structures preserved: literal estimator equations, inherited state/covariance,
physical probe record, gains, resets, frontend, scheduler and clocks. Relaxations
introduced: no new relaxation; existing off-orbit finite secants remain diagnostics,
not derivative bounds. This is a source-binding failure, not a failed stability
inequality or a genuine shipping counterexample. Existing local results and open
theorem flags are retained. The fresh native receipt is checked against the unchanged source verifier;
no motion/service qualification is inferred from
the new moving-boat dataset. The initial native command lacked clang++; the
available g++ compiler is used with the same diagnostic protocol.

Tangent research quality-gate repair (class E): Ruff found an unused import,
and CI cppcheck found three inherited group-operation name collisions. The
import is removed; tangent geometry now composes the existing physical storage
and two-frame algebra instead of permitting conversion to a different group law.
No operation formula or gate threshold changes. Structures preserved: existing
shipping OU3 source closure and all physical filter equations. Relaxations
introduced: none. The ten-filter 900-s seed-4101 traces are byte-identical and
all 40 score rows agree except host runtime. Latest local cppcheck 2.22 reports
the same 47 pre-existing findings on this tree and producer main; no finding is
added. CI's older cppcheck remains the deployment gate. No mathematical or
shipping instability conclusion follows. Next check: the unchanged hosted gate
and full repository regressions using rebuilt binaries.

Tangent-CARMA integration evidence recovery (class E): the PR evidence-contract
and OU-validation jobs reject seven stale replay dependency records predating the
current timing fixes. The completed full validation, robustness and TFG artifacts
from main run 43501064300 are recovered by their GitHub SHA256 archive digests.
OU manifests are normalized by the existing contract in the original 480d78f
producer checkout, then their immutable replay closure is checked unchanged here.
No fresh replay at this integration revision is claimed. Structures preserved:
all original full-cohort rows, physical histories, producer identities, covariance
and shipping operation results. Relaxations introduced: none. No shipping
counterexample is inferred from this publication failure. Existing theorem flags
and local results are retained; the limiter is evidence freshness. Next check:
full repository tests and hosted publication contracts on the recovered bundles.

Tangent-CARMA integration source-binding maintenance (class E): new CARMA and
Lie-group files invalidate the conservative whole-src digest of the native OU3
secant diagnostic. The 40-file compiler dependency closure is byte-identical to
main c897c92. Both complete inherited 252-direction words were rerun at epsilon
0.01 and 0.005; both native secant byte hashes match the previous report exactly.
Fresh Python/BLAS norm calculations differ only in final floating-point digits.
The invalidated hypothesis is source-digest freshness, not a stability inequality.
Structures preserved: the physical record, all 21 states/full covariance, gains,
Joseph updates, resets, projections, frontend/tuner, forcing, scheduler and clocks.
Relaxations introduced: none beyond the existing diagnostic off-orbit perturbations;
finite secants remain neither derivative enclosures nor reachable-history proofs.
The common-root SPD-ball sufficient-bound failure remains class D; no genuine
shipping counterexample is found. Existing operation identities and qualified
local results remain valid. The all-time physical-service and finite-error theorem
limiters remain open. Next falsifiable check: provenance/evidence tests and full
repository validation with this refreshed source binding and unchanged gates.

Post-merge CI source-binding repair (class E): run 43501064300 completed the
full OU/robustness and TFG studies, then failed publication on the native secant
whole-src digest. The only src change since a053a1d2 is W3dSimCommon.cpp, outside
the compiler -MM closure of planar_native_secant_probe.cpp. Every probe closure
file is byte-identical; the original whole-src digest was independently
reconstructed from that commit before refreshing identity metadata. No replay
is claimed. Native results/byte hashes and all theorem/qualification flags are
retained. Structures preserved: complete original diagnostic execution and
shipping operations. Relaxations introduced: none. This repairs source identity,
not a mathematical inequality; theorem completion and existing physical-service
limiters remain open. Next falsifiable check: proof provenance and all evidence
tests against the refreshed binding, followed by CI publication of fresh bundles.

## Gyro-sampling runtime rollout

Finite rollout selection: simultaneous frontend centering plus Schmidt holding
exceeded the unchanged shared-mode replay attitude/bias bars. This invalidated
the diagnostic assumption that integration accuracy alone preserves every RMS
bar (class E validation/rollout failure). It is not an OU-III theorem
counterexample. Separate core/frontend/holding experiments retain noiseless
rotation/covariance facts and support the narrower core-only default. The
frontend, held-bias policy and calibrated tuning remain native. No threshold
is enlarged. The final gated profiles and full-build result are recorded in
`reports/results/live_bias_consistency/validation.json`.

A separate native-API compatibility attempt also failed (class E): changing
raw CARMA holding defaults broke the historical dense-FD covariance fixture
(error up to 1.335 versus 1e-6) and a sinusoid bound (0.10019 versus 0.1).
Raw holding defaults are restored; the optional Schmidt policy retains its
independent dense-Joseph audit. The restored CARMA suite passes 74/74.

The broader 820-test source audit exposed two stale carried-history native records
and the whole-source-tree fingerprint on the native secant diagnostic (class E).
The invalidated hypothesis was completeness of the first provenance refresh, not
a mathematical bound. The 240-s and 1200-s histories are genuinely replayed;
all-time service/admission and theorem flags remain false. The complete native
secant replay is repeated after freezing source, before rerunning the audit.
A concurrent local build was killed in `bench_carma3_filters_ou.o` (class E
resource failure); serial compilation without other compilers is the next check.
The host warning gate additionally found three missing aggregate member initializers
for `Window::{an,mn,on}` in the test driver (class E build failure). Explicit empty
initializers preserve their existing zero-initialized `FiniteMean` semantics; the
same strict compile flags pass locally. No shipping operation or theorem premise
changes. The broad evidence fingerprint still requires a genuine full replay.
The native default replay also exceeds four legacy yaw fixtures by 0.003--0.007
degrees (class E numerical-contract mismatch); the measured increases are kept.
The historical integrator and every original limit now run under `--regression`.
A separate corrected-default fixture uses the existing ~0.5% margin rule; this
is not a claim of improved yaw or latent-bias recovery. Both profiles pass all
64 native record checks. The all-source publication test exposed three old
source-text guard assertions; the unchanged conditioned-accelerometer route
through the gyro helper is now checked, and all 30 article contract tests pass.
Clang-tidy exceeded its 45-minute infrastructure budget; every TU/check is
retained with a 90-minute budget. No shipping source or theorem premise changes
in these fixture/test repairs.


The live OU-II/OU-III/TFG and CARMA sample paths now default to measured-endpoint
gyro integration. OU-III prediction uses the same effective interval rate in
the mean, gyro-bias Jacobian and covariance kernel. The measurement-only
frontend retains its endpoint rule, with optional centering isolated behind
a separate configuration field. Bootstrap primes the live core;
endpoint lever-arm force terms retain endpoint rates. Invalid IMU packets are
rejected before the shared startup clock/history advances. This is an explicitly
authorized runtime change, not a source-only reorganization.

Structures preserved: OU/S state equations, covariance and gain chronology,
parameter adaptation, magnetic-reference acquisition, BA activation/mean
reversion/projection, scheduler and reset chronology. The prior instantaneous
endpoint-as-interval forcing map is replaced. Previous raw gyro, previous
interval and the sampling-mode/exposure configuration are persistent execution
memory and must be carried through proof-word boundaries. Internal attitude
resets retain physical sample history; a new acquisition clears it.

Relaxations introduced: none in runtime. The finite rotation/replay diagnostics
and the linearized white-noise covariance envelope do not qualify hardware
phase delays, SLOW+FAST physical envelopes, regional capture or end-to-end
stability. Prior source-bound native receipts are archived. The finite native/observer
records are regenerated against the changed complete execution, with the
cache and noise-map audit in `ou3-gyro-sampling-audit.md`. Classification E applies to stale
source bindings, not to the mathematical theorem. No shipping counterexample
was found. Existing algebraic identities for unchanged kernels remain
available; lifting their physical execution qualifications to this sampling
map is OPEN. Next falsifiable proof calculation: qualify the cache-augmented physical
forcing against SF3--SF9 and the still-open construction/release and all-time
magnetic-service obligations, without deleting inherited sample memory.

The final replay publisher exposed a class E archive-fixture failure:
`test_source_archive_provenance.py` copied estimator sources but omitted
`tests/common/GyroReplayContract.h`. The fixture now copies the shared test
headers, preserving the requirement to verify their recorded hashes in source
archives. All six fresh replay shards and the companion study completed on
81d16672; their publication failed at this fixture check, not a numerical gate.
Because the conservative replay fingerprint includes tests, the repaired
revision requires a new full replay before final publication. No source-bound
proof receipt or mathematical result is promoted by this test repair.

## Historical native secant source binding

The following receipt describes its archived source revision; the runtime
rollout above supersedes its current-tree binding. It does not certify the new
sampling map.

The finite native diagnostic was bound to that complete `src/` tree.
The sample-path covariance kernels of the experimental TFG-CARMA cores were
restructured to lower microcontroller processor load
(`src/kalman_carma/Kalman3D_Wave_TFG_CARMA3.h`, with plain-array access in the
discretization helpers of `src/kalman_carma/CarmaVerticalModel.h` and
`src/kalman_carma_shared_modes/CarmaSharedModesMotion.h`), and the series
coefficients of `Exp` and `J_l` were factored out of
`src/lie/SO3Jacobians.h` for both TFG families. That changed the
broad source fingerprint, so the unchanged shipping OU-III probe was rebuilt
and replayed at both registered perturbation sizes. Both raw native byte
hashes and every other field of the committed report are reproduced; only the
digest is rebound. See
`reports/results/tfg_shared_modes/shipping_faithfulness.json` for the receipt.
The QMEKF-CARMA, QMEKF-shared, TFG and OU-III wrapper traces stay bit-identical
to the previous tree over 6,000 samples. The two TFG-CARMA wrapper traces change
by binary32 rounding only: their binary64 cores agree with the previous tree
to print precision over 300 s, and the binary32 differences are of the size
that compiler multiply-add contraction alone produces on the previous tree
(`reports/results/tfg_shared_modes/sample_path_cost.json`). A non-finite gyro
sample still leaves the shipping OU-III wrapper state non-finite; that open
implementation item for a shipping change is unchanged.

Structures preserved: the complete shipping OU-III dependency closure, inherited
execution, covariance, frontend/tuner, scheduler and diagnostic qualifications.
Relaxations introduced: none; the inherited off-orbit finite secants remain
non-promoting diagnostics, not reachable or uniform derivative certificates.
The invalidated source binding is a class E implementation/diagnostic failure;
it does not establish estimator instability. Historical source-identity
records are preserved in `planar-native-secants-source-binding-archive.json`
and the provenance references resolve to their respective archived records.
No shipping counterexample was found; established shipping-level results and
open theorem obligations remain unchanged. Device runtime and peak stack/heap
measurements still require physical hardware.

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

Analytical work only: exact algebra/rationals, rigorous compact inequalities
after derivation, and identity/source regressions. No new replay, secant,
sampled grid or empirical factor. Before substantial calculations state the
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

## Evidence

**E — QMEKF CARMA sample-path storage access rebound the whole-`src/` digest.**
The experimental QMEKF CARMA cores index plain storage in their sample-path
loops (`src/kalman_carma/Kalman3D_Wave_CARMA3.h`, `CarmaCoupledMotion.h`,
`src/kalman_carma_shared_modes/CarmaSharedModesMotion.h`) with every
arithmetic expression and its order unchanged; the fail-closed validator
rejected the native secants' `shipping_source_tree_sha256`. Failed quantity:
source provenance only; classification E. Evidence: none of the three files is
in the shipping OU-III dependency closure; the diagnostic regenerated with
g++ 13.3.0 reproduces both native byte hashes and every field except that
digest (run summaries to 4.4e-16, host library round-off); the QMEKF-CARMA,
QMEKF-shared, TFG-CARMA, TFG-shared and OU-III wrapper traces are
bit-identical to the previous tree, the two QMEKF wrappers also with compiler
contraction enabled. Structures preserved: every shipping operation, constant
and chronology, bit for bit. Relaxations introduced: none. Rebound: the
digest, the secant report's blob and this ledger's blob. No proof flag,
payload or OPEN dependency changes; the current limiter and next experiment
are unchanged.

**E — TFG-CARMA sample-path kernels rebound the whole-`src/` digest.** The
experimental TFG-CARMA covariance prediction, reset transport, process-noise
quadrature and group injection were restructured for processor load in
`src/kalman_carma/Kalman3D_Wave_TFG_CARMA3.h`, with plain-array access in
`src/kalman_carma/CarmaVerticalModel.h` and
`src/kalman_carma_shared_modes/CarmaSharedModesMotion.h` and the `Exp`/`J_l`
series coefficients factored out in `src/lie/SO3Jacobians.h`; the fail-closed
validator rejected the native secants' `shipping_source_tree_sha256`. Failed
quantity: source provenance only; classification E. Evidence: none of the
four files is in the shipping OU-III dependency closure; the diagnostic
regenerated with g++ 13.3.0 reproduces every field, both native byte hashes
included, except that digest. Structures preserved: every shipping operation,
constant and chronology, bit for bit. Relaxations introduced: none. Rebound:
the digest, the secant report's blob and this ledger's blob. No proof flag,
payload or OPEN dependency changes; the current limiter and next experiment
are unchanged.

**E — single 27-state shared-modes filter rebound the whole-`src/` digest.**
Keeping one shared-modes filter removed three `src/kalman_carma_shared_modes/`
headers and changed four there and two `src/kalman_carma/` headers
(experimental blocks removed, comments); the fail-closed validator rejected
the native secants' `shipping_source_tree_sha256`. Failed quantity: source
provenance only; classification E. Evidence: the probe's compiler include
closure (38 `src/` files) is byte-identical to main and contains none of the
nine changed files; the diagnostic regenerated with clang 18.1.3 reproduces
every field except that digest. Structures preserved: every shipping
operation, constant and chronology, bit for bit. Relaxations introduced:
none. Rebound: the digest, the secant report's blob and this ledger's blob.
No proof flag, payload or OPEN dependency changes; the current limiter and
next experiment are unchanged.

**E — 27-state shared-modes headers rebound the whole-`src/` digest.** The
experimental minimal realization of the shared-modes CARMA QMEKF added three
`src/kalman_carma_shared_modes/` headers and changed two there and
`src/kalman_carma/Kalman3D_Wave_CARMA3.h` (an optional per-axis
re-coordinatization hook); the fail-closed validator rejected the native
secants' `shipping_source_tree_sha256`. Failed quantity: source provenance
only; classification E. Evidence: the probe's compiler include closure (38
`src/` files) is byte-identical to main and contains none of the six changed
files; the diagnostic regenerated with clang 18.1.3 reproduces every field
except that digest. Structures preserved: every shipping operation, constant
and chronology, bit for bit. Relaxations introduced: none. Rebound: the
digest, the secant report's blob and this ledger's blob. No proof flag,
payload or OPEN dependency changes; the current limiter and next experiment
are unchanged.

**E — CARMA shared-modes headers rebound the whole-`src/` digest.** The
experimental CARMA work added `src/kalman_carma_shared_modes/` and changed
three `src/kalman_carma/` headers; the fail-closed validator rejected the
native secants' `shipping_source_tree_sha256`. Failed quantity: source
provenance only; classification E. Evidence: the probe's compiler include
closure (38 `src/` files) is byte-identical to main and contains none of the
seven changed files; the diagnostic regenerated with clang 18.1.3 reproduces
every field except that digest. Structures preserved: every shipping
operation, constant and chronology, bit for bit. Relaxations introduced:
none. Rebound: the digest, the secant report's blob and this ledger's blob.
No proof flag, payload or OPEN dependency changes; the current limiter and
next experiment are unchanged.

**E — header reorganization rebound; no shipping operator changed.**
The generic attitude, covariance-projection and schedule helpers moved
verbatim from `KalmanOUCoreMath.h` into `src/kalman_common/`
(`MekfAttitudeMath.h`, `CovarianceProjection.h`, `PeriodicSchedule.h`;
namespace `ocean_imu::kalman::filter_detail`, every name re-exported into
`ou_detail`). OU-III's tuning envelope, r_S law constants and `TuneState`
moved verbatim from `SeaStateFusionFilter_OU_III.h` into
`src/kalman_common/SeaStateTuning.h`, re-exported into `ocean_imu::ou3`. The
CARMA research filters therefore include no OU header.
`Kalman3D_Wave_OU_III.h` is byte-identical. The unchanged fail-closed
validator rejected seven bound blobs, the source hashes recorded in
`planar-service-frontend-binding.json` and
`planar-generated-port-structure.json`, and the native secants' whole-`src/`
digest. Failed quantity: source provenance only; classification E. Evidence:
at the repository's own flags (`-O3 -march=native -funroll-loops`),
normalized assembly (symbol renames and label numbering only) is identical
for all twelve OU-III, OU-II and TFG translation units compared, including
the OU-III simulator, shipping-contract, core-numerics and gyro-projection
tests, and at `-Os` for all eight ESP32-S3 size-probe variants, including
the deployed `SeaStateFusion_OU_III`. The native-secant diagnostic was
regenerated from the new tree and every field except
`shipping_source_tree_sha256` reproduced exactly. Structures preserved:
every shipping operation, constant and chronology, bit for bit. Relaxations
introduced: none. Rebound: the seven blobs, the two certificates' hash
fields, the secant digest and the committed secant report's blob; newly
bound: `SeaStateTuning.h` (shipping) and the three math headers (operation
lemmas). Readers of moved literal text (gyro-bias radius, quaternion
polynomial, `MAX_SIGMA_A` premise, SpectralMSE density) read the new files.
No proof flag, payload or OPEN dependency changes; the current limiter and
next experiment are unchanged.

**D — an inverse-domain estimate alone does not settle the linked receipt.**
The exact failed sufficient step is promoting CR45's covariance-column
payment to `alpha>=dbar` without checking the SAME prefix's threshold.
Cauchy with full Qr proves `d/v>=v/Qr_aa`; monotonicity of
`(4r+1)/(2r+1)^2` gives the allocation ceiling
`alpha*v^2<=Qr_aa*(Qr_aa+4v)/(Qr_aa+2v)^2`.
If `k=dbar*v^2` exceeds it, no upper energy estimate can cross the
allocation threshold. This is a conditional failure of that sufficient
method, not a proved occurrence or negative mode in the admitted family.
Retained facts: CR44--CR45 and every linked full-process/transverse,
conditional/regression/acc/S/magnetic term; the ceiling is NOT an upper
bound on their total. Current limiter: no quantified retained-action
comparison excludes all feasible CR30 chi>=1 modes, and no common root
reserve is established. Next falsifiable test: derive that comparison
from CR5/CR31--CR32 and the literal applied columns on the actual inherited
image, retaining first marginal, target lag and every prefix guard; only
then test the complete CR43 remainder at the same reserve. No sign,
rank, inverse, OPEN dependency or broader flag is promoted.

**E — dependency installation failed before compilation.** Exact head
`5ecbe504d0d9` push build run 37644598944 failed jobs 112872753415
(imu_calibrate) and 112872753700 (OU-III): `apt-get did not succeed after
3 attempts` at Install build dependencies. The Ubuntu mirror download
timed out; no compilation or proof check failed in those jobs. Retry only
failed jobs; retain successful checks and inspect each new exact head.
The PR build still running rejected individual job reruns; do not infer
a successful retry. A transient GitHub fetch 502 and an overwrite attempt
on read-only Git objects were recovered by fetching the complete clone
and copying only missing objects. The full automatic evidence commits
`070e7de3633f` and `944f7bdffa37` are preserved; no scientific payload is
manually replaced. At `944f7bdffa37`, build run 37657130742, OU-II job
112925147621 again exits before compilation: `apt-get did not succeed
after 3 attempts` while downloading texlive-fonts-extra. This is E,
not a proof failure; retry failed jobs after the workflow completes.
The mistaken test-file path was corrected by repository search to
`tests/validation/test_ou3_measurement_frame.py`; no result used that read.

**E — exact-head two-column render counter repaired.** At `17c53af4b5ac`,
PR run 37643337778 / job 112867577478 fails `lualatex` with
`LaTeX Error: Counter too large` at the final shipping-audit paragraph.
The added paragraph exceeds IEEEtran's alphabetic paragraph counter;
the article-class wrapper does not test that display qualification.
Use numeric paragraph reference/display forms inside the already scoped
appendix group and split the new displays to fit columns. The arguments,
equations, status and shipping operations stay unchanged. The next test
is the exact-head IEEEtran two-column render, not another proof estimate.

**CR44--CR45 — quantitative nonzero receipt payment; uniform sign OPEN.**
The SAME L_rem in CR37 pays alpha*mu^2 with
`alpha=(4d+v)/(v*(2d+v)^2)>0`, `X=F P F'`, `Qr=Q_ship-q_aw ee'`,
`v=X_aa`, `d=(X e)' Qr^-1 (X e)`. Its full retained noise is at least
Qr/2; the transverse Fisher square is retained. This prices every actual
nonzero receipt with full integrated Q, earlier deletions and actual prefix
covariance, without the coarse kernel reserve. Substitute L_rem once in Hbar;
all acc/S/magnetic action and complete endogenous cross cost stay linked.
The direct allocation crosses its individual threshold only if
`k=dbar*v^2<1` and `d/v<=sqrt(1-k)/(2*(1-sqrt(1-k)))`.
No uniform bound on the actual covariance-column energy
`d/v=(P e)' F' Qr^-1 F (P e)/P_aa` crossing that threshold is proved.
The remaining grouped word action may pay its residual charge. Pointwise
alpha>0 does not exclude ALL feasible chi>=1 modes: promoting it would be
D, insufficient domination. Receipt blocker NO; full Schur NO; homogeneous
margin NO; no OPEN dependency discharged or count/flag promotion. Structures
preserved: same complete process/face/acc/S/magnetic prefixes and common tangent.
Relaxations: none in this inequality; existing real partial-fibre/cone scope.


**E — merged-main binding and warning repair.** Fresh main is `a47da994`;
new work uses `proof/ou3-receipt-schur-continuation`, never the merged branch.
Exact reproduction confirms the two reported binding failures. The native
profile still includes the same 34 shipping files: 33 are byte-identical;
the Mahony constant declaration changes only whitespace. Compiler preprocessing
differs only by that exact line-wrap replacement, proving token identity under
the unchanged C++20/O1/Eigen profile. Refresh the two affected bindings with
all diagnostic payloads, false flags and original audits retained. Main's
calibration/compensation/source-selection behavior stays intact. The separate
warnings job also rejected uninitialized output arrays in the magnetometer
test; value-initialize those test buffers, retaining its success and accuracy
checks. These repairs change no shipping operator or proof margin. Next:
reproduce evidence and the affected warning/test commands, then attack the
same nonzero receipt inequality and full endogenous cost.

The incorporated main restores the original Mahony declaration wrapping.
Its native dependency closure and full exact-profile preprocessing match the
earlier valid baseline byte for byte; rebind the affected source fields with
the mechanical refresh and all earlier audits retained. No original finite
payload or upstream automatic evidence is replaced. Required parallel make
failed by launching tests before binaries existed. The first serial command
exited 2 with no final diagnostic in its saved log, whose last line reports
passing TFG Jacobians. Classification E; no proof or shipping behavior is
invalidated. The subsequent serial `make all` at `0779e7d2a61d` completed with live
exit 0 and 1123 validation cases; CMake configuration was skipped because
CMake is unavailable locally. Do not credit the incomplete invocation as
a pass or transfer that completed run to a later head.

**E — status-path lookup corrected during handoff preparation.**
The command searching `docs/ou3-theorem-status.md` failed with
`No such file or directory`. The maintained record is
`reports/results/ou3_stability/theorem-status.json`, with its generating
source `tools/stability/ou3_theorem/theorem_status.py`. This invalidates
only the guessed documentation path. The current status and manifest were
then read; all proof flags/counts and CR41--CR43's unpaid sign are retained.
Use those actual paths in the next focused receipt/Schur calculation.

**E — current-head workflow admission and source-binding repair.**
At `ca3bd78d255d1e1b1e5eaceadf4c3f405aa50946`, all six PR workflows
reported `action_required`, with no check runs; the sole commit status was
successful. Thus prior-head CI success does not certify this automated head.
The PR was open and conflict-free, but main had advanced to `8cecdf39`.
Incorporate that main without dropping either history or any of the newer
20-file evidence bundle. Its full-src digest is
`94f68621549fce3b83bd04652f3ce6cf9789a72173eccffba4cef27b04272222`;
the preceding diagnostic binding was
`9186fa4aa6133ad3822daaec94466e252a5dfb3d89a31be6a5091e871b22f592`.
The exact-profile compiler include closure still has 34 shipping files,
all byte-identical to immutable `108ce145`; the eight source paths changed
since that original are outside the closure. Refresh only the provenance
binding, preserving the previous audit and all finite diagnostic results.
The literal OU-III sketch now calls `configureAtomS3RMag`, as main requires;
update that audited sketch binding separately. No replay or gate change.
This is an infrastructure/provenance failure, not a shipping counterexample.
No mathematical premise or proof flag changes; the AW and complete-gap
blockers remain OPEN. Required next verification: `make all` and exact-head
CI after publication. The next mathematical calculation remains the SAME
feasible linked `R_W(c)` sign, with the complete CR43 Schur cost retained.

**E — source-identity repair of merge CI; no mathematical promotion.**
The PR-head proof checks on `45efd074` passed, but its merge checks failed
with `finite quotient/secant evidence verification failed: shipping source
tree changed`. Main added two AtomS3R diagnostic headers. The actual native
secant probe's compiler include closure has 34 shipping dependencies; every
dependency has identical bytes to the preserved PR source, and neither added
header is in that closure. The old report's whole-src digest also matches
that immutable PR source. Refresh only the full-src provenance binding,
retaining the old binding, exact dependency hashes and explicit no-replay
record in the existing diagnostic. All diagnostic payloads and false
qualification flags remain. The validator and its full-src gate are unchanged.
Resolve generated-evidence conflicts with the complete newer PR bundle,
including its nonconflicting transition file, and retain main's diagnostic
code. The source-identity argument is a CI repair, not an AW sign proof.
Main then advanced through PR #658. Its additional source edits and deletion
are also outside the same 34-file include closure. Recheck their bytes,
retain the earlier refresh in the audit, refresh the literal OU-III sketch
binding and retain main's low-noise calibration/runtime setup. This does not
verify its magnetic residual envelopes or strengthen MAGNETIC SERVICE.

**Full endogenous receipt/kernel elimination — CR41--CR43; sign OPEN.**
Insert the SAME reassembled Hbar and all remaining linked process/S/mag/
acc utilde action into the original full homogeneous block. With
`K_c=lambda*Hbar-c*J_r_OO`, independent actual receipt rows `L=E*L0`,
`Z_c=L0*K_c^-1*L0'`, `Y_c=K_c^-1*L0'*Z_c^-1` and
`K_c_sharp=K_c^-1-Y_c*Z_c*Y_c'`, eliminate the paid odd receipt kernel.
The single existing remainder is
`R_W(c)=[[S_c,Y_c'*X_c],[X_c'*Y_c,A_c-X_c'*K_c_sharp*X_c]]`,
`S_c=Z_c^-1-lambda*E'*Dbar*E`.
The quadratic identity is lossless on a justified linear inherited slice;
on scalar cones the minimizing kernel component must retain actual guards.
The first full process qualifies K_0 at each fixed word, not a uniform
positive-reserve inverse. No tiny restricted-kernel margin prices the cross.
If S_c>0, the actual full Schur remainder is exactly
`T_c-B_c'*S_c^-1*B_c`, with `B_c=Y_c'*X_c` and
`T_c=A_c-X_c'*K_c_sharp*X_c`. Metric cross entries, generated suffixes,
linked process/source and cross/square work, reset/projection, inherited
auxiliary/reference and endpoint gauge are all retained.

The remaining linked action has NOT decided the receipt corner. A possible
feasible S_0*t=0 mode reconstructs y=Y_0*t, mu=E*t, with positive first
process action balancing the receipt. It does not force transported service
columns to vanish. Even such an odd zero would require B_0'*t=0 for full
linear-slice positivity. No actual such mode is exhibited or excluded.
Existing zero-action and physical service premises do not supply that
quantitative sign. AW blocker: NO. Complete homogeneous gap: NO.
This is a full-gap reduction, not the requested positivity breakthrough.
No uniform compact forward domain, inherited-image closure or flag/count
promotion is assumed. One remaining comparison: feasible R_W(c)>=0 for
one common c>0, with its existing inverse and cone qualifications.

**Joint process/face/acc receipt payment — CR36--CR40; sign still OPEN.**
Complete the existing dC and entire same transported dT jointly, keeping
the linked acc utilde square and all process/S/magnetic work. For each
actually paired triple this gives the strictly lower negative weight
`dbar^-1=A^2-C_star^2+eta*A^2/[omega+eta*kappa/(2*gamma*A^2)]`,
where `gamma=1/C_star-1/A>0` and
`kappa=ell' B_face^-1 ell<=r`, with `ell=U_T B_acc k'` from the actual
intervening S/reset transport. There is no division by kappa, no restriction
to zero receipts and no independent minimization of dB or a loss coordinate.
Insert immediately as `G_OO/lambda=Hbar-L_pmu' Dbar L_pmu`,
`Hbar=H_pmu-L_pmu' (D_pmu-Dbar) L_pmu`; other weights stay unchanged.
The first retained full process still justifies a fixed-word zero-reserve
inverse. Symbolic expansion verifies the completion and strict weight
difference; no numerical contraction search was performed.

The price reduction does not exclude the balanced shipping direction:
`H y=L' D mu` iff `Hbar y=L' Dbar mu` at threshold one.
No particular eigenvector was previously proved feasible, and none is
exhibited or excluded here. Treating the smaller weight alone as a positive
margin would be D: the positive Gram has decreased on the same reader
directions. The ONE sign obligation remains the feasible CR30 eigenmode
with `(Hbar,L_pmu,Dbar)`, followed by the full endogenous CR17 cost.
AW blocker: NO. Complete homogeneous gap: NO. No flag/count promotion.

**Constrained dangerous direction — CR29--CR35; no margin promotion.**
The requested root-normalized minimization now has its actual constrained
Euler equations, including inherited linear constraints and scalar face
guards read from the SAME prefixes. On each justified closed fixed-word
history-tangent slice, the first retained full process makes H_pmu invertible
at zero reserve; this is inserted immediately into CR28. It does not give a
common positive reserve or assume full covariance/image compactness.

The sign test reduces losslessly to feasible receipt eigenmodes:
`A_F t=chi t`, `A_F=D^1/2 L_F H_F^-1 L_F' D^1/2`,
`y=H_F^-1 L_F' D^1/2 t/chi`, with the actual inherited image and all
branch guards enforced. A nonpositive gap exists on a closed slice iff
such a reconstructed feasible mode has `chi>=1`. All dependent receipts
remain. No actual eigenmode is asserted feasible or excluded.
CR31--CR32 give its simultaneous dB/dT/dC/mu and differentiated full
integrated g/Q_ao relations; actual transported S/acc/magnetic columns
are retained. At saturation the mode is the inverse of the EXISTING
H action applied to adjoint-transported AW receipt rows (CR35).

**D — the automatic near-kernel contradiction fails at a precise step.**
A zero signed gap gives `y'H y=mu'D mu>0`, not zero positive action.
Its first process loss is strictly positive; correction losses and the
completed square need not vanish. Existing zero-action/kernel theorems
exclude zero action with nonzero receipts, but do not exclude this balanced
receipt-driven mode. Treating its individual losses as zero would assume
the disputed domination. No shipping-negative tangent or impossibility
is established. The ONE immediate remaining blocker is to exclude/pay
the reconstructed feasible `chi>=1` mode in CR30; its saturation is CR35.
The full CR17 cost is preserved for substitution after a proved odd margin.
AW blocker: NO. Complete homogeneous gap: NO. Counts remain 25/4/1.

**Nonzero receipts with the literal preceding process — CR25--CR28.**
The full four-kernel integrated covariance, including Q_ao, has normalized
inverse AA entry 8. Its existing source-defect/Volterra comparison gives
Q_ship >= 2*q_aw*E_aw*E_aw', with a uniform allocated q_aw
`248188600375173268173354737/4768491208509685862400000000000000`
>5.2e-8 on the existing regular real isotropic profile. This is an
algebraic split of ONE actual prediction, not a new covariance recursion.
Combining its Fisher loss with the subsequent face proves CR26's joint
signed inequality for all actual nonzero receipts. With C_star=C-q_aw,
A=Cplus, the retained expression is
`L_rem+(dC/C_star)^2-((dC-mu)/A)^2+2*(1/C_star-1/A)*dT B dT'`.
Its exact completion keeps the shared square and has effective gap
Delta+q_aw>=q_aw; no division by Delta is required. Before a first face,
the real inherited receipt remains; queued/committed targets are distinct.

Nonzero scalar zero-gap directions now use their actual cones:
mu<=0 deletes AA, mu>=0 is the identity; their maps agree at mu=0.
Each guard uses that branch's actual previous prefixes. No receipt kernel,
nonzero rank or independent covariance origin is assumed. The joint
completion is inserted in the original whole odd matrix as
`H_pmu-L_pmu' D_pmu L_pmu` with each preceding process counted once.
Actual acc/S/magnetic losses and all other prefixes remain. The unsplit
Q_mu/L_mu rows and CR17a's signed matrix are unchanged.

The remaining negative receipt form still needs uniform same-root service
payment (CR28), and then the full endogenous cross Schur cost. In an
actually feasible dB=dC=0, dT=zeta*T prefix its allocated conditional
sign depends on `A*(Delta+q_aw)-2*beta*C_star`; residual full process
and actual transported service action must pay any deficit. Feasibility
is not assumed or excluded. No negative shipping word is produced.
The fresh amount is not c_O or c_complete. AW blocker: NO. Complete
homogeneous margin: NO. No existing full-word OPEN dependency is discharged;
all uniform margin/supply/activation/nonlinear/arithmetic flags stay false/null.

**Actual acc/floor coupled payment — CR18--CR24.** The literal same-cycle
floor--S--acc--mag order and the planar acc AW coefficient one supply a new
payment inside OF4. After actual S/block resets, C_acc=Cplus_face and
dC_acc=-d(beta_face) on the SAME covariance tangent. Exact completion of
the applied acc Fisher loss proves `L_P=utilde' K_u utilde+eta*(t-b)^2`,
with `1>eta>Cplus/s>0`, `b=d(beta_face)/Cplus`, and
`t=dT_acc B_acc (h_o+T_acc)'/Cplus`. Thus the combined face/acc work is
its original positive face action plus `utilde'K_u utilde+eta/omega*t^2`
minus `omega*(b+eta/omega*t)^2`, where omega=1-eta>0. No independent
receipt, released q, scalar maximum or new storage is used.

This pays a strictly positive fraction of the original b-square coefficient,
but retains its exact endogenous b/t cross cost. Paired acc loss is removed
once from H_W. The original full odd gap becomes
`G_OO/lambda=H_acc-R_acc' Omega R_acc`; actual process/S/magnetic action,
earlier AA deletions, unpaired faces, queued target lag and full integrated
Q remain. The unpaid threshold is
`Omega^-1-R_acc*(H_acc-c_O H_0)^-1*R_acc'>=0`, with its inverse qualified.
Uniform innovation/reader coverage on the actual inherited image is NOT
proved. At k=h_o+T=0 the surviving weight is [R/(R+Cplus)]^2, and the
complementary `T dB T'` reader remains. Actual S/magnetic covariance-column
payment, generated through the nonzero integrated g and linked Q_ao in
CR24, must cover it; the physical two-column service floor is insufficient
as a certificate. No shipping-negative direction is asserted.

CR23 retains the existing >1e-37 whole signed reserve on the ACTUAL
nuisance-supported image intersected with ker R_acc and zero-gap
compatibility, after 17 s default regular A21 with the first prediction.
Unlike a zero-receipt assertion this can pay nonzero mu, when the actual
linked-reader constraint holds. The intersection may be empty; its rank,
uniform parameterization and full cross Schur cost remain unproved.
The odd AW blocker and complete homogeneous gap remain NO; no OPEN
dependency, count, source supply, capture, radius or arithmetic flag is
promoted. Exact-rational tests check the whole matrix polarization using
the existing prefix generator, not an empirical contraction search.

**Current shared-tangent AW attack — the odd blocker and homogeneous gap
are NOT closed.** CR12 defines the actual construction-history image
intersected with the existing pure odd-covariance fibre. Its partial
zero-energy kernel test passes by differentiated determinism; a uniform
parameterization and nonzero rank are still unproved. No inherited
auxiliary column is deleted outside this fibre.

Literal complete prefixes, including every active AA deletion, give
`G_OO/lambda=H_mu+Q_mu' L_mu+L_mu' Q_mu-L_mu' L_mu` (CR13a).
All q, mu, dT and correction losses use the SAME root tangent, actual
integrated F/Q, applied acc/S/magnetic rows and queued target chronology.
The unsplit formula is valid without dividing by Delta. CR17a eliminates
these common-root coordinates exactly; its unpaid receipt matrix is
`Z^-1+QY+Y'Q'-E'E-Q K_sharp Q'`. Its inverse qualification and uniform
nonnegative sign are not proved by the physical two-column service floor.

The strongest strict result is (CR16):
`y'G_OO y >= gamma_n*lambda*y'H_0 y`,
`gamma_n in (1.37935,1.37936)*10^-37`, on the ACTUAL nuisance-supported
image intersected with zero target-relative receipts at every active or
zero-gap face, at pre-prediction
roots after 17 s of the existing regular default A21 qualification.
Process information loss pays this whole signed word because face work is
nonnegative on that intersection. Inactive/absent faces do not reset it.
Its actual dimension is unknown; an empty intersection proves no strict
sub-block. Other directions remain in the full matrix. Substitution gives
the exact remaining test `A_c-X_c' N_c^-1 X_c>=0` (CR17), with metric
cross entries retained, for one common `0<c<10^-37`. No affordability of
that cross cost or absorption of nonzero receipts is inferred. The nonempty
word includes the first prediction. At a zero-gap face the actual AA
directional map is `D_aa_plus=max(D_aa,0)`, so active receipt rows alone
do not qualify that crossing. All partial matrices use the same
zero-gap-compatible image. Other crossing directions remain unresolved.

**LR8/LR9 remains controlling.**
Appendix app:correlated-complete-gap constructs the exact block
[[G,H],[H',-S]] by polarizing the actual complete prefixes and pulling back
both endpoint gauge energies. Only actual admissible physical/source
columns may enter d. Inherited endogenous auxiliary columns are folded into
the justified graph M+R_a Lambda or retained in the endogenous augmented
root. A zero reverse locked-reference edge does not set that graph to zero.
An augmented root with zero stored energy must satisfy (CR1)'s terminal
kernel test; the uniform inherited origin/graph qualification is not proved.
The homogeneous test is the full relative Schur remainder (CR3), before
source supply. Completing with G alone would consume the entire contraction
margin; a later supply must reserve G-theta J_root.

The unpaid AW-face work is substituted into that same signed gap using the
literal marginal receipt (CR4)--(CR5). After an actual active fixed-target
face on the qualified pure covariance fibre, its target-relative tangent is
mu=-sum w_j dI_a,j. Each scalar acc/S/magnetic correction has the sharp
same-operation directional Fisher reader
(dI_cov)^2 <= I_a*(2P_aa-I_a)*L_P.
Row/noise derivatives remain signed causal ports. General inherited
tuner/target/source variation remains in (CR5); it is not set to zero or
called independent external forcing. The queued floor target may lag the
newly committed process target. Actual Q_aa includes its small-step
polynomial and PSD repair; exact stationary variance is not imposed.

Exact face completion (CR8) retains dC, dT and mu together. Freeing their
shared-tangent correlations gives the attempted Schur condition (CR9).
Its one-correction fixed-target specialization has relative charge
(2alpha-I)/(2alpha-I-2beta)>1 for beta>0. **D**, sufficient-bound failure:
this analytic released-coordinate test is not a reached shipping
eigen-direction or a service-admitted word. The sharp reader cannot be
tightened to fix that relaxation. Returning to the literal S/magnetic
equations (CR7) shows the omitted link: dI and dC depend on the SAME dB,dT.
That link is now inserted along with actual integrated process and service
directions, but its uniform positive Schur remainder remains unproved.

No numerical search, new replay or new storage is used. No full-word OPEN
dependency is discharged. The complete homogeneous margin, admissible
linked supply, finite-error/prefix/arithmetic transfer and activation are
still OPEN. Counts stay 25/4/1 and all promotion flags stay false/null.
Intervening automated evidence 0e8313729b5605b561e9a145c03e88fd25796478
was fetched and incorporated without changing the proof edits or shipping
source. The general service-conditioned result and special planar
all-time admission remain separate.

**Direct signed-gap attack — analytical, no new storage:** Appendix
`app:odd-covariance-signed-gap` writes the common-coordinate complete matrix
(OF1) and substitutes the combined process/source/generated matrix (OF1a).
Its `z'Lambda z=E_i+S_i-2C_i` stays inside the cross/square work. An independent
uniform E_L cap is not imposed as a prerequisite for this matrix calculation.
The existing complementary-base Schur absorption is inserted with coefficient
one. All even/source/reference/reset/AW/gauge cross terms remain.

On the qualified central-planar regular A21 **pure odd-covariance partial
derivative**, the actual causal equations give zero mean, score, coefficient,
reference and non-AW generated self-ports. This is not the claim that arbitrary
odd mean or even BG variations have zero coefficient work. The ambient block
has 45 symmetric coordinates; restrict it to the actual justified causal
image, preserving inherited root/auxiliary correlations. Full E_L may be
nonzero, but its reader on this partial block vanishes exactly.

At an actual active AW_y face, in the existing conditional coordinates,
`beta=T B T'`, `C+=C+Delta`, and `dC+=-d(beta)` for this fixed-target fibre.
The **same** Fisher decrement is

```
(W-W+)/lambda = (dC/C)^2
  +2*(1/C-1/C+)*dT B dT' - (d(beta)/C+)^2.
G_OO = lambda*(H_W-R_W'R_W).
```

H_W contains the complete transported process, acc, S and actually applied
magnetic Fisher losses plus the positive face terms. R_W stacks the actual
regression rows after every preceding face derivative. The directional
mag/S loss is `2*h D J D h'/S-(h D h')^2/S^2`, retaining actual whitening,
covariance and chronology; no weakest-eigenvalue service replacement is used.
The exact threshold is `I-R_W*(H_W-c H_0)^-1*R_W' >=0`, when that inverse
is legitimate. Neither its uniform margin nor the full cross-block Schur
margin is proved.

**D — exact obstruction to automatic face absorption:** conditional precision
nonincrease does not pay this regression square. Even a root with D_AA=0
can fail after an anisotropic non-AW correction and pending sync. The exact
formal word in (OF6)--(OF7), including a positive process Q=I/1000, has signed
gap `-691428110973/567390082009` and relative reader
`3600981887018337723830924642750000/2360546894727403698724996038046493 > 1`.
This tests the sufficiency of the available matrix premises. Its F=I,
three-state process and observation are explicitly **not shipping-reached
or service-admitted**. No shipping counterexample or instability follows.
The missing shipping comparison is now the actual AW regression reader
against the same integrated process and transported S/magnetic action.
No machine-readable full-word OPEN dependency is discharged; counts and
all promotion flags remain unchanged.

**PROVED — analytical, actual causal reference reduction:**
`app:locked-live-causal-reference` proves that on the magnetically locked Live
default-Complementary continuation, fixed inherited auxiliary state and
delivered history give identical reference/offset executions for finite MEKF
root variations. The active reverse reference port and its homogeneous word
work are exactly zero. The Live gravity gate still reads MEKF attitude, but
its acquisition guard is inactive while lock holds. Refinement uses private
Mahony; its yaw write remains an actual MEKF jump. Cold/reacquisition,
watchdog branch qualification and external reconfiguration are not covered
by deleting that gate. This is a proved causal reduction, not an assumption
that every reference is exogenous on every chronology.

The same-source affine canonical map gives
`|B_ref-B_anchor| <= ||K-K_anchor||*M_mag + |K*b-K_anchor*b_anchor|`.
Convex slew and the actual accepted-matrix/offset gates prove, in the default
real lift after acquisition, `|B_ref|, |K*(m-b)| <= (9/5)*M_mag`. This supplies reference
amplitude finiteness without assuming comparison retention. It does not supply
alignment, acquisition time, horizontal derivative denominators, finite
arithmetic qualification, or an affordable complete-word charge. The complete
source differential groups shared statistics before norms (LR6); its gain,
Joseph and generated suffix work and endpoint gauge terms remain in (LR8).
The actual homogeneous MEKF gap and grouped forced-work domination remain
OPEN. No radius or count promotion.
Inherited auxiliary/root correlations are retained: (LR9) uses
`M_lift=M+R_a*Lambda` when the actual causal graph is proved, otherwise the
full augmented lift. The zero reverse edge is not zero total derivative for
an auxiliary state inherited from a MEKF-dependent acquisition.

**PROVED — analytical, restoration audit:** `app:bordered-causal-restoration`
proves `W_lambda,0-W_lambda,N=(c_0/2)*(A_B-F_B-R_W)` for the full causal
derivative, with `lambda=c_0/2`. The negative budget-derivative square in
restoration cancels the square already stored at the terminal point; it is
not additional dissipation. Paying restoration and all causal work is exactly
the original fixed-weight gap, not a weaker obligation. The actual causal
lift admits the necessary-and-sufficient Gram criterion (BR2); no positive
uniform constant for shipping has been proved. A fixed-coefficient rational
example propagates one nonzero root tangent with `dP=dc=0` and metric ratio
`2/(c_0+9n/7)`, so causality alone cannot supply that constant. This is a
formal scope test, not an admitted shipping history. No radius calculation.

**PROVED — analytical, new storage test:** `app:bordered-comparison-storage`
lifts the actual comparison into `M=[[P,e],[e',kappa]]`, with positive slack
`c=kappa-e'P^-1 e`. Its Fisher storage is exactly
`tr(J dP J dP)+2 eta'J eta/c+(dc/c)^2`. Supported prediction is a PSD addition
and optimal correction is a joint-covariance Schur complement. With fixed
process/row/noise/comparison-defect operands both operations dissipate this
ONE storage, absorbing their full covariance-induced score work without a
separate residual or comparison maximum. This is a proof lift, not a runtime
observer. Causal dF/dQ/ds and dH/dR/d(delta), plus reset/projection/AW and gauge
terms, remain exact signed work. Held masks require the reached reduction.

**Not a retention theorem:** slack obeys `c_N=c_0+B_W`. An exact formal
process/correction word keeps `P=1/2,e=1,V=2` but adds `9/7` to slack each word;
the mean weight tends to zero. No shipping reachability is asserted for that
scope test. Restoring fixed slack costs exactly (BC8), including its negative
`(dB_W)^2` term; it is not a free proof reset. Uniform equivalence and actual
causal-port absorption remain OPEN. Do not replace the original Schur-gap
obligation with monotonicity of this degenerating metric.

**PROVED — analytical:** the continuous exact MOVING physical fibre has identical
complete records and preserves the stated MARINE/SLOW+FAST conditions. The pair
satisfies max(V+,V-) >= 1600(g sin(theta))^2 > 15.38649 on the stated SPD branch.
The physical angle cap and curvature bound
`g^2*(beta^6/36+beta^4/4)` are not precision-normalized nominal-error bounds.
Same-record physical variation leaves the entire deterministic estimator/P
execution identical; this is different from perturbing its nominal root.

**PROVED — analytical:** the corrected local operation calculus retains dR,
dK*r, same-model dr, reset differential and the exact finite gain remainder.
Held-BA accelerometer S uses BA uncertainty/cross terms, while PCt omits BA
columns and masks BA rows. Its literal Joseph differential is separate from
the unmasked optimal Riccati formula. Exact dual-number and rational identities
verify these noncommuting formulas without finite differences.

**PROVED — analytical, explicitly qualified implication:** the private frontend
now has source-derived binary32 error charges, rather than assumed caps.
The bit-seed norm argument extends to squared input norms in [1/4,128].
A rational guard/seed induction and the literal new-integral-before-proportional
operation order close the pitch/integral cell `||z||_G<=.006` at every prefix.
The corrected rate is below .0067; the derived pitch charge is below 46e-9 and
scaled integral charge below 2e-9. This uses stated RN/sqrt/exp qualifications;
actual target-toolchain libm qualification remains OPEN. It is not a MEKF gain.

**PROVED — analytical reference-cone implication:** within that private domain,
exact normalized tilt and nonnegative same-history reference averaging give
`a*=11/5000`, `Bx>=75*(1-u)*(1-a*^2/2)`, `|Bz|<=75*(1+u)*a*`,
`||B||<=75*(1+u)`, where `u=2^-24`. Delivered magnetic-sample rounding is
retained even before charging reference arithmetic. For the exact-real
unrounded record the earlier norm bound 75 is unchanged.
Exact planar hard-iron symmetry retains its zero information direction.
Float reference accumulation, continuous-offset arithmetic, finite acquisition
and the complete inherited nominal MEKF domain remain OPEN.

**PROVED — analytical:** on each unmasked optimal correction, the covariance
part of `dK*r` satisfies `||f_P||_(C^-1)^2 <= (NIS/2)*L_P`, where `L_P` is
that SAME operation's covariance Fisher-storage loss. The proof retains linked
P/H/r and has exact constant 1/2. It is not applied to the held-BA branch.
For the complete joint quotient/covariance word, the exact storage identity
retains the shared gauge/source forcing and all cross terms. If its matrix gap
is uniformly at least `epsilon*J0`, Young's inequality yields rho=1-epsilon/2
with an explicit linked supply matrix. **No positive uniform epsilon is proved.**

**PROVED — analytical identities and qualified first-boundary implication:**
`planar_innovation_storage.py` proves the exact variational NIS identity
`NIS=min_z z'P^-1 z+(r-Hz)'R^-1(r-Hz)`. The default accelerometer rotation
remainder gives `NIS<=V+(kappa V+d)^2/R_min` when the same actual P proves
`B<=kappa P^-1`. At the literal first-Live boundary, the central beta=0 fibre
has `V<21/1000`, retaining its small nonzero BA and inherited physical S.
The bound is phase-uniform exact rational algebra in the real-operation lift,
conditional on the retained private-angle/reference cone and handoff. The
boundary-operand acc/mag NIS bounds are `<21/1000` and `<13/200` respectively.
They do not bound the next prediction/correction, prove finite acquisition,
or qualify float reference/covariance arithmetic. They are not an A21 restart.

**PROVED — analytical reachable H18 structure:** default hold decouples BA;
prediction, masked corrections, reset and AW preserve the zero cross blocks
and constant held BA covariance. The active 18-state Joseph block is optimal
with `R_eff=R_acc+P_BA`, retaining `dR_eff=dR_acc+dP_BA` and the held BA
comparison/source port. This does not apply to arbitrary held cross blocks
or discard the three BA coordinates after release.

**PROVED — analytical coupled storage identity:** the invertible tangent shear
`eta=de-dP P^-1 e` cancels covariance-induced `dK*r` exactly for fixed-row,
fixed-noise additive corrections. Its full identity retains dH/dR and the
literal `dr+H de` mismatch. The comparison must be inherited; resetting it
would destroy the cancellation. Prediction has its own exact shear formula
including physical OU mismatch, dF and dQ; its fixed-coefficient zero-source
port has the linked charge `L_mean*L_cov/2`. Actual source terms remain.
Uniform coercivity requires a bound on the comparison precision energy.
The word balance below now composes reset/AW/projection and moving quotient
ports; their uniform absorption remains OPEN. **No positive gap follows.**

**PROVED — analytical:** the magnetic-information differential and finite
resolvent bound are derived for complete same-history probe words. For stacked
center rows bounded by b, stacked row difference bounded by d, and relative
innovation radius e<1, the information perturbation is at most
`(e*b*b + 2*b*d + d*d)/(1-e)`. These are conditional domain bounds, not point
norms. No finite service value is substituted for a uniform center floor.

The causal center is defined using the literal real-operation maps, current
center state, same source, inherited clocks and its own branch decisions.
It consumes no observed future H/G/K. On a regular planar branch the even
mean/P block drives the odd P block, which does not return through the nonzero
mean correction. Both covariance blocks and all odd error/gauge directions
remain. The one-way tuner theorem is not extended to delete reference/gate
feedback. The center's compact forward domain is **OPEN**.

**FINITE DIAGNOSTIC ONLY, retained unchanged:** the existing native stream and
untapped control agree, and 7800 finite service roots have minimum 7.0273018107.
The conservative carried floor 7.0247646056 leaves a finite Weyl budget
6.0247646056. The homogeneous factor word has Euclidean norm 122.6598857,
metric norm .93478409, physical quotient norm .87311969 and C_Q norm .0002290123.
It is not the full nonlinear derivative. Native secants are step-sensitive.
The observed-coefficient center defect .0006179245 is finite and uses future
observed H/G. None of these numbers is a source-uniform theorem constant.

## Current limiter

CR45 now pays an explicit part of each nonzero receipt from the already
retained full process. Its covariance-column energy and the remaining linked
receipt charge have no uniform sign certificate. Keep the original CR41--CR43
remainder; alpha is not a common reserve and supplies no full Schur sign.


The one remaining comparison is CR41--CR43's full feasible R_W(c)>=0
for one common c>0. Its receipt corner is exactly CR39/CR30 using
Hbar,Dbar,L_pmu and the actual inherited image/guards; the remaining
linked full process/S/magnetic and acc utilde action has not decided it.
No particular threshold-one mode is asserted feasible or excluded.
The full endogenous cost is retained in the same remainder and becomes
T_c-B_c'*S_c^-1*B_c after receipt positivity. Its common-reserve inverse
and sign are unproved; no compact forward domain or image closure is
assumed. Source completion cannot start before this comparison closes.

The first decisive test is now explicitly (CR3) for the ACTUAL correlated
homogeneous lift, with true exogenous variation set to zero and inherited
endogenous columns retained. No uniform positive remainder has been proved.
On the qualified odd covariance restriction the actual common-tangent
matrix (CR13a) is assembled; CR20--CR21 now pay a directional fraction
of its active regression square with the actual same-cycle acc correction.
CR25--CR27 also insert the preceding actual process and cover nonzero
scalar zero-gap branch cones. Their finite joint receipt weight and
retained square leave CR28's signed reader comparison unpaid; CR22 and
CR17a remain equivalent accounts of that correlated work.
The actual CR24 integrated process and S/magnetic covariance Grams have
not been proved to make the complete remainder nonnegative uniformly.
Neither the inherited image/rank nor the full cross remainder
(CR17) has been certified. No independent receipt bound repairs this gap.
Actual target/tuner and inherited auxiliary origin/graph qualifications,
other mean/covariance cross entries, source/reference and endpoint gauge
remain. A released-coordinate direction is not a shipping obstruction.
No finite-error, capture, arithmetic, radius or regime promotion follows.

On the qualified locked-Live lift, do not retain a fictitious reference-to-MEKF
feedback loop: (LR1) removes its reverse edge exactly. Source/reference
variations are still the linked (LR6)--(LR7) ports, including anchor history.
The unpaid target is (LR8)'s actual MEKF quadratic gap and grouped source/gauge
work; the new finite reference amplitude bound alone cannot make it positive.

**PROVED — analytical fixed-input even prediction reduction:**
`app:planar-even-process-work` proves `dF_E=dQ_E=0` for MEKF-root variations
on the exact planar stratum, using the literal rotation/integral axis kernel
and one-way private tuning. It does NOT prove `ds_E=0`: the normalized shipping
quaternion polynomial has angle derivative defect at most `x^6/7680` on
`|x|<.01`. Its remaining even auxiliary process charge is bounded by
`1e-27*(d BG_y)^2` on the stated real-operation profile. The process score
feedback and its signed cross with that residual remain grouped; the constant
is not a complete-word decrement. Odd blocks, varying source/private histories,
projection, reset, reference and AW work are not eliminated. In moving AW
coordinates the coefficient derivatives remain connection differences until
complete-word composition. This removes a fictitious independent even-process
parameter charge, not the requested inherited E_L or full active-work bound.

**PROVED — analytical uniform held-score reduction:**
`app:held-score-boundary` telescopes the unsupported physical held-BA score
to `P_b^-1 Delta b_s`, with energy at most
`j_b min(D_s^2 T^2,4 B_s^2)`. Its active/AW component and generated active
suffix contributions vanish on the reached zero-cross stratum. Residual,
`dr`, and `dR_eff=dR_acc+dP_b` ports remain. Held tangent cross/square work
also telescopes to one exact boundary term; no active contraction of held
coordinates is asserted and release/A21 are not covered.

**Precise remaining budget domain:** the same-operation identity
`B_W=E_L+S_q-2 C_W`, `C_W^2<=E_L S_q`, shows that with uniformly bounded
source-only score action, uniform B_W control is equivalent to control of
the actual process-loss comparison energy E_L. No general activated cap on
E_L is proved. This is less than a full state box, but is not supplied by
positive base loss, physical source amplitudes or first-Live energy alone.
The remaining unsupported AG/chart scores and grouped signed active
auxiliary/cross work against the exact Schur loss are still OPEN.

**PROVED — analytical combined process/source score:**
`app:process-source-score` gives `score' L^dagger score=V-V+ +s' Q^dagger s`
with exact support checks. Thus the covered complete score satisfies
`q_c q_c' <= B_W G`, where B_W telescopes root/terminal comparison energy,
actual physical process action, correction defect action minus NIS, and
signed uncovered energy changes. Unsupported held/chart/projection scores
remain explicitly in q_u; first-Live energy is not restarted. This replaces
separate loss/source-score bounds on the covered process substeps.

**PROVED — analytical causal generated-process square control:** an exact
algebraic factorization puts dP,dF,dQ,ds into one process Fisher loss and
proves `|f|^2 <= (chi/2)*ell_P`. Its complete storage identity retains signed
auxiliary work, mean loss, the feedback cross term and square. Derivatives
include actual BG, jointly generated tuning and auxiliary/source recurrence;
consistent AW-frame derivatives retain all endpoint connections. It is not
an independent source model or an extra decrement. Do not separately maximize
the auxiliary charges and lose their cancellation with Fisher/cross work.
The [exact formulas and scope](ou3-linked-word-score.md#combined-physical-process-score-and-causal-generated-work)
are controlling. Uniform B_W/uncovered-score bounds and affordable signed
work absorption against the exact Schur loss remain OPEN. No radius is solved.


**E — rendering corrected:** the new CP7 equation tag was initially inside
`split`; moved it outside. No formula, premise or validation gate changed.
The exact-rational checks and source-bound evidence remain controlling.

**PROVED — analytical complementary base-packet absorption:**
`app:conditional-packet-absorption` proves `P_self,b+P_mix,b <= W_complement`
from the exact full joint base Schur loss. It retains all covariance cross
entries and consistent AW coordinates. This removes independent pricing of
base complementary-root propagation; it does not absorb the actual score
or generated coefficient/source ports. The remaining signed work is exactly
`2<T_base v,R>+|R|^2`, with all suffix scores and endpoint gauge retained.

**D — sufficient-bound failure and formulation review:** the fixed half-charge
criterion is negative (`-1173/21200`) on an exact formal strictly contractive
base specialization whose actual loss is `159/400`; the exact base Schur
remainder is `3627/21200>0`. This is not a reached shipping word. Therefore
half-charge domination is optional, not a necessary milestone. Separately,
the linked anisotropic score `q=(I-T)e` has the exact directional threshold
`a^2(1-epsilon)/(4 epsilon(1+epsilon))` in a weak direction, despite fixed
loss-score energy `a^2/2`. The isotropic cancellation cannot be extrapolated.
No shipping counterexample, failed replay, or interval instability is inferred.
Retain the actual directional reader and complete signed matrix; do not launch
a scalar tightening campaign. Uniform score and remaining signed-work margins
are OPEN. See [the proof and activation table](ou3-linked-word-score.md#complementary-base-packets-are-already-paid-for).


**PROVED — analytical actual conditional mixed-block completion:**
`app:conditional-mixed-word` now evaluates the root mean/covariance score
coupling in the actual conditional metric. Its exact three-row test is
`lambda K^-1-E_q G_F^-1 E_q'>0`, with C,T,q from the same word. It does not
use the unknown full gap as an input or replace the score by an independent
maximum. All internally generated covariance terms, including U*q_suffix,
remain in aggregate terminal packets t,Z. Exact completion gives the bound
`gap >= (1/2) z' D z+W_complement-P_self-2 P_mix+endpoint_adjustment`.
The packet terms are endogenous functions of the same inherited variation,
not independent disturbances. The state-dependent AW frame preserves C,dC,
u,Y,T,q_a and this mixed reader, even with dL retained. Its endpoint work
is exactly in the transformed complementary storage and packet terms (CM12).
Use those consistently without adding a second independent frame penalty;
input-frame accounting still needs Phi_0-Phi_N. Gauge work remains signed.

**OPEN:** a uniform positive actual score margin for one lambda, and linked
packet absorption against the remaining process/S/magnetic action. The
proved c_aw makes the base inverses legitimate; it does not establish this
net margin. The exact isotropic relative-loss identity cancels the small
loss from the score threshold, but is not imposed on general reached
conditional covariances. No radius or numerical search is performed.
See [the explicit mixed-work bound](ou3-linked-word-score.md#actual-conditional-mixed-block-bound).


**PROVED — analytical conditional AW process loss:**
`app:aw-conditional-loss` splits joint storage into conditional AW and
complementary terms with all covariance cross terms retained. Optimal additive
mag/S substeps leave the conditional term invariant, even with non-AW row/noise
ports. Their endpoint-square coefficient decreases through the actual attitude
covariance loss. The three mean/six covariance directions missed by those rows
are not physical gauges: the full integrated OU noise gives a positive uniform
base-action loss `c_aw W_cond`, with exact rational `c_aw>10^-14` on the existing
qualified isotropic regular profile. This is not complete tangent contraction;
process-generated covariance scores and all signed mixed work remain.
The existing 16-second LIN path matrix supplies conditional precision
`J_aw,aw<15942618 I` on its qualified A21 domain. Finiteness there is no longer
the blocker. The exact process recurrence (CA6)--(CA7) retains signed
integrated-chain and attitude/BG covariance work.

**OPEN — actual endpoint Schur margin:** no positive net complementary block
or full-word c has been proved. The precise missing comparison is between
these same-history signed mixed blocks (including covariance suffix scores,
reference/source and moving gauge) and the actual process/S/magnetic action.
A positive base conditional block cannot replace the net signed block in a
Schur complement. This is D, insufficient proved domination; no violating
shipping history or numerical failure is inferred. The new block coefficient
is not a reason to launch a scalar refinement campaign. See the
[minimal activation table](ou3-measurement-frame-loss.md#conditional-aw-action-and-the-endpoint-schur-blocker).


**PROVED — analytical complete-word frame reduction:**
`app:aw-frame-word` composes the AW frame jump with all subsequent literal
operations and generated covariance suffix scores. Internal frame connections
cancel in the full action-minus-work balance, leaving `Phi_N-Phi_0` only.
The endpoint work is exactly `Phi=2 h' daw+daw' Q daw`, with `Q>0` and
`Phi=(daw+Q^-1 h)'Q(daw+Q^-1 h)-h'Q^-1 h >= -W`. This is a signed identity
and lower bound, not a free deletion or upper bound on terminal work.
[The exact three-row Schur test](ou3-measurement-frame-loss.md#complete-word-frame-cancellation-and-directional-endpoint-absorption)
retains the actual endpoint tangent maps and moving gauge. Its threshold is
one, with the complementary block required positive. No margin is certified.
The per-jump accumulation obligation is removed. In addition, the actual
PSD AW synchronization adds its increment exactly to the conditional AW
covariance. Conditional precision, conditional comparison energy and the
Fisher connection-square coefficient cannot increase at that operation.
The same-history endpoint-square bound is `(V_cond+2 lambda)*K(daw)`;
this monotone coefficient is not a second decrement or a bound on signed h. Remaining activation for this
term uses actual **conditional AW precision** `(P^-1)_aw,aw`, full attitude
covariance, and linked comparison/tangent rows. Neither the marginal AW floor
nor the pitch-only ceiling proves those bounds. The qualified LIN-path
conditional precision bound above supplies only its stated component. Uniform absorption together
with physical prediction/reset/reference/source work remains **OPEN**; no
radius, new gauge dimension or independent covariance box is introduced.


**PROVED — analytical AW shear reduction:** the full-state invertible frame
`L=I-E_aw[aw]x E_theta^T` makes the nominal-world accelerometer row constant.
The remaining AW-row covariance port is exactly
`A(Gamma P+P Gamma')A'-(Gamma C+C Gamma')`, so it cancels in the covariant
correction tangent. The fixed-noise additive substep retains the full Fisher
covariance loss, with exact joint balance
`Delta W=-||H eta-m||_S^-1^2-lambda L_P+||m||_R^-1^2`.
The mismatch m is linked and endogenous, not independent forcing.
[The derivation](ou3-measurement-frame-loss.md#exact-aw-shear-removes-the-remaining-nominal-accelerometer-row-derivative)
retains the actual post-correction jump
`B_L=I+N(K_aw r)`, `dB_L=N(dK_aw r+K_aw dr)` and its downstream covariance
score. Literal OU prediction retains its attitude-frame discrepancy, BG and
integrated LIN columns. The actual AW increment is unchanged, but its face
uses the recovered original AW marginal. Nothing is independently boxed.
The physical gauge is transformed, not deleted.

For the central planar physical comparison, the sensor equation cancels the
nominal BA term in the mean mismatch. Its remaining curvature obeys
`||m0||_R^-1^2 <= V*tr(R^-1 J_y B_v P B_v' J_y')`, where B_v uses the
actual physical-force divided rotation and nominal variation. The exact
coefficient retains signed pitch/AW cross covariance; no nominal AW/BA box
is required. Source cross-plus-square work remains linked. Uniform activation
of V and the joint covariance block, and word-level absorption, are OPEN.

The qualified magnetic 9/10 loss and its pitch/BG ceilings survive in this
same joint storage: parity, attitude/BG marginal and physical comparison
energy are unchanged. **OPEN:** absorb the combined mean mismatch, actual
frame-jump, prediction/reset and noise/reference/AW-face work against the
complete transverse action. Do not separately charge the now-cancelled AW-row
norm or add losses belonging to different joint storage coordinates.

**PROVED — analytical linked suffix-score normal form:**
[The word-score proof and minimal activation table](ou3-linked-word-score.md)
regroup the exact complete tangent before norms. Every optimal additive
correction has zero local information score for arbitrary residual. The word
retains all generated covariance ports, including the essential `U_i q_(i+1,N)`
term. For a score in the base-action range, `Gz=q` gives the exact linked bound
`||M dP q||_J_N^2 <= (q'z/2)*L_P`. Its process/AW loss-generated part satisfies
`q'z <= sum e_i'Q_i e_i` on the same history. Physical/source/chart score and
any component outside the action range remain forcing. No NIS cap or new
nominal p/v box is required for these identities.

**PROVED — analytical inherited comparison budget:** telescoping the literal
optimal correction and prediction/AW/reset/projection comparison energy gives
`chi_loss <= E_loss = V0-VN+Supply_W-sum NIS`. The supply retains signed
source/chart work and actual correction defect energy. At `S=0` its defect is
`-S_physical`, not zero. Thus the loss-generated score needs root comparison
energy plus linked source supply, not separate all-time bounds on intermediate
nominal coordinates. Realized NIS is not a transverse observability bound.
Uniform supply, source-score kernel forcing and generated-port absorption
remain OPEN. First-Live energy cannot be reseeded at later roots.

**CONDITIONAL root-block threshold:** if a justified base action has margin
`c>0` and combined score charge at most `chi_*`, then
`chi_* < 2*c*lambda` makes the root mean/covariance-loss block positive.
This does not cover the remaining generated ports or prove covariance
coercivity. Uniform `c`, the linked comparison-loss/source-score budget and
signed absorption of the complete `U,v` work remain **OPEN**. Apply this
regrouping to the existing word identity, not as a second charge for magnetic
feedback already absorbed below. The physical fibre is not equated with an
action kernel. An out-of-range score cannot be discarded by a pseudoinverse.

**New analytical progress, qualified uniform magnetic loss:**
[Moving-frame loss and its activation table](ou3-measurement-frame-loss.md)
derives the exact transformed measurement rows, retaining their moving-frame
connection. In the planar nominal stratum the coupled magnetic correction
loses at least `9/10*(L_eta+L_P)` when the qualified reference cone, retained
six-degree chart and comparison energy `E<=100` hold. The rational loss-matrix
margin is `4673/1250000>0`, including the covariance-induced square charge.
This is an operation loss fraction, **not** a complete-word contraction rate.

A prior-observation covariance reader proves `P_pitch,pitch<3/5000` and
`P_BG_y,BG_y<1/4000` after `753/250` seconds of regular planar applied service,
without an initial covariance ceiling. The windows precede the correction
being analyzed. Together with the qualified reference cone this gives `t<6`;
the retained chart gives `|kappa|<3/500`. Hence the covered coupled correction
loses at least `1125*eta_pitch^2+1260*(dP J dP)_pitch,pitch`. These are derived
uniform row constants, not finite observations. The wrapper's `b0=1e-10f`
comes through `initialize_ext`; do not use the different standalone default.
The S substep retains its exact coupled loss on the one-way fixed-noise path.

**OPEN activation:** the reader assumes existing service and the qualified
reference cone, and does not prove the planar witness admitted. Chart/energy
retention, float reference acquisition and arithmetic transfer remain open.
The nominal planar stratum is not the physical roll/BA compatibility fibre.
The improved word action absorbs the covered magnetic/S ports; accelerometer
AW-row feedback, prediction/reset/frame/AW/projection work and source/gauge
terms remain linked. The precise missing gap is the improved action minus
that remaining signed work and endpoint gauge transport, uniformly positive
on a justified transverse root lift. No unrelated p/v box is required.

**PROVED — analytical word balance and kernel characterization:**
`information_shear_word.py` extends the coupled shear through literal
reset discrepancy, default AW increment, bias projection and moving-projector
transport. With actual causal prefixes, the nonnegative action has kernel
`intersection ker(Q_i T_i)`. The exact word balance retains signed port work,
its square charge and endpoint gauge-energy difference. No endogenous port is
reclassified as independent noise. The physical same-record tangent is a
forced comparison, not automatically the homogeneous or action kernel.

**OPEN:** for every activated admitted word on a justified compatible stratum,
prove the linked transverse matrix gap in the strategy note with one common
`c>0`. Existing magnetic service controls its physical probes; a quantitative
comparison to the actual coupled shear action, process/S action and absorption
of signed work/projector terms is missing. Pointwise rank is insufficient.
The general compatible-set/stratum relation must also be justified, rather
than extrapolating the planar fibre. These are proof dependencies, not new
physical assumptions.

The strategy note's activation table limits domain work to the quantities used
by that comparison: required metric/innovation ratios, comparison energy and
row geometry, reference/gates, reset chart, AW/clock branch coverage and source/
gauge/arithmetic charges. The first-Live central-fibre bound is only a boundary
lemma. No all-time compactness is assumed to prove its own retention. Uniform
`c`, actual charges and every-prefix bounds must precede radius algebra.

The planar profile is sigma_a=.2 with adaptive S; AtomS3R uses sigma_a=.12 at
its reference rate with fixed S and differing settings. No profile transfer.
Physical S is never reset by the estimator pseudo-measurement.

## Failed approaches / DEAD_ENDS

* **D, receipt-corner decision using qualitative linked action:** the
  requested sign of CR41's `S_c` was NOT established. The integrated
  `g,Q_ao` recurrence and the actual transported `D h_S'`, `D h_mag'`
  and `D h_acc'` were checked on the same CR31--CR32 tangent. At the
  possible threshold mode, CR40 still gives
  `Hbar*y=L_pmu'*Dbar*mu`, `mu=L_pmu*y`, and
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

The shipping sketch's BMI270 accelerometer runs at the gyroscope's 200 Hz
(`configureAtomS3RAccel200Hz()`: `ACC_CONF` 0xA9, same filter). M5Unified had
left it at the 100 Hz reset rate and the sketch accepts a frame only with both
data-ready bits, so deployed IMU intervals were 10 ms, outside the 0.006 s
admitted by [sampling fidelity](ou3-sampling-fidelity.md). Register polling
against the sensor clock and loop overruns can still give an occasional longer
interval, so the interval bound stays a deployment qualification. The sketch
also stops before the estimator when the saved calibration cannot be applied
in full. Estimator source, constants, gates and physical assumptions are
unchanged; the sketch's provenance blob and the native-secant `src/` hash are
rebound.

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
domain checks. The heading correction adds a sensor-derived steering-force
model and gyro-relative affine magnetic calibration to the shipping path.
Native observer/control, carried-state, mean-action and planar secant records
are regenerated; finite records do not establish a source-uniform extension.
The world-frame geometry scope explicitly excludes nonzero modeled steering
speed, while the generic operation identities retain complete literal rows.
Current test/CI results and command receipts are recorded with the source
refresh and in the PR metadata.

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
