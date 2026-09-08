# OU-III proof research state

## Current hypothesis

BRMM is the primary physical motion source. The target is bounded
accelerometer-bias error plus practical ISS of the other 18 errors, using
the actual full 21-state estimator. Preserve the common physical history,
every actual anisotropic R_S update, full Q/floors/cross covariance, finite
reset/projection, physical -S_true forcing and the closed .4 bias estimate
ball. Bias-error decay and a zero practical floor are not requirements.
BRMM-0 requires a_m=dv_m/dt with uniformly bounded same-history velocity;
constant physical acceleration DC is excluded in Q, O and mixed histories.
A constant sensor offset remains BIAS0. Bounded velocity alone does not
bound displacement or S. Spectral models are reference specializations.

## Certified subcertificate

Main 85e040371812d1fa32508d521acded0a6328102f has rebuilt conditional BRMM
P3 in run 34173418232, job 101899163631. All 56 P3 tests passed. Both H18
and A21 have delta=1e-18; the H18 worst outward LDLT pivot is
4.987499868870966e-14 and the first active A21 bias margin is
1.2499987189052501e-9. Artifact 10036647100 has verified ZIP SHA256
03d4011e49f0c54224cc129709631a26d5ee8b5399a435e9e083f22026b8ab03.

This proves the matrix implication under the explicit configured
Normal-Live premises, including accepted-vector geometry/recurrence.
It does not establish physical admission, projected nonlinear contraction,
P4 or P5. The numerical configuration remains acceleration 4 m/s^2, rate
30 deg/s, Racc=.04 I and Rmag=.09 I. A broader proposed domain cannot
inherit this PASS. BRMM recurrence alone does not imply vector PE.

## Runtime evidence

The runtime audit executes all eight unchanged v1.1.3 reference records with
the ordinary simulator, noise/seeds, 25 Hz magnetometer and quality gates.
An insertion-only overlay records outer/inner Live before and after each
sample, held/active bias, actual guard action, accepted Jacobian vectors,
reference writes, resets and every S update with its actual full 3x3 R_S.
Runtime flags define the phase partition independently of physical caps.
No source sample is removed to manufacture admission.

All eight dual runs passed byte-for-byte ordinary output comparison with
identical simulator compiler flags including -ffp-contract=off. They retain
1,920,000 samples and 315,721 actual S corrections. Sample summaries and
source/output/trace identities are in docs/ou3-brmm-runtime-audit.json;
the workflow retains complete compressed traces and violating intervals.
The archive remains the same input used by the prior source audit.

| Reference | H_s m | Runtime Live samples | Acceleration >4 | Body rate >30 deg/s |
| --- | ---: | ---: | ---: | ---: |
| JONSWAP | .27 | 233632 | 0 | 249 |
| JONSWAP | 1.5 | 233632 | 1909 | 14401 |
| JONSWAP | 4 | 233632 | 11442 | 30076 |
| JONSWAP | 8.5 | 233632 | 29832 | 47120 |
| PM-Stokes | .27 | 233632 | 0 | 1670 |
| PM-Stokes | 1.5 | 233632 | 603 | 26812 |
| PM-Stokes | 4 | 233632 | 9622 | 49405 |
| PM-Stokes | 8.5 | 233502 | 32476 | 71050 |

These violations occur during runtime Live with the guard dormant, no
hard tilt reset and no rejected accelerometer update. They are not explained
by startup alone. Six records also violate the candidate 2 m/s impulse cap
on complete 10-second runtime Live windows. That candidate is still unfrozen.
The prior 1,904,008-window source lobe audit found only O windows, no Q
coverage. Sampled lobe/primitive statistics remain non-promoting.

Actual replay Racc=.000866618473 I and Rmag=3.68640018 I differ from the
configured P3 matrices. Therefore no complete reference history is admitted
by that exact configuration, independently of physical cap violations.
Every sample-aligned 1-second runtime Live window has an accepted mag event;
this is packet recurrence, not a transported PE certificate. Some same-sample
vector sines drop below .1; that does not disprove asynchronous window PE.
Continuous magnetic-reference writes are retained variation, not automatically
forbidden hard regauges.

## Separate phase envelopes and physical qualification

The requested phase separation gives sampled maxima of approximately
7.083 m/s^2 and 178.480 deg/s before Live, versus 16.217 m/s^2 and
198.355 deg/s in Live. A rounded envelope with at least 25 percent headroom
would be 9 m/s^2 and 225 deg/s before Live, and 21 m/s^2 and 250 deg/s in
Live. These are candidate qualification targets, not certified physical
bounds or permission to transfer the old P3 result.

Only the H_s=8.5 JONSWAP record exceeds g in non-gravitational acceleration:
845 samples, all during Live. At sample 54116 (270.58 s) its vector is
approximately (6.52644,-3.42066,14.4454) m/s^2. The release generator
v1.1.3/c5ddd8dd reproduces that point: linear vertical contribution 6.67798,
simplified second-order contribution 7.76747 m/s^2. Linear vertical motion
alone ranges from -5.91312 to 6.70258 over the record. The second-order
correction exceeding the linear term is a model-validity concern, not proof
of physical vessel acceleration. The generator follows surface slopes and
uses 128 components over .02--.8 Hz plus sum-frequency terms up to 1.6 Hz;
it is not a qualified vessel response model. Do not treat its sampled
extrema as a physical validation of wider BRMM constants.
The primary implementation is the release's
[Jonswap3dStokesWaves.h](https://github.com/bareboat-necessities/oceanography-waves-lib/blob/v1.1.3/src/Jonswap3dStokesWaves.h)
with parameters from
[waves_sim.cpp](https://github.com/bareboat-necessities/oceanography-waves-lib/blob/v1.1.3/data-sim/waves_sim.cpp).
The point decomposition uses that header's amplitudes/frequencies, its
mt19937 seed-42 paired phases, T_p=float(11.4), and t=54116*float(1/200).
The linear vertical term is -sum A_i*(2*pi*f_i)^2*sin(phi_i-2*pi*f_i*t);
subtract it from the header's returned vertical acceleration for the bound term.

## Failure analysis and current limiter

The hypothesis that startup exclusions reconcile the old caps with all
reference Live histories fails. Class: source/configuration admission
failure, not filter instability, numerical enclosure failure or a P3
counterexample. It invalidates coverage claims for these whole replays by
the existing deployment configuration. Narrowing sample boxes cannot fix it.

A candidate Live cap of 21 m/s^2 also removes the triangle-inequality force
floor g-4. Any wider theorem requires an independent nonvanishing specific
force/geometry premise or a different window observability argument.
The candidate rate 250 deg/s was tested through the existing windowed PE
producer (with positive candidate force floor .5). It fails closed with
"declared PE recurrence/rate box is too wide for the current two-occurrence
transport bound". At the unchanged one-second recurrence, the controlling
bracket 1 - omega_max*(3 s)/2 is -5.544984694978735. Class: proof-method
transport failure; it does not falsify observability using all actual
accepted vectors. No gain from tightening interval rounding can cross zero.

Critic pass: the two-occurrence scalar transport bound cannot support the
requested high-rate envelope. Alternatives are (1) use the full transported
accepted-vector information history, retaining all intermediate acc events;
(2) qualify an explicitly different recurrence premise from deployment,
without silently substituting hardware ODR for accepted-vector PE; (3)
qualify a vessel-response source and separate surface-follower stress tests.
The physical validity of the extreme second-order wave motion must be
resolved before choosing the enlarged deployment domain.

## Retained facts and dead ends

- The old independent-port common-gain P4 witness is quantitatively unusable.
  Point factors .9999569976489486/.9788191291615017 with gains 2^32/2^36
  imply storage bounds at least 1.4037e15/2.2313e13. These are neither
  nonlinear contraction factors nor actual filter error floors. Its discarded
  source correlations cannot be recovered by interval refinement.
- Retain the common primitive graph coupling latent increments, physical
  -S_true, bias and every actual R_S. Useful channel gains and every-prefix
  retention must precede source-uniform covering and P4/P5 promotion.
- A qualified true-bias bound plus the .4 estimate ball gives bounded bias
  error. BIAS0/1 need qualification; optional BIAS2 needs its actual nonlinear
  sector. Configured tau_b=5000 s is not sensor qualification.
- Native FMA overlay output differed at sample 6371, initially in gyro-bias
  last digits. Class: instrumentation/code-generation failure. The single
  refinement used identical explicit -ffp-contract=off in baseline and
  overlay and passed all eight output comparisons. No tolerance was relaxed.

## Next falsifiable experiment

Qualify the phase-specific physical source envelope and actual measurement
configuration, including the large second-order reference peak. Then test
full transported vector information on that declared domain before a wider
P3 rebuild at delta=1e-18. Keep the current narrow conditional certificate
identified separately. P4 next consumes same-history primitive equalities
and explicit channel budgets; it requires a useful practical bound and
prefix retention before any rigorous source covering. P5 remains blocked.

## Validation state

PR #505 revision 3d487972 passed the eight-case runtime CI job 101905861666,
conditional P3 job 101906259900 (56 tests), and quality gates. All runtime
counts agree with the local result. The eight attachment tests and eight raw
source tests pass. No nonlinear P4/P5 promotion follows from green tests.

Full validation/robustness replays completed on main 85e04037 in run
34173418479. Publication job 101903490256 timed out after the source-only
spectral-moment tests, before the next BRMM search completed. Class: CI routing
failure; the stage-prefix exclusion missed renamed BRMM modules. It invalidates
publication completion, not the regenerated raw rows. Increasing timeouts
again would leave the same coupling. Alternatives: route source-only tests to
their existing proof jobs, shard them, or reuse exact-revision proof artifacts.
The selected correction routes them to explicit CI owners, with a regression
checking ownership; the full test target and evidence contracts remain intact.

The two full artifacts were recovered, their ZIP identities verified, and
provenance initialized by the normal contract at the genuine replay commit.
The contract passes against the current replay dependency closure. Original
raw-row identity and replay commit are retained. Artifact identities and
release checks are in docs/ou3-brmm-main-handover.md. Main must still complete
its automatic full publication and conservative fingerprint check before a
release tag; do not restamp a fingerprint to avoid that gate.

The user owns the generator correction. Continue physical qualification only
after corrected reference records are supplied. Until then preserve the v1.1.3
audit as evidence of the existing records and leave all proposed caps unfrozen.
