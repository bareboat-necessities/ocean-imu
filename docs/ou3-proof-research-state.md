# OU-III proof research state

## Current hypothesis

BRMM is the primary physical source: bounded accelerometer-bias error and
practical ISS of the other 18 errors, using the full 21-state estimator.
Keep the common physical history, all actual anisotropic R_S corrections,
full Q/floors/cross covariance, finite reset/projection, physical -S_true,
and the closed .4 bias-estimate ball. BRMM-0 requires a_m=dv_m/dt and bounded
same-history velocity; fixed physical acceleration DC is excluded. Constant
sensor offset remains BIAS0. Bounded velocity alone does not bound p or S.

## Current dataset and evidence

All current replays use oceanography-waves-lib v1.2.1,
`sim-data-files-vessel-rao-28ft.zip`, SHA256
`6d6eb92db78e97f1d2456c6b92387993bf23be14a9a6cf26943f0a22fc6fee60`.
Its generator commit is e442150682f560384be427df4cc7815956a091c5.
It applies an estimated stationary 28 ft fin-keel sailboat RAO to incident
components and emits exact CG kinematics and body-rate derivatives. Incident
spectra and filenames retain their incident-sea meaning. This engineering
preset is not hull-specific measurement or a continuum BRMM certificate.

The source and selected-configuration runtime audits are current. Among 1,920,000
source samples, 10,205 exceed 4 m/s^2 and 31,937 exceed 30 deg/s.
All 1,904,008 sampled 10 s windows are diagnostic O windows. Runtime Live
contains 10,002 acceleration and 31,589 rate violations. No sample passes
all configured execution checks: actual vector covariances differ from P3.
Two of eight ordinary replays fail unchanged performance gates; the trace
observer is byte-identical to its matched baseline in all eight cases.
Current runtime summary: `docs/ou3-brmm-runtime-audit.json`; full interval audit:
`reports/results/rao_stability/runtime-audit.json`. Source audit run:
https://github.com/bareboat-necessities/ocean-imu/actions/runs/34186435547

The source-connected PM 1.5 m probe now uses this RAO and one common incident
root. Independent Python reconstruction agrees with the C++ source to
5.22e-15 (H18) and 1.42e-13 (A21). Both fixed 600-sample windows satisfy the
checked point caps/branch premises. The actual forced endpoint ratios are
7.2745272436 (H18) and .4360473133 (A21); the A21 prefix maximum is 1.2968761874.
These ratios include nonzero physical forcing and are not homogeneous P4
contraction tests or counterexamples. No source or direction search was used.

## Conditional P3 scope

The established matrix implication uses delta=1e-18 for H18 and A21 under
explicit configured Normal-Live premises: acceleration 4 m/s^2, body rate
30 deg/s, Racc=.04 I, Rmag=.09 I and accepted-vector geometry/recurrence.
The retained H18 outward LDLT pivot is 4.987499868870966e-14 and first active
A21 bias margin 1.2499987189052501e-9. The selected-configuration rerun retains the conditional result; canonical
outputs and hashes are in `reports/results/rao_stability/`. Conditional P3 passes for both H18 and A21;
deployment/source admission remains false.
Canonical run: https://github.com/bareboat-necessities/ocean-imu/actions/runs/34186435497
These are conditional matrix facts, not full physical admission, nonlinear
projection stability, P4, or P5. BRMM recurrence alone does not imply PE.
Actual replay Racc=.000866618473 I and Rmag=58.9824028 I remain different.
Performance retuning is now user-authorized; the canonical proof premises
remain fixed. The runtime audit verifies passive-trace parity on the selected configuration.

## Failure analysis and current limiter

The old claim that startup exclusions alone reconcile surface-reference
histories with the configured caps failed source admission, not stability.
The new RAO source audit still observes cap/impulse violations in larger
incident seas. Qualification must use complete runtime phases, actual
accepted vectors and actual measurement matrices before widening P3.

The forced H18 ratio above one invalidates interpreting forced endpoint
energy as homogeneous contraction. It does not invalidate ISS with a supply
term, P3, or the filter. The limiter remains a useful same-history practical
supply/storage bound with every-prefix retention and the actual nonlinear
projection. A favorable forced A21 endpoint cannot discharge that obligation.

Canonical P4 reran and remains false: the full augmented endpoint inequality
`-(L_W + sum lambda_j Pi_j) > 0`, BIAS1 projection compatibility, BIAS2
corrected-error separation, and every-prefix finite gain/domain retention
are unclosed. Classification: proof-method obligations remain unproved,
not a theorem counterexample. The rerun validates the retained algebra but
invalidates claiming a complete nonlinear certificate. P4-motion is also
false; P5 and P5-motion may not start. The next experiment must quantify a
useful same-history supply/storage bound before any interval refinement.

The required `make all` compiles and reaches the existing frequency gate:
`freq-track.cpp:321: Assertion all_quality_ok && "freq_track quality gate failed"`.
This is a performance regression under the changed source, not a proof failure.
NLO also fails its unchanged vertical gate (for example 12.7763% > 7.13%).
Collect-all execution retains every record and a failing final exit status.
No performance threshold is changed. Parameter selection is evaluated
separately on paired default and additional sensor/initialization draws;
proof margins are not an optimization objective. Small NLO theta gains
produced large drift errors (up to 325.515% Hs), invalidating that tuning
direction, not the stability theorem. Retain the baseline unless a candidate
improves the full record set and separate validation draws.

The frequency failures occur on weak cnoidal CG acceleration with injected
DC bias and noise. The zero-crossing threshold is larger than the wave
amplitude and reports its fallback. Forcing estimator outputs together or
changing the 70% agreement threshold is not a valid tuning remedy; a
separately justified bias-rejection change would be an algorithm change.

The prior travel reference `azimuth + 180` fails on all four JONSWAP
RAO records by 166--174 degrees. This is a source-convention implementation
defect, not a stability failure. v1.2.1 uses `cos(k.d - omega*t + phase)`;
its azimuth is propagation-to. Independent RAO-matched orbital correlation passes on all eight records;
all 11 convention/law-mirror tests pass, with the existing 20-degree
convention-test tolerance unchanged. No performance threshold is relaxed.

## Retained facts / DEAD_ENDS

- Independent-port common-gain P4 storage was quantitatively unusable:
  point factors .9999569976489486/.9788191291615017 and gains 2^32/2^36 gave
  bounds at least 1.4037e15/2.2313e13. Do not refine that discarded coupling.
- Two-occurrence scalar PE transport cannot cover the proposed 250 deg/s,
  one-second recurrence: 1-omega_max*(3 s)/2=-5.544984695. Interval tightening
  cannot repair that sign. Those surface-derived cap proposals are unfrozen.
- Bias qualification plus the .4 estimate ball bounds bias error; BIAS0/1
  need qualification and optional BIAS2 needs its actual nonlinear sector.
- Passive runtime instrumentation must compare byte-for-byte to its baseline;
  matched -ffp-contract=off resolved the prior native FMA instrumentation defect.

## Critic pass, alternatives and next experiment

The strongest reason to abandon the old architecture is its loss of physical
source correlation and uninformative common gain. Retain full matrix/group
structure and compare: (1) common latent-increment/p/S primitives with explicit
channel budgets; (2) full transported accepted-vector information instead of
two-occurrence scalar transport; (3) a separately qualified recurrence/response
premise from deployment. None may be fitted to replay extrema for PASS.

Physical admission requires the actual covariance and vector-history premises,
not merely the runtime Live flag. A new P4 method needs useful practical
bounds and every-prefix retention before rigorous source covering. P5
remains blocked by the unclosed P4 obligation.
