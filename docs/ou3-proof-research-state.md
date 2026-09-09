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

The current limiting quantity is no longer the separated bias budget. The
declared-domain coordinate retention experiment
(`reports/results/rao_stability/domain-retention.md`) propagates the same
attached word in the declared physical coordinates instead of a scalar
storage. With the closed .4 m/s^2 bias ball and the same forcing template as
the only initial deviation, every motion coordinate is retained: 30-degree
chart .0703 (H18) and .0041 (A21) of its radius, velocity .1713 / .2382,
position .0444 / .0844, and latent acceleration .1431 / .1397. Route 1's
sufficient storages 6.4111 / 2.2693 are quadratic, so their linear equivalents
are 2.532 / 1.506 chart radii, a factor 36.0 / 364 above these. On the rows
that carry this conclusion the attained bound is within 1.3% of the certified
one, so it is a tight comparison, not a looser enclosure.

Classification: route 1's chart failure is a proof-method (scalarization)
failure. Charging a bias-driven velocity and displacement excursion to
attitude through the worst direction of the information metric manufactures
it. This invalidates treating the separated bias budget as the limiter and
invalidates any further refinement of that budget. It does not prove P4,
admit a physical source, change any endpoint or prefix ratio, or establish
anything uniform over the nonlinear source family.

The accelerometer-bias row is the one exception and a separate open premise:
H18 maps the closed .4 m/s^2 ball into itself exactly, but A21 reaches
.400126 m/s^2, a relative 3.146e-4 above the declared radius and far outside
roundoff. The ball is therefore invariant for one mode only, and a
bounded-bias statement over repeated words has to carry that growth. No motion
row depends on it.

The new limiter is the declared handoff set itself. Under the complete
declared product box the same word leaves the domain, limited by velocity at
35.358 (H18) and 9.357 (A21); the largest box of that shape it retains is
2.83% / 10.69% of the declared one. The dominant single source is the declared
300 m*s integral-displacement ball, which alone drives H18 velocity to
159.7 m/s against 5 m/s and position to 379.3 m against 20 m. That ball is not
independently reachable: the integral state is the running integral of the
position state, so 300 m*s with position at most 20 m needs a sustained 20 m
error for 15 s. Deleting it still leaves velocity at 4.033 / 4.921, with
position and attitude as the next sources, so a correlated integral/position
fact is necessary but not sufficient.

The 17 conditional P3 producer/promotion tests pass.

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

The frequency, TFG, PII and selected OU-II default records pass. OU-II's
vessel profile retains a 5 s post-statistical stillness decay and the paired
bias/drift retune in `reports/results/sigma_horizon/study.md`: all eight
defaults pass, and fresh failing records decrease from 28 to 24, with a 4.6%
mean pitch tradeoff. The moment averaging remains four periods / 35 s;
increasing it alone failed to improve low-wave accuracy. OU-III's two
high-sea Z-bias gates and NLO's low-wave PM-Stokes vertical gate remain
unresolved. `reports/results/rao_parameter_tuning/ou3_acc_z_bias/study.md`
screens the coefficients that move the OU-III metric: only the accelerometer
bias prior, its correlation time and its per-axis driving noise do, and the
one screened profile that clears all eight default-seed gates changes the
scored quantity by -0.4% over eight fresh paired seeds while costing 1.1% of
low-wave vertical accuracy. Over those seeds the committed profile itself
scores between 2.51% and 14.58% on the two failing records, so the 4.3% bar,
cut from one realization, is violated on 22 of 64 fresh cells before and after
every screened profile. No performance gate, deployed coefficient or canonical
proof premise is changed.

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

The point experiment freezes the attached finite coefficients and keeps
one common scalar multiplier on the entire physical forcing template,
including the latent increments and physical -S_true at every due update.
It tests the complete 21-state root/template quadratic form, including the
corrected bias trajectory, against the 18-error endpoint and every prefix.
This is a necessary feasibility screen for retaining source correlation,
not a certificate for varying physical roots or nonlinear coefficients.
This reduces endpoint gains to 65536 (H18) and 1048576 (A21), but gives
candidate composed storage bounds 1.6884e9 and 8.8012e7. H18's every-prefix
test at Gamma=2 fails at the first accelerometer correction with ratio
2.70996. This is a proof-method/supply failure, not a legal nonlinear
counterexample: the root direction and frozen coefficients are relaxed.
The common-template/common-gain restriction is also a DEAD_END pending a
new controlling fact. It does not invalidate correlated-source ISS itself.

BIAS2 is invoked in the dense endpoint and prefix master matrices, using
actual accelerometer covariance and the corrected 21-state history through
finite reset/projection. For A21, point energies are g=272.62755,
b=234.48615, cross=-240.71555 and y=25.68261; kappa=.9493553.
For Xi=sum b_error^T Racc^-1 b_error, the point sufficient mu is .1095272.
These are point values, not admitted source-uniform constants. H18 has zero
captured bias error and cannot identify a uniform mu from that point.

The explicit test is L + lambda*(Q_y - mu*Q_b) < 0, lambda>=0,
with each prefix using only its executed accelerometer events. The fixed
conditional mu/gain/multiplier grid makes no previously nonstrict master
strict. This is a coefficient-feasibility failure, not a proof that every
BIAS2 multiplier or the nonlinear theorem fails. It invalidates treating the
observed positive separation as the missing uniform certificate. The tested
common-gain bounds exceed the point 30-degree chart storage budgets by
1.8021e7 (H18) and 3856 (A21), even before uniform nonlinear enclosure.
All directions, 80-digit ratios, energy accounting and candidate multipliers
are in `reports/results/rao_stability/bias2-motion-gain.json`.

Critic alternatives are (1) separate bias and physical-template budgets,
(2) storage about the source-generated forced response, retaining its physical
output error explicitly, or (3) a physically scaled storage with a new
whole-word dissipativity argument instead of the ill-conditioned information
metric. The next falsifiable experiment is a complete-word forced-response
storage with separate physical output and initial-error budgets, tested on
H18 and A21 before any source covering. No interval refinement is justified
by the present gains. Physical admission, P4 and P5 remain unclosed.

## Three-route complete-word result

`reports/results/rao_stability/three-storage-routes.json` records all three
executed routes on the same attached H18/A21 source, with maximizing root
vectors, endpoint/prefix propagation at 80 digits and signed operation costs.
No source search, coefficient fitting, changed filter or reduced domain.

1. Source-centered **18-error motion** storage has endpoint ratios
   .999572658878 / .932948381611 and prefix maxima
   1.003245112488 / 1.000001212298. Full A21 error storage additionally gives
   .997124968091. Full corrected 21-state feedback remains in every map;
   initial bias is separately budgeted, not silently held at zero in the
   complete forcing inequality. Squared initial-bias gains reach
   3160.343 / 314671.332; particular-response point energies reach
   4.13762 / 16.50404. The .4 zero-true-bias capture ball plus the common
   template gives sufficient prefix bounds 600.63 / 51794.18, exceeding
   point 30-degree chart storage by 6.41 / 2.27. This invalidates using this
   unqualified separated budget for chart retention, not source-centered
   ISS or the nonlinear theorem. H18's limiting endpoint direction couples
   Y velocity, Y wave acceleration and Y displacement; A21 is dominated by
   integrated Z displacement. H18's smallest endpoint margin is .00042734.
2. Fixed SI storages give rho=1072.993 / 6659.686; fixed gravity/3-second
   storages give 11.18378 / 6.57608. These diagonal physical-metric formulations
   fail even at the coefficient point. Congruent rescaling of the information
   storage leaves rho unchanged to numerical precision. Classification:
   proof-method failure, not poor interval conditioning or a legal nonlinear
   counterexample. Freeze these metrics; finer arithmetic cannot repair them.
3. Full transported accepted-vector information is positive at this point:
   eta6 SI minimum eigenvalues 972.647 / 862.261. Keeping the signed whole-word
   matrix remainder reproduces positive motion margins .00042734 / .06705162.
   Separate scalar minima instead give negative bounds -.380703 / -.257014.
   The scalarized transport tactic is rejected; the full signed matrix route
   survives point feasibility. This is the corrected finite residual transport,
   not an automatic replacement for the canonical open-loop PE certificate.

## Declared-domain retention result

`reports/results/rao_stability/domain-retention.json` records the fourth
executed route on the same attached H18/A21 source, produced by
`tests/kalman_ou_iii/ou3_p4_domain_retention.py`. Each prefix carries a
certified subadditive upper bound and an attained maximizing functional, so a
group is retained only when its upper bound stays inside its own declared ball
and definitely violated only when its attained bound leaves it. Radii come
from `tools/stability/ou3_proof_operating_domain.json` with the same 30-degree
Cayley chart route 1 uses. No metric, no scalarization, no fitted radius, no
reduced domain.

| Reached / declared bound | H18 | A21 |
|---|---:|---:|
| Bias ball and template, worst motion group | .1713 | .2382 |
| Bias ball and template, 30-degree chart | .0703 | .0041 |
| Full declared box, worst group (velocity) | 35.358 | 9.357 |
| Full declared box without the integral ball | 4.033 | 4.921 |
| Largest retained box of the declared shape | 2.83% | 10.69% |

## Current critic and next falsifiable experiment

The strongest reason to reject promotion is unchanged in kind: every executed
route is a point diagnostic on one frozen capture, and none establishes a
uniform nonlinear matrix inequality, a BIAS2 mu, or physical source admission.
The retention result removes an obstruction rather than supplying a
certificate. The complete-word endpoint ratios .999572659 / .932948382 and
prefix maxima 1.003245112 / 1.000001212 are untouched by it, and a prefix
maximum above one still has to be discharged by retention rather than assumed
away.

The critic's own objection to this route is that it replaces one unproved
premise with another: the declared handoff box was never a reachable set, so
showing it is not invariant is weak evidence about admissible trajectories.
That objection is the next experiment. The three remaining architectures are a
correlated reachable set for the kinematically linked (p, S) and attitude
coordinates, a globally bounded particular solution from BRMM physical
primitives instead of a response restarted every word, and a full signed
vector-information enclosure with explicit source admission.

The quantitative next test is the first of those: whether the exact kinematic
identity relating the integral state to the position state, applied as a
correlation constraint on the initial set rather than as a smaller box, brings
the H18 velocity ratio 35.358 and the A21 ratio 9.357 below one. The integral
ball alone accounts for 35.358 -> 4.033 (H18) and 9.357 -> 4.921 (A21), so the
identity is necessary but cannot close either mode by itself; attitude and
position must supply the remainder. Refining the declared box without such an
identity is a dead end: it is a smaller box around a set that was already
shown not to be invariant. Physical admission, strict uniform P4 and P5 remain
unclosed.
