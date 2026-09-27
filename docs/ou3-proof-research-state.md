# OU-III proof research state

## Current hypothesis

One persistent MARINE MOTION / IMU BIAS / MAGNETIC SERVICE execution follows
construction -> capture -> informed H18 -> refinement/release -> recurring
A21 -> regional practical stability. STILL / TRANSITION / MOVING qualify
portions of this same route. No shipping mode switch is enabled.

Excitation means Delta_g([t,t+T_E])>=theta_E on **every complete window
contained in one maximal physical moving episode**, between nondegenerate
rest intervals. Isolated zero-rate instants do not split an episode. Crossing
windows do not owe excitation. T_E/theta_E remain symbolic. The old global
rolling stillness-or-span rule is removed, not repaired by more assumptions.

The stationary target must allow physical attitude/BA ambiguity and cover
all measurement-compatible histories before quiet evidence changes estimator
behavior. The existing assumptions admit exact rest/motion indistinguishability;
a sound exact-rest detector with finite entry and universal finite exit is
impossible on that admitted parameter family. This does not establish that
no estimator change is needed. Full stationary practical boundedness is OPEN.
The moving signed margins Delta_col/Delta_gyr, six-pivot uniform premises,
B_*, J_AG, rho_0 and nonlinear every-prefix retention remain OPEN.
The exact zero-residual quiet nominal subfamily now has a root-independent
historical action/full covariance ceiling and qualitative homogeneous linear
loss. For MOVING, inverse-frame transport and exact same-prediction-cell
groups sharpen the existing reader's sufficient premises; they do not supply
uniform positive geometry, injection or action bounds.

## Evidence

- `ou3-regime-design.md` states observable stationary information, required
  practical theorem, detector obligations, carried transition state and the
  obstruction before any shipping change. Terminal gyro-bias averaging has
  exact supply N_g+Omega_tube+D_g sum(weight*sample_age); coherent deterministic
  noise does not average away. Gravity/BA retain a one-dimensional kernel
  about the known magnetic axis. Existing projections alone do not close it.
- The sin-cubed rest/motion/rest witness uses alpha=1/1000, nu=1/40,
  p=v=a=0, B=75e_x, b_g=-phi' e_x and b_a=g(R'e_z-e_z). Its bias/rate
  upper bounds are .00980665, .00073549875, .000075 and .000005625 in
  their respective SI units. All lie strictly inside unchanged limits;
  joins are C2, jerk/primitive vanish and complete 80*pi windows have span
  .002. The measured packets are exactly level/north for arbitrary duration.
  Existing actually-applied stationary service retains a floor >1 after
  the exact rational lower factor (1-alpha^2/2)^2. This is a witness parameter
  instance, not a deployment qualification of symbolic T_E/theta_E.
- `regimes.py` implements complete-window requirements, exact conditional
  every-prefix bridge composition and deterministic stationary gyro supply.
  Its constant-memory necessary-evidence monitor uses physical/noise bounds,
  a supplied continuous dwell and immediate evidence-failure exit. Invalid
  packets/gaps restart dwell. It never returns a physical-STILL certificate.
  Estimated wave states cannot distinguish identical input histories.
- `ou3-moving-six-pivots.md` continues the same-history reader: chronological
  prediction/reset recursion bounds |B-tI|, then the actual two-group row
  budget gives s=c/[1+(a+1)/b0]-epsilon. If s>0, every one of six greedy
  row pivots is >=s/sqrt(m). All asynchronous/reset/reference defects remain.
  Applied acc/mag groups within a prediction cell now give E=0 exactly by
  retaining every intervening reset in C_i. No simultaneous-event assumption
  is made. Source-uniform c,b0 and factor/count bounds remain unproved.
  The exact supplied matrix audit has s=19/100; its 80-digit singular value
  is .35070823886858529042926658460119746953531050023679. This is algebraic
  feasibility, not a contraction construction or a shipping certificate.
- The exact zero-residual quiet nominal record supplies a special two-group
  floor: actually applied acc/mag pairs eight qualified predictions apart
  have C'C>=81 I, A=I, B=(sum h_i)I with sum h_i>=4/125, and E=0.
  Thus s=18/127 and six pivots >=9/254. This is a
  nominal row theorem, compatible with physical tilt/BA ambiguity. It is
  not a physical identifiability or nonlinear theorem. The same historical
  reader theta_read=z1, bg_read=(z1-z0)/T now cancels arbitrary AG root and
  cross covariance. Inherited AW/BA, correlated process and actual observation
  bounds give B_q, an every-operation full upper C_q, and qualitative
  homogeneous linear loss on this exact subcase. The full compatible-class
  nonlinear theorem remains OPEN; covariance retention is not state retention.
- Literal reset G=I+[d]/2 and both prediction R branches have inverse norm
  <=1. C=A^-1 B is unchanged by resets; later gyro injections retain A^-1.
  Exact prefix recurrences bound |C-tI| without amplifying its past defect
  by every reset. A terminal d=4e_x no longer erases an earlier .005-s floor.
  The 80-digit supplied noncommuting audit has floor .0396561974478783
  versus singular value .0399999504651045. The quiet correlated-root action
  audit ratio is .413929501952256 at root scales 1 and 10^12. These algebra
  audits precede enclosure and are non-promoting.
- A fresh carried-source audit reuses the existing read-only observer and
  unchanged driver, with exact untapped terminal parity on quiet and moving
  inputs. Each has eight complete same-cell groups. First/last group E=0
  is exact over exported rational operators; 80-digit six-column budgets
  are 1.204325434846729 and 1.074485913647320, below actual singular values
  1.922502992341225 and 1.717587283388896. These finite numbers are neither
  uniform certificates nor real-trajectory or magnetic-service enclosures.
- Same-cell geometry obeys c^2>=min(|f|^2,|b|^2)(1-|u'k|), where u=f/|f|
  and k is the normalized G_local^-1 b. Literal reset noncontraction proves
  the two magnetic nonzero singular values are >=|b|. This enters c in the
  same six-pivot inequality. The supplied 80-digit floor 8.54167423631316
  is below singular value 8.80825850549787; no physical uniformity is inferred.
- The current excited moving forced-balance word already has roll
  .02 sin(t/2) and displacement .4 sin(.6t)e_z. It satisfies the moving span
  premise for T_E>=4*pi, 0<theta_E<=.04. Its 80-digit balance residual is
  about 3.38446651750e-82, with exact observer/control terminal parity.
  Its failed rotation and variation budgets remain failed after removing
  indefinite stillness; see DEAD_ENDS 13--14. No scalar refinement is retried.
- The exact carried-adjoint compatibility criterion is ker C subset ker W.
  Otherwise both residual sums remain. The forced-data identity retains
  root action, joint signed rotation/reference, literal gains, all Euclidean
  means and reset/projection/rounding defects. In the frozen-attitude map,
  F E_BG=E_BG and Z E_BG=0; that cancellation is not gyro feedback.
  The quiet word's exact seven-column null witness has |Wv|>
  754573637/5000000000000000, while its actual quiet innovations vanish.
- Existing shipping projection evidence is retained unchanged: residual
  gyro-bias <=.5 rad/s, qualified prediction angle <.007, one-step transport
  floor .003999991833333333 s. A complete turn would require bias norm
  >1046.5466859583 rad/s on the qualified 4--6 ms family. None of these is
  a complete signed temporal gyro margin. No projection was active in the
  verified 840-case validation plus 18 calibration and 310-case robustness
  plus one calibration paired studies; those source/input-bound replays
  are not regenerated or restamped by this proof-only continuation.

## Current limiter

Exact physical rest is not identifiable from the current sensor/bias model.
A stationary practical theorem for the entire compatible class is required
before shipping detector/state/covariance changes can be justified. Conditional
finite bridge algebra does not supply detection liveness, a retained set or a
budget for arbitrarily repeated switches. These are explicit OPEN obligations.

On MOVING windows, physical span still does not control the nominal force,
reference and complete signed gyro transport. The six-pivot lemma supplies
a quantitative sufficient budget, not its source-uniform premises. Same-cell
grouping removes E for this selection rule; the inverse-frame formula still
needs actual injection control to give b0. Quiet-subcase homogeneous decay
does not control compatible physical mismatch. Preserve
actual chronological gains, resets, OU and bias histories. No independent
nominal boxes, unsigned energy, sampled Gramian or finite replay closes B_*.
Every new lemma's role in V_next<=rho V+c_d|d|^2 is stated in the design:
stationary information sets a supply/ambiguity radius, bridge products and
supplies compose rho/c_d, and six pivots feed the historical covariance/loss
comparison. Source-uniform numerical contraction enclosure is not yet justified.

## Failed approaches / DEAD_ENDS

1. **Unsigned cumulative mean-action energy: DEAD_END.** E=817885.0623259853
   while collinearity cost is about 1.0138313892; exact E_col-E lies in
   [-817884.048495,-817884.048494]. Invalidated: the correlated energy
   ellipsoid separates all nominal forces from the field. Retain its finite
   old-domain BG<1 and force/field sine>2/5 certificates, actual factors and
   arithmetic bounds. Do not refine precision/subdivision of this relaxation.
2. **Independent nominal coefficients: exact nullspace.** f parallel B gives
   historical AG rank four; the complete-turn nominal h=.005,
   bhat_g=-400pi e_z example also has rank four and information-floor margin
   -mu. Invalidated: physical bounds or innovation bounds alone control a
   free nominal root. Construction reachability and all-time service for
   these relaxed examples were unproved. The implemented gyro projection now
   excludes the full-turn bias state; force/field collinearity and full temporal
   rank remain separate. Retain actual-history linkage.
3. **Endpoint-free forced adjoint: missing compatibility.** The exact failed
   equations are Z_i=Z_(i+1)A_i and W_i=Z_(i+1)K_i together. Zero-mean
   projection generally destroys them; the finite carried witness above
   also fails with unrestricted terminal multiplier. Invalidated: zero mean
   or spline jets alone leave only physical bias increments and defects.
   Retain both residual sums and signed physical Abel/velocity supplies.
4. **Pairwise unsigned tilt transfer: nonpositive margin.** The subtraction
   min(2A_max,(J_max+Omega_max A_max)T_E) already exceeds the maximum available
   gravity chord for every feasible theta_E<=pi/2. Invalidated: a pairwise
   triangle bound excludes nominal collinearity. Retain complete-window span
   and the exact signed pair relation; use whole-window physical integrals.
5. **Two margins to six pivots: unsupported implication.** No rank theorem
   proves that force/field separation and nonaliasing exhaust all varying
   six-column loss mechanisms. Coefficient compactness is also missing.
   Invalidated: simply taking p=min(Delta_col,Delta_gyr). Retain exact
   historical factor pivoting, root cancellation and the conditional p>0 bound.
6. **Universal V<=36 capture: refuted by admitted stillness.** For true
   R=Rx(2 atan(1/100)), p=v=a=omega=0, B=75e_x and
   b_a=g(R' e_z-e_z), ||b_a|| approximately .1961232<B_a, the measured record
   is level/north. The stationary service floor remains >1 after multiplying
   by (9999/10001)^2. At applicable regular boundaries P_ba,ba<=I/1600 gives
   V>=38468153689/625062500>61.54289>36. Invalidated: universal entry into
   the convenient projection-inactive ball or exact attitude/bias convergence
   in indefinite stillness. Retain the literal bias-projection sector;
   no larger/shaped invariant region has yet been proved.
7. **Old constant-attitude moving diagnostics: domain exclusion.** The smooth
   diagonal wave has finite 400--600-s mean tilt 8.211611 degrees, nominal
   acceleration 9.776391 and endpoint V>11492.6752 (guard margin <-11456.6752).
   These numbers are retained as diagnostics of the older domain, not current
   physical counterexamples. The 200-Hz stationary-looking moving witness is
   also excluded by jerk and tilt span. Its old-domain service proof remains.
8. **Local arithmetic factors: signed defect, not PSD noise.** The exported
   complete sync/symmetry defect has Rayleigh quotient -2^-44. Retain the
   full signed rank-one/two factor enclosure. Invalidated: treating that
   defect as a PSD process increment or restricting it to three AW coordinates.

9. **Future-only AG loss with a free prior: analytic obstruction.**
   P_root=diag(t I6,I15) gives D_AG,AG<=I6/t. At t=10^12 the margin against
   10^-6 I6 is <=-9.99999e-7. Retain historical root cancellation; future
   excitation alone cannot establish absolute J for that relaxed prior.
10. **Entrywise Riccati and cross ceilings: DEAD_END.** The former prediction
    enclosure [-18.7907040,2502.41442] loses acc/mag inverse verification;
    the cross-ceiling/Gershgorin margin is about -5.1223e8. Retain full matrix
    factors. LDL pivots are not eigenvalue floors and more entrywise
    subdivision is not the active approach.
11. **Restricted information lifting and endpoint process sums: invalid.**
    Heading/bias restriction I2 can coexist with full loss
    [[1,0,1],[0,1,0],[1,0,1]], cancelled by nuisance vector (1,0,-1).
    Interleaved corrections also defeat uncorrected accumulated-process
    lower bounds. Retain nuisance elimination and corrected endpoint paths.
12. **Automatic vanishing reset remainder: invalid.** At fixed nonzero
    injection d, the reset derivative J_l(d) differs from I+[d]/2. Smoothness
    alone does not prove eta(r)->0 for that comparison. Retain the explicit
    injection-dependent remainder and require a complete strict word margin.
13. **Forced-data rotation triangle: failed relaxation.** On the carried
    moving word, charge 39.22565272637 exceeds recorded projected gravity
    8.77133455729: margin -30.45431816909 before other supplies. Invalidated:
    separately norming the two rotations can close this budget. Retain the
    exact jointly signed rotation/reference action and actual magnetic gains.
    This does not falsify Delta_col or the physical assumptions. Current
    limiter: no uniform coupled rotation/reference bound. Next: include
    actual attitude dynamics and magnetic information before taking norms.
14. **Physical acceleration variation norm: stop after one refinement.**
    The cellwise velocity+jerk charge is 239.97053019538. Summing signed
    weights over actual S intervals first gives 23.98805965312+3.82238280532
    =27.81044245844, still margin -19.03910790115 against the same recorded
    threshold. Invalidated: either variation-norm bound closes the budget;
    no precision/subdivision retry is justified. Retain both exact integral
    identities, bias-rate Abel supply and the actual signed acceleration
    action (norm about .005 on this word). Architecture review: the single
    construction/H18/A21 path remains; the missing joint physical-integral,
    rotation/reference and gain bound cannot be replaced by separate norms.
    Finite small actions do not supply source-uniform ceilings.

15. **Exact-STILL entry/exit detection: analytic identifiability failure.**
    The new sin-cubed witness has exactly the same full IMU/magnetic record
    as true rest, including across C2 joins. Any causal detector entering on
    rest follows the identical state sequence through arbitrarily long hidden
    motion. Invalidated: finite dwell, gyro/acc variation or estimated v/aw
    alone guarantee both sound exact-rest entry and finite exit. Retain the
    stationary gyro information and physical ambiguity tube; no instability
    or failure of signed nominal margins follows. Limiter: no practical
    theorem for the whole compatible class. Next: derive that quotient/set
    storage before considering stationary mean or covariance changes.
16. **Finite bridges alone: switching-composition failure.** A moving norm
    gain 1/2 and finite bridge norm gain 3 give repeated gain (3/2)^n.
    Invalidated: individual finite retention establishes recurring stability.
    Retain exact every-prefix product/supply composition. Limiter: uniform
    product and transported-supply control through all switches. Next: derive
    that budget from actual operations or prove eventual regime retention;
    do not add a convenient physical switching/dwell assumption.
17. **Unrestricted reset inverse nonexpansion: exact algebraic obstruction.**
    With zero corrected rate, predict h0=3/625, apply d=4e_z,4e_z,(8/3)e_z,
    then eight h=1/200 predictions. Every literal reset is nonsingular with
    inverse norm <=1, yet B=diag(0,0,28/625). Invalidated: inverse nonexpansion
    and each one-prediction floor alone force full inter-anchor gyro rank.
    Retain exact inverse-frame recurrence and same-cell rows. This relaxed
    sequence has no shipping reachability or magnetic-service proof, and
    intermediate rows could restore full historical rank. Limiter: actual
    injection/geometry bounds. Next: bound those coupled to applied service,
    not by deleting resets or shrinking independent nominal boxes.
    Likewise the relaxed same-cell f=(1,0,1), b=e_x, d=2e_y gives Cf=0
    although raw f,b are nonparallel. The pulled-back field direction,
    not raw simultaneous separation, is the required geometry. Neither
    relaxed witness establishes source reachability or complete-row rank loss.

## Retained facts

The sampled acceleration mean bound, constant-field joint physical vector
floor, full covariance-energy identity, LIN path action, full nuisance floor,
recurring nuisance upper covariance, historical reader factor/root algebra,
complete corrected-word matrix bootstrap, actual-gain finite-error composition,
reset remainder and projection-sector relations remain available.
The physical signed reserve is 27049050188592/625000000000000000, approximately
4.32784803017472e-5, under its stated physical field hypotheses. Its nominal
transfer ceiling remains unknown. The six-column-to-full contraction lemma
is conditional, not a certified numerical contraction for shipping histories.

The retained supplied-word 80-digit comparison has rho approximately
0.995037970315 and conditional decrement approximately 2.64851e-5; it is not
a source-uniform shipping result. Coarse existence quantities (nuisance ratio
approximately 1.02976e-18, scalar process floor approximately 2.49159e-26)
are not practical-radius certificates.

The theorem default profile uses g=9.80665, whereas the AtomS3R calibration-site
configuration uses default g=9.8025605 and its deployed S cadence. Other
configured profiles require their own source enclosure. No constants or
deployment behavior are changed to align them with the proof.

No H18 full-state contraction, independent nominal boxes, fitted/sampled rho,
scalar information lifting, scalar/Gershgorin process comparison, or finite
replay promotion is substituted for the retained proof path.

## Alternatives

For stationary operation, pursue a set/quotient practical storage retaining
physical tilt/BA ambiguity separately from estimator covariance. A gyro-only
observation must charge the entire compatible physical-rate tube. Do not hold
BA, force zero velocity, or shrink covariance solely from quiet-looking data.

For MOVING, couple the signed forced-data identity to actual attitude/BG
transport and magnetic reference action before taking norms. Use the explicit
inverse-frame budget and exact same-cell row factorization
to prove all six pivots. The source-uniform row/defect margin must exceed the
recorded failed budgets; a convenient supplied example cannot select it.

## Next falsifiable experiment

Extend the quiet-subcase reader to a rigorously delimited compatible tube,
charging physical tilt/BA ambiguity and actual nonlinear supplies. On MOVING,
derive coupled bounds for sigma_min([C_acc;C_mag G_local]) and inverse-frame
gyro transport from actual gains, innovations and applied service. The E=0
identity permits testing these two budgets directly on carried histories.
Evaluate any resulting complete finite-error ratio at high precision before
rigorous contraction enclosure; reject a construction if its tested ratio
exceeds one. Passing finite tests remain non-promoting. Uniform bounds are
required before instantiating B_*, J_AG, rho and nonlinear retained radius.
Do not repeat the failed pairwise, unsigned-energy or variation-norm tactics.

## Validation and infrastructure

Based on main `1c8776e7`, including the latest calibration fixes. The complete
command `EIGEN_INCLUDE_DIR=<Eigen-3.4.0> W3D_WRITE_TIMESERIES=0 make all
EIGEN_DIR=<Eigen-3.4.0>` returns **exit 0**, with all **589 shared tests**
and every native suite. GCC memory-management options
`--param ggc-min-expand=10 --param ggc-min-heapsize=4096` preserve O3 and
numerical behavior while compiling the updated-main calibration sources.
Only completed, ignored generated CSV outputs are discarded after scoring;
inputs, scenarios, gates and scalar logs are preserved.

The focused OU-III suite passes **294 tests**, including 18 new inverse-frame,
same-cell geometry/grouping, covariance and fail-closed scope regressions.
The fresh two-profile carried-source audit preserves observer/control terminal
parity and verifies E=0 over exported rational operators. The persistent
900-second construction/rest/hidden-rocking/rest regression, existing
stationary/rest-wave-rest and gyro-projection regressions pass in `make all`.
These finite checks do not establish an all-time stability theorem.

Exact evidence reproduction, the unchanged TFG publication manifest, Ruff
and `git diff --check` pass. OU-III estimator bytes and previous source replay
artifacts are unchanged. New/changed proof artifacts are generated from their
actual current inputs. The final IEEE article passes two LuaLaTeX runs and
all 15 pages render and are visually inspected, with no overfull box,
unresolved reference or balancing warning. The reported source hash matches
the verified PDF run.

Resolved findings and remaining validation limit:

- **Article review:** a literal tab had replaced the backslash in theta_E.
  Restore `\theta_E`; the new render displays the intended excited range.
- **CI clang-tidy on inherited calibration logging:**
  `tools/quality_gates.sh clang-tidy` failed at `AtomS3R_CalLog.h:23` with
  `bugprone-not-null-terminated-result`. The binary CRLF packet is consumed
  only by `write(pointer,len+2)`. A narrowly scoped annotation documents that
  intentional non-string use; no runtime operation or global gate changes.
  The full rebuilt native calibration/logging tests pass. Local clang-tidy
  installation failed with `Could not open file /var/cache/apt/archives/partial/
  ... - open (13: Permission denied)` (exit 100), so its rerun remains a CI
  obligation, not a claimed local success.
- **Prior resolved infrastructure:** a stale TFG fingerprint was fixed by
  restoring the original simulator Makefile and using the supplemental
  regression Makefile, without restamping the manifest. Earlier compiler
  OOM/session loss was addressed with the memory options above; both Eigen
  variables are required by the separate Python-compiled device test.
  None of those interrupted runs was counted as a successful full gate.
- **Publication reference resolution:** the named branch fetch returned
  main's tree `04ab4912` instead of the published tree `4d4c65ef`; the
  verification assertion stopped checkout synchronization without altering
  the validated index. Fetching immutable commit `3009b003` returned the
  exact expected tree. This was a ref-resolution failure, not changed proof
  inputs. Use the published commit ID for subsequent synchronization checks.
