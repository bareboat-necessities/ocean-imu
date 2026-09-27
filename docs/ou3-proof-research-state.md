# OU-III proof research state

## Current hypothesis

One persistent physical execution follows construction -> startup/capture ->
magnetically informed H18 -> refinement/release -> recurring magnetically
informed A21 -> regional practical stability. Shipping behavior, tuning,
simulations, physical constants and quality gates are unchanged.

On every physical T_E window, MARINE MOTION requires complete stillness or
complete-interval gravity-direction span >=theta_E>0. T_E/theta_E remain
symbolic; no pinned RAO qualification supplies their numerical values.
The diagonal-wave p=-(3/2)sin(2t)(1,0,1), R=I is now inadmissible. Quiet water
and quiet-bias ambiguity remain admissible. MAGNETIC SERVICE is still actual
applied normalized heading/axial-bias information, not cadence or six-state
observability. Its assumptions have not been strengthened.

The controlling moving-history targets remain positive source-uniform
Delta_col and Delta_gyr, using signed sensor identities and actual chronological
gains, OU transport, resets, S corrections, references and physical bias history.
Both targets remain OPEN; neither has a current-domain counterexample.

## Evidence

- PR #602's retained branch started at bad0ef86823f94f2a660ea25dc7748c7a888005e.
  GitHub already marked #602 merged before this continuation; this work does
  not merge anything or change main.
- Threshold hygiene: theta_E is checked finite and strictly positive before
  stillness and excitation. Tests cover zero, negative, NaN and both infinities
  under both branches, including invalid span. Initial focused validation: 11
  tests passed; expanded signed-temporal/physical validation: 25 passed.
- Physical sampled span is >=max(0,theta_E-2 Omega_max eta). The literal signed
  accelerometer pair identity retains physical acceleration, nominal force,
  bias updates and innovations. It does not assume R_hat follows R_true.
- Exact adjoint criterion: for C_i=Phi_(N,i+1)K_i, compatibility is equivalent
  to ker C subset ker W. Otherwise the state and innovation residual sums
  remain. Zero terminal multipliers are impossible for the first nonsingular
  regular S atom c0(I-K_SS).
- An actual quiet construction-to-release observer/control run has identical
  terminal outputs. Four actual S knots are (0.025,0.160,0.295,0.430) s relative
  to the 225-s root. Its exported rational full-mean functional has a
  seven-column witness ||v||_infinity=1, Cv=0 and
  |Wv|>754573637/5000000000000000. The 80-digit value is approximately
  3.0182945480072294e-7. This is exact incompatibility of that coefficient
  word, not an enclosure of all real-arithmetic histories or all-time service.
  The actual quiet innovations vanish; no instability follows.
- Signed physical bias summation uses the total and tail sums of weights,
  then bounds physical rate. Signed physical acceleration summation uses the
  same velocity integral, weight variation and the existing jerk sampling
  error. These are useful residual supplies conditional on actual multiplier
  bounds, not nominal mean clamps or separation certificates.
- Construction BG mean is zero and the first prediction increment is
  <=0.0039051914291880918 rad. A later complete turn requires nominal bias
  norm >1046.5466859583 rad/s for allowed h. A bound on the signed cumulative
  gain/innovation sum excluding approach to this barrier is still missing.
- The historical implication is corrected: SIX uniform residual-pivot floors
  plus operator/factor bounds imply exact L O=T_h and finite B_*. The former
  claim that two proposed temporal margins automatically supply these six
  floors was unproved and has been withdrawn. The backward action bound now
  includes its terminal selector and all observation subtractions.

A conditional smooth moving family also disproves mandatory nominal tilt
response: phi=(theta_E/2)sin(2pi t/T_E), p=0, b_g=-omega,
b_a=g(R'e_z-e_z), B=75e_x produces exactly quiet measurements. The sufficient
symbolic bias/rate/service conditions are proved in the signed-temporal note;
no numerical T_E/theta_E is selected. This family meets physical tilt span
when those conditions hold, yet nominal attitude is fixed. Its nominal force
and field remain separated, so it does not falsify either target margin.

## Current limiter

The physical span does not yet control the nominal force/field or gyro
transport on every carried execution. Compatible adjoint norms, or useful
bounds on both residual sums when incompatible, must be derived without
independent nominal boxes or unsigned innovation energy. Then a quantitative
six-pivot bridge and coefficient compactness are needed. No uniform historical
B_*, J_AG>0, full covariance upper bound or rho_0<1 is instantiated.

All new lemmas enter V_(j+1)<=rho V_j+c_d||d||^2 only through these unresolved
steps. No new contraction enclosure is attempted with missing inputs; the
required high-precision feasibility check must precede any future one.

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
   these relaxed examples are unproved. Retain actual-history linkage.
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

Use observation-forced multipliers with explicit state/innovation residuals,
then exploit their actual signed temporal weights against physical velocity,
primitive and bias histories. A homogeneous compatible terminal map is allowed
only if its kernel condition and norms are proved. Gyro control must retain
construction-origin level information; zero-mean projection alone cannot
exclude the constant bias offset needed for a complete turn. Stillness keeps
the literal projection sector pending a shaped-region proof.

## Next falsifiable experiment

Construct the forced balance for a carried moving word and compute the signed
residual coefficient action before taking norms. Test whether physical
velocity/primitive and bias-rate tail sums absorb it while preserving the
sampled tilt chord and actual two-column magnetic information. If not, record
the failed residual coefficient; do not refine the failed pairwise or energy
relaxations. Only a source-uniform positive margin and quantitative six-pivot
bridge justify a common historical action enclosure. Before any new contraction
enclosure evaluate the complete matrix ratio at high precision and stop if
it exceeds one.

## Validation and infrastructure

- Serial `make all` passed using the installed Eigen 3.4.0 headers and the
  normal verified pinned release archive; all 535 shared validation tests
  passed. No shipping sources, physical constants, simulation inputs or
  quality thresholds changed.
- All 246 OU-III proof tests, 39 focused publication/workflow/provenance-shape
  tests, the exact LIN path generator, and the standalone proof evidence gate
  passed. The latter reproduces the signed certificate and verifies the exact
  carried incompatibility witness and source/observer fingerprints.
- The final shared evidence/publication suite passed all 535 tests after the
  proof changes. Repository-wide Ruff and `git diff --check` passed.
- Historical AG source diagnostics reproduced observer/control parity for all
  four cases. The 80-digit action maxima range from 1199.060966781881 to
  1450.252861111712; exact root cancellation and full correlated exported-word
  enclosures pass for the quiet and moving words. Uniform coverage remains false.
- LuaLaTeX completed two passes; the 11-page article was rendered and all pages
  visually inspected, with the new mathematics inspected at full page size.
  No overfull boxes or undefined references remain. Initial IEEE class font
  substitutions and its existing amsmath `over` warning are nonfatal; final
  body/math fonts render correctly.
- Native construction diagnostics reproduced the existing storage and
  historical-domain results. The unchanged mean-action diagnostic reproduced
  89998 predictions and 122251 corrections at 80 digits; the 40-digit outward
  enclosure and exact recorded-force audit passed. Its collinearity margin is
  still in `[-817884.048495, -817884.048494]`, so these regressions do not reopen
  the failed energy tactic or certify current-domain source-uniform coverage.

The first apt command failed with setgroups/setegid/seteuid errors in the
container. Existing runtime Eigen headers resolved native builds. The first
article compile failed on missing IEEEtran/luaotfload; matching local TeX
packages, newtx fonts/maps and the generic binhex dependency resolved it.
The CTAN redirect returned HTTP 502 and a guessed binhex path returned 404;
verified Ubuntu package extraction supplied the dependencies instead. No
article-format or quality-gate change was used to hide these environment errors.
