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
- The forced data adjoint now rewrites both compatibility residual sums
  exactly: with r=y-Hbar u+epsilon, L=W+Z_next K and
  Z=Z_next A-L Hbar, sum W r=Z0 u0+sum L(y+epsilon)+sum Z_next d.
  All 18 Euclidean means, literal gains/OU, reference, resets and projection
  defects remain. No homogeneous compatibility is assumed. Its source-uniform
  root/rotation/multiplier action remains open.
- A moving construction-to-289-s observer/control pair has identical terminal
  outputs. Four actual S knots in the tail are (.065,21.365,42.665,63.96) s;
  the 80-digit forced-balance residual is about 3.38446651750e-82. This is a
  finite diagnostic, not an interval or all-time magnetic-service certificate.
  The exact frozen-attitude Euclidean map has F E_BG=E_BG and Z E_BG=0:
  removal of this functional's root BG coefficient supplies no gyro restoring
  feedback. Alias exclusion must include the actual attitude recurrence.
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
independent nominal boxes or unsigned innovation energy. The forced data
identity retains root action and joint signed rotation/reference action;
cellwise and actual-S-interval norm relaxations failed as quantified below.
Then a quantitative
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

Couple the forced data balance to the literal attitude and BG recurrences.
Derive a joint supply for the rotation/reference and physical integral terms
using sampled tilt and actual two-column magnetic information before separate
norms. First test that proposed coupling on a carried moving word, including
the construction-origin gyro level; its frozen Euclidean map alone has an
identity BG column. A source-uniform action bound must then control root
action and coefficient variation, with a quantitative margin that exceeds
the recorded failed budgets. Do not refine the failed variation, pairwise or energy
relaxations. Only a source-uniform positive margin and quantitative six-pivot
bridge justify a common historical action enclosure. Before any new contraction
enclosure evaluate the complete matrix ratio at high precision and stop if
it exceeds one.

## Validation and infrastructure

Current forced-data continuation:

- The exact signed tests pass (18), including inhomogeneous adjoint root,
  innovation/mean arithmetic, projection, the unit BG column, sharp interval
  jerk remainder and diagnostic promotion/tamper rejection. The complete
  OU-III proof suite passes 250 tests. The final shared evidence/proof/
  publication gate, `make -C tests/validation test`, passes all 539 tests.
  `build_evidence.py`, Ruff and `git diff --check` pass; theorem_closed=false.
- The moving native observer/control parity and 80-digit signed balance pass.
  Physical g=9.80665 and literal float g=9.8066501617431640625 are separated
  by an explicit signed representation defect; no physical constant changes.
- The initial native export had zero bytes and was rejected by the terminal
  mean-continuity check. Classification: export/infrastructure failure, with
  no numerical result accepted. A direct rerun produced the complete trace;
  an explicit stream-flush/write check and full observer/control rerun now
  guard export. The retained report concerns only that complete trace.
- The first full-build session disappeared before returning status
  (`write_stdin`: unknown process id 4731); its log ended in the frequency
  stage. No pass was inferred. Re-running the same `make all` command with
  continuous session polling returned exit 0. Its copied log was incomplete,
  so the complete shared gate was additionally captured through stdout and
  returned exit 0 with all 539 tests. Quality gates and fixtures are unchanged.
- LuaLaTeX completed two passes for the final 12-page article. All pages
  were rendered/decoded and visually inspected, with the new mathematics
  checked at full resolution. No overfull boxes or unresolved references
  remain. An intermediate PDF/PNG copy was truncated (PNG at 262144 bytes;
  PDF had an invalid xref); it was not delivered. Fresh output, completed
  rendering and decoding of every final PNG resolved the artifact failure.
  Existing class-font/amsmath warnings remain nonfatal. The article is saved.

Retained earlier validation evidence:

- Native historical-reader observer/control parity passed in all four cases;
  80-digit action maxima were 1199.060966781881--1450.252861111712. Exact root
  cancellation and correlated exported-word enclosures passed for quiet and
  moving words. Construction/mean-action regressions reproduced 89998
  predictions and 122251 corrections at 80 digits, with a 40-digit outward
  enclosure and exact recorded-force audit. None is source-uniform coverage.
- The earlier reconciled build passed native stages but four archive-copy
  tests hit `[Errno 28] No space left on device`: 508 MiB was free while the
  copied simulator trees required about 2.4 GiB. Removing this task's 3.6-GiB
  temporary TeX cache restored capacity, and all 535 then-current shared
  tests passed. This was infrastructure capacity, not a proof failure;
  fixture logic and quality gates remained unchanged.
- Container apt operations failed on setgroups/setegid/seteuid; a CTAN
  redirect and guessed binhex path also failed. Runtime Eigen and verified
  local Ubuntu/CTAN package extraction supplied the needed build/TeX inputs.
  No article format or source requirement was changed to hide those failures.
