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
  Source-uniform c,b0,epsilon and factor/count bounds have not been obtained.
  The exact supplied matrix audit has s=19/100; its 80-digit singular value
  is .35070823886858529042926658460119746953531050023679. This is algebraic
  feasibility, not a contraction construction or a shipping certificate.
- The exact zero-residual quiet nominal record supplies a special two-group
  floor: actually applied acc/mag pairs eight qualified predictions apart
  have C'C>=81 I, A=I, B=(sum h_i)I with sum h_i>=4/125, and E=0.
  Thus s=18/127 and six pivots >=9/254. This is a
  nominal row theorem, compatible with physical tilt/BA ambiguity. It is
  not a full stationary covariance/noise-action or nonlinear theorem.
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
reference and complete signed gyro transport. The new six-pivot lemma supplies
a quantitative sufficient budget, not its source-uniform premises. Preserve
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
chronological two-group budget, or a sharper same-history row factorization,
to prove all six pivots. The source-uniform row/defect margin must exceed the
recorded failed budgets; a convenient supplied example cannot select it.

## Next falsifiable experiment

Derive a proposed stationary compatible-class storage and a moving joint
rotation/reference/physical-integral budget with all supplied constants
explicit. Evaluate their full finite-error ratios at high precision on the
carried quiet and already-excited moving words before any rigorous contraction
enclosure. Reject a proposed contraction construction if its worst tested
ratio exceeds one; a passing finite test remains non-promoting. Establish
uniform reset and actual-row budgets before instantiating six pivots and B_*.
Do not repeat the failed pairwise, unsigned-energy or variation-norm tactics.

## Validation and infrastructure

Based on main ed5c020e36935df09e19d12b7e8210fff00006a7, including both
concurrent calibration fixes. The final complete command
`EIGEN_INCLUDE_DIR=<Eigen-3.4.0> W3D_WRITE_TIMESERIES=0 make all
EIGEN_DIR=<Eigen-3.4.0>` returns **exit 0**, including all **571 shared tests**
and every native suite. Inputs, scenarios, gates and source are unchanged.
Only reproducible ignored CSV outputs are discarded after scoring to keep
workspace capacity available; scalar logs and evidence are retained.

The focused new suites pass 21 tests, and the OU-III proof suite passes 276.
The native 900-s persistent construction/rest/hidden-rocking/rest regression
passes: Live sample 15843, reference refinement/BA release sample 24016,
20520 actually applied magnetic corrections, hidden tilt .001 rad. Physical
bias/rate bounds, zero nominal means, covariance positivity and no reseeding
are checked. Existing stationary/rest-wave-rest and gyro-projection regressions
also pass. These finite checks do not establish an all-time stability theorem.

Exact evidence reproduction, the unchanged TFG publication manifest, Ruff
and `git diff --check` pass. Shipping-source bytes and prior finite replay
artifacts are unchanged; new exact certificates and genuinely changed proof
bindings are reproduced from their current inputs. Source-fidelity review
charges the native gravity representation residual 1.617431640625e-7 m/s^2
and uses eight qualified prediction intervals, at least 4/125 s, for the
quiet nominal row floor rather than equating .005f with exact .005.

The final IEEE article completes two LuaLaTeX passes; all 14 pages are
rendered/decoded and visually inspected. No overfull box, unresolved reference
or balancing warning remains. Its source hash matches the verified PDF run.

Resolved failures and retained scope:

- **Publication fingerprint:** the first full run failed at
  `python3 ../../tools/tfg_comparison.py --check --require-ci` with
  `ValueError: TFG comparison source provenance is stale` (make exit 2).
  Exactly one recorded input differed: the OU-III simulator Makefile acquired
  the new auxiliary target. Restore that Makefile byte-for-byte and use the
  existing supplemental-Makefile pattern for the regression. The original
  manifest passes without restamping or changing its gate; native results
  and simulator behavior remain valid.
- **Execution/resources:** an updated-main retry lost handle 40461; an
  awaited retry then failed with `exec-server transport disconnected;
  executor key changed during session recovery` during calibration compilation.
  The local cgroup reported limit 8,589,934,592 bytes, peak 8,592,605,184
  bytes and an OOM kill. Neither interrupted run is counted as complete.
  Automatic command review rejected an `rm -f` cache-removal attempt;
  recovery instead forced the calibration target through Make, using GCC
  `ggc-min-expand=10` and `ggc-min-heapsize=4096`. O3 and all numerical gates
  remain unchanged. That rebuild and the final supervised full run pass;
  a durable exit-status record matches the complete log hash.
- **Eigen discovery:** the next full run passed 570/571 shared tests but
  `test_device_compass_startup.py:80` failed with `Install Eigen or set
  EIGEN_INCLUDE_DIR`. Make's EIGEN_DIR does not configure that Python-compiled
  device test. Point both variables to the same Eigen tree. The focused
  1800-orientation checks in all three variants and the final full run pass;
  no assertion or dependency requirement is relaxed.
- **Article output:** the initial command using only TEXMFHOME failed with
  `font txmiaX at 438 not found`; the existing extracted font maps via
  TEXMFVAR resolve it without changing IEEE/newtx. An unused alignment
  column caused a 4.27032pt overflow and was removed. Optional balancing
  then produced a second-column warning and a 16.61555pt page-13 overflow;
  restore standard IEEE flow with all equations retained. Generated native
  CSVs also filled scratch during an intermediate TeX run; its partial log
  is not accepted. Capacity cleanup, two complete passes, rendering and the
  source-hash check establish the final article result.
