# OU-III controlling proof obligations

The exact status is reproduced by `theorem_status.py` and committed in
`reports/results/ou3_stability/theorem-status.json`. Algebraic infrastructure
is distinct from a shipping source-uniform theorem certificate.

MARINE MOTION retains acceleration jerk <=100 m/s^3. Excitation is required only on every complete T_E window contained in one
maximal physical moving episode between nondegenerate rest intervals. Isolated
zero-rate instants do not restart an episode. T_E/theta_E remain symbolic.
Windows crossing rest boundaries carry finite-transition obligations; no
arbitrarily short departure window owes a full positive span. Read
`ou3-regime-design.md` for the exact quantifiers, stationary observability,
indistinguishable rest/motion histories and detector requirements. Quiet
packets cannot certify physical STILL with both finite entry and guaranteed
finite exit under the existing bias bounds. The shipping estimator is unchanged;
stationary practical robustness and certified transition retention remain OPEN.
`ou3-sampling-fidelity.md` proves the resulting sharp sampled-mean bound,
32-s exclusion of fixed-attitude stationary-sample aliases, joint 3-D vector
information and physical LIN prediction supply in the full covariance metric.
The full corrected state loss and general capture remain open.

| Obligation | Current state | Required certificate |
|---|---|---|
| Stationary physical observability | Partial information and ambiguity CLOSED; practical stability OPEN | Direct stationary gyro information has a deterministic noise/bias-rate supply; gravity/BA retain an irreducible ambiguity. A21 nominal detectability is not physical attitude/BA separation |
| Exact quiet nominal covariance subcase | Historical action, every-operation full upper and qualitative homogeneous linear loss CLOSED | `stationary_covariance.py`: two actual normalized acc/mag groups cancel arbitrary AG root and cross covariance; inherited nuisance/process bounds give C_q. This does not cover general quiet-compatible inputs or nonlinear physical stability |
| Certified STILL detector | Exact-rest entry plus universal finite exit impossible under current bounds | `regimes.py` gives smooth rest/motion/rest identical inputs; prove a practical theorem for the entire compatible class before enabling behavior changes |
| Regime transitions | Complete-window quantifiers and conditional every-prefix composition CLOSED; duration/retention OPEN | Carry the entire execution through both directions, including startup/H18/refinement/release; derive a cumulative storage budget for recurring bridges |
| Six historical pivots | Conditional chronological transport/two-group implication CLOSED; source budget OPEN | `ou3-moving-six-pivots.md`: s=c/(1+(a+1)/b0)-epsilon>0 implies every greedy row pivot >=s/sqrt(m). Same-prediction-cell actual acc/mag rows have E=0 exactly with all resets retained in C_i; uniform geometry and transport bounds remain required |
| Inverse-frame chronological gyro transport | Literal inverse nonexpansion and conditional prefix bound CLOSED | C=A^-1 B cancels resets exactly; inverse A carries their effect on later predictions. Large relaxed reset sequences can still cancel B; reachability and full historical rank loss are not established |
| World-frame historical rows | Factorization, reset Gram identity and attitude-invariant same-cell geometry CLOSED | `ou3-world-frame-rows.md`: rows are `-R_k[f_k]x[A~_k R_0',B~_k]`; attitude and its error enter only through world injections and the nominal rotation integral |
| Literal injection budget | Loewner lemma CLOSED; multi-second norm-summed use FAILED | `dd'<=NIS K S K'<=NIS P_theta,theta`. Summed norms overcharge 3-s transport about 460 times; the signed world injection sum is the retained refinement |
| Signed world-injection transport | Lemma I* and third-order reset factor CLOSED; perturbative 16-s charge FAILED (DEAD_END 21); half-angle remainder OPEN | `signed_injection.py`: `M_N M_0'=J K` gives `angle(J)<=alpha_0+alpha_N+int|omega_tilde|` without norm sums; `N=I-[x]/2+O(|x|^3)`; the rotating-frame Corollary A** charges a signed mean of relative rotations |
| Same-cell uniform floor | OPEN; not implied by motion/bias bounds alone | A 1-Hz collinear MARINE MOTION/IMU BIAS history degenerates every same-cell group, but its cadence fails MAGNETIC SERVICE (lambda_min 2.03e-5). Jerk forbids all-collinear cadences with length-weighted mean gap below .051 s (h=1/5, L=16 s). A floor needs that cadence coupling; aggregate rows avoid it |
| Aggregate world-frame six-column floor | Corollary A* (nominal signed mean) CLOSED; pointwise physical AW tracking REFUTED (DEAD_END 20); Lemma T and injection-free Theorem G0 CLOSED; nominal window statistics and injections OPEN | `aw_tracking.py`, `aggregate_floor.py`: transverse nominal mean `<1.96133 m/s^2` suffices for the attitude columns; 1-s service gaps and the nominal-rate bound make the field-axis gyro coordinate monotone; `s^2>=1.486786e-3` under `m_perp<=2/5`, `u1<=6/5` (carried worst .348, 1.091) with A~=I |
| Joint recurring lower covariance | CLOSED in real arithmetic at regular A21 post-prediction roots after a 16-s window | `root_covariance_certificate.py`: convex combination of fresh AG/BA injection and corrected LIN matrix action, with all cross covariance retained |
| Full A21 information/loss and rho0 | OPEN; reduced to one kernel-bounded word diameter; covariance-ceiling route existence-only (DEAD_END 22) | `ou3-corrected-word-proof.md` section 7: `kappa_W=lambda_max(J^-1 A)=lambda_max(Pi^-1 P_diff)` gives `rho_W<=tanh(log(kappa_W)/4)` for every root covariance. With the physical tilt/BA kernel `nu`, `A<=kappa J+lambda nu nu'` and the scalar `nu'P_0 nu<=c` suffice. Open: a source-uniform ceiling on the kernel-bounded diameter `kappa_nu` and a tilt ceiling about the body field axis |
| Estimated gyro bias / one-step gyro transport | CLOSED in real arithmetic on the qualified 4--6 ms family | Implemented norm <=.5 rad/s; angle <.007; transverse singular floor .003999991833333333 s. The API itself has no maximum positive dt; whole-word float32 and signed Delta_gyr are not closed |
| Historical AG readout action | Conditional matrix implication CLOSED; source-uniform action OPEN | `ou3-ag-readout-proof.md`: cancel the historical AG root exactly using a six-column reader, retain every process/nuisance correlation, and bound its 6x6 action uniformly over actual varying coefficients/resets. Carried quiet/moving probes and an exact exported-word enclosure do not provide this uniform bound; independent nominal force/field ranges still fail at rank-four arrays; the historical full-turn bias relaxation is excluded by the implemented gyro invariant |
| Word Riccati diameter and kernel reduction | Exact identity, composition, slow/fast factorization, rank-one kernel corollary, invariance implication and S-chain cancellation CLOSED; source-uniform gyro-bias persistence cap CLOSED; source-uniform `kappa_nu` OPEN | `word_diameter.py`: every word has `kappa_W>=sigma_g^2(1-delta)/((1+e)^2 b0 T^2)` (71.2, 17.8, 4.45 at 16, 32, 64 s). Carried words (`information-ratio-source-feasibility.json`) pass the joint-reader kill criterion on MOVING words and near steady state on quiet words |
| LIN/BA nuisance covariance upper bound | CLOSED for the regular default A21 profile after 17 s | `nuisance_upper_certificate.py` and `ou3-nuisance-upper-proof.md`: cancel the neutral root with three actual S observations; bound OU forcing, source Q defects and actual PSD sync; retain all nuisance cross covariance |
| Pre-prediction nuisance floor and coupled upper implication | CLOSED as stated algebra and source nuisance bounds | Propagate the embedded full nuisance floor through actual acc/S corrections and resets; Schur complement retains all cross covariance. The full upper bound still requires the open six-column J premise |
| Full covariance upper bound | OPEN; needed for coercivity, not for rho0 | Certify `I_eff>=mu` (hence `B_*`) under varying realized coefficients and actual corrections/resets for the nonlinear supplies. The word contraction needs only the scalar kernel variance, whose BA part is the proved `P_ba<=I/1600` |
| Actual-gain finite-error composition and finite-angle reset | CLOSED as operation/word inequalities | Freeze coefficients generated by the actual estimator, retain additive defects, source polynomial injection and the injection-dependent linear reset remainder; no two independently varying filters are equated |
| Explicit nonlinear retained radius | OPEN | Bound the complete nonlinear remainder, including projection/reset/tuner behavior, against the verified strict linear margin |
| Projection in the covariance storage | Prefix guard CLOSED; retention OPEN | The inherited BA marginal gives projection defect zero when the actual pre-projection sqrt(V)<=6; source margin is .02483339501604595 m/s^2. No entry or prefix invariance is inferred |
| Whole-word float32 supply | Composition and recorded sync-operation enclosure only | Literal operation counts, magnitude envelopes and certified prefix gains; real-arithmetic covariance positivity is not float32 totality |
| Finite startup/capture | OPEN under the jerk/tilt-span domain | Fixed-attitude stationary-looking aliases >=6 degrees are excluded; prove actual moving-attitude capture and retention through the shipping proxy, reference refinement and release |
| Finite H18 retention | Composition only | Actual entry set, history-dependent bridge duration and supply/gain bounds |
| Reference refinement and bias release | Conditional captured-domain completion only | Retain the actual refinement state machine, accepted-update count, one-second guard and covariance release |
| Release into the A21 region | OPEN | Compare the certified release set with the nonlinear retained region/projection-sector bound |
| Every-prefix tail retention | Composition only | Uniform bounds at every intermediate operation, including between recurring prediction roots |
| Recurring magnetic service | Continuation schema only | One-history every-window certificate from actually applied informative corrections and actual innovation covariance |
| Physical/sensor/bias qualification | Composition only | Simultaneous MARINE MOTION, IMU BIAS and MAGNETIC SERVICE continuations plus assembled sensor/mount/calibration qualification |
| Implementation/arithmetic totality | OPEN | Every finite branch, innovation solve, projection/reset, scheduler and arithmetic enclosure |

Do not revive the retired entrywise Riccati subdivision, endpoint-batch floor,
restricted-information lifting, or scalar normalization of an unproved full
Gramian. Local interval kernels remain available as conditional algebra; no
recurring interval box is currently certified. The active construction retains
matrix factors and the full covariance-energy loss identity.

The complete stability claim remains unproved under the revised domain. Consult `ou3-proof-research-state.md`
for the current failure classification and next falsifiable experiment.

The exact forward-only relaxation fails `D_AG,AG >= 10^-6 I6` at the supplied
prior `diag(10^12 I6,I15)`: its Rayleigh margin is at most `-9.99999e-7`.
This prior is not certified shipping-reachable or magnetically serviced.
It diagnoses a proof prerequisite, not physical instability. The next source
test is `L(W)O(W)=T_h(W)` with uniformly bounded historical action `B_W`,
followed by a matrix process/loss comparison. The measured-vector Gramian
cannot replace these nominal, reset-transported six-column rows.

The stationary A21 detectability proof in `ou3-stationary-detectability.md`
finds no nondecaying unobservable mode at rest with nonparallel gravity and
magnetic field. It is not a uniform varying-history theorem. Quiet water is
not excluded. The old sampling obstruction uses nonzero motion and is excluded by the
jerk condition, not by deleting quiet water.

The all-row noise-action Schur test also fails on the independent-coefficient
relaxation: B=(45,0,45), nominal a=(-g/2,0,g/2) gives AG rank four and
Rayleigh margin -1 against I6. This family has no proved nominal-mean or
magnetic-service reachability. The missing certificate is a quantitative
same-history exclusion of sustained nominal force/field collinearity and
positive signed temporal gyro margin, followed by the common full
matrix-action ceiling. The implemented gyro-bias ball already excludes
complete-turn bias aliases on the qualified prediction domain. Finer
pivots, precision or coefficient subdivision cannot remove an exact nullspace.

A historical, unprojected regular-root relaxation has quiet truth,
h=.005 s and nominal gyro bias -400 pi e_z. All innovations vanish, but the
literal full-turn bias transport is h e_z e_z': two gyro columns are invisible
and the Gram-floor margin is -mu. It violates the shipping .5 rad/s gyro-bias
invariant and is no longer an admissible nominal state. The one-prediction
transport floor is certified in `ou-gyro-bias-projection.md`; this does not
close the full chronological signed margin. The moving word's
complete sync/symmetry arithmetic is now enclosed with signed matrix factors;
its indefinite -2^-44 defect is charged, not treated as PSD process noise.
This local finite-word closure does not close whole-word float32 supply.

The full-construction stress history in `construction-history-feasibility.json`
obeys the older motion/bias envelopes, violates the current tilt-span
condition because R=I, and reaches the literal
refinement/release. Its finite 400--600 s tilt remains about 8.21 degrees;
no six-degree entry certificate follows from those stage flags. All-time
magnetic service, eventual capture and a uniform reachable-set enclosure
remain unproved. Construction-linked joint bounds, not free nominal boxes,
control the next step.

The full construction mean-action attempt is now executed and enclosed:
all recorded gyro prefixes have norm <1, and all 40000 pre-accelerometer
states in its 400--600 s tail have force/field sine >2/5. This still does not
cover all admitted histories. Its correlated cumulative-energy ellipsoid
fails the sufficient pointwise exclusion `E_col-E>0` by a certified margin
in [-817884.048495,-817884.048494]. No collinear source trajectory has been
established. The missing ingredient is the signed, time-linked innovation
recursion, not a new physical restriction or more numerical precision.

The current signed-temporal note proves the carried compatibility criterion
ker C subset ker W, retains both compatibility residual sums, and gives
physical tilt sampling and signed bias/velocity summation bounds. It does
not prove either temporal margin. The former assertion that two positive
margins automatically supply all six historical pivots is withdrawn; that
quantitative bridge and coefficient compactness remain open.
