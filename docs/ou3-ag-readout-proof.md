# Historical six-column action for the corrected AG loss

Regime scope: `ou3-regime-design.md` separates stationary observability and
finite transitions from complete excited moving windows on this same carried
proof path. `ou3-moving-six-pivots.md` supplies a conditional all-six-pivot
bound with chronological reset and actual-row defects. Its uniform premises
remain OPEN; neither this reduction nor the mode split certifies B_*, J_AG,
rho_0, a nonlinear retained radius or every-prefix retention.

The exact zero-residual quiet nominal subcase is now instantiated in
`ou3-stationary-detectability.md`: the same reader cancels arbitrary AG root
and cross covariance, yielding a uniform action, every-operation full upper
comparison and qualitative homogeneous linear loss on that restricted class.
General stationary physical robustness and the moving source-uniform bound
remain OPEN. Reset inverses and C=A^-1 B give a complementary chronological
transport budget without deleting any reset.
Actual acc/mag groups in one prediction cell factor with E=0. A conditional
geometry bound uses the magnetic direction pulled back through their resets;
its source-uniform positive premise remains OPEN. In world coordinates the
array O is attitude-free (`ou3-world-frame-rows.md`); same-cell geometry
depends on magnetic cadence, while aggregate rows need only the nominal
signed AW mean for their attitude columns; Theorem G0 there bounds the
injection-free array.


This result enters the existing finite-error inequality through a source-uniform
`J > 0`, then a matrix bound on `rho0`, in
`sqrt(V_next) <= sqrt(rho0) sqrt(V_root) + supply`.
It is a conditional construction within the single A21 proof. The existing
nuisance lower/upper comparisons and six-column Schur implication are retained.
No source-uniform J, useful decay rate, or retained region is certified here.

## Why future observations alone do not give an absolute J

The complete corrected energy identity necessarily gives

`0 <= D_W,hh <= (P_root^-1)_hh = (P_hh-P_hn P_nn^-1 P_nh)^-1`.

In particular, for `P_root=diag(t I6,I15)`, every future word satisfies
`D_W,hh <= I6/t`, regardless of its realized gains, resets or observations.
The nuisance upper/lower comparisons can remain fixed as t increases. At
`t=10^12`, a candidate `J=10^-6 I6` fails: every unit AG Rayleigh margin is
at most `10^-12-10^-6=-9.99999e-7`. This is an exact obstruction to removing
the AG prior from a *forward-only* absolute-loss argument. It is not a
shipping-reachable covariance family, an admitted magnetic-service execution,
or a stability counterexample. In particular, magnetic service may restrict
this family; no service claim is made for it.

The existing implication is valid. Its premise must use the inherited
history, not just positive future excitation. A positive future normalized
loss can coexist with an arbitrarily small absolute root loss.

## Backward readout action

Freeze the actual regular A21 coefficients on a historical window, retaining
each prediction F, process factor U (`Q=U U'`), applied observation H with
effective noise factor V (`R_eff=V V'`), and literal reset G. The window
starts after the established 17-s nuisance upper comparison applies. PSD
sync increments are additional process factors with identity transition.
Safety increments belong to the effective observation noise. A mean-only
projection has no covariance operation; its finite-error defect stays in
the separate nonlinear proof. Release/frame events require their actual maps
and are not silently omitted by this regular-window construction.

Realize the frozen covariance recursion as an auxiliary linear estimation
problem. Independent unit-covariance inputs are a matrix-algebra device;
they impose no stochastic law on physical marine or bias histories. At an
observation the latent state has identity transition. Optimal conditioning
produces exactly the actual Joseph covariance and realized gain; afterward
the actual G acts by congruence. This does not run a second adaptive filter.

Let O stack the six AG columns of its **raw auxiliary** observation maps,
including all preceding predictions and resets. Let T_h be the AG root
columns of its terminal AG map. Choose a block reader L with

`L O = T_h`.

These raw rows are not the physical vector Gramian or the corrected loss
rows. They are used only to construct a trial estimator of the terminal AG
state. The optimal-estimation comparison below accounts for every applied
correction, including corrections whose reader block is zero.

The action of that reader can be evaluated backward using only a 6x21
residual Y and a 6x6 action B. Initialize `Y=E_h'`, `B=0`. In reverse
operation order:

| Operation | Action addition | Residual predecessor |
|---|---|---|
| Prediction | `(Y U)(Y U)'` | `Y F` |
| Applied observation i | `(L_i V_i)(L_i V_i)'` | `Y-L_i H_i` |
| Literal reset | zero | `Y G_i` |

Require exact cancellation `Y E_h=0` at the historical root. Write the
remaining columns as `Y_n`. With the already established full nuisance
upper comparison `P_nn <= U_n`, set

`B_W = B + Y_n U_n Y_n'`.

**Lemma.** The actual terminal AG covariance satisfies `P_hh <= B_W`.

**Proof.** Expand the error of the trial estimate `sum_i L_i y_i` of the
terminal AG state. The recursion gives its coefficient on every root and
noise input. Its AG root coefficient is exactly zero, so both the unknown
AG covariance and all root AG/nuisance cross terms cancel. The remaining
root covariance contribution is at most `Y_n U_n Y_n'`. Each fresh input
contributes the corresponding full factor action. In particular, different
coordinates driven by one process factor are not split into independent
inputs. Optimal linear estimation has covariance no larger than this trial
estimator, which proves the claim. All gains and resets belong to the same
frozen actual execution. The inequality is between covariance matrices, not
between two nonlinear estimator trajectories.

The implementation selects rows of O by exact largest-residual factor
pivoting and solves a 6x6 system exactly. It does not invert a 21-state covariance or an observation normal
matrix. For a uniform proof the reader may depend on the realized coefficients;
one fixed reader is generally invalid under varying coefficients. Even a
nonzero root residual of size `10^-60` is rejected without a separate AG
root bound: its action can grow without bound as the root uncertainty grows.

## Joint minimum-action reader (exact)

The pivot reader is one feasible reader; the best one needs no selected
minor. Stack every action source of the trial estimator as one unit-covariance
column: the nuisance-root factor `Gamma` (`Gamma Gamma'=U_n`), each prediction
or sync factor `U_j` in operation order, and each applied noise factor `V_i`.
The frozen word is then exactly (`augmented_design` in `ag_readout.py`)

`y = O_h h0 + A s`,  terminal AG state `= T_h h0 + T s`,

with `h0` the unknown AG root and `s` the source vector. Every reader with
`L O_h = T_h` has action `B(L)=(T-L A)(T-L A)'`.

**Theorem (joint reader).** Let `Sigma=A A'>0` and let
`I_eff = O_h' Sigma^-1 O_h`.

1. **Minimizer.** If `I_eff>0`, the unique Loewner-minimal feasible reader is
   `L* = T A' Sigma^-1 + Tt I_eff^-1 O_h' Sigma^-1`, where
   `Tt = T_h - T A' Sigma^-1 O_h`.
2. **Minimum action.** Its action is `B* = Pi + Tt I_eff^-1 Tt'`, where
   `Pi = T (I - A' Sigma^-1 A) T'`.
3. **Every feasible reader.** For each `L` with `L O_h = T_h`,
   `B(L) = B* + (L-L*) Sigma (L-L*)'`.
4. **Diffuse limit.** For `h0 ~ N(0, t I)` the auxiliary posterior covariance
   of the terminal AG state is `Pi + Tt (I/t + I_eff)^-1 Tt'`. This equals the
   frozen-coefficient Riccati recursion started at `diag(t I6, U_n)`, and it
   increases to `B*` as `t` grows.
5. **Coercivity reduction.** For every `g>0`,
   `B* <= (1+1/g) T T' + (1+g) T_h I_eff^-1 T_h'`.
   Here `T T'` is the source-driven terminal AG covariance. Because shipping
   predictions have no nuisance-to-AG block, `T T'` is the AG process Gramian
   of the window.

*Proof.* Write `L=L*+Delta` with `Delta O_h=0`. Then
`(T-L* A)A' = -Tt I_eff^-1 O_h'`, so every cross term vanishes, which gives
item 3. `I - A' Sigma^-1 A` is the orthogonal projector onto `ker A`, which
gives item 2. Woodbury gives item 4; the Gaussian model is realized by the
Riccati recursion with realized coefficients frozen. For item 5,
`Sigma^-1/2 O_h I_eff^-1 O_h' Sigma^-1/2` is a projector, hence
`T A' Sigma^-1 O_h I_eff^-1 O_h' Sigma^-1 A T' <= T T'`; apply Young to `Tt`. ∎

Consequences:

- **Prior-independent AG bound.** `P_hh <= B*` holds for every actual root with
  `P_nn <= U_n`, retaining all root cross covariance. No AG prior, pivot or
  threshold enters.
- **Full 21×21 bound.** Riccati monotonicity gives `P_end <= lim_t Ric_W(diag(t I6, U_n'))`
  for every `U_n'` with `P_nn < U_n'`. This full-matrix upper bound needs no
  Young split and no `U_n` inflation of the terminal nuisance block.
- **Remaining premise.** The only data-dependent quantity is `I_eff`: the AG
  root information after marginalizing the nuisance root and every
  process/measurement source. A source-uniform `B_*` follows from a
  source-uniform floor `I_eff >= mu` together with the explicit process
  Gramian and a bound on `|T_h|`.
- **Role.** `B_*` and the full bound serve coercivity (`P <= C` in the
  nonlinear supplies). They are not needed for `rho_0`: on carried MOVING
  words the information-ratio bound with this `C` is attained at `k = 0`,
  and the word diameter of `ou3-corrected-word-proof.md` section 7 needs at
  most the scalar kernel variance `nu' P_0 nu`.
- **Measurement-only special case.** The earlier all-row noise reader
  (`minimum_noise_reader`) is the case with no process or nuisance sources.
- **Verification.** `joint_reader_audit` checks items 1–5 exactly on the
  supplied 21-state word; the pivot reader's action is strictly larger there.

**Carried feasibility (non-promoting).** The table below uses real-arithmetic
optimal-gain replays of the literal carried coefficients (corrected OU core)
and windows ending at 225.32 s. Each ratio is the largest generalized
eigenvalue against the replayed actual covariance.

| Window | Quiet `B*`/`P_hh` | Wave `B*`/`P_hh` | Quiet full/`P` | Wave full/`P` |
|---:|---:|---:|---:|---:|
| 0.32 s | 6.3e8 | 4.9e8 | 9.5e8 | 5.6e8 |
| 4 s | 1827 | 2560 | 1844 | 2569 |
| 16 s | 12.2 | 12.4 | 76.6 | 75.0 |
| 64 s | 5.07 | 5.40 | 75.4 | 55.1 |

On 0.32-s words `B*` has `lambda_max` 778 (quiet) and 193 (wave). The pivot
reader's committed values are 1199 and 1450. Short words cannot identify the
gyro-bias root, so no reader is commensurate there. From 16 s on, the joint
reader bounds the actual AG covariance within a factor of about 12, and
within about 5 at 64 s. The full bound stays within 55–77. These replays do
not certify `I_eff` source-uniformly.

## Conditional bootstrap to J and full contraction

If the historical action is uniformly bounded by a single `B_* > 0`, then
at its terminal roots block Cauchy--Schwarz gives, for every eta>0,

`P <= C_eta = diag((1+eta) B_*, (1+1/eta) U_n)`.

This preserves arbitrary cross covariance through the comparison, without
claiming that the actual covariance is block diagonal. If, uniformly at the
next prediction,

`Q - epsilon F C_eta F' >= 0`, `epsilon>0`,

then the complete corrected loss satisfies

`D_W >= delta P^-1`, `delta=epsilon/(1+epsilon)`,

and in particular

`D_W,hh >= J := delta B_*^-1/(1+eta) > 0`.

This J satisfies the retained six-column sufficient premise. The direct
historical upper comparison can be sharper than reconstructing it through
the small nuisance ratio alpha; the existing Schur implication remains
valid and available. The matrix comparison, rather than a scalar Q floor,
is verified with exact PSD elimination. Later corrections/resets cannot
undo this homogeneous loss. For a useful margin, however, retain and enclose
the entire actual word factor; the first-prediction implication is not
automatically a practical nonlinear margin.

**This route has a uniform ceiling.** `first_prediction_relative_ceiling` in
`corrected_word.py` proves it for every upper comparison `C >= P_root`,
whatever `B_*`, `U_n` or `eta`. Test the premise on one axis's S coordinate
(the displacement integral) with `y=e_S`:

`epsilon <= Q_SS/(F L F')_SS <= 3.3741e-10`

- **Numerator.** The OU triple-integrator impulse response in S is at most
  `t^3/6`. Hence one step has `Q_SS <= Sigma_aw (h/tau) h^6/126 (1+eps_q)`.
- **Denominator.** `L` is the certified pre-prediction LIN floor.

So `N` predictions certify at most `N epsilon`. That is at most `1.73e-4` per
2048-s proof word, against a carried word contraction of about `0.39` per
64 s. On carried 0.32-s words the ideal ratio with the literal root
covariance is `2e-23` (quiet) and `4e-22` (wave), and the best structured
chain from `(B*, U_n)` is about `1e-36`.

The implication itself stays valid algebra. The contraction is recast in
`ou3-corrected-word-proof.md` §6 as one information-ratio inequality and in
§7 as a word Riccati diameter, where no root covariance matrix enters.

## Executed checks and unresolved source inequality

The 80-digit supplied 21-state experiment includes varying transitions,
rank-three sensor updates, cross covariance, correlated process columns,
and nonorthogonal literal-form attitude resets. Its weakest absolute AG
loss decreases from `0.0281547113907` to `9.99999999965e-13` as the AG root
scale goes from one to `10^12`. The reader's maximum action eigenvalue is
`34.4659867861`; the actual terminal covariance stays below that same matrix.
The exact rational checker verifies root cancellation and matrix dominance,
and a conditional first-prediction decrement `1/10000000001`.
These are algebra audits; their rational LIN/process coefficients are not
shipping discretization, and their supplied roots are not reachability evidence.

The outstanding source inequality is a *uniform historical action bound*

`for every admitted realized window W: L(W) O(W)=T_h(W), B_W <= B_* < infinity`.

Neither source-uniform full AG rank nor a common action ceiling has been
established for these nominal coefficients and injection transports. The
physical 3-D Gramian cannot be substituted for O; the actual nominal AW,
gyro estimate and resets still need their same-history comparison.


## Executed source audit and the uniform-certificate failure

`ag_readout_source_diagnostic.py` compiles a read-only observer and an untapped
control from the shipping header. It carries construction, startup, reference
refinement, release, tuning, covariance and means to 225 s without a reseed.
The observed window ends at 225.32 s. Each control has exactly the same terminal
state, quaternion, covariance, stage times and accepted magnetic count.

The 80-digit results are:

| Input | Applied rows | Minimum singular value of full raw AG array | Maximum action eigenvalue |
|---|---:|---:|---:|
| Quiet, body field heading 0 | 225 | 7.1517206581 | 1199.0609668 |
| Quiet, heading .001 rad | 225 | 7.1517206581 | 1199.0622214 |
| Quiet, heading .000001 rad | 225 | 7.1517206581 | 1199.0609668 |
| Moving vessel | 228 | 6.2782239972 | 1450.2528611 |

For the moving record, `p_z=.4 sin(.6t)`, `roll=.02 sin(.5t)` and
`B_world=(60,0,30)`, with zero physical biases. These continuous truth formulas
obey the displacement/primitive, rate, acceleration and jerk limits; the
finite float samples and applied events do not certify all-time magnetic
service or arithmetic totality. The largest realized injection in this
window is about 3.65717e-6 rad. The diagnostic keeps actual time-varying
coefficients, accepted acc/S/mag observations, resets and sync increments.
It uses the established nuisance upper comparison, including all correlations.

An exact audit of the exported quiet coefficients follows the diagnostic.
Binary floating operands are interpreted as rational inputs, not as an
enclosure of the underlying real-arithmetic shipping trajectory. For each
correlated Q/R block, exact semidefinite LDL elimination gives `Q=L D L'`;
upward rational square roots of D give `U U' >= Q` by congruence. The reader
cancels the AG root **exactly**, and exact PSD elimination verifies a full
6x6 action ceiling with its off-diagonal entries retained. Nonzero source
sync increments are included. This supplied-sequence certificate does not
cover a neighborhood, all windows, startup reachability in exact arithmetic,
or accumulated float32 error.
The moving export has a raw sync skew of exactly 2^-49. The complete
operation is now observed before inflation and after the shipping addition
and full covariance symmetry pass. Write its literal boundaries as P_b,P_a,
and form E=P_a-(P_b+P_b')/2 exactly. This is not generally PSD: at event 96,
the normalized (e_0+e_2) direction has Rayleigh quotient -2^-44. The negative
quantity is a rounding defect, not evidence that P_a is indefinite.

Signed rank-one/two Schur elimination gives E=sum d_i v_i v_i'. Retaining
only positive terms and rounding their square roots upward gives a rational
factor U with U U' >= E, verified by exact full matrix PSD elimination.
This establishes P_a <= sym(P_b)+U U' for **each recorded complete operation**.
The antisymmetric part of the pre-state and the final averaging roundoff are
therefore not silently discarded. All 21 coordinates are retained; symmetry
roundoff is not confined to the three AW coordinates. Zero-diagonal,
nonzero-off-diagonal defects use two-column pivots, not a diagonal norm bound.

Both the quiet and moving exported words now have exact rational action
ceilings, exact AG-root cancellation, and exact comparisons with their
literal terminal AG covariance. This continues the factor calculation past
the old operand-symmetry rejection. It does not enclose prediction, Joseph,
solve, reset or state-update arithmetic, or provide a common all-history
factor ceiling. Those separate operations remain open. The exact action is
computed without rounding the reader or substituting a finite replay for
source-uniform coverage.

### Why changing the minor is necessary but insufficient

A first-independent-row rule can select a magnetic component
`B_y` tending to zero while ignoring a strong `B_x`. In a quiet two-epoch
heading/bias calculation its measurement-noise action alone is

`(R_m/B_y^2) [[1, 1/h], [1/h, 2/h^2]]`.

Thus no common finite ceiling exists for that selector as nonzero `B_y -> 0`,
even though the full observation array retains rank. The replacement performs
exact largest-residual factor pivoting across **all** rows, then solves a
six-coordinate system. There is no normal-equation rank threshold, dense
interval Riccati propagation, or scalar/Gershgorin contraction reduction.
The selector alone does not prove its action bounded.

To test uniform feasibility independently of any minor, retain the full
rank-three observation blocks and form

`I_W = sum_i O_i' R_i^-1 O_i`.

For every reader with `L O=T_h`, completing the square gives

`B_W >= L R L' >= T_h I_W^-1 T_h'`.

Equivalently, a proposed common ceiling must satisfy the full Schur condition
`[[B_*, T_h], [T_h', I_W]] >= 0`. If I_W is singular and T_h is invertible,
no exact reader exists. This is a **necessary noise-action test**, not a lift
of restricted information or an upper bound on nuisance/process action.

The attempt to enclose all words using independent nominal coefficient ranges
fails this test exactly. Put

`B=(45,0,45), a_hat=(-g/2,0,g/2), f_hat=(-g/2,0,-g/2), g=9.80665`.

Even `|a_hat|=g/sqrt(2)<8.8`; merely imposing the physical acceleration ceiling
on the nominal acceleration would not fix this relaxation. With zero corrected
rate and identity resets, the literal AG prediction is `[I,h I;0,I]` and both
sensor attitude blocks are cross products with parallel vectors. The complete
raw AG array has rank four. The two exact null columns are `(z,0)` and `(0,z)`,
where `z=(1,0,1)`. The endpoint map does not kill either column. Hence
`L O=T_h` is impossible, and the Rayleigh margin for any positive proposed
Gram floor `mu I6` is exactly `-mu` (reported at mu=1 as -1).
Nuisance columns and process correlations cannot repair this root cancellation.

This is a failure of the **independent-coefficient relaxation**, not a
shipping counterexample: the constant nominal AW values have not been linked
to the actual OU mean, innovations, pseudo-observations and resets. No all-time
MAGNETIC SERVICE or physical-history membership is asserted for that family.
The refinement from a fragile minor to all-row factors fixes the selector;
the all-row nullspace test then rules out further pivot, precision or interval
refinement as a cure for the relaxed domain. The next technique must use
same-history nominal dynamics to exclude sustained near-collinearity and
control the full signed temporal gyro margin before enclosing the residual
matrix action. The implemented bias ball supplies the separate one-prediction
gyro-transport floor below; true-vector sampling fidelity alone does not
establish the remaining historical reader.

In particular the literal accelerometer relation is

`f_hat = f_measured - b_hat_a - r_acc`

in this default, zero-lever-arm, reference-temperature profile. An attempted
transfer of the true-vector Gram bound must retain the actual innovation
`r_acc`, the nominal gyro-bias error and reset transport. There is no certified
all-window bound on those defects here. Replacing them with the physical
sensor/bias bounds would silently identify estimator innovations with sensor
noise. This is the precise remaining source-uniform certificate gap.

### Historical unprojected alias and implemented exclusion

This test enters the controlling tail inequality through the still-required
uniform historical AG action. Innovation bounds were the next candidate for
excluding the previous rank-four coefficient relaxation. They are insufficient
alone, even with nonparallel nominal force and field and zero innovations.

In real arithmetic take h=1/200 s, quiet truth, zero physical gyro bias,
measured gyro zero, nominal gyro bias -400 pi e_z, and zero nominal
v,p,S,a_w,b_a. Let the nominal attitude initially agree with truth and let
the committed field be (75,0,0). Use the source gravity constant in the quiet
accelerometer sample (any representational discrepancy fits the existing
sensor residual bound). The literal large-angle quaternion branch makes one
complete turn at each prediction: its quaternion changes sign, but its
rotation matrix returns to the same value. Every applied acc/S/mag innovation
is zero, for any realized finite gain. The nuisance zero mean is preserved
under arbitrary admitted OU coefficients in the unprojected relaxation.
The shipping gyro projection changes this bias immediately, so this is no
longer compatible with the implemented nominal mean recursion.

The literal Rodrigues/integral helper gives

R(h)=I,  B(h)=h e_z e_z',  F_AG=[[I,B(h)],[0,I]].

Consequently the initial gyro-bias columns e_bg,x and e_bg,y never reach an
attitude or nuisance column. Every observation annihilates them. Corrections
with any realized gain also leave them unchanged, since (I-KH)v=v when Hv=0.
Resets are identity. The full raw AG array has rank four, while T_h preserves
both missing columns. Hence LO=T_h is impossible and the Gram-floor margin
against mu I6 is exactly -mu. Correlated process factors cannot repair an
exact deterministic root-column annihilator.

For theta=|omega_hat|h, the two transverse singular values of B(h) are
2|sin(theta/2)|/|omega_hat|. An 80-digit diagnostic at theta=2pi+delta gives
7.9564805098e-7, 7.9577458881e-10 and 7.9577471533e-13 s for delta=10^-3,
10^-6 and 10^-9 respectively. The exact complete-turn nullspace, rather than
the numerical near-null values, certifies this failure.

**The historical relaxation is outside the implemented estimator family.**
OU-II/III now enforce `|b_hat_g|<=0.5 rad/s` as residual-bias protection;
this is separate from the unchanged physical `|b_g|<=0.02` qualification.
[The source-bound certificate](ou-gyro-bias-projection.md) charges physical
rate, calibrated residual, fast measurement residual and estimated bias.
On the qualified 4--6 ms family, the prediction angle is below .007 rad and
every transverse real-source transport singular value is at least
`.003999991833333333 s`. Both Rodrigues and small-rate polynomial branches
are included. This removes the exact and near-complete-turn bias relaxation.
It does not close the full signed temporal margin or the historical action:
actual chronological observations, resets and mean-projection defects still
need joint source-uniform bounds. Nominal force/field collinearity remains.
The device API's unrestricted positive timesteps are outside this timing
certificate. No all-history float32 totality is inferred.

### Construction and carried-entry attempt

The construction starts the MEKF means at zero. Before Live, `updateFrontEnd`
drives the proxy/tuner but leaves those means unchanged, and magnetic
corrections are withheld from the MEKF. `goLive` writes the proxy attitude and
its covariance; it does not import a proxy gyro-bias estimate. The initial
MEKF gyro bias is therefore zero, with marginal covariance 10^-6 I3. Any
magnetic corrections before the first prediction have zero gyro gain because
the handoff clears attitude/gyro cross covariance. With the existing physical
rate, bias and sensor envelopes, the first corrected angular increment is
at most .006*(.6108652381980153+.02+.02)=.0039051914291880918 rad.
This tighter construction bound applies at the first prediction. After the
subsequent coupled acc/S/mag corrections, the implemented .5 rad/s bias ball
provides the uniform qualified angle bound below .007 rad described above.

The attempted propagation must retain the true initial translation and
primitive, the front end, clocks, references and tuning, rather than replace
the reached root by a convenient mean/covariance box. The new reproducible
`construction_history_diagnostic.py` executes that full path for

p(t)=-(3/2)sin(2t)(1,0,1), R(t)=I, B=(45,0,45),

with zero true biases. Its continuous squared envelopes are p:9/2, v:18,
a:72, jerk:288, angular rate:0 and primitive span:9/2. These satisfy the older pointwise/jerk limits and have no displacement or
acceleration DC, but fail the current MARINE MOTION tilt-span condition.
This is an old-domain diagnostic, not an admitted moving history.
The native driver uses only `begin`, `update` and `updateMag`, with no reseed,
and retains the shipping front end, reference refinement and BA release.

It reaches Live at step 30002 and refinement/release at step 36008. On the
400--600 s diagnostic tail, mean tilt error is 8.211611 degrees and maximum
is 8.260681 degrees. The proposed six-degree entry comparison has negative
mean margin -2.211611 degrees over that finite interval. The largest nominal
acceleration from construction is 9.776391 m/s^2, illustrating why the old-domain proof could not substitute the physical
8.8 bound for a nominal coefficient bound. The current moving domain still
requires its own construction-linked derivation. The
largest gyro estimate is .000161120 rad/s and the tail minimum normalized
force/field cross magnitude is .4432715: this particular trajectory does not
approach either raw-rank degeneracy.

These are **finite float diagnostics**, not a proof of every-window magnetic
service, a real-arithmetic reachable-set enclosure, or a refutation of
history-dependent eventual capture. The actual six-degree entry set, the
uniform joint nominal-history bound, and the common B_* remain unproved.
This construction attempt does not justify moving a free-root relaxation
into the admitted set or promoting a sampled positive margin to a theorem.

The finite endpoint admits a stronger **exact** storage check. Interpret the
recorded float coefficients as rationals and verify the entire 21x21 P is
SPD by exact elimination. For the known zero true BA, the covariance Schur
identity gives, with all cross covariance retained,

V=e'P^-1e >= e_ba' P_ba,ba^-1 e_ba
  =41601169495079937572567/3619798568741048312 >11492.6752.

Thus the candidate V<=36 has margin at most
-41470856746605259833335/3619798568741048312 <-11456.6752 at that recorded
endpoint. This is an exact failure to lie in the sufficient projection guard,
not a claim that projection must be active, nor an all-time capture or
magnetic-service counterexample. It needs only a three-coordinate inverse
after checking full covariance positivity. The validator reproduces this
certificate from the committed endpoint and binds its generated native driver.

The construction-linked joint mean action is now executed and enclosed in
`ou3-construction-mean-action.md`. It proves finite recorded gyro and force
separation but its proposed energy-only uniform exclusion fails by an exact
margin in [-817884.035305,-817884.035304]. This is a deficient enclosure,
not a new physical counterexample. The signed whole-window source relations
and the common historical action ceiling remain unproved.
