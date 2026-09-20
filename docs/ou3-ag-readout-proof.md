# Historical six-column action for the corrected AG loss

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
The moving export fails the symmetry precondition for a process factor:
its largest sync skew entry is exactly 1/562949953421312 (2^-49). The
diagnostic records this rejection. Enclosing the literal addition and final
symmetry operation, rather than silently replacing its operand, remains open.

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
sampled gyro aliasing quantitatively, before enclosing the residual matrix
action. True-vector sampling fidelity alone does not establish this exclusion.

In particular the literal accelerometer relation is

`f_hat = f_measured - b_hat_a - r_acc`

in this default, zero-lever-arm, reference-temperature profile. An attempted
transfer of the true-vector Gram bound must retain the actual innovation
`r_acc`, the nominal gyro-bias error and reset transport. There is no certified
all-window bound on those defects here. Replacing them with the physical
sensor/bias bounds would silently identify estimator innovations with sensor
noise. This is the precise remaining source-uniform certificate gap.
