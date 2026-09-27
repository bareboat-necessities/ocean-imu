# Quiet water is not, by itself, an A21 detectability obstruction

Consider the literal A21 linear comparison at a stationary nominal solution:
zero corrected angular rate, zero wave acceleration, constant attitude and
finite constant OU parameters. Use a fixed sample step h>0 and a recurring
periodic correction word containing S, accelerometer and accepted magnetometer
observations. Measurement noise is positive definite. The magnetic field and
gravity are nonparallel. This proposition concerns that linear comparison,
not every filter trajectory driven by a quiet physical history, startup
capture, or the nonlinear finite-error theorem.

**Proposition.** The lifted stationary pair has no unobservable eigenvalue
on or outside the unit circle, including all nuisance coordinates.

**Proof.** The prediction matrix is block diagonal in (theta,b_g),
(v,p,S,a_w), and b_a. Its only eigenvalues are

`1, exp(-h/tau_aw), exp(-h/5000)`.

The latter two lie strictly inside the unit circle. For a lifted word of N
steps, the only possibly nondecaying eigenvalue remains 1. If
`F^N x=x`, its stable AW and BA components vanish. The attitude block is
`[[I,Nh I],[0,I]]`, so b_g=0. On the neutral translation block,
`F_N(Nh)-I` forces v=p=0. Thus every such eigenvector has the form

`x=(theta,0,0,0,S,0,0)`.

Prediction leaves this vector constant at all events. A zero S observation
forces S=0. The accelerometer and magnetic observations respectively force
`[g_b]_cross theta=0` and `[B_b]_cross theta=0`. Their kernels intersect only
at zero because g_b and B_b are nonparallel. Any lever-arm gyro-bias Jacobian
multiplies the already-zero b_g and cannot change this conclusion. Therefore
no nonzero unit-eigenvalue vector lies in the kernel of the lifted observation
matrix. The unobservable invariant subspace consequently has only the two
strictly stable OU eigenvalues. This proves detectability, including possible
Jordan chains: any nonzero unobservable generalized unit-eigenspace would
contain an unobservable unit-eigenvector. No restricted-information lifting
or independent nuisance-state assumption is used.

This calculation uses the source transition signs and Jacobians:
`Phi_theta,bg=+h I`, `J_acc,theta=-[R(aw-g)]_cross`,
`J_acc,aw=R`, `J_acc,ba=I`, `J_mag,theta=-[R B]_cross`, and `H_S` the S selector.
At the stationary zero-residual solution, attitude injection/reset is identity.
The small-x OU transition branch changes the AW-to-neutral column, but not
the above eigenvector argument because AW=0 at eigenvalue 1.

The active BA decay is essential to this argument. Replacing it by an identity
predictor would permit, for example, a constant tilt theta parallel to B and
a compensating constant accelerometer bias. That is an H18 comparison, not
the shipping A21 predictor. The physical bias is still an independent physical
history; its mismatch with the estimator OU predictor remains a disturbance
in the finite-error proof.

## Consequence for the marine-motion domain

This proposition does not establish a uniform contraction margin over the
current admissible histories. It does establish that zero wave amplitude is
not automatically an unobservable nondecaying mode of A21. There is therefore
no proven quiet-water instability here that would justify excluding it.

The moving proof uses complete excited windows inside a physical moving episode,
with finite-transition obligations at its boundaries (`ou3-regime-design.md`).
It does not impose a positive span on indefinite physical rest. Nor does this
nominal detectability proposition independently identify physical tilt and
accelerometer bias: BA OU decay is an estimator prior, not physical evidence.
The stationary sensor map has an attitude/BA ambiguity even with known magnetic
field. Direct gyro-bias information has a bounded-noise/rate supply. A full
stationary practical theorem for all compatible physical histories remains open;
its deterministic ambiguity must not disappear through nominal covariance decay.

The exact smooth rest/motion/rest construction has identical IMU and magnetic
packets throughout, with unchanged physical bias bounds and applied magnetic
service. Any detector's entire internal history is identical on both truths.
No finite dwell or estimated-wave-state threshold repairs that identifiability
obstruction. No stationary estimator change is enabled by this proposition.

## Uniform historical action on the exact quiet nominal subcase

This continues the same historical reader of `ou3-ag-readout-proof.md`.
It enters `V_next<=rho V+c_d|d|^2` through the full covariance upper
comparison/coercivity and the existing homogeneous process-loss implication.
It supplies neither the nonlinear disturbance coefficient nor a retained radius.

Use the already specified regular A21 zero-residual nominal record: identity
attitude, corrected gyro rate/AW/means zero, fixed reference B e_x with
g,B>=9, no frame/relock or reference change, and actually applied coincident
acc/mag groups every eight qualified 4--6 ms predictions. All mean projections
are then inactive and literal attitude resets are identity. This is a nominal
subfamily, not a conclusion about every stationary or quiet-compatible input.
The bounded default process profile and existing regular S service are retained.
The inherited auxiliary AW and BA standard-deviation ceilings are 156 and
1/40, even when subsequent observations are omitted. Use analytic ceilings
sigma_acc<=1/3, sigma_mag<=1, q_theta<=2e-6 and q_bg<=1e-9. These cover the
stated .00135/.8/default-bias profile without changing any noise setting.
The existing positive process lower bounds are retained separately.

At each anchor select actual acc-y, negative acc-x and mag-y rows and divide
by g,g,B. Their auxiliary observation is `z_i=theta_i+eta_i`, with

`eta_i=((aw_y+ba_y+n_a,y)/g, -(aw_x+ba_x+n_a,x)/g, n_m,y/B)`.

Here all symbols are auxiliary linear errors/noises, not physical mean
observations. They express the literal sensor Jacobians. Holding the actual
coefficients and sync additions fixed gives `Cov(eta_i)<=D`, where

`D=diag(((156+1/40+1/3)/9)^2, ((156+1/40+1/3)/9)^2, 1/81)`.

Covariance Cauchy--Schwarz retains AW/BA correlations; fresh observation noise
is independent only in the auxiliary covariance realization, never by a new
physical-noise assumption. The nuisance means are not replaced by boxes.
Unconditioned auxiliary AW stays below its inherited variance ceiling under
the actual frozen sync increments and stable OU recursion. BA obeys its own
inherited bound. Neutral v,p,S have no column in these selected rows.

For two anchors separated by T, use the terminal trial reader

`theta_read=z_1`, `bg_read=(z_1-z_0)/T`.

Its six AG root columns cancel exactly for arbitrary AG covariance and
arbitrary AG/nuisance cross covariance. This is the existing root-cancellation
condition LO=T_h. The gyro process residual
`xi=theta_1-theta_0-T bg_1` has covariance

`Cov(xi)<=(q_theta T+q_bg T^3/3) I`.

This follows by composing the full correlated source AG process blocks.
At zero rate the Simpson integral is exact, so nonuniform prediction steps
compose to the same continuous-time integral. Positive-definite real-arithmetic
LDL cleanup adds no floor. The two gyro/bias coordinates have not been split
within a step. This fresh AG process residual is independent of the
auxiliary root nuisance and AW/BA/measurement factors. The two eta epochs
can be fully correlated: `Cov(eta_1-eta_0)<=4D`.

With Tmin=4/125 and Tmax=6/125, the trial error and optimal-conditioning
comparison therefore give the uniform historical action ceiling

`P_hh<=B_q:=2 diag(D, 4D/Tmin^2+(q_theta/Tmin+q_bg Tmax/3) I)`.

The factor two retains the attitude/gyro cross covariance by a full Loewner
comparison. Every actual intermediate correction remains in the optimal
covariance; zero reader weights do not omit those corrections. This bound is
independent of the unknown AG root and of the length of the preceding rest.

For a target at delay d<=dmax=Tmax after the latest complete group, put
`B_q=diag(D_theta,D_bg)`. Ignoring additional optimal corrections is conservative.
The actual zero-rate prediction and full correlated Q give principal ceilings

`A_theta=D_theta+dmax^2 D_bg+(q_theta dmax+q_bg dmax^3/3) I`,
`A_bg=D_bg+q_bg dmax I`.

Thus `P_hh<=B_prefix:=2 diag(A_theta,A_bg)` at **every operation** after the
second complete group. After the existing 17-second nuisance comparison is
available, arbitrary AG/nuisance cross covariance is retained through

`P<=C_q:=2 diag(B_prefix,U_n)`.

`stationary_covariance.py` reproduces these rational matrices. The same
construction applies on the identical-input hidden-motion witness because it
has the same nominal execution, not because its physical bias follows the OU
prior. A supplied 80-digit correlated-root audit has normalized action ratio
about .4139295 at AG root scales 1 and 10^12, with zero root residual; exact
rational audits independently check cancellation and matrix dominance. These
finite calculations audit algebra and are not the proof of uniformity above.

### What linear loss now follows, and what does not

On this subcase the retained source has a positive uniform full process floor,
bounded F and the just-proved C_q. Hence there exists epsilon>0 satisfying the
**full matrix** inequality `Q>=epsilon F C_q F'`. The already established
covariance-energy identity gives homogeneous contraction
`V_next<=V/(1+epsilon)` at the next prediction; intervening optimal corrections
and the identity resets cannot increase that homogeneous storage. This is a
qualitative linear-comparison result within the existing proof. No sampled
spectrum or fitted rho is used. The coarse ceilings do not certify a useful
numerical margin, a nonlinear radius, capture, or a stationary physical theorem.

The physical quiet tilt/BA ambiguity persists, and the deterministic physical
mismatch supplies must still be charged. Uniform covariance of a nominal
comparison is not shrinking physical ambiguity. General quiet-compatible
histories, all moving words, nonzero injections/reference variation, finite
transition retention and float32 totality remain outside this subcase result.

## Startup gate check

Zero input from reset leaves the `WavePeriodEstimator` proxy states and
variances zero. Its variance guards therefore prevent `hasUsablePeriod()`
from becoming true, so the ordinary tuner-ready handoff cannot be inferred
from elapsed time alone. This does not prove permanent startup failure:
`maybeHandOffToMekf_()` has a separate `ready_by_timeout` branch requiring
`proxy_ready`, the timeout, and `mag_gravity_aligned_branch_`, without
requiring tuner readiness. Finite capture must prove these actual predicates
and the subsequent magnetic refinement/release. Neither the wave-period gate
alone nor the stationary tail proposition settles that capture obligation.
