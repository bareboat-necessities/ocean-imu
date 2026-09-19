# OU-III proof research state

## Current hypothesis

The productive proof route does not require H18 asymptotic contraction. The
shipping wrapper holds accelerometer-bias learning while its magnetic reference
is provisional and releases that external hold after reference refinement. The
internal gate then clears after the accepted-update threshold and one-second
guard. H18 is therefore a finite bridge; the recurring asymptotic tail is A21.

The single route is
`construction -> finite capture -> finite H18 bridge -> finite reference
refinement/release -> recurring magnetically informed A21 -> regional practical
stability`.

## Negative result retained from PR #557

Strict full-state H18 incremental contraction in the covariance metric is
obstructed by the held accelerometer-bias identity block on a feasible held
segment. This invalidates that proof tactic only. The finite-difference export,
candidate-rho search, and associated measurement scripts were research
diagnostics, not proof machinery, and have been removed.

## Productive lemmas

`tail_stability.py` contains theorem algebra only. A finite H18 recurrence
`V+ <= g V + d` stays bounded for any finite number of bridge steps, including
`g >= 1`. Under MAGNETIC SERVICE, every service window contains at least one
actually applied informative correction; once magnetic-reference refinement
finishes at finite time, the remaining accepted-update count and one-second
guard therefore clear in finite time.

For the recurring A21 tail, if the source-uniform linear word satisfies
`||F e||_W <= sqrt(rho0)||e||_W` and the finite nonlinear remainder difference
has storage-norm Lipschitz gain `eta`, then the complete finite-error word has

`rho = (sqrt(rho0) + eta)^2`.

Thus the decisive small-gain condition is `sqrt(rho0)+eta < 1`. This is the
route to a real finite-error theorem; no sampled trajectory or candidate
numerical rho is used.

## Current limiter

The fixed-coordinate long-word information certificate is constructive:
`mu_N >= 2.04e-3`. The controlling obligation is covariance normalization.

The literal 21-state entrywise midpoint-radius Riccati enclosure has now been
executed after its single permitted near-identity refinement. It failed by
interval dependency: predicted covariance spectral enclosure
`[-18.7907040215, 2502.41441854]`; S innovation `r/a=0.66666386445`
verified, while accelerometer and magnetometer innovation boxes did not.
Pointwise shipping covariance remains PSD.

Accordingly the next covariance representation must preserve PSD structure.
No further entrywise subdivision is permitted by the research protocol. Until
a recurring positive root covariance floor `p_min` is certified,
`mu_cov/rho0`, explicit nonlinear `r_*`, arithmetic practical radius,
capture/release-to-tail, and every-prefix tail retention remain open.

Independently retained: finite magnetic-reference refinement follows once
capture enters the <=7 degree tilt domain; H18 remains a finite bridge; the
literal accelerometer-bias projection sector is dissipative and projection-
inactive release is not required.

## Failed approaches / DEAD_ENDS

Do not reconstruct the PR #557 H18 full-state contraction search. H18 need not
be asymptotically contractive because it is a finite pre-tail bridge. Do not use
finite differences, a single simulated superword, or covariance-shaped point
ratios as theorem evidence. Do not narrow the physical motion/bias assumptions
to make a numerical contraction easier.

## Next proof step

Formalize the magnetic-reference refinement state machine and prove finite
release from recurring accepted informative corrections. In parallel, derive
the A21 linear error-transition energy inequality from the literal shipping
operations and the MAGNETIC SERVICE Gramian. The nonlinear small-gain remainder
is then bounded on the retained finite-error domain against that analytic
linear margin.


## A21 structural observability advance

For one scalar active-bias OU-III translation chain with state
`(v,p,S,a_w,b_a)`, exact constant-parameter OU transition, four successive
integral observations `H_S F^j`, `j=0..3`, and one attitude-normalized
accelerometer row `H_a=(0,0,0,1,1)`, a five-row observability minor has the
closed-form determinant

`det O = -Delta^3 tau^3 (1-exp(-Delta/tau))^3`.

It is nonzero for every finite `Delta,tau>0`. This removes a possible
structural rank obstruction: the integral chain plus active accelerometer-bias
state is observable per axis; the bias does not create an unobservable
translation mode after release.

This is not yet the shipping full-information certificate. The deployed tuner
can change `tau`, integral pseudo-updates occur on their actual scheduler, and
attitude is coupled rather than known. The next analytic step is to extend this
minor to bounded time-varying shipping transitions and combine its three-axis
floor with gravity/attitude and MAGNETIC SERVICE. That extension, rather than a
sampled matrix rank calculation, is the current route to the complete A21
information floor.


## Uniform time-varying A21 closure

The time-varying tuner/scheduler no longer requires freezing tau. The shipping
pseudo cadence is progress-preserving and satisfies
`T_S(tau)=clip((0.015/1.1)tau,0.005,0.15)`. On the admitted
`dt in [0.004,0.006]`, `tau in [0.02,12]` envelope, the accumulated OU
decay exponent between successive S updates is bounded by 0.55, including one
sample of scheduler overshoot. Across three gaps the latent-acceleration kernel
therefore retains at least `exp(-1.65)>0.19`.

For arbitrary positive time-varying `lambda(t)=1/tau(t)`, the S response to
initial acceleration has a kernel `K` with
`K'''(t)=exp(-integral lambda)>0`. Hence
`{1,t,t^2/2,K(t)}` is an extended complete Chebyshev system. Four successive
S observations have determinant equal to a Vandermonde factor times
`K'''(xi)/12`; the 4 ms sample-spacing floor makes this uniformly positive.
This proves time-varying translational observability for all three axes without
enumerating tuner histories.

The accelerometer attitude Jacobian has two nonzero singular values equal to
`|a_w-g|`, uniformly at least `9.80665-8.8=1.00665`. On the neutral
quotient (active accelerometer bias removed because its A21 OU predictor is
strictly stable), two separated gravity observations expose the two tilt/gyro
bias pairs. Known body rotation is orthogonal transport. The remaining
heading/axial-gyro-bias pair is exactly the subspace covered by MAGNETIC
SERVICE. Thus the non-decaying quotient is uniformly observable. With strict
inner guard margins, the finite hybrid branch cells form a compact union, so
the normalized quotient Gramian has an **existential uniform floor**

`mu_N := min_h lambda_min(J_N(h)) > 0`.

This is a proof of positivity, not yet a useful numerical enclosure of mu_N.

The A21 prediction is uniformly completely controllable on the same compact
envelope: the integrated-OU process covariance is SPD for positive dt, tau and
sigma; the exact attitude/gyro-bias covariance has positive gyro white-noise
and bias-RW densities; and active accelerometer bias has positive OU driving
density. The default periodic a_w covariance synchronization only queues a
bounded PSD increment, while bias release only raises finite P_ba diagonal
entries. These events preserve covariance compactness.

Uniform quotient observability/detectability plus uniform controllability gives
uniform upper/lower Riccati bounds and therefore a source-uniform linear A21
covariance-metric ratio `rho0<1`. No sampled rho is used.

The earlier attempt to require projection-inactive release was too strong and
has been removed. H18 holds the accelerometer-bias estimate while physical
truth may be anywhere in the all-time ball ||b_a||<=0.2251666 m/s^2. Even with
a zero held estimate, universal release error can therefore be 0.2251666
m/s^2, which exceeds the projection-inactive radius 0.1748334 m/s^2. The A21
proof now carries the literal projection sector
||e_ba^+||^2+||d_b||^2<=||e_ba^corr||^2 instead. This is dissipative and gives
the global post-projection bound ||e_ba||<=B_a+R_b=0.6251666 m/s^2; capture no
longer has to manufacture an unjustified 0.15 m/s^2 bias error.

The remaining work is constructive rather than architectural: rigorously
enclose a usable numerical `mu_N`/`rho0`, derive an explicit `r_*`, bound
the additive arithmetic supply, and prove finite capture/H18 release enters and
retains that explicit inner domain.


## Marine motion supplies the missing attitude diversity

A separate attitude-excitation assumption is not required. Let
`f=a-g` be the world specific-force vector and `B` the true geomagnetic
field. On every interval of length `T`,

`integral f x B dt = (v(T)-v(0)) x B - T g x B`.

The marine contract gives `||v||<=5.5 m/s`; the magnetic contract gives
`B_h>=15 uT` and `||B||<=75 uT`. Therefore

`(1/T)||integral f x B dt|| >= g B_h - 2 V_max B_max/T`.

The right side is positive for `T>5.60844 s`. Choosing an 8 s subwindow gives
a uniform diversity floor `43.97475 (m/s^2) uT`. Hence every 8 s physical
marine history contains an instant where accelerometer and magnetic vector
sensitivities are non-collinear. MAGNETIC SERVICE supplies an applied magnetic
observation in every 1 s interval; transporting that sensitivity through the
known attitude transition to the diversity instant preserves its norm. The two
rank-two vector observations are therefore jointly full rank for attitude.
Two consecutive 8 s diversity windows expose gyro bias through its attitude
injection. The A21 proof word is consequently taken as 16 s.

This closes the earlier gap in the attitude/gyro-bias quotient argument using
existing MARINE MOTION and MAGNETIC SERVICE assumptions, rather than adding a
persistent-excitation assumption. The 16 s word also improves the active
accelerometer-bias homogeneous norm factor to `exp(-16/5000)=0.996805...`
(`rho_b=exp(-32/5000)<0.994`).


## Quantitative enclosure status

A first fully analytic normalized translation enclosure is now constructive.
On the 16 s word, selecting S updates after times 0, 8 and 16 s with the
literal maximum scheduler delay 0.156 s, state scales
(V,p,S)=(5.5,8.1,1100), and worst declared S-noise standard deviation 100,
the determinant/Frobenius certificate gives

`mu_trans >= 2.04734e-3`.

This is a real lower bound, not a sampled singular value.

A deliberately sparse two-epoch attitude/gyro calculation gives a much weaker
candidate scale (~6.18e-7) and therefore is **not promoted as the final
shipping mu_N certificate**. The reason is important: MAGNETIC SERVICE is an
already-transported two-coordinate heading/axial-bias Gramian over a window,
not an instantaneous pure-heading row. A tight full certificate must compose
that actual 2-D service Gramian directly with the many accelerometer rows in
the same 16 s word; replacing it by a fictitious instantaneous attitude
measurement would be an invalid shortcut. The next quantitative proof step is
therefore the aggregate Schur/Gramian bound on the literal partition, using all
recurring accelerometer and S information rather than two sparse rows.

The float32 arithmetic path is likewise separated correctly. A straight-line
kernel with n rounded operations has the standard gamma_n bound
`gamma_n=n*u/(1-n*u)`, u=2^-24. This is implemented as a certificate
primitive, but no whole-word arithmetic supply is claimed until literal kernel
operation counts and magnitude envelopes are composed. Roundoff remains
additive supply and is not allowed to consume the nonlinear derivative margin.


## Constructive long-word neutral information

The weak gyro-bias scale is handled by accumulation, not by strengthening the
physical assumptions. The proof word is now 2048 s. At every disjoint one-second
MAGNETIC SERVICE root, the service inequality J_hb>=I permits extraction of one
unit of **pure root-heading information** (subtract diag(1,0); the remainder is
PSD). The simultaneous root accelerometer sample has no preceding gyro-bias
transport. With the commissioned detector-band residual bound 0.3 m/s^2,
shipping vibration gain 0.75, and nominal accelerometer std 0.2 m/s^2, the
literal effective accelerometer std is bounded by

`sigma_acc,eff <= hypot(0.2,0.75*0.3)=0.3010399 m/s^2`.

Together with `|f_z|>=g-A_max=1.00665`, `|f_h|<=8.8`, and one unit of
heading service, each service root supplies a pure attitude information floor
greater than 0.0129 in the declared proof coordinates.

Let `y_k=theta_0+C_k b_g` be those extracted root-attitude observations.
Bounded angular rate gives

`sigma_min(C_(k+1)-C_k) >=
  s_bg * 2 sin(Omega_max/2)/Omega_max`

with `s_bg=0.02`, so the normalized increment is at least about 0.01969.
The path-graph inequality
`sum |y_(k+1)-y_k|^2 <= 4 sum |y_k|^2`, combined with the root observation
`y_0=theta_0`, gives the joint attitude/gyro-bias floor

`mu_ag = j*c/(1+c),  c=(N-1) beta^2/4`.

For N=2048 this exceeds 2.1e-3. The independent 16-s S certificate remains
`mu_trans>=2.04734e-3`; therefore the complete **fixed-coordinate neutral
Gramian** satisfies

`mu_N >= 2.04e-3`.

This construction never integrates attitude over 2048 s: it uses only
one-second transport increments, so bounded rotations cannot cancel the
gyro-bias information.

### Important normalization correction

The number above is not yet a covariance-metric contraction factor. The
information-form identity `rho0<=1/(1+mu)` requires the Gramian to be whitened
by the actual root covariance. If, in the same proof coordinates,
`P_root>=p_min I`, then

`mu_cov >= p_min mu_N`,
`rho0 <= 1/(1+p_min mu_N)`.

Using the fixed-coordinate 2.04e-3 directly in the latter formula would be a
coordinate-dependent and invalid shortcut. The next quantitative certificate
is therefore a recurring lower enclosure `p_min>0` for the shipping A21 root
covariance. Only after that enclosure is established will rho0, the explicit
L2/r* radius, and the arithmetic/practical-radius bounds be promoted.


## Riccati enclosure route

The covariance-normalization step now has two layers. First, the corrected
classical UCO/UCC bounds are implemented as a fail-closed cross-check. For
`alpha1 I<=O<=alpha2 I` and `beta1 I<=C<=beta2 I` over N discrete steps,

`P_max <= (alpha1 + N alpha2^2 beta2)/alpha1^2`,

`P_min >= beta1^2/(beta1 + N alpha2 beta2^2)`.

The corrected form is important: the simpler historical upper bound
`O^{-1}+C` is not generally valid. These formulas prove a positive recurring
covariance floor once all four Gramian constants are enclosed, but on the raw
IMU-step horizon they are expected to be too conservative to leave a useful
nonlinear margin.

Therefore they are a verification fallback, not the primary numerical route.
The primary route is a **verified interval Riccati enclosure** of the bounded
shipping schedule, using midpoint-radius arithmetic and a residual/Krawczyk
inverse check for each 3x3 innovation solve. A floating-point midpoint trajectory
may propose the box; only outward-rounded residual inclusion can certify it.
The enclosure must cover the full dt/tau/R_S/R_acc/service ranges and the
literal covariance sync/release events. Its resulting p_min/p_max are then
cross-checked against the classical UCO/UCC bounds before rho0 is promoted.


## Verified interval Riccati certificate implementation

The covariance-normalization blocker now has executable fail-closed certificate
machinery. A candidate covariance trajectory may propose a spectral box
`p_min I <= P <= p_max I`, but it cannot promote it.

For every accelerometer, integral, and magnetic innovation interval, write
`S=S_0+E`. The certificate requires
`lambda_min(S_0)>=a`, `||E||_2<=r<a`. Then the midpoint inverse residual is
bounded by `r/a<1`, and Neumann/Krawczyk inclusion gives
`||S^{-1}||<=1/(a-r)`. If this strict inequality fails, the certificate fails
closed.

The complete interval Riccati image must then satisfy, with outward-rounding
slack included,

`[R(P,parameters)] subseteq [p_min I,p_max I]`

for every covariance in the proposed box and every admitted shipping parameter
box. The certificate also requires explicit inclusion of the periodic
`a_w` covariance synchronization, H18-to-A21 bias-release covariance floor,
the complete bounded tuner/scheduler envelope, and float32 rounding. Only after
all flags and inclusions verify does the machinery compute

`mu_cov=p_min*mu_N`, `rho0=1/(1+mu_cov)`.

A committed status artifact currently remains OPEN. This is deliberate: the
residual/inclusion checker is now present, but the full outward-rounded shipping
Riccati image has not yet been generated. No midpoint-only p_min is accepted.


## Interval Riccati enclosure attempt

The first full 21-state interval enclosure exposed a dependency failure rather
than a covariance instability. The original attitude prediction constructor
used the valid but useless entrywise box R_step in [-1,1]. At the shipping
bound Omega <= 0.610865 rad/s and dt <= 0.006 s, the actual step rotation is
within 0.003666 rad of identity, so that box discards almost all structure
before the first Riccati multiplication.

Classification: **interval-conditioning failure**, not a mathematical
certificate failure. The fixed-coordinate information floor, UCO/UCC result,
literal 21-state covariance/Joseph maps, and hard-event semantics remain valid.

One mathematically motivated refinement has been applied, as permitted by the
research protocol: Rodrigues bounds now retain
`|R_ii-1| <= 1-cos(Omega dt)` and
`|R_ij| <= sin(Omega dt)+1-cos(Omega dt)`, while the gyro-bias injection
integral keeps its near `-dt I` diagonal and O(Omega dt^2) off-diagonal
structure. At the declared ceiling this reduces the attitude diagonal interval
radius from 1 to below 7e-6 and off-diagonal radius below 0.0037.

The finite branch-subdivision machinery is also in place. It exactly covers the
unsplit parameter/H box, verifies every cell independently, and forms the hull
of every branch image. It may be used only once to test whether the tightened
near-identity transition crosses the innovation/self-inclusion threshold.

**Next falsifiable experiment:** run the recurring-box certificate with the
Rodrigues-tight F interval and one finite subdivision of the dominant
accelerometer/magnetic attitude-Jacobian entries. Record the first failed
innovation residual `r/a`, or the recurring spectral inclusion margins. If
this single refinement does not cross the threshold, stop interval subdivision
and review the covariance-normalization architecture rather than subdividing
again.


## Interval Riccati enclosure attempt

The first full 21-state interval enclosure exposed a dependency failure rather than a covariance instability. The original attitude prediction constructor used the valid but useless entrywise box R_step in [-1,1]. At the shipping bound Omega <= 0.610865 rad/s and dt <= 0.006 s, the actual step rotation is within 0.003666 rad of identity, so that box discards almost all structure before the first Riccati multiplication.

Classification: **interval-conditioning failure**, not a mathematical certificate failure. The fixed-coordinate information floor, UCO/UCC result, literal 21-state covariance/Joseph maps, and hard-event semantics remain valid.

One mathematically motivated refinement has been applied, as permitted by the research protocol: Rodrigues bounds now retain |R_ii-1| <= 1-cos(Omega dt) and |R_ij| <= sin(Omega dt)+1-cos(Omega dt), while the gyro-bias injection integral keeps its near -dt I diagonal and O(Omega dt^2) off-diagonal structure. At the declared ceiling this reduces the attitude diagonal interval radius from 1 to below 7e-6 and off-diagonal radius below 0.0037.

The finite branch-subdivision machinery is also in place. It exactly covers the unsplit parameter/H box, verifies every cell independently, and forms the hull of every branch image. It may be used only once to test whether the tightened near-identity transition crosses the innovation/self-inclusion threshold.

**Next falsifiable experiment:** run the recurring-box certificate with the Rodrigues-tight F interval and one finite subdivision of the dominant accelerometer/magnetic attitude-Jacobian entries. Record the first failed innovation residual r/a, or the recurring spectral inclusion margins. If this single refinement does not cross the threshold, stop interval subdivision and review the covariance-normalization architecture rather than subdividing again.


## Refined interval-Riccati experiment result

The required non-promoting feasibility experiment was run in CI on the literal 21-state interval predictor after the Rodrigues near-identity refinement. It did **not** cross the certification threshold. From the constructor covariance seed, the predicted interval covariance had spectral enclosure `[-18.7907040215, 2502.41441854]`. The integral-S innovation remained verifiable with residual ratio `r/a = 0.66666386445 < 1`, but the accelerometer and magnetometer innovation boxes had no strict inverse certificate (`r/a = infinity`).

Classification: **second interval-dependency/conditioning failure of the entrywise midpoint-radius covariance mechanism**. This is not evidence that the shipping covariance is indefinite: every pointwise Riccati covariance remains PSD. It invalidates using further blind entrywise subdivision as the route to `p_min`. The retained facts are the fixed-coordinate `mu_N >= 2.04e-3`, source-uniform UCO/UCC, literal shipping covariance/Joseph semantics, and the verified 3x3 inverse criterion.

Per AGENTS.md, no second round of tighter scalar boxes or deeper subdivision is allowed. The covariance-normalization architecture must now preserve PSD structure during enclosure. The next falsifiable construction is a structure-preserving spectral/square-root Riccati enclosure: represent the covariance set by a positive midpoint/factor plus an operator-norm radius, propagate prediction as a factorized PSD sum, and certify measurement updates through information/Joseph equivalence with verified inverse residuals. It must reproduce the pointwise shipping map and produce a positive recurring spectral margin before it can replace the failed entrywise box.


## PSD-preserving covariance feasibility result

The architecture-review probe succeeded at one shipping sample without manufacturing a negative covariance eigenvalue. With the same constructor covariance seed and source-uniform F/H/R bounds, the PSD-preserving spectral representation certified `sigma_min(F) >= 0.64408710396`, a post-prediction covariance floor `4.14848197494e-7`, and post-correction floors `4.14825252226e-7` (S), `3.01357720533e-7` (accelerometer), and `1.49402788135e-7` (magnetometer). This contrasts with the failed entrywise covariance box on the identical source ranges.

This is a feasibility result, not recurring `p_min`: it currently omits a constructive full-rank lower floor for the literal process covariance Q and therefore cannot be iterated over the 2048-s word without collapsing conservatively. The next falsifiable step is to derive a structure-preserving lower factor/floor for the exact shipping `Q_AA`, integrated-OU `Q_LL`, and active-bias OU `Q_BA`, then iterate the PSD spectral recurrence with the actual correction cadence. If that recurring floor is numerically too small to leave a usable nonlinear margin, covariance normalization must use a block/factor metric rather than the global scalar eigenvalue.


## Multi-step process-controllability probe

The first structure-preserving process-noise feasibility construction formed explicit integrated-OU impulse columns over a 0.156-s scheduler-scale window, but certified them with an unscaled Gershgorin lower bound on the resulting 4x4 Gramian. That lower bound was zero at tau = 0.02, 0.2, 2, and 12 s. Separately, iterating the PSD spectral recurrence with q_floor=0 drove a 1e-6 root floor to 1.46713043e-19 after only 32 maximal-correction samples.

Classification: **conditioning failure of the scalar/Gershgorin controllability reduction**, not loss of controllability. The exact integrated-OU Gramian is SPD for positive window, tau and process density, and uniform controllability was already proved structurally. The failed hypothesis is that an unscaled scalar eigenvalue/Gershgorin reduction can preserve enough of that structure for covariance normalization.

This is the first implementation of the new PSD-preserving mechanism, so one mathematically motivated refinement remains available. The quantitative reason is explicit: the state units differ by powers of time, making the raw Gramian strongly scaled. The next falsifiable experiment is therefore a **proof-coordinate-scaled factor/LDLT enclosure** using the already declared (v,p,S,a) scales rather than another scalar norm. It must certify positive pivots of the multi-step controllability Gramian uniformly over tau/dt; otherwise the PSD-preserving covariance-normalization route must be reviewed again.


## Block/factor covariance metric

The covariance-normalization architecture has moved to a block/factor metric after the second conditioning failure of the entrywise midpoint-radius Riccati box. The state is partitioned as AG=(attitude,gyro bias), LIN=(v,p,S,a_w), and BA=(active accelerometer bias). Each diagonal covariance block is represented by a certified factor floor P_ii >= L_i L_i^T, while cross-block covariance is carried by operator-norm bounds.

After scaling by the factor floors, block Gershgorin/Schur coercivity requires gamma=1-max_i sum_(j!=i) ||P_ij||/(ell_i ell_j)>0. Only then can fixed-coordinate information be covariance-normalized. For block information floors mu_i, the certificate uses mu_cov >= gamma*min_i(mu_i ell_i^2), followed by rho0<=1/(1+mu_cov). This prevents a poorly conditioned global scalar covariance eigenvalue from destroying useful block information while still accounting for every cross-covariance.

The old entrywise interval Riccati representation is retired and must not be subdivided again. A non-promoting scaled integrated-OU factor probe is now in CI. The next source-uniform certificate must enclose tau/process-noise variation and the within-slice controllability remainder, then establish AG and BA factor floors and recurring cross-block operator bounds under the literal shipping correction/event cadence.


## Block/factor coercivity refinement

The first block/factor implementation attempted to recover coercivity from
absolute cross-block bounds through
`1-max_i sum_j ||P_ij||/(ell_i ell_j)`.  The universal PSD cross bound is far
too weak because the process-factor floors are small while covariance ceilings
carry physical units and long-integrator scales.  This is a conditioning
failure of that sufficient Gershgorin reduction, not a covariance failure.

The literal shipping prediction supplies a stronger identity that does not
discard cross structure:
[
 P^- = F P^+ F^T + Q,qquad
 Q=operatorname{diag}(Q_{AG},Q_{LIN},Q_{BA}).
]
If source-uniform process factors satisfy
`Q_i>=L_i L_i^T`, then exactly
[
 P^- - operatorname{diag}(L_iL_i^T)
 = F P^+ F^T + (Q-operatorname{diag}(L_iL_i^T)) >=0.
]
Thus the post-prediction block metric has coercivity **gamma=1 independent of
the carried cross covariance**.  Cross-block operator bounds remain required
for finite-prefix magnitude/totality, but they do not need to establish root
coercivity.

For a correction, information form gives the exact sufficient recurrence
`gamma^+=gamma^-/(1+gamma^- j_D)`, where `j_D` is an upper bound on the
measurement information in the block-factor coordinates.  Finite measurement
noise therefore preserves a strictly positive prefix metric margin.  The next
controlling numeric task is consequently the source-uniform factorization of
the three literal process blocks, especially Q_LIN; AG and BA have closed-form
process factors.


## Block-factor source-range feasibility result

The source-range diagnostic was executed in CI. It gives an exact scaled AG one-second process factor ell_AG=4.999974644e-4; the active BA recurrence gives ell_BA=5.618273739e-4. The 16-s scaled LIN controllability probe is positive at every tested tau, with its smallest sampled pivot 1.194309698e-7 at tau=0.02 s. These values support the factor architecture but the LIN number is not promoted until the continuous tau interval is enclosed.

The first cross-block mechanism, the universal PSD ceiling bound, decisively fails the requested scalar block-Gershgorin test. Its normalized cross bounds were (4.092217619e8, 2.013744150e6, 1.030074464e8), giving gamma=-5.122292072e8. Classification: **conditioning failure of absolute cross-ceiling / scalar-factor Gershgorin**, not failure of covariance coercivity. The invalidated hypothesis is that generic PSD cross ceilings can be divided by the small process factor floors to certify the root metric.

The mathematically justified block/factor refinement uses the additive shipping process identity instead. At a post-prediction root, P-=F P+ F^T+Q and Q is block diagonal. Certified process factors Q_i>=L_iL_i^T imply exactly P->=diag(L_iL_i^T) regardless of the carried cross covariance, i.e. root metric gamma=1. Corrections preserve a positive factor-metric margin through gamma+=gamma-/(1+gamma- j_D) for finite normalized information ceiling j_D. Absolute cross-block bounds remain a prefix/totality obligation, but the failed generic cross-Gershgorin reduction is no longer used to establish root coercivity.

The remaining quantitative issue is LIN: the useful 16-s controllability factor must be propagated through the intervening literal corrections, or an equivalent verified 4x4 LIN Riccati factor must be enclosed. Using the one-sample scalar eigenfloor would be positive but numerically useless and is not promoted.


## LIN correction-survival falsification

The generic batch endpoint comparison has been retired fail-closed.  An
accumulated process Gramian floor together with per-event information ceilings
does not by itself survive corrections interleaved with process injection.  A
two-dimensional exact counterexample is

`Q1=[[4,2],[2,5/4]], Q2=diag(100,1/100), H=[1,0], R=1`.

Although `Q1+Q2 >= I`, the exact correction between `Q1` and `Q2` gives
endpoint covariance `[[504/5,2/5],[2/5,23/50]]`, whose second diagonal is
`23/50 < 1/2`.  Thus the former endpoint formula would overclaim the floor.

This does not invalidate the additive post-prediction identity, the
single-correction factor-metric inequality, UCO/UCC, or the fixed-coordinate
`mu_N >= 2.04e-3`.  It narrows the controlling LIN obligation: prove an
effective endpoint-information comparison through the literal shipping
correction/synchronization sequence, or propagate a structure-preserving LIN
factor recurrence event by event.  The theorem remains fail-closed until that
shipping comparison is numerical and source-uniform.
