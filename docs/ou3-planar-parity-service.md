# Planar covariance parity and service reduction

Status: analytical real-arithmetic factorization. The forward-invariant interval
cell and all-time service lower bound remain OPEN until their numerical
containment certificates are produced.

For the exact planar MOVING record, all nominal rotations are about body Y and
all delivered accelerometer/magnetometer vectors lie in the XZ plane. In the
shipping 21-state order

    (theta,b_g,v,p,S,a_w,b_a),

define the two coordinate sets

    E={theta_y,bg_y, v_x,v_z,p_x,p_z,S_x,S_z,aw_x,aw_z,ba_x,ba_z},
    O={theta_x,theta_z,bg_x,bg_z, v_y,p_y,S_y,aw_y,ba_y}.

They have dimensions 12 and 9.

Prediction preserves E/O: the exact attitude/BG transition for omega parallel
to e_y mixes x with z and leaves y alone; the LIN OU chain is axis-separable;
the BA OU factor is scalar by axis. Process Q has the same parity.

For any XZ predicted vector f=(f_x,0,f_z), -[f]x has the pattern

       [ 0  f_z  0]
       [-f_z 0  f_x]
       [ 0 -f_x  0].

Hence the X/Z accelerometer or magnetic residual rows couple only to theta_y
(and, for acceleration, X/Z AW/BA), while the Y row couples only to theta_x,z
(and Y AW/BA). The S=0 selector is axis-separable. Therefore H P H'+R, K,
the Joseph correction and the error-reset Jacobian for a Y-only injection all
preserve the E/O covariance block diagonal structure. The startup/relock
covariances are diagonal and hence start in this invariant set.

The four service probes are [h+,bg+,h-,bg-], with

    d_+ = Ry(-psi)(0,+sin(theta),cos(theta)),
    d_- = Ry(-psi)(0,-sin(theta),cos(theta)).

The orthogonal +/- transform splits the auxiliary four-column action into

    (d_+-d_-)/sqrt(2)  in E,
    (d_++d_-)/sqrt(2)  in O,

with the corresponding BG probes. This is useful for propagating covariance and
homogeneous columns without cross-parity dependency. The literal MAGNETIC
SERVICE criterion, however, is not replaced by a 4x4 eigenvalue: shipping checks
the two physical 2x2 heading/BG Gramians for d_+ and d_- separately. Each such
physical pair is reconstructed from its E and O components and must itself have
lambda_min >= mu_M. No cross-sign 4x4 criterion is promoted.

A useful exact geometry identity is also phase independent. With
B_b=75 Ry(-psi)e_x,

    ||[-B_b]x d_+|| = ||[-B_b]x d_-|| = 75,

because the common Ry rotation cancels in the cross product and d_+/- are unit
vectors orthogonal to e_x. This identity is retained symbolically before any
intervalization.

## Continuous hard-iron symmetry

On the exact planar record the continuous hard-iron tracker cannot create a spurious accepted 3-D offset. Its excitation matrix is M=I-A' A with A the weighted mean of tilt rotations. Every rotation is about Y, so A e_y=e_y and M e_y=0. Therefore lambda_min(M)=0 and the literal information gate in ContinuousMagHardIronEstimator::solve_ fails. The continuous estimate is not promoted and the applied hard-iron correction remains its zero startup value on the exact zero-hard-iron history. The corrected magnetic input consequently retains norm 75 uT and the learned planar reference remains finite/nonzero.

## All-time acceptance

After Live/refinement, the wrapper forwards every supplied 25-Hz magnetic call
to the MEKF. The MEKF magnetic path has no NIS rejection: for finite nonzero
magnetic input it rejects only if the innovation factorization fails. In exact
real arithmetic, P>=0 and Rmag=0.8^2 I imply

    S_m=H_m P H_m' + Rmag >= 0.64 I,

so the factorization is SPD and every call is accepted. Bounded tuner clamps,
finite process matrices and Joseph corrections preserve finite PSD covariance
at every finite prefix. This proves the actual accepted cadence on the planar
real-arithmetic execution once the reference/refinement stage has completed.
It does not by itself prove the required information floor.

## Literal accelerometer nuisance is a Schur complement

For one linearized accelerometer correction partition the error as target x=(attitude,BG) and nuisance n=(AW,BA and the other coordinates entering the row). With prior covariance [[Pxx,Pxn],[Pnx,Pnn]] and row y=A x+B n+v, exact Gaussian elimination gives n|x covariance N=Pnn-Pnx Pxx^-1 Pxn, effective target row A_eff=A+B Pnx Pxx^-1, and effective noise R_eff=R+B N B'. The target information increment is exactly

    J_x^+ - J_x^- = A_eff' R_eff^-1 A_eff.

This is the same update implemented by the literal PCt/S/Joseph chronology; no independent nuisance gain is introduced. The identity is regression-tested in accel_conditional_information.py. It establishes the correct comparison object for the direct-attitude oracle.

However, the desired service comparison is not a free monotonicity theorem. Stronger nonmag information shrinks target covariance but also changes the transported homogeneous heading/BG probes before the next magnetic sample. Therefore samplewise I_acc,literal <= I_acc,oracle alone does not yet imply an ordering of the later magnetic service Gramian. The remaining bridge must propagate the paired covariance+probe action (or prove an equivalent information-form variational inequality) through prediction, S, accel and reset on the same history. This prevents an invalid shortcut.

## Paired covariance--probe correction identity

For S=HPH'+R, K=PH'S^-1, Joseph Pplus and Phi_plus=(I-KH)Phi,

    Phi_plus' Pplus^-1 Phi_plus = Phi' P^-1 Phi - (H Phi)' S^-1 (H Phi).

Every accelerometer, S and magnetic correction consumes its PSD action. Reset
congruence preserves paired storage; prediction generally consumes additional
storage. A product of prediction-retention factors times a replay service
minimum is not a magnetic-service lower bound. All intervening correction
losses must be carried. The old equality without the loss was erroneous.

## Attitude reset metric

The shipping first-order attitude covariance/probe reset uses G=I+(1/2)[dtheta]x on the attitude block. Because [dtheta]x is real skew-symmetric, G'G=I+(1/4)(||dtheta||^2 I-dtheta dtheta'). Its singular values are exactly 1 and sqrt(1+||dtheta||^2/4) (twice). Thus the reset is invertible, has sigma_min=1, and its inverse has sigma_min=1/sqrt(1+||dtheta||^2/4). This closes the coordinate-conditioning formula needed by the comparison; a uniform numerical dtheta bound must still come from the same-history correction cell rather than an independent clamp.

## Scheduler coordinate is structurally invariant

The scheduler component of the augmented cell is not itself open: the literal retarget function preserves elapsed time when it is below the new period and otherwise parks it immediately below the new period; periodic_update_due then returns an elapsed value in [0,T_S). Hence the causal set {(elapsed,T_S): T_S>0, 0<=elapsed<T_S} is forward invariant under every shipping retarget/due operation. What remains open is covariance/service enclosure uniformly over that scheduler coordinate, not boundedness of the scheduler state.

## Exact-rational oracle margin

A deliberately more-informative two-state attitude/BG feasibility oracle has now been solved exactly over rational arithmetic. It uses dt=1/200 s, direct isotropic attitude information H_a=30 with R_a=1/25 at every IMU sample, magnetic H_m=75 with R_m=16/25 every eight samples, P_theta,0=1/2000 and P_bg,0=11/10000000, and the literal gyro/process densities used by the replay. Over one second the exact rational Riccati map satisfies P_0-P_1 positive definite, so this oracle cell is self-containing by Riccati monotonicity. The accumulated magnetic heading/BG information satisfies I_mag-I_2 positive definite exactly; det(I_mag-I_2) is about 0.98970433. Thus the magnetic geometry has strict margin even after intentionally strong direct attitude conditioning.

This is not yet a shipping service proof. Promotion requires a formal Schur/conditional-information comparison showing that the literal accelerometer+S nuisance chronology cannot suppress the magnetic incremental information more than this oracle, an all-time bound corresponding to H_a<=30 on the literal predicted specific force/attitude row, and treatment of the shipping reset coordinates. These dependencies are explicit in planar_oracle_rational.py and fail closed.

## Literal prediction and scheduler reuse

No independent reimplementation of the 12x12 LIN prediction is needed. The proof instrumentation already receives F_LL and Q_LL after the shipping calls to IntegratedOUChain::transition/process_covariance, together with F_AA,Q_AA, the BA phi/Q block and every correction/reset in operation order. The factor stream therefore treats those callback matrices as the literal coefficients. Likewise S=0 is observed only when the shipping time_update calls periodic_update_due and then applyIntegralZeroPseudoMeas; the factor stream consumes that actual correction event instead of approximating a cadence. This removes the surrogate F_LL/Q_LL and scheduler chronology from the intended certificate. A temporary non-promoting diagnostic exports the first-200 literal prediction retention and S-event count to verify plumbing before interval promotion.

## Complete phase-aware causal cell: OPEN

`planar_moving_center.py` implements a conditional moving covariance center
on the complete observed 40000->48000 word, inherited without reseeding at
44000. It recomputes each gain from the center's linked P,H,R, and retains
the exact observed F/Q, resets, AW targets and S/AW event placements.
**FINITE DIAGNOSTIC ONLY:** relative Frobenius center defects are .0003693294
at 44000 and .0006158214 at 48000; the maximum over all 8000 sample prefixes
is .0006179245. The corresponding fixed-initial-center drifts are .5378221
and .4473738. BA-y variance grows by factors 1.13510 and 1.26914. The AW ages
at these endpoints are 20 and 9 samples, with S elapsed times .0878590 and
.0428655 seconds. They are coordinates from one observed chronology, not
independent phase boxes or evidence of a 20-second return.

This reference still consumes the observed future H and reset G. Therefore
its small defect cannot be inserted as a uniform q_P or used to build an
autonomous all-time center. The calculation isolates finite transport defect
from center drift; the next missing object is a linked future coefficient and
arithmetic enclosure, plus the nonlinear mean/physical-gauge chart. The
physical arc and its quadratic bias remainder are derived in M2a of
`ou3-moving-quiet-compatibility.md`; its physical angle cap is not the required
precision-normalized amplitude. No proof count or admission status changes.

The same diagnostic evaluates the **curved** physical arc in each exported
actual-P metric. Over the finite 8000-prefix tail, the largest central-physical-
chart gauge amplitude is 16.6316254 (sample 40001), and the largest transverse
curvature is .184909934 (sample 47992). Both arc endpoints are tested as a
linked one-parameter family. The transverse extremum formula is justified
because t-sin(t) and 1-cos(t) increase for 0<=t<=theta; maximizing over the
sign of beta makes the mixed term nonnegative. A positive gauge-coordinate
derivative comparison validates its endpoint extrema. These floating results
are **FINITE DIAGNOSTIC ONLY**, and the central physical chart is not the
actual nominal-error chart. Neither supplies an all-time gauge/curvature bound.

The exact covariance decomposition is 12+9, but these blocks are NOT independent
causal histories. A forward-invariant cell must retain the planar mean, both
covariances, private raw Mahony quaternion/integral, guard, frequency/variance,
staged/applied joint tuner tuple, S phase, AW synchronization clock and pending
target, reference/refinement state and gates. The S interval alone is invariant;
that is not covariance or service containment.

Default AW synchronization is P -> P+Ew(Sigma-Pww)_+Ew'. It is not Loewner
monotone. `planar_service_cell.py` proves its fixed-target Frobenius bound and
the rank-4/even, rank-2/odd linked S/prediction commutator. The old scalar-norm
bound from P<=Pupper and the orientation-independent rank-three accelerometer
ceiling were false; corrected bounds and exact counterexamples are documented
in the appendix. These are proof-code failures, not shipping counterexamples.

`planar_service_stream.py` exports the true pre-prediction P, literal F/Q/R_S,
all actual correction H/R/S/K/PCt/residuals, pre/post reset, AW target/pending
state and sample state in operation order. `planar_service_audit.py` checks this
word before measuring local defects. Its S/prediction commutator is not the
complete scheduler-cell difference, which must also traverse acc/reset and
possibly mag/AW operations and future mean-dependent coefficients.

The entire covariance partial derivative, retaining all AW replacements, has
finite twenty-second relative Frobenius gains approximately 0.74261434 and
0.87382130. This is a frozen-coefficient Jacobian block in different root/end
metrics, not a self-containing shipping cell. The full same-history mean and
coefficient feedback remains open. Products of one-step/sync-local norms lose
this cancellation and are retired after the failed grouping refinement.

The literal private float Mahony observer is not exactly normalized. Its raw
quaternion is used to produce vertical acceleration. The old interval reference
also omitted gravity subtraction and lacks literal seeding/period/tuner/clock
binding; it must not feed a shipping HistoryCell. The default entry now fails
closed with an implementation-binding error. A conditional pitch/integral
quadratic tube has an exact positive margin, but its every-step error and
initialization caps still need source binding.

One all-time component has closed: `planar_service_guard.py` proves that the
exact planar record's seeded guard stays inactive in real arithmetic, since its
detector RMS is below 0.000531<0.03. Float32 transfer is not inferred from this.

For service, the orthogonal plus/minus transform gives parity information
matrices I_E and I_O with each physical pair I_+=I_-=(I_E+I_O)/2. Therefore the
required lower floor is for those physical 2x2 matrices, not the minimum of the
parity eigenvalues and not a 4x4 eigenvalue. Every placed window/root and all
future scheduler/adaptation phases must be covered after complete causal
self-inclusion. Current all-time admission and exclusion remain OPEN.

Structures preserved: full 21-state literal execution; parity is a lossless
permutation; actual gains, Joseph, resets, AW synchronization and both clocks;
physical bias/S histories are inherited and never reset by proof boundaries.

Relaxations introduced: fixed-coefficient covariance derivatives and local
covariance cells are explicitly subordinate blocks/outer sets, not reachable
shipping trajectories. Numerical spectra and reference oracles are finite
feasibility diagnostics only. No physical assumption, estimator parameter or
quality gate has changed.

## Quotient chart, metrics and complete native secants

Status: **FINITE DIAGNOSTIC ONLY** for the numerical word. The physical line is
fixed before any SVD. For R_bw=Rx(alpha)Ry(psi), the shipping left W->B error has
dtheta=-Ry(-psi)e_x d(alpha), whereas db_a=g e_y d(alpha). Orienting the line
with positive roll therefore gives r=(Ry(-psi)e_x,0,...,-g e_y). The positive
BA-y handoff seed must not be used unchanged in this chart. At the central
physical family the literal acc and mag rows annihilate this first-order line;
this does not identify a nonlinear fibre around the finite +/- pair.

For P=L L', whiten by L^-1 and normalize u=L^-1 r/||L^-1 r||. Check uu' and
I-uu' at the original dimensionless tolerance. This avoids the false absolute
precision-scale failure in Pi'J-JPi; it does not relax SPD. If U spans u's
orthogonal complement, T=L_N^-1 M L_0 gives
M_Q=U_N' T U_0, C_Q=U_N' T u_0, and b_Q=U_N' L_N^-1 b. Retain
xi_N=M_Q xi_0+b_Q+C_Q alpha_0 and a physical bound on alpha.

On 40000 -> 44000 the F/(I-KH)/G homogeneous product has Euclidean gain 122.66,
full covariance-metric gain .93478409 and quotient gain .87311969, with
||C_Q||=.00022901228. The legacy positive-BA seed gives .87160337 and .48225195.
A 60-digit endpoint check agrees; it does not enclose the rounded product.
The homogeneous product omits nonlinear dK*r and prediction/injection/reset
terms, so it is not the complete mean/P Jacobian. The private Mahony gain
.58270576 is also not a MEKF mean gain; its old .05265365 coupling budget is
withdrawn.

The new native secant probe forks the complete inherited execution and runs the
unchanged wrapper for all 21 mean and 231 covariance perturbations, at epsilon
.01 and .005. Endpoint-metric quotient gains are about .87312274; the b_Q
secant changes from .00253837 to .00425304, c_Q is about .0585573, and the full
matrix step-halving difference is .01714895. Float secants are not derivatives
or uniform bounds. The null fork is bitwise identical and the forcing/clock
hashes match for these finite runs.

A fixed-root covariance metric changes the comparison materially: rho_P about
.991 and point drift q_P=.537822 require an optimistic zero-gauge relative
radius above 60. A two-sided relative SPD ball requires radius below one.
This **D_SUFFICIENT_BOUND_FAILURE** rejects that candidate representation,
not the estimator or all possible phase/history-dependent tubes. Both gauge
injections must be retained where present (native root-metric C_P about .00133,
C_Q about .00022527). The scalar radius solver never promotes its algebra to
forward invariance or service. Uniform derivatives, physical gauge amplitude,
future forcing/center drift, exact S/AW phase coverage and every-window
Delta I<6.024764605642485 remain OPEN.

The numerical profile above is the existing planar probe configuration of the
shipping wrapper (sigma_a=.2, adaptive S cadence). The AtomS3R sketch explicitly
uses sigma_a=.12 and fixed S cadence, as well as its own gravity/magnetic-start
settings. No transfer of the carried floor to that distinct profile is proved.
