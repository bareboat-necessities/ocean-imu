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

The correct same-history comparison variable is now explicit. For any positive-noise linear correction with S=HPH'+R, K=PH'S^-1, Joseph covariance P+=(I-KH)P(I-KH)'+KRK', and homogeneous probe Phi+=(I-KH)Phi, the information-form identity P+^-1=P^-1+H'R^-1H and I-KH=P+P^-1 imply exactly

    Phi+' P+^-1 Phi+ = Phi-' P^-1 Phi-.

Thus literal accelerometer and S=0 corrections do not consume covariance-metric heading/BG probe storage when P and Phi are propagated together. This is stronger and cleaner than trying to order covariance or gains separately. The Schur identity above explains how nuisance enters the correction, but no additional nuisance-loss charge is needed for this paired storage.

This does not by itself prove the next magnetic service summand: prediction adds process covariance and therefore can reduce Phi'P^-1Phi, the magnetic row rotates with the nominal attitude, and reset changes coordinates. Those are now the only comparison losses that must be bounded between magnetic events. The reset singular values are given below. The remaining prediction loss should be evaluated directly from the literal Q/F factors on each 5-ms step, preserving the two parity blocks and scheduler phase.

## Attitude reset metric

The shipping first-order attitude covariance/probe reset uses G=I-(1/2)[dtheta]x on the attitude block. Because [dtheta]x is real skew-symmetric, G'G=I+(1/4)(||dtheta||^2 I-dtheta dtheta'). Its singular values are exactly 1 and sqrt(1+||dtheta||^2/4) (twice). Thus the reset is invertible, has sigma_min=1, and its inverse has sigma_min=1/sqrt(1+||dtheta||^2/4). This closes the coordinate-conditioning formula needed by the comparison; a uniform numerical dtheta bound must still come from the same-history correction cell rather than an independent clamp.

## Scheduler coordinate is structurally invariant

The scheduler component of the augmented cell is not itself open: the literal retarget function preserves elapsed time when it is below the new period and otherwise parks it immediately below the new period; periodic_update_due then returns an elapsed value in [0,T_S). Hence the causal set {(elapsed,T_S): T_S>0, 0<=elapsed<T_S} is forward invariant under every shipping retarget/due operation. What remains open is covariance/service enclosure uniformly over that scheduler coordinate, not boundedness of the scheduler state.

## Exact-rational oracle margin

A deliberately more-informative two-state attitude/BG feasibility oracle has now been solved exactly over rational arithmetic. It uses dt=1/200 s, direct isotropic attitude information H_a=30 with R_a=1/25 at every IMU sample, magnetic H_m=75 with R_m=16/25 every eight samples, P_theta,0=1/2000 and P_bg,0=11/10000000, and the literal gyro/process densities used by the replay. Over one second the exact rational Riccati map satisfies P_0-P_1 positive definite, so this oracle cell is self-containing by Riccati monotonicity. The accumulated magnetic heading/BG information satisfies I_mag-I_2 positive definite exactly; det(I_mag-I_2) is about 0.98970433. Thus the magnetic geometry has strict margin even after intentionally strong direct attitude conditioning.

This is not yet a shipping service proof. Promotion requires a formal Schur/conditional-information comparison showing that the literal accelerometer+S nuisance chronology cannot suppress the magnetic incremental information more than this oracle, an all-time bound corresponding to H_a<=30 on the literal predicted specific force/attitude row, and treatment of the shipping reset coordinates. These dependencies are explicit in planar_oracle_rational.py and fail closed.

## Literal prediction and scheduler reuse

No independent reimplementation of the 12x12 LIN prediction is needed. The proof instrumentation already receives F_LL and Q_LL after the shipping calls to IntegratedOUChain::transition/process_covariance, together with F_AA,Q_AA, the BA phi/Q block and every correction/reset in operation order. The factor stream therefore treats those callback matrices as the literal coefficients. Likewise S=0 is observed only when the shipping time_update calls periodic_update_due and then applyIntegralZeroPseudoMeas; the factor stream consumes that actual correction event instead of approximating a cadence. This removes the surrogate F_LL/Q_LL and scheduler chronology from the intended certificate. A temporary non-promoting diagnostic exports the first-200 literal prediction retention and S-event count to verify plumbing before interval promotion.

## Remaining Poincare certificate

The same-history enclosure now needs two independent parity cells rather than a
generic 21x21 box. It must also carry the filtered adaptive parameters and
pseudo-measurement phase because those are causal scheduling variables. For one
20-s source period define the exact literal maps F_E,F_O on those cells.

A valid certificate must provide outward-rounded cells C_E,C_O,C_T such that

    F_E(C_E,C_T) subset int(C_E),
    F_O(C_O,C_T) subset int(C_O),
    F_T(C_T)     subset int(C_T),

and then propagate each 2x2 service probe through every sample root in a
one-second placed window. The phase-uniform result is

    min_phi min(lambda_min(I_E(phi)),lambda_min(I_O(phi))) > 1.

The 1200-s replay is only a seed-selection diagnostic. Tail sample-root sweeping
and 20-s Poincare drift are recorded separately and are not theorem evidence
until the interval inclusion above succeeds.

Structures preserved: literal 21-state chronology through an exact invariant
permutation; actual covariance/gain/Joseph/reset; actual accepted magnetic
events; coupled tuner/scheduler.

Relaxations introduced: none in the parity factorization. Any interval hull used
later must be recorded as an explicit outer enclosure.
