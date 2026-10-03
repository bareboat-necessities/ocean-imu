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
