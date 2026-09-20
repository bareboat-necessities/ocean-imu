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

If the varying-history proof needs motion excitation, the added marine-motion premise
must state a positive lower information bound over every finite window, with
actual sample/event timing and nuisance cancellation included. Merely requiring
nonzero motion is inadequate: amplitudes can tend to zero, and no uniform
positive margin follows from strict nonzero amplitude. A physical excitation
condition must be derived into the required full-state inequality; assuming
that inequality outright would conceal the outstanding proof obligation.

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
