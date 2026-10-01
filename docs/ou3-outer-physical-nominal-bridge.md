# Outer same-history physical-to-nominal accelerometer bridge

The correct outer object is not a pointwise |a_hat-a| bound. The earlier
17-s sync-slab derivation already produced the right signed functional:
prediction plus actual S forcing telescopes physical p/v/S boundaries, while
the accelerometer reader jump remains in the SAME Joseph port. Summing slabs
cancels all internal physical boundaries before norms. The result is one
endpoint+accelerometer functional.

That route was previously abandoned because FAST gyro was only amplitude
bounded. The admitted witness n_g=.02 sin(t) b generated about .055 m/s^2 of
signed field-axis supply, exceeding the old .0338 margin. Under the candidate
new temporal profile it is no longer admissible: its signed accumulation over
placed 60-s windows is O(.04 rad), whereas C_g=.0042 rad. Thus the historical
negative result does not transfer to the SLOW+FAST model.

For a complete candidate-profile proof word, define Phi_W as the already
derived augmented full-root endpoint+accelerometer functional, with physical
p/v/a/S boundaries telescoped before norms and all acc/S/mag correction jumps
retained through their actual gains. Its difference between physical and
nominal histories decomposes into:

  E_IMU = predecessor-linked slow accel/gyro Abel terms
        + fast accel/gyro carried-primitive terms,

  E_state = the TWO external full-root endpoint reader actions
          + nonlinear SO(3)/reset/arithmetic defects.

There is no interior pointwise AW error term. Every interior base innovation
is eliminated by the same-history adjoint/Joseph identity. This is the desired
outer bridge architecture.

The H18/A21 work on PR #640 now supplies what the old derivation lacked for
the endpoint terms: compact A21 release, held-H18 LIN BIBO, q^T v boundary
control, and compact covariance/coefficient families. Therefore E_state is
finite source-uniformly on any retained compact outer set. However FINITE is
not enough for the six-column contradiction: a numerical upper bound smaller
than the physical reserve is still required.

For the hypothetical profile the already-audited physical attitude reserve
after slow+fast IMU charges is about .00305612 rad (0.1751 deg). To close the
outer bridge one must evaluate the actual augmented endpoint-reader operator
on the compact release/outer class and prove its normalized state/nonlinear
charge is below this reserve. The repository currently has no source-uniform
numerical operator norm for that endpoint functional. Carried 17-s/100-s
values cannot be promoted.

Thus the new temporal model DOES remove the previous explicit fast-gyro
obstruction and reduces E_state to endpoint/nonlinear action rather than a
pointwise AW tube. But it does not yet provide the numerical E_state needed
for J_AG>0. The next falsifiable calculation is a source-uniform bound on the
existing augmented endpoint-reader operator over the compact A21 release/outer
class; do not reintroduce interior innovation or AW amplitude boxes.
