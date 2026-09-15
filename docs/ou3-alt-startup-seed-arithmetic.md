# Startup seed and scalar Mahony arithmetic

The scalar startup word has a conditional totality certificate in
`finite_mahony_prefix_totality.py`. For both commissioned IMU profiles, an
ordinary `FromTwoVectors` seed and the actual gains `two_kp=0.2f`,
`two_ki=0.02f` give the following bounds through the existing 30,602-step
timeout-plus-word horizon. They assume the named scalar binary32 profile,
correctly rounded square root, and continuation of the same sensor bounds.

| Quantity | Strict upper bound |
| --- | ---: |
| Ordinary seed squared quaternion norm | 1.06 |
| Squared quaternion norm after each update | 1.112 |
| Each Mahony feedback component | 1.3 |
| Each integral-feedback component | 4.1 |
| Squared-norm sum before quaternion normalization | 6 |
| Absolute vertical output, m/s² | 32 |

The ordinary-seed proof retains the cancellation in
`(v0.x²+v0.y²)/(2(1+v0.z))`. It uses the rounded vector's norm defect jointly
with the branch lower bound on `1+v0.z`; separate bounds on the cross product
and reciprocal would lose this result near the branch boundary. The prefix
proof bounds the actual integral recurrence, then the quaternion Euler update,
then applies the exhaustive fast-inverse-square-root normalization bound to
close the induction. It includes nonnegative zero/subnormal norm words. It
does not use tilt accuracy or covariance as an error bound.

The finite horizon is a conditional arithmetic domain. It does not establish
that startup reaches Live by the timeout. The discarded Mahony Euler-angle
outputs and their library calls are outside this state/output dependency
theorem. The near-antiparallel seed and actual compiler/library correspondence
remain required for totality over every startup branch.

## Near-antiparallel solver reduction

The pinned Arduino Eigen 0.3.2 package contains Eigen 3.4.0. The reviewed
`Quaternion.h`, `JacobiSVD.h`, `ColPivHouseholderQR.h`, `Householder.h` and
`HouseholderSequence.h` hashes are in `finite_seed_svd_axis_reduction.py`.
`audit_headers` rejects a different header payload even if its version label
matches.

For the actual real `JacobiSVD<Matrix<float,2,3>>(m, ComputeFullV)`:

1. `diagSize=min(2,3)=2`. The preconditioner factors the scaled matrix's
   3×2 adjoint with column-pivoted Householder QR and writes its full 3×3 Q to V.
2. The Jacobi loop can rotate only columns `(1,0)` of V. The real sign
   adjustment changes U only. Sorting can swap only columns `(0,1)` of V.
3. Consequently V column 2 is unchanged after QR, bit for bit, on every
   returning execution. This applies to rank-one and rank-two inputs and does
   not need a singular-value gap or a unique nullspace axis.
4. The axis consumed by `FromTwoVectors` therefore equals column 2 of the
   two-reflector QR Q. No Jacobi rotation-roundoff accumulation enters that
   axis.

This reduces the axis proof to the scaled source input, pivot choice, two
Householder constructions, and Q evaluation. It does not remove the need to
prove termination of the shipping Jacobi loop: the function must return before
the column is read. Replacing the solver by a cross product, an arbitrary
unit nullspace axis, or an unqualified native result would not prove shipping
correspondence.

The remaining limiting obligations are source-uniform QR arithmetic and
axis bounds, source-uniform termination of the 2×2 Jacobi iteration, and actual
target arithmetic correspondence. The next falsifiable step is to construct
the two-reflector target operation graph and establish those bounds together
with a finite iteration bound. A rank-deficient witness and a near-threshold
input must be retained. Failure of that attempt would require reconsidering
the solver proof architecture; subdivision or host replay alone would not
resolve it.
