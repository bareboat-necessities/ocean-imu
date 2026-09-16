# Startup seed and Mahony arithmetic

`finite_mahony_prefix_totality.py` proves conditional scalar arithmetic
bounds for both commissioned IMU profiles. An ordinary `FromTwoVectors` seed
and the actual gains `two_kp=0.2f`, `two_ki=0.02f` satisfy the following bounds
through the 30,602-step timeout-plus-word arithmetic horizon.

| Quantity | Strict upper bound |
| --- | ---: |
| Ordinary seed squared quaternion norm | 1.06 |
| Squared quaternion norm after each update | 1.112 |
| Each Mahony feedback component | 1.3 |
| Each integral-feedback component | 4.1 |
| Squared-norm sum before quaternion normalization | 6 |
| Absolute vertical output, m/s² | 32 |

The ordinary-seed proof retains the cancellation in
`(v0.x²+v0.y²)/(2(1+v0.z))`. Its rounded-vector norm defect and branch lower
bound on `1+v0.z` are used together. The prefix proof bounds the actual
integral recurrence, then the Euler update, then applies the exhaustive
fast-inverse-square-root normalization bound. Nonnegative zero/subnormal
norm words are included. The arithmetic horizon does not establish startup
timeout reachability.

## Near-antiparallel seed

`finite_seed_eigen_svd.py` materializes the pinned Eigen 3.4.0 scalar QR/Jacobi
producer. Both QR pivots, rank-one and rank-two inputs, tiny tails, the actual
Jacobi loop and its IEEE special-value arithmetic are retained. The axis is
computed from the normalized source vectors. A supplied legacy witness must
equal that computed axis.

The reviewed header hashes in `finite_seed_svd_axis_reduction.py` bind the
source deduction. All post-QR rotations and sorting swaps touch only V columns
0 and 1, so the requested third column is exactly the two-reflector QR column.
`finite_seed_svd_roundoff.py` proves its squared norm below 100 and proves
uniform return in at most two Jacobi sweeps. The complete operation-level
argument and local-FMA arithmetic-family scope are in
[the scalar SVD proof](ou3-alt-svd-scalar-termination.md).

`first_svd_seed_bridge()` handles the first Mahony event explicitly; it does
not assume that the coarsely bounded seed is already in the normalized shell.
With initial integral zero, the actual gains and 5ms dt, and raw stored
per-component gyro/accel caps of 35 rad/s and 160 m/s², it proves:

| Quantity | Strict upper bound |
| --- | ---: |
| SVD seed squared quaternion norm | 128 |
| First feedback component | 144 |
| First integral component | .015 |
| First corrected rate component | 64 |
| First squared-norm sum before normalization | 2048 |
| Successor squared quaternion norm | 1.112 |
| First absolute vertical output, m/s² | 322 |

The common normalization theorem covers the scalar and permitted fused
inverse-square-root correction. The first successor therefore meets the
quaternion and integral premises of `finite_live_input_contract`'s all-time
invariant. A too-small initial accelerometer norm preserves the default
observer state until a qualifying sample arrives.

These source arithmetic theorems are distinct from actual target compiler
correspondence and startup admission. Native comparisons are regression
checks. Unused Euler-angle library outputs are outside the state/output
dependency theorem. Target qualification and the complete-word proof must
consume the relevant source and arithmetic-family premises explicitly.
