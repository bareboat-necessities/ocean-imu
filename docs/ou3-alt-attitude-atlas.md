# Finite attitude atlas and ungauged Live continuation

The finite runtime carries `(chart,z24,P21,q_hat,Reference)`. Chart selection
changes proof coordinates only. Physical attitude, nominal attitude, all motion
and bias coordinates, full shipping covariance, and the one-time Live/S origin
are retained. No source heading is excluded.

## Universal cover

For any nonzero relative quaternion q choose the smallest index k maximizing
`abs(q[k])`. The chart coordinates are `c[j]=2*q[j]/q[k]`, in ascending order
of the three indices other than k. Its homogeneous representative inserts 2
at position k and c in the other positions.

Since `sum(q[i]^2) <= 4*q[k]^2`, the selected normalized component has magnitude
at least 1/2. Therefore each coordinate is at most 2 in magnitude and
`||c||^2 <= 12`. Selection is invariant under every nonzero scale, including
the quaternion sign. Its deterministic tie rule covers chart boundaries.
The bound is an algebraic consequence, not a sampled enclosure or new source
restriction. The south quaternion `(0,0,0,1)` is chart 3, coordinates zero.

For `h_k(c)` the chart representative, its rotation is the homogeneous
quaternion rotation numerator divided by `D=4+||c||^2`. Hence `D>=4` on every
chart. The coefficient tests prove `N(q)'N(q)=||q||^4 I`,
`N(q*r)=N(q)N(r)` and `||q*r||^2=||q||^2||r||^2` over the exact polynomial
ring. In particular nonzero quaternion products never lose every chart.

## Exact event and chart transport

Every represented attitude update has `q_next=L*h_k(c)*R`. Prediction uses
the same physical rotation increment on the left and the conjugated nominal
increment on the right. Measurement injection uses identity on the left and
the conjugate of the same `(w,k_inj*d)` injection quaternion on the right.
No physical increment is replaced by a sampled estimator increment.

Let A be the 4-by-4 matrix of that quaternion product. Select the largest
component l of `u=A*h_k(c)`. For each i different from l,

`c_next[i] = (4*A[i,k] + 2*sum_j A[i,j]*c[j]) / u[l]`.

The sum omits k. This is an exact finite affine/rational relation with its
same-state denominator; all 16 source/target chart combinations are covered.
No Jacobian is substituted. Consecutive maps compose by ordinary substitution,
including chart switches at a Cayley pole.

For sensors, the actual shipping H remains unchanged. In chart zero the existing
Cayley secant supplies the finite residual identity. In charts 1..3 the attitude
term `(R(h_k(c))-I)*f_hat` is retained as a rational function of the SAME c.
The remaining terms include `R_true*e_aw`, held/active `e_ba`, and the original
physical residual. This attitude term is never an independent disturbance.
The inverse-free descriptor retains `S*solve = Hbar*z + offset(c) + nu`,
`d=N_theta*solve`, the affine chart transport and the true-bias projection.
The gain, Joseph covariance and local-error covariance reset reuse the original
21-state algebra. Shipping covariance is not a covariance of atlas coordinates.

Tilt resets and magnetic yaw writes derive a new chart from their SAME updated
nominal quaternion. H18/A21/hold edges and covariance-only events preserve the
chart index. The legacy local Cayley lemmas remain valid on their stated chart.

## Ungauged runtime branch

`finite_live_interleave.from_startup` accepts the normalized ungauged proxy seed
and keeps the unready magnetic accumulator, gravity gate, private observer,
continuous-calibration statistics and bias lock. Ungauged IMU events continue
the gravity gate from the same raw packet and the post-update MEKF attitude.
Source qualification and all required norm/LPF arithmetic remain explicit.

The dual-clock magnetic composer accumulates continuous statistics once, then
executes initial north acquisition. Shipping selects these distinct frames:

| Operation | Frame source |
| --- | --- |
| Initial acquisition before Live | Private startup observer |
| Initial acquisition after ungauged Live | Current MEKF boat quaternion |
| Gravity gate after Live | Current MEKF boat quaternion |
| Continuous hard-iron accumulation | Private startup observer |
| Refinement | Private startup observer |

While north is absent, no inner magnetic measurement or bias-unlock count is
invented. On the actual north-ready event, the composer writes the reference
and absolute yaw, executes any due refinement/application, and passes that same
packet to the MEKF. There is no second continuous accumulation. The post-gauged
service clock starts at this event; the physical Live/S origin does not move.

The executable regressions cover waiting -> IMU -> initial north -> next IMU,
simultaneous initial north/refinement, distinct MEKF/private tilts, and rejection
of premature measurement operands. These are conditional event identities,
not proof that every admitted physical history reaches north.

## Remaining qualification boundary

The atlas removes the single-chart representation obstruction. It does not
prove universal finite startup, magnetic accuracy, target arithmetic, counter
safety, or a complete source-uniform 600-transition deployment word.

Chart k>0 has its coordinate origin at a 180-degree rotation, not zero error.
Consequently a common `z'Mz` with every chart origin treated as zero would not
be a coercive attitude storage. Any future storage must respect the actual
rotation and compatible chart transport (for example using chart-dependent
homogeneous/affine terms). No storage or rho search is authorized by this lemma.

The magnetic control graph now rejects the undefined signed-int32 successor
at INT_MAX. This is faithful partial program semantics, not a proof that all
admitted histories avoid it. MAG-CALL-SCHEDULE-v1 still supplies no uniform
upper call count. The AtomS3R sketch makes at most one updateMag call per
updateFilter_ invocation, but using that bound requires attachment of that
specific caller and its actual variable-dt/configuration/build profile. It is
not silently substituted for the broader current asynchronous theorem language.
