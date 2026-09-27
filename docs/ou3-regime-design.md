# STILL, TRANSITION and MOVING on the carried OU-III execution

## Design decision and theorem scope

The single route remains construction -> capture -> magnetically informed
H18 -> reference refinement/release -> recurring A21 -> regional practical
stability. Physical regimes qualify portions of that route; they do not
replace H18 by an indefinitely held-bias estimator or restart A21.

The proposed automatic **certified** STILL switch is blocked by an exact
identifiability result below. No shipping state update, detector, zero-velocity
observation, bias hold, covariance reset or noise change is justified by quiet
packets alone. This is not a proof that the existing estimator needs no change:
stationary robustness of the full nonlinear shipping execution remains OPEN.
The executable changes implement the regime contracts, conditional bridge
algebra, a necessary-evidence monitor for falsification, and counterexample
regressions. The monitor is proof tooling, not an enabled shipping detector.

Every result here enters the one target
`V_next <= rho V + c_d |d|^2` as indicated. No finite replay supplies a uniform
certificate. STILL may ultimately use distance to an observable consistency
class plus a deterministic ambiguity radius, rather than zero full-state error.

## What is observable at physical rest

For true rest, in the fixed deployed frame and after accounting for temperature,
lever-arm and calibration terms, the observations are

`gyro=b_g+n_g`, `f=-g Q e_z+b_a+n_a`, `m=Q B+b_m+n_m`.

Thus zero physical rate gives direct gyro-bias information, but not exact bias
in bounded noise. A time average over a dwell D estimates the *terminal* bias
with error at most `N_g+D_g D/2`; a certified rate tube `|omega|<=eps_w` adds
`eps_w`. Deterministic bounded noise does not decrease as `1/sqrt(N)`. The
discrete version charges the weighted sample ages. Longer dwell alone cannot
make this bound tend to zero.

Zero physical v/a would support stationary observations **if** certified.
Neither velocity nor absolute wave displacement is observed by an IMU packet.
On an indefinite rest continuation, the bounded primitive excludes nonzero
constant wave displacement, so p=0. The already carried physical S is constant,
not necessarily zero: its capture origin is never reset.

Even with known nonvertical B and noiseless measurements, an attitude rotation
about B can be traded against b_a. Locally the stationary sensor map on
`(theta,b_a)` has a one-dimensional kernel:
`theta parallel Q B`, `delta b_a=g delta(Q e_z)`. Additional unknown magnetic
reference/residual terms only enlarge the ambiguity. The existing quiet-bias
example has two indistinguishable truths separated by about .02 rad and
.196 m/s^2. Ordinary moving six-column observability does not separate them.
The estimator BA OU prior is not physical evidence removing this kernel.

The required STILL theorem is boundedness/practical stability relative to this
consistency class, with an explicit ambiguity tube for physical attitude/BA,
bounded wave states, covariance comparisons, magnetic-reference uncertainty,
and every-prefix retention for arbitrary duration. A covariance decrease along
an estimator OU mode cannot be reported as shrinking the deterministic physical
ambiguity. Keep that uncertainty as a separate set/supply unless a justified
covariance rule is proved. The existing projections give `|e_bg|<=.52` and,
on completed BA projections, `|e_ba|<=B_a+.4`; they do not prove the remaining
state or covariance theorem. This observability result sets an irreducible
practical radius, not a value of rho.

## Exact rest/motion indistinguishability

Let alpha=1/1000, nu=1/40, and on a moving episode [a,b] of an integer number
of periods put `phi(t)=alpha sin^3(nu(t-a))`. Outside it put phi=0. Use
`R=Rx(phi)`, p=v=a=0, B=75 e_x,
`b_g=-phi' e_x`, `b_a=g(R' e_z-e_z)`, and zero fast residuals.
Both phi' and phi'' vanish at the joins. The physical bias histories are
locally absolutely continuous across rest/motion joins and satisfy

`|b_a|<=g alpha`, `|dot b_a|<=3g alpha nu`,
`|b_g|<=3 alpha nu`, `|dot b_g|<=9 alpha nu^2`.

All four are strictly inside the unchanged qualification. Translation, jerk
and primitive are zero. Every complete `T_E=2pi/nu=80pi` window **contained in
the moving episode** has gravity span `2alpha=.002`. No numerical deployment
qualification of these symbolic parameters is asserted. The measured record
is identically `gyro=0`, `f=-g e_z`, `m=75 e_x`, before, during and after motion.
The same construction works for a nonempty small-alpha/slow-nu parameter
family; the existing assumptions do not exclude that family.

For any deterministic causal detector (including all its dwell/hysteresis,
wave states, gains, covariance, clocks and magnetic reference), equality of
input prefixes and initial state implies equality of its entire internal
history by induction. Randomization cannot give a sure guarantee either.
Consequently finite entry on indefinite rest and guaranteed exit on every
admitted physical departure cannot both hold. Motion may start after any
finite rest entry time and continue indefinitely without a detectable packet
change. False entry and delayed exit have no general finite upper bound.

The real-arithmetic witness uses g=9.80665. The native regression separately
charges the constant residual `(g-g_float)e_z`, of norm
.0000001617431640625 m/s^2, to produce the literal float32 quiet packet.
It does not silently redefine physical gravity or a physical bias bound.

The actual nominal execution is the existing quiet construction, so the
stationary applied MAGNETIC SERVICE argument is retained. Its true-axis
projection costs at most `cos(alpha)^2` at each root, bounded below by
`(1-alpha^2/2)^2`. This leaves the already proved stationary service reserve
above one; this is not information from attempted callbacks. The construction
does not claim instability or refute the signed nominal margins: its nominal
force and field remain separated and the gyro transport is ordinary.

## Detector requirements derived from that obstruction

A sound implementation must distinguish necessary quiet evidence from physical
certification. Under the stated bounds, exact rest necessarily satisfies

`|gyro_i|<=B_g+N_g`, `||f_i|-g|<=B_a+N_a`,
`|gyro_j-gyro_i|<=min(2B_g,D_g |t_j-t_i|)+2N_g`,
`|f_j-f_i|<=min(2B_a,D_a |t_j-t_i|)+2N_a`.

The reference monitor checks these against the first sample of a continuous
dwell. Invalid packets, nonincreasing time or excessive gaps invalidate that
dwell. An entry dwell must be supplied explicitly; it is not fitted to traces
or claimed to improve deterministic observability. Entry is delayed and exit
on inconsistent evidence is immediate (temporal hysteresis). No amplitude
threshold tighter than the physical/noise envelope is called a necessary
rest condition. Consistent samples are only STILL-compatible. They never
certify STILL, zero velocity, excited MOVING, or informative magnetic service.
Estimated v/a cannot resolve the exact identical-input example.

Before enabling an estimator-changing switch, either prove a stationary
practical theorem for the **entire measurement-compatible class**, including
hidden motion, or provide independently justified stationary evidence. Neither
is added as a new physical assumption in this continuation. Numeric qualification
and an embedded implementation of such a certified switch remain OPEN.

## Complete-window excitation and finite bridges

A physical rest episode is a nondegenerate interval of zero physical v/a/rate.
Use maximal moving episodes I=(a,b) between those intervals (including the
initial or final unbounded episode). Isolated zero-rate instants do not split
a moving episode; periodic rocking retains one episode.
They are properties of one continuous history, not labels that can be restarted
at each proof word. The revised excitation quantifier is precisely

`for every moving episode I and every t with [t,t+T_E] subset closure(I):`
`Delta_g([t,t+T_E]) >= theta_E > 0`.

No excitation requirement is imposed on windows crossing a rest boundary.
A departure at a first permits a complete excited window at a+T_E. Shorter
moving episodes have no full excited window and are wholly transition for this
proof. A long constant-attitude translation cannot evade the requirement by
relabeling successive transition intervals. The all-time episode/window
certificate is explicit and cannot be inferred from a finite trace.

The theorem regimes are STILL -> TRANSITION -> MOVING and the reverse.
The certification layer retains a physical episode's boundary separately from
the detector's delayed decisions. MOVING entry additionally requires the
existing A21 release, applied magnetic information, historical-reader horizon,
retained-region entry and covariance premises. The 17-s regular-A21 clock is
not reset just because a physical label changes, and no prior history is
silently discarded. H18/refinement during physical rest still need their own
finite completion proof; STILL is not an external BA hold.

At every bridge operation let `sqrt(V_i)<=g_i sqrt(V_(i-1))+s_i`, including
the actual reset/projection/hard-event supplies and any storage comparison at
a regime change. Induction gives

`sqrt(V_k)<=prod_(i<=k)g_i sqrt(V_0)+sum_(i<=k)(prod_(i<j<=k)g_j)s_i`.

This is finite for each finite completed word. A retained bridge needs these
bounds at **every prefix**, inside the domain used to bound its operations.
No general finite detection duration or source-uniform prefix bounds have
been proved, so this algebra alone is not a certified bridge. For any eta>0,
its endpoint G,S gives `V_end<=(1+eta)G^2 V_0+(1+1/eta)S^2`.
This identifies its contribution to rho and the supply in the common target.

Repeated finite bridges need an additional *derived storage budget*, not an
invented motion assumption: accumulated bridge gains and storage-comparison
factors must be outweighed by contraction on certified words, with bounded
transported supplies and prefix gains. For example, a cycle with moving norm
factor q=.5 and bridge norm factor 3 has factor 1.5 and grows indefinitely.
Every bridge is finite; the repeated composition is nevertheless unstable.
Thus eventual single-regime retention, or a proved recurring product bound,
is an OPEN obligation. Arbitrarily many short motion episodes cannot simply
be declared harmless.

## Moving proof continuation

Removing indefinite stillness fixes the window quantifier and separates the
physical BA ambiguity. It supplies no new control of actual innovations,
nominal AW, reference action or reset sums on an excited episode. In particular
the existing moving forced-balance diagnostic already obeys excitation; its
failed rotation and signed-velocity norm budgets remain failed. Do not repeat
that relaxation. The rest/motion witness also invalidates a mandatory nominal
tilt response on the revised moving branch.

`ou3-moving-six-pivots.md` provides a new conditional bridge from chronological
gyro transport and two actual sensor-row groups to **all six** historical
pivots. It retains an explicit row-defect budget, reset action and coordinate
normalization, then invokes the existing historical action comparison. Its
uniform premises remain open; neither positive one-step gyro transport nor
physical span discharges them. Thus B_*, J_AG, full covariance upper bound,
rho_0, the nonlinear radius and every-prefix retention remain OPEN.

The same reader also has an explicit quiet nominal subcase: two actual
coincident acc/mag groups eight qualified predictions apart at the zero-residual
stationary nominal solution give six-column floor 18/127 and all-six-pivot
floor 9/254 in raw AG
scaling. This analytic result reinforces the distinction from physical
attitude/BA ambiguity; it does not prove stationary nonlinear robustness.
