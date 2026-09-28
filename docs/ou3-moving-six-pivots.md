# Chronological gyro transport and six historical pivots

These are conditional lemmas on the existing historical AG reader, subordinate
to `B_* -> J_AG -> full covariance upper -> rho_0 < 1` in
`V_next <= rho V + c_d |d|^2`. They do not certify the source-uniform premises,
bound nonlinear supplies, or replace the signed temporal balance. The declared
complete-window excitation applies only inside physical moving episodes.

## 1. The complete chronological gyro block

In the fixed default deployment frame, the historical auxiliary AG transition
has form `T=[[A,B],[0,I]]`. Predictions have literal blocks `(R,D)` and the
attitude covariance reset is `G=I+[delta_theta]/2`. Applied observations do
not change the *latent auxiliary state*: their actual gains enter the optimal
covariance/readout comparison. They are never omitted from that comparison.
Mean projections retain their separate finite-error defects.

Maintain exact nonnegative majorants
`|A|<=a`, `|B-t I|<=b`, starting at `(a,b,t)=(1,0,0)` at an anchor.
For a prediction of duration h, suppose
`|R|<=r`, `|R-I|<=u`, `|D-h I|<=v`. Then

`a_next=r a`, `b_next=r b+u t+v`, `t_next=t+h`.

For a reset, `|G|<=r`, `|G-I|<=u` gives

`a_next=r a`, `b_next=r b+u t`, `t_next=t`.

These follow by expanding `RB+D-(t+h)I` or `GB-tI`. Consequently
`sigma_min(B)>=max(0,t-b)`. The bounds permit nonorthogonal resets and retain
their actual order. Frame changes require their literal maps and are outside
this particular fixed-frame formula. No independent box is substituted for
the nominal history: the *uniform* bounds on its reset/injection action are
still missing. For normalized gyro coordinates, replace hI and D by the
same scaled injection throughout; do not mix radians/s and normalized bias.

### Inverse-frame refinement retaining the literal resets

There is a stronger complementary bound that does not repeatedly amplify a
triangle defect by reset norms. It enters the same six-pivot premise below,
then B_*, J_AG and rho; it is not a new comparison filter or a nonlinear bound.
Put C=A^-1 B. For every literal reset G=I+[d]/2,

`G'G=I+(|d|^2 I-dd')/4 >= I`, `|G^-1|<=1`,
`|G^-1-I|=(|d|/2)/sqrt(1+|d|^2/4)<=min(1,|d|/2)`.

The Rodrigues prediction R is orthogonal. In the small-rate source branch,
its transverse singular values are sqrt(1+theta^4/4), theta=h|omega_hat|,
and its axial value is one. Thus both branches have `|R^-1|<=1` and every
chronological A has sigma_min(A)>=1, independently of the injection sizes.
In these fixed raw gyro coordinates the exact updates are

`prediction: C_next=C+A^-1 R^-1 D`, `reset: C_next=C`.

The reset still changes A^-1 and therefore every subsequent injection.
This cancellation does not remove any reset from the actual comparison.
Starting at `(beta,e,t,a)=(0,0,0,1)`, maintain
`|A^-1-I|<=beta`, `|C-tI|<=e`, `|A|<=a`. With theta an actual upper
bound on the prediction argument, a source-valid recurrence is

`e_next=e+h beta+h(theta/2+theta^2/3)`, `t_next=t+h`,
`beta_next=min(2,beta+theta+theta^2/2)`, `a_next=a(1+theta^4/8)`.

For a reset with |d|<=s, use
`beta_next=min(2,beta+min(1,s/2))`, `a_next=a(1+s^2/8)`, leaving e,t fixed.
For Rodrigues, `|R^-1 D-hI|<=h theta/2` by its exact integral. For the
small-rate branch, expand `D-hR=h^2 W/2-h^3 W^2/3`, then use inverse
nonexpansion. The common bound on R^-1-I follows from R-I and the same
inverse inequality. Product norms give the displayed beta bounds; since
|A^-1|<=1, its distance to I never exceeds two. Consequently, at every prefix,

`sigma_min(B)>=sigma_min(C)>=max(0,t-e)`.

The source .5-rad/s projection supplies the qualified theta bound, not a
bound on the preceding reset injections. One may take the maximum of this
floor and the direct recurrence's floor. Neither dominates in all cases.
A terminal reset cannot consume the new gyro floor: h=.005, omega_hat=0,
followed by d=4e_x still gives floor .005; the direct triangle bound gives
zero. This is an exact supplied algebra audit, not a claimed source history.

Inverse nonexpansion alone is insufficient. Start with one qualified
prediction h0=3/625, zero corrected rate, then the literal resets
d=(4,4,8/3)e_z and eight predictions h=1/200. On the transverse plane the
three reset multipliers have product
`(1+2i)^2(1+4i/3)=-25/3`. Therefore the inter-anchor gyro block is exactly
`B=diag(0,0,28/625)`, although every prediction D is a positive scalar identity
and every reset is nonsingular. This rejects an unrestricted-reset implication
from the one-step floor. It is **not** a shipping-reachable or magnetic-service
counterexample; intermediate observation rows could restore full historical
rank. The necessary source-uniform injection/row budget remains OPEN.

The implemented .5 rad/s gyro invariant bounds each prediction's R,D errors
using the existing qualified rate and source-polynomial allowance. It does
not bound the intervening injection sum. In particular, positive singular
values of individual D do not imply a positive singular value of their
signed chronological sum. This recurrence identifies the missing budget
without claiming the single-step certificate settles it.

## 2. Two groups of actual rows give a six-column bound

All matrices below use one fixed AG coordinate scaling. Suppose two groups
of historical applied sensor rows, transported to their own anchors, give
three-column matrices C0,C1 with
`sigma_min(C0)>=c`, `sigma_min(C1)>=c`, c>0. They may each have more than
three rows; individual magnetic/accelerometer skew blocks remain rank two.
Let A,B be their complete chronological inter-anchor map, with
`|A|<=a` and `sigma_min(B)>=b0>0` from the preceding lemma or a sharper bound.
Write their actual selected historical rows as

`O = diag(C0,C1) [[I,0],[A,B]] + E`, `|E|<=epsilon`.

E is an explicit *actual-row* defect: asynchronous sensor timing, within-group
gyro transport, resets, reference changes, and force variation must be charged
there or retained directly in C0,C1,A,B. This is not an assertion that events
are simultaneous, that the nominal force equals physical force, or that
innovations are sensor noise. Remaining historical rows only add information.

The inverse of the square block matrix is
`[[I,0],[-B^-1 A,B^-1]]`. Its norm is at most
`K=1+(a+1)/b0` (a deliberately conservative block triangle bound). Thus

`sigma_min(O)>=s:=c/K-epsilon`.

When s>0 this proves **all six** columns are independent, including the
gyro directions. This is a corrected sufficient premise, stronger than just
the proposed Delta_col/Delta_gyr margins. There is no proof that physical
span, the signed identities and two-column magnetic service imply uniform
c,b0,epsilon with s>0 on all realized moving histories. This budget is the
remaining falsifiable target, rather than the unsupported choice
`p=min(Delta_col,Delta_gyr)`.

### Exact grouping within one prediction cell

For regular zero-lever-arm operation, an **applied** accelerometer correction
followed by an **applied** magnetic correction with no prediction between them
gives an exact group. The literal AG Jacobians are `[C_acc,0]` and `[C_mag,0]`.
Let G_local be the chronological product of every intervening attitude reset,
including the accelerometer correction's reset. The group's matrix at its
accelerometer anchor is

`C=[C_acc; C_mag G_local]`.

Sync and observations themselves are identity maps on the latent AG state;
the corrections still enter the actual optimal covariance comparison. Thus
the group's raw historical rows are exactly `C [A_anchor,B_anchor]`.
For two such groups, in the first anchor's coordinates, the displayed
two-group factorization has **E=0**, even with nonzero injections. This is
an operation-order identity, not a simultaneity approximation. The wrapper
performs time prediction (and due S correction), then acc; a subsequently
forwarded accepted mag before the next prediction has this ordering. No
rejected packet or mere call counts as a row. A prediction separates cells;
frame/relock/reference hard events require their own retained maps and are
outside this regular grouping lemma.

`same_prediction_cell_groups` checks the exact identity on exported rational
operands, rejects direct BG columns and unsupported hard events, and retains
all noncommuting resets. This removes the generic epsilon charge for this
selection rule. Uniform positive c and b0, bounds on A and complete readout
action, and coverage by actual informative service remain OPEN. In the
common tail inequality this lemma supplies E=0 to the same six-column budget;
it does not supply rho or every-prefix nonlinear retention.

The remaining geometry has a precise reset-aware form. Write the literal
group as `C=[skew(f);skew(b)G]`, with signs immaterial, G the local reset
product, and nonzero actual nominal vectors f,b. Put `u=f/|f|` and
`k=G^-1 b/|G^-1 b|`. Since sigma_min(G)>=1, the magnetic block's two
nonzero singular values are at least |b| and its kernel is span(k). Hence

`C'C >= |f|^2(I-uu')+|b|^2(I-kk')`,
`sigma_min(C)^2 >= min(|f|^2,|b|^2)(1-|u'k|)`.

The last step uses eigenvalues 2,1+|u'k|,1-|u'k| of the sum of the two
projectors. Supplied lower vector norms and an upper absolute cosine therefore
supply c in the same six-column budget. These are unresolved actual-history
premises, not additional physical assumptions or independently chosen boxes.
The relevant field is **pulled back through all local resets**. Indeed the
relaxed group f=(1,0,1), b=e_x, d=2e_y has nonparallel raw f,b but Gf=2e_x,
so Cf=0. This is not a reachable source or whole-history rank counterexample.
The independent 80-digit noncommuting example gives conditional floor
8.54167423631316 below actual singular value 8.80825850549787; exact rational
tests audit the matrix inequality and the relaxed null vector separately.

### World-frame form and the limit of same-cell groups

`ou3-world-frame-rows.md` proves that each group satisfies
`C R=diag(-R,-R_2)[[f]x;[B_w]x N]`, with f the world nominal force, B_w the
committed reference and `N=R_2' G_local R`. The geometry is therefore
attitude-free, and the projector bound above sharpens to the exact least
eigenvalue `(F+B)/2-sqrt((F-B)^2/4+FB kappa^2)`. Literal injections obey
`dd'<=NIS P_theta,theta`, which bounds a and the inverse-frame b0 on short
words. A 1-Hz collinear MARINE MOTION/IMU BIAS history degenerates every
same-cell group while the aggregate array stays full rank; its cadence fails
MAGNETIC SERVICE. A jerk lemma excludes collinearity at every instant of a
cadence whose length-weighted mean gap is below about .051 s (h=1/5, 16-s
windows). A uniform c for this selection rule therefore needs a
magnetic-cadence coupling; the aggregate world-frame rows do not.

## 3. Quantitative bridge to the existing factor selector

If an m-by-6 array has `sigma_min(O)>=s>0`, then at each of its six
largest-residual row-selection stages the selected residual norm is at least
`s/sqrt(m)`. Indeed let Pi project onto the complement of the previously
selected row span, of dimension d>=1. Then

`sum_i |Pi row_i(O)|^2 = trace(Pi O' O Pi) >= d s^2`.

At least one row has squared residual >=d s^2/m>=s^2/m. Already selected
rows have zero residual, so the next row is distinct. This proof includes
every stage, not just the determinant of one convenient minor. Adding other
rows preserves the singular floor; use their actual total row count in m.

Taking a rational p>0 with `m p^2<=s^2` therefore supplies the six-pivot
premise of `separated_reader_action_implication` in `signed_temporal.py`.
Uniform bounds on operation count, actual observation/transport norms,
complete noise factors, terminal map and inherited nuisance covariance then
give finite B_*. Keep the full matrix process comparison for rho_0; the coarse
scalar finiteness estimate is not a usable nonlinear decay margin.

## 4. Exact quiet nominal subcase

There is a quantitative stationary subcase of the same reader. On a regular
zero-residual nominal execution with corrected rate zero, AW=0, identity
attitude, reference B e_x and g>=9, B>=9, select two **actually applied**
acc/mag pairs at the same sample phases, separated by eight qualified
4--6 ms predictions. This covers the 200-Hz/every-eighth-mag quiet record
after reference establishment without equating .005f to exact .005. The actual
zero corrections have identity resets; projections leave the zero means fixed.
No physical bias is forced to follow the estimator prior.

At each anchor C consists of the gravity and magnetic skew blocks, so
`C'C=diag(g^2,g^2+B^2,B^2)>=81 I3`. Across the anchors `A=I`,
`B=(sum h_i)I`, with `sum h_i>=8*.004=4/125`; the groups are at common
latent AG times and `E=0`. The preceding lemma gives `s=18/127` and all
six pivots at least `9/254` for these twelve raw rows.
These are analytic bounds on this precisely stated nominal subcase, not a
sampled Gramian or a uniform result for all quiet-driven trajectories. Applied
gains and all nuisance/measurement/process correlations still enter B_W.
Full covariance/noise-action does not follow from this row bound alone.
The explicit same-reader construction in `ou3-stationary-detectability.md`
now supplies a uniform historical action, every-operation AG/full covariance
upper comparison and qualitative homogeneous linear loss on this precise
quiet nominal subcase. It retains inherited nuisance covariance and all actual
conditioning. Nonlinear stationary robustness remains OPEN.

The hidden-motion witness generates the same nominal history. Hence this
six-column nominal floor is fully compatible with physical attitude/BA
nonidentifiability. Removing physical stillness from MOVING does not, by
itself, improve that nominal floor or repair the failed signed-action bounds.
This subcase enters the common target only through the historical-reader
premise; it supplies no physical ambiguity removal or numerical rho.
By the world-frame form its six-column floor does not need identity attitude
or a horizontal field: any constant attitude and any committed reference with
|B|>=20 uT and horizontal fraction >=1/5 give `c^2>=1296/481`, s=16/635 and
six pivots >=4/635. Its covariance action is not regenerated for that class.

## Scope of the executable certificate

`moving_pivots.py` verifies the scalar bounds with exact fractions and an
independent supplied rational matrix audit. Its nonsimultaneous/reset defect
budget is conditional; there is no replay-derived universal value. The
80-digit feasibility check evaluates this supplied algebra before the exact
check. It is not a contraction construction or a shipping feasibility claim.
The two-margin implication stays OPEN. The same moving signed-action norm
relaxations that failed before the physical split still fail. No source-uniform
B_*, J_AG, rho_0, retained radius or every-prefix certificate is promoted.
