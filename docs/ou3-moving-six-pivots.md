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
Full covariance/noise-action and nonlinear stationary robustness do not follow
from this row bound alone.

The hidden-motion witness generates the same nominal history. Hence this
six-column nominal floor is fully compatible with physical attitude/BA
nonidentifiability. Removing physical stillness from MOVING does not, by
itself, improve that nominal floor or repair the failed signed-action bounds.
This subcase enters the common target only through the historical-reader
premise; it supplies no physical ambiguity removal or numerical rho.

## Scope of the executable certificate

`moving_pivots.py` verifies the scalar bounds with exact fractions and an
independent supplied rational matrix audit. Its nonsimultaneous/reset defect
budget is conditional; there is no replay-derived universal value. The
80-digit feasibility check evaluates this supplied algebra before the exact
check. It is not a contraction construction or a shipping feasibility claim.
The two-margin implication stays OPEN. The same moving signed-action norm
relaxations that failed before the physical split still fail. No source-uniform
B_*, J_AG, rho_0, retained radius or every-prefix certificate is promoted.
