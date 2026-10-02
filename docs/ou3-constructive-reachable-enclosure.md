# Constructive dependency-preserving shaped-supply enclosure

This is the constructive continuation of `ou3-source-uniform-reachable-graph.md`.
It obeys the permanent shipping-faithfulness protocol.

## C1. Shared-history cells only

A cover cell parameterizes only independent theorem inputs on one persistent
history: physical initial/boundary coordinates and admissible MARINE trajectory
degrees of freedom; SLOW bias root/rate degrees of freedom; FAST signed
primitive increments; delivered measurement/model/arithmetic uncertainty where
the theorem explicitly admits it; and magnetic-service delivered rows/events.

The following are NEVER independent cell axes: Mahony state, frequency,
variance, tau, sigma_aw, R_S, T_S, covariance, gains, scheduler phase, magnetic
reference or BA gate. They are propagated outputs of the literal shipping map.

## C2. One propagated dependency token

For one cell H, propagate a single enclosure X_j(H) through the literal
chronology. At every operation construct simultaneously

    L_j(H) = Q_j - A_j' Q_{j+1} A_j,
    d_j(H),

with Q_j generated from the same covariance/tuner state. Do not first bound L
and d separately. Accumulate the exact quadratic form

    q_j(e,u) =
      e' L_j e - 2 e' A_j' Q_{j+1} d_j - d_j' Q_{j+1} d_j.

All occurrences of a history generator retain the same affine/noise symbol.
Prediction/correction cancellation and tuner coboundaries therefore remain
visible. Summing gives one joint form q_W whose positive part is dissipation
minus physical supply.

## C3. Kernel-aware quotient

The complete superword homogeneous zero set is already excluded on qualified
MOVING histories. Use that theorem before division. On a history cell prove a
positive lower bound d_H on the homogeneous restriction of q_W after quotienting
the word-dependent compatibility line by the actual later-word action. Bound
the linked physical remainder on the SAME form by s_H. The cell certificate is

    S_W / D_W <= s_H/d_H.

No global lambda_min unrelated to the forcing direction is used.

## C4. Temporal physical coordinates

SLOW coordinates are propagated from one root using their derivative bounds;
they are not reselected samplewise. FAST coordinates use primitive increments
with all-placed-window constraints and are not reset at cell/word boundaries.
Physical acceleration contributions are represented by v,p,S boundary
variables plus the jerk/sampling remainder so the integrated chain telescopes
before interval magnitude is taken. MARINE 30-s span and displacement
conditions are imposed on the same trajectory parameters. Magnetic service is
imposed on actually applied informative rows.

## C5. Adaptive shipping propagation

For each history cell, the enclosure order is literally

    delivered sensor cell
      -> conditioning + measurement-only Mahony/frontend
      -> period/frequency/variance/bandpass
      -> staged tuner candidate
      -> next-sample commit of tau,sigma,R_S,T_S
      -> covariance/process/S scheduler
      -> literal gains and mean factors.

Clamps and gates create strata in HISTORY INPUT space. A branch is split only
when the same input cell can reach multiple literal branches. Generated outputs
are never bisected independently.

## C6. Constructive target

For 60/100-s superwords obtain finite leaves H_i satisfying

    D_W >= d_i W_0, d_i>0,
    S_W <= s_i,

with linked dependence retained. Then

    rho = 1-min_i d_i,
    E_phys = max_i s_i,

is allowed only if d_i and s_i came from the same cell construction. Prefer the
sharper cellwise absorbing bound max_i s_i/d_i rather than
(max s_i)/(min d_i).

The target is

    max_i s_i/d_i < W_entry,

together with every-prefix retention and a proved W_entry -> V<.15^2
comparison. If a cell is too wide, subdivide its shared history coordinate with
largest contribution to the JOINT quotient, not a generated tuner/gain output.

## C7. Current implementation boundary

`reachable_history_enclosure.py` now enforces the admissible subdivision and
shared-dependency contract and fails closed if dissipation has not been proved
positive after zero-set exclusion. Existing verified interval matrix arithmetic
may be reused underneath it, but the old independent operation-box admission
is not the theorem interface.

The remaining substantive implementation is the literal causal interval
propagator for frontend/tuner/covariance/gains and the joint quadratic-form
enclosure. This is a constructive numerical proof obligation, not an existence
lemma and not a replay experiment.

Structures preserved: all shipping-faithfulness protocol items.
Relaxations introduced: interval/affine enclosure over theorem input histories
is an outer numerical representation, but generated quantities remain causally
linked; failure of a coarse enclosure is class D, implementation error class E.
