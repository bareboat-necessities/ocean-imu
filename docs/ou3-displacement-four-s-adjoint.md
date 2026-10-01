# Boundary-completed displacement / four-S adjoint

The first contraction found the missing endpoint velocity row.  The existing
historical AG reader cannot supply it: that reader targets the six
attitude/gyro-bias root coordinates, whereas v(t0) is a LIN root coordinate.
Do not enlarge the AG theorem by silently treating v as nuisance with zero
action.

The exact completion uses the nonzero-boundary polynomial

    lambda(t)=T-t,  0<=t<=T.

It is a degree-one member of the allowed piecewise-cubic class.  Its jets are

    lambda(0)=T, lambda(T)=0,
    lambda'(0)=lambda'(T)=-1,
    lambda''=0.

Twice integrating the physical chain p'=v, v'=a gives exactly

    integral_0^T lambda(t) q^T a(t) dt
       = q^T[p(T)-p(0)] - T q^T v(0),                       (D1)

or

    q^T Delta p
       = T q^T v(0) + integral_0^T lambda q^T a dt.         (D2)

No norm has been taken.  There are no interior distributional atoms in lambda.

Now let psi be the existing literal four-S spline.  Its exterior value, first
derivative and second derivative all vanish.  Therefore for ANY scalar alpha,

    Lambda(t)=lambda(t)+alpha psi(t)

has exactly the same displacement boundary coefficients as lambda.  Its
interior S/accelerometer weights are precisely alpha times the existing
four-S weights, plus the smooth lambda accelerometer kernel.  Thus the
boundary completion does not change the scheduler, invent S observations or
alter the shipping estimator.

Schematically the literal forced-data identity can now be written

    q^T Delta p
      = T q^T v_0
        + sum_acc C_a r_a + sum_S C_S r_S
        + F_slow + F_fast + D_literal,                      (D3)

where C_a contains the sampled lambda*q physical-acceleration row together
with alpha times the existing four-S row, and C_S is alpha times the actual
four-S atom row.  The physical acceleration term is substituted through the
actual measurement equation BEFORE taking norms.  F_slow is then reduced by
the predecessor-linked slow Abel identity; F_fast by the carried U primitive
and all-placed-window cap.  Rotation/reference and arithmetic defects remain
in D_literal.  No independent residual box is reintroduced.

D3 is an exact algebraic completion, but it does NOT yet bound the boundary
velocity action.  The crude |T q^T v_0|<=T Vmax would usually consume P_E and
is therefore not the intended next inequality.  The boundary term must be
kept signed and linked to either (a) a second displacement chord/window, or
(b) the literal LIN root/terminal reader/covariance action.  This is now a
one-coordinate LIN boundary problem rather than an uncontrolled acceleration
functional.

Consequently the requested boundary-velocity COMPLETION IDENTITY is closed.
The source-uniform outer compatibility set is not yet proved empty.  The next
falsifiable calculation is to carry the q^T v root row through the existing
12-state LIN historical reader on the same interleaved S+accelerometer word
and evaluate its action jointly with the completed displacement functional.
