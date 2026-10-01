# Displacement-selected moment versus the literal four-S adjoint

This is the requested contraction before taking norms.  It finds a structural
boundary term; it does not promote a failed interior cancellation.

Choose t0<t1 in a complete displacement-excited MOVING window and
q=(p(t1)-p(t0))/|p(t1)-p(t0)|.  Put T=t1-t0.  Exactly,

    q^T[p(t1)-p(t0)]
      = T q^T v(t0) + integral_t0^t1 (t1-t) q^T a(t) dt             (D1)
      = T q^T v(t1) - integral_t0^t1 (t-t0) q^T a(t) dt.            (D2)

The left side is at least P_E.  These are signed identities on the SAME
physical history.

The literal four-S multiplier uses four actually applied S events s0<...<s3,

    c_j=-6/prod_(l!=j)(s_j-s_l),
    psi(t)=1/2 sum_j c_j (t-s_j)_+^2.

Its exterior jets satisfy

    psi=psi'=psi''=0

on both sides of the spline support.  Consequently its exact integration by
parts identity,

    integral psi a_hat
      = sum_j c_j r_S,j
        -sum_i(psi_i K_v,i-psi'_i K_p,i+psi''_i K_S,i) r_i
        -d_psi,                                                     (D3)

contains no physical/nominal endpoint v,p,S term.  That was the purpose of
the four-S construction.

Therefore D1/D2 are NOT in the scalar row space of D3 merely because both
contain acceleration.  The displacement chord requires the nonzero boundary
functional T q^T v(t0) (or T q^T v(t1)); the standard four-S spline annihilates
it.  Contracting the displacement-selected acceleration kernel against D3 and
then dropping the velocity boundary would be invalid.

This is useful: it identifies the minimal missing object.  The next exact
multiplier must augment the four-S/accelerometer adjoint by ONE boundary
velocity row in the selected chord direction, or use a nonzero-boundary
piecewise-cubic multiplier whose endpoint jet reproduces D1 while its interior
jumps remain the literal S/accelerometer weights.  Only after that completion
may the SLOW and FAST Abel bounds be applied.

The required augmented identity has the schematic exact form

    q^T Delta p
      = B_v(q;v0,v1)
        + sum_(actual S) C_S r_S
        + sum_(actual acc) C_a r_a
        + F_slow + F_fast + D_literal,                              (D4)

where B_v is RETAINED until it is either represented by the root/terminal
reader or cancelled by a second displacement chord.  F_slow uses predecessor
increments; F_fast uses the carried primitive U and all placed-window caps.
No independent residual box is allowed.

Hence displacement excitation has removed the zero-translation counterexample,
but the existing four-S interior adjoint by itself does not yet prove the
outer compatibility set empty.  The obstruction is now algebraically narrower:
construct and bound the boundary-velocity completion of D4, then intersect it
with gyro transport and actual MAGNETIC SERVICE.  No additional physical
assumption is indicated by this calculation.
