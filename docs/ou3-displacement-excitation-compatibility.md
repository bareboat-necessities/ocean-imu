# Displacement excitation in the moving compatibility proof

The amended MARINE MOVING class requires both gravity-direction span and
wave-displacement span.  This removes every p=v=a=0 rocking construction from
MOVING, but it does not manufacture a pointwise acceleration lower bound.

Choose s<u in a complete displacement window with
|p(u)-p(s)|>=P_E and put T=u-s, q=(p(u)-p(s))/|p(u)-p(s)|.  The exact
kinematic identities are

    q.(p(u)-p(s)) = integral_s^u q.v(t) dt,

and

    q.(p(u)-p(s))
      = T q.v(s) + integral_s^u (u-t) q.a(t) dt
      = T q.v(u) - integral_s^u (t-s) q.a(t) dt.

Therefore some time has |q.v|>=P_E/T.  Also each signed acceleration moment
obeys the conditional lower bound

    |integral_s^u (u-t) q.a(t) dt| >= max(0,P_E-T V_max),

and similarly from the right endpoint.  The latter can be zero for ordinary
slow waves; it is not an acceleration-excitation assumption.

This is nevertheless useful in the existing proof because the literal
physical chain is p'=v, v'=a, S'=p with one persistent origin.  The four-S
spline and the two Abel identities already act on signed polynomial moments.
The correct moving compatibility problem is therefore the intersection of:

1. Delta_p(W)>=P_E and Delta_g(W)>=theta_E on their complete moving windows;
2. the exact p/v/a/S continuation and globally bounded displacement primitive;
3. SLOW predecessor-linked bias increments and FAST signed primitives;
4. the actual accelerometer and gyro equations;
5. actually applied MAGNETIC SERVICE and the same generated tuner/covariance.

A zero-dissipation or outer-entry candidate cannot now choose p independently
of the acceleration history used to drive the frontend.  In particular the
old p=v=a=0 sin^3 family is outside MOVING.  But no source-uniform contradiction
is yet numerical because T_P,P_E and both FAST profiles are OPEN.  Even after
qualification, the proof must compare the displacement-selected signed moments
with the literal S/accelerometer adjoint weights; P_E alone need not make
max(0,P_E-T_P V_max) positive.

The local field-axis theorem is unaffected and remains stronger inside
sqrt(V)<=.15: bounded velocity plus gravity/field noncollinearity already gives
the 17-s physical force/field separation.  Displacement excitation primarily
repairs the outer physical class and removes the zero-translation loophole.
The remaining controlling task is source-uniform outer-to-inner entry using
the same-history signed p/v/a/S plus SLOW+FAST accel/gyro and magnetic action.
