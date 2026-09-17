# Exact startup docking and candidate invariant tail

The scalar binary32 reset-to-docking calculation reaches the seven stored
Mahony components at sample 30,002 without installing estimator state.
`finite_startup_trap.py` executes 3,660,244 literal rounded operations from
construction and audits every packet against the same analytic physical
circle/yaw source used by `finite_startup_timeout_alignment_obstruction.py`.
The source retains both commissioned sensor profiles and zero BIAS0/1/2 driver.
The report is `reports/results/ou3_alt_storage/startup-trap.json`.

This is a reached finite state and a candidate infinite continuation, not a
proof of indefinite startup failure. All promotion flags remain false.

## Candidate all-time argument

Write h(q)=R(q/||q||)^T e_z, h0=h(q0), I0 for the attained integral,
r=1e-4, k=16 and gamma=10.2264552. The source shadow recursively generates

- a = RN32(gamma h);
- e = the literal rounded Mahony half-error on (q,a);
- mu = (I0 + (0,0,0.6)) dot h;
- g = RN32(mu h + k(h0 cross h) - I0 - Kp e).

The accelerometer is antialigned with the observer's measured gravity.
These are sensor packets; the shadow never writes native filter state.
This recursive tail differs from the older double-precision docking controller.

For exact normalized Euler evolution with angular rate
w=mu h+k(h0 cross h), c=h0 dot h and d=||h-h0||, direct quaternion algebra gives

`d_next²/d² = 1 - (1+c)(dt k - dt² k² c/2)/(1+dt² ||w||²/4)`.

On d<=r, ||w||<1 this is at most 0.846401, below 0.93². The executable
checks the rational inequality and exhaustive covering cells for the actual
fast inverse square root. The resulting quaternion squared-norm enclosure
is inside [0.99,1.01].

The proposed rounding budgets are ||e||<=16u and direction defect<=64u,
u=2^-24. Conditional on those budgets, the integral increment is below
9.54e-11, inside every attained integral rounding half-cell; hence I stays
bitwise fixed. The radius successor is bounded by
`0.93r+64u = 9.6814697265625e-5 < r`.
The candidate gyro residual bound is below 0.002 rad/s and accelerometer
residual below 0.008124 m/s², inside the existing source contracts.
The centre-to-physical acceleration bound is checked with rational Taylor
enclosures, independently of native floating trigonometry.

**The two rounding budgets still require a complete operation-by-operation
majorant audit.** Parallel exact acceleration/half-gravity cancellation,
rounded cross products, gyro cancellation, Euler products/additions and final
normalization must all be accounted for, including subnormal additive errors.
The constants are not certified merely because the final inequalities pass.
`rounding_majorants_audit_closed` and
`conditional_observer_tail_invariant_closed` therefore remain false.

After that audit, composition must prove indefinitely that the actual guard
stays transparent, the public quaternion/world-acceleration/low-pass gate
stays on the non-Live side, and wrapper clocks/operations remain total.
Target firmware correspondence is another independent obligation.
The native 600-sample extension has zero guard engagement, remains non-Live,
and has minimum tail gate z >2.53; it proves only this finite prefix.

## Reproduction

From the repository root, with Eigen available at ../eigen:

```sh
g++ -std=c++20 -O1 -ffp-contract=off -fno-fast-math -DEIGEN_DONT_VECTORIZE -I../eigen -Isrc tests/ou3_alt_contraction/startup_antialigned_docking.cpp -o /tmp/ou3-docking
/tmp/ou3-docking /tmp/ou3-prefix.txt 30602
PYTHONPATH=.:tools/stability python3 tools/stability/ou3_alt_contraction/finite_startup_trap.py --prefix /tmp/ou3-prefix.txt --tail-output /tmp/ou3-tail.txt --output /tmp/ou3-startup-trap.json
/tmp/ou3-docking --verify-packets /tmp/ou3-tail.txt
```

The exact Python replay takes several minutes. Native verification requires
all seven docking components to match and verifies finite wrapper guard/gate
behavior. Compare its final quaternion with the report as an independent
continuation check. Neither native execution nor this scalar model establishes
the whole target's compiler/FCR/library correspondence.
