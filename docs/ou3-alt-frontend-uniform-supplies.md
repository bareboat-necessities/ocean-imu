# Configured frontend finite supplies

`finite_frontend_uniform_bounds.py` proves an exact rational induction for the
frontend arithmetic consumed by a complete startup or Live IMU event. It uses
the existing physical MEMS/Mahony vertical bound, the shipping 6 Hz tracker
LPF, the default adaptive band, and the default statistics configuration.
Every finite history from the literal frontend reset is covered. No startup
deadline, source phase, or numeric ISS `W` cap is assumed.

| Quantity | Uniform bound |
|---|---:|
| Mahony vertical input | absolute value ≤512 |
| Tracker LPF and band lowpass state | absolute value ≤1024 |
| Band output, consumed by the same sample's statistics | absolute value ≤2048 |
| Unit-noise covariance | `p00≤2`, `abs(p01)≤2`, `p11≤8` |
| Stillness energy | ≤32768 |
| Statistics mean and square numerators | absolute mean ≤4096; square ≤8388608 |
| Both statistics weights | [0,2] |
| Band noise sigma | <0.5 |
| Configured sigma target | ≤4 |

All scalar no-reassociation contraction alternatives are included. The
exponential calls retain their actual rounded arguments and allow relative
error `2^-20`, plus the final rounding cell. The pinned target library proof
covers all frontend arguments in `[-60,0]` within that allowance. The runtime
band, statistics, LPF, stillness and sigma producers use this error relation.
The square-root witness relation remains correctly rounded, as supported by
the pinned target scalar/library profile.

The covariance induction retains `q=1-alpha` exactly. For the reached binary32
exponential values, Sterbenz subtraction makes the complement exact. Its
absolute monomial sums satisfy

```
p00' ≤ 2 - 4*alpha_low + 3*alpha_low²
abs(p01') ≤ 2*(1-alpha_low)
p11' ≤ 8 - 12*alpha_high + 7*alpha_high².
```

The coefficient intervals cover every configured rounded frequency/horizon
argument, including arbitrary switching and the distinction between lagged
band frequency and current statistics frequency. A conservative 32-operation
rounding factor and additive gradual-underflow allowance leave strict positive
margins in every state bound. This is an induction, not an endpoint simulation.

The variance readout is evaluated only when both actual weights exceed the
shipping `1e-6` readiness threshold. The eager squared debiased mean is the
largest bounded ordinary intermediate, below `1.678e19`; it is included even
when the subsequent variance is floored to zero. The bound therefore supports
execution before the sigma clamp, rather than inferring finite arithmetic from
the clamp alone.

Bounded startup and strongest Live wrappers check the configured constants and
the retained predecessor/successor states. Each event binds the same Mahony
vertical to the band, the same new band output to statistics, and the resulting
variance/noise/stillness to sigma. No history is reset or replaced.

This supplies the frontend part of complete-event totality. Existing finite
runtime constructor count limits remain separate limitations; the induction
does not justify extending those limits or prove a startup deadline. Startup
seed capture, full frontend compiler correspondence, MEKF covariance/Joseph
totality, and startup reaching Live are not proved here. Under the dormant
guard and pre-Live RAO branch, Racc performs no inflation and either preserves
the current value or restores nominal Racc; the nominal covariance and its
MEKF use retain their own source/totality obligation.
