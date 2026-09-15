# ALT startup direction sampling

The commissioned continuous condition `r=m+d(xi)/dt`, with `|m|<=0.10`
and `|xi|<=1.5 s`, does not imply the same numerical condition for the
5 ms sample sequence. It therefore cannot be substituted into a discrete
Mahony comparison without a sampling argument. This is a proof-premise gap,
not a proof that the shipping observer fails to reach Live.

`finite_startup_direction_sampling.py` supplies an exact enclosure witness.
In the true world frame, let gravity be `(0,0,g)`, attitude identity, and

```
p_x(t) = (6.98/pi^2) sin(pi*t)
v_x(t) = (6.98/pi) cos(pi*t)
a_x(t) = -6.98 sin(pi*t).
```

The other translation coordinates, true gyro, and all biases are zero. This
is a 0.5 Hz physical wave with a bounded centered position primitive. It uses
no DC displacement, and zero bias/driver solves BIAS0, BIAS1 and BIAS2. Add
horizontal accelerometer residual `-0.12 cos(400*pi*t)`, a continuous 200 Hz
signal. This satisfies even the smaller declared fast-error allowance. No
noise frequency or derivative bound is imposed by the current source contract.

The true normalized measured down direction is consequently

```
u(t) = (6.98 sin(pi*t)+0.12 cos(400*pi*t), 0, g) / norm(...)
r(t) = u(t) - (0,0,1).
```

Take `m` as the continuous two-second mean of this periodic `r`. Its
zero-mean remainder has a continuous periodic primitive. Integrating from
the midpoint of the period bounds that primitive uniformly. At every sample,
`cos(400*pi*t_k)=1`, so the sample sequence sees a smooth sinusoid plus a constant
0.12 m/s^2 offset.

| Certified quantity | Bound, rounded outward for display |
| --- | ---: |
| Continuous mean direction norm | < 0.099858 |
| Continuous direction primitive norm | < 0.823858 s |
| Sampled mean x component | > 0.0090690 |
| Sampled mean z component | < -0.0997823 |
| Squared sampled mean norm | > 0.0100387 |
| Physical position norm | < 0.70723 m |
| Physical velocity norm | < 2.22181 m/s |
| Physical acceleration norm | 6.98 m/s^2 |
| Centered position primitive | < 0.45024 m s |

Machin's identity and alternating rational series bound pi. Rational Taylor
remainders bound each sine, and integer square-root enclosures bound the
normalization. A periodic midpoint-rule remainder bounds the continuous
baseline mean. For `u(x)=(x,0,g)/sqrt(x*x+g*g)`, differentiation gives
`|u''(x)|<=1/g^2`. Taylor's integral remainder and integration by parts bound
the continuous mean change from the 200 Hz residual by
`N*A/(400*g^2)+N^2/(2*g^2)<0.000097`. This is an
interval certificate, not a sampled search for a universal positive result.

If the sampled sequence had `r_k=m_k+(xi_(k+1)-xi_k)/dt` with the same
0.10 mean and uniformly bounded 1.5 s primitive, averaging successively more
whole periods would imply `|mean(r_k)|<=0.10`. The strict certified excess
contradicts that implication. The excess in this witness is small: it does
not by itself exclude a different decomposition on the 150-second startup
prefix, and it establishes neither startup noncapture nor an all-time
shipping trajectory that remains in startup.

The exact-real two-pole guard detector is bounded by 0.02445 m/s^2, including
its initial transient, below its 0.03 lower rail. Its coefficient bound and
difference recurrence are checked. The certificate does not promote that
real calculation to complete binary32 target/guard correspondence.

A startup proof must derive the needed sampled forcing bound from the actual
physical and residual histories, with its quadrature error included, or work
directly with the continuous source and discrete observer. The existing source
must not silently gain a sampled-observability requirement. Both universal
startup prerequisites and all rho/PASS gates remain open.
