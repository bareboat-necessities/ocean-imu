# OU-III sampling-map audit

The runtime change is relative to main `063017d44e3c422621a7f7b260f950cacb43a616`.
It changes live prediction forcing and finite-packet rejection. It does not
close the regional stability theorem or qualify hardware timing.

The frontend, parameter adaptation, magnetic-reference acquisition, S scheduler,
held/active bias equations, projections, measurement gains and reset chronology
retain their existing defaults. Prediction uses the average of the previous and
current raw gyro samples. The native mean, gyro-bias Jacobian and gyro process
noise use the same effective rate; lever-arm force terms retain endpoint rates.

Previous raw gyro, previous interval, mode and exposure configuration are
persistent execution memory. Bootstrap primes that memory before live prediction.
On a same-history comparison the cache is inherited and identical, rather than
reset or independently extremized. Averaging is a deterministic sensor-history
map, not an independent acceleration or attitude observation. Its norm does not
exceed the maximum of the endpoint norms. SF3--SF9's raw, averaged and weighted
forcing qualifications remain necessary; this audit does not discharge their
open physical SLOW+FAST envelopes.

For unspecified exposure the local noise-envelope multiplier is
`max(1, (1 + dt/dt_previous)/2)`. On the existing 4--6 ms qualified grid it lies
between 1 and 1.25. With the same effective rate, the added gyro covariance is
PSD; native lower-process comparisons remain conservative. Existing formulas
must not be interpreted as a newly proved recurring upper covariance bound.
Other exposure settings need their own bound. Correlated endpoint noise and
quadrature error remain separate from a nonlinear consistency theorem.

The first unprimed native sample does not invent a preceding rotation. Qualified
live tails inherit a primed cache. Nonfinite sensor packets are outside the
finite, bounded forcing premises; the wrapper rejects them before changing its
clock or histories. The focused rejection control compares startup and live
state/covariance with an identical valid-packet sequence. It does not establish
hardware gap timing or delayed-measurement compensation.

Native secants fork the complete inherited execution, including the cache, then
perturb the existing mean/covariance coordinates on the same delivered record.
Their unperturbed control and applied magnetic cadence must reproduce. These
are finite diagnostics with off-orbit perturbations, not derivative enclosures,
reachable-cell certificates, or all-time magnetic service.

Structures preserved: the existing OU/S operators, gain/reset algebra, physical
history and acquisition/release state. Relaxations introduced: none in runtime;
the secant experiment retains its explicit off-orbit diagnostic relaxation.
Stale source receipts are class E failures. A failed sufficient storage test
remains class D. No admitted shipping counterexample or completed theorem is
claimed. Fresh source-bound receipts require regeneration, not hash substitution.
