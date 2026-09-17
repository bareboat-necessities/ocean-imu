# History-dependent capture and nonuniform ALT contraction

## Answer and quantifiers

A numerical maximum capture time is not necessary for an eventual stability
result. The proposed weaker ALT target explicitly assumes that, for each
admitted history h and allowed initial state, a finite time T(h,x0) exists at
which the actual trajectory enters a specified retained tail domain. History
includes the actual physical source, sensor disturbances and control/events.
There is no assertion that sup_h,x0 T(h,x0) is finite.

The conclusion below is **conditional history-wise eventual practical
boundedness**, with convergence in the zero-supply case. It is weaker than
uniform ISS, uniform exponential stability, a uniform capture deadline, or a
uniform accuracy guarantee. It is a meaningful theorem if the tail premises
are independently established. Assuming entry into an unspecified "stable
trajectory" would merely restate the desired conclusion and is not allowed.
The tail domain and storage law must be specified before the capture premise.

Existence and calculation are different: one may prove that each T is finite
without calculating it. Here finite capture is instead explicitly assumed at
the user's selected scope. Neither the old missed deadline nor the finite
antialigned witness disproves that assumption; neither proves it either.

## Conditional theorem A: one strict tail factor per history

Fix one history and initial condition satisfying capture. After capture index
service-word boundaries by j=0,1,..., at finite times t_j tending to infinity.
The same unchanged shipping trajectory and all carried memory continue through
these boundaries. No state/covariance/tuner or Live origin is reseeded.

Assume the following independently checkable tail premises, with constants
that may depend on this history and initial condition but are fixed along its
entire tail:

1. All required maps are well posed and the trajectory is forward complete.
   Capture occurs in a predeclared domain; word endpoints and all intervening
   prefixes remain in the certified domain. No infinite accumulation of events
   in finite time is allowed.
2. A declared storage law, such as M(P)=diag(P21^-1,I3), gives nonnegative
   V_j, with 0<m_h<=M_h<infinity and
   m_h ||e_j||² <= V_j <= M_h ||e_j||² on the tail. The same m_h must work
   for every j; pointwise positivity alone is insufficient.
3. For every tail word and its actual finite error, not merely its tangent,
   V_(j+1) <= r_h V_j + c_h D_h², where 0<=r_h<1,
   0<=c_h<infinity and D_h is a finite bound on the declared disturbance/supply.
   Canonical supplied bias coordinates and physical forcing must be accounted
   for; zero sensor noise alone does not imply D_h=0 for a forced wave model.
4. At all intermediate times t in [t_j,t_(j+1)), a storage W(t) satisfies
   m_h ||e(t)||² <= W(t) <= K_h V_j + ell_h D_h²,
   with finite K_h>=1 and ell_h>=0, including all hybrid jumps.

Then

V_j <= r_h^j V_0 + c_h D_h² (1-r_h^j)/(1-r_h),

and

limsup_(t->infinity) ||e(t)||²
  <= [K_h c_h/(1-r_h) + ell_h] D_h² / m_h.

In particular, if D_h=0 then e(t)->0. If the actual finite pre-capture prefix
is bounded, the entire trajectory is bounded, with a possibly history-dependent
bound. No common rate, common ultimate radius, or maximum T is asserted.

**Proof.** Repeated substitution proves the first inequality by induction:
the homogeneous term gains one factor r_h and the forced terms form a finite
geometric sum. Since r_h<1, r_h^j->0. Apply premise 4, divide by m_h>0, and
use t_j->infinity to obtain the limsup statement. The union of a bounded
pre-capture prefix and a bounded tail is bounded. The proof uses only finiteness
of capture; it never uses its numerical value or a supremum over histories.

Premise 1 must not be assumed from a word inequality valid only inside a domain.
For a domain containing a whole storage sublevel, a useful retention check is
V_0<=R_h and c_h D_h²<=(1-r_h)R_h: the boundary recurrence then preserves
V_j<=R_h. Prefix retention and arithmetic totality still require their own
checks. This avoids circularly applying a local estimate after leaving its
validity region.

The constants may exist without supplied numerical values, but their existence
still needs a mathematical proof. "The measured words are below one" proves
neither a tail supremum below one nor history-wise coercivity.

## Conditional theorem B: genuinely varying factors

More generally let V_(j+1)<=r_j V_j+q_j with r_j,q_j>=0, and define
P(n,k)=product_(i=k)^(n-1) r_i, with P(n,n)=1. Exact substitution gives

V_n <= P(n,0)V_0 + sum_(k=0)^(n-1) P(n,k+1)q_k.

This identity even allows some r_j>=1. It shows the two separate obligations:
forget the initial condition, and control the accumulated forcing.

For q_k<=c_h D_h², sufficient tail conditions are

P(n,0)->0, and G_h = sup_n sum_(k=0)^(n-1) P(n,k+1) < infinity.

They give limsup V_n<=c_h G_h D_h², and premise 4 of theorem A gives the
corresponding state bound. A simple sufficient certificate is
P(n,k)<=C_h lambda_h^(n-k) for every n>=k, with C_h finite and 0<lambda_h<1.
Then G_h<=C_h/(1-lambda_h). These constants need not be uniform over histories.

Alternatively, if 0<=r_j<=1 and the *proved* supply satisfies
q_j<=(1-r_j)B_h, telescoping gives

V_n <= P(n,0)V_0 + [1-P(n,0)]B_h.

Thus P(n,0)->0 suffices for limsup V_n<=B_h under this stronger supply law.
For example r_j=(j+1)/(j+2) has no strict upper bound below one, yet its product
is 1/(n+1). Divergence of sum_j (1-r_j) is sufficient for product decay,
using log r<=-(1-r); a zero factor also annihilates the initial product.
A fixed nonzero arithmetic/noise floor cannot silently be replaced by a
supply shrinking with (1-r_j). That replacement is a separate proof obligation.

**Proof.** The displayed product-sum bound follows by induction. Bound each
q_k by c_h D_h² for the first corollary. For the second, use the exact identity
sum P(n,k+1)(1-r_k)=1-P(n,0). The same prefix/coercivity argument as theorem A
then converts storage bounds to state bounds. These are proofs of abstract
implications, not assertions that the shipping filter meets their premises.

## Exact counterexamples to insufficient weakenings

1. **Every factor below one is insufficient.** With
   r_j=1-1/(j+2)² and q_j=0, the product through n words is
   (n+2)/(2(n+1))->1/2. Strict decrease at every word does not force convergence.
2. **Product decay alone is insufficient under persistent noise.** With
   r_j=(j+1)/(j+2), q_j=1 and V_0=0, exact induction gives
   V_n=n(n+3)/(2(n+1)), which diverges although P(n,0)=1/(n+1)->0.
3. **Storage decay without tail coercivity is insufficient.** Let e_j=1 and
   M_j=2^(-j). Then V_j=2^(-j) halves at every word while the state error does
   not decay. This is precisely the logical risk highlighted by the native
   no-magnetic-service inverse-covariance diagnostic.
4. **Finite capture alone does not bound startup amplification uniformly.**
   Consider histories indexed by N with x_(j+1)=2x_j for j<N and x_(j+1)=x_j/2
   thereafter. Each has finite capture time N and a contracting tail. Taking
   x_0=2^(-N) gives x_N=1, so no neighborhood of zero controls the transient
   uniformly over this family. Each fixed N has a finite but different gain.
   To claim Lyapunov stability from startup, one must prove suitable transient
   control and constants valid across the relevant nearby initial conditions.

`nonuniform_capture_theorem.py` and its tests verify exact rational finite
identities and these comparison counterexamples. The infinite conclusions
follow from the proofs above, not from extrapolation of finite tests.

## What changes for OU-III

The conditional target may now assume finite history-dependent capture, and
may seek r_h<1 and m_h>0 on each tail instead of one global capture deadline,
one global rho and one covariance enclosure covering every history. This
removes uniformity requirements from this weaker target. It does not discard
the remaining existence proofs or modify the original uniform proof track.

The captured domain must include actual Live/heading/service ancestry and the
finite-error region where the word estimate is valid. Informative service
alone does not imply capture, retention, or contraction. Without heading
service, use the explicit transverse/centre target, not full-heading decay.
Axial bias is still represented, and BIAS0/1/2 and all H18/A21 edges survive.

For the actual binary32 covariance law, `ou3-alt-binary32-coercivity.md`
now proves a format-wide coercivity bound on the finite symmetric SPD domain.
Thus the qualitative coercivity requirement reduces to preservation of that
domain, already part of the outstanding machine-arithmetic obligations. It no
longer independently requires a useful history-specific covariance envelope.
The existence lemma gives no useful numerical accuracy bound.

What is proved here is the conditional comparison theorem. The shipping
finite-error word inequality, history-wise tail coercivity and prefix
retention/totality remain open. Capture is an assumption, not a newly proved
startup result. All existing ALT PASS and storage-search gates remain false.
A future certificate can report this weaker conditional theorem separately;
it must not relabel it as the original uniform ISS or P4/P5 result.

For a zero-supply asymptotic-stability claim rather than eventual convergence,
local stability must also be established. For a useful noisy estimator result,
existence of an arbitrarily large history-specific radius is mathematically
meaningful but does not guarantee any specified physical accuracy.
