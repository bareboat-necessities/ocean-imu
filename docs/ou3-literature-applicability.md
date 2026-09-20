# OU-III literature applicability and realized-observer contract

This map belongs to the single proof in PR #560. A cited theorem does not
discharge a shipping hypothesis. The numbering below refers to the linked
primary-source versions, so it remains unambiguous across publication revisions.

| Exact result | Required hypotheses and conclusion | Shipping correspondence | Established and remaining |
|---|---|---|---|
| [Bonnabel–Slotine, arXiv:1211.6624v2](https://arxiv.org/pdf/1211.6624v2), Theorem 1 and Corollary 2; Assumptions 1–2, virtual system (5) | Continuous-time smooth EKF; uniformly positive and bounded P, positive Q; a local differential inequality, with Hessian bounds for an explicit radius. Gives local virtual-observer contraction. | Use the realized filter's gain as a coefficient of an auxiliary observer. Covariance energy is compatible with this viewpoint. | The recurring regular-root lower covariance is closed. Upper covariance, discrete/hybrid remainder, retained region and entry are open. The continuous-time conclusion is not imported into a discrete Joseph/reset word. |
| [Barrau–Bonnabel, arXiv:1410.1465v4](https://arxiv.org/pdf/1410.1465v4), Theorem 3, conditions (i)–(v), presenting the Deyst–Price linear criterion | Nonsingular bounded transition, positive observation noise and driving covariance, and full finite-window controllability/observability comparisons. The sampled-data linear result gives bounded covariance and decay in covariance energy. | Group the realized predictions and actual applied observations; preserve all state columns. Uniform full Gramians and reset/PSD-sync transfer must be checked for the implemented discrete coefficients. | Positive process mechanisms and regular-root lower covariance are available. Full observability, upper covariance and uniform hybrid transfer remain open. A two-column magnetic restriction is insufficient for condition (v). |
| Same source, Theorems 1–2 and 4 | Special group dynamics/output structure yields autonomous and log-linear invariant errors; Theorem 4 also needs Theorem 3's hypotheses. | A quaternion representation alone supplies none of those structural identities. | No equivalence of shipping OU-III to this IEKF has been proved. Its nonlinear stability theorem is not applied. |
| [D'Souza–Zanetti, author manuscript](https://www.researchgate.net/publication/325301963_Information_Formulation_of_the_UDU_Kalman_Filter), equations (1)–(6) and Section III; [published DOI](https://doi.org/10.1109/TAES.2018.2850379) | Factorized information calculations allow vector observations with non-diagonal noise. This is a numerical representation result. | Shipping three-dimensional innovation solves and the proof's retained measurement factors use the same small observation dimension. | Exact real-arithmetic equivalence is checked independently. Numerical factor stability does not imply observer convergence or bound the shipping float32 Joseph error. The new proof engine uses orthogonal square-root arrays, not an unimplemented UDU estimator. |
| [Aslam–Haydar–Akhtar, arXiv:2503.08255](https://arxiv.org/pdf/2503.08255), Theorem 1, Assumption 1, equations (36)–(37) | Continuous generalized SO(3)-MEKF with curvature-dependent terms, bounded positive P, exact matrix attitude measurements, zero process disturbance, and exclusion of the pi-error set. | Those modified gain/gain-evolution equations and measurement model are not the shipping observer. | The almost-global conclusion is inapplicable to unchanged OU-III and is not used for capture. |

The primary [Zhang–Zhang covariance-boundedness paper](https://doi.org/10.1016/j.ifacol.2021.08.381)
is identified for a discrete Riccati-bounds audit. Its full text was not
retrievable from the publisher/HAL during this review; no theorem number,
bound or discharged hypothesis is inferred from its abstract or search excerpts.
The same rule applies to any newer differentiation/factorization paper: it must
resolve a named hypothesis before being used as stability evidence.

## What is fixed in the auxiliary observer

Let zeta_k denote the complete *realized* shipping record: nominal estimate,
covariance, gain, measurement, applied/rejected event, tuner parameters,
reference state, projection/reset operation and clocks. For that realization,
define the auxiliary operation using its realized gain and event decision.
For a Euclidean correction the map is

`z^+ = z^- + K_k [y_k - h_k(z^-)]`.

The realized estimate is one trajectory of this map. Physical truth is a
forced trajectory, with measurement/model mismatch explicitly retained. A
second independently running EKF is not this auxiliary trajectory: its gains,
tuner and gates generally differ. No assertion equates those executions.

For quaternion injection, bias projection and other hybrid operations, an
auxiliary map must reproduce the actual retraction/projection, and its error
must be expressed in the same reset chart. Defining a frozen gain does not
prove that the covariance reset is the exact finite-error differential.
That discrepancy, sensor forcing and model mismatch must be included in

`e_(k+1) = A_k(zeta_k) e_k + r_k(e_k,zeta_k,d_k)`.

The linear comparison uses the actual ordered F, I-KH and congruent reset G.
The open nonlinear task is to bound the complete composed r in that storage
uniformly over admitted realized records. It does not require differentiating
the gains of every independently perturbed covariance trajectory. Event
discontinuities, projection sectors and retention still require their own
finite-error argument; they cannot be declared smooth by freezing the record.

## Exact discrete result currently used

The proof derives, rather than imports from a continuous theorem,

`M_W^T P_W^-1 M_W + D_W = P_0^-1`.

Its hypotheses are finite same-dimensional real-arithmetic operations,
P_0 positive definite, Q positive semidefinite, effective R positive definite,
the optimal gain for the actual innovation matrix, and nonsingular congruent
resets. No full-state observability premise is needed for this identity.

If the shipping safety path adds E>=0 to S, set R_eff=R+E. The shipping
gain and expanded Joseph covariance both use that same adjusted S. The
identity therefore applies to that real-arithmetic comparison. A failed
factorization, nonoptimal/frozen gain row, noncongruent projection or arithmetic
defect is not silently admitted as a regular A21 correction.

The required next implication is a *uniform* lower singular bound for the
full root-whitened loss factor, including all nuisance columns, followed by
the same-storage nonlinear/supply bound and capture-to-retention composition.
Neither a cited covariance theorem with unchecked hypotheses nor the
positive quiet-water diagnostic is a substitute.
