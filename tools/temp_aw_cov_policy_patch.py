from pathlib import Path


def replace_once(path, old, new):
    p = Path(path)
    text = p.read_text()
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"{path}: expected one match, found {n}")
    p.write_text(text.replace(old, new, 1))


def replace_between(path, start, end, new):
    p = Path(path)
    text = p.read_text()
    i = text.find(start)
    if i < 0:
        raise SystemExit(f"{path}: start marker not found: {start}")
    j = text.find(end, i)
    if j < 0:
        raise SystemExit(f"{path}: end marker not found: {end}")
    p.write_text(text[:i] + new + text[j:])


# OU-II: prediction-time PSD covariance inflation by default; one flag restores
# the historical immediate marginal replacement.
p = "src/kalman_ou_ii/Kalman3D_Wave_OU_II.h"
old = r'''    // Synchronize only the a_w marginal with the stationary process model while
    // preserving cross-covariances learned by the running filter.
    void synchronize_aw_covariance_to_stationary() {
        Pext.template block<3,3>(OFF_AW, OFF_AW) = Sigma_aw_stat;
        symmetrize_Pext_();
    }
'''
new = r'''    // Select the pre-PSD-inflation covariance policy. false (default) keeps
    // synchronization inside the next Kalman prediction as a positive-
    // semidefinite process-covariance increment. true restores the historical
    // immediate P_awaw block replacement exactly.
    void set_legacy_aw_covariance_replacement(bool on) {
        legacy_aw_covariance_replacement_ = on;
        if (on) aw_covariance_floor_pending_ = false;
    }
    [[nodiscard]] bool legacy_aw_covariance_replacement() const noexcept {
        return legacy_aw_covariance_replacement_;
    }

    // Request an a_w covariance floor at the current stationary OU covariance.
    //
    // Default policy: queue the request and apply it inside the next prediction
    // as
    //   Delta = Pi_+(Sigma_aw_stat - P_awaw^-),
    //   P^-   <- P^- + E_a Delta E_a^T.
    // This can only add covariance, preserves PSD, leaves cross-covariances
    // untouched, and is algebraically an additional bounded process-noise term.
    //
    // Legacy policy: replace P_awaw immediately while retaining the old cross-
    // covariances. This is kept only for regression/rollback comparisons.
    void synchronize_aw_covariance_to_stationary() {
        if (legacy_aw_covariance_replacement_) {
            aw_covariance_floor_pending_ = false;
            Pext.template block<3,3>(OFF_AW, OFF_AW) = Sigma_aw_stat;
            symmetrize_Pext_();
            return;
        }
        aw_covariance_floor_target_ = T(0.5) * (Sigma_aw_stat + Sigma_aw_stat.transpose());
        aw_covariance_floor_pending_ = true;
    }
'''
replace_once(p, old, new)
replace_once(
    p,
    '''    bool linear_block_enabled_ = true;\n    bool acc_bias_updates_enabled_ = true;\n\n    bool aw_process_correlated_ = false;\n''',
    '''    bool linear_block_enabled_ = true;\n    bool acc_bias_updates_enabled_ = true;\n\n    // false is the default/proof-compatible policy. true restores the old\n    // immediate marginal replacement for exact regression comparisons.\n    bool legacy_aw_covariance_replacement_ = false;\n    bool aw_covariance_floor_pending_ = false;\n    Matrix3 aw_covariance_floor_target_ = Matrix3::Zero();\n\n    bool aw_process_correlated_ = false;\n''',
)
replace_once(
    p,
    '''        Pext.template block<3,3>(OFF_AW, OFF_AW) = Sigma_aw_stat;\n        symmetrize_Pext_();\n    }\n\n    // Select the pre-PSD-inflation covariance policy.''',
    '''        Pext.template block<3,3>(OFF_AW, OFF_AW) = Sigma_aw_stat;\n        aw_covariance_floor_pending_ = false;\n        symmetrize_Pext_();\n    }\n\n    // Select the pre-PSD-inflation covariance policy.''',
)
marker = '''    EIGEN_STRONG_INLINE void symmetrize_Pext_() {\n        for (int i = 0; i < NX; ++i) {\n            for (int j = i + 1; j < NX; ++j) {\n                const T v = T(0.5) * (Pext(i,j) + Pext(j,i));\n                Pext(i,j) = v;\n                Pext(j,i) = v;\n            }\n        }\n    }\n'''
helper = marker + r'''

    // Apply a queued covariance floor as part of the prediction covariance.
    // The positive spectral projection is the Frobenius-nearest PSD matrix to
    // Sigma_aw_stat - P_awaw^-. Negative directions are deliberately ignored:
    // uncertainty may be added here, never deleted.
    EIGEN_STRONG_INLINE void apply_pending_aw_covariance_inflation_() {
        if (!aw_covariance_floor_pending_ || !linear_block_enabled_) return;
        aw_covariance_floor_pending_ = false;

        Matrix3 P_aw = Pext.template block<3,3>(OFF_AW, OFF_AW);
        P_aw = T(0.5) * (P_aw + P_aw.transpose());
        Matrix3 Delta = aw_covariance_floor_target_ - P_aw;
        Delta = T(0.5) * (Delta + Delta.transpose());

        Eigen::SelfAdjointEigenSolver<Matrix3> es(Delta);
        if (es.info() != Eigen::Success) return;
        Vector3 evals = es.eigenvalues();
        for (int i = 0; i < 3; ++i) evals(i) = std::max(T(0), evals(i));
        Delta = es.eigenvectors() * evals.asDiagonal() * es.eigenvectors().transpose();
        Delta = T(0.5) * (Delta + Delta.transpose());

        Pext.template block<3,3>(OFF_AW, OFF_AW) += Delta;
    }
'''
replace_once(p, marker, helper)
replace_once(
    p,
    '''    symmetrize_Pext_();   // Symmetry hygiene\n\n    if (linear_block_enabled_ &&\n''',
    '''    if (linear_block_enabled_) {\n        apply_pending_aw_covariance_inflation_();\n    }\n    symmetrize_Pext_();   // Symmetry hygiene\n\n    if (linear_block_enabled_ &&\n''',
)

# OU-III: same policy, retaining the separate congruence experiment.
p = "src/kalman_ou_iii/Kalman3D_Wave_OU_III.h"
old = r'''    // Synchronize only the a_w marginal with the stationary process model,
    // leaving the cross-covariances the running filter has learned untouched.
    //
    // This is the deployed policy. It is not a consistent posterior operation:
    // the implied correlation coefficients between a_w and the rest of the
    // state are rescaled by the square root of the marginal change, which in a
    // developed sea is a factor of ten or more. See
    // synchronize_aw_covariance_to_stationary_congruent() for the consistent
    // alternative and docs/ou-iii-adaptation-and-travel-sense.md for what each
    // one costs.
    void synchronize_aw_covariance_to_stationary() {
        Pext.template block<3,3>(OFF_AW, OFF_AW) = Sigma_aw_stat;
        symmetrize_Pext_();
    }
'''
new = r'''    // Select the pre-PSD-inflation covariance policy. false (default) keeps
    // synchronization inside the next Kalman prediction as a positive-
    // semidefinite process-covariance increment. true restores the historical
    // immediate P_awaw block replacement exactly.
    void set_legacy_aw_covariance_replacement(bool on) {
        legacy_aw_covariance_replacement_ = on;
        if (on) aw_covariance_floor_pending_ = false;
    }
    [[nodiscard]] bool legacy_aw_covariance_replacement() const noexcept {
        return legacy_aw_covariance_replacement_;
    }

    // Request an a_w covariance floor at the current stationary OU covariance.
    // The default path queues the request and applies it inside the next
    // prediction as Delta = Pi_+(Sigma_aw_stat - P_awaw^-). Because the full
    // covariance increment is E_a Delta E_a^T >= 0, PSD and all existing cross-
    // covariances are preserved. The legacy flag restores the previous raw
    // posterior block replacement for exact rollback/regression comparisons.
    void synchronize_aw_covariance_to_stationary() {
        if (legacy_aw_covariance_replacement_) {
            aw_covariance_floor_pending_ = false;
            Pext.template block<3,3>(OFF_AW, OFF_AW) = Sigma_aw_stat;
            symmetrize_Pext_();
            return;
        }
        aw_covariance_floor_target_ = T(0.5) * (Sigma_aw_stat + Sigma_aw_stat.transpose());
        aw_covariance_floor_pending_ = true;
    }
'''
replace_once(p, old, new)
replace_once(
    p,
    '''    void synchronize_aw_covariance_to_stationary_congruent() {\n        Matrix3 P_aw = Pext.template block<3,3>(OFF_AW, OFF_AW);\n''',
    '''    void synchronize_aw_covariance_to_stationary_congruent() {\n        aw_covariance_floor_pending_ = false;\n        Matrix3 P_aw = Pext.template block<3,3>(OFF_AW, OFF_AW);\n''',
)
replace_once(
    p,
    '''    bool linear_block_enabled_ = true;\n    bool acc_bias_updates_enabled_ = true;\n\n    bool aw_process_correlated_ = false;\n''',
    '''    bool linear_block_enabled_ = true;\n    bool acc_bias_updates_enabled_ = true;\n\n    // false is the default/proof-compatible policy. true restores the old\n    // immediate marginal replacement for exact regression comparisons.\n    bool legacy_aw_covariance_replacement_ = false;\n    bool aw_covariance_floor_pending_ = false;\n    Matrix3 aw_covariance_floor_target_ = Matrix3::Zero();\n\n    bool aw_process_correlated_ = false;\n''',
)
replace_once(
    p,
    '''        Pext.template block<3,3>(OFF_AW, OFF_AW) = Sigma_aw_stat;\n        symmetrize_Pext_();\n    }\n\n    // Select the pre-PSD-inflation covariance policy.''',
    '''        Pext.template block<3,3>(OFF_AW, OFF_AW) = Sigma_aw_stat;\n        aw_covariance_floor_pending_ = false;\n        symmetrize_Pext_();\n    }\n\n    // Select the pre-PSD-inflation covariance policy.''',
)
marker = '''    EIGEN_STRONG_INLINE void symmetrize_Pext_() {\n        for (int i = 0; i < NX; ++i) {\n            for (int j = i + 1; j < NX; ++j) {\n                const T v = T(0.5) * (Pext(i,j) + Pext(j,i));\n                Pext(i,j) = v;\n                Pext(j,i) = v;\n            }\n        }\n    }\n'''
helper = marker + r'''

    EIGEN_STRONG_INLINE void apply_pending_aw_covariance_inflation_() {
        if (!aw_covariance_floor_pending_ || !linear_block_enabled_) return;
        aw_covariance_floor_pending_ = false;

        Matrix3 P_aw = Pext.template block<3,3>(OFF_AW, OFF_AW);
        P_aw = T(0.5) * (P_aw + P_aw.transpose());
        Matrix3 Delta = aw_covariance_floor_target_ - P_aw;
        Delta = T(0.5) * (Delta + Delta.transpose());

        Eigen::SelfAdjointEigenSolver<Matrix3> es(Delta);
        if (es.info() != Eigen::Success) return;
        Vector3 evals = es.eigenvalues();
        for (int i = 0; i < 3; ++i) evals(i) = std::max(T(0), evals(i));
        Delta = es.eigenvectors() * evals.asDiagonal() * es.eigenvectors().transpose();
        Delta = T(0.5) * (Delta + Delta.transpose());

        Pext.template block<3,3>(OFF_AW, OFF_AW) += Delta;
    }
'''
replace_once(p, marker, helper)
replace_once(
    p,
    '''    symmetrize_Pext_();   // Symmetry hygiene\n\n    if (linear_block_enabled_ &&\n''',
    '''    if (linear_block_enabled_) {\n        apply_pending_aw_covariance_inflation_();\n    }\n    symmetrize_Pext_();   // Symmetry hygiene\n\n    if (linear_block_enabled_ &&\n''',
)

# Stability appendix: make the new default part of the UES theorem and clearly
# separate the legacy raw replacement.
p = "doc/kalman_ou_iii/w3d-iss-stability.tex-part"
theorem = r'''\begin{theorem}[UES of the 18-state Live configuration with PSD covariance inflation]
\label{thm:iss-18-ues}
Consider the reduced Live filter with accelerometer-bias estimation disabled,
no hard reset inside the analysis interval, and the normal Kalman
predict/update Riccati recursion. At an optional covariance-maintenance event,
let the implemented default policy augment the predicted covariance by
\begin{equation}
\mat P_k^-\leftarrow \mat P_k^-
+\mat E_a\mat\Delta_k\mat E_a^\top,
\qquad
\mat\Delta_k
=\Pi_+\!\left(
\mat\Sigma_{aw,k}^{\mathrm{stat}}-\mat P_{aa,k}^-
\right),
\label{eq:iss-psd-inflation}
\end{equation}
where $\mat E_a$ selects the three $\vct a_w$ coordinates and $\Pi_+$ replaces
negative eigenvalues of a symmetric matrix by zero. Assume in addition to the
standing conditions that
\begin{equation}
\mat 0\preceq\mat\Sigma_{aw,k}^{\mathrm{stat}}
\preceq \bar\sigma_P\mat I_3<\infty
\label{eq:iss-Sigma-bound}
\end{equation}
uniformly. Then there exist finite constants
$0<p_-\le p_+<\infty$, $M\ge1$, and $0<\rho<1$ such that
\begin{equation}
p_-\mat I_{18}\preceq\mat P_k\preceq p_+\mat I_{18}
\label{eq:iss-P-bounds}
\end{equation}
and the homogeneous closed-loop estimation-error transition satisfies
\begin{equation}
\|\mat\Psi_{18}(k,j)\|
\le M\rho^{k-j},
\qquad k\ge j.
\label{eq:iss-ues}
\end{equation}
Thus the nominal linearized 18-state Live configuration using the default PSD
inflation policy is uniformly exponentially stable.
\end{theorem}

\begin{proof}
Lemma~\ref{lem:iss-18-uco} establishes uniform complete observability and
Eq.~\eqref{eq:iss-Q18-bounds} gives the positive bounded process covariance of
the ordinary IMU prediction. The new maintenance operation is not a covariance
replacement. It is exactly an additional process-covariance term
\begin{equation}
\mat Q_{\mathrm{infl},k}
=\mat E_a\mat\Delta_k\mat E_a^\top\succeq\mat0.
\label{eq:iss-Qinfl}
\end{equation}
Because $\mat P_{aa,k}^-\succeq\mat0$ and
Eq.~\eqref{eq:iss-Sigma-bound} holds, the positive part obeys
\begin{equation}
\mat0\preceq\mat\Delta_k\preceq
\bar\sigma_P\mat I_3 .
\label{eq:iss-Delta-bound}
\end{equation}
Consequently the effective prediction covariance satisfies
\begin{equation}
q_-\mat I_{18}
\preceq
\mat Q_k^{\mathrm{eff}}
=\mat Q_k+\mat Q_{\mathrm{infl},k}
\preceq
(q_++\bar\sigma_P)\mat I_{18}.
\label{eq:iss-Qeff-bounds}
\end{equation}
The inflation therefore preserves the uniform positive lower process-noise
bound and only changes its finite upper bound. It does not alter the state
transition or any observation matrix, so the UCO result is unchanged.
Standard discrete time-varying Kalman-filter results then give uniform upper and
lower Riccati bounds and UES of the homogeneous filter error dynamics
\cite{Jazwinski1970_Filtering,Maybeck1979_StochasticModels}. The argument is
uniform over every realized sequence of bounded PSD increments satisfying
Eq.~\eqref{eq:iss-Delta-bound}; the fact that the implementation computes the
increment from the current predicted marginal therefore does not weaken the
bound.

The MEKF covariance reset after a small quaternion correction does not alter
this first-order conclusion: its reset Jacobian equals the identity at zero
attitude error, so the linearized homogeneous map is the same map used above.
\end{proof}

The theorem covers both the new default PSD-inflation maintenance and the case
where periodic maintenance is disabled entirely. It does not cover the optional
legacy raw block-replacement mode, and it does not claim UES of the unrestricted
21-state estimator with active random-walk accelerometer-bias estimation.

'''
replace_between(
    p,
    r'\begin{theorem}[UES of the standard-Riccati 18-state Live proof configuration]',
    r'\subsection{Why unconditional 21-state UES is false}',
    theorem,
)

covsec = r'''\subsection{Default PSD covariance inflation and the optional legacy replacement}
\label{app:iss-cov-sync}

The covariance-maintenance purpose is to prevent the latent-acceleration
marginal from becoming unrealistically overconfident when the stationary OU
scale changes. The previous implementation achieved that by replacing the
posterior $\vct a_w$ block while retaining learned cross-covariances. The new
default instead requests the floor at the same adaptation cadence and applies
it during the next Kalman prediction through
Eq.~\eqref{eq:iss-psd-inflation}. Only covariance is added: directions in which
$\mat P_{aa}^-$ already exceeds the stationary target receive zero increment.
This is why the change was made. It preserves positive semidefiniteness by
construction, retains the learned cross-covariances, and keeps the operation
inside the bounded process-covariance class used by
Theorem~\ref{thm:iss-18-ues}.

The old behavior remains available through the optional runtime flag
\texttt{set\_legacy\_aw\_covariance\_replacement(true)}. Its default is
\texttt{false}. With the flag enabled, synchronization again performs the
historical immediate replacement
\begin{equation}
\mat P^+=
\begin{bmatrix}
\mat A&\mat C\\
\mat C^\top&\mat\Sigma_{aw}^{\mathrm{stat}}
\end{bmatrix}
\label{eq:iss-P-raw-sync}
\end{equation}
from a pre-sync partition
\begin{equation}
\mat P=
\begin{bmatrix}
\mat A&\mat C\\
\mat C^\top&\mat D
\end{bmatrix}\succeq0,
\qquad \mat A\succ0.
\label{eq:iss-P-partition}
\end{equation}
The Schur complement requires
\begin{equation}
\mat\Sigma_{aw}^{\mathrm{stat}}
-\mat C^\top\mat A^{-1}\mat C\succeq0
\label{eq:iss-sync-schur}
\end{equation}
for the replacement to preserve PSD. The pre-sync covariance only guarantees
$\mat D-\mat C^\top\mat A^{-1}\mat C\succeq0$, so a downward change in the
stationary target need not satisfy Eq.~\eqref{eq:iss-sync-schur}. The article
therefore makes no UES claim for this legacy mode; it is retained so previous
results can be reproduced exactly and the policy can be rolled back with one
flag.

OU--III also retains the separate congruence re-alignment experiment. That map
preserves PSD, but it is an additional covariance coordinate transformation and
is not the default policy or part of the theorem above. Setting periodic
covariance synchronization off altogether removes both maintenance policies and
returns to the ordinary Riccati recursion, which is also covered by the theorem
when its other conditions hold.

'''
replace_between(
    p,
    r'\subsection{Periodic covariance re-alignment remains outside the theorem}',
    r'\subsection{Scope, remaining conditions, and proof-oriented implementation choices}',
    covsec,
)
replace_once(
    p,
    r'''\paragraph{What is established.}
The OU--III translational block has an explicit uniform four-update
observability bound for arbitrary bounded measurable $\tau(t)$. Under bounded
noncollinear accelerometer/magnetometer geometry and bounded update gaps, the
attitude--gyro-bias block has a uniform finite-horizon observability bound. In
the reduced 18-state standard-Riccati configuration, these combine with
uniformly positive process covariance to give UCO, bounded Riccati solutions,
and UES by standard time-varying Kalman theory.
''',
    r'''\paragraph{What is established.}
The OU--III translational block has an explicit uniform four-update
observability bound for arbitrary bounded measurable $\tau(t)$. Under bounded
noncollinear accelerometer/magnetometer geometry and bounded update gaps, the
attitude--gyro-bias block has a uniform finite-horizon observability bound. In
the reduced 18-state Live configuration these combine with uniformly positive
process covariance and the bounded default PSD $\vct a_w$ inflation to give UCO,
bounded Riccati solutions, and UES by standard time-varying Kalman theory.
''',
)
replace_once(
    p,
    r'''\paragraph{What remains outside the theorem.}
The default periodic $\vct a_w$ covariance synchronization, active
accelerometer-bias estimation, hard initialization, tilt re-lock, gating and
fault recovery, and arbitrary loss of magnetometer updates are not covered by
Theorem~\ref{thm:iss-18-ues}. The physical $S=0$ pseudo-measurement is a
regularizer rather than a true observation, so stochastic physical mismatch
must be handled separately from nominal filter UES.
''',
    r'''\paragraph{What remains outside the theorem.}
The optional legacy raw $\vct a_w$ covariance replacement, active
accelerometer-bias estimation, hard initialization, tilt re-lock, gating and
fault recovery, and arbitrary loss of magnetometer updates are not covered by
Theorem~\ref{thm:iss-18-ues}. The physical $S=0$ pseudo-measurement is a
regularizer rather than a true observation, so stochastic physical mismatch
must be handled separately from nominal filter UES.
''',
)
replace_once(
    p,
    r'''\item \textbf{Disable periodic raw $\vct a_w$ covariance synchronization.}
This restores the standard Riccati recursion used by
Theorem~\ref{thm:iss-18-ues}. If synchronization is retained, prefer the
positive-semidefinite congruence form and prove uniform bounds for the resulting
periodic covariance map.
''',
    r'''\item \textbf{Use the default PSD $\vct a_w$ covariance inflation.}
The periodic maintenance call now queues
$\mat\Delta=\Pi_+(\mat\Sigma_{aw}^{\mathrm{stat}}-\mat P_{aa}^-)$ and adds it
inside the next prediction. This is the default and is covered by
Theorem~\ref{thm:iss-18-ues}. Set
\texttt{set\_legacy\_aw\_covariance\_replacement(true)} only to reproduce the
previous immediate block replacement; UES is not proved here for that optional
legacy behavior.
''',
)
replace_once(
    p,
    r'''The resulting hierarchy is therefore precise: the translational and
attitude--gyro-bias observability pieces can be closed analytically; a genuine
18-state UES theorem follows for the standard-Riccati reduced Live
configuration; the corresponding local ISS-type perturbation bound then follows
for bounded residual bias and model mismatch; but unconditional 21-state UES is
ruled out by the static attitude--accelerometer-bias ambiguity rather than left
as an unspecified future proof obligation.''',
    r'''The resulting hierarchy is therefore precise: the translational and
attitude--gyro-bias observability pieces can be closed analytically; the new
default bounded PSD covariance inflation remains inside the process-covariance
class used by the Riccati argument; and a genuine 18-state UES theorem follows
for the reduced Live configuration under the stated geometry, cadence, noise,
and no-reset conditions. The corresponding local ISS-type perturbation bound
then follows for bounded residual bias and model mismatch. The optional legacy
raw covariance replacement is deliberately retained but is not covered by that
proof, while unconditional 21-state UES is ruled out by the static
attitude--accelerometer-bias ambiguity.''',
)

# Introduction: synchronize claim summary.
p = "doc/kalman_ou_iii/w3d-intro.tex-part"
old = r'''Appendix~\ref{app:iss-stability} is intentionally limited to analytical
observability and stability conditions; it is not presented as an ISS proof of
the implemented 21-state estimator. Because the tuner is driven by raw IMU data
and a private Mahony attitude solution, its schedule is exogenous with respect
to the MEKF state. The appendix proves a uniform finite-horizon observability
bound for each OU--III translational chain under arbitrary bounded variation of
$\tau$, and shows how $r_S\propto\sigma_{aw}\tau^3$ preserves the dimensionless
strength of the integral pseudo-measurement. It also records a standard
conditional perturbation lemma: if uniform exponential stability of the
complete Live linearization is established independently, then sufficiently
small nonlinear remainder and bounded disturbances give a local ISS-type
bound. The appendix does \emph{not} establish that required full-system UES
property. A coupled attitude--gyro-bias--translation--accelerometer-bias
finite-horizon bound, uniform scheduled stabilizability, bounded stabilizing
Riccati behavior, treatment of covariance re-alignment, the stochastic meaning
of the $S=0$ regularizer, and hybrid reset logic all remain open. Consequently,
this is not a conditional sufficient-condition result for the implemented
filter and not an unconditional EKF convergence claim.
'''
new = r'''Appendix~\ref{app:iss-stability} gives an analytical stability result with a
carefully delimited state and operating region; it is not an ISS or convergence
proof of the unrestricted 21-state estimator. Because the tuner is driven by
raw IMU data and a private Mahony attitude solution, its schedule is exogenous
with respect to the MEKF state. The appendix proves uniform finite-horizon
observability of the OU--III translational block, finite-horizon observability
of attitude and gyro bias under bounded noncollinear accel/mag geometry, and
uniform exponential stability of the resulting 18-state Live performance
filter when accelerometer-bias estimation is frozen, process and measurement
covariances have positive finite bounds, magnetic updates have bounded gaps, and
no hard reset occurs in the analysis interval. The default periodic
$\vct a_w$ covariance maintenance is changed to prediction-time PSD inflation,
so it is an additional bounded process-covariance term and is included in that
UES proof. A single optional legacy flag restores the previous raw marginal
replacement for reproducibility, but UES is not proved for that legacy mode.
Unconditional 21-state UES is in fact impossible over a motion class containing
static trajectories with the present random-walk accelerometer-bias state,
because the appendix constructs a nondecaying attitude--bias ambiguity. The
$S=0$ regularizer's physical model mismatch and hybrid reset/recovery logic also
remain outside the local nominal theorem.
'''
replace_once(p, old, new)

# Conclusion: synchronize claim boundaries.
p = "doc/kalman_ou_iii/w3d-conclusion-summary.tex-part"
replace_between(
    p,
    r'Appendix~\ref{app:iss-stability} provides a partial analytical observability and',
    r'The results support drift-controlled wave-state reconstruction for the tested',
    r'''Appendix~\ref{app:iss-stability} now closes the nominal reduced-state stability
chain more tightly than the empirical results alone. It derives an explicit
uniform finite-horizon observability bound for every OU--III translational axis,
proves finite-horizon attitude--gyro-bias observability under bounded
noncollinear accelerometer/magnetometer geometry and bounded magnetic-update
gaps, and combines those pieces with positive bounded process and measurement
covariances to obtain UES of the 18-state Live performance filter. The default
periodic $\vct a_w$ covariance maintenance is part of that theorem: instead of
replacing the posterior marginal, the implementation now adds the PSD positive
part of the stationary-covariance deficit inside the next prediction. This
preserves PSD and can be absorbed into a uniformly bounded effective process
covariance.

The theorem has explicit boundaries. Accelerometer-bias estimation is frozen
and its residual error is treated as an input; no hard reset occurs inside the
analysis interval; the reference vectors remain bounded and uniformly
noncollinear; accepted magnetometer updates have a bounded gap; and process and
measurement covariances retain positive finite bounds. With the present
random-walk accelerometer-bias augmentation, unconditional 21-state UES over a
motion class that contains static trajectories is not merely unproved but
false: a static attitude--accelerometer-bias ambiguity gives a nondecaying
unobservable mode. The $S=0$ channel also remains a stochastic model
regularizer rather than a physical observation, so the UES theorem is a nominal
linearized result rather than an infinite-horizon deterministic bound on ocean
model mismatch.

For reproducibility the prior covariance policy remains available as an
optional rollback. Setting
\texttt{set\_legacy\_aw\_covariance\_replacement(true)} restores the previous
immediate $P_{a_wa_w}\leftarrow\Sigma_{aw}^{\mathrm{stat}}$ block replacement
while retaining cross-covariances. The default is the new PSD inflation policy;
UES is proved for that default under the conditions above, but is not claimed in
this article for the optional legacy replacement. Disabling periodic
maintenance entirely also returns to the ordinary Riccati case covered by the
same reduced-state theorem when its other assumptions hold. Hybrid startup,
re-lock, gating, and recovery remain outside the local Live result.

''',
)

# Regression test for both OU generations.
test = r'''#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>

#include "kalman_ou_ii/Kalman3D_Wave_OU_II.h"
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"

namespace {

template <typename MatrixT>
bool psd(const MatrixT& P, double tol = 1e-9) {
    Eigen::MatrixXd S = 0.5 * (P.template cast<double>() + P.transpose().template cast<double>());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
    return es.info() == Eigen::Success && es.eigenvalues().minCoeff() >= -tol;
}

template <typename Filter>
bool exercise(Filter& f, int off_aw, const char* name) {
    if (f.legacy_aw_covariance_replacement()) {
        std::cerr << name << ": legacy replacement must default to false\n";
        return false;
    }

    const auto before = f.covariance_full();
    Eigen::Vector3d std_hi(5.0, 6.0, 7.0);
    f.set_aw_stationary_std(std_hi);
    f.synchronize_aw_covariance_to_stationary();

    const auto queued = f.covariance_full();
    if ((queued - before).norm() > 1e-12) {
        std::cerr << name << ": default synchronization rewrote posterior covariance\n";
        return false;
    }

    f.time_update(Eigen::Vector3d::Zero(), 0.005);
    const auto after = f.covariance_full();
    Eigen::Matrix3d target = std_hi.array().square().matrix().asDiagonal();
    Eigen::Matrix3d margin = after.template block<3,3>(off_aw, off_aw) - target;
    if (!psd(margin, 1e-8) || !psd(after, 1e-8)) {
        std::cerr << name << ": prediction-time PSD inflation failed floor/PSD test\n";
        return false;
    }

    f.set_legacy_aw_covariance_replacement(true);
    if (!f.legacy_aw_covariance_replacement()) {
        std::cerr << name << ": legacy flag did not latch\n";
        return false;
    }
    Eigen::Vector3d std_legacy(1.1, 1.2, 1.3);
    f.set_aw_stationary_std(std_legacy);
    f.synchronize_aw_covariance_to_stationary();
    const auto legacy = f.covariance_full();
    Eigen::Matrix3d legacy_target = std_legacy.array().square().matrix().asDiagonal();
    if ((legacy.template block<3,3>(off_aw, off_aw) - legacy_target).norm() > 1e-10) {
        std::cerr << name << ": legacy mode did not reproduce block replacement\n";
        return false;
    }
    return true;
}

} // namespace

int main() {
    const Eigen::Vector3d sigma_a = Eigen::Vector3d::Constant(0.02);
    const Eigen::Vector3d sigma_g = Eigen::Vector3d::Constant(0.001);
    const Eigen::Vector3d sigma_m = Eigen::Vector3d::Constant(0.1);

    Kalman3D_Wave_OU_II<double> ou2(sigma_a, sigma_g, sigma_m);
    Kalman3D_Wave_OU_III<double> ou3(sigma_a, sigma_g, sigma_m);

    if (!exercise(ou2, 12, "OU-II")) return 1;
    if (!exercise(ou3, 15, "OU-III")) return 1;

    std::cout << "aw covariance policy: PASS\n";
    return 0;
}
'''
Path("tests/kalman_ou_iii/aw_covariance_policy-test.cpp").write_text(test)

# Register the compile/runtime regression target.
p = "CMakeLists.txt"
replace_once(
    p,
    'add_executable(kalman_ou_common-test "${OCEAN_IMU_TESTS_DIR}/kalman_ou_iii/kalman_ou_common-test.cpp")\n',
    'add_executable(kalman_ou_common-test "${OCEAN_IMU_TESTS_DIR}/kalman_ou_iii/kalman_ou_common-test.cpp")\nadd_executable(aw_covariance_policy-test "${OCEAN_IMU_TESTS_DIR}/kalman_ou_iii/aw_covariance_policy-test.cpp")\n',
)
