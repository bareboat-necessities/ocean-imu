#!/usr/bin/env python3
from pathlib import Path

path = Path("doc/kalman_ou_iii/w3d-baseline-comparison-study.tex-part")
text = path.read_text(encoding="utf-8")

replacements = [
    (
        """OU--II is the immediate stochastic predecessor of OU--III. It retains the
quaternion error-state filtering, colored latent-acceleration prior, adaptive
sea-state tuning, and integral-state regularization, but uses the shorter
second-order translational chain. Its comparison with OU--III is therefore an
internal ablation of the additional integral-displacement state and associated
covariance structure rather than an independent external benchmark.""",
        """OU--II is the immediate stochastic predecessor of OU--III. It retains the
quaternion error-state filtering, colored latent-acceleration prior, adaptive
sea-state tuning, and low-frequency drift regularization, but uses the shorter
second-order translational chain with position/velocity pseudo-constraints
rather than the additional integral-displacement state. Its comparison with
OU--III is therefore an internal architectural ablation rather than an
independent external benchmark.""",
    ),
    (
        """For that ablation to isolate the state, the two filters have to be tuned from
the same place. Both therefore derive the operating point from the wave-band
zero-crossing period of Sec.~\\ref{sec:wave-band-period}, with
$\\tau=\\widehat{T}_z/2$ and the same tuning-frequency floor; they differ only in
the regularization law that period feeds, which is $r_S\\propto\\sigma_{aw}\\tau^3$
for the integral state of OU--III and the pair
$r_{p0}\\propto\\sigma_{aw}\\tau^2$, $r_{v0}\\propto\\sigma_{aw}\\tau$ for the
second-order chain of OU--II. Both sets of coefficients and clamps were
re-fitted for the wave-band period on the same four stationary records; neither
was fitted against the paired study reported below.""",
        """For that ablation to isolate the architecture, the two filters have to be
tuned from the same wave statistic. Both therefore derive the operating point
from the wave-band zero-crossing period of Sec.~\\ref{sec:wave-band-period}, with
$\\tau=\\widehat{T}_z/2$ and the same tuning-frequency floor. OU--III then uses
the integral-state schedule of Sec.~\\ref{sec:adaptation}, including its
pseudo-update cadence and information-rate renormalization; away from clamps its
filter-input standard deviation scales as
$r_{S,\\mathrm{filter}}\\propto\\sigma_{aw}\\tau^{5/2}$. OU--II uses the
second-order pair $r_{p0}\\propto\\sigma_{aw}\\tau^2$ and
$r_{v0}\\propto\\sigma_{aw}\\tau$. Both sets of coefficients and clamps were
re-fitted for the wave-band period on the same four stationary records; neither
was fitted against the paired study reported below.""",
    ),
    (
        """    \\caption{Paired vertical-displacement reconstruction and error for the
    medium JONSWAP case. The reference, adaptive PII observer, TVG--NLO,
    OU--II, and OU--III estimates are evaluated from the same motion and
    inertial record.}""",
        """    \\caption{Common-record vertical-displacement reconstruction and error for
    the medium JONSWAP case. The reference, adaptive PII observer, TVG--NLO,
    OU--II, and OU--III estimates are evaluated from the same motion and
    inertial record.}""",
    ),
    (
        """The paired result records, generated tables, and comparison plots are produced
by the automated deterministic evaluation.""",
        """The common-record outputs, generated tables, and comparison plots are produced
by the automated deterministic evaluation.""",
    ),
]

for old, new in replacements:
    if old not in text:
        raise SystemExit("expected wording anchor not found:\n" + old[:180])
    text = text.replace(old, new, 1)

path.write_text(text, encoding="utf-8")
