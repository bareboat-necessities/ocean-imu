#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "doc" / "kalman_ou_iii"


def read(name):
    return (DOC / name).read_text(encoding="utf-8")


def write(name, text):
    (DOC / name).write_text(text, encoding="utf-8")


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{label}: expected one occurrence, found {count}")
    return text.replace(old, new, 1)


def replace_from(text, marker, replacement, label):
    count = text.count(marker)
    if count != 1:
        raise SystemExit(f"{label}: expected one marker, found {count}")
    return text.split(marker, 1)[0] + replacement


def replace_between(text, start, end, replacement, label):
    if text.count(start) != 1 or text.count(end) != 1:
        raise SystemExit(f"{label}: section markers are not unique")
    before, rest = text.split(start, 1)
    _, after = rest.split(end, 1)
    return before + replacement + end + after


# 1) Core adaptation narrative: keep the former cubic law as a documented
# comparator, but make the reduced-MSE SpectralMSE schedule the deployed default.
name = "w3d-adaptation-motivation.tex-part"
text = read(name)
marker = r"\subsection{Applied cubic base and cadence normalization}"
replacement = r"""\subsection{Legacy cubic reference and deployed reduced-MSE schedule}
\label{sec:riccati-similarity-design-guidance}
\label{sec:riccati-deployed-cadence}

The fixed-corner analysis above historically led to an amplitude-coupled cubic
schedule.  It remains useful as the comparator used in the earlier coefficient
and cadence studies, but it is no longer the deployed OU--III default.  In that
former schedule the dimensional base was
\begin{equation}
\boxed{
r_{S,\star}^{\mathrm{legacy,base}}
=c_R\sigma_{aw,\star}\tau_\star^3.}
\label{eq:adapt-rs-base}
\end{equation}
Together with the self-similar pseudo-update period
\begin{equation}
T_S=\clip(c_T\tau_{\mathrm{appl}},T_{S,\min},T_{S,\max})
\label{eq:engineering-pseudo-cadence}
\end{equation}
and the reference-cadence information normalization
\begin{equation}
r_{S,\mathrm{filter}}^{\mathrm{legacy}}
=r_{S,\mathrm{appl}}^{\mathrm{legacy,base}}
\sqrt{\frac{T_{S,0}}{T_S}},
\label{eq:engineering-rs-filter}
\end{equation}
the unclipped former schedule has
\begin{equation}
\boxed{
r_{S,\mathrm{filter}}^{\mathrm{legacy}}
\propto\sigma_{aw}\tau^{5/2}.}
\label{eq:deployed-effective-law}
\end{equation}
The $5/2$ exponent remains a useful fixed-corner/cadence reference, but it is
not the displacement-MSE optimum.

For that historical comparator, combining the reduced corner relation with
Eqs.~\eqref{eq:adapt-rs-base} and \eqref{eq:engineering-rs-filter} gives
\begin{equation}
\boxed{
(\omega_R\tau)^6
=\frac{q_{\mathrm{eff}}}
 {c_R^2\sigma_{aw}^2T_{S,0}}.}
\label{eq:adapt-cR-spectral-meaning}
\end{equation}
This relation explains why the earlier amplitude-law experiments were
informative about fixed normalized pole placement; it does not define the
current target.

The current implementation instead uses the reduced physical-MSE optimum
derived in Sec.~\ref{sec:adaptation-reduced-mse-envelope} and fixed
analytically in Sec.~\ref{sec:adaptation-reduced-mse-validation}.  For the
realized target cadence
\begin{equation}
T_{S,\star}=\clip(c_T\tau_\star,T_{S,\min},T_{S,\max}),
\end{equation}
the default SpectralMSE target is
\begin{equation}
\boxed{
r_{S,\star}^{\mathrm{MSE}}
=C_J q_{\mathrm{eff}}^{1/14}
 \widehat\sigma_{a,B,\star}^{\,6/7}
 \tau_\star^{24/7}T_{S,\star}^{-1/2},
\qquad C_J=0.0538.}
\label{eq:adapt-default-rs-law}
\end{equation}
Here the acceleration amplitude is the physical band-limited wave statistic,
not the OU-prior standard deviation.  Under the strong-observation reduction
used by the deployed coefficient, $q_{\mathrm{eff}}\simeq2R_a h$ to leading
order.  The factor $T_S^{-1/2}$ is already part of
Eq.~\eqref{eq:adapt-default-rs-law}; therefore the default SpectralMSE target is
\emph{not} multiplied again by the reference-cadence factor of
Eq.~\eqref{eq:engineering-rs-filter}.  Away from cadence clamps,
$T_S\propto\tau$ and
\begin{equation}
\boxed{
r_{S,\star}^{\mathrm{MSE}}
\propto\widehat\sigma_{a,B}^{\,6/7}\tau^{41/14}.}
\end{equation}
Thus the deployed period exponent is $41/14\simeq2.929$, rather than $5/2$,
and the deployed acceleration exponent is $6/7\simeq0.857$, rather than one.

The resulting current schedule is summarized by
\begin{equation}
\boxed{
\begin{aligned}
\tau_\star&=c_\tau T_{\rm sea},\\
\sigma_{aw,\star}&=c_\sigma\widehat\sigma_{a,B},\\
r_{S,\star}&=C_J q_{\mathrm{eff}}^{1/14}
 \widehat\sigma_{a,B}^{\,6/7}\tau_\star^{24/7}T_{S,\star}^{-1/2}.
\end{aligned}}
\label{eq:adapt-three-layer-summary}
\end{equation}
The first two mappings set the OU prior; the third separately chooses the
integral-regularization strength by balancing drift rejection against distortion
of genuine wave displacement.  The cubic schedule remains a supported low-cost
comparison/fallback, but it is not the default adaptation law.
"""
text = replace_from(text, marker, replacement, name)
write(name, text)


# 2) Engineering implementation subsection: state the exact default target and
# avoid applying the cubic cadence normalization twice.
name = "w3d-fus-methods.tex-part"
text = read(name)
start = r"\subsection{Targets, smoothing, and application cadence}"
end = r"\subsection{Three-axis application}"
replacement = r"""\subsection{Targets, smoothing, and application cadence}
\label{sec:adaptation-target-implementation}

The OU targets and the target pseudo-update period are
\begin{equation}
\begin{aligned}
\tau_\star&=c_\tau\frac{\widehat T_z}{2},\\
\widehat\sigma_{a,B,\star}
 &=\sqrt{\max(\widehat\sigma_{\mathrm{wave}}^2,\epsilon)},\\
\sigma_{aw,\star}&=c_\sigma\widehat\sigma_{a,B,\star},\\
T_{S,\star}&=\clip(c_T\tau_\star,T_{S,\min},T_{S,\max}),\\
r_{S,\star}&=C_Jq_{\mathrm{eff}}^{1/14}
 \widehat\sigma_{a,B,\star}^{\,6/7}
 \tau_\star^{24/7}T_{S,\star}^{-1/2}.
\end{aligned}
\label{eq:tuner-targets}
\end{equation}
The last line is the deployed SpectralMSE law of
Eq.~\eqref{eq:adapt-default-rs-law}, with $C_J=0.0538$.  The implementation
uses the strong-observation reduced intensity $q_{\mathrm{eff}}\simeq2R_a h$
and caches its $1/14$ power.  The physical wave statistic
$\widehat\sigma_{a,B}$, rather than $\sigma_{aw}$, enters the distortion term;
in code the tuner supplies $\sigma_{aw}=c_\sigma\widehat\sigma_{a,B}$ and the
factor $c_\sigma$ is divided back out when evaluating the regularizer target.

Applied $\tau$ and $\sigma_{aw}$ use a common EMA,
\begin{equation}
\begin{aligned}
\alpha&=1-e^{-\Delta t/\tau_{\mathrm{adapt}}},\\
\tau_{\mathrm{appl}}&\leftarrow\tau_{\mathrm{appl}}
 +\alpha(\tau_\star-\tau_{\mathrm{appl}}),\\
\sigma_{aw,\mathrm{appl}}&\leftarrow\sigma_{aw,\mathrm{appl}}
 +\alpha(\sigma_{aw,\star}-\sigma_{aw,\mathrm{appl}}),
\end{aligned}
\label{eq:adapt-ewma}
\end{equation}
with the sea-scaled horizon
\begin{equation}
\tau_{\mathrm{adapt}}=c_{\mathrm{adapt}}T_{\mathrm{sea}},\qquad
c_{\mathrm{adapt}}=0.40.
\label{eq:adapt-ema-sea-time}
\end{equation}
The regularizer target is smoothed separately,
\begin{equation}
\begin{aligned}
\tau_{r_S}&=m_{r_S}\tau_\star,\qquad
\alpha_{r_S}=1-e^{-\Delta t/\tau_{r_S}},\\
r_{S,\mathrm{appl}}&\leftarrow
 r_{S,\mathrm{appl}}
 +\alpha_{r_S}(r_{S,\star}-r_{S,\mathrm{appl}}),
\end{aligned}
\label{eq:adapt-RS}
\end{equation}
where the measured default is $m_{r_S}=1.5$.  Every valid physical sample
contributes to the tuner and to all three EMAs.  The smoothed candidate is
activated at approximately \SI{0.1}{s} intervals and staged by one IMU sample;
this activation cadence does not decimate the measurement stream.  Parameters
used at step $k$ are therefore measurable with respect to data through $k-1$.

The actual pseudo-update cadence follows the applied OU time scale,
\begin{equation}
T_S=\clip(c_T\tau_{\mathrm{appl}},T_{S,\min},T_{S,\max}),
\qquad
c_T=\frac{\SI{15}{ms}}{\SI{1.1}{s}},
\label{eq:engineering-pseudo-cadence-implementation}
\end{equation}
with $T_{S,\min}=\SI{5}{ms}$ and $T_{S,\max}=\SI{250}{ms}$.  For the deployed
SpectralMSE law, the $T_S^{-1/2}$ dependence is already contained in the target,
so the filter receives
\begin{equation}
r_{S,\mathrm{filter}}=r_{S,\mathrm{appl}}
\qquad\text{(SpectralMSE default).}
\label{eq:engineering-rs-filter-implementation}
\end{equation}
No additional $\sqrt{T_{S,0}/T_S}$ factor is applied.  That reference-cadence
renormalization belongs only to the supported cubic fallback described by
Eq.~\eqref{eq:engineering-rs-filter}.  Away from cadence clamps, the deployed
schedule therefore has
\begin{equation}
r_{S,\mathrm{filter}}
\propto\widehat\sigma_{a,B}^{\,6/7}\tau^{41/14},
\end{equation}
which is the reduced physical-MSE prediction rather than the former
$\sigma_{aw}\tau^{5/2}$ schedule.

"""
text = replace_between(text, start, end, replacement, name)
# Update the implementation-default table that drifted with the old schedule.
repls = [
    (r"adaptation commit / EMA horizon & $0.1/1.8$ s",
     r"adaptation activation / common EMA horizon & $0.1$ s / $0.40T_{\mathrm{sea}}$"),
    (r"& coefficients $(c_\tau,c_\sigma,c_R)$ & $(1.0,0.9,0.35)$",
     r"& default coefficients $(c_\tau,c_\sigma,C_J)$ & $(1.0,0.9,0.0538)$"),
    (r"$r_S$ smoothing horizon & $3\tau_\star$",
     r"$r_S$ smoothing horizon & $1.5\tau_\star$"),
    (r"& frequency / variance horizons & \SI{1.0}{s} / $2T_{\mathrm{eff}}$",
     r"& frequency / variance horizons & $0.10T_{\mathrm{sea}}$ / $2T_{\mathrm{eff}}$"),
]
for old, new in repls:
    text = replace_once(text, old, new, f"{name}: {old[:35]}")
write(name, text)


# 3) Cadence-only discussion: clearly mark it as the former cubic comparator.
name = "w3d-adaptation-cadence-interpretation.tex-part"
text = read(name)
old = r"""For the applied OU--III schedule, $a_R=\sigma_{aw}$, so this argument gives
$r_{S,\mathrm{filter}}\propto\sigma_{aw}\tau^{5/2}$.  The amplitude-free
experiment replaced $a_R$ by the constant sensor-noise scale $\sqrt{R_a}$ and
confirmed the same $5/2$ period behavior but degraded the multi-seed endpoint.
The subsequent sensor-floor-plus-leakage interpolation returned to
$a_R\propto\sigma_{aw}$ over the entire tested envelope and tied the applied law
under paired ten-seed evaluation.  The cadence argument is therefore supported
as an explanation of the period exponent while remaining deliberately agnostic
about the correct steady-state amplitude scale.  The remaining observed
stationary--transition tradeoff concerns how that amplitude target is followed
in time, not the square-root cadence normalization itself."""
new = r"""For the former amplitude-coupled cubic OU--III schedule,
$a_R=\sigma_{aw}$, so this argument gives
$r_{S,\mathrm{filter}}\propto\sigma_{aw}\tau^{5/2}$.  The amplitude-free
experiment replaced $a_R$ by the constant sensor-noise scale $\sqrt{R_a}$ and
confirmed the same $5/2$ period behavior but degraded the multi-seed endpoint.
The subsequent sensor-floor-plus-leakage interpolation returned to
$a_R\propto\sigma_{aw}$ over the tested envelope and tied that legacy cubic
comparator under paired ten-seed evaluation.  These results remain useful for
interpreting the fixed-corner/cadence family, but they do not define the current
default.  SpectralMSE instead contains $T_S^{-1/2}$ directly in its reduced-MSE
target and is not multiplied by this reference-cadence normalization a second
time."""
text = replace_once(text, old, new, name)
write(name, text)


# 4) Reduced-MSE envelope: relabel the former applied law as the legacy cubic
# comparator. Numerical values and equation labels are intentionally unchanged.
name = "w3d-reduced-mse-envelope.tex-part"
text = read(name)
for old, new in [
    ("The result is numerically close to the applied schedule", "The result is numerically close to the legacy amplitude-coupled cubic schedule"),
    (r"r_{S,\rm appl}\propto\sigma_{aw}\tau^{5/2}", r"r_{S,\rm cubic}\propto\sigma_{aw}\tau^{5/2}"),
    (r"\frac{r_{S,\rm MSE}^\star}{r_{S,\rm appl}}", r"\frac{r_{S,\rm MSE}^\star}{r_{S,\rm cubic}}"),
    ("Applied and reduced physical-MSE schedule shapes", "Legacy cubic and reduced physical-MSE schedule shapes"),
    (r"applied $\sigma_{aw}\tau^{5/2}$", r"legacy cubic $\sigma_{aw}\tau^{5/2}$"),
    ("of the applied $r_S$ schedule", "of the legacy cubic $r_S$ schedule"),
    (r"r_{S,\rm appl}\propto\tau^{3.607}", r"r_{S,\rm cubic}\propto\tau^{3.607}"),
    ("the applied and reduced-MSE laws", "the legacy cubic and reduced-MSE laws"),
    ("the applied $(1,5/2)$ pair", "the legacy cubic $(1,5/2)$ pair"),
    ("why the empirically selected law lies near the observed optimum", "why the former empirical cubic law lay near the observed optimum"),
]:
    if old in text:
        text = text.replace(old, new)
# Add a closing sentence making the now-completed deployment decision explicit.
text += r"""

The paired evaluation reported in
Sec.~\ref{sec:adaptation-reduced-mse-validation} subsequently favored the
reduced-MSE schedule and the implementation promoted SpectralMSE to the default.
Accordingly, the cubic curve in this subsection is a historical comparator, not
the current deployed adaptation law.
"""
write(name, text)


# 5) Coefficient-law investigation predates promotion of SpectralMSE. Preserve
# the experiment, but remove current-tense claims that the old cubic law is
# deployed.
name = "w3d-adaptation-coefficient-investigation.tex-part"
text = read(name)
anchor = r"\label{sec:adaptation-coefficient-investigation}"
intro = anchor + "\n\n" + r"""This subsection records coefficient experiments performed against the former
amplitude-coupled cubic schedule.  That schedule is retained here as the
\emph{legacy cubic comparator}; the later reduced-MSE study promoted SpectralMSE
to the deployed default.  The numerical comparisons below are therefore
historical evidence about the route to the present law, not a statement that
$\sigma_{aw}\tau^3$ is currently applied."""
text = replace_once(text, anchor, intro, name + " intro")
for old, new in [
    ("replacing the applied\namplitude-dependent law", "replacing the legacy cubic\namplitude-dependent law"),
    ("compared with the applied schedule", "compared with the legacy cubic schedule"),
    ("The applied law supplies", "The legacy cubic law supplies"),
    ("Applied & 4.181", "Legacy cubic & 4.181"),
    ("compared with the applied coefficient $0.35$", "compared with the legacy cubic coefficient $0.35$"),
    ("the applied linear-$\\sigma_{aw}$ law", "the legacy linear-$\\sigma_{aw}$ cubic law"),
    ("the applied and fitted laws give", "the legacy cubic and fitted laws give"),
    ("candidate minus applied normalized", "candidate minus legacy-cubic normalized"),
    ("tied the applied law", "tied the legacy cubic comparator"),
    ("improve the applied steady-state amplitude law", "improve the legacy cubic steady-state amplitude law"),
]:
    text = text.replace(old, new)
write(name, text)


# 6) Reduced-MSE validation: make clear that the winning law is now default and
# the other arm is the historical cubic comparator.
name = "w3d-reduced-mse-validation.tex-part"
text = read(name)
for old, new in [
    ("a direct paired comparison against the applied law,", "a direct paired comparison against the legacy amplitude-coupled cubic law,"),
    ("against the applied schedule. Vertical-displacement", "against the legacy cubic comparator. Vertical-displacement"),
    ("Scenario & Applied & MSE", "Scenario & Cubic & MSE"),
    ("against none for the applied cubic\nlaw", "against none for the cubic\ncomparator"),
    ("The applied cubic\nschedule is therefore retained as the low-cost configuration for embedded\ndeployment", "The cubic\nschedule is therefore retained as an explicitly selectable low-cost configuration"),
    ("and Eq.~\\eqref{eq:adapt-mse-applied-form} is the default elsewhere.", "while Eq.~\\eqref{eq:adapt-mse-applied-form} is the deployed/default law."),
]:
    text = text.replace(old, new)
# Strengthen deployment subsection opening without altering measured numbers.
old = "Both schedules are retained as supported configurations."
new = "SpectralMSE is the deployed/default configuration; the cubic family is retained as a supported low-cost comparator/fallback."
text = replace_once(text, old, new, name + " deployment")
write(name, text)


# 7) Historical robustness subsection: its coupled sweep really used the old
# 0.35*sigma_aw*tau^3 relation, so label it as such instead of calling it current.
name = "w3d-ou-robustness.tex-part"
text = read(name)
old = r"""A ten-seed robustness study perturbs the nominal OU--III operating point by
multipliers $0.5$, $0.75$, $1.0$, $1.25$, and $1.5$.  One-factor sweeps vary
$\tau$, $\sigma_{aw}$, or $r_S$ with the other parameters fixed; coupled sweeps
also move $r_S$ according to the applied relation
$r_S=\clip(0.35\,\sigma_{aw}\tau^3,0.15,400)$.  Before clipping, a coupled scale
change therefore gives $r_S\to c^{3}r_S$."""
new = r"""A ten-seed robustness study, performed before SpectralMSE became the default,
perturbs the then-nominal OU--III operating point by multipliers $0.5$, $0.75$,
$1.0$, $1.25$, and $1.5$.  One-factor sweeps vary $\tau$, $\sigma_{aw}$, or
$r_S$ with the other parameters fixed; the historical coupled sweeps also move
$r_S$ according to the former amplitude-coupled cubic relation
$r_S=\clip(0.35\,\sigma_{aw}\tau^3,0.15,400)$.  Before clipping, a coupled scale
change therefore gives $r_S\to c^{3}r_S$.  These coupled rows characterize the
legacy cubic schedule and must not be read as the current SpectralMSE law."""
text = replace_once(text, old, new, name)
text = text.replace("the sensitivity of the applied tuning law", "the sensitivity of that legacy cubic tuning law")
write(name, text)


# 8) Conclusion already has the correct law; make the default/fallback wording
# match the literal code default rather than imply an automatic target switch.
name = "w3d-conclusion-summary.tex-part"
text = read(name)
old = r"""SpectralMSE is therefore the
default adaptation law outside the low-cost embedded configuration, while the
cubic schedule remains a supported implementation option when transcendental
cost dominates a sub-percent displacement difference."""
new = r"""SpectralMSE is therefore the deployed/default adaptation law, while the cubic
schedule remains an explicitly selectable low-cost implementation option when
transcendental cost dominates a sub-percent displacement difference."""
text = replace_once(text, old, new, name)
write(name, text)

print("OU-III manuscript synchronized to SpectralMSE default")
