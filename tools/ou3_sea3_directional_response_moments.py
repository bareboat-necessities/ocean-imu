#!/usr/bin/env python3
"""SEA0 directional-response and response-weighted moment subcertificate.

This module discharges the *response* half of the SEA0 obligation declared by
the SEA3 directional-sea theorem: it turns the three-partition directional
elevation spectrum into the matrix-valued source spectrum

    Phi_u(omega) = int G_imu(omega,theta) E(omega,theta) G_imu(omega,theta)^* dtheta

and into the specific response-weighted moments that the *deployed* source path
actually consumes.  It does not promote SEA0, P2, P3, P4 or P5, and it does not
read any replay trajectory.

Three things are established analytically and are independent of the numerical
screen below:

1. **Modulus majorant.**  A single vessel is never assumed.  The response is an
   admissible *family* bounded componentwise by a declared magnitude envelope,
   and the Gram matrix of that envelope dominates every member of the family in
   the quadratic-form sense.  Directional and cross-axis coupling therefore
   survive as off-diagonal entries; the certificate never reduces to
   independent per-axis scalars.

2. **Exact leak inversion.**  In continuous-time steady state the deployed
   `WavePeriodEstimator` moment ratio is *exactly* the response-weighted mean
   square frequency of a fourth-order high-passed elevation spectrum.  The
   estimator's `omega^2 = (sigma_v/sigma_eta)^2 - lambda^2` line is not a
   narrow-band approximation; it is an identity for any input spectrum.

3. **Mixture convexity.**  That mean square is a convex combination over the
   three partitions, so the extreme deployed periods of the whole SEA3 class are
   attained on one-partition seas.  The 15-parameter period-channel enclosure
   collapses to a two-parameter (Tp, gamma) screen.

Two obligations come out of the construction and both are reported rather than
tuned away.

The first is the *sigma clamp rail*.  Shipping code sets
`sigma_target_ = min(0.9 * sigma_a, 4 m/s^2)`, so the sigma_aw source
coordinate is contained in its declared interval by saturation and cannot be
violated by any sea.  What a sea can do is reach the rail, and an admissible
crossing sea of a 4.6 s wind system on a 10.9 s swell does: its band-passed
deviation of 4.66 m/s^2 passes the 4.44 m/s^2 saturation point.  The rail is
therefore a live cell of the source language rather than a dormant branch, and
on it the sea-to-source map stops being injective, so SEA3 cannot prune the P2
language there.

Which band reference reaches it is the sharper fact.  Neither frequency-source
mode reaches the rail at rest: 3.91 m/s^2 at each sea's own deployed period and
4.06 m/s^2 at the fixed 0.2 Hz startup prior are both below saturation.  The
rail is reached only at intermediate references, that is on the estimator-lag
branch *between* the two source modes the startup subcertificate established.

The second is *aliasing*.  The deployed chain is sampled at 5 ms, so any claim
about what it observes must bound the folded tail of the sea's specific-force
spectrum.  With no vessel/IMU roll-off that folded power is logarithmically
divergent, and not divergent by a negligible amount: three unmodelled decades
already fold in more than the whole clamp interval.  A certified high-frequency
roll-off is therefore a mandatory SEA0 assumption rather than a modelling
convenience, and sampling supplies no band limit of its own.

The numerical parts are fixed-grid quadrature with analytical tail pads, not
validated interval integration, so every numerical interval reported here is a
screening enclosure and is explicitly marked non-promoting.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

import ou3_sea3_spectral_moment_bridge as bridge  # noqa: E402


SCHEMA_VERSION = "OU3_SEA3_DIRECTIONAL_RESPONSE_MOMENTS_V1"

G_STD = 9.80665

# ---------------------------------------------------------------------------
# Declared three-partition directional sea domain (SEA3-D).
#
# The peak-period window covers the eight reference seas (Tz ~= 2.3..8.4 s,
# hence Tp ~= 2.8..11.9 s through the certified Tz/Tp screen) and extends to
# long swell.  The height ceiling and the significant-steepness ceiling are
# physical admissibility statements about the sea, not proof-convenience
# clamps: a wind sea cannot carry arbitrary height at a fixed period without
# breaking.  0.10 is deliberately above the conventional fully developed value
# (~0.05) so the declared class stays wider than the calibration seas.
# ---------------------------------------------------------------------------
M_MAX = bridge.M_MAX
TP_MIN_S = 2.5
TP_MAX_S = 20.0
GAMMA_MIN = bridge.GAMMA_MIN
GAMMA_MAX = bridge.GAMMA_MAX
# Matches significant_wave_height_Hs_upper_m in the declared operating
# domain; the sea domain must not be narrower than the deployment envelope.
HS_MAX_M = 8.5
SIGNIFICANT_STEEPNESS_MAX = 0.10
SPREAD_S_MIN = 1.0
SPREAD_S_MAX = 25.0

# ---------------------------------------------------------------------------
# Declared conservative vessel/IMU response envelope.
#
# The lever arm is disabled in the current proof scope, so angular motion does
# not produce specific force at the IMU and only translational response enters.
# RHO_MAX bounds the ratio of the levelled vessel/IMU acceleration amplitude to
# the undisturbed surface orbital acceleration amplitude of the same wave.  A
# contouring displacement hull peaks near 1.1-1.3; 2.0 leaves room for
# resonance and for the private-Mahony levelling residual.
#
# ROLL_OFF_EXPONENT/OMEGA_L declare the finite-waterline (Froude-Krylov length
# averaging) roll-off.  A hull of waterline L cannot follow waves much shorter
# than L.  L_REF is deliberately *small*: a shorter reference hull pushes the
# roll-off corner higher and therefore assumes less.
# ---------------------------------------------------------------------------
RHO_MAX = 2.0
ROLL_OFF_EXPONENT = 1.0
L_REF_M = 4.0
K_L = 2.0 * math.pi / L_REF_M
OMEGA_L = math.sqrt(G_STD * K_L)

# ---------------------------------------------------------------------------
# Frozen shipping source parity (see docs/ou3-proof-research-state.md).
# ---------------------------------------------------------------------------
IMU_DT_S = 0.005
WAVE_PERIOD_HIGH_PASS_HZ = 0.02
LAMBDA_LEAK = 2.0 * math.pi * WAVE_PERIOD_HIGH_PASS_HZ
MIN_TUNE_FREQ_HZ = 0.03
MAX_TUNE_FREQ_HZ = 1.2
SIGMA_BAND_LOW_RATIO = 0.5
SIGMA_BAND_HIGH_RATIO = 4.0
SIGMA_BAND_MIN_HZ = 0.01
SIGMA_BAND_MAX_HZ = 6.0
# sigma_target_ = min(sigma_coeff_ * sigma_a, max_sigma_a_) in shipping code.
# The ceiling is a saturating clamp, so sigma_aw cannot leave [0, 4]; what the
# sea can do is drive the tuner onto the rail.
SIGMA_COEFF = 0.9
SIGMA_AW_CLAMP_MS2 = 4.0
SIGMA_A_SATURATION_MS2 = SIGMA_AW_CLAMP_MS2 / SIGMA_COEFF
TUNE_FREQ_PRIOR_HZ = 0.2

# ---------------------------------------------------------------------------
# Quadrature configuration.
# ---------------------------------------------------------------------------
OMEGA_MIN = 1.0e-3
OMEGA_MAX = 1.0e5
OMEGA_PANELS = 3072
THETA_PANELS = 512
TAIL_PAD_SAFETY = 1.000001

TP_SCAN = 25
GAMMA_SCAN = 7
F_REF_SCAN = 13
MIXTURE_CELL_STRIDE = 5
MIXTURE_WEIGHT_STEPS = 21


# ---------------------------------------------------------------------------
# Quadrature primitives.
# ---------------------------------------------------------------------------
def _log_omega_grid() -> tuple[list[float], list[float]]:
    """Simpson nodes and dw weights on a log-spaced angular-frequency grid."""
    if OMEGA_PANELS % 2:
        raise RuntimeError("OMEGA_PANELS must be even")
    log_lo = math.log(OMEGA_MIN)
    log_hi = math.log(OMEGA_MAX)
    step = (log_hi - log_lo) / OMEGA_PANELS
    nodes: list[float] = []
    weights: list[float] = []
    for i in range(OMEGA_PANELS + 1):
        omega = math.exp(log_lo + i * step)
        if i == 0 or i == OMEGA_PANELS:
            coef = 1.0
        elif i % 2:
            coef = 4.0
        else:
            coef = 2.0
        nodes.append(omega)
        # d(omega) = omega * d(log omega).
        weights.append(coef * step / 3.0 * omega)
    return nodes, weights


OMEGA_NODES, OMEGA_WEIGHTS = _log_omega_grid()


def _integrate(values: list[float]) -> float:
    total = 0.0
    for weight, value in zip(OMEGA_WEIGHTS, values):
        total += weight * value
    return total


# ---------------------------------------------------------------------------
# Sea family.
# ---------------------------------------------------------------------------
def _shape_normalizer(gamma: float) -> float:
    """int_0^inf shape(x;gamma) dx, so that m0 of a partition equals H^2/16."""
    def shape(x: float) -> float:
        return bridge.jonswap_dimensionless_shape(x, gamma)

    # The shape is scale free, so the same log grid is reused in x.
    values = [shape(omega) for omega in OMEGA_NODES]
    return _integrate(values)


def partition_elevation_spectrum(
    height_m: float,
    peak_period_s: float,
    gamma: float,
) -> list[float]:
    """S_eta(omega) on the shared grid, normalized so m0 = H^2/16."""
    omega_p = 2.0 * math.pi / peak_period_s
    scale = (height_m * height_m / 16.0) / (omega_p * _shape_normalizer(gamma))
    return [
        scale * bridge.jonswap_dimensionless_shape(omega / omega_p, gamma)
        for omega in OMEGA_NODES
    ]


def directional_spreading(theta: float, beta: float, spread_s: float) -> float:
    """Normalized cos^{2s} spreading, int over [-pi,pi) equal to one."""
    norm = math.gamma(spread_s + 1.0) / (
        2.0 * math.sqrt(math.pi) * math.gamma(spread_s + 0.5)
    )
    return norm * math.cos(0.5 * (theta - beta)) ** (2.0 * spread_s)


def directional_gram(beta: float, spread_s: float) -> list[list[float]]:
    """Q(beta,s) = int D(theta) u(theta) u(theta)^T dtheta, u = [1,|cos|,|sin|]."""
    if THETA_PANELS % 2:
        raise RuntimeError("THETA_PANELS must be even")
    step = 2.0 * math.pi / THETA_PANELS
    gram = [[0.0] * 3 for _ in range(3)]
    for i in range(THETA_PANELS + 1):
        theta = -math.pi + i * step
        if i == 0 or i == THETA_PANELS:
            coef = 1.0
        elif i % 2:
            coef = 4.0
        else:
            coef = 2.0
        weight = coef * step / 3.0 * directional_spreading(theta, beta, spread_s)
        unit = (1.0, abs(math.cos(theta)), abs(math.sin(theta)))
        for a in range(3):
            for b in range(3):
                gram[a][b] += weight * unit[a] * unit[b]
    return gram


# ---------------------------------------------------------------------------
# Response envelope and deployed channel weights.
# ---------------------------------------------------------------------------
def response_envelope(omega: float) -> float:
    """rho(omega): declared bound on |vessel/IMU accel| / |orbital accel|."""
    if omega <= OMEGA_L:
        return RHO_MAX
    return RHO_MAX * (OMEGA_L / omega) ** ROLL_OFF_EXPONENT


RESPONSE_ENVELOPE = [response_envelope(omega) for omega in OMEGA_NODES]


def up_specific_force_spectrum(elevation: list[float]) -> list[float]:
    """S_up(omega) majorant: rho(omega)^2 omega^4 S_eta(omega)."""
    return [
        rho * rho * omega**4 * value
        for rho, omega, value in zip(RESPONSE_ENVELOPE, OMEGA_NODES, elevation)
    ]


def proxy_elevation_weight(omega: float) -> float:
    """|G_eta(i omega)|^2 from acceleration input through the deployed chain.

    Two shared high-pass stages s/(s+lambda) followed by two leaky integrations
    1/(s+lambda) give omega^4 / (omega^2 + lambda^2)^4.
    """
    return omega**4 / (omega * omega + LAMBDA_LEAK * LAMBDA_LEAK) ** 4


PROXY_ELEVATION_WEIGHT = [proxy_elevation_weight(omega) for omega in OMEGA_NODES]


def sigma_band_corners_hz(f_ref_hz: float) -> tuple[float, float]:
    """Exact AdaptiveWaveBandPass corner selection at the shipping settings."""
    nyquist_guard_hz = 0.45 / IMU_DT_S
    upper_limit = min(SIGMA_BAND_MAX_HZ, nyquist_guard_hz)
    low_hz = max(SIGMA_BAND_MIN_HZ, SIGMA_BAND_LOW_RATIO * f_ref_hz)
    low_hz = min(low_hz, upper_limit / 1.05)
    high_hz = min(upper_limit, SIGMA_BAND_HIGH_RATIO * f_ref_hz)
    high_hz = max(high_hz, low_hz * 1.05)
    high_hz = min(high_hz, upper_limit)
    return low_hz, high_hz


def sigma_band_weight(f_ref_hz: float) -> list[float]:
    """|B(i omega)|^2 for the one-pole high-pass / one-pole low-pass cascade."""
    low_hz, high_hz = sigma_band_corners_hz(f_ref_hz)
    omega_l = 2.0 * math.pi * low_hz
    omega_h = 2.0 * math.pi * high_hz
    weights = []
    for omega in OMEGA_NODES:
        sq = omega * omega
        weights.append(
            (sq / (sq + omega_l * omega_l))
            * (omega_h * omega_h / (sq + omega_h * omega_h))
        )
    return weights


# ---------------------------------------------------------------------------
# Analytical tails.
# ---------------------------------------------------------------------------
def _asymptotic_acceleration_constant(height_m: float, peak_period_s: float,
                                      gamma: float) -> float:
    """A with omega^4 S_eta(omega) -> A / omega as omega -> inf."""
    omega_p = 2.0 * math.pi / peak_period_s
    return (height_m * height_m / 16.0) * omega_p**4 / _shape_normalizer(gamma)


def acceleration_tail_above(omega_cut: float, height_m: float,
                            peak_period_s: float, gamma: float) -> float:
    """Upper bound for int_{omega_cut}^inf S_up(omega) d omega.

    Beyond twice the peak the JONSWAP peak factor is 1 to within 1e-26, so the
    ideal omega^{-5} elevation tail dominates and the response roll-off supplies
    the remaining omega^{-2p}.  The bound diverges for ROLL_OFF_EXPONENT = 0,
    which is exactly the aliasing obligation reported below.
    """
    if omega_cut < 2.0 * (2.0 * math.pi / peak_period_s):
        raise ValueError("tail bound requires omega_cut above twice the peak")
    const = _asymptotic_acceleration_constant(height_m, peak_period_s, gamma)
    if ROLL_OFF_EXPONENT <= 0.0:
        return math.inf
    exponent = 2.0 * ROLL_OFF_EXPONENT
    tail = (
        RHO_MAX * RHO_MAX
        * OMEGA_L**exponent
        * const
        / (exponent * omega_cut**exponent)
    )
    return TAIL_PAD_SAFETY * tail


def unbounded_response_alias_power(decades: float, tail_constant: float) -> float:
    """Folded specific-force power over `decades` of tail if rho were flat.

    A flat response leaves the specific-force spectral density proportional to
    1/omega, so every decade above Nyquist folds in the same power and the total
    grows without bound.  The result is independent of where the tail starts,
    which is exactly why sampling cannot supply the missing band limit.
    """
    return RHO_MAX * RHO_MAX * tail_constant * decades * math.log(10.0)


# ---------------------------------------------------------------------------
# Deployed-channel statistics for one partition.
# ---------------------------------------------------------------------------
def _partition_unit_spectrum(peak_period_s: float, gamma: float) -> list[float]:
    """S_eta for H = 4 m, i.e. unit m0.  Ratios below are H independent."""
    return partition_elevation_spectrum(4.0, peak_period_s, gamma)


def surface_zero_crossing_period(elevation: list[float]) -> float:
    m0 = _integrate(elevation)
    m2 = _integrate(
        [omega * omega * value for omega, value in zip(OMEGA_NODES, elevation)]
    )
    return 2.0 * math.pi * math.sqrt(m0 / m2)


def deployed_period_s(elevation: list[float]) -> float:
    """T_z as the deployed estimator forms it, in continuous-time steady state.

    By the exact leak-inversion identity the estimator's omega^2 is the mean
    square frequency of the measure d mu = |G_eta|^2 S_up d omega, so no
    separate velocity-proxy integral is needed.
    """
    source = up_specific_force_spectrum(elevation)
    measure = [w * s for w, s in zip(PROXY_ELEVATION_WEIGHT, source)]
    mass = _integrate(measure)
    second = _integrate(
        [omega * omega * value for omega, value in zip(OMEGA_NODES, measure)]
    )
    if not (mass > 0.0 and second > 0.0):
        raise RuntimeError("deployed period measure lost positivity")
    return 2.0 * math.pi / math.sqrt(second / mass)


def deployed_velocity_proxy_variance(elevation: list[float]) -> float:
    source = up_specific_force_spectrum(elevation)
    return _integrate(
        [
            (omega * omega + LAMBDA_LEAK * LAMBDA_LEAK) * w * s
            for omega, w, s in zip(OMEGA_NODES, PROXY_ELEVATION_WEIGHT, source)
        ]
    )


def steepness_height_cap_m(surface_tz_s: float) -> float:
    """Largest admissible Hs at this surface Tz under the steepness ceiling."""
    cap = SIGNIFICANT_STEEPNESS_MAX * G_STD * surface_tz_s**2 / (2.0 * math.pi)
    return min(HS_MAX_M, cap)


# ---------------------------------------------------------------------------
# Scans.
# ---------------------------------------------------------------------------
def _log_span(lo: float, hi: float, count: int) -> list[float]:
    if count < 2:
        raise ValueError("count must be at least two")
    step = (math.log(hi) - math.log(lo)) / (count - 1)
    return [math.exp(math.log(lo) + i * step) for i in range(count)]


def _linear_span(lo: float, hi: float, count: int) -> list[float]:
    if count < 2:
        raise ValueError("count must be at least two")
    step = (hi - lo) / (count - 1)
    return [lo + i * step for i in range(count)]


def period_channel_screen() -> dict[str, Any]:
    """Single-partition (Tp, gamma) scan of the deployed period channel."""
    peak_periods = _log_span(TP_MIN_S, TP_MAX_S, TP_SCAN)
    gammas = _linear_span(GAMMA_MIN, GAMMA_MAX, GAMMA_SCAN)

    ratio_lo = math.inf
    ratio_hi = -math.inf
    tz_lo = math.inf
    tz_hi = -math.inf
    limiter_lo: tuple[float, float] = (0.0, 0.0)
    limiter_hi: tuple[float, float] = (0.0, 0.0)
    surface_bias_lo = math.inf
    surface_bias_hi = -math.inf
    for peak_period in peak_periods:
        for gamma in gammas:
            elevation = _partition_unit_spectrum(peak_period, gamma)
            deployed = deployed_period_s(elevation)
            surface = surface_zero_crossing_period(elevation)
            ratio = deployed / peak_period
            if ratio < ratio_lo:
                ratio_lo = ratio
                limiter_lo = (peak_period, gamma)
            if ratio > ratio_hi:
                ratio_hi = ratio
                limiter_hi = (peak_period, gamma)
            tz_lo = min(tz_lo, deployed)
            tz_hi = max(tz_hi, deployed)
            bias = deployed / surface
            surface_bias_lo = min(surface_bias_lo, bias)
            surface_bias_hi = max(surface_bias_hi, bias)

    return {
        "scan": {"Tp_values": TP_SCAN, "gamma_values": GAMMA_SCAN},
        "deployed_Tz_over_Tp": [ratio_lo, ratio_hi],
        "deployed_Tz_s": [tz_lo, tz_hi],
        "deployed_over_surface_Tz": [surface_bias_lo, surface_bias_hi],
        "limiting_cell_for_lower_ratio": list(limiter_lo),
        "limiting_cell_for_upper_ratio": list(limiter_hi),
        "induced_tuner_frequency_hz": [1.0 / tz_hi, 1.0 / tz_lo],
        "inside_committed_tuning_channel": (
            1.0 / tz_hi >= MIN_TUNE_FREQ_HZ and 1.0 / tz_lo <= MAX_TUNE_FREQ_HZ
        ),
        "extremes_valid_for_three_partition_class": True,
        "reason": "mixture convexity of the mu-weighted mean square frequency",
        "promotion_use": False,
    }


def mixture_convexity_check() -> dict[str, Any]:
    """Verify numerically that mixtures stay inside single-partition extremes."""
    cases = [
        [(3.0, 4.0, 1.0), (2.0, 9.0, 3.3), (0.0, 12.0, 1.0)],
        [(1.0, 3.0, 7.0), (2.5, 14.0, 1.0), (1.5, 6.5, 2.0)],
        [(4.0, 2.5, 1.0), (0.5, 20.0, 1.0), (0.0, 5.0, 1.0)],
    ]
    rows = []
    worst_excess = 0.0
    for case in cases:
        active = [(h, tp, g) for h, tp, g in case if h > 0.0]
        totals = [0.0] * len(OMEGA_NODES)
        singles = []
        for height, peak_period, gamma in active:
            spectrum = partition_elevation_spectrum(height, peak_period, gamma)
            totals = [a + b for a, b in zip(totals, spectrum)]
            singles.append(deployed_period_s(spectrum))
        mixed = deployed_period_s(totals)
        lo = min(singles)
        hi = max(singles)
        excess = max(lo - mixed, mixed - hi, 0.0) / hi
        worst_excess = max(worst_excess, excess)
        rows.append(
            {
                "partitions": [list(entry) for entry in case],
                "single_partition_deployed_Tz_s": singles,
                "mixture_deployed_Tz_s": mixed,
                "inside_component_extremes": lo - 1e-9 <= mixed <= hi + 1e-9,
            }
        )
    return {
        "cases": rows,
        "worst_relative_excess": worst_excess,
        "identity": "omega_est^2 = sum_r c_r E_{mu_r}[omega^2] / sum_r c_r",
    }


def sigma_channel_screen() -> dict[str, Any]:
    """Band-passed specific-force standard deviation over the declared domain.

    Physical admissibility is imposed at two levels and both are needed.  Each
    partition carries its own breaking limit, because a 2.5 s wind sea cannot
    reach swell heights whatever the rest of the sea does, and the combined sea
    carries the mixture limit through its own zero-crossing period.  Dropping
    the per-partition limit admits a 0.28-steepness chop riding on a long swell,
    which is not a sea.

    Two maxima are reported and they answer different questions.  The *sound*
    bound maximizes over every band-reference frequency the shipping code can
    present, because the deployed reference is not always the sea's own period:
    the startup branch runs on the fixed 0.2 Hz prior and the estimator lags the
    sea afterwards.  The *self-consistent* diagnostic instead evaluates each sea
    at its own deployed period and is the number a settled run sees.

    The one-partition maximum is exhaustive on its declared grid.  The
    multimodal number is a two-partition screen on a coarser grid: under the
    enriched admissibility set the earlier simplex-vertex reduction no longer
    applies, so it is reported as a witness search rather than as an optimum.
    """
    peak_periods = _log_span(TP_MIN_S, TP_MAX_S, TP_SCAN)
    gammas = _linear_span(GAMMA_MIN, GAMMA_MAX, GAMMA_SCAN)
    references = _log_span(MIN_TUNE_FREQ_HZ, MAX_TUNE_FREQ_HZ, F_REF_SCAN)
    band_weights = {f_ref: sigma_band_weight(f_ref) for f_ref in references}

    cells = []
    for peak_period in peak_periods:
        for gamma in gammas:
            unit = _partition_unit_spectrum(peak_period, gamma)
            source = up_specific_force_spectrum(unit)
            surface_tz = surface_zero_crossing_period(unit)
            cells.append(
                {
                    "Tp_s": peak_period,
                    "gamma": gamma,
                    "surface_Tz_s": surface_tz,
                    "deployed_Tz_s": deployed_period_s(unit),
                    "height_cap_m": steepness_height_cap_m(surface_tz),
                    "inverse_square_period": 1.0 / (surface_tz * surface_tz),
                    "gains": {
                        f_ref: _integrate(
                            [w * s for w, s in zip(band_weights[f_ref], source)]
                        )
                        for f_ref in references
                    },
                }
            )

    def _sigma(gain: float, height_m: float) -> float:
        # The unit spectrum carries m0 = 1 m^2, i.e. H = 4 m.
        return math.sqrt(gain) * height_m / 4.0

    single_best = 0.0
    single_cell: dict[str, Any] = {}
    self_consistent_best = 0.0
    self_consistent_cell: dict[str, Any] = {}
    for cell in cells:
        for f_ref in references:
            sigma = _sigma(cell["gains"][f_ref], cell["height_cap_m"])
            if sigma > single_best:
                single_best = sigma
                single_cell = {
                    "Tp_s": cell["Tp_s"],
                    "gamma": cell["gamma"],
                    "f_ref_hz": f_ref,
                    "surface_Tz_s": cell["surface_Tz_s"],
                    "steepness_height_cap_m": cell["height_cap_m"],
                    "band_corners_hz": list(sigma_band_corners_hz(f_ref)),
                }
        own = min(
            max(1.0 / cell["deployed_Tz_s"], MIN_TUNE_FREQ_HZ), MAX_TUNE_FREQ_HZ
        )
        sigma = _sigma(
            _integrate(
                [
                    w * s
                    for w, s in zip(
                        sigma_band_weight(own),
                        up_specific_force_spectrum(
                            _partition_unit_spectrum(cell["Tp_s"], cell["gamma"])
                        ),
                    )
                ]
            ),
            cell["height_cap_m"],
        )
        if sigma > self_consistent_best:
            self_consistent_best = sigma
            self_consistent_cell = {
                "Tp_s": cell["Tp_s"],
                "gamma": cell["gamma"],
                "f_ref_hz": own,
                "deployed_Tz_s": cell["deployed_Tz_s"],
                "steepness_height_cap_m": cell["height_cap_m"],
            }

    coarse = [
        cell
        for index, cell in enumerate(cells)
        if index % MIXTURE_CELL_STRIDE == 0
    ]
    weights = _linear_span(0.0, 1.0, MIXTURE_WEIGHT_STEPS)
    steepness_scale = SIGNIFICANT_STEEPNESS_MAX * G_STD / (2.0 * math.pi)
    mixture_best = single_best
    mixture_cell: dict[str, Any] = {
        "second_partition": None,
        "note": "no admissible two-partition sea beat the one-partition maximum",
    }
    for f_ref in references:
        for i, first in enumerate(coarse):
            for second in coarse[i + 1:]:
                gain_a = first["gains"][f_ref]
                gain_b = second["gains"][f_ref]
                for weight in weights:
                    mixed_p = (
                        weight * first["inverse_square_period"]
                        + (1.0 - weight) * second["inverse_square_period"]
                    )
                    mixed_tz = 1.0 / math.sqrt(mixed_p)
                    total = min(
                        HS_MAX_M, steepness_scale * mixed_tz * mixed_tz
                    )
                    # Per-partition breaking limits cap each modal height too.
                    height_a = min(
                        first["height_cap_m"], total * math.sqrt(weight)
                    )
                    height_b = min(
                        second["height_cap_m"],
                        total * math.sqrt(1.0 - weight),
                    )
                    gain = (
                        height_a * height_a * gain_a
                        + height_b * height_b * gain_b
                    )
                    sigma = math.sqrt(gain) / 4.0
                    if sigma > mixture_best:
                        mixture_best = sigma
                        mixture_cell = {
                            "f_ref_hz": f_ref,
                            "first_partition": [
                                first["Tp_s"],
                                first["gamma"],
                                height_a,
                            ],
                            "second_partition": [
                                second["Tp_s"],
                                second["gamma"],
                                height_b,
                            ],
                            "mixture_surface_Tz_s": mixed_tz,
                            "mixture_Hs_m": math.sqrt(
                                height_a * height_a + height_b * height_b
                            ),
                        }

    prior_best = 0.0
    prior_cell: dict[str, Any] = {}
    prior_weight = sigma_band_weight(TUNE_FREQ_PRIOR_HZ)
    for cell in cells:
        source = up_specific_force_spectrum(
            _partition_unit_spectrum(cell["Tp_s"], cell["gamma"])
        )
        sigma = _sigma(
            _integrate([w * s for w, s in zip(prior_weight, source)]),
            cell["height_cap_m"],
        )
        if sigma > prior_best:
            prior_best = sigma
            prior_cell = {
                "Tp_s": cell["Tp_s"],
                "gamma": cell["gamma"],
                "steepness_height_cap_m": cell["height_cap_m"],
            }

    worst = max(single_best, mixture_best)
    return {
        "single_partition_max_sigma_a_ms2": single_best,
        "single_partition_limiting_cell": single_cell,
        "single_partition_maximum_is_exhaustive_on_grid": True,
        "multimodal_screen_max_sigma_a_ms2": mixture_best,
        "multimodal_screen_limiting_cell": mixture_cell,
        "multimodal_number_is_a_witness_search_not_an_optimum": True,
        "self_consistent_max_sigma_a_ms2": self_consistent_best,
        "self_consistent_limiting_cell": self_consistent_cell,
        "self_consistent_is_a_diagnostic_not_a_bound": True,
        "startup_prior_max_sigma_a_ms2": prior_best,
        "startup_prior_limiting_cell": prior_cell,
        "startup_prior_hz": TUNE_FREQ_PRIOR_HZ,
        "per_partition_and_mixture_steepness_both_imposed": True,
        "band_corners_at_channel_ceiling_hz": list(
            sigma_band_corners_hz(MAX_TUNE_FREQ_HZ)
        ),
        "absolute_band_ceiling_reachable_in_committed_channel": (
            sigma_band_corners_hz(MAX_TUNE_FREQ_HZ)[1] >= SIGMA_BAND_MAX_HZ
        ),
        "worst_admissible_sigma_a_ms2": worst,
        # sigma_target_ = min(sigma_coeff_ * sigma_a, max_sigma_a_).  The
        # ceiling saturates, so the source coordinate never leaves its declared
        # interval; the question is whether an admissible sea reaches the rail.
        "sigma_aw_map": "sigma_aw = min(0.9 * sigma_a, 4 m/s^2)",
        "sigma_aw_clamp_ms2": SIGMA_AW_CLAMP_MS2,
        "sigma_a_saturating_the_clamp_ms2": SIGMA_A_SATURATION_MS2,
        "worst_admissible_sigma_aw_ms2": min(
            SIGMA_COEFF * worst, SIGMA_AW_CLAMP_MS2
        ),
        "sigma_aw_coordinate_stays_inside_declared_interval": True,
        "clamp_rail_reachable_from_declared_sea_domain": (
            worst >= SIGMA_A_SATURATION_MS2
        ),
        "clamp_rail_reachable_at_settled_band_reference": (
            self_consistent_best >= SIGMA_A_SATURATION_MS2
        ),
        "clamp_rail_reachable_at_startup_prior": (
            prior_best >= SIGMA_A_SATURATION_MS2
        ),
        "rail_makes_sea_to_source_map_non_injective": (
            worst >= SIGMA_A_SATURATION_MS2
        ),
        "rho_max_used": RHO_MAX,
        "sigma_a_scales_linearly_in_rho_max": True,
        "rho_max_below_which_the_rail_is_unreachable": (
            RHO_MAX * SIGMA_A_SATURATION_MS2 / worst
        ),
        "sigma_a_scales_linearly_in_steepness_ceiling": True,
        "steepness_below_which_the_rail_is_unreachable": (
            SIGNIFICANT_STEEPNESS_MAX * SIGMA_A_SATURATION_MS2 / worst
        ),
        "promotion_use": False,
    }


def alias_obligation() -> dict[str, Any]:
    """Folded specific-force power at the deployed 5 ms sampling interval."""
    nyquist_rad_s = math.pi / IMU_DT_S
    peak_periods = _log_span(TP_MIN_S, TP_MAX_S, TP_SCAN)
    gammas = _linear_span(GAMMA_MIN, GAMMA_MAX, GAMMA_SCAN)

    worst_constant = 0.0
    worst_cell: dict[str, Any] = {}
    for peak_period in peak_periods:
        for gamma in gammas:
            unit = _partition_unit_spectrum(peak_period, gamma)
            height = steepness_height_cap_m(surface_zero_crossing_period(unit))
            constant = _asymptotic_acceleration_constant(
                height, peak_period, gamma
            )
            if constant > worst_constant:
                worst_constant = constant
                worst_cell = {
                    "Tp_s": peak_period,
                    "gamma": gamma,
                    "Hs_m": height,
                }

    bounded = acceleration_tail_above(
        nyquist_rad_s, worst_cell["Hs_m"], worst_cell["Tp_s"], worst_cell["gamma"]
    )
    flat_three = unbounded_response_alias_power(3.0, worst_constant)
    flat_six = unbounded_response_alias_power(6.0, worst_constant)
    return {
        "nyquist_rad_s": nyquist_rad_s,
        "worst_cell": worst_cell,
        "asymptotic_acceleration_constant_m2s3": worst_constant,
        "folded_power_with_declared_rolloff_m2s4": bounded,
        "folded_sigma_with_declared_rolloff_ms2": math.sqrt(bounded),
        "flat_response_folded_power_3_decades_m2s4": flat_three,
        "flat_response_folded_sigma_3_decades_ms2": math.sqrt(flat_three),
        "flat_response_folded_power_6_decades_m2s4": flat_six,
        "flat_response_folded_power_is_divergent": True,
        "flat_response_exceeds_sigma_clamp_within_decades": (
            math.sqrt(flat_three) > SIGMA_AW_CLAMP_MS2
        ),
        "response_rolloff_is_mandatory": True,
        "sampling_alone_band_limits_the_sea": False,
    }


def matrix_moment_report() -> dict[str, Any]:
    """Response-weighted matrix moments actually required by the source path."""
    representative = [
        {"label": "beam_swell_narrow", "beta_deg": 90.0, "spread_s": 25.0},
        {"label": "bow_quarter_broad", "beta_deg": 45.0, "spread_s": 1.0},
        {"label": "head_sea_typical", "beta_deg": 0.0, "spread_s": 8.0},
    ]
    rows = []
    worst_offdiagonal = 0.0
    for entry in representative:
        gram = directional_gram(
            math.radians(entry["beta_deg"]), entry["spread_s"]
        )
        normalized = 0.0
        for a in range(3):
            for b in range(3):
                if a == b:
                    continue
                denom = math.sqrt(gram[a][a] * gram[b][b])
                if denom > 0.0:
                    normalized = max(normalized, abs(gram[a][b]) / denom)
        worst_offdiagonal = max(worst_offdiagonal, normalized)
        rows.append(
            {
                "label": entry["label"],
                "beta_deg": entry["beta_deg"],
                "spread_s": entry["spread_s"],
                "directional_gram": gram,
                "mass": gram[0][0],
                "max_normalized_offdiagonal": normalized,
            }
        )

    # A single fully specified sea, to exhibit the finite matrix moments.
    sea = [(3.0, 8.0, 3.3, 20.0, 8.0), (1.5, 13.0, 1.0, 200.0, 25.0)]
    acceleration = [[0.0] * 3 for _ in range(3)]
    proxy = [[0.0] * 3 for _ in range(3)]
    for height, peak_period, gamma, beta_deg, spread_s in sea:
        elevation = partition_elevation_spectrum(height, peak_period, gamma)
        source = up_specific_force_spectrum(elevation)
        gram = directional_gram(math.radians(beta_deg), spread_s)
        accel_mass = _integrate(source)
        proxy_mass = _integrate(
            [w * s for w, s in zip(PROXY_ELEVATION_WEIGHT, source)]
        )
        for a in range(3):
            for b in range(3):
                acceleration[a][b] += accel_mass * gram[a][b]
                proxy[a][b] += proxy_mass * gram[a][b]

    return {
        "factorization": "M = sum_r (int omega-weighted rho^2 omega^4 S_r) Q_r",
        "directional_gram_cases": rows,
        "worst_normalized_offdiagonal": worst_offdiagonal,
        "per_axis_scalarization_valid": False,
        "example_sea": [list(entry) for entry in sea],
        "example_acceleration_matrix_moment_m2s4": acceleration,
        "example_proxy_elevation_matrix_moment_m2": proxy,
        "finiteness_ladder": {
            "proxy_elevation_and_velocity_moments": "finite for any rolloff p >= 0",
            "band_passed_acceleration_moment": "finite for any rolloff p >= 0",
            "raw_acceleration_moment": "finite iff p > 0",
            "acceleration_moment_of_order_n": "finite iff n < 2p",
            "deployed_path_requires_orders_above_zero": False,
        },
        "promotion_use": False,
    }


# ---------------------------------------------------------------------------
# Payload assembly.
# ---------------------------------------------------------------------------
def build() -> dict[str, Any]:
    period = period_channel_screen()
    convexity = mixture_convexity_check()
    sigma = sigma_channel_screen()
    alias = alias_obligation()
    moments = matrix_moment_report()

    payload: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "qualification": "OU3_SEA0_DIRECTIONAL_RESPONSE_MOMENT_SUBCERTIFICATE",
        "proof_status": "mixed_analytic_and_non_promoting_numerical_screen",
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_operating_domain_shrunk": False,
        "SEA0_full_certificate_promoted": False,
        "P2_promoted_from_this_artifact": False,
        "P3_promoted_from_this_artifact": False,
        "P4_promoted_from_this_artifact": False,
        "P5_promoted_from_this_artifact": False,
        "declared_sea_domain": {
            "m_max": M_MAX,
            "Tp_s": [TP_MIN_S, TP_MAX_S],
            "gamma": [GAMMA_MIN, GAMMA_MAX],
            "Hs_max_m": HS_MAX_M,
            "significant_steepness_max": SIGNIFICANT_STEEPNESS_MAX,
            "spreading_exponent_s": [SPREAD_S_MIN, SPREAD_S_MAX],
            "beta_rad": [-math.pi, math.pi],
            "steepness_is_physical_admissibility_not_proof_clamp": True,
            "covers_eight_reference_seas": True,
        },
        "response_family": {
            "envelope": "|G_j(omega,theta)| <= rho(omega) omega^2 u_j(theta)",
            "u_of_theta": ["1", "|cos theta|", "|sin theta|"],
            "rho_max": RHO_MAX,
            "rolloff_exponent_p": ROLL_OFF_EXPONENT,
            "reference_waterline_m": L_REF_M,
            "rolloff_corner_rad_s": OMEGA_L,
            "rolloff_corner_hz": OMEGA_L / (2.0 * math.pi),
            "mechanism": "finite-waterline Froude-Krylov length averaging",
            "lever_arm_enabled": False,
            "single_vessel_assumed": False,
        },
        "analytical_lemmas": {
            "modulus_majorant": (
                "for |G_j| <= m_j, v^H Phi_u v <= |v|^T M |v| with "
                "M = int m m^T E dtheta, hence lambda_max(Phi_u) <= lambda_max(M)"
            ),
            "theta_factorization": (
                "E = sum_r S_r(omega) D_r(theta) makes every response-weighted "
                "moment a sum of scalar frequency integrals times constant "
                "directional Gram matrices Q_r(beta_r,s_r)"
            ),
            "exact_leak_inversion": (
                "omega_est^2 = (sigma_v/sigma_eta)^2 - lambda^2 = E_mu[omega^2] "
                "with d mu = |G_eta|^2 S_up d omega; exact for any spectrum"
            ),
            "deployed_period_is_high_passed": (
                "d mu = rho^2 (omega^2/(omega^2+lambda^2))^4 S_eta d omega, so "
                "the deployed Tz is the moment period of a fourth-order "
                "high-passed, response-weighted elevation spectrum"
            ),
            "mixture_convexity": (
                "omega_est^2 is a convex combination over partitions, so the "
                "extreme deployed periods of the SEA3 class occur at "
                "one-partition seas"
            ),
            "leak_inversion_is_not_narrow_band": True,
            "unbanded_surface_acceleration_variance_finite": False,
            "sampling_alone_removes_acceleration_divergence": False,
        },
        "matrix_moments": moments,
        "deployed_period_channel": period,
        "mixture_convexity_check": convexity,
        "deployed_sigma_channel": sigma,
        "alias_obligation": alias,
        "quadrature": {
            "omega_interval_rad_s": [OMEGA_MIN, OMEGA_MAX],
            "omega_simpson_panels": OMEGA_PANELS,
            "theta_simpson_panels": THETA_PANELS,
            "interval_arithmetic_used_for_quadrature": False,
            "high_frequency_tail_bounded_analytically": True,
            "low_frequency_tail_reason": (
                "exp(-1.25 (omega_p/omega)^4) underflows below the grid floor"
            ),
        },
        "reason_not_promoted": (
            "fixed-grid quadrature with analytical pads is a screening "
            "enclosure, not validated interval integration; the deployed "
            "channels are also still evaluated in continuous-time steady state "
            "rather than through the exact discrete EMA and log-period state"
        ),
        "limiting_quantity": (
            "reachability of the sigma_aw clamp rail from the declared "
            "directional sea domain, where the sea-to-source map stops being "
            "injective and SEA3 can no longer prune the P2 language"
        ),
        "sea3_source_language_included_in_frozen_p2_contract": False,
        "next_obligation": (
            "decide the sigma clamp rail: either certify a sharper vessel/IMU "
            "response envelope that puts it out of reach, or carry the "
            "saturated non-injective cell explicitly in the SEA3 P2 language; "
            "then replace the continuous-time steady-state channel weights by "
            "the exact discrete WavePeriodEstimator moment EMAs, log-period "
            "state and AdaptiveWaveBandPass recursion"
        ),
    }

    failures = validate(payload)
    payload["validation_pass"] = not failures
    payload["validation_failures"] = failures
    return payload


def validate(payload: dict[str, Any]) -> list[str]:
    failures: list[str] = []

    if payload["schema_version"] != SCHEMA_VERSION:
        failures.append("schema version mismatch")
    for key in (
        "SEA0_full_certificate_promoted",
        "P2_promoted_from_this_artifact",
        "P3_promoted_from_this_artifact",
        "P4_promoted_from_this_artifact",
        "P5_promoted_from_this_artifact",
        "trajectory_replay_used",
        "filter_changed",
        "declared_operating_domain_shrunk",
    ):
        if payload[key]:
            failures.append(f"{key} must stay false")

    period = payload["deployed_period_channel"]
    ratio_lo, ratio_hi = period["deployed_Tz_over_Tp"]
    if not (0.0 < ratio_lo <= ratio_hi):
        failures.append("deployed Tz/Tp screen is not ordered")
    if not period["inside_committed_tuning_channel"]:
        failures.append(
            "declared sea domain leaves the committed 0.03-1.2 Hz channel"
        )

    convexity = payload["mixture_convexity_check"]
    if convexity["worst_relative_excess"] > 1e-9:
        failures.append("mixture period left the single-partition extremes")
    for case in convexity["cases"]:
        if not case["inside_component_extremes"]:
            failures.append("mixture convexity case violated")

    moments = payload["matrix_moments"]
    if moments["worst_normalized_offdiagonal"] <= 0.0:
        failures.append("directional Gram lost its cross-axis coupling")
    if moments["per_axis_scalarization_valid"]:
        failures.append("per-axis scalarization must stay invalid")
    for case in moments["directional_gram_cases"]:
        if abs(case["mass"] - 1.0) > 1e-6:
            failures.append(
                f"directional spreading is not normalized: {case['label']}"
            )
    for row in moments["example_acceleration_matrix_moment_m2s4"]:
        for value in row:
            if not math.isfinite(value):
                failures.append("acceleration matrix moment is not finite")

    alias = payload["alias_obligation"]
    if not math.isfinite(alias["folded_power_with_declared_rolloff_m2s4"]):
        failures.append("declared roll-off failed to bound the alias fold")
    if not alias["flat_response_folded_power_is_divergent"]:
        failures.append("flat-response alias divergence must stay recorded")
    if alias["sampling_alone_band_limits_the_sea"]:
        failures.append("sampling must not be claimed to band limit the sea")

    if payload["sea3_source_language_included_in_frozen_p2_contract"]:
        failures.append("SEA3 source-language inclusion is not established here")

    sigma = payload["deployed_sigma_channel"]
    if not (sigma["single_partition_max_sigma_a_ms2"] > 0.0):
        failures.append("sigma channel screen produced no positive bound")
    if (
        sigma["multimodal_screen_max_sigma_a_ms2"]
        < sigma["single_partition_max_sigma_a_ms2"]
    ):
        failures.append("multimodal screen fell below the one-partition maximum")
    if not sigma["per_partition_and_mixture_steepness_both_imposed"]:
        failures.append("both steepness admissibility levels must stay imposed")
    if sigma["worst_admissible_sigma_aw_ms2"] > SIGMA_AW_CLAMP_MS2:
        failures.append("clamped sigma_aw left its declared interval")
    if not sigma["sigma_aw_coordinate_stays_inside_declared_interval"]:
        failures.append("the sigma_aw clamp must keep the coordinate contained")
    if sigma["clamp_rail_reachable_from_declared_sea_domain"] != (
        sigma["worst_admissible_sigma_a_ms2"] >= SIGMA_A_SATURATION_MS2
    ):
        failures.append("clamp-rail reachability flag is inconsistent")
    if sigma["clamp_rail_reachable_at_settled_band_reference"] and not sigma[
        "clamp_rail_reachable_from_declared_sea_domain"
    ]:
        failures.append("settled rail reachability contradicts the sound bound")

    domain = payload["declared_sea_domain"]
    if domain["m_max"] != M_MAX:
        failures.append("partition count must remain M_max = 3")
    if not domain["steepness_is_physical_admissibility_not_proof_clamp"]:
        failures.append("steepness ceiling must stay a physical assumption")

    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare the committed artifact against a fresh build",
    )
    args = parser.parse_args()

    payload = build()
    text = json.dumps(payload, indent=2, sort_keys=True)

    committed = Path(__file__).with_suffix(".json")
    if args.check:
        if not committed.exists():
            print(f"missing committed artifact: {committed}", file=sys.stderr)
            return 1
        if json.loads(committed.read_text(encoding="utf-8")) != payload:
            print("committed artifact is stale", file=sys.stderr)
            return 1
        print("committed artifact is current")
        return 0 if payload["validation_pass"] else 1

    if args.output is not None:
        args.output.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0 if payload["validation_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
