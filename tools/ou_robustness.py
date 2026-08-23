#!/usr/bin/env python3
"""Current OU-III robustness study using the deployed analytical SpectralMSE coupling.

The historical implementation is kept in ``ou_robustness_core`` so committed
legacy evidence can still be restated byte-for-byte. New simulator runs use the
same OFAT and degradation protocols, but the two coupled sensitivity directions
follow the deployed SpectralMSE law rather than the retired cubic relation.
"""
from __future__ import annotations

import json
import math
import re
from dataclasses import replace
from pathlib import Path
from typing import Any, Mapping, Sequence

import ou_robustness_core as _core
from ou_robustness_core import *  # noqa: F401,F403

_ORIGINAL_SCALED_TUNING_POINT = _core.scaled_tuning_point
_ORIGINAL_WRITE_SENSITIVITY_PLOT = _core.write_sensitivity_plot
_ORIGINAL_WRITE_PUBLICATION_TABLE = _core.write_publication_table
_ORIGINAL_CODE_SOURCE_PATHS = tuple(_core.CODE_SOURCE_PATHS)

# The source uses T_S = clip((0.015/1.1) tau, 0.005, 0.25) at the default
# 200-Hz schedule. Keep the exact cadence in the coupled tau perturbation so
# the study remains correct when a scale crosses a cadence clamp.
_PSEUDO_TAU_RATIO = 0.015 / 1.1
_PSEUDO_MIN_S = 0.005
_PSEUDO_MAX_S = 0.25


def _pseudo_period(tau_s: float) -> float:
    return min(max(_PSEUDO_TAU_RATIO * tau_s, _PSEUDO_MIN_S), _PSEUDO_MAX_S)


def spectral_mse_tau_ratio(old_tau_s: float, new_tau_s: float) -> float:
    """Exact r_S(new)/r_S(old) for a tau-only SpectralMSE perturbation."""
    if not (old_tau_s > 0.0 and new_tau_s > 0.0):
        raise ValueError("tau values must be positive")
    return (
        (new_tau_s / old_tau_s) ** (24.0 / 7.0)
        * math.sqrt(_pseudo_period(old_tau_s) / _pseudo_period(new_tau_s))
    )


def scaled_tuning_point(
    baseline: Any,
    parameter: str,
    scale: float,
) -> Any:
    """Apply OFAT or deployed-SpectralMSE-coupled sensitivity scaling."""
    if parameter not in SENSITIVITY_PARAMETERS:
        raise ValueError(f"unknown sensitivity parameter: {parameter}")
    if not (math.isfinite(scale) and scale > 0.0):
        raise ValueError("sensitivity scale must be positive and finite")
    if baseline.RS_ms is None:
        raise ValueError("OU-III sensitivity requires an r_S tuning value")

    if parameter == "tau":
        point = replace(baseline, tau_s=baseline.tau_s * scale)
    elif parameter == "sigma_aw":
        point = replace(baseline, sigma_a_mps2=baseline.sigma_a_mps2 * scale)
    elif parameter == "r_s":
        point = replace(baseline, RS_ms=baseline.RS_ms * scale)
    elif parameter == "sigma_aw_rs":
        # c_sigma is fixed in this controlled study, so a multiplier on sigma_aw
        # is the same multiplier on the physical sigma_a,B entering SpectralMSE.
        sigma_new = min(
            max(baseline.sigma_a_mps2 * scale, SIGMA_AW_BOUNDS_MPS2[0]),
            SIGMA_AW_BOUNDS_MPS2[1],
        )
        realized = sigma_new / baseline.sigma_a_mps2
        point = replace(
            baseline,
            sigma_a_mps2=sigma_new,
            RS_ms=baseline.RS_ms * realized ** (6.0 / 7.0),
        )
    else:
        # SpectralMSE: r_S ~ tau^(24/7) / sqrt(T_S(tau)). Away from cadence
        # clamps this reduces to tau^(41/14), but the study uses the exact ratio.
        tau_new = min(
            max(baseline.tau_s * scale, TAU_BOUNDS_S[0]),
            TAU_BOUNDS_S[1],
        )
        point = replace(
            baseline,
            tau_s=tau_new,
            RS_ms=baseline.RS_ms * spectral_mse_tau_ratio(baseline.tau_s, tau_new),
        )
    return _core._bound_tuning_point(point)


def write_sensitivity_plot(path: Path, summary: Sequence[Mapping[str, Any]]) -> None:
    labels = {
        "tau": r"$\tau$ only",
        "sigma_aw": r"$\sigma_{aw}$ only",
        "r_s": r"$r_S$ only ($R_S=\mathrm{diag}(r_S^2)$)",
        "sigma_aw_rs": r"$\sigma_{aw}$ + SpectralMSE $r_S$",
        "tau_rs": r"$\tau$ + SpectralMSE $r_S$",
    }
    colors = {
        "tau": "#0072B2",
        "sigma_aw": "#D55E00",
        "r_s": "#009E73",
        "sigma_aw_rs": "#CC79A7",
        "tau_rs": "#56B4E9",
    }
    styles = {
        "tau": "-",
        "sigma_aw": "-",
        "r_s": "-",
        "sigma_aw_rs": "--",
        "tau_rs": "--",
    }
    metrics = (
        ("disp_z_pct_hs", "Vertical reconstruction", r"RMS error (\% $H_s$)"),
        ("disp_3d_rms_m", "Full 3D displacement", "RMS error (m)"),
        ("pitch_rms_deg", "Pitch", "RMS error (deg)"),
    )
    fig, axes = _core.plt.subplots(1, 3, figsize=(10.2, 3.35))
    handles = []
    for axis, (metric, title, ylabel) in zip(axes, metrics):
        for parameter in SENSITIVITY_PARAMETERS:
            rows = sorted(
                (
                    row
                    for row in summary
                    if row["experiment"] == "sensitivity"
                    and row["parameter"] == parameter
                    and row["metric"] == metric
                ),
                key=lambda row: float(row["scale_label"]),
            )
            x = _core.np.asarray([float(row["scale_label"]) for row in rows])
            y = _core.np.asarray([float(row["mean"]) for row in rows])
            low = y - _core.np.asarray(
                [float(row["bootstrap_ci95_low"]) for row in rows]
            )
            high = _core.np.asarray(
                [float(row["bootstrap_ci95_high"]) for row in rows]
            ) - y
            handle = axis.errorbar(
                x,
                y,
                yerr=_core.np.vstack((low, high)),
                marker="o",
                linewidth=1.6,
                capsize=2.5,
                color=colors[parameter],
                linestyle=styles[parameter],
                label=labels[parameter],
            )
            if axis is axes[0]:
                handles.append(handle)
        axis.axvline(1.0, color="#555555", linewidth=0.9, linestyle="--")
        axis.set_xlabel("Nominal multiplier")
        axis.set_ylabel(ylabel)
        axis.set_title(title)
        axis.grid(True, alpha=0.25)
    fig.legend(
        handles,
        [labels[parameter] for parameter in SENSITIVITY_PARAMETERS],
        frameon=False,
        ncol=3,
        fontsize=8,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.01),
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.89))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format="svg", metadata={"Date": None})
    _core.plt.close(fig)


def write_publication_table(
    path: Path,
    summary: Sequence[Mapping[str, Any]],
    effects: Sequence[Mapping[str, Any]],
    macros_path: Path | None = None,
) -> None:
    """Write the standard tables, then replace only the retired coupling wording."""
    _ORIGINAL_WRITE_PUBLICATION_TABLE(path, summary, effects, macros_path=macros_path)
    text = path.read_text(encoding="utf-8")
    caption = (
        r"\caption{OU--III tuning sensitivity at the nominal $H_s=1.5$ m sea. "
        r"Each entry is vertical-displacement RMS error in percent of $H_s$ "
        r"(mean $\pm$ sample standard deviation, $n=\OURobustnessPairs$ paired "
        r"seed triplets). The first three columns hold the other two parameters "
        r"fixed. The last two couple $r_S$ to the perturbed parameter through "
        r"the deployed SpectralMSE law. A $\sigma_{aw}$ multiplier $c$ gives "
        r"$r_S\to c^{6/7}r_S$; a $\tau$ perturbation uses "
        r"$r_S\to c^{24/7}[T_S(\tau)/T_S(c\tau)]^{1/2}r_S$, which reduces to "
        r"$c^{41/14}r_S$ away from cadence clamps. The $r_S$-only multiplier "
        r"acts on the pseudo-measurement standard deviation, so its covariance "
        r"scales quadratically.}"
    )
    pattern = (
        r"\\caption\{OU--III tuning sensitivity.*?\}\n"
        r"\s*\\label\{tab:ou_robustness_sensitivity\}"
    )
    replacement = caption + "\n  " + r"\label{tab:ou_robustness_sensitivity}"
    text, count = re.subn(
        pattern,
        lambda _match: replacement,
        text,
        count=1,
        flags=re.S,
    )
    if count != 1:
        raise RuntimeError("robustness sensitivity caption anchor moved")
    text = text.replace(
        r"$\tau,r_S\propto\tau^3$",
        r"$\tau,r_S$ (SpectralMSE)",
    )
    path.write_text(text, encoding="utf-8")


def _activate_current_study() -> None:
    _core.scaled_tuning_point = scaled_tuning_point
    _core.write_sensitivity_plot = write_sensitivity_plot
    _core.write_publication_table = write_publication_table
    paths = [Path(__file__).resolve(), Path(_core.__file__).resolve()]
    for item in _ORIGINAL_CODE_SOURCE_PATHS:
        resolved = Path(item).resolve()
        if resolved not in paths:
            paths.append(resolved)
    _core.CODE_SOURCE_PATHS = tuple(paths)


def _activate_legacy_study() -> None:
    """Restore the historical core hooks for a legacy archived ensemble."""
    _core.scaled_tuning_point = _ORIGINAL_SCALED_TUNING_POINT
    _core.write_sensitivity_plot = _ORIGINAL_WRITE_SENSITIVITY_PLOT
    _core.write_publication_table = _ORIGINAL_WRITE_PUBLICATION_TABLE
    _core.CODE_SOURCE_PATHS = _ORIGINAL_CODE_SOURCE_PATHS


def _rows_use_current_spectral_mse(rows: Sequence[Mapping[str, Any]]) -> bool:
    """Identify the coupling from replayed tuning values, not editorial text.

    The legacy study coupled ``r_S`` linearly to ``sigma_aw``. The deployed
    SpectralMSE study uses the realized ``sigma_aw`` ratio to the 6/7 power.
    A non-unit coupled sensitivity row therefore distinguishes the two archived
    ensembles without trusting a caption or the current source revision.
    """
    by_repetition: dict[int, list[Mapping[str, Any]]] = {}
    for row in rows:
        if row.get("experiment") != "sensitivity" or row.get("parameter") != "sigma_aw_rs":
            continue
        by_repetition.setdefault(int(row["repetition"]), []).append(row)

    for repetition_rows in by_repetition.values():
        reference = next(
            (
                row
                for row in repetition_rows
                if math.isclose(float(row["scale_label"]), 1.0, abs_tol=1e-12)
            ),
            None,
        )
        if reference is None:
            continue
        sigma_ref = float(reference["configured_sigma_aw_mps2"])
        rs_ref = float(reference["configured_r_s_ms"])
        if not (sigma_ref > 0.0 and rs_ref > 0.0):
            continue
        for row in repetition_rows:
            if math.isclose(float(row["scale_label"]), 1.0, abs_tol=1e-12):
                continue
            sigma_ratio = float(row["configured_sigma_aw_mps2"]) / sigma_ref
            rs_ratio = float(row["configured_r_s_ms"]) / rs_ref
            current_ratio = sigma_ratio ** (6.0 / 7.0)
            legacy_ratio = sigma_ratio
            separation = abs(current_ratio - legacy_ratio)
            if separation <= 1e-10:
                continue
            tolerance = 1e-7 * max(1.0, abs(rs_ratio), abs(current_ratio))
            if abs(rs_ratio - current_ratio) <= tolerance:
                return True
            if abs(rs_ratio - legacy_ratio) <= tolerance:
                return False
            raise ValueError(
                "archived sigma_aw/r_S sensitivity rows match neither the "
                "legacy cubic coupling nor the deployed SpectralMSE coupling"
            )
    raise ValueError("cannot identify robustness coupling from archived sensitivity rows")


def restat_bundle(
    source: Path,
    output_dir: Path,
    bootstrap_resamples: int,
    stats_seed: int,
) -> int:
    """Restate using the coupling that produced the archived replay rows."""
    with source.open(encoding="utf-8") as stream:
        bundle = json.load(stream)
    rows = bundle.get("raw_runs", [])
    if not rows:
        # Preserve the core provenance/error ordering. In particular, a bare
        # empty JSON file must fail as an unprovenanced source bundle before any
        # coupling inference is attempted.
        return _core.restat_bundle(
            source,
            output_dir,
            bootstrap_resamples=bootstrap_resamples,
            stats_seed=stats_seed,
        )
    if _rows_use_current_spectral_mse(rows):
        _activate_current_study()
    else:
        _activate_legacy_study()
    return _core.restat_bundle(
        source,
        output_dir,
        bootstrap_resamples=bootstrap_resamples,
        stats_seed=stats_seed,
    )


def main(argv: list[str] | None = None) -> int:
    _activate_current_study()
    return _core.main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
