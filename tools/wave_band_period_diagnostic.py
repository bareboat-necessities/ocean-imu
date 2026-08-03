#!/usr/bin/env python3
"""Compare acceleration-band and wave-band period proxies on released records.

This is an evidence-only diagnostic. It does not alter the deployed filter.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

G_STD = 9.80665
TAU_NUMERATOR = 1.38 * 0.5
ORACLE_TAU_JONSWAP = {
    "0.270": 1.081,
    "1.500": 1.565,
    "4.000": 2.350,
    "8.500": 2.330,
}


@dataclass
class Row:
    height_m: float
    record: str
    samples: int
    duration_sec: float
    elevation_peak_hz: float
    acceleration_peak_hz: float
    body_proxy_peak_hz: float
    elevation_mean_hz: float
    elevation_zero_cross_hz: float
    acceleration_mean_hz: float
    acceleration_zero_cross_hz: float
    body_proxy_mean_hz: float
    body_proxy_zero_cross_hz: float
    reconstructed_elevation_mean_hz: float
    reconstructed_elevation_zero_cross_hz: float
    body_noise_psd: float
    oracle_tau_sec: float
    tau_from_acceleration_peak_sec: float
    tau_from_body_peak_sec: float
    tau_from_body_mean_sec: float
    tau_from_body_zero_cross_sec: float
    tau_from_reconstructed_elevation_mean_sec: float
    tau_from_reconstructed_elevation_zero_cross_sec: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("plots/kalman_ou_ii"))
    parser.add_argument("--output-dir", type=Path, default=Path("reports/results/wave_band_period"))
    parser.add_argument("--window-sec", type=float, default=900.0)
    parser.add_argument("--decimated-hz", type=float, default=4.0)
    parser.add_argument("--welch-sec", type=float, default=128.0)
    return parser.parse_args()


def find_records(data_dir: Path) -> list[tuple[str, Path]]:
    found: list[tuple[str, Path]] = []
    for token in ORACLE_TAU_JONSWAP:
        matches = sorted(data_dir.glob(f"wave_data_jonswap_H{token}_*.csv"))
        if not matches:
            raise FileNotFoundError(f"missing JONSWAP H{token} record in {data_dir}")
        found.append((token, matches[0]))
    return found


def load_columns(path: Path, names: Iterable[str]) -> tuple[dict[str, np.ndarray], float]:
    with path.open("r", encoding="utf-8") as stream:
        header = stream.readline().strip().split(",")
    names = tuple(names)
    indices = [header.index(name) for name in names]
    values = np.loadtxt(path, delimiter=",", skiprows=1, usecols=indices, dtype=np.float64)
    if values.ndim == 1:
        values = values.reshape(1, -1)
    columns = {name: values[:, index] for index, name in enumerate(names)}
    return columns, float(np.median(np.diff(columns["time"])))


def block_decimate(values: np.ndarray, source_dt: float, target_hz: float) -> tuple[np.ndarray, float]:
    factor = max(1, int(round(1.0 / (source_dt * target_hz))))
    count = values.size // factor
    if count < 8:
        raise ValueError("record too short after decimation")
    return values[: count * factor].reshape(count, factor).mean(axis=1), source_dt * factor


def detrend_linear(values: np.ndarray) -> np.ndarray:
    x = np.linspace(-1.0, 1.0, values.size)
    design = np.column_stack((np.ones(values.size), x))
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    return values - design @ coefficients


def welch_psd(values: np.ndarray, fs: float, segment_samples: int) -> tuple[np.ndarray, np.ndarray]:
    segment_samples = min(segment_samples, values.size)
    segment_samples = max(64, segment_samples)
    if segment_samples > values.size:
        raise ValueError("record too short for Welch estimate")
    hop = max(1, segment_samples // 2)
    window = np.hanning(segment_samples)
    scale = fs * float(np.sum(window * window))
    spectra: list[np.ndarray] = []
    for start in range(0, values.size - segment_samples + 1, hop):
        segment = detrend_linear(values[start : start + segment_samples])
        fft = np.fft.rfft(segment * window)
        psd = (np.abs(fft) ** 2) / scale
        if psd.size > 2:
            psd[1:-1] *= 2.0
        spectra.append(psd)
    if not spectra:
        raise ValueError("no Welch segments")
    return np.fft.rfftfreq(segment_samples, 1.0 / fs), np.mean(spectra, axis=0)


def band_metrics(frequency: np.ndarray, psd: np.ndarray, f_min: float = 0.05, f_max: float = 0.65) -> tuple[float, float, float]:
    mask = (frequency >= f_min) & (frequency <= f_max) & np.isfinite(psd) & (psd > 0.0)
    f = frequency[mask]
    p = psd[mask]
    if f.size < 3 or float(np.sum(p)) <= 0.0:
        return math.nan, math.nan, math.nan
    peak = float(f[int(np.argmax(p))])
    m0 = float(np.trapz(p, f))
    m1 = float(np.trapz(f * p, f))
    m2 = float(np.trapz(f * f * p, f))
    return peak, m1 / m0, math.sqrt(m2 / m0)


def safe_tau(frequency_hz: float) -> float:
    return TAU_NUMERATOR / frequency_hz if frequency_hz > 0.0 and math.isfinite(frequency_hz) else math.nan


def evaluate_record(token: str, path: Path, window_sec: float, decimated_hz: float, welch_sec: float) -> Row:
    columns, source_dt = load_columns(path, ("time", "disp_z", "acc_z", "acc_bz"))
    count = min(columns["time"].size, int(round(window_sec / source_dt)))
    sl = slice(columns["time"].size - count, columns["time"].size)

    elevation, dt = block_decimate(columns["disp_z"][sl], source_dt, decimated_hz)
    acceleration, dt_acc = block_decimate(columns["acc_z"][sl], source_dt, decimated_hz)
    body_proxy, dt_body = block_decimate(-(columns["acc_bz"][sl] + G_STD), source_dt, decimated_hz)
    if abs(dt - dt_acc) > 1e-12 or abs(dt - dt_body) > 1e-12:
        raise AssertionError("decimation cadence mismatch")

    fs = 1.0 / dt
    segment_samples = int(round(welch_sec * fs))
    f_eta, p_eta = welch_psd(elevation, fs, segment_samples)
    f_acc, p_acc = welch_psd(acceleration, fs, segment_samples)
    f_body, p_body = welch_psd(body_proxy, fs, segment_samples)
    if not (np.array_equal(f_eta, f_acc) and np.array_equal(f_eta, f_body)):
        raise AssertionError("frequency grids differ")

    noise_mask = (f_body >= 0.8) & (f_body <= min(1.8, 0.9 * fs / 2.0))
    noise_psd = float(np.median(p_body[noise_mask])) if np.any(noise_mask) else 0.0
    body_signal_psd = np.maximum(p_body - noise_psd, 0.0)
    omega4 = np.maximum((2.0 * math.pi * f_body) ** 4, np.finfo(np.float64).tiny)
    reconstructed_eta_psd = body_signal_psd / omega4
    reconstructed_eta_psd[f_body < 0.05] = 0.0

    eta_peak, eta_mean, eta_zc = band_metrics(f_eta, p_eta)
    acc_peak, acc_mean, acc_zc = band_metrics(f_acc, p_acc)
    body_peak, body_mean, body_zc = band_metrics(f_body, body_signal_psd)
    _, reconstructed_mean, reconstructed_zc = band_metrics(f_body, reconstructed_eta_psd)

    return Row(
        height_m=float(token), record=path.name, samples=int(count), duration_sec=count * source_dt,
        elevation_peak_hz=eta_peak, acceleration_peak_hz=acc_peak, body_proxy_peak_hz=body_peak,
        elevation_mean_hz=eta_mean, elevation_zero_cross_hz=eta_zc,
        acceleration_mean_hz=acc_mean, acceleration_zero_cross_hz=acc_zc,
        body_proxy_mean_hz=body_mean, body_proxy_zero_cross_hz=body_zc,
        reconstructed_elevation_mean_hz=reconstructed_mean,
        reconstructed_elevation_zero_cross_hz=reconstructed_zc, body_noise_psd=noise_psd,
        oracle_tau_sec=ORACLE_TAU_JONSWAP[token],
        tau_from_acceleration_peak_sec=safe_tau(acc_peak),
        tau_from_body_peak_sec=safe_tau(body_peak),
        tau_from_body_mean_sec=safe_tau(body_mean),
        tau_from_body_zero_cross_sec=safe_tau(body_zc),
        tau_from_reconstructed_elevation_mean_sec=safe_tau(reconstructed_mean),
        tau_from_reconstructed_elevation_zero_cross_sec=safe_tau(reconstructed_zc),
    )


def summarize(rows: list[Row]) -> dict[str, object]:
    candidates = [name for name in asdict(rows[0]) if name.startswith("tau_from_")]
    oracle = np.array([row.oracle_tau_sec for row in rows])
    errors: dict[str, dict[str, float]] = {}
    for name in candidates:
        estimate = np.array([getattr(row, name) for row in rows])
        delta = estimate - oracle
        errors[name] = {
            "mae_sec": float(np.mean(np.abs(delta))),
            "rmse_sec": float(np.sqrt(np.mean(delta * delta))),
            "max_abs_sec": float(np.max(np.abs(delta))),
            "correlation": float(np.corrcoef(estimate, oracle)[0, 1]),
        }
    return {"tau_numerator": TAU_NUMERATOR, "candidate_errors": errors}


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = [evaluate_record(token, path, args.window_sec, args.decimated_hz, args.welch_sec)
            for token, path in find_records(args.data_dir)]
    summary = summarize(rows)

    csv_path = args.output_dir / "wave_band_period_diagnostic.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)
    json_path = args.output_dir / "wave_band_period_diagnostic.json"
    json_path.write_text(json.dumps({"rows": [asdict(row) for row in rows], "summary": summary},
                                    indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(csv_path.read_text(encoding="utf-8"), end="")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
