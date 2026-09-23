#!/usr/bin/env python3
"""Independent double-precision oracle for the scalar detrender default fixture.

The input timestamps/samples are retained verbatim. Only reviewed algorithm or
configuration changes justify ``--write``; ordinary CI only checks the fixture.
This model does not import, compile, or read the implementation under test.
"""

import argparse
import csv
import io
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/detrend/detrend-basic-data.csv"
FIELDS = ("x_axis", "original_cm", "baseline_slow", "wave_raw", "wave_clean",
          "wave_freq_hz", "wave_period_s", "baseline_cutoff_hz", "baseline_tau_s",
          "cleanup_cutoff_hz", "slope_rms", "slope_threshold", "freq_valid", "schmitt_state")


def reference_rows(rows):
    """Replay default scalar equations with log-frequency EMA in binary64."""
    if not rows:
        raise ValueError("Empty detrender fixture")
    frequency, baseline = 0.12, float(rows[0]["original_cm"])
    elapsed = slope = rms2 = cleanup_lp = 0.0
    state, accepted = 0, 0
    last_cross = {1: None, -1: None}
    last_valid = -math.inf
    prev_x, prev_time = baseline, float(rows[0]["x_axis"])
    outputs = []
    for index, row in enumerate(rows):
        x, timestamp = float(row["original_cm"]), float(row["x_axis"])
        if not (math.isfinite(x) and math.isfinite(timestamp)):
            raise ValueError(f"Non-finite input at row {index}")
        wave_raw = wave_clean = 0.0
        if index:
            if timestamp <= prev_time:
                raise ValueError(f"Non-increasing timestamp at row {index}")
            dt = min(2.0, max(1e-4, timestamp - prev_time))
            previous_slope = slope
            elapsed += dt
            alpha = math.exp(-dt / 0.20)
            slope = alpha * slope + (1.0 - alpha) * (x - prev_x) / dt
            alpha = math.exp(-dt / 8.0)
            rms2 = alpha * rms2 + (1.0 - alpha) * slope * slope
            threshold = min(1e9, max(0.001, 0.30 * math.sqrt(rms2)))
            # Same-polarity slope crossings provide measured period T. The
            # smoothing weight is exp(-T/13), not the previous held period.
            if elapsed >= 2.0:
                for sign in (1, -1):
                    if state != sign and sign * previous_slope < threshold <= sign * slope:
                        dy = slope - previous_slope
                        fraction = (sign * threshold - previous_slope) / dy if abs(dy) >= 1e-12 else 1.0
                        crossing = elapsed - dt + min(1.0, max(0.0, fraction)) * dt
                        previous = last_cross[sign]
                        if previous is not None:
                            period = crossing - previous
                            if 1.0 / 1.20 <= period <= 1.0 / 0.02:
                                measured = min(1.20, max(0.02, 1.0 / period))
                                alpha = math.exp(-max(period, 1e-4) / 13.0)
                                frequency = math.exp(alpha * math.log(frequency)
                                                     + (1.0 - alpha) * math.log(measured))
                                accepted += 1
                                last_valid = elapsed
                        last_cross[sign], state = crossing, sign
                        break
            cutoff = min(0.25, max(0.003, 0.35 * frequency))
            alpha = math.exp(-2.0 * math.pi * cutoff * dt)
            baseline = alpha * baseline + (1.0 - alpha) * x
            wave_raw = x - baseline
            # One cleanup high-pass stage at the baseline cutoff.
            cleanup_lp = alpha * cleanup_lp + (1.0 - alpha) * wave_raw
            wave_clean = wave_raw - cleanup_lp
        cutoff = min(0.25, max(0.003, 0.35 * frequency))
        rms = math.sqrt(rms2)
        values = (baseline, wave_raw, wave_clean, frequency, 1.0 / frequency,
                  cutoff, 1.0 / (2.0 * math.pi * cutoff), cutoff, rms,
                  min(1e9, max(0.001, 0.30 * rms)),
                  accepted >= 2 and elapsed - last_valid <= 3.0 / frequency, state)
        output = {"x_axis": row["x_axis"], "original_cm": row["original_cm"]}
        output.update(zip(FIELDS[2:], values))
        outputs.append(output)
        prev_time, prev_x = timestamp, x
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Explicitly replace expected outputs")
    args = parser.parse_args()
    with FIXTURE.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ValueError("Unexpected detrender fixture schema")
        rows = list(reader)
    expected = reference_rows(rows)
    if args.write:
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(buffer, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(expected)
        FIXTURE.write_text(buffer.getvalue(), encoding="utf-8")
    else:
        for index, (actual, reference) in enumerate(zip(rows, expected)):
            for key in FIELDS[2:]:
                value = reference[key]
                matches = (actual[key] == str(value) if isinstance(value, bool) else
                           math.isclose(float(actual[key]), value, rel_tol=1e-10, abs_tol=1e-10))
                if not matches:
                    raise ValueError(f"Stale detrender reference: row {index}, {key}")
    print(f"Scalar detrender reference: {len(expected)} rows verified; inputs preserved")


if __name__ == "__main__":
    main()
