#!/usr/bin/env python3
"""Summarise an AtomS3R INS diagnostic capture (SEA_STATE_SERIAL_DIAG=1).

    python3 tools/ins_diag_report.py bench.log [--csv state.csv]

Reads the line types described in src/AtomS3R/AtomS3R_InsDiag.h and reports,
per boot:

* sensor configuration as read back from the BMI270/BMM150 (I line) and the
  calibration the sketch loaded (C line);
* accepted-frame timing from the R lines: interval statistics, gaps longer
  than the magnetic rotation history accepts (20 ms) and missing R lines;
* magnetometer acquisition cadence from the R lines' acquisition records;
* loop health and estimator state from the D lines: time to Live / north lock,
  startup gravity-gate trust, attitude-covariance re-seeds, published heading
  against the tilt-compensated compass on the same tilt (north/south swaps),
  large-tilt episodes, and heave or envelope values that stop changing while
  frames keep arriving.

--csv writes every D line as one CSV row (union of keys) for plotting.
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
import struct
import sys
from dataclasses import dataclass, field

FREEZE_MIN_S = 3.0          # value unchanged this long while frames arrive
HEADING_SWAP_DEG = 90.0     # published heading vs compass on the same tilt
HEADING_WARN_DEG = 20.0
LARGE_TILT_DEG = 60.0
ROTATION_HISTORY_GAP_US = 20000  # src/util/MagneticRotation.h MaxGapUs


def bits_to_float(text: str) -> float:
    return struct.unpack("<f", struct.pack("<I", int(text, 16)))[0]


def wrap180(deg: float) -> float:
    return (deg + 180.0) % 360.0 - 180.0


def as_float(value: str | None) -> float:
    try:
        return float(value) if value is not None else math.nan
    except ValueError:
        return math.nan


@dataclass
class Boot:
    seqs: list[int] = field(default_factory=list)
    sample_us: list[int] = field(default_factory=list)      # unwrapped
    gyro_dps: list[float] = field(default_factory=list)
    mag_frame_us: list[int] = field(default_factory=list)   # one per new acquisition
    mag_read_us: list[int] = field(default_factory=list)
    mag_nan_frames: int = 0
    state: list[dict[str, str]] = field(default_factory=list)
    events: list[str] = field(default_factory=list)
    registers: dict[str, str] = field(default_factory=dict)
    calibration: tuple[str, str] | None = None
    sketch: str = ""
    fifo_setup: dict[str, str] = field(default_factory=dict)   # F line (FIFO timing sketch)
    passes: list[dict[str, str]] = field(default_factory=list)  # P lines
    fits: list[dict[str, str]] = field(default_factory=list)    # X lines


def parse(lines) -> list[Boot]:
    boots: list[Boot] = []
    pre = Boot()          # H/I/C/E lines seen before the first frame of a boot
    last_mag_seq = None
    for raw in lines:
        line = raw.rstrip("\r\n")
        if len(line) < 2 or line[1] != ",":
            continue
        tag, tokens = line[0], line.split(",")[1:]
        if tag == "R":
            if len(tokens) != 19:
                continue
            seq = int(tokens[0], 16)
            if seq == 0 or not boots:
                boots.append(pre)
                pre = Boot()
                last_mag_seq = None
            boot = boots[-1]
            us = int(tokens[1], 16)
            if boot.sample_us:
                us = boot.sample_us[-1] + ((us - boot.sample_us[-1]) & 0xFFFFFFFF)
            boot.seqs.append(seq)
            boot.sample_us.append(us)
            w = [bits_to_float(t) for t in tokens[6:9]]
            boot.gyro_dps.append(math.degrees(math.sqrt(sum(x * x for x in w))))
            m = [bits_to_float(t) for t in tokens[9:12]]
            if any(math.isnan(x) for x in m):
                boot.mag_nan_frames += 1
            have, mag_seq = int(tokens[13], 16), int(tokens[15], 16)
            if have and mag_seq != last_mag_seq:
                last_mag_seq = mag_seq
                boot.mag_frame_us.append(int(tokens[16], 16))
                boot.mag_read_us.append((int(tokens[18], 16) - int(tokens[17], 16)) & 0xFFFFFFFF)
            continue
        target = boots[-1] if boots else pre
        if tag == "E" and "seq=0" in tokens and boots and boots[-1].seqs:
            target = pre  # logged by the next boot before its first frame
        if tag == "D":
            target.state.append(dict(t.split("=", 1) for t in tokens if "=" in t))
        elif tag == "I":
            target.registers = dict(t.split("=", 1) for t in tokens if "=" in t)
        elif tag == "C" and len(tokens) == 3 and target.calibration is None:
            target.calibration = (tokens[0], tokens[1])
        elif tag == "H":
            target.sketch = dict(t.split("=", 1) for t in tokens if "=" in t).get("sketch", "")
        elif tag == "E":
            target.events.append(",".join(tokens))
        elif tag in "FPX":
            values = dict(t.split("=", 1) for t in tokens if "=" in t)
            if tag == "F" and not target.fifo_setup:
                target.fifo_setup = values
            elif tag == "P":
                target.passes.append(values)
            elif tag == "X":
                target.fits.append(values)
    if not boots and (pre.state or pre.passes):
        boots.append(pre)  # a replay's output: state lines without R frames
    return boots


def episodes(times: list[float], flags: list[bool], min_s: float = 0.0) -> list[tuple[float, float]]:
    out, start = [], None
    for t, f in zip(times, flags):
        if f and start is None:
            start = t
        elif not f and start is not None:
            if t - start >= min_s:
                out.append((start, t))
            start = None
    if start is not None and times and times[-1] - start >= min_s:
        out.append((start, times[-1]))
    return out


def frozen(times: list[float], values: list[float], frames: list[float], min_s: float) -> list[tuple[float, float]]:
    """Spans where a value repeats exactly while frames keep arriving."""
    flags = [False] * len(values)
    for i in range(1, len(values)):
        same = values[i] == values[i - 1] and not math.isnan(values[i])
        flags[i] = same and frames[i] > 0
    return episodes(times, flags, min_s)


def fmt_spans(spans: list[tuple[float, float]], limit: int = 6) -> str:
    text = ", ".join(f"{a:.1f}-{b:.1f}s" for a, b in spans[:limit])
    return text + (f" (+{len(spans) - limit} more)" if len(spans) > limit else "")


def report_boot(index: int, boot: Boot, out) -> None:
    w = out.write
    w(f"\n=== boot {index}: {boot.sketch or 'unknown sketch'}\n")
    if boot.registers and boot.registers.get("bmi270") != "1":
        w("sensor: BMI270 registers not read (host replay or another IMU)\n")
    elif boot.registers:
        r = boot.registers
        w("sensor: acc_odr {} Hz (ACC_CONF {}), gyr_odr {} Hz (GYR_CONF {}), aux_odr {} Hz, bmm_odr {} Hz,"
          " bmm reps xy={} z={}, mag_source {}\n".format(
              r.get("acc_odr_hz", "?"), r.get("acc_conf", "?"), r.get("gyr_odr_hz", "?"),
              r.get("gyr_conf", "?"), r.get("aux_odr_hz", "?"), r.get("bmm_odr_hz", "?"),
              r.get("bmm_rep_xy", "?"), r.get("bmm_rep_z", "?"), r.get("mag_source", "?")))
    if boot.calibration:
        w(f"calibration: stored blob {'yes' if boot.calibration[0] == '1' else 'NO (raw sensor values)'}\n")
    for e in boot.events:
        w(f"event: {e}\n")
    if boot.fifo_setup:
        f = boot.fifo_setup
        w(f"FIFO acquisition setup: {'ok' if f.get('ok') == '1' else 'FAILED (' + f.get('failure', '?') + ')'};"
          f" ACC_CONF {f.get('acc_conf', '?')}, GYR_CONF {f.get('gyr_conf', '?')}, AUX_CONF {f.get('aux_conf', '?')},"
          f" trigger-to-readout {f.get('aux_offset_ms', '?')} ms, magnetometer calibration"
          f" {'applied' if f.get('mag_cal') == '1' else 'not applied'}\n")
    if boot.passes:
        p = boot.passes[-1]
        def flag(key: str) -> str:
            return "PASS" if p.get(key) == "1" else "fail"
        w(f"pass criteria after {p.get('window_s', '?')} s:\n"
          f"  rate     {flag('rate_ok')}  {p.get('rate_hz')} Hz, lost {p.get('lost')}\n"
          f"  no loss  {flag('no_loss_ok')}  {p.get('stalls')} display bursts, longest {p.get('stall_max_ms')} ms,"
          f" FIFO peak {p.get('fifo_max_bytes')} B\n"
          f"  spacing  {flag('spacing_ok')}  {p.get('spacing_min_us')}..{p.get('spacing_max_us')} us,"
          f" sensor clock {p.get('clock_ppm')} ppm\n"
          f"  mag      {flag('mag_ok')}  {p.get('aux_hz')} Hz, data-ready failures {p.get('drdy_fail')},"
          f" repeats {p.get('mag_repeat')}, off-grid intervals {p.get('aux_interval_bad')}\n"
          f"  delay    {flag('delay_ok')}  {p.get('delay_ms')} ms over {p.get('delay_fits')} fits,"
          f" spread {p.get('delay_spread_ms')} ms\n"
          f"  overall  {flag('all_ok')}\n")
    if boot.fits:
        ok = [as_float(x.get("delay_ms")) for x in boot.fits if x.get("ok") == "1"]
        reasons = sorted({x.get("reason", "?") for x in boot.fits if x.get("ok") != "1"})
        w(f"delay fits: {len(ok)} of {len(boot.fits)} valid"
          f"{', ' + ' '.join(f'{d:.0f}' for d in ok[:12]) + ' ms' if ok else ''}"
          f"{'; rejected: ' + ', '.join(reasons) if reasons else ''}\n")

    n = len(boot.seqs)
    if n:
        missing = sum(max(0, b - a - 1) for a, b in zip(boot.seqs, boot.seqs[1:]))
        start = "boot" if boot.seqs[0] == 0 else f"seq {boot.seqs[0]} (capture began after boot)"
        dts = [b - a for a, b in zip(boot.sample_us, boot.sample_us[1:])]
        span_s = (boot.sample_us[-1] - boot.sample_us[0]) * 1e-6
        w(f"frames: {n} from {start}, {span_s:.1f} s, {n / span_s if span_s > 0 else math.nan:.1f} frames/s,"
          f" missing R lines {missing}\n")
        if dts:
            q = statistics.quantiles(dts, n=100) if len(dts) >= 100 else [min(dts)] * 99
            w(f"frame interval us: min {min(dts)}  p1 {q[0]:.0f}  median {statistics.median(dts):.0f}"
              f"  p99 {q[98]:.0f}  max {max(dts)};  >{ROTATION_HISTORY_GAP_US // 1000} ms: "
              f"{sum(d > ROTATION_HISTORY_GAP_US for d in dts)},  >50 ms: {sum(d > 50000 for d in dts)}\n")
        fast = sum(g > 30.0 for g in boot.gyro_dps)
        w(f"gyro: max {max(boot.gyro_dps):.1f} dps, frames above 30 dps {100.0 * fast / n:.1f}%\n")
    if len(boot.mag_frame_us) > 1:
        mdt = [((b - a) & 0xFFFFFFFF) * 1e-3 for a, b in zip(boot.mag_frame_us, boot.mag_frame_us[1:])]
        w(f"magnetometer: {len(boot.mag_frame_us)} acquisitions, interval ms median {statistics.median(mdt):.1f}"
          f" max {max(mdt):.1f}, read ms median {statistics.median(boot.mag_read_us) * 1e-3:.2f},"
          f" frames without field {boot.mag_nan_frames}\n")

    st = boot.state
    if not st:
        w("no D lines\n")
        return
    t = [as_float(d.get("t_ms")) * 1e-3 for d in st]
    frames = [as_float(d.get("frames")) for d in st]
    w(f"loop: overruns {sum(as_float(d.get('overruns')) for d in st):.0f},"
      f" max tick {max(as_float(d.get('tick_max_us')) for d in st):.0f} us,"
      f" max frame gap {max(as_float(d.get('dt_max_us')) for d in st):.0f} us,"
      f" diagnostic lines dropped {st[-1].get('lines_dropped', '?')}\n")

    ready_key = next((k for k in ("live", "init", "seeded") if k in st[0]), None)
    first = None
    if ready_key:
        ready = [d.get(ready_key) == "1" for d in st]
        first = next((ti for ti, r in zip(t, ready) if r), None)
        lost = [b for a, b in episodes(t, [not r for r in ready]) if first is not None and a > first]
        w(f"{ready_key}: {'never' if first is None else f'at {first:.1f} s'}"
          f"{f', lost again at {fmt_spans([(x, x) for x in lost])}' if lost else ''}\n")
        if "north" in st[0]:
            north = next((ti for ti, d in zip(t, st) if d.get("north") == "1"), None)
            w(f"north lock: {'never' if north is None else f'at {north:.1f} s'}\n")
        if "gate_trusted" in st[0]:
            pre_live = [d for d, r in zip(st, ready) if not r] or st
            trusted = sum(d.get("gate_trusted") == "1" for d in pre_live)
            fast = sum(as_float(d.get("gyro_fast")) for d in pre_live)
            nfr = sum(as_float(d.get("frames")) for d in pre_live)
            w(f"before {ready_key}: gravity gate trusted {100.0 * trusted / len(pre_live):.0f}% of the time,"
              f" frames above 30 dps {100.0 * fast / nfr if nfr else math.nan:.0f}%,"
              f" max gate hold {max(as_float(d.get('gate_good_s')) for d in pre_live):.1f} s\n")
    if "reseeds" in st[0]:
        counts = [as_float(d.get("reseeds")) for d in st]
        times = [ti for ti, a, b in zip(t[1:], counts, counts[1:]) if b > a]
        w(f"attitude covariance re-seeds (yaw sigma jump > 30 deg): {len(times)}"
          f"{' at ' + ', '.join(f'{x:.1f} s' for x in times[:10]) if times else ''}\n")
    if "mag_applied" in st[0]:
        tries = sum(as_float(d.get("mag_tries")) for d in st)
        applied = sum(as_float(d.get("mag_applied")) for d in st)
        w(f"magnetic corrections: {applied:.0f} applied of {tries:.0f} offered\n")

    hdg = [as_float(d.get("hdg")) for d in st]
    cmp_ = [as_float(d.get("hdg_cmp")) for d in st]
    ok = [d.get("hdg_ok") == "1" for d in st]
    err = [abs(wrap180(a - b)) if o and not (math.isnan(a) or math.isnan(b)) else math.nan
           for a, b, o in zip(hdg, cmp_, ok)]
    finite = [e for e in err if not math.isnan(e)]
    if finite:
        swaps = episodes(t, [e > HEADING_SWAP_DEG for e in err])
        warn = sum(e > HEADING_WARN_DEG for e in finite)
        w(f"heading vs tilt-compensated compass: median {statistics.median(finite):.1f} deg, max {max(finite):.1f}"
          f" deg, >{HEADING_WARN_DEG:.0f} deg {100.0 * warn / len(finite):.1f}% of valid samples\n")
        if swaps:
            w(f"  NORTH/SOUTH-SCALE disagreement (>{HEADING_SWAP_DEG:.0f} deg): {fmt_spans(swaps)}\n")
    roll = [as_float(d.get("roll")) for d in st]
    pitch = [as_float(d.get("pitch")) for d in st]
    tilt = episodes(t, [abs(r) > LARGE_TILT_DEG or abs(p) > LARGE_TILT_DEG for r, p in zip(roll, pitch)])
    if tilt:
        w(f"large tilt (>{LARGE_TILT_DEG:.0f} deg roll or pitch): {fmt_spans(tilt)}\n")
        for _, b in tilt[:6]:
            after = [e for ti, e in zip(t, err) if b + 5.0 <= ti <= b + 30.0 and not math.isnan(e)]
            if after:
                w(f"  heading vs compass 5-30 s after {b:.1f} s: max {max(after):.1f} deg\n")
    # Outputs are held before the filter is ready and its envelope starts at
    # zero; judge freezes only from ten seconds after it is ready.
    since = first + 10.0 if first is not None else math.inf
    for key in ("heave", "heave_clean", "env", "a_vert", "p"):
        if key in st[0]:
            spans = [s_ for s_ in frozen(t, [as_float(d.get(key)) for d in st], frames, FREEZE_MIN_S)
                     if s_[0] >= since]
            if spans:
                w(f"{key} unchanged for >= {FREEZE_MIN_S:.0f} s while frames arrive: {fmt_spans(spans)}\n")
    for key, limit in (("heave", 19.9), ("p", 19.9)):
        if key in st[0]:
            clamp = episodes(t, [abs(as_float(d.get(key))) >= limit for d in st])
            if clamp:
                w(f"{key} at the 20 m clamp: {fmt_spans(clamp)}\n")


def write_csv(boots: list[Boot], path: str) -> None:
    keys: list[str] = ["boot"]
    for b in boots:
        for d in b.state:
            keys.extend(k for k in d if k not in keys)
    with open(path, "w", newline="") as f:
        out = csv.DictWriter(f, fieldnames=keys)
        out.writeheader()
        for i, b in enumerate(boots):
            for d in b.state:
                out.writerow({"boot": i, **d})


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("capture")
    parser.add_argument("--csv", help="write the D lines to this CSV file")
    args = parser.parse_args(argv)
    with open(args.capture, encoding="utf-8", errors="replace") as f:
        boots = parse(f)
    if not boots:
        print("no R frames in capture (was the sketch built with SEA_STATE_SERIAL_DIAG=1?)", file=sys.stderr)
        return 1
    for i, b in enumerate(boots):
        report_boot(i, b, sys.stdout)
    if args.csv:
        write_csv(boots, args.csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
