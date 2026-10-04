"""Literal planar operation stream; finite diagnostics, never all-time admission.

The exporter observes the same full shipping execution as the MOVING witness.
Only covariance arrays in the file are losslessly parity-compressed. The mean,
frontend, tuner, scheduler, gains, AW maintenance and reset run unmodified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import tempfile
from collections import Counter
import numpy as np

from .moving_compatibility_diagnostic import HEADER, REPO, instrument as base_instrument
from .planar_parity import EVEN, ODD

PROBE = Path(__file__).with_name("planar_service_probe.cpp")
MAGIC = b"OU3EVT1\0"
# prediction; acc, mag, S; pre-reset; post-AW-sync; sample; post-reset
SIZES = {1: 912, 2: 435, 3: 435, 4: 435, 5: 228, 6: 235, 7: 293, 8: 225}


def instrument(text: str) -> str:
    try:
        s = base_instrument(text)
    except AssertionError as exc:
        raise ValueError("base shipping tap anchors changed") from exc
    anchor = "    // Always propagate attitude(+gyro-bias) covariance"
    if s.count(anchor) != 1:
        raise ValueError("pre-prediction tap anchor changed")
    s = s.replace(anchor, "    if (moving_recording) planar_before_prediction(Pext.data());\n" + anchor)
    old = "Qba_diag.data(),Pext.data());"
    if s.count(old) != 2:
        raise ValueError("literal prediction tap changed")
    s = s.replace(old, "Qba_diag.data(),Pext.data(),R_S.data(),pseudo_update_elapsed_s_,pseudo_update_period_s_,Ts);")
    for kind, noise in (("acc", "Racc"), ("mag", "Rmag"), ("S", "R_S")):
        old = f'moving_correction("{kind}",H.data(),S_mat.data(),K.data());'
        if s.count(old) != 1:
            raise ValueError(f"{kind} correction tap changed")
        s = s.replace(old, f'moving_correction("{kind}",H.data(),S_mat.data(),K.data(),Pext.data(),{noise}.data(),PCt.data(),r.data());')
    old = "        apply_pending_aw_covariance_inflation_();"
    if s.count(old) != 1:
        raise ValueError("AW inflation tap changed")
    s = s.replace(old, """        const bool planar_sync_pending = aw_covariance_floor_pending_;
        apply_pending_aw_covariance_inflation_();
        if (moving_recording) planar_after_sync(Pext.data(),aw_covariance_floor_target_.data(),planar_sync_pending);""")
    old = "moving_reset(dtheta_injected.data());"
    if s.count(old) != 1:
        raise ValueError("reset tap changed")
    s = s.replace(old, "moving_reset(dtheta_injected.data(),Pext.data());")
    old = "        ocean_imu::kalman::ou_detail::apply_left_error_reset<T, NX>(Pext, dtheta_injected);"
    if s.count(old) != 1:
        raise ValueError("post-reset tap changed")
    s = s.replace(old, old + "\n        if (moving_recording) planar_after_reset(Pext.data());")
    return s


def expand(values: np.ndarray) -> np.ndarray:
    """Restore all 21 coordinates without dropping covariance cross terms."""
    if values.size != 225:
        raise ValueError("12+9 covariance block payload required")
    out = np.zeros((21, 21), dtype=np.float64)
    out[np.ix_(EVEN, EVEN)] = values[:144].reshape(12, 12)
    out[np.ix_(ODD, ODD)] = values[144:].reshape(9, 9)
    return out


def records(path: Path):
    """Stream strict, finite, chronological records; reject truncated evidence."""
    with Path(path).open("rb") as f:
        if f.read(8) != MAGIC:
            raise ValueError("bad operation-stream magic")
        previous = 0
        while True:
            head = f.read(12)
            if not head:
                break
            if len(head) != 12:
                raise ValueError("truncated operation header")
            kind, sample, count = struct.unpack("<III", head)
            if kind not in SIZES or count != SIZES[kind]:
                raise ValueError(f"unknown operation layout: {kind}, {count}")
            if sample < previous:
                raise ValueError("nonchronological operation stream")
            previous = sample
            payload = f.read(count * 4)
            if len(payload) != count * 4:
                raise ValueError("truncated operation payload")
            a = np.frombuffer(payload, dtype="<f4").astype(np.float64)
            if not np.isfinite(a).all():
                raise ValueError("nonfinite exported value")
            yield kind, sample, a


def run_native(eigen: Path, duration: float, tail: float, output: Path,
               cxx: str = "g++", verify_untapped: bool = True) -> dict:
    if duration < 220 or duration > 3600 or tail < 1 or tail > duration - 180:
        raise ValueError("require 220<=duration<=3600, 1<=tail<=duration-180")
    eigen, output = Path(eigen), Path(output).resolve()
    if not (eigen / "Eigen/Dense").is_file():
        raise ValueError("Eigen/Dense missing")
    tapped = instrument(HEADER.read_text())
    compiler = subprocess.run([cxx, "--version"], check=True, capture_output=True,
                              text=True, timeout=10).stdout.splitlines()[0]
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ou3-planar-service-") as d:
        work = Path(d)
        header = work / "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
        header.parent.mkdir(parents=True)
        header.write_text(tapped)
        command = [cxx, "-std=c++20", "-O1", "-DEIGEN_UNROLLING_LIMIT=0", "-DEIGEN_NON_ARDUINO",
                   "-I" + str(work), "-I" + str(REPO / "src"), "-I" + str(eigen), str(PROBE),
                   "-o", str(work / "probe")]
        p = subprocess.run(command, capture_output=True, text=True, timeout=150)
        if p.returncode:
            raise RuntimeError("compile failed:\n" + p.stderr)
        p = subprocess.run([str(work / "probe"), str(duration), str(tail), str(output)],
                           check=True, capture_output=True, text=True, timeout=max(240, int(duration)))
        native = json.loads(p.stdout)
        control_verified = False
        if verify_untapped:
            # Same probe and inputs, but the original shipping header resolves
            # from src. No observation callback is inserted into the estimator.
            control_command = [x for x in command if x != "-I" + str(work)]
            control_command[-1] = str(work / "control")
            p = subprocess.run(control_command, capture_output=True, text=True, timeout=150)
            if p.returncode:
                raise RuntimeError("untapped control compile failed:\n" + p.stderr)
            control_stream = work / "control.bin"
            p = subprocess.run([str(work / "control"), str(duration), str(tail), str(control_stream)],
                check=True, capture_output=True, text=True, timeout=max(240, int(duration)))
            control = json.loads(p.stdout)
            for key in ("terminal_state", "live_step", "refined_step", "active_step", "accepted_mag"):
                if control[key] != native[key]:
                    raise ValueError("read-only tap changed the shipping execution: " + key)
            from itertools import zip_longest
            original = ((k,a) for kind,k,a in records(control_stream) if kind == 7)
            observed = ((k,a) for kind,k,a in records(output) if kind == 7)
            for left,right in zip_longest(original,observed):
                if left is None or right is None or left[0] != right[0] or not np.array_equal(left[1],right[1]):
                    raise ValueError("read-only tap changed a complete tail sample")
            control_verified = True
    counts = Counter(kind for kind, _, _ in records(output))
    if counts[1] != round(tail * 200) or counts[7] != counts[1]:
        raise ValueError("incomplete sample/prediction coverage")
    if native["parity_max_abs"] != 0:
        raise ValueError("parity compression is not lossless on this replay")
    return {
        "qualification": "OU3_PLANAR_SERVICE_OPERATION_STREAM_V1",
        "result_type": "FINITE DIAGNOSTIC ONLY",
        "native": native,
        "compiler": compiler,
        "untapped_complete_tail_samples_bitwise_equal": control_verified,
        "record_counts": dict(sorted(counts.items())),
        "stream_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "probe_sha256": hashlib.sha256(PROBE.read_bytes()).hexdigest(),
        "driver_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "shipping_header_sha256": hashlib.sha256(HEADER.read_bytes()).hexdigest(),
        "instrumented_header_sha256": hashlib.sha256(tapped.encode()).hexdigest(),
        "all_time_magnetic_service_verified": False,
        "full_shipping_counterexample_admitted": False,
        "shipping_behavior_changed": False,
        "theorem_closed": False,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--eigen", type=Path, default=Path("/usr/include/eigen3"))
    p.add_argument("--cxx", default="g++")
    p.add_argument("--skip-untapped-control", action="store_true")
    p.add_argument("--duration", type=float, default=240.)
    p.add_argument("--tail", type=float, default=40.)
    p.add_argument("--stream", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    result = run_native(a.eigen, a.duration, a.tail, a.stream, a.cxx, not a.skip_untapped_control)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
