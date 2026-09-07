#!/usr/bin/env python3
"""Project host-observer physical payload H back to exact shipping structure.

The single shipping observer reconstructs H from P^{-1} PCt after the executed
Joseph update because production code intentionally does not expose a stored H.
In binary32 this introduces tiny symmetric/skew and structurally-zero-column
noise even when P*H^T reproduces PCt accurately.  That numerical reconstruction
must not become fake physical geometry in the non-promoting finite-map probe.

This tool changes only the diagnostic payload representation:

* S=0 H is the exact S selector;
* vector H keeps only the closest skew attitude block;
* accelerometer H keeps the closest skew attitude block, projects the recovered
  a_w block to the nearest SO(3) matrix, and restores the exact active-b_a I
  block in A21;
* measurement R is symmetrized and projected to its configured diagonal form.

P, Q, shipping linear maps, event order, actual R_S S events, schedule values,
and every source timestamp are copied bit-for-bit.  The downstream zero-state
event/whole-word parity checks remain the authority that this host-only
projection is consistent with the executed shipping word.  This is diagnostic
plumbing only and cannot promote P4.
"""
from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path

import numpy as np

MAGIC = b"OU3PHY1\0"
HEADER = struct.Struct("<8sIIIdd")
PREFIX = struct.Struct("<IIdfffff")
NX = 21
EV_PRED = 1
EV_FLOOR = 2
EV_S = 3
EV_ACC = 4
EV_VEC = 5
FLOATS_PER_RECORD = 6 + 4 * NX * NX + 3 * NX + 9 + 9
BYTES_PER_RECORD = PREFIX.size + 4 * FLOATS_PER_RECORD


def _nearest_so3(R: np.ndarray) -> np.ndarray:
    U, _, Vt = np.linalg.svd(R.astype(np.float64), full_matrices=False)
    Q = U @ Vt
    if np.linalg.det(Q) < 0.0:
        U[:, -1] *= -1.0
        Q = U @ Vt
    return Q.astype(np.float32)


def _structural_H(H: np.ndarray, event_type: int, mode_dim: int) -> np.ndarray:
    out = np.zeros((3, NX), dtype=np.float32)
    if event_type == EV_S:
        out[:, 12:15] = np.eye(3, dtype=np.float32)
        return out
    if event_type not in (EV_ACC, EV_VEC):
        return H.copy()

    att = 0.5 * (H[:, 0:3].astype(np.float64) - H[:, 0:3].T.astype(np.float64))
    out[:, 0:3] = att.astype(np.float32)
    if event_type == EV_ACC:
        out[:, 15:18] = _nearest_so3(H[:, 15:18])
        if mode_dim == 21:
            out[:, 18:21] = np.eye(3, dtype=np.float32)
    return out


def _structural_R(R: np.ndarray, event_type: int) -> np.ndarray:
    if event_type == EV_S:
        # Actual applied anisotropic R_S is already the direct source value.
        return R.copy()
    if event_type in (EV_ACC, EV_VEC):
        S = 0.5 * (R.astype(np.float64) + R.T.astype(np.float64))
        return np.diag(np.diag(S)).astype(np.float32)
    return R.copy()


def _take(vals: np.ndarray, start: int, count: int) -> tuple[np.ndarray, int]:
    """Return one contiguous record slice and the next parser offset."""
    end = start + count
    return vals[start:end], end


def canonicalize(src: Path, dst: Path) -> dict:
    raw = src.read_bytes()
    if len(raw) < HEADER.size:
        raise ValueError("physical payload is truncated")
    magic, version, nx, mode_dim, t0, t1 = HEADER.unpack_from(raw, 0)
    if magic != MAGIC or version != 1 or nx != NX or mode_dim not in (18, 21):
        raise ValueError("unsupported physical payload header")
    if (len(raw) - HEADER.size) % BYTES_PER_RECORD:
        raise ValueError("physical payload record size mismatch")

    out = bytearray(raw[:HEADER.size])
    offset = HEADER.size
    event_counts: dict[str, int] = {}
    max_H_relative_projection = 0.0
    max_H_supported_relative_projection = 0.0
    max_R_offdiag_removed = 0.0
    max_so3_orthogonality_before = 0.0
    max_so3_orthogonality_after = 0.0

    names = {EV_PRED: "prediction", EV_FLOOR: "aw_floor", EV_S: "S_zero", EV_ACC: "accelerometer", EV_VEC: "vector"}

    while offset < len(raw):
        prefix = raw[offset:offset + PREFIX.size]
        event_type, record_dim, *_ = PREFIX.unpack(prefix)
        offset += PREFIX.size
        if record_dim != mode_dim or event_type not in names:
            raise ValueError("physical payload event header mismatch")
        vals = np.frombuffer(raw, dtype="<f4", count=FLOATS_PER_RECORD, offset=offset).copy()
        offset += 4 * FLOATS_PER_RECORD

        k = 0
        _, k = _take(vals, k, 3)  # omega
        _, k = _take(vals, k, 3)  # dtheta
        _, k = _take(vals, k, NX * NX)  # Pbefore
        _, k = _take(vals, k, NX * NX)  # Pafter
        _, k = _take(vals, k, NX * NX)  # linear map
        _, k = _take(vals, k, NX * NX)  # Q
        H_start = k
        Hflat, k = _take(vals, k, 3 * NX)
        H = Hflat.reshape(3, NX)
        R_start = k
        Rflat, k = _take(vals, k, 9)
        R = Rflat.reshape(3, 3)
        _, k = _take(vals, k, 9)  # aux
        if k != FLOATS_PER_RECORD:
            raise RuntimeError("physical payload parser drifted")

        Hnew = _structural_H(H, event_type, mode_dim)
        Rnew = _structural_R(R, event_type)
        if event_type in (EV_S, EV_ACC, EV_VEC):
            denom = max(1.0, float(np.linalg.norm(H[:, :mode_dim])))
            rel = float(np.linalg.norm(Hnew[:, :mode_dim].astype(float) - H[:, :mode_dim].astype(float)) / denom)
            max_H_relative_projection = max(max_H_relative_projection, rel)

            if event_type == EV_ACC:
                cols = list(range(3)) + list(range(15, 18)) + (list(range(18, 21)) if mode_dim == 21 else [])
            elif event_type == EV_VEC:
                cols = list(range(3))
            else:
                cols = list(range(12, 15))
            supported = float(np.linalg.norm(Hnew[:, cols].astype(float) - H[:, cols].astype(float)) /
                              max(1.0, float(np.linalg.norm(H[:, cols]))))
            max_H_supported_relative_projection = max(max_H_supported_relative_projection, supported)

            if event_type == EV_ACC:
                Rwb = H[:, 15:18].astype(float)
                Qwb = Hnew[:, 15:18].astype(float)
                max_so3_orthogonality_before = max(
                    max_so3_orthogonality_before,
                    float(np.linalg.norm(Rwb.T @ Rwb - np.eye(3))),
                )
                max_so3_orthogonality_after = max(
                    max_so3_orthogonality_after,
                    float(np.linalg.norm(Qwb.T @ Qwb - np.eye(3))),
                )

            offdiag = R.astype(float) - np.diag(np.diag(R.astype(float)))
            max_R_offdiag_removed = max(max_R_offdiag_removed, float(np.max(np.abs(offdiag))))

        vals[H_start:H_start + 3 * NX] = Hnew.reshape(-1)
        vals[R_start:R_start + 9] = Rnew.reshape(-1)
        out.extend(prefix)
        out.extend(vals.astype("<f4", copy=False).tobytes())
        event_counts[names[event_type]] = event_counts.get(names[event_type], 0) + 1

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_bytes(bytes(out))
    if len(out) != len(raw):
        raise RuntimeError("canonical payload changed binary size")

    return {
        "qualification": "NON_PROMOTING_PHYSICAL_PAYLOAD_STRUCTURAL_H_PROJECTION",
        "input": str(src),
        "output": str(dst),
        "mode_dim": int(mode_dim),
        "t0": float(t0),
        "t1": float(t1),
        "event_counts": event_counts,
        "max_H_relative_projection_all_columns": max_H_relative_projection,
        "max_H_relative_projection_supported_columns": max_H_supported_relative_projection,
        "max_R_offdiag_removed": max_R_offdiag_removed,
        "max_accelerometer_Rwb_orthogonality_error_before": max_so3_orthogonality_before,
        "max_accelerometer_Rwb_orthogonality_error_after": max_so3_orthogonality_after,
        "P_Q_linear_map_schedule_timestamps_copied_bit_for_bit": True,
        "actual_RS_S_events_copied_bit_for_bit": True,
        "finite_source_or_theorem_word_generated_here": False,
        "P4_promoted": False,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    args = ap.parse_args()
    d = canonicalize(args.input, args.output)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(d, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
