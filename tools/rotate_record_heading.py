#!/usr/bin/env python3
"""Rotate the vessel heading of a wave record about the world vertical.

The sea is left untouched: world displacement, velocity and acceleration
columns, and therefore the wave field the platform sits in, are identical in the
output.  Only the platform's orientation changes, from ``R_old(t)`` to

    R_new(t) = Rz(psi) R_old(t).

Body specific force is re-resolved in the new body frame, because the specific
force is a world-frame field at an unchanged location.  Body angular rate is
*not* transformed: a constant world-frame rotation of the attitude gives
``omega_world_new = Rz(psi) omega_world_old``, so
``omega_body = R_new^T omega_world_new`` is unchanged.  The magnetometer is
synthesised by the simulator from the reference attitude, so it follows
automatically.

This is the controlled experiment for the wave travel-sense gauge question: the
same sea seen from a different heading must leave heave, 3D displacement and the
propagation axis unchanged, and must leave the *directed* propagation angle
unchanged once the heading is removed, while any output defined against the
estimator's own axis representative will invert.

Usage:
    python3 tools/rotate_record_heading.py IN.csv OUT.csv 180

Name the output so the simulator still parses the same generator azimuth (the
``_A<deg>_`` field), because the sea itself has not moved; vary another field
(for example the trailing ``_P<phase>_``) to keep the two records distinct.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

QUATERNION_COLUMNS = ("q_wb_zu_w", "q_wb_zu_x", "q_wb_zu_y", "q_wb_zu_z")
SPECIFIC_FORCE_COLUMNS = ("acc_bx", "acc_by", "acc_bz")


def quaternion_multiply(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    aw, ax, ay, az = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    bw, bx, by, bz = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack(
        (
            aw * bw - ax * bx - ay * by - az * bz,
            aw * bx + ax * bw + ay * bz - az * by,
            aw * by - ax * bz + ay * bw + az * bx,
            aw * bz + ax * by - ay * bx + az * bw,
        ),
        axis=-1,
    )


def quaternion_conjugate(q: np.ndarray) -> np.ndarray:
    out = q.copy()
    out[..., 1:] *= -1.0
    return out


def quaternion_rotate(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Apply the rotation represented by ``q`` to ``v``."""
    q_vector = q[..., 1:]
    t = 2.0 * np.cross(q_vector, v)
    return v + q[..., :1] * t + np.cross(q_vector, t)


def rotate_record(source: Path, destination: Path, heading_offset_deg: float) -> None:
    columns = source.read_text(encoding="utf-8").split("\n", 1)[0].strip().split(",")
    index = {name: position for position, name in enumerate(columns)}
    missing = [
        name
        for name in (*QUATERNION_COLUMNS, *SPECIFIC_FORCE_COLUMNS, "yaw_deg")
        if name not in index
    ]
    if missing:
        raise ValueError(f"record is missing columns: {', '.join(missing)}")

    data = np.loadtxt(source, delimiter=",", skiprows=1)

    quaternion_positions = [index[name] for name in QUATERNION_COLUMNS]
    # The record stores world -> body: rotating a body vector by its conjugate
    # reproduces the world-frame acceleration columns.
    q_world_to_body = data[:, quaternion_positions]
    q_world_to_body /= np.linalg.norm(q_world_to_body, axis=1, keepdims=True)
    q_body_to_world = quaternion_conjugate(q_world_to_body)

    half = np.deg2rad(heading_offset_deg) * 0.5
    q_heading = np.tile(
        np.array([np.cos(half), 0.0, 0.0, np.sin(half)]), (len(data), 1)
    )

    q_body_to_world_new = quaternion_multiply(q_heading, q_body_to_world)
    q_world_to_body_new = quaternion_conjugate(q_body_to_world_new)

    force_positions = [index[name] for name in SPECIFIC_FORCE_COLUMNS]
    specific_force_world = quaternion_rotate(
        q_body_to_world, data[:, force_positions]
    )
    data[:, force_positions] = quaternion_rotate(
        q_world_to_body_new, specific_force_world
    )

    data[:, quaternion_positions] = q_world_to_body_new
    yaw = data[:, index["yaw_deg"]] + heading_offset_deg
    data[:, index["yaw_deg"]] = (yaw + 180.0) % 360.0 - 180.0

    with destination.open("w", encoding="utf-8") as stream:
        stream.write(",".join(columns) + "\n")
        np.savetxt(stream, data, delimiter=",", fmt="%.9g")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("heading_offset_deg", type=float)
    args = parser.parse_args()
    rotate_record(args.source, args.destination, args.heading_offset_deg)


if __name__ == "__main__":
    main()
