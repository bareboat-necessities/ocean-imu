"""Exact paired covariance/probe dissipation for one Kalman correction.

For S=HPH'+R, K=PH'S^-1, A=I-KH and the Joseph covariance P+,
P+=AP=PA'. Hence A'(P+)^-1 A=P^-1-H'S^-1H, and propagation
Phi+=A Phi subtracts (H Phi)'S^-1(H Phi) from paired storage.
Covariance and homogeneous columns must be carried together; their storage
is not invariant across informative accelerometer, S or magnetic corrections.
"""
from __future__ import annotations
import numpy as np


def correction(P, Phi, H, R):
    P = np.asarray(P, float)
    Phi = np.asarray(Phi, float)
    H = np.asarray(H, float)
    R = np.asarray(R, float)
    S = H @ P @ H.T + R
    K = np.linalg.solve(S, H @ P).T
    A = np.eye(P.shape[0]) - K @ H
    Pp = A @ P @ A.T + K @ R @ K.T
    return Pp, A @ Phi


def storage(P, Phi):
    return Phi.T @ np.linalg.solve(P, Phi)


def measurement_loss(Phi, H, S):
    Y = H @ Phi
    return Y.T @ np.linalg.solve(S, Y)


def prior_to_posterior_drop(P, Phi, H, R):
    S = H @ P @ H.T + R
    Pp, Xp = correction(P, Phi, H, R)
    drop = storage(P, Phi) - storage(Pp, Xp)
    service = measurement_loss(Phi, H, S)
    return {"storage_drop": drop, "service": service,
            "storage_defect": float(np.linalg.norm(drop - service))}


def identity_defect(P, Phi, H, R):
    return prior_to_posterior_drop(P, Phi, H, R)["storage_defect"]


def certificate():
    rng = np.random.default_rng(7)
    A = rng.normal(size=(6, 6))
    P = A @ A.T + .3 * np.eye(6)
    X = rng.normal(size=(6, 2))
    H = rng.normal(size=(3, 6))
    R = .4 * np.eye(3)
    result = prior_to_posterior_drop(P, X, H, R)
    return {
        "qualification": "OU3_PAIRED_COVARIANCE_PROBE_DISSIPATION_V2",
        "identity_defect": result["storage_defect"],
        "measurement_loss_norm": float(np.linalg.norm(result["service"])),
        "identity": "Phi_plus^T P_plus^-1 Phi_plus = Phi_minus^T P_minus^-1 Phi_minus - (H Phi_minus)^T (H P_minus H^T + R)^-1 (H Phi_minus)",
        "consequence": "each informative accel/S/mag correction consumes exactly its own measurement action in paired storage; nonmagnetic losses cannot be omitted from a magnetic-service comparison",
        "warning": "total paired-storage dissipation alone does not lower-bound the magnetic portion; prediction/process noise and every nonmagnetic correction must also be retained",
        "shipping_service_proved": False,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
