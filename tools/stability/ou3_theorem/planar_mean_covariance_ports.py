"""Cell-dependent inequalities; physical acceleration is NOT nominal aw.

No point-orbit innovation floor or magnetic reference norm is a default proof
constant. Bounds must be proved on the same causal cell before application.
"""
from __future__ import annotations


def certificate(*, g=None, nominal_aw_bound=None,
                reference_norm_bound=None, measurement_noise_floor=None):
    values = (g, nominal_aw_bound, reference_norm_bound, measurement_noise_floor)
    supplied = all(x is not None for x in values)
    if any(x is not None and x < 0 for x in values):
        raise ValueError("cell bounds must be nonnegative")
    if measurement_noise_floor is not None and measurement_noise_floor <= 0:
        raise ValueError("a positive noise floor is required")
    # H_acc=[-[R(aw-g)]x,0,...,R,...,I]. Both the skew and R blocks vary.
    # ||dH|| <= sqrt((g+Ahat)^2+1)||dtheta||+||daw||.
    force = None if g is None or nominal_aw_bound is None else g + nominal_aw_bound
    return {
        "qualification": "OU3_PLANAR_MEAN_COVARIANCE_PORT_BOUNDS_V2",
        "result_type": "CONDITIONAL",
        "inequalities_status": "PROVED — analytical",
        "acc_H_attitude_coefficient_squared": None if force is None else str(force**2 + 1),
        "acc_H_lipschitz_aw": "1",
        "mag_H_lipschitz_attitude": None if reference_norm_bound is None else str(reference_norm_bound),
        "innovation_inverse_bound": None if measurement_noise_floor is None else str(1 / measurement_noise_floor),
        "symbolic_acc_H_bound": "||dH_acc|| <= sqrt((g+Ahat)^2+1)||dtheta||+||daw||",
        "symbolic_mag_H_bound": "||dH_mag|| <= Bref||dtheta||+||dBref||",
        "symbolic_inverse_bound": "P>=0 and R>=rI>0 imply ||S^-1||<=1/r",
        "cell_bounds_supplied": supplied,
        "cell_bounds_certified": False,
        "physical_acceleration_substituted_for_nominal_aw": False,
        "carried_point_innovation_floor_used": False,
        "same_history_required": True,
        "joint_cell_forward_invariant": False,
        "all_time_magnetic_service_verified": False,
        "theorem_closed": False,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
