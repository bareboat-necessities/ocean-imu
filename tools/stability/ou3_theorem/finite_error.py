"""Finite-error service-superword point audit, not a theorem certificate.

The supplied numbers and flags describe one diagnostic. They do not prove
source coverage, storage coercivity or the assertions represented by flags.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Sequence


@dataclass(frozen=True)
class ServiceSuperwordWitness:
    history_id: str
    V_start: float
    V_end: float
    disturbance_energy: float
    rho_candidate: float
    disturbance_gain: float
    prefix_values: tuple[float, ...]
    prefix_bounds: tuple[float, ...] | None
    magnetic_information_min_eigenvalue: float
    magnetic_information_floor: float
    complete_shipping_map: bool
    inherited_state: bool
    same_history_motion: bool
    same_history_bias: bool
    applied_magnetic_information: bool
    finite_error_map: bool
    arithmetic_enclosed: bool

    def __post_init__(self) -> None:
        if not isinstance(self.history_id, str) or not self.history_id.strip():
            raise ValueError("nonempty history_id required")
        scalars = (
            self.V_start, self.V_end, self.disturbance_energy,
            self.rho_candidate, self.disturbance_gain,
            self.magnetic_information_min_eigenvalue,
            self.magnetic_information_floor,
        )
        if not all(math.isfinite(v) for v in scalars):
            raise ValueError("finite witness scalars required")
        if min(self.V_start, self.V_end, self.disturbance_energy,
               self.disturbance_gain) < 0.0 or self.magnetic_information_floor <= 0.0:
            raise ValueError("nonnegative storage/energy/gain and positive information floor required")
        if not self.prefix_values:
            raise ValueError("prefix data must include the start and end values")
        if not all(math.isfinite(v) and v >= 0.0 for v in self.prefix_values):
            raise ValueError("finite nonnegative prefix values required")
        if self.prefix_values[0] != self.V_start or self.prefix_values[-1] != self.V_end:
            raise ValueError("prefix endpoints must match the storage endpoints")
        if self.prefix_bounds is not None:
            if len(self.prefix_values) != len(self.prefix_bounds):
                raise ValueError("prefix values and bounds must have equal length")
            if not all(math.isfinite(v) and v >= 0.0 for v in self.prefix_bounds):
                raise ValueError("finite nonnegative prefix bounds required")
        for value in self.structural_flags():
            if type(value) is not bool:
                raise ValueError("structural assertions must be literal booleans")

    def structural_flags(self) -> tuple[bool, ...]:
        return (
            self.complete_shipping_map, self.inherited_state,
            self.same_history_motion, self.same_history_bias,
            self.applied_magnetic_information, self.finite_error_map,
            self.arithmetic_enclosed,
        )


def audit_service_superword(w: ServiceSuperwordWitness) -> dict:
    """Check supplied point values without promoting asserted premises.

    The target is V_end <= rho*V_start + gain*energy with rho<1, plus
    every-prefix retention. Even a consistent point does not provide a
    source-uniform finite-error certificate or validate the structural flags.
    """
    strict_rho = 0.0 <= w.rho_candidate < 1.0
    rhs = w.rho_candidate*w.V_start + w.disturbance_gain*w.disturbance_energy
    finite_rhs = math.isfinite(rhs)
    dissipative = finite_rhs and w.V_end <= rhs
    prefixes = w.prefix_bounds is not None and all(
        v <= b for v, b in zip(w.prefix_values, w.prefix_bounds)
    )
    information = w.magnetic_information_min_eigenvalue >= w.magnetic_information_floor
    structural = all(w.structural_flags())
    consistent = bool(strict_rho and dissipative and prefixes and information and structural)
    return {
        "history_id": w.history_id,
        "strict_rho": strict_rho,
        "dissipativity": dissipative,
        "dissipativity_rhs": rhs if finite_rhs else None,
        "arithmetic_finite": finite_rhs,
        "prefix_retention": prefixes,
        "magnetic_information": information,
        "structural_premises": structural,
        "structural_premises_verified_here": False,
        "point_consistency_pass": consistent,
        "certificate_complete": False,
        "point_witness_promoted_to_source_uniform_theorem": False,
    }


def incomplete_diagnostic(history_id: str, storage_values: Sequence[float],
                          information_floor: float) -> ServiceSuperwordWitness:
    """Represent missing prefix bounds explicitly, rather than invent infinity."""
    values = tuple(float(v) for v in storage_values)
    if not values:
        raise ValueError("at least one storage value required")
    return ServiceSuperwordWitness(
        history_id=history_id, V_start=values[0], V_end=values[-1],
        disturbance_energy=0.0,
        rho_candidate=values[-1]/values[0] if values[0] > 0.0 else 1.0,
        disturbance_gain=0.0, prefix_values=values, prefix_bounds=None,
        magnetic_information_min_eigenvalue=0.0,
        magnetic_information_floor=information_floor,
        complete_shipping_map=False, inherited_state=False,
        same_history_motion=False, same_history_bias=False,
        applied_magnetic_information=False, finite_error_map=False,
        arithmetic_enclosed=False,
    )


@dataclass(frozen=True)
class HeldBiasSuperwordBlock:
    """Literal shipping facts about the held accelerometer bias on one superword.

    While H18 holds the accelerometer bias, the shipping estimator applies no
    bias mean dynamics, freezes the bias rows of every gain, and keeps the bias
    cross-covariances at the zero they were set to when the hold was taken. The
    three fields record exactly that, measured on the execution rather than
    assumed: the held-bias block of a same-history incremental superword map, the
    largest surviving bias cross-covariance at either endpoint, and the change
    in the held-bias covariance block across the superword.
    """
    history_id: str
    map_block: tuple[tuple[float, float, float], ...]
    cross_covariance_max: float
    covariance_block_change: float

    def __post_init__(self) -> None:
        if not isinstance(self.history_id, str) or not self.history_id.strip():
            raise ValueError("nonempty history_id required")
        if len(self.map_block) != 3 or any(len(row) != 3 for row in self.map_block):
            raise ValueError("the held accelerometer-bias block is 3x3")
        values = [v for row in self.map_block for v in row]
        values += [self.cross_covariance_max, self.covariance_block_change]
        if not all(math.isfinite(v) for v in values):
            raise ValueError("finite measured block entries required")
        if min(self.cross_covariance_max, self.covariance_block_change) < 0.0:
            raise ValueError("measured magnitudes are nonnegative")


def held_bias_non_contraction(block: HeldBiasSuperwordBlock,
                              tolerance: float = 0.0) -> dict:
    """Held-bias obstruction for the incremental covariance-storage map.

    On a feasible held segment with no projection-changing or frame/relock
    event, the same-history differential is [[A,B],[0,I]]. Physical bias
    increments cancel between those executions, not from absolute error:
    e_b_plus=e_b+w still holds. This is not a global linear finite-error map.

    With decoupled frozen P_bb, a direction (0,u) has endpoint energy
    (B u)^T P_oo,end^-1 (B u)+u^T P_bb^-1 u, at least its initial energy.
    The H18 complement therefore carries bounded held bias as an input in the
    single theorem path. Projection, hard events and release remain separate.

    A tolerance may diagnose proximity, but only exact supplied premises
    exclude contraction by this argument. Neither case certifies the source.
    Missing premises decide nothing about availability of other contractions.
    """
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("finite nonnegative diagnostic tolerance required")
    identity_defect = max(
        abs(block.map_block[i][j] - (1.0 if i == j else 0.0))
        for i in range(3) for j in range(3)
    )
    reproduces = identity_defect == 0.0
    decoupled = block.cross_covariance_max == 0.0
    frozen = block.covariance_block_change == 0.0
    obstructed = bool(reproduces and decoupled and frozen)
    return {
        "history_id": block.history_id,
        "held_bias_identity_defect": identity_defect,
        "held_bias_reproduces_itself": reproduces,
        "bias_covariance_decoupled": decoupled,
        "bias_covariance_frozen": frozen,
        "non_contraction_obstruction": obstructed,
        "within_diagnostic_tolerance": max(identity_defect, block.cross_covariance_max,
                                           block.covariance_block_change) <= tolerance,
        # Failing to establish this one obstruction says nothing about whether
        # rho < 1 is available: some other obstruction, or none, may apply. The
        # only two outcomes here are "excluded by this argument" and "this
        # argument does not decide it".
        "full_state_rho_status": "excluded" if obstructed else "undecided_here",
        "obstruction_is_an_instability_claim": False,
        "certificate_complete": False,
    }
