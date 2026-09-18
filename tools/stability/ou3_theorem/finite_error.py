"""Finite-error H18 service-superword proof target.

This module states the complete inequality that a future rigorous certificate
must close.  It deliberately refuses to promote point diagnostics to theorem
status.  All quantities belong to one persistent execution.
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
    prefix_bounds: tuple[float, ...]
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
        if not self.history_id:
            raise ValueError("nonempty history_id required")
        scalars=(
            self.V_start,self.V_end,self.disturbance_energy,self.rho_candidate,
            self.disturbance_gain,self.magnetic_information_min_eigenvalue,
            self.magnetic_information_floor,
        )
        if not all(math.isfinite(v) for v in scalars):
            raise ValueError("finite witness scalars required")
        if self.V_start < 0.0 or self.V_end < 0.0 or self.disturbance_energy < 0.0:
            raise ValueError("storage and disturbance energy must be nonnegative")
        if self.disturbance_gain < 0.0 or self.magnetic_information_floor <= 0.0:
            raise ValueError("nonnegative gain and positive information floor required")
        if len(self.prefix_values)!=len(self.prefix_bounds):
            raise ValueError("prefix values and bounds must have equal length")
        if not all(math.isfinite(v) and v>=0.0 for v in self.prefix_values+self.prefix_bounds):
            raise ValueError("finite nonnegative prefix data required")


def audit_service_superword(w: ServiceSuperwordWitness) -> dict:
    """Fail-closed audit of the intended H18 finite-error inequality.

    A complete certificate requires
      V_end <= rho V_start + c_d ||d||^2, rho < 1,
    every-prefix retention, the required applied magnetic information, and the
    literal inherited shipping map with arithmetic enclosure.  This function
    checks supplied evidence; it does not manufacture any missing premise.
    """
    strict_rho=0.0 <= w.rho_candidate < 1.0
    rhs=w.rho_candidate*w.V_start+w.disturbance_gain*w.disturbance_energy
    dissipative=w.V_end <= rhs
    prefixes=all(v<=b for v,b in zip(w.prefix_values,w.prefix_bounds))
    information=w.magnetic_information_min_eigenvalue>=w.magnetic_information_floor
    structural=all((
        w.complete_shipping_map,
        w.inherited_state,
        w.same_history_motion,
        w.same_history_bias,
        w.applied_magnetic_information,
        w.finite_error_map,
        w.arithmetic_enclosed,
    ))
    complete=bool(strict_rho and dissipative and prefixes and information and structural)
    return {
        "history_id":w.history_id,
        "strict_rho":strict_rho,
        "dissipativity":dissipative,
        "dissipativity_rhs":rhs,
        "prefix_retention":prefixes,
        "magnetic_information":information,
        "structural_premises":structural,
        "certificate_complete":complete,
        "point_witness_promoted_to_source_uniform_theorem":False,
    }


def incomplete_diagnostic(history_id: str, storage_values: Sequence[float],
                          information_floor: float) -> ServiceSuperwordWitness:
    """Represent a numerical diagnostic without falsely closing the theorem."""
    values=tuple(float(v) for v in storage_values)
    if not values:
        raise ValueError("at least one storage value required")
    return ServiceSuperwordWitness(
        history_id=history_id,
        V_start=values[0],
        V_end=values[-1],
        disturbance_energy=0.0,
        rho_candidate=(values[-1]/values[0] if values[0]>0.0 else 1.0),
        disturbance_gain=0.0,
        prefix_values=values,
        prefix_bounds=tuple(math.inf for _ in values),
        magnetic_information_min_eigenvalue=0.0,
        magnetic_information_floor=information_floor,
        complete_shipping_map=False,
        inherited_state=False,
        same_history_motion=False,
        same_history_bias=False,
        applied_magnetic_information=False,
        finite_error_map=False,
        arithmetic_enclosed=False,
    )
