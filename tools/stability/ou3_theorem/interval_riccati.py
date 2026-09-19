"""Verified interval-Riccati certificate schema for the OU-III A21 proof.

This module is intentionally fail-closed.  A floating-point Riccati trajectory
may propose midpoint/radius boxes, but only outward residual inclusion and
verified innovation inverses can promote a recurring covariance enclosure.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Iterable

from tools.stability.ou3_theorem.tail_stability import (
    covariance_floor_to_rho,
    riccati_box_inclusion,
    verified_interval_innovation_inverse,
)


@dataclass(frozen=True)
class InnovationCertificate:
    name: str
    midpoint_min_eigenvalue: float
    spectral_radius_bound: float

    def verify(self) -> dict:
        out=verified_interval_innovation_inverse(
            midpoint_min_eigenvalue=self.midpoint_min_eigenvalue,
            spectral_radius_bound=self.spectral_radius_bound)
        return {"name":self.name,**out}


@dataclass(frozen=True)
class RiccatiImageCertificate:
    proposed_lower: float
    proposed_upper: float
    image_lower: float
    image_upper: float
    outward_rounding_slack: float

    def verify(self) -> dict:
        return riccati_box_inclusion(
            proposed_lower=self.proposed_lower,
            proposed_upper=self.proposed_upper,
            image_lower=self.image_lower,
            image_upper=self.image_upper,
            outward_rounding_slack=self.outward_rounding_slack)


@dataclass(frozen=True)
class IntervalRiccatiCertificate:
    qualification: str
    fixed_coordinate_mu: float
    innovations: tuple[InnovationCertificate,...]
    recurring_image: RiccatiImageCertificate
    covariance_hard_events_included: bool
    full_shipping_schedule_included: bool
    float32_rounding_included: bool

    def verify(self) -> dict:
        if self.qualification != "OU3_A21_INTERVAL_RICCATI_V1":
            raise ValueError("wrong interval-Riccati qualification")
        if not math.isfinite(self.fixed_coordinate_mu) or self.fixed_coordinate_mu <= 0:
            raise ValueError("positive fixed-coordinate information floor required")
        flags=(self.covariance_hard_events_included,
               self.full_shipping_schedule_included,
               self.float32_rounding_included)
        if any(type(x) is not bool for x in flags):
            raise ValueError("literal certificate flags required")
        inv=[x.verify() for x in self.innovations]
        image=self.recurring_image.verify()
        complete=bool(all(x["verified"] for x in inv) and image["verified"] and all(flags))
        rho=None
        if complete:
            rho=covariance_floor_to_rho(
                fixed_coordinate_mu=self.fixed_coordinate_mu,
                covariance_floor=self.recurring_image.proposed_lower)
        return {
            "innovation_inverses":inv,
            "recurring_image":image,
            "covariance_hard_events_included":self.covariance_hard_events_included,
            "full_shipping_schedule_included":self.full_shipping_schedule_included,
            "float32_rounding_included":self.float32_rounding_included,
            "certificate_complete":complete,
            "rho_chain":rho,
        }


def load_certificate(path: str | Path) -> IntervalRiccatiCertificate:
    raw=json.loads(Path(path).read_text())
    innovations=tuple(InnovationCertificate(**x) for x in raw["innovations"])
    image=RiccatiImageCertificate(**raw["recurring_image"])
    return IntervalRiccatiCertificate(
        qualification=raw["qualification"],
        fixed_coordinate_mu=float(raw["fixed_coordinate_mu"]),
        innovations=innovations,
        recurring_image=image,
        covariance_hard_events_included=raw["covariance_hard_events_included"],
        full_shipping_schedule_included=raw["full_shipping_schedule_included"],
        float32_rounding_included=raw["float32_rounding_included"],
    )
