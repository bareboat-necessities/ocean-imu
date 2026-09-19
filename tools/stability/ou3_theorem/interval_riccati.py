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


def outward_dot_interval(a_mid: tuple[float,...], a_rad: tuple[float,...],
                         b_mid: tuple[float,...], b_rad: tuple[float,...]) -> tuple[float,float]:
    """Midpoint-radius dot product with one-ulp outward expansion.

    This primitive is deliberately elementary so the certificate does not
    depend on host BLAS rounding mode. The analytic radius includes every
    midpoint/radius product; nextafter then covers the final binary64
    evaluation. A float32 shipping kernel adds its separate gamma_n budget.
    """
    n=len(a_mid)
    if not (n==len(a_rad)==len(b_mid)==len(b_rad)):
        raise ValueError("equal interval-vector lengths required")
    mid=0.0; rad=0.0
    for am,ar,bm,br in zip(a_mid,a_rad,b_mid,b_rad):
        if not all(math.isfinite(x) for x in (am,ar,bm,br)) or ar<0 or br<0:
            raise ValueError("finite midpoint/radius entries required")
        mid += am*bm
        rad += abs(am)*br + abs(bm)*ar + ar*br
    lo=math.nextafter(mid-rad,-math.inf)
    hi=math.nextafter(mid+rad, math.inf)
    return lo,hi


def symmetric_interval_gershgorin(mid: tuple[tuple[float,...],...],
                                  rad: tuple[tuple[float,...],...]) -> tuple[float,float]:
    """Outward spectral enclosure for a symmetric midpoint-radius matrix."""
    n=len(mid)
    if n==0 or len(rad)!=n or any(len(r)!=n for r in mid) or any(len(r)!=n for r in rad):
        raise ValueError("square midpoint/radius matrices required")
    lower=math.inf; upper=-math.inf
    for i in range(n):
        if rad[i][i] < 0: raise ValueError("nonnegative radii required")
        center_lo=mid[i][i]-rad[i][i]
        center_hi=mid[i][i]+rad[i][i]
        row=0.0
        for j in range(n):
            if i==j: continue
            if rad[i][j] < 0: raise ValueError("nonnegative radii required")
            row += abs(mid[i][j])+rad[i][j]
        lower=min(lower,center_lo-row)
        upper=max(upper,center_hi+row)
    return math.nextafter(lower,-math.inf),math.nextafter(upper,math.inf)
