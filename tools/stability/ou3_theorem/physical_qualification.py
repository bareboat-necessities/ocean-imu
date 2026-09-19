"""Single-history physical qualification for the OU-III theorem."""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_theorem.marine_motion import (
    MarineContinuationCertificate,MarineLimits,continuation_admitted as marine_admitted,
)
from tools.stability.ou3_theorem.imu_bias import (
    BiasContinuationCertificate,BiasLimits,continuation_admitted as bias_admitted,
)
from tools.stability.ou3_theorem.magnetic_service import (
    MagneticServiceContinuationCertificate,continuation_admitted as magnetic_admitted,
)


@dataclass(frozen=True)
class PhysicalQualification:
    marine: MarineContinuationCertificate
    bias: BiasContinuationCertificate
    magnetic: MagneticServiceContinuationCertificate
    sensor_datasheet_and_commissioning_qualified: bool
    assembled_mounting_and_calibration_qualified: bool


def qualification_admitted(q: PhysicalQualification, *,
                            marine_limits: MarineLimits,
                            bias_limits: BiasLimits,
                            magnetic_window_s: float,
                            magnetic_information_floor: float) -> dict:
    ids=(q.marine.history_id,q.bias.history_id,q.magnetic.history_id)
    same_history=len(set(ids))==1
    pieces={
        "same_history":same_history,
        "marine":marine_admitted(q.marine,marine_limits),
        "bias":bias_admitted(q.bias,bias_limits),
        "magnetic":magnetic_admitted(
            q.magnetic,required_window_s=magnetic_window_s,
            required_information_floor=magnetic_information_floor),
        "sensor_datasheet_and_commissioning":q.sensor_datasheet_and_commissioning_qualified,
        "assembled_mounting_and_calibration":q.assembled_mounting_and_calibration_qualified,
    }
    pieces["qualified"]=all(pieces.values())
    return pieces
