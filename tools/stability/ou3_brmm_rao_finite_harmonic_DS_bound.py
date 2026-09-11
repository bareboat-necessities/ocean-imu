#!/usr/bin/env python3
"""Hard centered-S bound for the shipped 28-ft finite-harmonic RAO data source.

This is a *subset qualification*, not COMPLETE-BRMM closure.  It derives a
pathwise bound from the actual oceanography-waves-lib v1.2.1 construction used
by the repository's 28-ft RAO dataset:

* N=128 incident harmonics;
* f in [0.02,0.8] Hz;
* Hs <= 8.5 m;
* sum a_i^2 <= Hs^2/8 from the deterministic Hs normalization;
* VesselRao translation gain is bounded analytically from its transfer model.

No PSD/high-probability argument is used.  The broader COMPLETE-BRMM source also
admits separately certified shaping generators, so this number MUST NOT be
promoted to the whole family unless those generators are numerically bounded by
an equal-or-stronger certificate.
"""
from __future__ import annotations

import json
import math

QUALIFICATION="OU3_BRMM_28FT_RAO_FINITE_HARMONIC_DS_SUBSET_V1"
N=128
F_MIN_HZ=0.02
F_MAX_HZ=0.8
HS_MAX_M=8.5
HEAVE_DAMPING=0.45
LEGACY_P4_S_RADIUS_M_S=300.0


def _up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def build() -> dict:
    z=HEAVE_DAMPING
    if not 0.0 < z < 1.0/math.sqrt(2.0):
        raise RuntimeError("expected underdamped heave response")
    # max_r |1/(1-r^2-i 2 z r)| = 1/(2 z sqrt(1-z^2)).
    heave_gain=_up(1.0/(2.0*z*math.sqrt(1.0-z*z)))
    # Surge/sway share direction cos/sin and each lowpass/footprint magnitude <=1,
    # hence their joint Euclidean gain is <=1, not sqrt(2).
    translation_gain=_up(math.sqrt(1.0+heave_gain*heave_gain))
    # variance = Hs^2/16 = 0.5 sum a_i^2 (or lower after the bisection guard).
    l2_amp=_up(HS_MAX_M/math.sqrt(8.0))
    l1_amp=_up(math.sqrt(float(N))*l2_amp)
    disp_amplitude_mass=_up(translation_gain*l1_amp)
    omega_min=math.nextafter(2.0*math.pi*F_MIN_HZ, -math.inf)
    D_S=_up(2.0*disp_amplitude_mass/omega_min)
    return {
      "qualification":QUALIFICATION,
      "source":"oceanography-waves-lib v1.2.1 28-ft VesselRao over incidentHarmonics",
      "deterministic_pathwise":True,
      "probabilistic_or_PSD_argument_used":False,
      "N_harmonics":N,
      "frequency_support_hz":[F_MIN_HZ,F_MAX_HZ],
      "omega_min_rad_s_lower":omega_min,
      "Hs_max_m":HS_MAX_M,
      "incident_amplitude_l2_upper_m":l2_amp,
      "incident_amplitude_l1_upper_m":l1_amp,
      "heave_transfer_gain_upper":heave_gain,
      "translation_vector_transfer_gain_upper":translation_gain,
      "displacement_vector_amplitude_mass_M0_upper_m":disp_amplitude_mass,
      "centered_primitive_D_S_upper_m_s":D_S,
      "formula":"D_S <= 2*M0/omega_min; M0 <= G_trans*sqrt(N)*Hs/sqrt(8)",
      "legacy_300_m_s_radius_is_proved_sufficient_for_this_bound":D_S <= LEGACY_P4_S_RADIUS_M_S,
      "legacy_300_m_s_radius_m_s":LEGACY_P4_S_RADIUS_M_S,
      "complete_BRMM_numeric_qualification_closed":False,
      "reason_not_complete":"COMPLETE-BRMM also admits bounded shaping-generator realizations whose numerical potential diameter is not yet uniformly qualified",
      "P4_PASS":False,
      "P5_MAY_START":False,
    }


def validate(d:dict)->list[str]:
    f=[]
    for k in ("deterministic_pathwise",):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("probabilistic_or_PSD_argument_used","complete_BRMM_numeric_qualification_closed","P4_PASS","P5_MAY_START"):
        if d.get(k) is not False:f.append(k+" not false")
    x=float(d.get("centered_primitive_D_S_upper_m_s",math.nan))
    if not (math.isfinite(x) and 850.0 < x < 880.0):f.append("unexpected 28-ft RAO D_S bound")
    if d.get("legacy_300_m_s_radius_is_proved_sufficient_for_this_bound") is not False:f.append("300 m*s was incorrectly promoted")
    return f


def main()->int:
    d=build();f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f
    print(json.dumps(d,indent=2,sort_keys=True));return int(bool(f))

if __name__=="__main__":raise SystemExit(main())
