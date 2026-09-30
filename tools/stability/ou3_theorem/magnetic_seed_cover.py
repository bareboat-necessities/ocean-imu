"""Theorem-domain seed audit for the exhaustive magnetic certificate.

A seed may be emitted only from CLOSED theorem bounds.  This module records
which continuous/discrete coordinates have such bounds and refuses to invent
a full covariance box from carried histories.
"""
from __future__ import annotations

def magnetic_seed_audit():
    items={
      "dt":{"closed":True,"range":[0.004,0.006],"source":"IMU BIAS / regular A21"},
      "tau":{"closed":True,"range":[0.02,12.0],"source":"shipping tuner clamps"},
      "reset_dtheta":{"closed":True,"norm_upper_rad":0.10471975511965977,
                      "source":"retained 6-degree local domain"},
      "magnetic_service":{"closed":True,"window_s":1.0,"mu_M":1.0,
                          "source":"MAGNETIC SERVICE"},
      "aw_covariance":{"closed":True,"spectral_upper":16.0,
                       "source":"shipping AW synchronization target"},
      "ba_covariance":{"closed":True,"spectral_upper":1.0/1600.0,
                       "source":"proved BA ceiling"},
      "full_covariance_P":{"closed":False,
          "source":"AG6/full recurring upper covariance remains open",
          "needed_for":["S_m_actual=H_m P H_m'+Rmag","K=P H_m' S^-1","I-KH"]},
      "predicted_field_vector":{"closed":True,
          "source":"committed magnetic reference rotated by attitude",
          "note":"norm fixed by committed field; orientation covered by retained attitude chart"},
      "event_schedule":{"closed":False,
          "source":"MAGNETIC SERVICE constrains aggregate accepted information, not a finite callback pattern",
          "note":"arbitrarily many 4-6ms callback patterns exist; must parameterize by service rows, not enumerate schedules"},
    }
    blockers=[k for k,v in items.items() if not v["closed"]]
    return {"qualification":"OU3_MAGNETIC_SEED_AUDIT_V1",
            "verified":False,"seed_cover_constructible":False,
            "blockers":blockers,"items":items,
            "reason":"current theorem assumptions do not supply a compact full-P box or finite event-pattern cover"}

def theorem_seed_cover():
    """Fail closed until the two blockers in magnetic_seed_audit are resolved."""
    return (),magnetic_seed_audit()
