from __future__ import annotations
import json
import sys
import tempfile
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))

from tools.stability.ou3_theorem.interval_riccati import (
    InnovationCertificate, IntervalRiccatiCertificate,
    RiccatiImageCertificate, load_certificate,
)


class IntervalRiccatiCertificateTests(unittest.TestCase):
    def certificate(self, **kw):
        d=dict(
            qualification="OU3_A21_INTERVAL_RICCATI_V1",
            fixed_coordinate_mu=.00204,
            innovations=(InnovationCertificate("acc",.2,.01),
                         InnovationCertificate("S",.3,.01),
                         InnovationCertificate("mag",.4,.01)),
            recurring_image=RiccatiImageCertificate(.1,10.0,.2,9.0,.01),
            covariance_hard_events_included=True,
            full_shipping_schedule_included=True,
            float32_rounding_included=True,
        )
        d.update(kw)
        return IntervalRiccatiCertificate(**d)

    def test_complete_certificate_promotes_rho(self):
        r=self.certificate().verify()
        self.assertTrue(r["certificate_complete"])
        self.assertGreater(r["rho_chain"]["mu_cov"],0)
        self.assertLess(r["rho_chain"]["rho0"],1)

    def test_any_unverified_inverse_fails_closed(self):
        c=self.certificate(innovations=(
            InnovationCertificate("acc",.1,.1),))
        self.assertFalse(c.verify()["certificate_complete"])

    def test_missing_schedule_or_float32_scope_fails_closed(self):
        self.assertFalse(self.certificate(full_shipping_schedule_included=False).verify()["certificate_complete"])
        self.assertFalse(self.certificate(float32_rounding_included=False).verify()["certificate_complete"])

    def test_noninvariant_box_fails_closed(self):
        image=RiccatiImageCertificate(.1,10.0,.05,9.0,.01)
        self.assertFalse(self.certificate(recurring_image=image).verify()["certificate_complete"])

    def test_json_loader(self):
        raw={
            "qualification":"OU3_A21_INTERVAL_RICCATI_V1",
            "fixed_coordinate_mu":.00204,
            "innovations":[{"name":"acc","midpoint_min_eigenvalue":.2,
                            "spectral_radius_bound":.01}],
            "recurring_image":{"proposed_lower":.1,"proposed_upper":10.0,
                               "image_lower":.2,"image_upper":9.0,
                               "outward_rounding_slack":.01},
            "covariance_hard_events_included":True,
            "full_shipping_schedule_included":True,
            "float32_rounding_included":True,
        }
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"c.json"; p.write_text(json.dumps(raw))
            self.assertTrue(load_certificate(p).verify()["certificate_complete"])


if __name__=="__main__":
    unittest.main()
