from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.block_factor_metric import (
    BlockFactorFloor,BlockFactorMetric,active_bias_factor_floor,
    block_information_normalization,block_scaled_schur_floor,
    ba_ou_recurring_factor_floor,cross_block_psd_bound,
    ag_one_second_factor_floor,ba_process_variance,block_factor_feasibility,
    additive_process_metric_certificate,block_metric_correction_margin,
    normalized_measurement_information_ceiling,
    uniform_lin_jensen_factor_certificate,
)
from tools.stability.ou3_theorem.interval_riccati_21 import (
    integrated_ou_scaled_factor_probe,
)

class BlockFactorMetricTests(unittest.TestCase):
    def test_cross_coupling_fails_closed(self):
        m=BlockFactorMetric(BlockFactorFloor(1,2),BlockFactorFloor(1,2),
                            BlockFactorFloor(1,2),.6,.6,.1)
        self.assertFalse(block_scaled_schur_floor(m)["verified"])

    def test_positive_block_margin_normalizes_information(self):
        m=BlockFactorMetric(BlockFactorFloor(.2,2),BlockFactorFloor(.3,3),
                            BlockFactorFloor(.1,1),.005,.002,.002)
        r=block_information_normalization(
            information_floors=(.0021,.00204734,.003),
            metric=m)
        self.assertTrue(r["verified"])
        self.assertGreater(r["mu_cov"],0)
        self.assertLess(r["rho0"],1)

    def test_active_bias_factor_uses_literal_release_or_ou_floor(self):
        self.assertEqual(active_bias_factor_floor(
            stationary_sigma_min=.05,release_std_floor=.004),.004)

    def test_ba_recurring_factor_is_positive(self):
        r=ba_ou_recurring_factor_floor(
            phi_max=math.exp(-.004/5000.0),
            process_variance_floor=1e-9,
            release_variance_floor=1.6e-5,
            corrections_information_ceiling=1e4)
        self.assertGreater(r["factor_floor"],0)

    def test_generic_psd_cross_bound_is_too_weak_for_small_factors(self):
        c=cross_block_psd_bound(2.0,3.0)
        self.assertAlmostEqual(c,math.sqrt(6.0))

    def test_additive_process_factor_gives_root_gamma_one_despite_cross_terms(self):
        r=additive_process_metric_certificate(
            ag_q_factor=5e-4,lin_q_factor=1e-9,ba_q_factor=3e-5)
        self.assertTrue(r["verified"])
        self.assertEqual(r["root_metric_gamma"],1.0)
        j=normalized_measurement_information_ceiling(
            h_norm=10.0,max_factor=5e-4,noise_variance_floor=.0025)
        g=block_metric_correction_margin(gamma_in=1.0,
                                         normalized_information_ceiling=j)
        self.assertGreater(g,0.0);self.assertLessEqual(g,1.0)

    def test_source_uniform_lin_factor_certificate(self):
        r=uniform_lin_jensen_factor_certificate(
            tau_min=.02,tau_max=12.0,sigma_min=.05,window_s=16.0,
            scales=(5.5,8.1,1100.0,4.0),tau_cells=4096)
        print("LIN_UNIFORM_FACTOR",r)
        self.assertTrue(r["verified"])
        self.assertGreater(r["singular_factor_floor"],0)

    def test_source_range_factor_feasibility_diagnostic(self):
        ag=ag_one_second_factor_floor(
            gyro_white_density=.00157,gyro_bias_rw_density=1e-5,
            attitude_scale=1.0,gyro_bias_scale=.02)
        self.assertTrue(ag["verified"])
        lin=[]
        for tau in (.02,.05,.1,.2,.5,1,2,4,8,12):
            r=integrated_ou_scaled_factor_probe(
                tau=tau,sigma=.05,window_s=16.0,
                scales=(5.5,8.1,1100.0,4.0),slices=256)
            lin.append((tau,r["pivot_floor"]))
        qba=ba_process_variance(dt=.004,tau_bacc=5000.0,drive_density=5e-4)
        ba=ba_ou_recurring_factor_floor(
            phi_max=math.exp(-.004/5000.0),
            process_variance_floor=qba,
            release_variance_floor=1.6e-5,
            corrections_information_ceiling=1e4)
        probe=block_factor_feasibility(
            ag_factor=ag["factor_floor"],
            lin_factor=math.sqrt(min(x[1] for x in lin)),
            ba_factor=ba["factor_floor"],
            ag_ceiling=2.0,lin_ceiling=2500.0,ba_ceiling=.16)
        print("BLOCK_FACTOR_SOURCE_PROBE",{"ag":ag,"lin":lin,"qba":qba,"ba":ba,"generic_cross":probe})
        self.assertGreater(min(x[1] for x in lin),0)
        self.assertGreater(ba["factor_floor"],0)

    def test_translation_scaled_factor_probe_is_positive_for_representative_tau(self):
        # Non-promoting feasibility diagnostic for the new metric. The actual
        # theorem must interval-enclose tau and the within-slice remainder.
        r=integrated_ou_scaled_factor_probe(
            tau=1.0,sigma=.05,window_s=16.0,
            scales=(5.5,8.1,1100.0,4.0),slices=128)
        print("BLOCK_FACTOR_PROBE",r)
        self.assertTrue(r["verified"])
        self.assertGreater(r["pivot_floor"],0)

if __name__=="__main__":unittest.main()
