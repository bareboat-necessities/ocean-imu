from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.interval_riccati_21 import (
    IMat,N,accel_bias_release_event,add,aw_covariance_floor_event,
    innovation_covariance,innovation_inverse_spectral_certificate,
    joseph_covariance,matmul,predict_covariance,spectral_box,
    shipping_acc_update_intervals,shipping_integral_update_intervals,
    shipping_mag_update_intervals,shipping_prediction_intervals,
    verified_gain_interval,verified_joseph_update,iterate_recurring_box,
    shipping_max_correction_step,shipping_word_map,split_interval_matrix,
    adaptive_verified_update,recurring_box_over_cells,exact_midpoint_seed,
    interval_failure_metrics,
)


def exact(rows):
    return IMat(tuple(tuple(float(x) for x in r) for r in rows),
                tuple(tuple(0.0 for _ in r) for r in rows))


def diag(n,x):
    return exact([[x if i==j else 0.0 for j in range(n)] for i in range(n)])


class IntervalRiccati21Tests(unittest.TestCase):
    def test_prediction_keeps_cross_covariance(self):
        p=diag(N,1.0)
        frows=[[1.0 if i==j else 0.0 for j in range(N)] for i in range(N)]
        frows[0][3]=.5
        q=diag(N,.1)
        out=predict_covariance(p,exact(frows),q)
        self.assertAlmostEqual(out.mid[0][3],.5)
        self.assertAlmostEqual(out.mid[3][0],.5)
        self.assertAlmostEqual(out.mid[0][0],1.35)

    def test_integral_innovation_and_joseph(self):
        p=diag(N,2.0)
        hrows=[[0.0]*N for _ in range(3)]
        for a in range(3): hrows[a][12+a]=1.0
        h=exact(hrows); r=diag(3,1.0)
        s=innovation_covariance(p,h,r)
        cert=innovation_inverse_spectral_certificate(s)
        self.assertTrue(cert["verified"])
        krows=[[0.0]*3 for _ in range(N)]
        for a in range(3): krows[12+a][a]=2.0/3.0
        out=joseph_covariance(p,exact(krows),h,r)
        self.assertAlmostEqual(out.mid[12][12],2.0/3.0,places=12)

    def test_interval_innovation_can_fail_closed(self):
        p=diag(N,1.0)
        hrows=[[0.0]*N for _ in range(3)]
        for a in range(3): hrows[a][a]=1.0
        h=exact(hrows)
        r=IMat(((.01,0,0),(0,.01,0),(0,0,.01)),
               ((2.0,0,0),(0,2.0,0),(0,0,2.0)))
        cert=innovation_inverse_spectral_certificate(innovation_covariance(p,h,r))
        self.assertFalse(cert["verified"])

    def test_hard_covariance_events_are_enclosed(self):
        p=diag(N,.01)
        aw=aw_covariance_floor_event(p,.25)
        for i in range(15,18):
            self.assertGreaterEqual(aw.mid[i][i]+aw.rad[i][i],.25)
        ba=accel_bias_release_event(p,.04)
        for i in range(18,21):
            self.assertGreaterEqual(ba.mid[i][i]-ba.rad[i][i],.04)

    def test_shipping_prediction_interval_has_literal_blocks(self):
        f,q=shipping_prediction_intervals(
            dt_min=.004,dt_max=.006,tau_min=.02,tau_max=12.0,
            omega_max=.6108652382,tau_bacc=5000.0,
            gyro_white_density=.00157,gyro_bias_rw_density=1e-5,
            aw_sigma_max=4.0,accel_bias_drive_density=5e-4)
        self.assertEqual(f.shape,(N,N)); self.assertEqual(q.shape,(N,N))
        self.assertEqual(f.mid[3][3],1.0)
        self.assertLess(f.rad[0][0],1e-4)
        self.assertLess(f.rad[0][1],.004)
        self.assertLess(f.rad[0][3],.0021)
        self.assertGreater(f.mid[18][18]-f.rad[18][18],0.999)
        self.assertGreater(q.mid[15][15]+q.rad[15][15],0.0)

    def test_verified_gain_and_joseph_update(self):
        p=diag(N,2.0)
        hrows=[[0.0]*N for _ in range(3)]
        for a in range(3): hrows[a][12+a]=1.0
        h=exact(hrows);r=diag(3,1.0)
        k,cert=verified_gain_interval(p,h,r)
        self.assertTrue(cert["verified"])
        self.assertAlmostEqual(k.mid[12][0],2/3,places=12)
        out,cert2=verified_joseph_update(p,h,r)
        self.assertTrue(cert2["verified"])
        self.assertAlmostEqual(out.mid[12][12],2/3,places=12)

    def test_literal_measurement_interval_constructors(self):
        hs,rs=shipping_integral_update_intervals(.01,100.0)
        self.assertEqual(hs.mid[0][12],1.0)
        self.assertGreater(rs.mid[0][0]+rs.rad[0][0],9999.0)
        hm,rm=shipping_mag_update_intervals(75.0,.1,2.0)
        self.assertEqual(hm.rad[0][1],75.0)
        self.assertEqual(hm.rad[0][3],0.0)
        ha,ra=shipping_acc_update_intervals(18.7,.05,.31)
        self.assertEqual(ha.mid[0][18],1.0)
        self.assertEqual(ha.rad[0][15],1.0)
        self.assertGreater(ra.mid[0][0],0.0)

    def test_literal_event_order_prediction_S_acc_mag(self):
        p=diag(N,2.0);fmat=diag(N,1.0);q=diag(N,.01)
        hsrows=[[0.0]*N for _ in range(3)]
        harows=[[0.0]*N for _ in range(3)]
        hmrows=[[0.0]*N for _ in range(3)]
        for a in range(3):
            hsrows[a][12+a]=1.0;harows[a][15+a]=1.0;hmrows[a][a]=1.0
        out,certs=shipping_max_correction_step(
            p,f=fmat,q=q,h_s=exact(hsrows),r_s=diag(3,1.0),
            h_acc=exact(harows),r_acc=diag(3,1.0),
            h_mag=exact(hmrows),r_mag=diag(3,1.0))
        self.assertEqual([x["name"] for x in certs],["S","acc","mag"])
        self.assertTrue(all(x["verified"] for x in certs))
        self.assertLess(out.mid[0][0],p.mid[0][0])

    def test_recurring_box_requires_positive_self_inclusion(self):
        seed=diag(N,1.0)
        def word(box):
            return box,[{"verified":True,"name":"synthetic"}]
        out=iterate_recurring_box(seed,word,max_iterations=2)
        self.assertTrue(out["verified"])
        self.assertGreater(out["spectral_lower"],0.0)

    def test_recurring_box_fails_on_unverified_innovation(self):
        seed=diag(N,1.0)
        def word(box):
            return box,[{"verified":False,"name":"synthetic"}]
        self.assertFalse(iterate_recurring_box(seed,word)["verified"])

    def test_split_cells_cover_original_entry(self):
        a=diag(N,1.0)
        m=[list(r) for r in a.mid];q=[list(r) for r in a.rad]
        m[0][1]=0.0;q[0][1]=2.0
        x=IMat(tuple(tuple(r) for r in m),tuple(tuple(r) for r in q))
        cells=split_interval_matrix(x,((0,1),))
        self.assertEqual(len(cells),2)
        lo=min(z.mid[0][1]-z.rad[0][1] for z in cells)
        hi=max(z.mid[0][1]+z.rad[0][1] for z in cells)
        self.assertLessEqual(lo,-2.0);self.assertGreaterEqual(hi,2.0)

    def test_recurring_box_requires_all_branch_cells(self):
        seed=diag(N,1.0)
        def good(p): return p,[{"verified":True}]
        def bad(p):
            out=diag(N,2.0)
            return out,[{"verified":True}]
        r=recurring_box_over_cells(seed,(good,bad),max_iterations=1)
        self.assertFalse(r["verified"])
        self.assertEqual(r["branch_count"],2)

    def test_adaptive_update_accepts_exact_cell_without_split(self):
        p=diag(N,1.0)
        hrows=[[0.0]*N for _ in range(3)]
        for a in range(3): hrows[a][12+a]=1.0
        out,certs=adaptive_verified_update(
            p,exact(hrows),diag(3,1.0),split_entries=((0,12),),max_depth=1)
        self.assertEqual(out.shape,(N,N))
        self.assertEqual(certs[0]["depth"],0)

    def test_failure_metrics_are_nonpromoting_and_numeric(self):
        p=exact_midpoint_seed((1.0,)*N)
        hrows=[[0.0]*N for _ in range(3)]
        for a in range(3): hrows[a][12+a]=1.0
        m=interval_failure_metrics(p,exact(hrows),diag(3,1.0))
        self.assertTrue(m["innovation_verified"])
        self.assertAlmostEqual(m["covariance_spectral_lower"],1.0,places=14)

    def test_spectral_box_is_finite(self):
        lo,hi=spectral_box(diag(N,1.0))
        self.assertLessEqual(lo,1.0); self.assertGreaterEqual(hi,1.0)


if __name__=="__main__": unittest.main()
