"""Fail-closed tests for the actual-word experiment, never theorem gates."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import mpmath as mp
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/stability/ou3_alt_contraction'))
import shipping_probe as probe
import word_diagnostic as diagnostic


class ShippingWordTests(unittest.TestCase):
    def test_overlay_is_only_added_observations(self):
        original=(ROOT/probe.HEADER).read_text()
        observed=probe.overlay(original)
        self.assertEqual(probe.strip_observations(observed),original)
        self.assertEqual(observed.count('    alt_solve('),3)
        self.assertEqual(observed.count('    AltScope alt_scope('),4)
        for tag in ('prediction','aw_floor','joseph','reset','projection'):
            self.assertIn('alt_state("'+tag+'"',observed)

    def test_changed_anchor_fails_closed(self):
        original=(ROOT/probe.HEADER).read_text()
        with self.assertRaisesRegex(ValueError,'anchor changed'):
            probe.overlay(original.replace('Vector3 const& gyr_body, T Ts)',
                                           'Vector3 const& gyr_body, T changed)'))

    def test_joint_metric_is_fixed_coercive_and_keeps_bias_cross(self):
        M=diagnostic.METRIC
        self.assertEqual(M.shape,(24,24))
        self.assertGreater(np.linalg.eigvalsh(M).min(),0)
        self.assertEqual(M[18,21],25)
        self.assertAlmostEqual(np.linalg.eigvalsh(M).min(),.01)

    def test_physical_probe_has_all_time_kinematic_envelope(self):
        d=diagnostic.analytic_probe_envelope()
        self.assertTrue(d['kinematic_caps_satisfied'])
        self.assertFalse(d['full_source_cover'])
        self.assertFalse(d['full_BRMM_PE_sensor_startup_admission'])
        self.assertLess(d['BIAS2_derivative_bound'],.0002)
        self.assertLess(d['BIAS2_component_bound'],.13)

    def test_pending_attitude_extension_and_high_precision_chart(self):
        a=np.zeros((2,37));a[:,0]=[1000,1001];a[:,2]=1
        a[:,6]=[1e-4,2e-4]
        t=diagnostic.Tape(np.array(['word_start','word_end']),a,2,-1,.001,1.2,1000,800,0)
        e=diagnostic.errors(t)
        self.assertTrue(np.isfinite(e).all())
        for j in (0,1):
            with mp.workdps(80):
                ref=np.array(list(diagnostic.mp_energy_of_endpoint(t,j)),dtype=float)
            np.testing.assert_allclose(ref,e[j],atol=3e-14,rtol=3e-14)
        self.assertFalse(np.allclose(diagnostic.effective_quaternion(a),a[:,2:6]))

    def test_event_mismatch_is_not_silently_discarded(self):
        a=np.zeros((2,37));a[:,0]=[1000,1600];a[:,2]=1
        b=a.copy();b[:,12]=[-.001,-.002]
        tags=np.array(['word_start','word_end'])
        base=diagnostic.Tape(tags,a,2,-1,.001,1.2,1000,800,-1)
        other=diagnostic.Tape(np.array(['word_start','different_event']),b,2,-1,.001,1.2,1000,800,6)
        result=diagnostic.compare(base,other)
        self.assertFalse(result['same_event_word'])
        self.assertIsNone(result['signed_delta_by_operation'])
        self.assertFalse(result['actual_full_joint24_map_certified'])

    def test_solve_graph_retains_endogenous_nominal_terms(self):
        with tempfile.TemporaryDirectory() as td:
            paths=[]
            for j in range(2):
                path=Path(td)/str(j)
                r=np.array([1.,2.,3.])+j*np.array([.1,.2,-.1])
                N=np.arange(63,dtype=float).reshape(21,3)/100+j/50
                S=np.eye(3)*(2+j/10)
                K=N@np.linalg.inv(S)
                # An observed floating correction defect is NOT discarded.
                if j:K[0,0]+=1e-6
                P=np.eye(21)*(1+j/10)
                row=['acc','1','0']+[repr(float(x)) for v in (r,N,S,K,P,np.zeros(21),(K.astype(np.float32)@r.astype(np.float32)).astype(float)) for x in v.flatten()]
                path.write_text(','.join(row)+'\n');paths.append(path)
            d=diagnostic.solve_graph(*paths)
            self.assertTrue(d['aligned'])
            self.assertGreater(d['maxima']['delta_N_times_q0_norm'],0)
            self.assertGreater(d['maxima']['delta_S_times_q0_norm'],0)
            self.assertGreater(d['maxima']['actual_mean_minus_ideal_Nq_increment_norm'],1e-7)
            self.assertLess(d['maxima']['descriptor_residual_norm'],1e-70)
            self.assertFalse(d['deployment_roundoff_enclosed'])

    def test_no_final_gate_can_be_promoted(self):
        d=diagnostic.open_obligations()
        for key in ('ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS','P4_promoted','P5_promoted'):
            self.assertIs(d[key],False)

if __name__=='__main__':unittest.main()
