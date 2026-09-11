"""Analytic operand binding and anti-shortcut regressions, not stability gates."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction.shipping_graph import (
    measurement_graph, quaternion_matrix, same_driver_bias_graph, binding_defects,
)
from tools.stability.ou3_alt_contraction.shipping_word import (
    DTYPE, HEADER, analyze_pair, gain_pair_identity, high_precision_pair,
    instrument_header, instrument_fusion, read_trace, source_contract, bias_contracts, FUSION, verify_transparency, reference_error, state_difference,
)


def record(active=False,kind=10):
    row=np.zeros(1,dtype=DTYPE)[0]
    row['kind']=kind;row['active']=active;row['q']=[1,0,0,0]
    row['P']=np.eye(21,dtype=np.float32).ravel()
    row['R']=np.eye(3,dtype=np.float32).ravel()
    row['parameters']=[9.80665,5000,.4,0,0,.0034641]
    row['measured']=[.1,.2,-9.7];row['mag_reference']=[22,0,43]
    row['source'][15:18]=row['measured']
    graph=measurement_graph(kind,row['P'],row['q'],row['x'],row['measured'],
                            row['R'],active,row['mag_reference'],row['parameters'][0])
    for name in ('N','S','r'):row[name]=graph[name].ravel()
    row['K']=np.linalg.solve(graph['S'],graph['N'].T).T.ravel()
    return row


class ShippingAttachmentTests(unittest.TestCase):
    def test_H18_keeps_bias_uncertainty_but_masks_actual_gain(self):
        a=record(False);b=record(True)
        np.testing.assert_array_equal(a['S'],b['S'])
        self.assertGreater(np.linalg.norm(b['N'].reshape(21,3)[18:]),0)
        np.testing.assert_array_equal(a['N'].reshape(21,3)[18:],0)
        h=measurement_graph(10,a['P'],a['q'],a['x'],a['measured'],a['R'],False,
                            a['mag_reference'],a['parameters'][0])
        self.assertFalse(np.array_equal(h['H_residual'],h['H_gain']))

    def test_dense_same_P_H_R_relation_with_nonzero_cross_blocks(self):
        rng=np.random.default_rng(76)
        b=rng.normal(size=(21,21));p=b@b.T+np.eye(21)
        q=np.array([.98,.1,.04,-.15]);q/=np.linalg.norm(q)
        for active in (False,True):
            h=measurement_graph(10,p,q,np.linspace(.01,.21,21),[1,2,-9],np.eye(3),active,[22,0,43],9.80665)
            np.testing.assert_allclose(h['S'],h['H_residual']@p@h['H_residual'].T+np.eye(3),atol=1e-10)
            expected=p@h['H_gain'].T
            if not active:expected[18:]=0
            np.testing.assert_allclose(h['N'],expected)

    def test_mag_and_S_graphs_use_actual_reference_and_anisotropic_R(self):
        p=np.diag(np.arange(1.,22.));noise=np.diag([.2,.3,.7])
        for kind in (11,12):
            g=measurement_graph(kind,p,[1,0,0,0],np.arange(21.),[1,2,3],noise,True,[20,1,40],9.8)
            np.testing.assert_allclose(g['S'],g['H_residual']@p@g['H_residual'].T+noise)
        g=measurement_graph(12,p,[1,0,0,0],np.arange(21.),[1,2,3],noise,False,[20,1,40],9.8)
        np.testing.assert_array_equal(g['r'],[-12,-13,-14])

    def test_nonzero_nominal_residual_exposes_covariance_memory(self):
        a=record();b=a.copy()
        p=b['P'].reshape(21,21);p[0,0]+=np.float32(.01)
        g=measurement_graph(10,p,b['q'],b['x'],b['measured'],b['R'],False,b['mag_reference'],b['parameters'][0])
        for name in ('N','S','r'):b[name]=g[name].ravel()
        b['K']=np.linalg.solve(g['S'],g['N'].T).T.ravel()
        s=gain_pair_identity(a,b)
        self.assertEqual(s['frozen_correction_norm'],0.)
        self.assertGreater(s['endogenous_correction_norm'],1e-9)
        self.assertLess(s['equality_max_abs'],1e-12)
        self.assertLess(s['correction_identity_max_abs'],1e-14)
        for digits in (80,120):
            h=high_precision_pair(a,b,digits)
            self.assertLess(float(h['row_graph_residual']),1e-60)
            self.assertIn('NOT high-precision nonlinear',h['scope'])

    def test_lower_triangle_solve_does_not_replace_raw_Joseph_S(self):
        a=record();b=a.copy();b['S'][1]+=np.float32(.25)
        s=gain_pair_identity(a,b)
        self.assertEqual(s['dS_norm'],0.)
        self.assertGreater(s['raw_S_asymmetry_max_abs'],0.)
        self.assertLess(s['equality_max_abs'],1e-12)

    def test_binding_uses_pre_state_not_reconstructed_H(self):
        a=record();d=binding_defects(a)
        self.assertLess(max(d[k+'_scaled_defect'] for k in ('N','S','r')),1e-5)
        a['N'][0]+=.5
        self.assertGreater(binding_defects(a)['N_absolute_defect'],.49)

    def test_nonzero_thermal_coefficient_is_allowed_at_reference_temperature(self):
        a=record();a['parameters'][5]=.0034641018
        self.assertLess(binding_defects(a)['r_absolute_defect'],1e-5)

    def test_out_of_scope_lever_branch_is_not_silently_accepted(self):
        a=record();a['parameters'][3]=1
        with self.assertRaises(ValueError):binding_defects(a)

    def test_joint_bias_prediction_retains_same_driver_and_mismatch(self):
        b0=np.array([.04,-.02,.03]);b1=b0+np.array([.001,0,-.001])
        bh0=np.array([.01,-.01,.02])
        for pt in (.99,1.):
            for ph in (.98,1.):
                w,defect=same_driver_bias_graph(b0,b1,bh0,ph*bh0,pt,ph)
                np.testing.assert_allclose(defect,0,atol=1e-17)
                np.testing.assert_allclose(pt*b0+w,b1)

    def test_quaternion_polynomial_does_not_normalize_away_roundoff(self):
        q=np.array([.99,.02,.03,.04])
        self.assertGreater(np.max(np.abs(quaternion_matrix(q)-quaternion_matrix(q/np.linalg.norm(q)))),1e-6)

    def test_mismatched_hybrid_words_are_not_zipped_or_truncated(self):
        a=np.array([record()],dtype=DTYPE);b=np.array([record(kind=11)],dtype=DTYPE)
        r=analyze_pair(a,b)
        self.assertFalse(r['event_paths_equal'])
        self.assertIn('MISMATCH',r['classification'])

    def test_hooks_are_temporary_and_fail_closed_on_anchor_changes(self):
        original=(ROOT/HEADER).read_text();patched=instrument_header(original)
        self.assertEqual((ROOT/HEADER).read_text(),original)
        self.assertIn('alt::solve(10',patched)
        self.assertIn('alt::ExitGuard alt_guard(120',patched)
        self.assertIn('alt::point(70',patched)
        with self.assertRaises(ValueError):
            instrument_header(original.replace('xext.noalias() += K * r;','xext += K * r;'))
        with self.assertRaises(ValueError):
            instrument_header(original.replace('tempC_ref = T(35.0)','tempC_ref = T(25.0)'))

    def test_all_three_bias_definitions_are_invoked_separately(self):
        r=bias_contracts()
        self.assertEqual(set(r),{'BIAS0','BIAS1','BIAS2'})
        for v in r.values():
            self.assertTrue(v['definition_validator_pass'])
            self.assertTrue(v['analytic_example_parameter_membership'])
            self.assertFalse(v['fresh_Live_admission'])

    def test_actual_continuous_hard_iron_is_kept_in_input_binding(self):
        a=record(kind=11);a['hard_iron']=[.02,-.01,.03]
        a['source'][21:24]=a['measured']+a['hard_iron']
        self.assertLess(binding_defects(a)['physical_input_max_abs_defect'],1e-6)
        a['hard_iron']=0
        self.assertGreater(binding_defects(a)['physical_input_max_abs_defect'],.029)
        original=(ROOT/FUSION).read_text()
        self.assertIn('alt::set_hard_iron',instrument_fusion(original))
        self.assertEqual(original,(ROOT/FUSION).read_text())

    def test_transparency_comparison_detects_one_bit_of_state_change(self):
        import json
        with tempfile.TemporaryDirectory() as d:
            left,right=Path(d)/'left',Path(d)/'right'
            left.mkdir();right.mkdir()
            manifest={'words':[{'name':'sample'}]}
            for folder in (left,right):
                (folder/'manifest.json').write_text(json.dumps(manifest))
                for direction in ('base','tilt','held_bias','covariance'):
                    a=np.array([record()],dtype=DTYPE);a['kind']=40
                    a.tofile(folder/('sample-'+direction+'.bin'))
            self.assertTrue(verify_transparency(left,right)['all_sample_records_bit_identical'])
            a['x'][0,3]=1e-9;a.tofile(right/'sample-base.bin')
            with self.assertRaises(ValueError):verify_transparency(left,right)

    def test_reference_error_keeps_physical_truth_and_one_S_origin(self):
        a=record();a['kind']=40
        a['source'][:12]=np.arange(12,dtype=np.float32)
        a['source'][18:21]=[.03,-.02,.01]
        e=reference_error(a)
        np.testing.assert_array_equal(e[6:9],[3,4,5])
        np.testing.assert_array_equal(e[9:12],[0,1,2])
        np.testing.assert_array_equal(e[12:15],[9,10,11])
        np.testing.assert_array_equal(e[18:21],e[21:24])
        b=a.copy();b['x'][18]+=.01
        self.assertAlmostEqual(state_difference(a,b)[18],-.01,places=8)
        np.testing.assert_array_equal(state_difference(a,b)[21:24],0)
        b['kind']=50
        with self.assertRaises(ValueError):reference_error(b)

    def test_truncated_trace_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'broken.bin';p.write_bytes(b'bad')
            with self.assertRaises(ValueError):read_trace(p)

    def test_analytic_source_is_not_promoted_to_COMPLETE_BRMM(self):
        c=source_contract()
        self.assertTrue(c['kinematic_envelope_sufficient'])
        self.assertFalse(c['BRMM_full_admission'])
        self.assertFalse(c['PE_admission'])
        self.assertTrue(c['probes_are_not_source_cover'])


if __name__=='__main__':unittest.main()
