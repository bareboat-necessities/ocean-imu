import sys,unittest
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.signed_temporal import (
    atoms,
    adjoint_compatibility,
    balanced_gyro_weights,
    certificate,
    endpoint_annihilating_multiplier_constraints,
    endpoint_cancelled_source_bound,
    endpoint_jets_vanish,
    forced_adjoint_source_bound,
    forced_data_adjoint,
    gyro_zero_mean_companion,
    gyro_construction_barrier,
    literal_signed_functional_bound,
    margin_to_Bstar_theorem,
    physical_bounds,
    physical_interval_sampling_error,
    projection_sector_gap,
    separated_reader_action_implication,
    source_margin_attempt,
    sampled_tilt_span_lower,
    signed_bias_terms,
    signed_acceleration_supply,
    zero_terminal_homogeneous_adjoint_impossible,
)

class SignedTemporalTests(unittest.TestCase):
    def test_spline_atoms_and_endpoint_jets(self):
        k=(F(0),F(1),F(2),F(3))
        self.assertEqual(sum(atoms(k)),0)
        self.assertTrue(endpoint_jets_vanish(k))
    def test_homogeneous_adjoint_failure_is_exact(self):
        self.assertTrue(zero_terminal_homogeneous_adjoint_impossible(3))
    def test_physical_margin_positive(self):
        b=physical_bounds()
        self.assertGreater(b["physical_floor_margin"],0)
    def test_projection_gap_positive(self):
        self.assertGreater(projection_sector_gap(),0)
    def test_margin_to_Bstar_implication(self):
        t=margin_to_Bstar_theorem()
        self.assertFalse(t["implication_closed"])
        self.assertTrue(t["six_pivot_floor_conditional_implication_closed"])
        self.assertFalse(t["premise_margins_source_uniformly_certified"])
        b=separated_reader_action_implication(F(1,200),F(2),F(3),F(5),F(7),10)
        self.assertTrue(b["B_star_finite"])
        self.assertGreater(b["B_star_scalar_ceiling"],0)
        for bad in (0,-1,float('nan'),float('inf')):
            with self.subTest(pivot=bad), self.assertRaises((ValueError,OverflowError)):
                separated_reader_action_implication(bad,2,3,5,7,10)
        with self.assertRaises(TypeError):
            separated_reader_action_implication(F(1,200),2,3,5,7,1.5)

    def test_source_margin_attempt_records_exact_gap(self):
        a=source_margin_attempt()
        self.assertGreater(a["physical_collinearity_reserve"],0)
        self.assertIsNone(a["forced_adjoint_signed_transfer_ceiling"])
        self.assertIsNone(a["Delta_col_lower"])
        self.assertFalse(a["new_physical_assumption_needed"])

    def test_literal_adjoint_eliminates_innovation_energy(self):
        b=forced_adjoint_source_bound()
        self.assertFalse(b["innovation_energy_needed"])
        self.assertTrue(b["BA_endpoint_bounded"])
        self.assertFalse(b["BG_endpoint_bounded"])
        self.assertFalse(b["AW_endpoint_bounded"])
        self.assertIsNone(b["finite_numeric_ceiling"])
        self.assertEqual(literal_signed_functional_bound(F(2),F(3),F(4),F(5)),F(32))

    def test_endpoint_annihilation(self):
        a=endpoint_annihilating_multiplier_constraints()
        self.assertTrue(a["four_S_spline_satisfies_boundary_jets"])
        self.assertFalse(a["coupled_estimator_endpoint_cancellation_certified"])
        b=balanced_gyro_weights([1,2,4,8])
        self.assertEqual(sum(b),0)
        g=gyro_zero_mean_companion(b)
        self.assertTrue(g["endpoint_zero"])
        self.assertEqual(endpoint_cancelled_source_bound(2,3,4,F(1,100000),64,
                         state_residual_bound=5,innovation_residual_bound=7),
                         F(18)+F(256,100000))

    def test_exact_compatibility_and_zero_mean_failure(self):
        from tools.stability.ou3_theorem.matrix_certificates import matmul
        # Noncommuting chronological transports are already included in C.
        c=[[1,2,0,3],[0,1,1,2]]
        z=[[2,-3],[1,4]]
        w=matmul(z,c)
        found=adjoint_compatibility(c,w)
        self.assertTrue(found['compatible'])
        self.assertEqual(matmul(found['terminal_multiplier'],c),w)
        # Zero temporal mean by itself destroys a compatible adjoint.
        bad=adjoint_compatibility([[1,2]],[balanced_gyro_weights([1,2])])
        self.assertFalse(bad['compatible'])
        v=[[x] for x in bad['kernel_witness']]
        self.assertEqual(matmul([[1,2]],v),[[0]])
        self.assertNotEqual(matmul([balanced_gyro_weights([1,2])],v),[[0]])

    def test_signed_bias_cancellation_keeps_increment_history(self):
        c=[2,-3,1]; b=[F(1,50),F(2001,100000),F(1999,100000)]
        s=signed_bias_terms(c,b)
        self.assertEqual(s['boundary'],0)
        self.assertEqual(s['increments'],sum(x*y for x,y in zip(c,b)))
        shifted=signed_bias_terms(c,[x+100 for x in b])
        self.assertEqual(s,shifted)
        self.assertEqual(s['tail_weights'],[-2,1])

    def test_forced_data_adjoint_keeps_root_rounding_and_projection(self):
        from tools.stability.ou3_theorem.matrix_certificates import add, matmul
        # Noncommuting operations, incompatible arbitrary desired weights,
        # physical input, innovation arithmetic and a separate projection.
        ops=[{'A':[[1,1],[0,F(1,2)]], 'K':[[0],[0]], 'H':[[0,0]]},
             {'A':[[1,0],[0,1]], 'K':[[2],[-1]], 'H':[[1,3]]},
             {'A':[[1,0],[0,1]], 'K':[[0],[0]], 'H':[[0,0]]},
             {'A':[[1,0],[0,1]], 'K':[[1],[2]], 'H':[[-2,1]]}]
        weights=[[[0]],[[3]],[[0]],[[-2]]]
        y=[[[0]],[[5]],[[0]],[[-4]]]
        eps=[[[0]],[[F(1,100)]],[[0]],[[F(-1,200)]]]
        defects=[[[F(1,50)],[0]],[[0],[F(-1,1000)]],
                 [[0],[F(1,3)]],[[F(1,700)],[0]]]
        state=[[F(2)],[F(-1)]]; root=state
        direct=[[F(0)]]
        for op,w,yi,ei,di in zip(ops,weights,y,eps,defects):
            r=add(add(yi,ei),matmul(op['H'],state),-1)
            direct=add(direct,matmul(w,r))
            state=add(add(matmul(op['A'],state),matmul(op['K'],r)),di)
        adj=forced_data_adjoint(ops,weights)
        result=matmul(adj['Z'][0],root)
        for i in range(len(ops)):
            result=add(result,matmul(adj['L'][i],add(y[i],eps[i])))
            result=add(result,matmul(adj['Z'][i+1],defects[i]))
        self.assertEqual(result,direct)
        self.assertNotEqual(matmul(adj['Z'][0],root),[[0]])
        self.assertEqual(adj['Z'][-1],[[0,0]])

    def test_forced_data_adjoint_does_not_supply_gyro_feedback(self):
        # First column is the carried BG mean. Literal Euclidean sensor maps
        # have zero BG column when the shipping lever arm is disabled.
        ops=[{'A':[[1,0],[0,F(9,10)]], 'K':[[0],[0]], 'H':[[0,0]]},
             {'A':[[1,0],[0,1]], 'K':[[2],[3]], 'H':[[0,1]]}]
        result=forced_data_adjoint(ops,[[[0]],[[1]]])
        self.assertTrue(all(z[0][0]==0 for z in result['Z']))
        # Even a large gyro gain leaves (A-KH)e_BG=e_BG exactly.
        from tools.stability.ou3_theorem.matrix_certificates import add,matmul
        f=add(ops[1]['A'],matmul(ops[1]['K'],ops[1]['H']),-1)
        self.assertEqual(matmul(f,[[1],[0]]),[[1],[0]])

    def test_physical_span_sampling_and_construction_scope(self):
        self.assertEqual(sampled_tilt_span_lower(F(1,100),F(3,5),F(3,500)),F(7,2500))
        self.assertEqual(sampled_tilt_span_lower(F(1,1000),F(3,5),F(3,500)),0)
        g=gyro_construction_barrier()
        self.assertLess(g['first_prediction_increment_ceiling'],F(391,100000))
        self.assertGreater(g['complete_turn_requires_bias_norm_at_least'],1046)
        self.assertIsNone(g['all_time_signed_gain_innovation_sum_ceiling'])
        self.assertFalse(g['construction_unreachable_certified'])

    def test_signed_acceleration_supply_against_integrated_linear_history(self):
        # a(t)=2t, v=t^2 on [0,2]; C=(1,-1), h=1.
        # Left sampled action=-2, endpoint/variation action=-4+2=-2;
        # the two signed quadrature defects cancel exactly. The bound permits
        # their worst-case sum without discarding the physical velocity chain.
        c=[F(1),F(-1)]; v=[F(0),F(1),F(4)]
        telescoped=c[-1]*v[-1]-c[0]*v[0]+(c[0]-c[1])*v[1]
        self.assertEqual(telescoped,F(-2))
        self.assertEqual(signed_acceleration_supply(2,2,2,4,2),18)

    def test_interval_jerk_remainder_is_sharp_even_with_zero_total_weight(self):
        # C=(1,-1), a(t)=t, v=t^2/2 on [0,1]. D=0 cancels the
        # velocity term; the literal signed sample action is -1, not zero.
        self.assertEqual(physical_interval_sampling_error([1,1],[0,1],0,1,1),1)
        # A constant coefficient sampled at the midpoint has the sharp
        # Lipschitz mean-error factor H/4, not H/2.
        self.assertEqual(physical_interval_sampling_error([2],[F(1,2)],0,1,3),F(3,2))
        with self.assertRaises(ValueError):
            physical_interval_sampling_error([1],[2],0,1,1)

    def test_carried_word_witness_and_fail_closed_tamper(self):
        import json
        from copy import deepcopy
        from tools.stability.ou3_theorem.signed_temporal_diagnostic import verify_diagnostic
        p=ROOT/'reports/results/ou3_stability/signed-adjoint-diagnostic.json'
        report=json.loads(p.read_text())
        self.assertTrue(verify_diagnostic(report))
        bad=deepcopy(report); bad['kernel_witness'][0][0]='0'
        with self.assertRaises(ValueError): verify_diagnostic(bad)
        bad=deepcopy(report); bad['source_uniform_verified']=True
        with self.assertRaises(ValueError): verify_diagnostic(bad)

    def test_moving_forced_balance_stays_non_promoting(self):
        import json
        from copy import deepcopy
        from tools.stability.ou3_theorem.signed_temporal_balance import verify_diagnostic
        p=ROOT/'reports/results/ou3_stability/signed-balance-diagnostic.json'
        report=json.loads(p.read_text())
        self.assertTrue(verify_diagnostic(report))
        self.assertLess(F(report['S_interval_acceleration_margin_before_other_supplies']),0)
        self.assertLess(F(report['rotation_triangle_margin_before_other_supplies']),0)
        bad=deepcopy(report); bad['source_uniform_verified']=True
        with self.assertRaises(ValueError): verify_diagnostic(bad)
        bad=deepcopy(report); bad['root_term'][0]='1'
        with self.assertRaises(ValueError): verify_diagnostic(bad)

    def test_fail_closed(self):
        c=certificate()
        self.assertFalse(c["source_uniform_nominal_force_field_temporal_margin"])
        self.assertFalse(c["source_uniform_nominal_gyro_alias_temporal_margin"])
        self.assertFalse(c["uniform_historical_AG_readout_action"])
        self.assertFalse(c["full_21_covariance_upper"])
        self.assertFalse(c["rho0_certified"])
        self.assertFalse(c["temporal_margins_imply_finite_B_star"])
        self.assertFalse(c["zero_mean_projection_preserves_adjoint"])
        self.assertFalse(c["forced_data_adjoint_source_uniform_action_bound"])

if __name__=="__main__": unittest.main()
