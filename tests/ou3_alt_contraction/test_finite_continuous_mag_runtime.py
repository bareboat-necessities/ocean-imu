"""Exact calibration algebra/guard regressions, not admitted BRMM trajectories."""
from dataclasses import replace
from fractions import Fraction as F
import unittest
from tools.stability.ou3_alt_contraction import finite_continuous_mag_runtime as X

I=X.IDENTITY
ROTATIONS=(I, ((1,0,0),(0,-1,0),(0,0,-1)),
           ((-1,0,0),(0,1,0),(0,0,-1)), ((-1,0,0),(0,-1,0),(0,0,1)))


def fit(extra_identity=False):
    rotations=ROTATIONS+((I,) if extra_identity else ())
    n=len(rotations)
    eigen=F(24,25) if extra_identity else F(1)
    cfg=X.Config(memory=0,solve_period=n,min_weight=n,min_information=1,
                 ridge=eigen,relative_ridge=0,max_residual=3)
    s=X.State()
    for j,rot in enumerate(rotations):
        raw=X.add(X.mv(X.transpose(rot),(30,0,0)),(2,0,0))
        kw=dict(eigen_success=True,spectrum=X.Spectrum(I,(eigen,)*3)) if j==n-1 else {}
        out=X.update(s,cfg,dt=1,rotation=rot,raw_body=raw,**kw)
        s=out.state
    return s,cfg


class Tests(unittest.TestCase):
    def test_exact_unridged_matrix_regularized_fit_and_residual(self):
        s,cfg=fit()
        self.assertTrue(s.estimate.valid)
        self.assertEqual(s.estimate.bias,(1,0,0))
        self.assertEqual(s.estimate.field,(30,0,0))
        self.assertEqual(s.estimate.residual_squared,1)
        self.assertEqual(s.elapsed,0)
        self.assertEqual(s.estimate.information,4)

    def test_nonzero_mean_rotation_kept_in_both_fit_and_reference(self):
        s,_=fit(True)
        self.assertEqual(s.estimate.bias,(1,0,0))
        self.assertEqual(s.estimate.field,(F(151,5),0,0))
        self.assertEqual(X.level_reference(s,(0,0,0),F(1,1000)),(F(152,5),0,0))
        self.assertEqual(X.level_reference(s,(1,0,0),F(1,1000)),s.estimate.field)

    def test_raw_not_precompensated_packet_enters_all_sufficient_statistics(self):
        cfg=X.Config(memory=0,solve_period=10)
        s=X.update(X.State(),cfg,dt=F(1,100),rotation=ROTATIONS[1],raw_body=(3,4,5)).state
        self.assertEqual(s.body_sum,(3,4,5)); self.assertEqual(s.level_sum,(3,-4,-5))
        self.assertEqual(s.square_sum,50); self.assertEqual(s.rotation_sum,ROTATIONS[1])

    def test_same_decay_scales_every_statistic(self):
        cfg=X.Config(memory=2,solve_period=10)
        s=X.update(X.State(),cfg,dt=1,rotation=I,raw_body=(30,0,0),decay=X.Decay(1,2,F(1,2))).state
        s=X.update(s,cfg,dt=1,rotation=ROTATIONS[1],raw_body=(0,30,0),decay=X.Decay(1,2,F(1,2))).state
        self.assertEqual(s.weight,F(3,2)); self.assertEqual(s.square_sum,F(1350))
        self.assertEqual(s.body_sum,(15,30,0)); self.assertEqual(s.level_sum,(15,-30,0))
        self.assertEqual(s.rotation_sum,((F(3,2),0,0),(0,F(-1,2),0),(0,0,F(-1,2))))

    def test_decay_cannot_detach_clock_or_memory(self):
        with self.assertRaisesRegex(ValueError,'detached'):
            X.update(X.State(),X.Config(),dt=1,rotation=I,raw_body=(30,0,0),decay=X.Decay(2,600,1))

    def test_rejected_sample_does_not_advance_solve_elapsed(self):
        s=replace(X.State(),elapsed=F(1,2))
        out=X.update(s,X.Config(),dt=1,rotation=None,raw_body=(0,0,0))
        self.assertIs(out.state,s); self.assertFalse(out.sample_accepted)
        with self.assertRaisesRegex(ValueError,'consumes no'):
            X.update(s,X.Config(),dt=1,rotation=None,raw_body=(0,0,0),decay=X.Decay(1,600,1))

    def test_not_due_consumes_no_eigenpair(self):
        with self.assertRaisesRegex(ValueError,'not-due'):
            X.update(X.State(),X.Config(memory=0),dt=F(1,100),rotation=I,raw_body=(30,0,0),eigen_success=False)

    def test_underweight_due_replaces_stale_valid_estimate(self):
        s,cfg=fit(True)
        cfg=replace(cfg,memory=600)
        out=X.update(s,cfg,dt=5,rotation=I,raw_body=(30,0,0),decay=X.Decay(5,600,F(1,100)))
        self.assertEqual(out.branch,'weight'); self.assertFalse(out.state.estimate.valid)
        self.assertEqual(out.state.estimate.bias,(0,0,0)); self.assertEqual(out.state.elapsed,0)

    def test_eigensolver_rejection_and_detached_spectrum(self):
        cfg=X.Config(memory=0,min_weight=1,min_information=0)
        out=X.update(X.State(),cfg,dt=1,rotation=I,raw_body=(30,0,0),eigen_success=False)
        self.assertEqual(out.branch,'eigensolver_rejected')
        with self.assertRaisesRegex(ValueError,'same sufficient statistics'):
            X.update(X.State(),cfg,dt=1,rotation=I,raw_body=(30,0,0),eigen_success=True,
                     spectrum=X.Spectrum(I,(1,1,1)))

    def test_information_gate_not_inferred_from_regularization(self):
        cfg=X.Config(memory=0,min_weight=1)
        out=X.update(X.State(),cfg,dt=1,rotation=I,raw_body=(30,0,0),eigen_success=True,
                     spectrum=X.Spectrum(I,(0,0,0)))
        self.assertEqual(out.branch,'information'); self.assertFalse(out.state.estimate.valid)

    def test_fit_bias_gate_and_residual_gate(self):
        for field,bias,max_fraction,max_residual,branch in (
            ((0,0,30),(30,0,0),F(1,10),100,'bias_norm'),
            ((0,0,30),(2,0,0),1,F(1,2),'residual')):
            cfg=X.Config(memory=0,solve_period=4,min_weight=4,min_information=1,
                         ridge=1,relative_ridge=0,max_bias_fraction=max_fraction,max_residual=max_residual)
            s=X.State()
            for j,rot in enumerate(ROTATIONS):
                raw=X.add(X.mv(X.transpose(rot),field),bias)
                kw=dict(eigen_success=True,spectrum=X.Spectrum(I,(1,1,1))) if j==3 else {}
                out=X.update(s,cfg,dt=1,rotation=rot,raw_body=raw,**kw); s=out.state
            self.assertEqual(out.branch,branch); self.assertFalse(s.estimate.valid)

    def test_refinement_blocks_application_even_with_valid_fit(self):
        s,_=fit(True); a=X.Applied()
        out=X.apply(a,s,(30,0,0),time=90,refinement_done=False)
        self.assertEqual(out.state,a); self.assertFalse(out.wrote_reference)
        self.assertIsNone(out.state.last_time); self.assertIsNone(out.state.anchor_reference)

    def test_reference_and_offset_change_together_from_same_current_statistics(self):
        s,_=fit(True)
        out=X.apply(X.Applied(),s,(30,0,0),time=120,refinement_done=True,slew_tau=0,
                    new_norm=X.HorizontalNorm((F(151,5),0),F(151,5)),
                    anchor_norm=X.HorizontalNorm((F(152,5),0),F(152,5)))
        self.assertTrue(out.wrote_reference); self.assertEqual(out.state.applied,(1,0,0))
        self.assertEqual(out.reference,(F(149,5),0,0))
        self.assertEqual(out.state.anchor_reference,(30,0,0)); self.assertEqual(out.state.anchor_bias,(0,0,0))

    def test_same_timestamp_apply_uses_fallback_not_zero_dt(self):
        s,_=fit()
        a=X.Applied(last_time=120)
        out=X.apply(a,s,(30,0,0),time=120,refinement_done=True,
                    decay=X.Decay(F(1,200),45,F(1,2)),
                    new_norm=X.HorizontalNorm((30,0),30),anchor_norm=X.HorizontalNorm((30,0),30))
        self.assertEqual(out.state.applied,(F(1,2),0,0))
        with self.assertRaisesRegex(ValueError,'detached'):
            X.apply(a,s,(30,0,0),time=120,refinement_done=True,
                    decay=X.Decay(1,45,F(1,2)))

    def test_anchor_and_clock_survive_level_reference_failure(self):
        s,_=fit(True)
        out=X.apply(X.Applied(),s,(30,0,0),time=120,refinement_done=True,
                    fraction=152,slew_tau=0)
        self.assertEqual(out.branch,'level_reference'); self.assertFalse(out.wrote_reference)
        self.assertEqual(out.state.applied,(0,0,0)); self.assertEqual(out.state.last_time,120)
        self.assertEqual(out.state.anchor_reference,(30,0,0))

    def test_anchor_and_clock_survive_horizontal_reference_failure(self):
        s,_=fit(True)
        out=X.apply(X.Applied(),s,(F(1,1000),0,30),time=120,refinement_done=True,slew_tau=0,
                    new_norm=X.HorizontalNorm((F(151,5),0),F(151,5)),
                    anchor_norm=X.HorizontalNorm((F(152,5),0),F(152,5)))
        self.assertEqual(out.branch,'horizontal_reference'); self.assertEqual(out.state.applied,(0,0,0))
        self.assertEqual(out.state.last_time,120); self.assertIsNotNone(out.state.anchor_bias)

    def test_detached_horizontal_norm_rejected(self):
        s,_=fit(True)
        with self.assertRaisesRegex(ValueError,'detached'):
            X.apply(X.Applied(),s,(30,0,0),time=120,refinement_done=True,slew_tau=0,
                    new_norm=X.HorizontalNorm((30,0),30),anchor_norm=X.HorizontalNorm((30,0),30))

    def test_deterministic_bound_constants_and_no_theorem_promotion(self):
        b=X.deterministic_bounds()
        self.assertEqual(b['raw_norm_uT'],82)
        self.assertEqual(b['accepted_fitted_bias_norm_uT'],F(287,10))
        self.assertEqual(b['active_reference_norm_uT'],F(1107,10))
        self.assertEqual(b['physical_to_active_residual_norm_uT'],F(1107,5))
        r=X.readiness()
        for name in ('ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS',
                     'complete_word_finite_identity','deployment_exp_cast_solver_nonfinite_branches_closed'):
            self.assertFalse(r[name])

    def test_moment_matrix_is_weighted_rotation_variance_not_independent_box(self):
        s,_=fit(True)
        rotations=ROTATIONS+(I,)
        a=tuple(X.scale(row,1/s.weight) for row in s.rotation_sum)
        normal=tuple(tuple(I[i][j]-X.mm(X.transpose(a),a)[i][j] for j in range(3)) for i in range(3))
        variance=[[F(0)]*3 for _ in range(3)]
        for rot in rotations:
            diff=tuple(X.sub(row,avg) for row,avg in zip(rot,a))
            term=X.mm(X.transpose(diff),diff)
            for i in range(3):
                for j in range(3): variance[i][j]+=term[i][j]/s.weight
        self.assertEqual(tuple(map(tuple,variance)),normal)


if __name__=='__main__': unittest.main()
