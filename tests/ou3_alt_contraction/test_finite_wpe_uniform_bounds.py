from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_wpe_uniform_bounds as U
from tools.stability.ou3_alt_contraction import finite_wpe_machine_binary32 as W


def config(): return W.SHADOW.WPEConfig(U.LAMBDA,4,F(1,20),20,180)
def rounded_exp(x):
    lo,hi=U.exp_interval(F(x)); return U.B.rn32((lo+hi)/2)


class Tests(unittest.TestCase):
    def test_exact_induction_covers_every_alpha_and_nonproducing_branches(self):
        r=U.build(); self.assertEqual(U.validate(r),[])
        self.assertTrue(all(x>0 for x in r['induction_margins'].values()))
        self.assertGreater(r['log_induction_margin'],0)
        self.assertGreater(r['raw_period_derived_lower'],U.RAW_LO)
        self.assertLess(r['raw_period_derived_upper'],U.RAW_HI)
        self.assertTrue(r['no_positive_variance_or_period_production_assumed'])
        self.assertFalse(r['upstream_observer_totality_or_startup_reachability_closed'])
        self.assertFalse(r['target_libm_and_compiler_correspondence_closed'])
        r['state_bounds']['hp1']=F(1)
        self.assertTrue(U.validate(r))
        self.assertEqual(U.build()['state_bounds']['hp1'],U.HP1)

    def test_eager_frequency_exp_is_checked_even_before_usable_latch(self):
        from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as WF
        from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as ST
        s=W.initial(config(),bounded_profile=True)
        logs=replace(s.logs,separate=replace(s.logs.separate,log_period=F(1)),
                     fma=replace(s.logs.fma,log_period=F(1)))
        s=replace(s,logs=logs)
        exact=W.SHADOW.WPEState(log_period=F(1))
        good=WF.FrequencyExp(F(-1),rounded_exp(-1))
        kw=dict(logs=logs,fma_getter=good,shadow_frequency=None,
                stats_cfg=ST.StatsConfig(4,F(3,10),60,F(1,20),5),
                exact_min_hz=F(3,100),exact_max_hz=F(6,5))
        out=WF.machine_frequencies(exact,s,separate_getter=good,**kw)
        self.assertEqual(out[0].external.exact_shadow_frequency,WF.PRIOR_EXACT)
        legacy=WF.getters(WF.bind_log_state(exact,F(1)),period_exp=rounded_exp(1),frequency_exp=rounded_exp(-1))
        WF.machine_frequencies(exact,s,separate_getter=legacy,**kw)
        for bad in (WF.FrequencyExp(F(-1),None),WF.FrequencyExp(F(-1),F(1))):
            with self.assertRaisesRegex(ValueError,'exp witness detached'):
                WF.machine_frequencies(exact,s,separate_getter=bad,**kw)

    def test_transcendental_relations_reject_detached_and_nonfinite_witnesses(self):
        for x in (F(-48),F(-1),F(0),F(1),F(48)):
            U.check_exp(x,rounded_exp(x))
            with self.assertRaisesRegex(ValueError,'exp witness detached'):
                U.check_exp(x,None)
        for x in (U.RAW_LO,F(1,2),F(1),F(2),U.RAW_HI):
            lo,hi=U.log_interval(x); y=U.B.rn32((lo+hi)/2)
            U.check_log(x,y)
            with self.assertRaisesRegex(ValueError,'log witness detached'):
                U.check_log(x,U.B.add(y,1))

    def test_first_valid_log_really_comes_from_machine_moments(self):
        # Component induction predecessor, not asserted to be a reachable
        # physical startup fixture. All bounded predecessors must be covered.
        from test_finite_admitted_wpe_machine_clock_interleaved_prefix import general_witness
        s=W.initial(config(),bounded_profile=True)
        moment=W.MOM.State(elapsed=U.B.rn32(30),weight=F(1),velocity_sq=F(9),elevation_sq=F(1))
        s=replace(s,separate=moment,fma=moment)
        ws=[]
        for mode in ('separate','fma'):
            w=general_witness(s,F(0),U.DT,mode)
            raw=w.raw_log.raw_period; lo,hi=U.log_interval(raw)
            lr=U.B.rn32((lo+hi)/2)
            w=replace(w,raw_log=W.RawLogBinding(raw,lr),log=W.LOG.InitWitness(lr),
                      usable=W.USABLE.PeriodWitness(lr,rounded_exp(lr)))
            ws.append(w)
        out,*_=W.step(s,dt=U.DT,vertical_accel=0,separate=ws[0],fma=ws[1])
        self.assertTrue(out.bounded_profile)
        U.check_state(out)
        bad=replace(ws[0],raw_log=replace(ws[0].raw_log,log_raw=F(0)),log=W.LOG.InitWitness(F(0)))
        with self.assertRaisesRegex(ValueError,'log witness detached'):
            W.step(s,dt=U.DT,vertical_accel=0,separate=bad,fma=ws[1])

    def test_bounds_are_required_by_the_qualified_runtime_state(self):
        s=W.initial(config(),bounded_profile=True)
        with self.assertRaisesRegex(ValueError,'source input exceeds'):
            W.step(s,dt=U.DT,vertical_accel=F(33),separate=W.ModeWitnesses({}),fma=W.ModeWitnesses({}))
        with self.assertRaisesRegex(ValueError,'uniform hp1 bound'):
            replace(s,separate=replace(s.separate,hp1=2*U.HP1))
        with self.assertRaisesRegex(ValueError,'default construction constants'):
            W.initial(replace(config(),min_horizon_sec=F(1)),bounded_profile=True)

    def test_scaled_envelopes_prove_margins_and_explicit_totality_limit(self):
        for exponent in range(6):
            r=U.build(exponent)
            self.assertEqual(U.validate(r),[])
            self.assertTrue(all(m>0 for m in r['induction_margins'].values()))
            self.assertGreater(r['log_induction_margin'],0)
            self.assertGreater(r['finite_binary32_intermediate_margin'],0)
            self.assertEqual(r['state_bounds']['velocity'],U.V*2**exponent)
            self.assertEqual(r['state_bounds']['velocity_sq'],2*U.V**2*4**exponent)
            lo=r['raw_period_interval'][0]
            l,h=U.log_interval(lo,exponent)
            U.check_log(lo,U.B.rn32((l+h)/2),exponent)
        with self.assertRaisesRegex(ValueError,'cannot prove binary32 totality'):
            U.build(6)
        with self.assertRaisesRegex(ValueError,'outside proved period bounds'):
            U.check_log(U.RAW_LO/2,F(-40))

    def test_widening_preserves_nontrivial_machine_history_and_advances(self):
        from test_finite_admitted_wpe_machine_clock_interleaved_prefix import general_witness
        before=W.initial(config(),bounded_profile=True)
        witnesses=[general_witness(before,F(16),U.DT,mode) for mode in ('separate','fma')]
        before,*_=W.step(before,dt=U.DT,vertical_accel=16,separate=witnesses[0],fma=witnesses[1])
        self.assertIs(W.widen_supply(before,input_abs_upper=16),before)
        wide=W.widen_supply(before,input_abs_upper=512)
        self.assertTrue(W.same_machine_history(before,wide))
        for field in ('cfg','separate','fma','logs'):
            self.assertIs(getattr(wide,field),getattr(before,field))
        self.assertFalse(W.same_machine_history(wide,before))
        detached=replace(wide,logs=replace(wide.logs,separate=replace(wide.logs.separate,log_period=F(1))))
        self.assertFalse(W.same_machine_history(before,detached))
        witnesses=[general_witness(wide,F(512),U.DT,mode) for mode in ('separate','fma')]
        out,*_=W.step(wide,dt=U.DT,vertical_accel=512,separate=witnesses[0],fma=witnesses[1])
        self.assertEqual(out.supply_scale_exponent,4)
        self.assertEqual(out.separate.samples,before.separate.samples+1)
        U.check_state(out)

    def test_live_iss_uses_its_own_norm_instead_of_startup_sensor_profile(self):
        small=U.live_iss_certificate(0); larger=U.live_iss_certificate(500)
        self.assertEqual(small['scale_exponent'],0)
        self.assertEqual(larger['scale_exponent'],5)
        self.assertEqual(larger['raw_accel_residual_norm_upper'],750)
        self.assertFalse(larger['startup_sensor_residual_bound_reused_for_Live'])
        self.assertFalse(larger['unrestricted_Live_ISS_amplitude_domain_closed'])
        self.assertTrue(larger['observer_totality_assumed_not_proved'])
        larger['scale_exponent']=0
        self.assertEqual(U.live_iss_certificate(500)['scale_exponent'],5)
        with self.assertRaisesRegex(ValueError,'cannot prove binary32 totality'):
            U.live_iss_certificate(1000)

    def test_error_inclusive_libm_relation_has_strict_positive_uniform_margins(self):
        for exponent in (0,4,5):
            r=U.build(exponent,U.ERROR_PROFILE)
            self.assertEqual(U.validate(r),[])
            self.assertTrue(all(m>0 for m in r['induction_margins'].values()))
            self.assertGreater(r['log_induction_margin'],0)
            self.assertGreater(r['finite_binary32_intermediate_margin'],0)
            self.assertFalse(r['target_libm_and_compiler_correspondence_closed'])
        for x in (-48,-1,0,1,48):
            exact=rounded_exp(x)
            shifted=U.B.rn32(exact*(1+U.EXP_SQRT_REL_ERROR/2))
            U.check_exp(x,shifted,U.ERROR_PROFILE)
            with self.assertRaisesRegex(ValueError,'exp witness detached'):
                U.check_exp(x,shifted)
        raw=U.RAW_LO/16
        lo,hi=U.log_interval(raw,4); nominal=U.B.rn32((lo+hi)/2)
        shifted=U.B.add(nominal,U.LOG_ABS_ERROR/2)
        U.check_log(raw,shifted,4,U.ERROR_PROFILE)
        with self.assertRaisesRegex(ValueError,'log witness detached'):
            U.check_log(raw,shifted,4)
        root=U.M.SQRT.sqrt32(F(2)); shifted=U.B.rn32(root*(1+U.EXP_SQRT_REL_ERROR/2))
        U.check_sqrt(2,shifted,U.ERROR_PROFILE)
        with self.assertRaisesRegex(ValueError,'sqrt result detached'):
            U.check_sqrt(2,shifted)

    def test_error_profile_is_carried_into_actual_moment_and_log_product(self):
        state=W.initial(config(),bounded_profile=True,libm_profile=U.ERROR_PROFILE)
        x=F(32); d=U.exp_result_interval(-U.B.mul(U.LAMBDA,U.DT),U.ERROR_PROFILE)[0]
        gain=U.B.div(1-d,U.LAMBDA)
        hp1=U.B.mul(d,x); hp2=U.B.mul(d,hp1)
        vel=U.B.mul(gain,hp2); elev=U.B.mul(gain,vel)
        witness=W.ModeWitnesses(dict(decay_exp=d,velocity_successor=vel,elevation_successor=elev))
        out,*_=W.step(state,dt=U.DT,vertical_accel=x,separate=witness,fma=witness)
        self.assertEqual(out.libm_profile,U.ERROR_PROFILE)
        self.assertEqual(out.separate.samples,1)
        with self.assertRaisesRegex(ValueError,'exp argument RNE cell'):
            W.step(replace(state,libm_profile=U.RNE_PROFILE),dt=U.DT,vertical_accel=x,
                   separate=witness,fma=witness)

    def test_MEMS_domain_produces_uniform_512_supply_without_ISS_amplitude_limit(self):
        for profile in (U.RNE_PROFILE,U.ERROR_PROFILE):
            r=U.physical_input_certificate(profile)
            self.assertLess(r['vertical_abs_upper'],322)
            self.assertEqual(r['vertical_input_abs_upper'],512)
            self.assertEqual(r['scale_exponent'],4)
            self.assertEqual(r['raw_accel_component_upper'],160)
            self.assertEqual(r['raw_gyro_component_upper'],35)
            self.assertTrue(r['supply_does_not_require_a_numeric_ISS_W_limit'])
            self.assertTrue(r['source_uniform_WPE_supplies_under_declared_MEMS_prefix_closed'])
            self.assertTrue(r['all_initialized_finite_prefixes_totality_closed'])
            self.assertFalse(r['startup_deadline_required_for_WPE_supply_bound'])
            self.assertFalse(r['startup_seed_capture_and_target_qualification_proved_here'])

    def test_actual_target_exp_log_proofs_fit_same_WPE_argument_and_error_domains(self):
        r=U.target_error_profile_certificate()
        self.assertTrue(r['actual_WPE_arguments_covered_by_target_library_proof'])
        self.assertTrue(r['target_exp_log_satisfy_WPE_error_profile_under_scalar_and_link_premises'])
        self.assertLess(r['exp_relative_error_upper'],U.EXP_SQRT_REL_ERROR)
        self.assertLess(r['log_absolute_error_upper'],U.LOG_ABS_ERROR)
        self.assertFalse(r['transcendental_correct_rounding_assumed'])
        self.assertFalse(r['scalar_sqrt_ROM_and_firmware_link_qualification_supplied_here'])

    def test_per_mode_fused_omega_flows_through_moment_raw_log_and_usable_product(self):
        from test_finite_admitted_wpe_machine_clock_interleaved_prefix import general_witness
        state=W.initial(config(),bounded_profile=True,libm_profile=U.ERROR_PROFILE)
        near=U.B.add(U.B.mul(U.LAMBDA,U.LAMBDA),F(8,2**29))
        moment=W.MOM.State(elapsed=U.B.rn32(30),weight=F(1),velocity_sq=near,elevation_sq=F(1))
        state=replace(state,separate=moment,fma=moment)
        witnesses=[]; omegas=[]
        for mode,which in (('separate',-1),('fma',0)):
            original=general_witness(state,F(0),U.DT,mode)
            args=dict(original.moment)
            ratio=U.B.div(args['velocity_var_successor'],args['elevation_var_successor'])
            omega=W.MOM.omega_sq_choices(ratio,U.LAMBDA)[which]
            self.assertGreater(omega,W.MOM.OMEGA_GATE)
            root=W.MOM.SQRT.sqrt32(omega); raw=U.B.div(W.MOM.TWO_PI,root)
            lo,hi=U.log_interval(raw); log=U.B.rn32((lo+hi)/2)
            args.update(omega_sq_successor=omega,sqrt_omega=root)
            witnesses.append(W.ModeWitnesses(args,raw_log=W.RawLogBinding(raw,log),
                log=W.LOG.InitWitness(log),usable=W.USABLE.PeriodWitness(log,rounded_exp(log))))
            omegas.append(omega)
        out,separate,fma,_=W.step(state,dt=U.DT,vertical_accel=0,
                                separate=witnesses[0],fma=witnesses[1])
        self.assertNotEqual(omegas[0],omegas[1])
        self.assertEqual((separate.omega_sq,fma.omega_sq),tuple(omegas))
        self.assertNotEqual(separate.raw_period,fma.raw_period)
        self.assertEqual(out.logs.separate.log_period,witnesses[0].raw_log.log_raw)
        self.assertEqual(out.logs.fma.log_period,witnesses[1].raw_log.log_raw)


if __name__=='__main__': unittest.main()
