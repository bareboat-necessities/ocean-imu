"""Coherent dual-compiler WPE-log -> frequency -> tau Live regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_tau_interleaved_prefix as TAUJOIN
from tools.stability.ou3_alt_contraction import finite_admitted_wpe_tau_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as LEDGER
from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as WPEF
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as WPELOG
import test_finite_admitted_source_imu_word as IBASE
import test_finite_admitted_tau_interleaved_prefix as TBASE
import test_finite_source_bound_live_word as LBASE
import test_finite_source_bound_prediction_word as PBASE


def base_state(*,usable=False):
    p=TBASE.usable_half_prefix() if usable else TBASE.qualified_prefix()
    return TAUJOIN.begin(p,LEDGER.State(updates=17))


def machine_wpe_for(base):
    exact=TAUJOIN._entry_wpe(base)
    if exact.log_period is None: return WPELOG.initial()
    q=B.rn32(exact.log_period); t=WPELOG.Track(q,1)
    return WPELOG.State(t,t,0)


def half_getters(state):
    exact=TAUJOIN._entry_wpe(state.base)
    sl=WPEF.StoredLogPeriod(F(exact.log_period),state.wpe.separate.log_period,
                            state.wpe.separate.log_period-F(exact.log_period))
    fl=WPEF.StoredLogPeriod(F(exact.log_period),state.wpe.fma.log_period,
                            state.wpe.fma.log_period-F(exact.log_period))
    return (WPEF.getters(sl,period_exp=B.rn32(2),frequency_exp=B.rn32(F(1,2))),
            WPEF.getters(fl,period_exp=B.rn32(2),frequency_exp=B.rn32(F(1,2))))


class Tests(unittest.TestCase):
    def test_usable_entry_orders_frequency_tau_then_current_WPE_update(self):
        # Component-only rational-root cell.  The carried WPE state owns the
        # output; these preupdate operands are the exp witnesses consumed by it.
        b=base_state(usable=True); s=X.begin(b,machine_wpe_for(b))
        witness,segment,raw,r,br,dynamic=IBASE.operands(s.base.prefix.prefix.live)
        dynamic=dict(dynamic)
        dynamic['preupdate_period']=F(2); dynamic['preupdate_frequency']=F(1,2)
        dynamic['spectral']=CAND.SpectralWitness(1,1)
        sg,fg=half_getters(s)
        e=TBASE.machine_decay(F(1,2),s.base.prefix.prefix.live.live_word.runtime.candidate_cfg)
        before_wpe=s.wpe
        out=X.imu_step(s,separate_getter=sg,fma_getter=fg,shadow_frequency=F(1,2),
            separate_tau_exp_decay=e,fma_tau_exp_decay=e,
            restricted=r,bias_restricted=br,witness=witness,raw=raw,
            packet_id='imu-wpe-tau',**PBASE.root_args(),**dynamic)
        self.assertEqual(out.state.base.tau.updates,18)
        self.assertEqual(out.state.wpe.samples,before_wpe.samples+1)
        self.assertEqual(out.tau_step.separate_step.frequency,out.tau_step.fma_step.frequency)
        self.assertEqual(out.separate_supply.frequency,0); self.assertEqual(out.fma_supply.frequency,0)
        self.assertFalse(out.wpe_step.produced_period)

    def test_preusable_frequency_source_is_prior_without_claiming_full_spectral_step(self):
        b=base_state(); s=X.begin(b,machine_wpe_for(b)); exact=TAUJOIN._entry_wpe(s.base)
        self.assertFalse(exact.usable_period)
        sf,ff=X._frequency_sources(s,separate_getter=None,fma_getter=None,shadow_frequency=None)
        self.assertEqual(sf.external.branch,'prior'); self.assertEqual(ff.external.branch,'prior')
        self.assertEqual(sf.stored.stored_hz,WPEF.PRIOR); self.assertEqual(ff.stored.stored_hz,WPEF.PRIOR)
        # Full deployment SpectralMSE now uses an interval exact-root relation;
        # this legacy exact-rational component fixture intentionally does not.

    def test_mag_and_hold_preserve_both_machine_ledgers(self):
        b=base_state(); s=X.begin(b,machine_wpe_for(b)); w=s.wpe; t=s.base.tau
        s,_=X.mag_step(s,**LBASE.mag_kwargs(s.base.prefix.prefix.live.live_word))
        self.assertIs(s.wpe,w); self.assertIs(s.base.tau,t)
        s,_=X.set_hold(s,hold=False)
        self.assertIs(s.wpe,w); self.assertIs(s.base.tau,t)

    def test_machine_WPE_initialization_must_match_exact_carried_state(self):
        b=base_state(); exact=TAUJOIN._entry_wpe(b)
        if exact.log_period is None:
            bad=WPELOG.State(WPELOG.Track(B.rn32(F(1,2)),1),WPELOG.Track(B.rn32(F(1,2)),1),0)
        else:
            bad=WPELOG.initial()
        with self.assertRaisesRegex(ValueError,'initialization detached'):
            X.begin(b,bad)

    def test_readiness_closes_global_history_topology_not_numerics(self):
        r=X.readiness()
        for k in ('global_separate_and_FMA_WPE_log_tracks_carried_in_Live_product',
                  'sample_entry_WPE_frequency_precedes_current_sample_WPE_update',
                  'separate_WPE_frequency_advances_only_separate_tau_track',
                  'FMA_WPE_frequency_advances_only_FMA_tau_track',
                  'global_compiler_track_coherence_includes_WPE_log_and_tau_states',
                  'exact_vs_machine_frequency_target_decay_supplies_retained_per_mode',
                  'MAG_and_HOLD_preserve_both_WPE_and_tau_ledgers'):
            self.assertTrue(r[k])
        self.assertFalse(r['startup_WPE_machine_ledger_provenance_closed'])
        self.assertFalse(r['WPE_log_std_log_target_libm_correspondence_closed'])
        self.assertFalse(r['source_uniform_WPE_frequency_supply_bound_closed'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])


if __name__=='__main__': unittest.main()
