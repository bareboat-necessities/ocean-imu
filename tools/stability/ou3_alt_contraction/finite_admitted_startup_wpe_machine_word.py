"""Admitted startup/goLive constructor retaining the full WPE machine history.

This closes a composition gap between the startup-rooted full binary32 WPE
moment/log product and the admitted COMPLETE-BRMM + BIAS + disturbance Live
word.  No Live WPE snapshot is accepted: the only WPE state is the one advanced
from literal startup reset and preserved by the same goLive control edge.

This is ancestry/attachment only.  Universal startup reachability, target-libm
correspondence, source-uniform arithmetic supplies and the literal 600-edge
finite master remain fail-closed.
"""
from __future__ import annotations

from tools.stability.ou3_alt_contraction import finite_startup_wpe_machine_history as START
from tools.stability.ou3_alt_contraction import finite_admitted_wpe_machine_clock_interleaved_prefix as LIVE
from tools.stability.ou3_alt_contraction import admitted_startup_mahony_invariant as MAHONY

QUALIFICATION='OU3_ALT_ADMITTED_STARTUP_WPE_MACHINE_WORD_V1'


def build(go:START.GoLive,magnetic,origin,bias_history,**source_witnesses):
    if not isinstance(go,START.GoLive):
        raise TypeError('startup-rooted full-WPE goLive result required')
    out=LIVE.begin_from_startup(go,magnetic,origin,bias_history,**source_witnesses)
    if out.wpe is not go.wpe:
        raise ValueError('admitted Live word reseeded full WPE machine history')
    if out.entry_wpe_samples!=go.wpe.separate.samples or out.entry_wpe_samples!=go.wpe.fma.samples:
        raise ValueError('admitted Live WPE entry count detached from startup history')
    if out.wpe_steps!=0 or out.base.clock_steps!=0:
        raise ValueError('admitted constructor executed a hidden Live IMU event')
    return out


def readiness():
    s=START.readiness(); l=LIVE.readiness(); mah=MAHONY.build()
    mf=MAHONY.validate(mah)
    if mf: raise RuntimeError('admitted startup Mahony invariant invalid: '+repr(mf))
    return {
      'literal_reset_full_WPE_machine_history_attached_through_startup':s['full_WPE_machine_histories_rooted_with_literal_startup_reset'],
      'same_private_Mahony_vertical_drives_startup_WPE_machine':s['same_machine_Mahony_vertical_output_drives_full_WPE_moment_product'],
      'full_WPE_machine_history_preserved_by_same_goLive_edge':s['goLive_preserves_full_WPE_moment_log_machine_history_by_identity'],
      'admitted_Live_constructor_accepts_only_startup_produced_WPE_state':l['startup_produced_full_WPE_state_constructor_available'],
      'admitted_source_private_Mahony_startup_to_Live_invariant_closed':mah['admitted_source_private_Mahony_startup_to_Live_invariant_closed'],
      'startup_guard_Mahony_WPE_machine_to_admitted_Live_attachment_closed':True,
      'every_admitted_startup_history_reaches_this_product':False,
      'source_uniform_startup_WPE_supply_bounds_closed':False,
      'target_WPE_libm_and_compiler_profile_correspondence_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,'ALT_END_TO_END_PASS':False,
    }
