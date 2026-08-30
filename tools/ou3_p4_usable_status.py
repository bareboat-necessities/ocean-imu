#!/usr/bin/env python3
"""Machine-readable status for the usable OU-III P4 replacement route.

This consumer intentionally separates a successful *partial* full-word result
from a usable nonlinear P4 certificate.  It consumes already-produced JSON so
CI does not rerun the expensive outward-rounded translation calculation.

A useful P4 may be promoted only after the complete source-reachable fixed-mode
word, including attitude/gyro-bias (and active accelerometer bias in A mode), is
validated and the exact finite-angle return map closes on a non-microscopic
set.  The old scalar B*W frontier is never used for promotion here.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path


def build_from_payloads(translation: dict, design: dict, path: dict) -> dict:
    failures=[]; modes={}
    for mode in ('H','A'):
        t=translation.get('modes',{}).get(mode,{})
        delta=float(t.get('complete_word_translation_margin_lower',0.0))
        old=float(t.get('old_single_seed_translation_margin_lower',0.0))
        factor=float(t.get('margin_widening_factor_lower',0.0))
        horizons=design.get('modes',{}).get(mode,{}).get('horizons',{})
        d8=float(horizons.get('8.0',{}).get('worst_grid_point',{}).get('translation_complete_word_generalized_margin_design',0.0))
        if not (delta>0 and old>0 and factor>1 and d8>0): failures.append(f'{mode}: incomplete full-word evidence')
        modes[mode]={
            'validated_one_second_translation_margin_lower':delta,
            'old_single_seed_translation_margin_lower':old,
            'validated_translation_widening_factor_lower':factor,
            'eight_second_design_translation_margin':d8,
            'translation_complete_word_progress':'PASS' if delta>old else 'NOT_ESTABLISHED',
            'full_state_complete_word_validated':False,
            'exact_finite_angle_complete_return_map_validated':False,
            'usable_certificate_status':'NOT_ESTABLISHED',
        }
    recurrent=bool(path.get('old_worst_corner_has_internal_recurrent_cycle',False))
    if path.get('path_graph_ready') is not True: failures.append('source path graph not ready')
    return {
        'qualification':'OU3_P4_USABLE_CERTIFICATE_STATUS',
        'source_only':True,
        'trajectory_replay_used_for_promotion':False,
        'old_scalar_frontier_used_for_promotion':False,
        'source_path_graph_ready':path.get('path_graph_ready') is True,
        'old_worst_corner_has_internal_recurrent_cycle':recurrent,
        'translation_complete_word_validated_on_old_worst_cell':translation.get('P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS')=='PASS',
        'meaningful_linear_improvement_established':all(modes[m]['validated_translation_widening_factor_lower']>1e5 for m in modes),
        'P4_USABLE_CERTIFICATE_STATUS':'NOT_ESTABLISHED',
        'next_blocker':(
            'validate full-state complete-word Phi/Omega on the recurrent worst source cell and then validate the exact finite-angle complete return map; '
            'the source-path graph cannot eliminate this corner because it contains an internal recurrent cycle'
        ),
        'modes':modes,
        'failures':failures,
    }


def validate(d: dict) -> list[str]:
    f=list(d.get('failures',[]))
    if d.get('source_only') is not True or d.get('trajectory_replay_used_for_promotion') is not False: f.append('status is not source-only')
    if d.get('old_scalar_frontier_used_for_promotion') is not False: f.append('retired scalar frontier was promoted')
    if d.get('translation_complete_word_validated_on_old_worst_cell') is not True: f.append('validated translation progress missing')
    if d.get('meaningful_linear_improvement_established') is not True: f.append('no meaningful validated linear improvement')
    if d.get('P4_USABLE_CERTIFICATE_STATUS')!='NOT_ESTABLISHED': f.append('partial evidence prematurely promoted usable P4')
    for mode in ('H','A'):
        m=d.get('modes',{}).get(mode,{})
        if m.get('translation_complete_word_progress')!='PASS': f.append(f'{mode}: translation progress did not pass')
        if m.get('full_state_complete_word_validated') is not False: f.append(f'{mode}: full-state word incorrectly claimed')
        if m.get('exact_finite_angle_complete_return_map_validated') is not False: f.append(f'{mode}: nonlinear return map incorrectly claimed')
    return f


def main() -> int:
    ap=argparse.ArgumentParser();ap.add_argument('--translation',type=Path,required=True);ap.add_argument('--design',type=Path,required=True);ap.add_argument('--path',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build_from_payloads(json.loads(a.translation.read_text()),json.loads(a.design.read_text()),json.loads(a.path.read_text()));f=validate(d);d['validation_failures']=f;a.output.write_text(json.dumps(d,indent=2,sort_keys=True))
    print(json.dumps({'status':d['P4_USABLE_CERTIFICATE_STATUS'],'meaningful_linear_improvement':d['meaningful_linear_improvement_established'],'worst_corner_recurrent':d['old_worst_corner_has_internal_recurrent_cycle'],'modes':{m:{'validated_delta':d['modes'][m]['validated_one_second_translation_margin_lower'],'factor':d['modes'][m]['validated_translation_widening_factor_lower'],'design_8s_delta':d['modes'][m]['eight_second_design_translation_margin']} for m in ('H','A')},'next_blocker':d['next_blocker'],'failures':f},indent=2))
    return 0 if not f else 2
if __name__=='__main__':raise SystemExit(main())
