"""Point diagnostic: qualified fresh-Live entry through every literal P4 prefix.

Unlike the legacy declared-working-box diagnostic, the entrance radii here come
from the deterministic qualified fresh-Live source relation:
  * e_S is exactly zero by shared-origin re-anchoring;
  * v/p are bounded by the uniformly qualified BRMM physical primitives;
  * b_a uses the admitted physical true-bias bound at fresh zero estimate;
  * attitude, gyro-bias and a_w use their qualified startup bounds.

The output is compared with the unchanged full P4 working tube.  The attached
coefficient/source word is still one binary64 point trace, so this file is a
falsifiable feasibility screen only and cannot promote P4.
"""
from __future__ import annotations
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'tools'/'stability'))

import ou3_p4_bounded_bias_cocycle as C
import ou3_p4_motion_gain as G
import ou3_p4_qualified_fresh_entry as QENTRY

GROUPS=("attitude","gyro_bias","velocity","position","integral_displacement","latent_acceleration")
ENTRIES=GROUPS+("accelerometer_bias",)

def prefix_retention(samples,entry_radii,working_radii):
    transition=np.eye(21);response=np.zeros(21);worst=np.zeros(len(GROUPS));where=[None]*len(GROUPS);count=0
    for sample in samples:
        for full,forcing,step in sample['prefixes']:
            phi=full@transition;r=full@response+forcing;count+=1
            for g in range(len(GROUPS)):
                sl=slice(3*g,3*g+3)
                reach=sum(float(np.linalg.norm(phi[sl,3*j:3*j+3],2))*entry_radii[j] for j in range(len(ENTRIES)))
                total=reach+float(np.linalg.norm(r[sl]));ratio=total/working_radii[g]
                if ratio>worst[g]:
                    worst[g]=ratio;where[g]={'index':step['index'],'stage':step['stage'],'kind':step['kind'],'ratio':ratio,'bound':total}
        transition=sample['full']@transition
        response=sample['full']@response+sample['forcing']
    return count,worst,where

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--prefix',required=True,type=Path);ap.add_argument('--attachment',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
    attachment=json.loads(a.attachment.read_text());paths={s:Path(str(a.prefix)+s) for s in ('.root.json','.inputs.csv','.prefixes.jsonl','.events.jsonl')};hashes={s:hashlib.sha256(p.read_bytes()).hexdigest() for s,p in paths.items()}
    if not attachment['read_only_trace_recovers_baseline_bit_for_bit']:raise ValueError('passive trace parity required')
    if hashes['.events.jsonl']!=attachment['event_trace_sha256'] or any(hashes[s]!=attachment['baseline_capture_sha256'][s] for s in ('.root.json','.inputs.csv','.prefixes.jsonl')):raise ValueError('detached source/word capture')
    root=json.loads(paths['.root.json'].read_text());rows=[json.loads(x) for x in paths['.events.jsonl'].read_text().splitlines()];points=[json.loads(x) for x in paths['.prefixes.jsonl'].read_text().splitlines()]
    q=QENTRY.build();qf=QENTRY.validate(q)
    if qf:raise ValueError('qualified fresh entry invalid: '+repr(qf))
    sb=q['source_bounds'];entry=np.array([sb['cayley_norm'],sb['gyro_bias_norm_rad_s'],sb['velocity_norm_mps'],sb['position_norm_m'],sb['centered_S_norm_m_s'],sb['latent_acceleration_norm_mps2'],sb['true_accelerometer_bias_norm_mps2']],dtype=float)
    working=np.array([q['entry_radii']['attitude_cayley_norm'],q['entry_radii']['gyro_bias_norm_rad_s'],q['entry_radii']['velocity_norm_mps'],q['entry_radii']['position_norm_m'],q['entry_radii']['integral_displacement_norm_m_s'],q['entry_radii']['latent_acceleration_norm_mps2']],dtype=float)
    if entry[4]!=0.0:raise ValueError('fresh centered S must be exact zero')
    report={'experiment':'QUALIFIED_FRESH_ENTRY_POINT_PREFIX_RETENTION','capture_sha256':hashes,'entry_radii':dict(zip(ENTRIES,entry.tolist())),'working_radii':dict(zip(GROUPS,working.tolist())),'fresh_centered_S_exact_zero':True,'covariance_membership_used':False,'legacy_300_m_s_entry_factor_used':False,'source_or_coefficient_search_used':False,'trajectory_replay_promoted':False,'binary64_point_diagnostic_only':True,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,'modes':{}}
    for mode in ('H18','A21'):
        if attachment['modes'][mode]['decision']!='CONNECTED_POINT_ATTACHMENT_PASS':raise ValueError('unattached '+mode)
        steps,_,_,defects,counts=G.build_word(root,rows,points,mode)
        if defects.failures:raise ValueError('finite coefficient attachment failed: '+repr(defects.failures[:1]))
        samples=C.sample_lifts(steps);n,worst,where=prefix_retention(samples,entry,working)
        report['modes'][mode]={'sample_count':len(samples),'completed_prefix_count':n,'maximum_working_tube_ratios':dict(zip(GROUPS,map(float,worst))),'limiting_prefix':dict(zip(GROUPS,where)),'all_coordinates_inside_working_tube':bool(np.all(worst<=1.0)),'attitude_inside_70deg_hard_reinit_guard':bool(worst[0]*working[0] < 2*np.tan(np.deg2rad(70.0)/2.0)),'counts':dict(counts)}
        print('QUALIFIED_FRESH_POINT',mode,'inside=',bool(np.all(worst<=1.0)),'max=',float(np.max(worst)),'att=',float(worst[0]))
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
