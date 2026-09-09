#!/usr/bin/env python3
"""Relevant-subspace signed Joseph diagnostic for exact finite-angle vector chords.

For isotropic accelerometer/magnetometer R=r I, write A=R^-1/2 H P^1/2.
The exact chord obeys ||q||^2>=k||p||^2 and ||q-p||^2=||p||^2-||q||^2
<= (1-k)||p||^2.  Therefore

 q^T S^-1 q - (q-p)^T R^-1(q-p)
 >= [ k/(1+a2) - (1-k) ] p^T R^-1 p,

where a2>=||A||^2.  Unlike the rejected whole-state Ptrace/Smax bound, a2 is
formed only from covariance coordinates actually touched by that measurement H.
This is still a diagnostic scalarization; a negative factor means the production
proof must retain the dense same-history S/K direction in the augmented master.
A positive factor is a valid source-uniform lower that can be fed back into the
universal information composition.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_brmm_riccati_tube_smallx_scaled as TUBE
import ou3_p4_complete_brmm_exact_chord_joint_coordinate as CHORD


def up(x):return math.nextafter(float(x),math.inf)
def down(x):return math.nextafter(float(x),-math.inf)

def build():
    t=TUBE.build_base();c=CHORD.build();vf=CHORD.validate(c)
    if vf:raise RuntimeError('chord invalid: '+repr(vf))
    k=float(c['information_retention_factor_lower_full_entry'])
    # Config/runtime vector bounds.
    fmax=13.80665;mmax=200.0;ra=down(0.2**2);rm=down(0.3**2)
    modes={}
    for mode,key in (('H18','H'),('A21','A')):
        p=list(map(float,t['modes'][key]['Pbar_diagonal_variance_upper']))
        patt=max(p[0:3]);paw=sum(p[15:18]);pba=sum(p[18:21]) if mode=='A21' else 0.0
        # ||[v]x||_F^2=2||v||^2; weighting each attitude column by max Pii.
        acc_hph=up(2.0*fmax*fmax*patt + paw + pba)
        mag_hph=up(2.0*mmax*mmax*patt)
        events={}
        for ev,hph,r in (('accelerometer',acc_hph,ra),('magnetometer',mag_hph,rm)):
            a2=up(hph/r)
            factor=down(k/up(1.0+a2) - (1.0-k))
            events[ev]={'HPHt_lambda_max_upper_relevant_subspace':hph,'normalized_A_norm_squared_upper':a2,
              'signed_chord_information_factor_lower':factor,'nonnegative':factor>=0}
        modes[mode]={'attitude_Pii_max':patt,'aw_trace_upper':paw,'ba_trace_upper':pba,'events':events}
    return {'qualification':'OU3_P4_EXACT_CHORD_SIGNED_JOSEPH_RELEVANT_SUBSPACE_V1','canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'full_entry_chord_k_lower':k,'whole_state_Ptrace_not_used':True,'relevant_H_subspace_only':True,
      'dense_same_history_master_still_preferred_if_negative':True,'modes':modes,
      'all_event_scalar_factors_nonnegative':all(e['nonnegative'] for m in modes.values() for e in m['events'].values()),
      'P4_promoted_here':False}
def validate(d):
    f=[]
    for k in ('whole_state_Ptrace_not_used','relevant_H_subspace_only','dense_same_history_master_still_preferred_if_negative'):
        if d.get(k) is not True:f.append(k+' not true')
    if d.get('P4_promoted_here') is not False:f.append('P4 promoted')
    for mode,m in d['modes'].items():
        for ev,e in m['events'].items():
            for k in ('HPHt_lambda_max_upper_relevant_subspace','normalized_A_norm_squared_upper','signed_chord_information_factor_lower'):
                if not math.isfinite(float(e[k])):f.append(mode+' '+ev+' '+k+' nonfinite')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())