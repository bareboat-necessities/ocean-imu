#!/usr/bin/env python3
"""Non-promoting spacing study for the limiting eta6/a_w H18 information block.

The current four-S lemma selects actual shipping S=0 firings in windows
[0,g],[2g,3g],[4g,5g],[6g,7g].  The canonical BRMM word is 3 s long, so this
study asks whether the same scheduler guarantee can select four firings farther
apart while remaining inside the SAME word and source.

For integer spacing m>=2 use windows
  [0,g], [m g,(m+1)g], [2m g,(2m+1)g], [3m g,(3m+1)g].
All other due S firings remain in the literal word.  Every bound below retains
worst-case OU decay, process nuisance, actual applied SpectralMSE R_S and the
same frozen 1e-18 gate.  This file is diagnostic only and cannot promote P3/P4.
"""
from __future__ import annotations
import json, math
from fractions import Fraction

import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_brmm_complete_source as COMPLETE
import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_brmm_windowed_vector_pe as PE
import ou3_brmm_four_s_translation_information_tight as BASE

GATE=1.0e-18
WORD_S=3.0

def down(x): return math.nextafter(float(x),-math.inf)
def up(x): return math.nextafter(float(x),math.inf)

def dd_rows(m:int):
    # Guaranteed scaled pair separations for windows [jm,jm+1].
    d01=Fraction(m-1); d02=Fraction(2*m-1); d03=Fraction(3*m-1)
    d12=Fraction(m-1); d13=Fraction(2*m-1); d23=Fraction(m-1)
    q0=Fraction(1)
    q1=Fraction(2,1)/d01
    q2=sum((Fraction(1,1)/(d01*d02),Fraction(1,1)/(d01*d12),Fraction(1,1)/(d02*d12)),Fraction(0))
    q3=sum((Fraction(1,1)/(d01*d02*d03),Fraction(1,1)/(d01*d12*d13),Fraction(1,1)/(d02*d12*d23),Fraction(1,1)/(d03*d13*d23)),Fraction(0))
    return list(map(float,(q0,q1,q2,q3)))

def candidate(m:int,dyn:dict,comp:dict,sched:dict,pe:dict):
    inv=dyn['dynamic_invariant']; g=up(float(sched['certified_uniform_max_gap_s']))
    T=up((3*m+1)*g)
    if T>WORD_S: return None
    tau_lo=float(inv['tau_applied_s'][0]); sigma_hi=float(inv['sigma_aw_filter_mps2'][1]); rs_hi=float(inv['R_S_applied'][1])
    q0,q1,q2,q3=dd_rows(m)
    # c'''(u)>=exp(-T/tau_lo), hence third divided difference >= /6.
    third=down(math.exp(-up(T/tau_lo))/6.0)
    if third<=0: return None
    c0=up(1/6)
    c01=up(((m+1)**2)/2.0)      # max c'(u), u<=m+1
    c012=up((2*m+1)/2.0)        # half max c''(u), u<=2m+1
    halfsum=up((m+2)/2.0)       # (u0+u1)/2 upper
    rowA=up(q3/third)
    rowV=up(2*up(q2+up(c012*rowA)))
    rowP=up(q1+up(halfsum*rowV)+up(c01*rowA))
    rowS=up(q0+rowP+up(.5*rowV)+up(c0*rowA))
    rows={'S':rowS,'g*p':rowP,'g^2*v':rowV,'g^3*a_w':rowA}

    rs=comp['R_S_regularizer']; factors=list(map(float,rs['axis_std_factors']))
    meas=max(up(up(rs_hi*f)**2) for f in factors)
    qc=up(2*sigma_hi*sigma_hi/tau_lo)
    proc=up(qc*up(T**7)/252.0)
    lam=up(meas+up(4*proc))
    directional={k:down(1.0/up(4*up(lam*up(v*v)))) for k,v in rows.items()}

    alpha=float(pe['eta6_information']['alpha_6_information_lower'])
    ra=float(pe['measurement_runtime']['accelerometer_variance_upper'])
    # Current proof's conservative cross, for apples-to-apples comparison.
    cross_current=up(2.0/(ra*(g**6)))
    # Same-history OU tightening: second required PE occurrence starts >=2W.
    W=float(pe['declared_normal_live_PE']['recurrence_window_s'])
    tau_hi=float(inv['tau_applied_s'][1])
    second_sq=up(math.exp(-down(4.0*W/tau_hi)))
    cross_decay=up((1.0+second_sq)/(ra*(g**6)))

    d=directional['g^3*a_w']
    def coupled(c2): return down((alpha*d)/up(alpha+d+c2))
    old=coupled(cross_current); improved=coupled(cross_decay)
    return {'m':m,'windows_scaled':[[0,1],[m,m+1],[2*m,2*m+1],[3*m,3*m+1]],'last_window_end_s':T,
            'third_divided_difference_lower':third,'inverse_row_l1_upper':rows,
            'measurement_covariance_lambda_max_upper':meas,'process_four_record_covariance_trace_upper':up(4*proc),'total_covariance_lambda_max_upper':lam,
            'directional_information_lower':directional,'aw_cross_current':cross_current,'aw_cross_with_second_PE_OU_decay':cross_decay,
            'coupled_eta6_aw_current_cross':old,'coupled_eta6_aw_with_decay':improved,
            'gate_ratio_with_decay':improved/GATE,'gate_pass_with_decay':improved>=GATE}

def build():
    dyn=DYNAMIC.build(); comp=COMPLETE.build(); sched=SCHED.build(); pe=PE.build(); base=BASE.build()
    bad={'dynamic':DYNAMIC.validate(dyn),'complete':COMPLETE.validate(comp),'scheduler':SCHED.validate(sched),'PE':PE.validate(pe),'base':BASE.validate(base)}
    bad={k:v for k,v in bad.items() if v}
    if bad: raise RuntimeError('spacing sweep prerequisites failed: '+repr(bad))
    rows=[]
    for m in range(2,20):
        c=candidate(m,dyn,comp,sched,pe)
        if c is not None: rows.append(c)
    if not rows: raise RuntimeError('no spacing candidates fit canonical word')
    best=max(rows,key=lambda x:x['coupled_eta6_aw_with_decay'])
    return {'qualification':'OU3_P4_ETA6_AW_FOUR_S_SPACING_SWEEP_V1','canonical_word_s':WORD_S,'gate':GATE,
            'filter_changed':False,'gate_changed':False,'source_changed':False,'diagnostic_only':True,
            'current_m':2,'candidates':rows,'best':best,
            'improvement_over_current_7p092e19':best['coupled_eta6_aw_with_decay']/7.092471820569812e-19,
            'P3_promoted':False,'P4_PASS':False,'P5_MAY_START':False}

def main():
    d=build(); print(json.dumps(d,indent=2,sort_keys=True)); return 0
if __name__=='__main__': raise SystemExit(main())
