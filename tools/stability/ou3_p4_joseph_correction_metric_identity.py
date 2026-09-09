#!/usr/bin/env python3
"""Exact same-event correction-metric identity for OU-III P4 reset transport.

For one accepted shipping measurement with prior covariance P, residual
Jacobian H, noise R, innovation S=H P H^T+R, Kalman gain K=P H^T S^-1 and
Joseph posterior P+, the information-form identities give

    P+^-1 = P^-1 + H^T R^-1 H,
    K      = P+ H^T R^-1,
    S^-1  = R^-1 - R^-1 H P+ H^T R^-1.

Consequently

    K^T P+^-1 K = R^-1 - S^-1.                         (1)

For the actual full-state correction k=K y,

    k^T P+^-1 k = y^T (R^-1-S^-1) y.                   (2)

This is the complementary innovation information paid by the SAME Joseph event.
It is exactly what finite reset transport needs: the correction coordinate can
be retained in the posterior metric without a source-uniform marginal P ceiling
or an independent norm box for K.

NumPy is imported lazily inside ``build`` so the retained stability package can
be imported on the minimal source-foundation runner.  The numerical smoke still
requires NumPy when executed; the theorem identity itself is algebraic.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

SCHEMA=1
QUALIFICATION='OU3_P4_JOSEPH_CORRECTION_POSTERIOR_METRIC_IDENTITY_V1'

def _check(np,P,H,R):
    P=np.asarray(P,dtype=float);H=np.asarray(H,dtype=float);R=np.asarray(R,dtype=float)
    S=H@P@H.T+R
    K=P@H.T@np.linalg.inv(S)
    Pp=P-K@S@K.T
    lhs=K.T@np.linalg.inv(Pp)@K
    rhs=np.linalg.inv(R)-np.linalg.inv(S)
    err=float(np.max(np.abs(lhs-rhs)))
    scale=max(1.0,float(np.max(np.abs(rhs))))
    return err,err<=5e-12*scale

def build():
    import numpy as np
    cases=[]
    for P,H,R in (
      ([[2.,0.,0.],[0.,3.,0.],[0.,0.,4.]],[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],[[.2,0.,0.],[0.,.3,0.],[0.,0.,.4]]),
      ([[2.,.2,0,.1],[.2,1.4,.1,0],[0,.1,1.8,.3],[.1,0,.3,2.2]],
       [[1.,.2,0,.1],[0,.7,.3,-.2],[.1,0,.8,.4]],
       [[.4,.03,0],[.03,.5,.02],[0,.02,.6]]),
      ([[1.5,.1,.05],[.1,2.0,-.08],[.05,-.08,1.2]],
       [[.3,-.4,.8],[.7,.2,-.1],[-.2,.5,.6]],
       [[.25,.01,.02],[.01,.32,0],[.02,0,.29]]),
    ):
        err,ok=_check(np,P,H,R);cases.append({'max_abs_identity_error':err,'pass':ok})
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'identity':'K^T(Pplus)^-1K=R^-1-S^-1',
      'correction_energy_identity':'(Ky)^T(Pplus)^-1(Ky)=y^T(R^-1-S^-1)y',
      'same_event_P_H_R_S_K_required':True,
      'actual_K_eliminated_only_by_exact_same_event_identity':True,
      'independent_K_box_used':False,'marginal_Ptheta_ceiling_used':False,
      'complementary_innovation_information_retained':True,
      'usable_as_reset_correction_supply_coordinate':True,
      'numpy_import_is_lazy_execution_only':True,
      'smoke_cases':cases,'smoke_identity_closed':all(x['pass'] for x in cases),
      'filter_changed':False,'declared_domain_changed':False,'P4_promoted_here':False,
    }
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_event_P_H_R_S_K_required','actual_K_eliminated_only_by_exact_same_event_identity','complementary_innovation_information_retained','usable_as_reset_correction_supply_coordinate','numpy_import_is_lazy_execution_only','smoke_identity_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_K_box_used','marginal_Ptheta_ceiling_used','filter_changed','declared_domain_changed','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'closed':d['smoke_identity_closed'],'errors':[x['max_abs_identity_error'] for x in d['smoke_cases']],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
