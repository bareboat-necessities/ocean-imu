"""Exhaustive/fail-closed 1-s magnetic-service stratum certificate driver.

The driver is intentionally generic about the source-domain seed.  A seed is
a finite collection of CLOSED operation strata whose union has been proved by
the caller to cover the theorem domain.  Each leaf is verified with the
literal operation-box exporter and E_hb/non-E_hb Schur certificate.

No carried history is accepted as a coverage proof.
"""
from __future__ import annotations
from dataclasses import dataclass,replace
import json,math
from pathlib import Path
from .interval_riccati_21 import IMat
from .rank_loss_literal_boxes import magnetic_certificate_from_operation_stratum
from .shipping_operation_interval_boxes import split_interval_matrix,widest_entry
from .magnetic_literal_box_export import PredictionBox,CorrectionBox,ResetBox,SyncBox


@dataclass(frozen=True)
class ClosedMagneticStratum:
    stratum_id: str
    root_dim: int
    ops: tuple
    ehb_indices: tuple[int,...]
    nuisance_indices: tuple[int,...]
    coverage_tag: str
    admissible: bool=True


def _imat_width(a:IMat):
    i,j,r=widest_entry(a);return r,(i,j)


def _splittable_fields(op):
    out=[]
    if isinstance(op,PredictionBox):
        out.append(("F",op.F))
    elif isinstance(op,CorrectionBox):
        if op.H is not None: out.append(("H",op.H))
        if op.S_actual is not None: out.append(("S_actual",op.S_actual))
        out.append(("A",op.A))
    elif isinstance(op,ResetBox):
        out.append(("G",op.G))
    return out


def widest_operation_entry(stratum:ClosedMagneticStratum):
    best=None
    for k,op in enumerate(stratum.ops):
        for name,a in _splittable_fields(op):
            r,ij=_imat_width(a)
            # relative scale prevents large deterministic entries dominating.
            m=abs(a.mid[ij[0]][ij[1]])
            score=r/max(1.0,m)
            if best is None or score>best[0]:
                best=(score,k,name,ij,r)
    return best


def _replace_matrix(op,name,a):
    if isinstance(op,PredictionBox) and name=="F": return PredictionBox(a)
    if isinstance(op,ResetBox) and name=="G": return ResetBox(a)
    if isinstance(op,CorrectionBox):
        kw=dict(A=op.A,sensor=op.sensor,applied=op.applied,H=op.H,
                S_actual=op.S_actual,t_lo=op.t_lo,t_hi=op.t_hi)
        kw[name]=a
        return CorrectionBox(**kw)
    raise ValueError("unsupported split target")


def split_stratum(s:ClosedMagneticStratum):
    w=widest_operation_entry(s)
    if w is None or w[4]<=0: return ()
    _,k,name,(i,j),_=w
    a=dict(_splittable_fields(s.ops[k]))[name]
    parts=split_interval_matrix(a,i,j)
    out=[]
    for n,p in enumerate(parts):
        ops=list(s.ops);ops[k]=_replace_matrix(ops[k],name,p)
        out.append(replace(s,stratum_id=f"{s.stratum_id}.{n}",ops=tuple(ops)))
    return tuple(out)


def certify_stratum(s:ClosedMagneticStratum):
    if not s.admissible:
        return {"verified":True,"excluded":True,"gamma_M_lower":math.inf,
                "stratum_id":s.stratum_id,"coverage_tag":s.coverage_tag}
    z=magnetic_certificate_from_operation_stratum(
        s.root_dim,s.ops,s.ehb_indices,s.nuisance_indices)
    z.update({"stratum_id":s.stratum_id,"coverage_tag":s.coverage_tag,
              "excluded":False})
    return z


def exhaustive_certificate(seeds,*,max_depth=14,max_leaves=200000):
    """Branch-and-bound all supplied closed strata.

    Positive completion requires BOTH all leaves verified and a caller-supplied
    nonempty coverage_tag on every seed.  This prevents a finite list of
    diagnostic strata from masquerading as source-domain coverage.
    """
    seeds=tuple(seeds)
    if not seeds:
        return {"verified":False,"gamma_M_lower":0.0,
                "reason":"no source-domain strata supplied","unresolved":[]}
    if any(not s.coverage_tag for s in seeds):
        return {"verified":False,"gamma_M_lower":0.0,
                "reason":"seed missing source-domain coverage proof tag","unresolved":[]}
    stack=[(s,0) for s in seeds];accepted=[];unresolved=[]
    while stack:
        if len(accepted)+len(unresolved)+len(stack)>max_leaves:
            unresolved.append({"stratum_id":stack[-1][0].stratum_id,
                               "reason":"leaf budget exceeded"})
            break
        s,d=stack.pop();z=certify_stratum(s)
        if z.get("verified",False):
            accepted.append(z);continue
        if d>=max_depth:
            unresolved.append(z);continue
        kids=split_stratum(s)
        if not kids:
            unresolved.append(z);continue
        stack.extend((x,d+1) for x in kids)
    finite=[x["gamma_M_lower"] for x in accepted
            if not x.get("excluded") and math.isfinite(x["gamma_M_lower"])]
    lower=min(finite) if finite else 0.0
    ok=not unresolved and bool(finite) and lower>0.0
    return {"qualification":"OU3_MAGNETIC_1S_EXHAUSTIVE_V1",
            "verified":ok,"gamma_M_lower":lower if ok else 0.0,
            "accepted_leaf_count":len(accepted),
            "unresolved_leaf_count":len(unresolved),
            "unresolved":unresolved[:100],
            "coverage_tags":sorted({s.coverage_tag for s in seeds}),
            "max_depth":max_depth,"max_leaves":max_leaves}


def emit(result,path):
    Path(path).write_text(json.dumps(result,indent=2,sort_keys=True,default=str)+"\n")
