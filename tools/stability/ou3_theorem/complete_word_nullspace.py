"""Exact block-nullspace lemma for literal complete-word zero action.

No observability shortcut is encoded here. Callers stack the authoritative
literal homogeneous zero-action rows (process/sync, S, magnetic, accelerometer,
gyro/transport, BA compatibility). Exact rational rank then certifies that the
common kernel is exactly the supplied compatibility line.
"""
from fractions import Fraction as F
from .linked_supply import matrix

def rref(a):
    a=[row[:] for row in matrix(a)]; m=len(a); n=len(a[0]); piv=[]; i=0
    for j in range(n):
        p=next((k for k in range(i,m) if a[k][j]),None)
        if p is None: continue
        a[i],a[p]=a[p],a[i]; q=a[i][j]; a[i]=[x/q for x in a[i]]
        for k in range(m):
            if k!=i and a[k][j]:
                q=a[k][j]; a[k]=[x-q*y for x,y in zip(a[k],a[i])]
        piv.append(j); i+=1
        if i==m: break
    return a,piv

def matvec(a,x):
    return [sum((F(v)*F(w) for v,w in zip(row,x)),F(0)) for row in a]

def stack_blocks(blocks):
    rows=[]; names=[]
    n=None
    for name,b in blocks:
        bb=matrix(b)
        if n is None:n=len(bb[0])
        if len(bb[0])!=n: raise ValueError("common column dimension required")
        rows.extend(bb); names.extend([name]*len(bb))
    if not rows: raise ValueError("nonempty block stack required")
    return rows,names

def compatibility_line_certificate(blocks,line):
    C,names=stack_blocks(blocks); r=[F(x) for x in line]; n=len(C[0])
    if len(r)!=n or all(x==0 for x in r): raise ValueError("nonzero line in stack coordinates required")
    annihilated=all(x==0 for x in matvec(C,r))
    _,piv=rref(C); rank=len(piv); nullity=n-rank
    # rank n-1 plus Cr=0 proves ker C=span(r), exactly.
    return {"rows":len(C),"columns":n,"rank":rank,"nullity":nullity,
            "pivot_columns":piv,"compatibility_line_annihilated":annihilated,
            "kernel_equals_compatibility_line":annihilated and nullity==1,
            "block_row_counts":{name:names.count(name) for name in sorted(set(names))}}

def service_row_stratum_certificate(block_builder, strata):
    """Exact rank audit over a theorem-supplied finite SERVICE-ROW stratum cover.

    This intentionally refuses callback-pattern enumeration. Each stratum must
    provide literal root-coordinate blocks and its declared compatibility line.
    """
    out=[]; all_ok=True
    for st in strata:
        blocks,line=block_builder(st)
        c=compatibility_line_certificate(blocks,line)
        c["stratum"]=st
        out.append(c); all_ok &= c["kernel_equals_compatibility_line"]
    return {"strata":out,"all_strata_kernel_line":all_ok,
            "callback_pattern_enumeration_used":False}

def positive_attitude_then_BA_graph_status(attitude_constraint):
    """Exact rank status of the aggregate attitude block under A*.

    A* strictness means the 3-D attitude Gram is positive definite, hence its
    row stack has rank 3/nullity 0. The physical tilt/BA compatibility line is
    therefore NOT the kernel of an attitude-only constraint; it arises only
    after nuisance/BA elimination in the joint equations.
    """
    C=matrix(attitude_constraint); n=len(C[0]); _,piv=rref(C)
    return {"attitude_rank":len(piv),"attitude_nullity":n-len(piv),
            "Astar_attitude_PD_equivalent":len(piv)==n,
            "attitude_only_nullity_one_claim_valid":False}

def graph_kernel_certificate(attitude_constraint,A):
    """Exact graph-kernel lemma for word-dependent BA compatibility.

    C_theta must have nullity one. A maps attitude coordinates to BA. The
    stacked canonical equations [C_theta 0; A I] then have exactly the graph
    kernel {(a,-Aa): a in ker C_theta}; no fixed word-independent generator is
    assumed.
    """
    Ct=matrix(attitude_constraint); AA=matrix(A)
    p=len(Ct[0]); q=len(AA)
    if len(AA[0])!=p or any(len(row)!=q for row in AA):
        raise ValueError("A must be q by p")
    _,piv=rref(Ct); ntheta=p-len(piv)
    C=[row+[F(0)]*q for row in Ct]
    C += [AA[i]+[F(int(i==j)) for j in range(q)] for i in range(q)]
    _,allp=rref(C); nullity=p+q-len(allp)
    return {"attitude_nullity":ntheta,"graph_stack_nullity":nullity,
            "kernel_is_word_dependent_BA_graph":ntheta==1 and nullity==1,
            "fixed_global_compatibility_vector_required":False}

def certificate():
    blocks=[("four_S_LIN",[[1,0,0,0],[0,1,0,0]]),
            ("magnetic_AG",[[0,0,1,-1]]),
            ("accelerometer_aggregate",[[0,0,2,-2]]),
            ("gyro_process",[[1,1,0,0]])]
    c=compatibility_line_certificate(blocks,[0,0,1,1])
    return {"qualification":"OU3_COMPLETE_WORD_NULLSPACE_V1",**c,
            "literal_event_strata_regressed":False,
            "finite_callback_pattern_cover_exists":False,
            "required_stratum_parameterization":"aggregate MAGNETIC SERVICE rows + literal continuous factor parameters",
            "Astar_persistence_exclusion_separate":True,
            "word_dependent_graph_kernel_required":True,
            "theorem_closed":False}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
