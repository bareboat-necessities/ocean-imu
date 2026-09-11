"""Exact rank-three signed Joseph/reset algebra for the actual shipping numerator.

Unlike the ordinary Kalman information identity this does not assume N=P H^T,
which is false for the H18 held-bias measurement path. Let J=P^-1, K=N S^-1,
and P+=P-N S^-1 N^T. Define the 3D core

    L=N^T J,   C=S-N^T J N.

If P+ is SPD (equivalently C is SPD), then

    J+=J+L^T C^-1 L,
    J+ K=L^T C^-1,
    K^T J+ K=C^-1-S^-1.

For physical residual q and state error e, xi=q-L e gives the exact signed
energy identity

    V_J-V = -q^T S^-1 q + xi^T C^-1 xi.

If the finite reset is e_r=G t+rho, b=G^-1 rho, t=e-Kq, then

    V_r-V = -q^T S^-1 q + xi^T C^-1 xi
            +2 e^T J b -2 xi^T C^-1 Lb
            +b^T J b +(Lb)^T C^-1(Lb).

Only S and C are 3x3. The same identities are inverse-free with S u=q,
C v=xi and C w=Lb. This module is algebra only; source-uniform xi/reset sectors
and complete-word dissipation remain separate obligations.
"""
from __future__ import annotations
from fractions import Fraction as F

QUALIFICATION="OU3_ALT_ACTUAL_NUMERATOR_RANK3_SIGNED_JOSEPH_V1"

def shape(A):
    r=len(A);c=len(A[0]) if r else 0
    if any(len(x)!=c for x in A): raise ValueError("ragged matrix")
    return r,c

def mt(A): return [list(x) for x in zip(*A)]
def mm(A,B):
    if shape(A)[1]!=shape(B)[0]: raise ValueError("dimension mismatch")
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))),F(0)) for j in range(len(B[0]))] for i in range(len(A))]
def mv(A,x): return [sum((A[i][j]*x[j] for j in range(len(x))),F(0)) for i in range(len(A))]
def add(A,B): return [[A[i][j]+B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def sub(A,B): return [[A[i][j]-B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def vadd(a,b): return [x+y for x,y in zip(a,b)]
def vsub(a,b): return [x-y for x,y in zip(a,b)]
def dot(a,b): return sum((x*y for x,y in zip(a,b)),F(0))
def quad(A,x): return dot(x,mv(A,x))
def eye(n): return [[F(1 if i==j else 0) for j in range(n)] for i in range(n)]
def zero(A): return all(x==0 for row in A for x in row)

def inverse(A):
    n,m=shape(A)
    if n!=m or n==0: raise ValueError("square matrix required")
    aug=[list(A[i])+eye(n)[i] for i in range(n)]
    for j in range(n):
        p=next((i for i in range(j,n) if aug[i][j]!=0),None)
        if p is None: raise ValueError("singular matrix")
        aug[j],aug[p]=aug[p],aug[j]
        d=aug[j][j]; aug[j]=[x/d for x in aug[j]]
        for i in range(n):
            if i==j: continue
            a=aug[i][j]
            if a: aug[i]=[x-a*y for x,y in zip(aug[i],aug[j])]
    return [r[n:] for r in aug]

def objects(P,S,N):
    n=shape(P)[0]
    if shape(P)!=(n,n) or shape(S)!=(3,3) or shape(N)!=(n,3): raise ValueError("expected P nxn, S 3x3, N nx3")
    J=inverse(P);Si=inverse(S);K=mm(N,Si);Pp=sub(P,mm(mm(N,Si),mt(N)));Jp=inverse(Pp)
    L=mm(mt(N),J);C=sub(S,mm(L,N));Ci=inverse(C)
    return dict(J=J,Si=Si,K=K,Pp=Pp,Jp=Jp,L=L,C=C,Ci=Ci)

def signed_terms(P,S,N,e,q,b=None):
    o=objects(P,S,N);J=o['J'];Jp=o['Jp'];K=o['K'];L=o['L'];Ci=o['Ci'];Si=o['Si']
    xi=vsub(q,mv(L,e));t=vsub(e,mv(K,q));before=quad(J,e);after=quad(Jp,t)
    joseph=-quad(Si,q)+quad(Ci,xi)
    out=dict(xi=xi,t=t,joseph_identity_residual=(after-before)-joseph,joseph_signed=joseph)
    if b is not None:
        Lb=mv(L,b);reduced=joseph+2*dot(e,mv(J,b))-2*dot(xi,mv(Ci,Lb))+quad(J,b)+quad(Ci,Lb)
        direct=quad(Jp,vadd(t,b))-before
        out.update(reset_identity_residual=direct-reduced,reset_reduced_signed=reduced,Lb=Lb)
    return out

def descriptor_terms(P,S,N,e,q,b=None):
    o=objects(P,S,N);L=o['L'];C=o['C'];J=o['J'];xi=vsub(q,mv(L,e))
    u=mv(o['Si'],q);v=mv(o['Ci'],xi);out=dict(u=u,v=v,xi=xi,S_u_minus_q=vsub(mv(S,u),q),C_v_minus_xi=vsub(mv(C,v),xi),joseph_signed=-dot(q,u)+dot(xi,v))
    if b is not None:
        Lb=mv(L,b);w=mv(o['Ci'],Lb);out.update(w=w,C_w_minus_Lb=vsub(mv(C,w),Lb),reset_reduced_signed=out['joseph_signed']+2*dot(e,mv(J,b))-2*dot(xi,w)+quad(J,b)+dot(Lb,w))
    return out

def masked_regression():
    P=[[F(3),F(1,2),F(1,5)],[F(1,2),F(2),F(1,4)],[F(1,5),F(1,4),F(5,2)]]
    H=[[F(1),F(2),F(-1)],[F(-1),F(1),F(1,2)],[F(1,3),F(-2,3),F(1)]]
    R=[[F(4),F(1,5),F(0)],[F(1,5),F(5),F(1,7)],[F(0),F(1,7),F(6)]]
    S=add(mm(mm(H,P),mt(H)),R);Hg=[r[:] for r in H]
    for i in range(3): Hg[i][2]=F(0)
    N=mm(P,mt(Hg));N[2]=[F(0),F(0),F(0)]
    e=[F(2,7),F(-1,5),F(1,9)];eta=[F(1,13),F(-2,17),F(1,19)];q=vadd(mv(H,e),eta);b=[F(1,101),F(-1,97),F(1,89)]
    o=objects(P,S,N);info=add(o['J'],mm(mm(mt(o['L']),o['Ci']),o['L']))
    return dict(Jplus=zero(sub(o['Jp'],info)),JplusK=zero(sub(mm(o['Jp'],o['K']),mm(mt(o['L']),o['Ci']))),KtJplusK=zero(sub(mm(mm(mt(o['K']),o['Jp']),o['K']),sub(o['Ci'],o['Si']))),Ldiff=o['L']!=H,Cdiff=o['C']!=R,signed=signed_terms(P,S,N,e,q,b),descriptor=descriptor_terms(P,S,N,e,q,b))

def unmasked_regression():
    P=[[F(3),F(1,2)],[F(1,2),F(2)]];H=[[F(1),F(2)],[F(-1),F(1)],[F(1,2),F(-1,3)]];R=[[F(3),F(1,5),F(0)],[F(1,5),F(4),F(1,7)],[F(0),F(1,7),F(5)]]
    S=add(mm(mm(H,P),mt(H)),R);o=objects(P,S,mm(P,mt(H)));return o['L']==H,o['C']==R

def build():
    m=masked_regression();ul,uc=unmasked_regression();s=m['signed'];d=m['descriptor']
    return {'qualification':QUALIFICATION,'actual_shipping_numerator_is_primitive':True,'unqualified_N_equals_P_Htranspose_required':False,'rank_three_core_dimension':3,'posterior_information_identity_exact':m['Jplus'],'posterior_gain_cross_identity_exact':m['JplusK'],'posterior_gain_energy_identity_exact':m['KtJplusK'],'masked_regression_L_differs_from_physical_H':m['Ldiff'],'masked_regression_C_differs_from_physical_R':m['Cdiff'],'unmasked_limit_recovers_L_equals_H':ul,'unmasked_limit_recovers_C_equals_R':uc,'exact_masked_Joseph_signed_identity_closed':s['joseph_identity_residual']==0,'exact_masked_Joseph_reset_signed_identity_closed':s['reset_identity_residual']==0,'inverse_free_S_descriptor_closed':all(x==0 for x in d['S_u_minus_q']),'inverse_free_C_xi_descriptor_closed':all(x==0 for x in d['C_v_minus_xi']),'inverse_free_C_reset_descriptor_closed':all(x==0 for x in d['C_w_minus_Lb']),'independent_K_box_used':False,'physical_residual_Jacobian_substituted_for_masked_L':False,'source_uniform_C_positive_definite_closed_here':False,'source_uniform_xi_sector_closed_here':False,'source_uniform_reset_absorption_closed_here':False,'source_uniform_complete_word_dissipation_closed_here':False,'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False}
def validate(d):
    f=[]
    for k in ('actual_shipping_numerator_is_primitive','posterior_information_identity_exact','posterior_gain_cross_identity_exact','posterior_gain_energy_identity_exact','masked_regression_L_differs_from_physical_H','masked_regression_C_differs_from_physical_R','unmasked_limit_recovers_L_equals_H','unmasked_limit_recovers_C_equals_R','exact_masked_Joseph_signed_identity_closed','exact_masked_Joseph_reset_signed_identity_closed','inverse_free_S_descriptor_closed','inverse_free_C_xi_descriptor_closed','inverse_free_C_reset_descriptor_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('unqualified_N_equals_P_Htranspose_required','independent_K_box_used','physical_residual_Jacobian_substituted_for_masked_L','source_uniform_C_positive_definite_closed_here','source_uniform_xi_sector_closed_here','source_uniform_reset_absorption_closed_here','source_uniform_complete_word_dissipation_closed_here','ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    if d.get('rank_three_core_dimension')!=3:f.append('rank-three dimension changed')
    return f
