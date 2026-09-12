"""Finite IntegratedOUChain<3> process covariance used by shipping OU-III.

This is the real-arithmetic formula graph behind QdAxis4x1_analytic. It carries
the nested 3x3 marginal and final 4x4 PSD-hygiene branches explicitly. Runtime
exp values, machine epsilon and Eigen LDLT/eigensolver outcomes remain declared
witnesses until deployment finite precision is enclosed.

The accepted-LDLT branch is no longer a free boolean: before a witness may claim
that shipping returned from ``regularize_psd_if_needed`` after LDLT, the SAME
rational matrix must pass an exact symmetric LDL^T inertia check with the same
``-tol`` acceptance threshold. This is deliberately conservative with respect
to Eigen pivoting and does not claim binary32/Eigen correspondence.

The coefficient formula may additionally be selected by an externally derived
``small_branch`` Boolean. The strongest source-bound path supplies that Boolean
from the exact binary32 h/tau branch graph rather than re-deciding it from the
real-rational h/tau quotient. Legacy/local formula tests may omit it.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


def diag(v):
    out=M.zeros(len(v),len(v))
    for i,x in enumerate(v): out[i][i]=x
    return out


def _exact_ldlt_min_d(S):
    n=len(S); A=M.mat(S,n,n)
    if A != M.transpose(A): raise ValueError('exact LDLT requires symmetric input')
    L=M.eye(n); D=[F(0)]*n
    for j in range(n):
        d=A[j][j]-sum((L[j][k]*L[j][k]*D[k] for k in range(j)),F(0))
        D[j]=d
        for i in range(j+1,n):
            r=A[i][j]-sum((L[i][k]*L[j][k]*D[k] for k in range(j)),F(0))
            if d==0:
                if r!=0: return None
                L[i][j]=F(0)
            else:
                L[i][j]=r/d
    return min(D) if D else F(0)


@dataclass(frozen=True)
class PSDWitness:
    ldlt_accepts: bool
    eigensolver_success: bool | None = None
    eigenvectors: tuple | None = None
    eigenvalues: tuple | None = None
    def __post_init__(self):
        if not isinstance(self.ldlt_accepts,bool): raise TypeError('literal LDLT acceptance branch required')
        if self.ldlt_accepts:
            if any(x is not None for x in (self.eigensolver_success,self.eigenvectors,self.eigenvalues)):
                raise ValueError('accepted LDLT branch consumes no eigensolver witness')
        else:
            if not isinstance(self.eigensolver_success,bool): raise TypeError('rejected LDLT branch requires eigensolver outcome')
            if not self.eigensolver_success and (self.eigenvectors is not None or self.eigenvalues is not None):
                raise ValueError('failed eigensolver branch consumes no eigensystem')
            if self.eigensolver_success and (self.eigenvectors is None or self.eigenvalues is None):
                raise ValueError('successful eigensolver requires eigensystem')


def regularize_psd(S,witness:PSDWitness,*,machine_epsilon):
    n=len(S); S=M.mat(S,n,n); eps=P.rational(machine_epsilon)
    if n not in (3,4) or eps <= 0: raise ValueError('N=3/4 and positive machine epsilon required')
    S=M.scaled(M.plus(S,M.transpose(S)),F(1,2))
    scale=max(F(1),max(abs(x) for row in S for x in row)); tol=64*eps*scale
    if witness.ldlt_accepts:
        min_d=_exact_ldlt_min_d(S)
        if min_d is None or min_d < -tol:
            raise ValueError('LDLT-accept witness detached from SAME PSD-hygiene matrix/tolerance')
        return S
    if not witness.eigensolver_success:
        out=[r[:] for r in S]
        for i in range(n): out[i][i]+=tol
        return M.scaled(M.plus(out,M.transpose(out)),F(1,2))
    U=M.mat(witness.eigenvectors,n,n); lam=P.vec(witness.eigenvalues,n)
    if M.mm(M.transpose(U),U) != M.eye(n): raise ValueError('orthonormal eigensystem required')
    if M.mm(S,U) != M.mm(U,diag(lam)): raise ValueError('eigensystem detached from SAME PSD-hygiene input')
    clipped=[max(F(0),x) for x in lam]
    out=M.mm(M.mm(U,diag(clipped)),M.transpose(U))
    return M.scaled(M.plus(out,M.transpose(out)),F(1,2))


def _branch(x, small_branch):
    if small_branch is None: return abs(x) < F(1,100)
    if not isinstance(small_branch,bool): raise TypeError('literal Qaxis coefficient branch required')
    return small_branch


def marginal_raw(tau,h,sigma2,alpha,*,small_branch=None):
    """Literal IntegratedOUChain<T,2> formula under the selected deployed branch."""
    tau,h,sigma2,alpha=map(P.rational,(tau,h,sigma2,alpha))
    if tau < F(1,10**7) or h <= 0 or sigma2 < 0 or not 0 < alpha <= 1: raise ValueError('post-clamp tau, positive h, sigma2>=0, valid exp root required')
    inv=1/tau; x=h*inv; i2=inv*inv; i3=i2*inv; i4=i3*inv; i5=i4*inv; i6=i5*inv; i7=i6*inv; i8=i7*inv; i9=i8*inv
    h2=h*h; h3=h2*h; h4=h3*h; h5=h4*h; h6=h5*h; h7=h6*h; h8=h7*h; h9=h8*h
    Q=M.zeros(3,3)
    if _branch(x,small_branch):
        Q[0][0]=sigma2*(F(2,3)*h3*inv-F(1,2)*h4*i2+F(7,30)*h5*i3-F(1,12)*h6*i4+F(31,1260)*h7*i5-F(1,160)*h8*i6+F(127,90720)*h9*i7)
        Q[0][1]=sigma2*(F(1,4)*h4*inv-F(1,6)*h5*i2+F(5,72)*h6*i3-F(1,45)*h7*i4+F(17,2880)*h8*i5-F(41,30240)*h9*i6)
        Q[0][2]=sigma2*(h2*inv-h3*i2+F(7,12)*h4*i3-F(1,4)*h5*i4+F(31,360)*h6*i5-F(1,40)*h7*i6+F(127,20160)*h8*i7-F(17,12096)*h9*i8)
        Q[1][1]=sigma2*(F(1,10)*h5*inv-F(1,18)*h6*i2+F(5,252)*h7*i3-F(1,180)*h8*i4+F(17,12960)*h9*i5)
        Q[1][2]=sigma2*(F(1,3)*h3*inv-F(1,3)*h4*i2+F(11,60)*h5*i3-F(13,180)*h6*i4+F(19,840)*h7*i5-F(1,168)*h8*i6+F(247,181440)*h9*i7)
        Q[2][2]=sigma2*(2*h*inv-2*h2*i2+F(4,3)*h3*i3-F(2,3)*h4*i4+F(4,15)*h5*i5-F(4,45)*h6*i6+F(8,315)*h7*i7-F(2,315)*h8*i8+F(4,2835)*h9*i9)
    else:
        a=alpha; a2=a*a; qc=2*sigma2*inv; t2=tau*tau; t3=t2*tau; t4=t3*tau; t5=t4*tau; x2=x*x; x3=x2*x
        K00=t3*(-a2+4*a+2*x-3)/2
        K01=t4*(a2+2*a*(x-1)+x2-2*x+1)/2
        K02=t2*(a2-2*a+1)/2
        K11=t5*(-a2/2-2*a*x+x3/3-x2+x+F(1,2))
        K12=t3*(-a2-2*a*x+1)/2
        K22=tau*(1-a2)/2
        Q[0][0]=qc*K00; Q[0][1]=qc*K01; Q[0][2]=qc*K02; Q[1][1]=qc*K11; Q[1][2]=qc*K12; Q[2][2]=qc*K22
    Q[1][0]=Q[0][1]; Q[2][0]=Q[0][2]; Q[2][1]=Q[1][2]
    return Q


def qaxis4(tau,h,sigma2,alpha,*,marginal_psd:PSDWitness,final_psd:PSDWitness,machine_epsilon,small_branch=None):
    """Literal IntegratedOUChain<T,3> graph under one shared formula branch."""
    tau,h,sigma2,alpha=map(P.rational,(tau,h,sigma2,alpha)); inv=1/tau; x=h*inv
    small=_branch(x,small_branch)
    marginal=regularize_psd(marginal_raw(tau,h,sigma2,alpha,small_branch=small),marginal_psd,machine_epsilon=machine_epsilon)
    Q=M.zeros(4,4); idx=(0,1,3)
    for i in range(3):
        for j in range(3): Q[idx[i]][idx[j]]=marginal[i][j]
    if small:
        i2=inv*inv; i3=i2*inv; i4=i3*inv; i5=i4*inv; i6=i5*inv
        h2=h*h; h3=h2*h; h4=h3*h; h5=h4*h; h6=h5*h; h7=h6*h; h8=h7*h; h9=h8*h
        qvS=sigma2*(F(1,15)*h5*inv-F(1,24)*h6*i2+F(41,2520)*h7*i3-F(7,1440)*h8*i4+F(109,90720)*h9*i5)
        qpS=sigma2*(F(1,36)*h6*inv-F(1,72)*h7*i2+F(13,2880)*h8*i3-F(1,864)*h9*i4)
        qSS=sigma2*(F(1,126)*h7*inv-F(1,288)*h8*i2+F(13,12960)*h9*i3)
        qSa=sigma2*(F(1,12)*h4*inv-F(1,12)*h5*i2+F(2,45)*h6*i3-F(1,60)*h7*i4+F(11,2240)*h8*i5-F(73,60480)*h9*i6)
    else:
        a=alpha; a2=a*a; qc=2*sigma2*inv; t4=tau**4; t5=tau**5; t6=tau**6; t7=tau**7; x2=x*x; x3=x2*x; x4=x3*x; x5=x4*x
        K02=t5*(-3*a2+3*a*(x2+4)+x3-3*x2+6*x-9)/6
        K12=t6*(a2/2+a*(-x2+2*x-2)/2+x4/8-x3/2+x2-x+F(1,2))
        K22=t7*(-a2/2+a*x2+2*a+x5/20-x4/4+F(2,3)*x3-x2+x-F(3,2))
        K23=t4*(a2-a*(x2+2)+1)/2
        qvS=qc*K02; qpS=qc*K12; qSS=qc*K22; qSa=qc*K23
    Q[0][2]=Q[2][0]=qvS; Q[1][2]=Q[2][1]=qpS; Q[2][2]=qSS; Q[2][3]=Q[3][2]=qSa
    return regularize_psd(Q,final_psd,machine_epsilon=machine_epsilon)


def readiness():
    return {
      'Qaxis_small_branch_formula_materialized':True,
      'Qaxis_general_branch_formula_materialized':True,
      'nested_marginal_and_final_formula_branch_can_be_forced_from_one_deployed_decision':True,
      'nested_marginal_psd_hygiene_materialized':True,
      'final_Qaxis_psd_hygiene_materialized':True,
      'free_Qaxis_matrix_removed_by_this_lemma':True,
      'LDLT_accept_branch_has_same_matrix_exact_inertia_guard':True,
      'arbitrary_LDLT_accept_boolean_can_bypass_matrix_relation':False,
      'alpha_exp_runtime_source_attached':False,
      'machine_epsilon_deployment_attached':False,
      'Eigen_LDLT_eigensolver_outcomes_attached':False,
      'nonfinite_replacement_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
