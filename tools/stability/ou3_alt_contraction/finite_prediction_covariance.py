"""Exact block structure of the shipping OU-III covariance prediction.

This closes only the algebra between one prediction predecessor and the point
immediately before pending a_w covariance inflation, symmetrization hygiene and
periodic S=0 service. It does not qualify the primitives that generate the
attitude/noise blocks and it is not a storage certificate.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P

N, NA, NL, NB, OFF_L, OFF_B = 21, 6, 12, 3, 6, 18


def linear_transition(axis_coefficients):
    """Shipping [v,p,S,a_w] 12x12 transition from the SAME mean coefficients."""
    if len(axis_coefficients) != 3:
        raise ValueError('three axis coefficient tuples required')
    out = M.zeros(NL, NL); idx = (0, 3, 6, 9)
    for axis, raw in enumerate(axis_coefficients):
        va, pa, Sa, alpha, h = P.vec(raw, 5)
        phi = [[1,0,0,va], [h,1,0,pa], [h*h/F(2),h,1,Sa], [0,0,0,alpha]]
        for i in range(4):
            for j in range(4): out[idx[i]+axis][idx[j]+axis] = phi[i][j]
    return out


def correlated_linear_process(qaxis_unit, sigma_aw):
    """Literal correlated-Q assembly: each group block is Sigma_aw*qaxis[g,h]."""
    q, sig = M.mat(qaxis_unit,4,4), M.mat(sigma_aw,3,3)
    if q != M.transpose(q) or sig != M.transpose(sig):
        raise ValueError('symmetric qaxis and Sigma_aw required')
    out=M.zeros(NL,NL); goff=(0,3,6,9)
    for g in range(4):
        for h in range(4):
            for i in range(3):
                for j in range(3): out[goff[g]+i][goff[h]+j]=sig[i][j]*q[g][h]
    return out


def independent_linear_process(qaxis):
    """Literal independent-axis Q assembly used by the shipping fallback branch."""
    if len(qaxis) != 3: raise ValueError('three per-axis Qd matrices required')
    out=M.zeros(NL,NL); idx=(0,3,6,9)
    for axis, qa in enumerate(qaxis):
        qa=M.mat(qa,4,4)
        if qa != M.transpose(qa): raise ValueError('symmetric per-axis Qd required')
        for i in range(4):
            for j in range(4): out[idx[i]+axis][idx[j]+axis]=qa[i][j]
    return out


@dataclass(frozen=True)
class Blocks:
    F_AA: tuple; Q_AA: tuple; F_LL: tuple; Q_LL: tuple
    phi_b: F; Q_BB: tuple; active_bias: bool
    def __post_init__(self):
        fa=M.mat(self.F_AA,NA,NA); qa=M.mat(self.Q_AA,NA,NA)
        fl=M.mat(self.F_LL,NL,NL); ql=M.mat(self.Q_LL,NL,NL)
        qb=M.mat(self.Q_BB,NB,NB); ph=P.rational(self.phi_b)
        if qa != M.transpose(qa) or ql != M.transpose(ql) or qb != M.transpose(qb):
            raise ValueError('symmetric process covariance blocks required')
        if not isinstance(self.active_bias,bool): raise TypeError('literal held/active bias branch required')
        if self.active_bias:
            if not 0 < ph <= 1: raise ValueError('active bias factor must lie in (0,1]')
        elif ph != 1 or any(any(x for x in row) for row in qb):
            raise ValueError('held H18 BA covariance must use phi_b=1 and Q_BB=0')
        object.__setattr__(self,'F_AA',tuple(map(tuple,fa))); object.__setattr__(self,'Q_AA',tuple(map(tuple,qa)))
        object.__setattr__(self,'F_LL',tuple(map(tuple,fl))); object.__setattr__(self,'Q_LL',tuple(map(tuple,ql)))
        object.__setattr__(self,'Q_BB',tuple(map(tuple,qb))); object.__setattr__(self,'phi_b',ph)


def dense_FQ(blocks):
    Fm,Qm=M.zeros(N,N),M.zeros(N,N)
    for i in range(NA):
        for j in range(NA): Fm[i][j]=blocks.F_AA[i][j]; Qm[i][j]=blocks.Q_AA[i][j]
    for i in range(NL):
        for j in range(NL): Fm[OFF_L+i][OFF_L+j]=blocks.F_LL[i][j]; Qm[OFF_L+i][OFF_L+j]=blocks.Q_LL[i][j]
    for i in range(NB):
        Fm[OFF_B+i][OFF_B+i]=blocks.phi_b
        for j in range(NB): Qm[OFF_B+i][OFF_B+j]=blocks.Q_BB[i][j]
    return Fm,Qm


def shipping_block_successor(covariance, blocks):
    """Shipping block order, before pending-aw inflation/symmetrize/S service."""
    P0=M.mat(covariance,N,N)
    if P0 != M.transpose(P0): raise ValueError('symmetric predecessor covariance required')
    FA,QA=list(map(list,blocks.F_AA)),list(map(list,blocks.Q_AA)); FL,QL=list(map(list,blocks.F_LL)),list(map(list,blocks.Q_LL))
    QB,ph=list(map(list,blocks.Q_BB)),blocks.phi_b; out=[r[:] for r in P0]
    aa=M.plus(M.mm(M.mm(FA,[r[:NA] for r in P0[:NA]]),M.transpose(FA)),QA)
    for i in range(NA): out[i][:NA]=aa[i]
    ll=M.plus(M.mm(M.mm(FL,[r[OFF_L:OFF_L+NL] for r in P0[OFF_L:OFF_L+NL]]),M.transpose(FL)),QL)
    for i in range(NL): out[OFF_L+i][OFF_L:OFF_L+NL]=ll[i]
    al=M.mm(M.mm(FA,[r[OFF_L:OFF_L+NL] for r in P0[:NA]]),M.transpose(FL))
    for i in range(NA):
        for j in range(NL): out[i][OFF_L+j]=al[i][j]; out[OFF_L+j][i]=al[i][j]
    bb=M.plus(M.scaled([r[OFF_B:] for r in P0[OFF_B:]],ph*ph),QB)
    for i in range(NB): out[OFF_B+i][OFF_B:]=bb[i]
    ab=M.scaled(M.mm(FA,[r[OFF_B:] for r in P0[:NA]]),ph); lb=M.scaled(M.mm(FL,[r[OFF_B:] for r in P0[OFF_L:OFF_L+NL]]),ph)
    for i in range(NA):
        for j in range(NB): out[i][OFF_B+j]=ab[i][j]; out[OFF_B+j][i]=ab[i][j]
    for i in range(NL):
        for j in range(NB): out[OFF_L+i][OFF_B+j]=lb[i][j]; out[OFF_B+j][OFF_L+i]=lb[i][j]
    return out


def dense_successor(covariance, blocks):
    Fm,Qm=dense_FQ(blocks)
    return M.plus(M.mm(M.mm(Fm,covariance),M.transpose(Fm)),Qm)


def validate_shipping_identity(covariance, blocks):
    a=shipping_block_successor(covariance,blocks); b=dense_successor(covariance,blocks)
    if a != b: raise AssertionError('shipping block predictor differs from block-diagonal congruence')
    return a


def readiness():
    return {'shipping_covariance_block_identity':True,'full_F21_Q21_free_operands_removed_by_this_lemma':True,
            'linear_F_from_same_mean_coefficients':True,'correlated_and_independent_Q_LL_assemblies':True,
            'held_BA_branch_exactly_identity':True,'attitude_F_Q_primitives_source_attached':False,
            'Qaxis_analytic_primitives_source_attached':False,'active_BA_qd_scale_source_attached':False,
            'pending_aw_covariance_inflation_attached':False,'periodic_S_service_attached':False,
            'finite_precision_attached':False,'complete_word_finite_identity':False,'ALT_LIVE_PASS':False}
