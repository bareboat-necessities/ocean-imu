"""Square-root full-word feasibility engine; never a uniform certificate.

Enters V_next <= rho V + supply through D >= delta P_root^-1. If P0=C0 C0',
carry T=M C0 and compress root-whitened loss rows by QR, not normal equations.
The covariance factors represent the exact-real Joseph comparison. They do not
replace the shipping estimator's floating-point implementation.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import cholesky, qr, solve_triangular, svd


def _matrix(value):
    a=np.asarray(value,dtype=float)
    if a.ndim!=2 or not np.isfinite(a).all():
        raise ValueError('finite matrix required')
    return a


def _lower(a):
    """Triangular factor of a a', without forming that product."""
    return qr(a.T,mode='economic',check_finite=False)[1].T


def _symmetric(value):
    a=_matrix(value)
    if a.shape[0]!=a.shape[1] or not np.array_equal(a,a.T):
        raise ValueError('exactly symmetric supplied covariance required')
    return a


def psd_factor(q):
    """Explicit diagnostic adapter; never clip a negative eigenvalue to zero.

    Prefer supplied block factors. Eigenanalysis is only a fallback for test
    inputs and singular increments, and does not prove interval positivity.
    """
    q=_symmetric(q)
    if not np.any(q):
        return np.zeros((len(q),0))
    try:
        return cholesky(q,lower=True)
    except np.linalg.LinAlgError:
        w,v=np.linalg.eigh(q)
        if np.min(w)<0:
            raise ValueError('negative process spectrum; supply a verified factor') from None
        keep=w>0
        return v[:,keep]*np.sqrt(w[keep])


class FactorWord:
    """Stream one supplied regular A21 linear comparison in root coordinates.

    No event is accepted on behalf of the estimator. Rejected events are
    explicitly skipped. Covariance/gain/event schedules are realized inputs,
    not schedules of independently perturbed filters.
    """
    def __init__(self,root_covariance):
        p=_symmetric(root_covariance)
        self.c=cholesky(p,lower=True)
        self.root_factor=self.c.copy()
        self.t=self.c.copy()
        self.n=len(p)
        self.loss=np.zeros((0,self.n))
        self.counts={'prediction':0,'correction':0,'reset':0,'rejected':0}
        self.max_measurement_rows=0

    def _append(self,rows):
        if len(rows):
            self.loss=qr(np.vstack((self.loss,rows)),mode='economic',
                         check_finite=False)[1]

    def prediction(self,f,process_factor):
        f,u=_matrix(f),_matrix(process_factor)
        if f.shape!=(self.n,self.n) or u.shape[0]!=self.n:
            raise ValueError('incompatible prediction dimensions')
        # [F C,U]' = O [L';0]. The complement's first n columns factor
        # I-(L^-1 F C)'(L^-1 F C), without a precision subtraction.
        orth,r=qr(np.hstack((f@self.c,u)).T,mode='full',check_finite=False)
        cnew=r[:self.n,:].T
        if np.any(np.diag(cnew)==0):
            raise ValueError('prediction covariance is singular')
        z=solve_triangular(self.c,self.t,lower=True,check_finite=False)
        self._append(orth[:self.n,self.n:].T@z)
        self.t=f@self.t
        self.c=cnew
        self.counts['prediction']+=1

    def correction(self,h,r,*,accepted=True,innovation_increment=None):
        if not accepted:
            self.counts['rejected']+=1
            return None
        h,r=_matrix(h),_symmetric(r)
        m=len(h)
        if not 0<m<=3 or h.shape[1]!=self.n or r.shape!=(m,m):
            raise ValueError('regular correction requires one to three rows')
        if innovation_increment is not None:
            inc=_matrix(innovation_increment)
            if inc.shape!=(m,m):
                raise ValueError('innovation increment must match the observation')
            psd_factor(inc)  # refuse an indefinite safety adjustment
            r=r+inc
        vr=cholesky(r,lower=True)
        # Orthogonal array gives Ls, K Ls and posterior C simultaneously.
        # Its covariance is [[S,H P],[P H',P]], so the bottom Schur factor
        # is exactly the optimal-gain Joseph posterior with the effective R.
        a=np.block([[vr,h@self.c],[np.zeros((self.n,m)),self.c]])
        post=_lower(a)
        ls=post[:m,:m]
        gain=solve_triangular(ls.T,post[m:,:m].T,lower=False,
                              check_finite=False).T
        hm=h@self.t
        self._append(solve_triangular(ls,hm,lower=True,check_finite=False))
        self.t=self.t-gain@hm
        self.c=post[m:,m:]
        self.counts['correction']+=1
        self.max_measurement_rows=max(self.max_measurement_rows,m)
        return {'gain':gain,'innovation':ls@ls.T}

    def reset(self,g):
        g=_matrix(g)
        if g.shape!=(self.n,self.n) or np.linalg.matrix_rank(g)<self.n:
            raise ValueError('nonsingular same-dimension congruence required')
        self.c=_lower(g@self.c)
        self.t=g@self.t
        self.counts['reset']+=1

    def snapshot(self):
        z=solve_triangular(self.c,self.t,lower=True,check_finite=False)
        singular=svd(self.loss,compute_uv=False,check_finite=False)
        delta=float(singular[-1]**2) if len(singular)==self.n else 0.0
        rho=float(svd(z,compute_uv=False,check_finite=False)[0]**2)
        # Products below are diagnostic endpoints only, never the accumulation.
        residual=z.T@z+self.loss.T@self.loss-np.eye(self.n)
        return {'delta_diagnostic':delta,'rho_diagnostic':rho,
                'identity_residual_fro':float(np.linalg.norm(residual,'fro')),
                'compressed_loss_rows':len(self.loss),
                'max_measurement_rows':self.max_measurement_rows,
                'events':dict(self.counts),'source_uniform_verified':False,
                'arithmetic_enclosure_verified':False,'theorem_closed':False}


def run_word(root_covariance,events):
    word=FactorWord(root_covariance)
    for event in events:
        kind=event['kind']
        if kind=='prediction':
            u=event['U'] if 'U' in event else psd_factor(event['Q'])
            word.prediction(event['F'],u)
        elif kind=='correction':
            word.correction(event['H'],event['R'],accepted=event.get('accepted',True),
                            innovation_increment=event.get('innovation_increment'))
        elif kind=='reset':
            word.reset(event['G'])
        else:
            raise ValueError('unsupported operation; do not silently omit it')
    return word


def nuisance_residual(factor,heading_columns,*,rank_tolerance):
    """Diagnostic range projection; numerical rank is explicitly unverified.

    The exact companion handles singular nuisance columns without tolerances.
    """
    a=_matrix(factor)
    h=list(heading_columns)
    if len(set(h))!=len(h) or any(i<0 or i>=a.shape[1] for i in h):
        raise ValueError('distinct valid heading columns required')
    if not np.isfinite(rank_tolerance) or rank_tolerance<0:
        raise ValueError('nonnegative finite rank tolerance required')
    n=[i for i in range(a.shape[1]) if i not in h]
    u,s,_=svd(a[:,n],full_matrices=True,check_finite=False)
    rank=int(np.sum(s>rank_tolerance))
    residual=u[:,rank:].T@a[:,h]
    return {'residual_factor':residual,'nuisance_rank_diagnostic':rank,
            'rank_tolerance':rank_tolerance,'rank_certified':False}
